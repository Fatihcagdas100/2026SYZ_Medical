"""Varyant-düzeyi ve gene-disjoint ayrımlar arasındaki rol geçişlerini denetle."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    from scripts.analyze_oof_calibration import sha256_file
except ModuleNotFoundError:  # Doğrudan `python scripts/...` çalıştırması için.
    from analyze_oof_calibration import sha256_file  # type: ignore[no-redef]


DEFAULT_OLD_SPLIT = Path("data/processed/baseline_stratified_split.tsv")
DEFAULT_NEW_SPLIT = Path("data/processed/baseline_gene_disjoint_split.tsv")
DEFAULT_REPORT_DIR = Path("reports/generated/split_transition_audit")
SPLIT_ORDER = ("train", "validation", "test")


def load_split(path: Path, *, label: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{label} ayrım dosyası bulunamadı: {path}")
    table = pd.read_csv(
        path,
        sep="\t",
        usecols=["VariationID", "target", "split"],
        dtype={"VariationID": "string", "target": "Int8", "split": "string"},
    )
    if table.isna().any().any():
        raise ValueError(f"{label} ayrımın gerekli alanlarında boş değer bulundu.")
    if table["VariationID"].duplicated().any():
        raise ValueError(f"{label} ayrımda tekrarlanan VariationID bulundu.")
    if set(table["split"].unique()) != set(SPLIT_ORDER):
        raise ValueError(f"{label} ayrım train, validation ve test içermelidir.")
    if set(table["target"].astype(int).unique()) != {0, 1}:
        raise ValueError(f"{label} ayrım iki hedef sınıfını da içermelidir.")
    return table


def merge_splits(old_split_path: Path, new_split_path: Path) -> pd.DataFrame:
    old = load_split(old_split_path, label="Eski").rename(
        columns={"target": "old_target", "split": "old_split"}
    )
    new = load_split(new_split_path, label="Yeni").rename(
        columns={"target": "new_target", "split": "new_split"}
    )
    old_ids = set(old["VariationID"])
    new_ids = set(new["VariationID"])
    if old_ids != new_ids:
        difference = old_ids.symmetric_difference(new_ids)
        raise ValueError(f"Ayrımların VariationID kümeleri eşleşmiyor: {min(difference)}")
    merged = old.merge(new, on="VariationID", how="inner", validate="one_to_one")
    mismatch = merged["old_target"].astype(int) != merged["new_target"].astype(int)
    if mismatch.any():
        variation_id = merged.loc[mismatch, "VariationID"].iloc[0]
        raise ValueError(f"Ayrımlar arasında hedef uyuşmazlığı: {variation_id}")
    return merged.drop(columns="new_target").rename(columns={"old_target": "target"})


def build_transition_table(merged: pd.DataFrame) -> pd.DataFrame:
    counts = (
        merged.groupby(["new_split", "old_split"], observed=True)
        .size()
        .rename("count")
        .reset_index()
    )
    counts["new_split_total"] = counts.groupby("new_split")["count"].transform("sum")
    counts["percentage_within_new_split"] = (
        100 * counts["count"] / counts["new_split_total"]
    )
    counts["new_split"] = pd.Categorical(
        counts["new_split"], categories=SPLIT_ORDER, ordered=True
    )
    counts["old_split"] = pd.Categorical(
        counts["old_split"], categories=SPLIT_ORDER, ordered=True
    )
    return counts.sort_values(["new_split", "old_split"]).reset_index(drop=True)


def build_final_test_roles(merged: pd.DataFrame) -> pd.DataFrame:
    final_test = merged.loc[merged["new_split"] == "test"]
    counts = final_test["old_split"].value_counts().reindex(SPLIT_ORDER, fill_value=0)
    result = counts.rename_axis("previous_split").rename("count").reset_index()
    result["percentage_of_final_test"] = 100 * result["count"] / len(final_test)
    result["previous_development_role"] = result["previous_split"].isin(
        ["train", "validation"]
    )
    return result


def build_transition_by_target(merged: pd.DataFrame) -> pd.DataFrame:
    result = (
        merged.groupby(["target", "new_split", "old_split"], observed=True)
        .size()
        .rename("count")
        .reset_index()
    )
    result["target_label"] = result["target"].map({0: "Benign", 1: "Pathogenic"})
    totals = result.groupby(["target", "new_split"])["count"].transform("sum")
    result["percentage_within_target_and_new_split"] = 100 * result["count"] / totals
    return result.sort_values(["target", "new_split", "old_split"]).reset_index(drop=True)


def plot_transition_heatmap(merged: pd.DataFrame, output_path: Path) -> None:
    matrix = pd.crosstab(merged["new_split"], merged["old_split"]).reindex(
        index=SPLIT_ORDER, columns=SPLIT_ORDER, fill_value=0
    )
    figure, axis = plt.subplots(figsize=(7.5, 5.5))
    image = axis.imshow(matrix.to_numpy(), cmap="Blues")
    for row in range(len(SPLIT_ORDER)):
        for column in range(len(SPLIT_ORDER)):
            count = int(matrix.iloc[row, column])
            percentage = 100 * count / int(matrix.iloc[row].sum())
            axis.text(
                column,
                row,
                f"{count:,}\n%{percentage:.2f}",
                ha="center",
                va="center",
                color="white" if image.norm(count) > 0.55 else "black",
            )
    axis.set_xticks(np.arange(3), SPLIT_ORDER)
    axis.set_yticks(np.arange(3), SPLIT_ORDER)
    axis.set_xlabel("Önceki varyant-düzeyi rol")
    axis.set_ylabel("Son gene-disjoint rol")
    axis.set_title("Veri ayrımları arasında VariationID rol geçişi")
    figure.colorbar(image, ax=axis, label="Varyant sayısı")
    figure.tight_layout()
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)


def format_report(manifest: dict[str, object]) -> str:
    exposure = manifest["final_test_prior_exposure"]
    assert isinstance(exposure, dict)
    return f"""# Veri ayrımı geçiş denetimi

## Ana bulgu

Son gene-disjoint test bölümündeki {exposure['final_test_rows']:,} varyantın
{exposure['previous_train_rows']:,} tanesi önceki varyant-düzeyi train,
{exposure['previous_validation_rows']:,} tanesi önceki validation ve
{exposure['previous_test_rows']:,} tanesi önceki test rolündeydi.

Önceki train veya validation bölümlerinde yer alan toplam
{exposure['previous_development_rows']:,} varyant, son testin
`%{exposure['previous_development_percentage']:.2f}`'sini oluşturur.

## Metodolojik sonuç

Son gene-disjoint model bu test varyantlarıyla doğrudan eğitilmemiş ve test
gen grupları son modelin train/validation gruplarıyla çakışmamıştır. Bununla
birlikte varyantların büyük bölümü daha önceki deneylerin train veya validation
rollerinde bulunduğu için bu değerlendirme geliştirme süreci açısından tamamen
görülmemiş, bağımsız bir nihai test olarak sunulamaz.

Mevcut sonuç **iç gene-disjoint değerlendirme** olarak raporlanmalıdır.
Bağımsız nihai performans iddiası için geliştirme sürecinde hiç kullanılmamış
harici bir veri kaynağı veya prospektif holdout gerekir.
"""


def audit_split_transition(
    old_split_path: Path,
    new_split_path: Path,
    report_dir: Path,
) -> dict[str, object]:
    merged = merge_splits(old_split_path, new_split_path)
    transition = build_transition_table(merged)
    final_test_roles = build_final_test_roles(merged)
    by_target = build_transition_by_target(merged)

    report_dir.mkdir(parents=True, exist_ok=True)
    transition_path = report_dir / "split_transition_counts.csv"
    roles_path = report_dir / "final_test_previous_roles.csv"
    by_target_path = report_dir / "split_transition_by_target.csv"
    plot_path = report_dir / "split_transition_heatmap.png"
    report_path = report_dir / "split_transition_report.md"
    transition.to_csv(transition_path, index=False, float_format="%.10f")
    final_test_roles.to_csv(roles_path, index=False, float_format="%.10f")
    by_target.to_csv(by_target_path, index=False, float_format="%.10f")
    plot_transition_heatmap(merged, plot_path)

    role_counts = final_test_roles.set_index("previous_split")["count"]
    previous_development_rows = int(role_counts["train"] + role_counts["validation"])
    final_test_rows = int((merged["new_split"] == "test").sum())
    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "old_variant_level_split": {
                "path": str(old_split_path.resolve()),
                "sha256": sha256_file(old_split_path),
            },
            "new_gene_disjoint_split": {
                "path": str(new_split_path.resolve()),
                "sha256": sha256_file(new_split_path),
            },
        },
        "integrity": {
            "rows": len(merged),
            "unique_variation_ids": int(merged["VariationID"].nunique()),
            "id_sets_match": True,
            "target_mismatches": 0,
        },
        "final_test_prior_exposure": {
            "final_test_rows": final_test_rows,
            "previous_train_rows": int(role_counts["train"]),
            "previous_validation_rows": int(role_counts["validation"]),
            "previous_test_rows": int(role_counts["test"]),
            "previous_development_rows": previous_development_rows,
            "previous_development_percentage": (
                100 * previous_development_rows / final_test_rows
            ),
        },
        "interpretation": {
            "direct_training_leakage_into_final_gene_disjoint_model": False,
            "fully_process_unseen_final_test": False,
            "recommended_label": "internal_gene_disjoint_evaluation",
            "external_evaluation_required": True,
        },
        "outputs": {
            "transition_counts": str(transition_path.resolve()),
            "final_test_previous_roles": str(roles_path.resolve()),
            "transition_by_target": str(by_target_path.resolve()),
            "heatmap": str(plot_path.resolve()),
            "report": str(report_path.resolve()),
        },
    }
    report_path.write_text(format_report(manifest), encoding="utf-8")
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("Veri ayrımı geçiş denetimi tamamlandı")
    print(f"Son gene-disjoint test: {final_test_rows:,}")
    print(
        f"Önceki train/validation rolü: {previous_development_rows:,} "
        f"(%{100 * previous_development_rows / final_test_rows:.2f})"
    )
    print(f"Önceki test rolü: {int(role_counts['test']):,}")
    print("Önerilen adlandırma: iç gene-disjoint değerlendirme")
    print(f"Rapor: {report_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Eski ve yeni veri ayrımlarındaki VariationID rollerini denetler."
    )
    parser.add_argument("--old-split", type=Path, default=DEFAULT_OLD_SPLIT)
    parser.add_argument("--new-split", type=Path, default=DEFAULT_NEW_SPLIT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    audit_split_transition(args.old_split, args.new_split, args.report_dir)


if __name__ == "__main__":
    main()
