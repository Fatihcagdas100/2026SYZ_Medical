"""Missense özellik kazanımını train-only gene-disjoint OOF ile doğrula."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

try:
    from scripts.cross_validate_gene_disjoint import cross_validate
    from scripts.run_hgvs_ablation import (
        calculate_subgroup_metrics,
        load_validation_annotations,
    )
    from scripts.train_logistic_baseline import sha256_file
except ModuleNotFoundError:  # Doğrudan `python scripts/...` çalıştırması için.
    from cross_validate_gene_disjoint import cross_validate  # type: ignore[no-redef]
    from run_hgvs_ablation import (  # type: ignore[no-redef]
        calculate_subgroup_metrics,
        load_validation_annotations,
    )
    from train_logistic_baseline import sha256_file  # type: ignore[no-redef]


DEFAULT_FEATURES = Path("data/processed/clinvar_missense_features.tsv.gz")
DEFAULT_SPLIT = Path("data/processed/baseline_gene_disjoint_split.tsv")
DEFAULT_REPORT_DIR = Path("reports/generated/missense_gene_disjoint_cv")
DEFAULT_SEED = 2027
DEFAULT_N_SPLITS = 5
DEFAULT_THRESHOLD = 0.5
DEFAULT_CHUNK_SIZE = 100_000
PROFILES = (
    ("hgvs_reference", "hgvs", "Tam HGVS referans"),
    ("hgvs_plus_missense_biology", "hgvs_missense", "HGVS + missense biyoloji"),
)
METRICS = (
    "roc_auc",
    "average_precision",
    "precision",
    "recall_sensitivity",
    "specificity",
    "f1",
)


def _assert_identical_folds(reference_path: Path, enriched_path: Path) -> None:
    reference = pd.read_csv(reference_path, sep="\t", usecols=["VariationID", "fold"])
    enriched = pd.read_csv(enriched_path, sep="\t", usecols=["VariationID", "fold"])
    merged = reference.merge(
        enriched,
        on="VariationID",
        how="outer",
        suffixes=("_reference", "_enriched"),
        indicator=True,
        validate="one_to_one",
    )
    if not merged["_merge"].eq("both").all() or not (
        merged["fold_reference"] == merged["fold_enriched"]
    ).all():
        raise RuntimeError("İki profil aynı OOF satırları ve fold atamalarını kullanmadı.")


def build_oof_subgroup_table(
    predictions_path: Path,
    features_path: Path,
    *,
    threshold: float,
    chunk_size: int,
) -> pd.DataFrame:
    """Sadece train kaynaklı OOF tahminlerinde HGVS alt gruplarını ölç."""
    predictions = pd.read_csv(
        predictions_path,
        sep="\t",
        dtype={"VariationID": "string", "source_partition": "string"},
    )
    if set(predictions["source_partition"].dropna()) != {"train"}:
        raise ValueError("Alt grup analizi yalnızca train kaynaklı OOF tahminlerini kabul eder.")
    train_ids = set(predictions["VariationID"].astype(str))
    annotations = load_validation_annotations(
        features_path, train_ids, chunk_size=chunk_size
    )
    merged = predictions.merge(
        annotations, on="VariationID", how="left", validate="one_to_one"
    )
    groups = (
        ("all_variants", pd.Series(True, index=merged.index)),
        ("missense", merged["protein_effect"].eq("missense")),
        ("synonymous", merged["protein_effect"].eq("synonymous")),
        ("frameshift", merged["protein_effect"].eq("frameshift")),
        (
            "stop_gained_or_changed",
            merged["protein_effect"].eq("stop_gained_or_changed"),
        ),
        (
            "protein_annotation_present",
            merged["has_protein_annotation"].eq(1),
        ),
        (
            "protein_annotation_absent",
            merged["has_protein_annotation"].eq(0),
        ),
    )
    return pd.DataFrame(
        [
            calculate_subgroup_metrics(
                merged.loc[mask], group_name=name, threshold=threshold
            )
            for name, mask in groups
        ]
    )


def run_missense_gene_disjoint_cv(
    features_path: Path,
    split_path: Path,
    report_dir: Path,
    *,
    seed: int = DEFAULT_SEED,
    n_splits: int = DEFAULT_N_SPLITS,
    threshold: float = DEFAULT_THRESHOLD,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, object]:
    report_dir.mkdir(parents=True, exist_ok=True)
    manifests: dict[str, dict[str, object]] = {}
    subgroup_tables: list[pd.DataFrame] = []
    fold_tables: list[pd.DataFrame] = []

    for key, feature_profile, label in PROFILES:
        profile_dir = report_dir / key
        manifest = cross_validate(
            features_path,
            split_path,
            profile_dir,
            seed=seed,
            n_splits=n_splits,
            threshold=threshold,
            chunk_size=chunk_size,
            feature_profile=feature_profile,
        )
        manifests[key] = manifest
        subgroup = build_oof_subgroup_table(
            profile_dir / "train_out_of_fold_predictions.tsv.gz",
            features_path,
            threshold=threshold,
            chunk_size=chunk_size,
        )
        subgroup.insert(0, "profile", key)
        subgroup.insert(1, "label", label)
        subgroup_tables.append(subgroup)

        fold_metrics = pd.read_csv(profile_dir / "fold_metrics.csv")
        fold_metrics.insert(0, "profile", key)
        fold_tables.append(fold_metrics)

    reference_key, _, _ = PROFILES[0]
    enriched_key, _, _ = PROFILES[1]
    _assert_identical_folds(
        report_dir / reference_key / "train_out_of_fold_predictions.tsv.gz",
        report_dir / enriched_key / "train_out_of_fold_predictions.tsv.gz",
    )

    subgroups = pd.concat(subgroup_tables, ignore_index=True)
    folds = pd.concat(fold_tables, ignore_index=True)
    subgroup_path = report_dir / "oof_subgroup_comparison.csv"
    fold_path = report_dir / "fold_model_comparison.csv"
    subgroups.to_csv(subgroup_path, index=False, float_format="%.10f")
    folds.to_csv(fold_path, index=False, float_format="%.10f")

    reference = subgroups.loc[subgroups["profile"] == reference_key].set_index("group")
    enriched = subgroups.loc[subgroups["profile"] == enriched_key].set_index("group")
    delta_rows: list[dict[str, object]] = []
    for group in reference.index:
        for metric in METRICS:
            delta_rows.append(
                {
                    "group": group,
                    "metric": metric,
                    "enriched_minus_reference": float(
                        enriched.loc[group, metric] - reference.loc[group, metric]
                    ),
                }
            )
    deltas = pd.DataFrame(delta_rows)
    delta_path = report_dir / "oof_metric_deltas.csv"
    deltas.to_csv(delta_path, index=False, float_format="%.10f")

    figure, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for axis, group in zip(axes, ("all_variants", "missense"), strict=True):
        subset = subgroups.loc[subgroups["group"] == group].set_index("profile")
        metrics = ("roc_auc", "average_precision", "recall_sensitivity")
        positions = range(len(metrics))
        width = 0.36
        axis.bar(
            [position - width / 2 for position in positions],
            subset.loc[reference_key, list(metrics)],
            width,
            label="HGVS referans",
        )
        axis.bar(
            [position + width / 2 for position in positions],
            subset.loc[enriched_key, list(metrics)],
            width,
            label="+ missense biyoloji",
        )
        axis.set_xticks(list(positions), ["ROC-AUC", "AP", "Recall"])
        axis.set_title(group)
        axis.set_ylim(0, 1.02)
        axis.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Birleştirilmiş OOF skoru")
    axes[1].legend(loc="lower right")
    figure.suptitle("Train-only gene-disjoint OOF missense karşılaştırması")
    figure.tight_layout()
    plot_path = report_dir / "oof_missense_comparison.png"
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)

    report_path = report_dir / "missense_gene_disjoint_cv.md"
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write("# Train-only gene-disjoint OOF missense karşılaştırması\n\n")
        handle.write(
            "Her iki profil yalnızca sabit `train` bölümünde aynı "
            "StratifiedGroupKFold atamalarıyla çalıştırıldı. Sabit validation ve "
            "kilitli test yüklenmedi.\n\n"
        )
        handle.write("## Birleştirilmiş OOF sonuçları\n\n")
        handle.write("| Grup | Profil | n | ROC-AUC | AP | Precision | Recall | Özgüllük | F1 |\n")
        handle.write("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for group in ("all_variants", "missense"):
            for key, _, label in PROFILES:
                row = subgroups.loc[
                    (subgroups["group"] == group) & (subgroups["profile"] == key)
                ].iloc[0]
                handle.write(
                    f"| {group} | {label} | {int(row['row_count']):,} | "
                    f"{row['roc_auc']:.4f} | {row['average_precision']:.4f} | "
                    f"{row['precision']:.4f} | {row['recall_sensitivity']:.4f} | "
                    f"{row['specificity']:.4f} | {row['f1']:.4f} |\n"
                )
        handle.write("\nBu iç geliştirme doğrulamasıdır; dış veri veya klinik geçerlilik kanıtı değildir.\n")

    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "features": {"path": str(features_path.resolve()), "sha256": sha256_file(features_path)},
            "gene_disjoint_split": {"path": str(split_path.resolve()), "sha256": sha256_file(split_path)},
        },
        "data_policy": {
            "source_partition": "train",
            "fixed_validation_used": False,
            "test_used": False,
            "same_seed_n_splits_threshold_and_fold_assignments": True,
        },
        "profiles": manifests,
        "outputs": {
            "subgroup_comparison": str(subgroup_path.resolve()),
            "fold_comparison": str(fold_path.resolve()),
            "metric_deltas": str(delta_path.resolve()),
            "plot": str(plot_path.resolve()),
            "report": str(report_path.resolve()),
        },
    }
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    missense_delta = deltas.loc[
        deltas["group"].eq("missense") & deltas["metric"].eq("roc_auc"),
        "enriched_minus_reference",
    ].iloc[0]
    print("Train-only gene-disjoint missense karşılaştırması tamamlandı")
    print(f"Missense OOF ROC-AUC farkı: {missense_delta:+.4f}")
    print(f"Rapor: {report_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="HGVS ve HGVS + missense profillerini train-only OOF ile karşılaştırır."
    )
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--split", type=Path, default=DEFAULT_SPLIT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--n-splits", type=int, default=DEFAULT_N_SPLITS)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    run_missense_gene_disjoint_cv(
        args.features,
        args.split,
        args.report_dir,
        seed=args.seed,
        n_splits=args.n_splits,
        threshold=args.threshold,
        chunk_size=args.chunk_size,
    )


if __name__ == "__main__":
    main()
