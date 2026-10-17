"""Missense biyolojik özelliklerinin kontrollü validation ablasyonunu çalıştır."""

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
    from scripts.run_hgvs_ablation import build_subgroup_table
    from scripts.train_logistic_baseline import (
        MISSENSE_CATEGORICAL_FEATURES,
        MISSENSE_NUMERIC_FEATURES,
        sha256_file,
        train_baselines,
    )
except ModuleNotFoundError:  # Doğrudan ``python scripts/...`` çalıştırması.
    from run_hgvs_ablation import build_subgroup_table
    from train_logistic_baseline import (
        MISSENSE_CATEGORICAL_FEATURES,
        MISSENSE_NUMERIC_FEATURES,
        sha256_file,
        train_baselines,
    )


DEFAULT_FEATURES = Path("data/processed/clinvar_missense_features.tsv.gz")
DEFAULT_SPLIT = Path("data/processed/baseline_gene_disjoint_split.tsv")
DEFAULT_REPORT_DIR = Path("reports/generated/missense_feature_ablation")
DEFAULT_ARTIFACT_DIR = Path("artifacts/missense_feature_ablation")
DEFAULT_SEED = 2027
DEFAULT_THRESHOLD = 0.5
DEFAULT_CHUNK_SIZE = 100_000
PRIMARY_MODEL = "without_missingness_flag"
PROFILES = (
    ("hgvs_reference", "hgvs", "Tam HGVS referans"),
    (
        "hgvs_plus_missense_biology",
        "hgvs_missense",
        "HGVS + missense biyolojik özellikler",
    ),
)
DELTA_METRICS = (
    "roc_auc",
    "average_precision",
    "precision",
    "recall_sensitivity",
    "specificity",
    "f1",
)


def run_missense_feature_ablation(
    features_path: Path,
    split_path: Path,
    report_dir: Path,
    artifact_dir: Path,
    *,
    seed: int = DEFAULT_SEED,
    threshold: float = DEFAULT_THRESHOLD,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, object]:
    report_dir.mkdir(parents=True, exist_ok=True)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    overall_rows: list[dict[str, object]] = []
    subgroup_tables: list[pd.DataFrame] = []
    manifests: dict[str, object] = {}

    for profile_key, training_profile, label in PROFILES:
        profile_report_dir = report_dir / "models" / profile_key
        profile_artifact_dir = artifact_dir / profile_key
        profile_manifest = train_baselines(
            features_path,
            split_path,
            profile_report_dir,
            profile_artifact_dir,
            seed=seed,
            threshold=threshold,
            chunk_size=chunk_size,
            feature_profile=training_profile,
        )
        if profile_manifest["data_policy"]["test_partition"] != (
            "locked_not_loaded_into_development_dataframe"
        ):
            raise ValueError("Missense ablasyonu test bölümünü kilitli tutmadı.")
        metrics = profile_manifest["validation_metrics"][PRIMARY_MODEL]
        overall_rows.append(
            {
                "profile": profile_key,
                "label": label,
                "training_feature_profile": training_profile,
                **metrics,
            }
        )
        subgroup = build_subgroup_table(
            profile_report_dir / "validation_predictions.tsv.gz",
            features_path,
            threshold=threshold,
            chunk_size=chunk_size,
        )
        subgroup.insert(0, "profile", profile_key)
        subgroup_tables.append(subgroup)
        manifests[profile_key] = profile_manifest

    overall = pd.DataFrame(overall_rows)
    subgroups = pd.concat(subgroup_tables, ignore_index=True)
    overall_path = report_dir / "validation_model_comparison.csv"
    subgroup_path = report_dir / "validation_subgroup_comparison.csv"
    overall.to_csv(overall_path, index=False, float_format="%.10f")
    subgroups.to_csv(subgroup_path, index=False, float_format="%.10f")

    reference = subgroups.loc[subgroups["profile"] == PROFILES[0][0]].set_index(
        "group"
    )
    enriched = subgroups.loc[subgroups["profile"] == PROFILES[1][0]].set_index(
        "group"
    )
    delta_rows: list[dict[str, object]] = []
    for group_name in reference.index:
        for metric in DELTA_METRICS:
            delta_rows.append(
                {
                    "group": group_name,
                    "metric": metric,
                    "enriched_minus_reference": float(
                        enriched.loc[group_name, metric]
                        - reference.loc[group_name, metric]
                    ),
                }
            )
    deltas = pd.DataFrame(delta_rows)
    delta_path = report_dir / "validation_metric_deltas.csv"
    deltas.to_csv(delta_path, index=False, float_format="%.10f")

    enriched_coefficients = pd.read_csv(
        report_dir / "models" / PROFILES[1][0] / "coefficients.csv"
    )
    feature_prefixes = (*MISSENSE_CATEGORICAL_FEATURES, *MISSENSE_NUMERIC_FEATURES)
    targeted_coefficients = enriched_coefficients.loc[
        enriched_coefficients["model_name"].eq(PRIMARY_MODEL)
        & enriched_coefficients["transformed_feature"].str.startswith(
            feature_prefixes
        )
    ].copy()
    targeted_coefficients = targeted_coefficients.sort_values(
        "absolute_coefficient", ascending=False
    )
    coefficient_path = report_dir / "missense_feature_coefficients.csv"
    targeted_coefficients.to_csv(coefficient_path, index=False, float_format="%.10f")

    plot_groups = ("all_variants", "missense")
    plot_metrics = ("roc_auc", "average_precision", "recall_sensitivity")
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for axis, group_name in zip(axes, plot_groups, strict=True):
        subset = subgroups.loc[subgroups["group"] == group_name].set_index("profile")
        positions = range(len(plot_metrics))
        width = 0.36
        axis.bar(
            [position - width / 2 for position in positions],
            subset.loc[PROFILES[0][0], list(plot_metrics)],
            width,
            label="HGVS referans",
        )
        axis.bar(
            [position + width / 2 for position in positions],
            subset.loc[PROFILES[1][0], list(plot_metrics)],
            width,
            label="+ missense biyoloji",
        )
        axis.set_xticks(list(positions), ["ROC-AUC", "AP", "Recall"])
        axis.set_title(group_name)
        axis.set_ylim(0, 1.02)
        axis.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Validation skoru")
    axes[1].legend(loc="lower right")
    figure.suptitle("Missense biyolojik özellik ablasyonu")
    figure.tight_layout()
    plot_path = report_dir / "missense_feature_ablation.png"
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)

    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "features": {
                "path": str(features_path.resolve()),
                "sha256": sha256_file(features_path),
            },
            "split": {
                "path": str(split_path.resolve()),
                "sha256": sha256_file(split_path),
            },
        },
        "data_policy": {
            "training_partition": "train",
            "evaluation_partition": "validation",
            "test_partition": "locked_not_loaded_into_development_dataframe",
            "threshold": threshold,
            "same_estimator_split_seed_and_threshold": True,
            "hyperparameter_search": False,
        },
        "profiles": [
            {
                "profile": key,
                "training_feature_profile": training_profile,
                "label": label,
                "feature_columns": manifests[key]["models"][PRIMARY_MODEL][
                    "feature_columns"
                ],
            }
            for key, training_profile, label in PROFILES
        ],
        "reports": {
            "model_comparison": str(overall_path.resolve()),
            "subgroup_comparison": str(subgroup_path.resolve()),
            "metric_deltas": str(delta_path.resolve()),
            "missense_feature_coefficients": str(coefficient_path.resolve()),
            "plot": str(plot_path.resolve()),
        },
    }

    report_path = report_dir / "missense_feature_ablation.md"
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write("# Missense biyolojik özellik ablasyonu\n\n")
        handle.write(
            "İki lojistik model aynı gene-disjoint train/validation ayrımı, seed "
            "ve sabit 0.5 eşikle karşılaştırıldı. Test yüklenmedi.\n\n"
        )
        handle.write("| Profil | ROC-AUC | AP | Precision | Recall | Özgüllük | F1 |\n")
        handle.write("| --- | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for row in overall.to_dict("records"):
            handle.write(
                f"| {row['label']} | {row['roc_auc']:.4f} | "
                f"{row['average_precision']:.4f} | {row['precision']:.4f} | "
                f"{row['recall_sensitivity']:.4f} | {row['specificity']:.4f} | "
                f"{row['f1']:.4f} |\n"
            )
        handle.write("\n## Missense alt grubu\n\n")
        handle.write("| Profil | n | ROC-AUC | AP | Precision | Recall | Özgüllük | F1 |\n")
        handle.write("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for row in subgroups.loc[subgroups["group"] == "missense"].to_dict(
            "records"
        ):
            handle.write(
                f"| {row['profile']} | {row['row_count']:,} | {row['roc_auc']:.4f} | "
                f"{row['average_precision']:.4f} | {row['precision']:.4f} | "
                f"{row['recall_sensitivity']:.4f} | {row['specificity']:.4f} | "
                f"{row['f1']:.4f} |\n"
            )
        handle.write(
            "\nBu sonuçlar validation bulgusudur; klinik geçerlilik veya dış veri "
            "genellemesi anlamına gelmez.\n"
        )
    manifest["reports"]["report"] = str(report_path.resolve())
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("Missense biyolojik özellik ablasyonu tamamlandı")
    for row in overall.to_dict("records"):
        print(
            f"{row['profile']}: ROC-AUC={row['roc_auc']:.4f}, "
            f"AP={row['average_precision']:.4f}, F1={row['f1']:.4f}"
        )
    for row in subgroups.loc[subgroups["group"] == "missense"].to_dict("records"):
        print(
            f"missense/{row['profile']}: ROC-AUC={row['roc_auc']:.4f}, "
            f"AP={row['average_precision']:.4f}, recall={row['recall_sensitivity']:.4f}"
        )
    print(f"Rapor: {report_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Missense biyolojik özelliklerinin validation katkısını ölçer."
    )
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--split", type=Path, default=DEFAULT_SPLIT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--artifact-dir", type=Path, default=DEFAULT_ARTIFACT_DIR)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    run_missense_feature_ablation(
        args.features,
        args.split,
        args.report_dir,
        args.artifact_dir,
        seed=args.seed,
        threshold=args.threshold,
        chunk_size=args.chunk_size,
    )


if __name__ == "__main__":
    main()
