"""Gene-disjoint validation üzerinde HGVS özellik ablasyonu ve alt grup analizi."""

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
from sklearn.metrics import (
    average_precision_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)

try:
    from scripts.train_logistic_baseline import sha256_file, train_baselines
except ModuleNotFoundError:  # Doğrudan ``python scripts/...`` çalıştırması.
    from train_logistic_baseline import sha256_file, train_baselines


DEFAULT_FEATURES = Path("data/processed/clinvar_hgvs_features.tsv.gz")
DEFAULT_SPLIT = Path("data/processed/baseline_gene_disjoint_split.tsv")
DEFAULT_REPORT_DIR = Path("reports/generated/hgvs_ablation")
DEFAULT_ARTIFACT_DIR = Path("artifacts/hgvs_ablation")
DEFAULT_SEED = 2027
DEFAULT_CHUNK_SIZE = 100_000
PRIMARY_MODEL = "without_missingness_flag"
PROFILES = (
    ("structural", "baseline", "Yapısal"),
    (
        "structural_plus_hgvs_biological",
        "hgvs_biological",
        "Yapısal + HGVS biyolojik sonuç",
    ),
    (
        "structural_plus_hgvs_biological_and_representation",
        "hgvs",
        "Yapısal + HGVS biyolojik sonuç + gösterim",
    ),
)


def calculate_subgroup_metrics(
    frame: pd.DataFrame,
    *,
    group_name: str,
    threshold: float,
) -> dict[str, object]:
    targets = frame["target"].astype(int)
    probabilities = frame["probability_pathogenic"].astype(float)
    if frame.empty:
        return {
            "group": group_name,
            "row_count": 0,
            "negative_count": 0,
            "positive_count": 0,
            "positive_prevalence": float("nan"),
            "precision": float("nan"),
            "recall_sensitivity": float("nan"),
            "specificity": float("nan"),
            "f1": float("nan"),
            "roc_auc": float("nan"),
            "average_precision": float("nan"),
            "true_negatives": 0,
            "false_positives": 0,
            "false_negatives": 0,
            "true_positives": 0,
        }
    predictions = (probabilities >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(targets, predictions, labels=[0, 1]).ravel()
    both_classes = targets.nunique() == 2
    return {
        "group": group_name,
        "row_count": len(frame),
        "negative_count": int((targets == 0).sum()),
        "positive_count": int((targets == 1).sum()),
        "positive_prevalence": float(targets.mean()) if len(targets) else float("nan"),
        "precision": precision_score(targets, predictions, zero_division=0),
        "recall_sensitivity": recall_score(targets, predictions, zero_division=0),
        "specificity": tn / (tn + fp) if tn + fp else float("nan"),
        "f1": f1_score(targets, predictions, zero_division=0),
        "roc_auc": roc_auc_score(targets, probabilities) if both_classes else float("nan"),
        "average_precision": (
            average_precision_score(targets, probabilities)
            if int((targets == 1).sum()) > 0
            else float("nan")
        ),
        "true_negatives": int(tn),
        "false_positives": int(fp),
        "false_negatives": int(fn),
        "true_positives": int(tp),
    }


def load_validation_annotations(
    features_path: Path,
    validation_ids: set[str],
    *,
    chunk_size: int,
) -> pd.DataFrame:
    retained: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        features_path,
        sep="\t",
        usecols=["VariationID", "protein_effect", "has_protein_annotation"],
        dtype={
            "VariationID": "string",
            "protein_effect": "string",
            "has_protein_annotation": "Int8",
        },
        chunksize=chunk_size,
    ):
        selected = chunk.loc[chunk["VariationID"].isin(validation_ids)]
        if not selected.empty:
            retained.append(selected)
    annotations = pd.concat(retained, ignore_index=True)
    if annotations["VariationID"].duplicated().any():
        raise ValueError("Validation anotasyonlarında tekrarlanan VariationID bulundu.")
    if set(annotations["VariationID"].astype(str)) != validation_ids:
        raise ValueError("Tüm validation kimlikleri için HGVS anotasyonu bulunamadı.")
    return annotations


def build_subgroup_table(
    predictions_path: Path,
    features_path: Path,
    *,
    threshold: float,
    chunk_size: int,
) -> pd.DataFrame:
    predictions = pd.read_csv(
        predictions_path,
        sep="\t",
        dtype={"VariationID": "string", "model_name": "string", "split": "string"},
    )
    predictions = predictions.loc[predictions["model_name"] == PRIMARY_MODEL].copy()
    if set(predictions["split"].dropna()) != {"validation"}:
        raise ValueError("Alt grup analizi yalnızca validation tahminlerini kabul eder.")
    validation_ids = set(predictions["VariationID"].astype(str))
    annotations = load_validation_annotations(
        features_path, validation_ids, chunk_size=chunk_size
    )
    merged = predictions.merge(annotations, on="VariationID", how="left", validate="one_to_one")
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


def run_hgvs_ablation(
    features_path: Path,
    split_path: Path,
    report_dir: Path,
    artifact_dir: Path,
    *,
    seed: int = DEFAULT_SEED,
    threshold: float = 0.5,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, object]:
    report_dir.mkdir(parents=True, exist_ok=True)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    metrics_rows: list[dict[str, object]] = []
    profile_manifests: dict[str, object] = {}

    for profile_key, training_profile, label in PROFILES:
        profile_report_dir = report_dir / "models" / profile_key
        profile_artifact_dir = artifact_dir / profile_key
        manifest = train_baselines(
            features_path,
            split_path,
            profile_report_dir,
            profile_artifact_dir,
            seed=seed,
            threshold=threshold,
            chunk_size=chunk_size,
            feature_profile=training_profile,
        )
        if manifest["data_policy"]["test_partition"] != (
            "locked_not_loaded_into_development_dataframe"
        ):
            raise ValueError("Ablasyon profili test bölümünü kilitli tutmadı.")
        primary_metrics = manifest["validation_metrics"][PRIMARY_MODEL]
        metrics_rows.append(
            {
                "profile": profile_key,
                "label": label,
                "training_feature_profile": training_profile,
                **primary_metrics,
            }
        )
        profile_manifests[profile_key] = manifest

    metrics = pd.DataFrame(metrics_rows)
    metrics_path = report_dir / "ablation_metrics.csv"
    metrics.to_csv(metrics_path, index=False, float_format="%.10f")

    gain_metrics = ("roc_auc", "average_precision", "balanced_accuracy", "f1")
    gain_rows: list[dict[str, object]] = []
    for previous_index, current_index in ((0, 1), (1, 2), (0, 2)):
        previous = metrics.iloc[previous_index]
        current = metrics.iloc[current_index]
        for metric in gain_metrics:
            gain_rows.append(
                {
                    "from_profile": previous["profile"],
                    "to_profile": current["profile"],
                    "metric": metric,
                    "absolute_gain": float(current[metric] - previous[metric]),
                }
            )
    gains = pd.DataFrame(gain_rows)
    gains_path = report_dir / "incremental_gains.csv"
    gains.to_csv(gains_path, index=False, float_format="%.10f")

    full_profile_key = PROFILES[-1][0]
    predictions_path = (
        report_dir / "models" / full_profile_key / "validation_predictions.tsv.gz"
    )
    subgroup_table = build_subgroup_table(
        predictions_path,
        features_path,
        threshold=threshold,
        chunk_size=chunk_size,
    )
    subgroup_path = report_dir / "full_model_validation_subgroups.csv"
    subgroup_table.to_csv(subgroup_path, index=False, float_format="%.10f")

    figure, axis = plt.subplots(figsize=(10, 5.5))
    x_positions = range(len(metrics))
    axis.plot(x_positions, metrics["roc_auc"], marker="o", label="ROC-AUC")
    axis.plot(
        x_positions,
        metrics["average_precision"],
        marker="o",
        label="Average precision",
    )
    axis.set_xticks(list(x_positions), ["Yapısal", "+ HGVS biyolojik", "+ gösterim"])
    axis.set_ylabel("Validation skoru")
    axis.set_ylim(max(0.0, float(metrics[["roc_auc", "average_precision"]].min().min()) - 0.03), 1.0)
    axis.set_title("Gene-disjoint validation HGVS özellik ablasyonu")
    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    plot_path = report_dir / "ablation_roc_ap.png"
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)

    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "features": {"path": str(features_path.resolve()), "sha256": sha256_file(features_path)},
            "split": {"path": str(split_path.resolve()), "sha256": sha256_file(split_path)},
        },
        "data_policy": {
            "training_partition": "train",
            "evaluation_partition": "validation",
            "test_partition": "locked_not_loaded_into_development_dataframe",
            "primary_model": PRIMARY_MODEL,
            "threshold": threshold,
            "same_split_seed_and_estimator_across_profiles": True,
        },
        "profiles": [
            {
                "profile": key,
                "training_feature_profile": training_profile,
                "label": label,
                "feature_columns": profile_manifests[key]["models"][PRIMARY_MODEL][
                    "feature_columns"
                ],
            }
            for key, training_profile, label in PROFILES
        ],
        "interpretation_limit": (
            "The biological profile includes protein_effect; its "
            "no_protein_annotation category can implicitly encode annotation absence. "
            "The final step therefore measures the incremental value of explicit "
            "representation fields, not all annotation-availability signal."
        ),
        "reports": {
            "ablation_metrics": str(metrics_path.resolve()),
            "incremental_gains": str(gains_path.resolve()),
            "full_model_validation_subgroups": str(subgroup_path.resolve()),
            "plot": str(plot_path.resolve()),
        },
    }
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    report_path = report_dir / "hgvs_ablation_report.md"
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write("# HGVS özellik ablasyonu\n\n")
        handle.write(
            "Bütün modeller aynı gene-disjoint train/validation ayrımında ve aynı "
            "lojistik regresyon ayarlarıyla çalıştırıldı. Test bölümü yüklenmedi.\n\n"
        )
        handle.write("| Profil | ROC-AUC | AP | Duyarlılık | Özgüllük | F1 |\n")
        handle.write("| --- | ---: | ---: | ---: | ---: | ---: |\n")
        for row in metrics.to_dict("records"):
            handle.write(
                f"| {row['label']} | {row['roc_auc']:.4f} | "
                f"{row['average_precision']:.4f} | {row['recall_sensitivity']:.4f} | "
                f"{row['specificity']:.4f} | {row['f1']:.4f} |\n"
            )
        handle.write("\n## Tam model validation alt grupları\n\n")
        handle.write("| Grup | n | Negatif | Pozitif | Precision | Recall | ROC-AUC | AP | Özgüllük |\n")
        handle.write("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for row in subgroup_table.to_dict("records"):
            handle.write(
                f"| {row['group']} | {row['row_count']:,} | {row['negative_count']:,} | "
                f"{row['positive_count']:,} | {row['precision']:.4f} | "
                f"{row['recall_sensitivity']:.4f} | {row['roc_auc']:.4f} | "
                f"{row['average_precision']:.4f} | "
                f"{row['specificity']:.4f} |\n"
            )
        handle.write("\n## Yorum sınırı\n\n")
        handle.write(
            "Biyolojik profil `protein_effect` özelliğini içerir; "
            "`no_protein_annotation` kategorisi anotasyon yokluğunu dolaylı biçimde "
            "taşıyabilir. Bu nedenle son basamak tüm anotasyon-varlığı etkisini değil, "
            "açık gösterim alanlarının (`hgvs_sequence_type`, "
            "`has_protein_annotation`) ek katkısını ölçer.\n"
        )
    manifest["reports"]["report"] = str(report_path.resolve())
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("HGVS özellik ablasyonu tamamlandı")
    for row in metrics.to_dict("records"):
        print(
            f"{row['profile']}: ROC-AUC={row['roc_auc']:.4f}, "
            f"AP={row['average_precision']:.4f}, F1={row['f1']:.4f}"
        )
    print(f"Rapor: {report_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="HGVS özellik profillerini gene-disjoint validation üzerinde karşılaştırır."
    )
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--split", type=Path, default=DEFAULT_SPLIT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--artifact-dir", type=Path, default=DEFAULT_ARTIFACT_DIR)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    run_hgvs_ablation(
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
