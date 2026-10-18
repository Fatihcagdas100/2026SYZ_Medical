"""Yalnızca train bölümünde gen-grubu çakışmasız HGVS çapraz doğrulaması yap."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

try:
    from scripts.train_logistic_baseline import (
        FEATURE_PROFILES,
        build_pipeline,
        calculate_metrics,
    )
except ModuleNotFoundError:  # Doğrudan `python scripts/...` çalıştırması için.
    from train_logistic_baseline import (  # type: ignore[no-redef]
        FEATURE_PROFILES,
        build_pipeline,
        calculate_metrics,
    )


DEFAULT_FEATURES = Path("data/processed/clinvar_hgvs_features.tsv.gz")
DEFAULT_SPLIT = Path("data/processed/baseline_gene_disjoint_split.tsv")
DEFAULT_REPORT_DIR = Path("reports/generated/hgvs_gene_disjoint_cv")
DEFAULT_SEED = 2027
DEFAULT_N_SPLITS = 5
DEFAULT_THRESHOLD = 0.5
DEFAULT_CHUNK_SIZE = 100_000
FEATURE_PROFILE = "hgvs"
SUMMARY_METRICS = (
    "roc_auc",
    "average_precision",
    "log_loss",
    "brier_score",
    "accuracy",
    "balanced_accuracy",
    "precision",
    "recall_sensitivity",
    "specificity",
    "f1",
    "matthews_correlation_coefficient",
)


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def load_train_group_assignments(
    split_path: Path,
) -> tuple[pd.DataFrame, dict[str, int]]:
    if not split_path.is_file():
        raise FileNotFoundError(f"Ayrım dosyası bulunamadı: {split_path}")
    split_table = pd.read_csv(
        split_path,
        sep="\t",
        usecols=["VariationID", "split", "gene_group"],
        dtype={"VariationID": "string", "split": "string", "gene_group": "string"},
    )
    if split_table.isna().any().any():
        raise ValueError("Ayrım dosyasının gerekli alanlarında boş değer bulundu.")
    if split_table["VariationID"].duplicated().any():
        raise ValueError("Ayrım dosyasında tekrarlanan VariationID bulundu.")
    observed_splits = set(split_table["split"].unique())
    if observed_splits != {"train", "validation", "test"}:
        raise ValueError("Ayrım dosyası train, validation ve test içermelidir.")
    partition_counts = {
        split_name: int((split_table["split"] == split_name).sum())
        for split_name in ("train", "validation", "test")
    }
    train_groups = split_table.loc[
        split_table["split"] == "train", ["VariationID", "gene_group"]
    ].copy()
    return train_groups, partition_counts


def load_train_features(
    features_path: Path,
    train_groups: pd.DataFrame,
    *,
    feature_names: tuple[str, ...],
    categorical_features: tuple[str, ...],
    chunk_size: int,
) -> pd.DataFrame:
    if not features_path.is_file():
        raise FileNotFoundError(f"Özellik dosyası bulunamadı: {features_path}")
    if chunk_size <= 0:
        raise ValueError("chunk_size pozitif olmalıdır.")

    train_ids = set(train_groups["VariationID"].tolist())
    dtypes: dict[str, str] = {
        "VariationID": "string",
        "target": "Int8",
        "is_canonical_snv": "Int8",
        "has_protein_annotation": "Int8",
    }
    dtypes.update({feature: "string" for feature in categorical_features})
    retained: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        features_path,
        sep="\t",
        usecols=["VariationID", *feature_names, "target"],
        dtype=dtypes,
        chunksize=chunk_size,
    ):
        mask = chunk["VariationID"].isin(train_ids)
        if mask.any():
            retained.append(chunk.loc[mask].copy())
    if not retained:
        raise ValueError("Train kimliklerinin hiçbiri özellik dosyasında bulunamadı.")

    train_data = pd.concat(retained, ignore_index=True)
    if train_data["VariationID"].duplicated().any():
        raise ValueError("Train özelliklerinde tekrarlanan VariationID bulundu.")
    found_ids = set(train_data["VariationID"].tolist())
    if found_ids != train_ids:
        missing = train_ids.difference(found_ids)
        extra = found_ids.difference(train_ids)
        detail = min(missing or extra)
        raise ValueError(f"Train kimliği eşleşmiyor: {detail}")
    if train_data["target"].isna().any() or set(
        train_data["target"].astype(int).unique()
    ) != {0, 1}:
        raise ValueError("Train özelliklerinde iki hedef sınıfı da bulunmalıdır.")

    train_data = train_data.merge(
        train_groups,
        on="VariationID",
        how="inner",
        validate="one_to_one",
    )
    return train_data


def plot_fold_stability(fold_metrics: pd.DataFrame, output_path: Path) -> None:
    metrics = ["roc_auc", "f1", "balanced_accuracy", "recall_sensitivity", "specificity"]
    positions = np.arange(len(fold_metrics))
    width = 0.15
    figure, axis = plt.subplots(figsize=(11, 5.5))
    offsets = np.linspace(-2 * width, 2 * width, len(metrics))
    for offset, metric in zip(offsets, metrics):
        axis.bar(positions + offset, fold_metrics[metric], width=width, label=metric)
    axis.set_xticks(positions, [f"Fold {fold}" for fold in fold_metrics["fold"]])
    axis.set_ylim(0, 1)
    axis.set_ylabel("Metrik değeri")
    axis.set_title("Train-only gene-disjoint çapraz doğrulama kararlılığı")
    axis.grid(axis="y", alpha=0.25)
    axis.legend(fontsize=8, ncol=3)
    figure.tight_layout()
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)


def cross_validate(
    features_path: Path,
    split_path: Path,
    report_dir: Path,
    *,
    seed: int = DEFAULT_SEED,
    n_splits: int = DEFAULT_N_SPLITS,
    threshold: float = DEFAULT_THRESHOLD,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    feature_profile: str = FEATURE_PROFILE,
) -> dict[str, object]:
    if n_splits < 2:
        raise ValueError("n_splits en az 2 olmalıdır.")
    if not 0 < threshold < 1:
        raise ValueError("Eşik 0 ile 1 arasında olmalıdır.")
    if feature_profile not in FEATURE_PROFILES:
        raise ValueError(
            f"Bilinmeyen özellik profili: {feature_profile}. "
            f"Seçenekler: {', '.join(FEATURE_PROFILES)}"
        )

    profile = FEATURE_PROFILES[feature_profile]
    categorical_features = tuple(profile["categorical"])
    numeric_features = tuple(profile["numeric"])
    feature_names = (*categorical_features, *numeric_features)
    train_groups, partition_counts = load_train_group_assignments(split_path)
    train_data = load_train_features(
        features_path,
        train_groups,
        feature_names=feature_names,
        categorical_features=categorical_features,
        chunk_size=chunk_size,
    )
    if train_data["gene_group"].nunique() < n_splits:
        raise ValueError("Gen grubu sayısı fold sayısından az olamaz.")

    targets = train_data["target"].astype(int).reset_index(drop=True)
    groups = train_data["gene_group"].reset_index(drop=True)
    features = train_data.loc[:, feature_names].reset_index(drop=True)
    train_data = train_data.reset_index(drop=True)
    splitter = StratifiedGroupKFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=seed,
    )
    oof_probabilities = np.full(len(train_data), np.nan, dtype=float)
    oof_folds = np.full(len(train_data), -1, dtype=np.int16)
    fold_metric_rows: list[dict[str, object]] = []
    coefficient_tables: list[pd.DataFrame] = []

    for fold, (fit_indices, holdout_indices) in enumerate(
        splitter.split(features, targets, groups=groups),
        start=1,
    ):
        fit_groups = set(groups.iloc[fit_indices])
        holdout_groups = set(groups.iloc[holdout_indices])
        if fit_groups.intersection(holdout_groups):
            raise RuntimeError(f"Fold {fold} içinde gen grubu çakışması bulundu.")

        pipeline = build_pipeline(
            feature_names,
            seed=seed,
            categorical_features=categorical_features,
        )
        pipeline.fit(features.iloc[fit_indices], targets.iloc[fit_indices])
        probabilities = pipeline.predict_proba(features.iloc[holdout_indices])[:, 1]
        oof_probabilities[holdout_indices] = probabilities
        oof_folds[holdout_indices] = fold
        metrics = calculate_metrics(
            targets.iloc[holdout_indices],
            pd.Series(probabilities, index=holdout_indices),
            threshold=threshold,
        )
        classifier = pipeline.named_steps["classifier"]
        fold_metric_rows.append(
            {
                "fold": fold,
                "fit_rows": len(fit_indices),
                "holdout_rows": len(holdout_indices),
                "fit_gene_groups": len(fit_groups),
                "holdout_gene_groups": len(holdout_groups),
                "iterations": int(classifier.n_iter_[0]),
                "converged": bool(classifier.n_iter_[0] < classifier.max_iter),
                **metrics,
            }
        )

        transformed_features = pipeline.named_steps[
            "preprocessing"
        ].get_feature_names_out()
        coefficients = classifier.coef_[0]
        coefficient_tables.append(
            pd.DataFrame(
                {
                    "fold": fold,
                    "transformed_feature": transformed_features,
                    "coefficient": coefficients,
                    "absolute_coefficient": abs(coefficients),
                }
            )
        )

    if np.isnan(oof_probabilities).any() or (oof_folds < 0).any():
        raise RuntimeError("Bazı train satırları out-of-fold tahmin almadı.")
    if len(set(train_data["VariationID"])) != len(train_data):
        raise RuntimeError("Out-of-fold tahminlerde VariationID tekilliği bozuldu.")

    fold_metrics = pd.DataFrame(fold_metric_rows)
    metric_summary_rows = []
    for metric in SUMMARY_METRICS:
        values = fold_metrics[metric]
        metric_summary_rows.append(
            {
                "metric": metric,
                "mean": values.mean(),
                "sample_std": values.std(ddof=1),
                "minimum": values.min(),
                "maximum": values.max(),
            }
        )
    metric_summary = pd.DataFrame(metric_summary_rows)
    oof_metrics = calculate_metrics(
        targets,
        pd.Series(oof_probabilities),
        threshold=threshold,
    )
    oof_predictions = pd.DataFrame(
        {
            "VariationID": train_data["VariationID"],
            "gene_group": groups,
            "source_partition": "train",
            "fold": oof_folds,
            "target": targets,
            "probability_pathogenic": oof_probabilities,
            "prediction": (oof_probabilities >= threshold).astype(int),
        }
    )
    if oof_predictions.groupby("gene_group")["fold"].nunique().max() != 1:
        raise RuntimeError("Bir gen grubu birden çok holdout fold'una atandı.")

    fold_coefficients = pd.concat(coefficient_tables, ignore_index=True)
    coefficient_stability = (
        fold_coefficients.groupby("transformed_feature")["coefficient"]
        .agg(
            folds_present="count",
            mean="mean",
            sample_std="std",
            minimum="min",
            maximum="max",
        )
        .reset_index()
    )
    sign_stability = (
        fold_coefficients.assign(positive=fold_coefficients["coefficient"] > 0)
        .groupby("transformed_feature")["positive"]
        .mean()
        .rename("positive_fraction")
        .reset_index()
    )
    coefficient_stability = coefficient_stability.merge(
        sign_stability,
        on="transformed_feature",
        validate="one_to_one",
    )
    coefficient_stability["absolute_mean"] = coefficient_stability["mean"].abs()
    coefficient_stability = coefficient_stability.sort_values(
        "absolute_mean", ascending=False
    )

    report_dir.mkdir(parents=True, exist_ok=True)
    fold_metrics_path = report_dir / "fold_metrics.csv"
    summary_path = report_dir / "metric_summary.csv"
    predictions_path = report_dir / "train_out_of_fold_predictions.tsv.gz"
    fold_coefficients_path = report_dir / "fold_coefficients.csv"
    coefficient_stability_path = report_dir / "coefficient_stability.csv"
    plot_path = report_dir / "fold_metric_stability.png"
    fold_metrics.to_csv(fold_metrics_path, index=False, float_format="%.10f")
    metric_summary.to_csv(summary_path, index=False, float_format="%.10f")
    oof_predictions.to_csv(
        predictions_path,
        sep="\t",
        index=False,
        compression="gzip",
        float_format="%.10f",
    )
    fold_coefficients.to_csv(
        fold_coefficients_path, index=False, float_format="%.10f"
    )
    coefficient_stability.to_csv(
        coefficient_stability_path, index=False, float_format="%.10f"
    )
    plot_fold_stability(fold_metrics, plot_path)

    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "features": {
                "path": str(features_path.resolve()),
                "sha256": sha256_file(features_path),
            },
            "gene_disjoint_split": {
                "path": str(split_path.resolve()),
                "sha256": sha256_file(split_path),
            },
        },
        "data_policy": {
            "source_partition": "train",
            "fixed_validation_used": False,
            "test_used": False,
            "feature_profile": feature_profile,
            "feature_columns": list(feature_names),
            "vcf_missingness_flag_used": False,
            "threshold": threshold,
            "threshold_policy": "fixed_in_advance_not_tuned",
        },
        "cross_validation": {
            "method": "StratifiedGroupKFold",
            "n_splits": n_splits,
            "shuffle": True,
            "seed": seed,
            "group_column": "gene_group",
        },
        "counts": {
            "train_rows": len(train_data),
            "train_gene_groups": int(groups.nunique()),
            "excluded_fixed_validation_rows": partition_counts["validation"],
            "excluded_test_rows": partition_counts["test"],
            "out_of_fold_prediction_rows": len(oof_predictions),
            "gene_groups_in_multiple_holdout_folds": 0,
        },
        "out_of_fold_metrics": oof_metrics,
        "outputs": {
            "fold_metrics": str(fold_metrics_path.resolve()),
            "metric_summary": str(summary_path.resolve()),
            "out_of_fold_predictions": str(predictions_path.resolve()),
            "fold_coefficients": str(fold_coefficients_path.resolve()),
            "coefficient_stability": str(coefficient_stability_path.resolve()),
            "stability_plot": str(plot_path.resolve()),
        },
    }
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    roc_summary = metric_summary.loc[metric_summary["metric"] == "roc_auc"].iloc[0]
    print("Train-only gene-disjoint çapraz doğrulama tamamlandı")
    print(f"Özellik profili: {feature_profile}")
    print(f"Train satırı: {len(train_data):,}; gen grubu: {groups.nunique():,}")
    print(
        f"Fold ROC-AUC: {roc_summary['mean']:.4f} ± "
        f"{roc_summary['sample_std']:.4f} "
        f"(min={roc_summary['minimum']:.4f}, max={roc_summary['maximum']:.4f})"
    )
    print(f"Birleşik OOF ROC-AUC: {oof_metrics['roc_auc']:.4f}")
    print("Sabit validation ve test verisi kullanılmadı.")
    print(f"Rapor: {manifest_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Yalnızca gene-disjoint train bölümünde HGVS modeli için "
            "grup-korumalı çapraz doğrulama yapar."
        )
    )
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--split", type=Path, default=DEFAULT_SPLIT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--n-splits", type=int, default=DEFAULT_N_SPLITS)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    parser.add_argument(
        "--feature-profile",
        choices=tuple(FEATURE_PROFILES),
        default=FEATURE_PROFILE,
        help="Eğitimde kullanılacak özellik profili.",
    )
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    cross_validate(
        args.features,
        args.split,
        args.report_dir,
        seed=args.seed,
        n_splits=args.n_splits,
        threshold=args.threshold,
        chunk_size=args.chunk_size,
        feature_profile=args.feature_profile,
    )


if __name__ == "__main__":
    main()
