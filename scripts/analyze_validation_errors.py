"""Ana baseline modelinin validation hatalarını alt gruplarda incele."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    precision_score,
    roc_auc_score,
)


DEFAULT_PREDICTIONS = Path(
    "reports/generated/logistic_baseline/validation_predictions.tsv.gz"
)
DEFAULT_FEATURES = Path("data/processed/clinvar_baseline_features.tsv.gz")
DEFAULT_REPORT_DIR = Path("reports/generated/logistic_baseline/error_analysis")
DEFAULT_MODEL_NAME = "without_missingness_flag"
DEFAULT_THRESHOLD = 0.5
DEFAULT_MIN_GROUP_SIZE = 100
DEFAULT_TOP_ERRORS = 50
DEFAULT_CHUNK_SIZE = 100_000

FEATURE_COLUMNS = (
    "variant_type",
    "chromosome",
    "mutation_class",
    "substitution",
    "is_canonical_snv",
    "vcf_alleles_missing",
    "reference_length",
    "alternate_length",
    "length_delta",
    "absolute_length_delta",
)
GROUP_FEATURES = (
    "variant_type",
    "chromosome",
    "mutation_class",
    "substitution",
    "is_canonical_snv",
    "vcf_alleles_missing",
    "absolute_length_delta_bucket",
)
PREDICTION_COLUMNS = {
    "VariationID",
    "model_name",
    "split",
    "target",
    "probability_pathogenic",
}


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def load_model_validation_predictions(
    predictions_path: Path,
    model_name: str,
) -> pd.DataFrame:
    if not predictions_path.is_file():
        raise FileNotFoundError(f"Tahmin dosyası bulunamadı: {predictions_path}")
    predictions = pd.read_csv(
        predictions_path,
        sep="\t",
        dtype={
            "VariationID": "string",
            "model_name": "string",
            "split": "string",
            "target": "int8",
            "probability_pathogenic": "float64",
        },
    )
    missing = sorted(PREDICTION_COLUMNS.difference(predictions.columns))
    if missing:
        raise ValueError("Eksik tahmin sütunları: " + ", ".join(missing))
    if set(predictions["split"].dropna().unique()) != {"validation"}:
        raise ValueError("Hata analizi yalnızca validation tahminlerini kabul eder.")

    selected = predictions.loc[predictions["model_name"] == model_name].copy()
    if selected.empty:
        available = ", ".join(sorted(predictions["model_name"].dropna().unique()))
        raise ValueError(f"Model bulunamadı: {model_name}. Mevcut modeller: {available}")
    if selected[list(PREDICTION_COLUMNS)].isna().any().any():
        raise ValueError("Seçilen modelin gerekli tahmin alanlarında boş değer var.")
    if selected["VariationID"].duplicated().any():
        raise ValueError("Seçilen modelde tekrarlanan VariationID bulundu.")
    if set(selected["target"].unique()) != {0, 1}:
        raise ValueError("Validation tahminlerinde iki hedef sınıfı da bulunmalıdır.")
    if not np.isfinite(selected["probability_pathogenic"]).all():
        raise ValueError("Sonlu olmayan olasılık değeri bulundu.")
    if not selected["probability_pathogenic"].between(0, 1).all():
        raise ValueError("Olasılıklar 0 ile 1 arasında olmalıdır.")
    return selected


def load_validation_features(
    features_path: Path,
    validation_ids: set[str],
    *,
    chunk_size: int,
) -> pd.DataFrame:
    """Özellik dosyasından yalnızca verilen validation kimliklerini belleğe al."""

    if not features_path.is_file():
        raise FileNotFoundError(f"Özellik dosyası bulunamadı: {features_path}")
    if chunk_size <= 0:
        raise ValueError("chunk_size pozitif olmalıdır.")

    retained: list[pd.DataFrame] = []
    dtypes = {
        "VariationID": "string",
        "variant_type": "string",
        "chromosome": "string",
        "mutation_class": "string",
        "substitution": "string",
        "is_canonical_snv": "Int8",
        "vcf_alleles_missing": "Int8",
    }
    for chunk in pd.read_csv(
        features_path,
        sep="\t",
        usecols=["VariationID", *FEATURE_COLUMNS],
        dtype=dtypes,
        chunksize=chunk_size,
    ):
        mask = chunk["VariationID"].isin(validation_ids)
        if mask.any():
            retained.append(chunk.loc[mask].copy())

    if not retained:
        raise ValueError("Validation kimliklerinin hiçbiri özellik dosyasında bulunamadı.")
    features = pd.concat(retained, ignore_index=True)
    if features["VariationID"].duplicated().any():
        duplicate = features.loc[
            features["VariationID"].duplicated(), "VariationID"
        ].iloc[0]
        raise ValueError(f"Özellik dosyasında tekrarlanan VariationID: {duplicate}")
    found_ids = set(features["VariationID"].tolist())
    missing_ids = validation_ids.difference(found_ids)
    if missing_ids:
        raise ValueError(
            f"Özellik dosyasında bulunmayan validation VariationID: {min(missing_ids)}"
        )
    return features


def add_analysis_columns(data: pd.DataFrame, *, threshold: float) -> pd.DataFrame:
    result = data.copy()
    result["prediction"] = (
        result["probability_pathogenic"] >= threshold
    ).astype("int8")
    result["error_type"] = "correct"
    result.loc[
        (result["target"] == 1) & (result["prediction"] == 0), "error_type"
    ] = "false_negative"
    result.loc[
        (result["target"] == 0) & (result["prediction"] == 1), "error_type"
    ] = "false_positive"
    length_bucket = pd.cut(
        result["absolute_length_delta"],
        bins=[-np.inf, 0, 1, 5, 10, 50, 100, np.inf],
        labels=["0", "1", "2-5", "6-10", "11-50", "51-100", "101+"],
    ).astype("string")
    result["absolute_length_delta_bucket"] = length_bucket.fillna("missing")
    return result


def calculate_group_metrics(
    group: pd.DataFrame,
    *,
    min_group_size: int,
) -> dict[str, object]:
    targets = group["target"].to_numpy()
    predictions = group["prediction"].to_numpy()
    probabilities = group["probability_pathogenic"].to_numpy()
    true_negatives, false_positives, false_negatives, true_positives = confusion_matrix(
        targets, predictions, labels=[0, 1]
    ).ravel()
    positive_count = int((targets == 1).sum())
    negative_count = int((targets == 0).sum())
    both_classes = positive_count > 0 and negative_count > 0
    return {
        "row_count": len(group),
        "below_minimum_group_size": len(group) < min_group_size,
        "negative_count": negative_count,
        "positive_count": positive_count,
        "pathogenic_prevalence": float(targets.mean()),
        "mean_predicted_probability": float(probabilities.mean()),
        "roc_auc": roc_auc_score(targets, probabilities) if both_classes else np.nan,
        "brier_score": brier_score_loss(targets, probabilities),
        "accuracy": accuracy_score(targets, predictions),
        "balanced_accuracy": (
            balanced_accuracy_score(targets, predictions) if both_classes else np.nan
        ),
        "precision_positive_predictive_value": precision_score(
            targets, predictions, zero_division=0
        ),
        "sensitivity": true_positives / positive_count if positive_count else np.nan,
        "specificity": true_negatives / negative_count if negative_count else np.nan,
        "false_negative_rate": false_negatives / positive_count if positive_count else np.nan,
        "false_positive_rate": false_positives / negative_count if negative_count else np.nan,
        "f1": f1_score(targets, predictions, zero_division=0),
        "true_negatives": int(true_negatives),
        "false_positives": int(false_positives),
        "false_negatives": int(false_negatives),
        "true_positives": int(true_positives),
    }


def build_subgroup_metrics(
    data: pd.DataFrame,
    *,
    min_group_size: int,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = [
        {
            "group_feature": "overall",
            "group_value": "all_validation_rows",
            **calculate_group_metrics(data, min_group_size=min_group_size),
        }
    ]
    for feature in GROUP_FEATURES:
        grouped = data.groupby(feature, dropna=False, observed=True, sort=False)
        for group_value, group in grouped:
            value = "missing" if pd.isna(group_value) else str(group_value)
            rows.append(
                {
                    "group_feature": feature,
                    "group_value": value,
                    **calculate_group_metrics(group, min_group_size=min_group_size),
                }
            )
    result = pd.DataFrame(rows)
    return result.sort_values(
        ["group_feature", "row_count", "group_value"],
        ascending=[True, False, True],
        ignore_index=True,
    )


def build_confident_errors(data: pd.DataFrame, *, top_errors: int) -> pd.DataFrame:
    false_negatives = data.loc[data["error_type"] == "false_negative"].copy()
    false_negatives["wrong_class_confidence"] = (
        1 - false_negatives["probability_pathogenic"]
    )
    false_negatives = false_negatives.sort_values(
        ["wrong_class_confidence", "VariationID"], ascending=[False, True]
    ).head(top_errors)

    false_positives = data.loc[data["error_type"] == "false_positive"].copy()
    false_positives["wrong_class_confidence"] = false_positives[
        "probability_pathogenic"
    ]
    false_positives = false_positives.sort_values(
        ["wrong_class_confidence", "VariationID"], ascending=[False, True]
    ).head(top_errors)

    errors = pd.concat([false_negatives, false_positives], ignore_index=True)
    errors["rank_within_error_type"] = errors.groupby("error_type").cumcount() + 1
    columns = [
        "error_type",
        "rank_within_error_type",
        "VariationID",
        "target",
        "prediction",
        "probability_pathogenic",
        "wrong_class_confidence",
        *FEATURE_COLUMNS,
        "absolute_length_delta_bucket",
    ]
    return errors.loc[:, columns]


def plot_error_rates(
    subgroup_metrics: pd.DataFrame,
    group_feature: str,
    output_path: Path,
    *,
    min_group_size: int,
) -> None:
    rows = subgroup_metrics.loc[
        (subgroup_metrics["group_feature"] == group_feature)
        & (subgroup_metrics["row_count"] >= min_group_size)
    ].copy()
    rows = rows.sort_values("false_negative_rate", ascending=True)
    if rows.empty:
        return

    positions = np.arange(len(rows))
    bar_height = 0.38
    figure_height = max(4.5, 0.45 * len(rows) + 1.5)
    figure, axis = plt.subplots(figsize=(10, figure_height))
    axis.barh(
        positions - bar_height / 2,
        rows["false_negative_rate"],
        height=bar_height,
        label="Yanlış negatif oranı",
    )
    axis.barh(
        positions + bar_height / 2,
        rows["false_positive_rate"],
        height=bar_height,
        label="Yanlış pozitif oranı",
    )
    labels = [
        f"{value} (n={count:,})"
        for value, count in zip(rows["group_value"], rows["row_count"])
    ]
    axis.set_yticks(positions, labels)
    axis.set_xlim(0, 1)
    axis.set_xlabel("Hata oranı")
    axis.set_title(f"Validation hata oranları: {group_feature}")
    axis.grid(axis="x", alpha=0.25)
    axis.legend()
    figure.tight_layout()
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)


def analyze_errors(
    predictions_path: Path,
    features_path: Path,
    report_dir: Path,
    *,
    model_name: str = DEFAULT_MODEL_NAME,
    threshold: float = DEFAULT_THRESHOLD,
    min_group_size: int = DEFAULT_MIN_GROUP_SIZE,
    top_errors: int = DEFAULT_TOP_ERRORS,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, object]:
    if not 0 < threshold < 1:
        raise ValueError("Eşik 0 ile 1 arasında olmalıdır.")
    if min_group_size <= 0 or top_errors <= 0:
        raise ValueError("Grup boyutu ve hata örneği sayısı pozitif olmalıdır.")

    predictions = load_model_validation_predictions(predictions_path, model_name)
    validation_ids = set(predictions["VariationID"].tolist())
    features = load_validation_features(
        features_path,
        validation_ids,
        chunk_size=chunk_size,
    )
    data = predictions.merge(
        features,
        on="VariationID",
        how="inner",
        validate="one_to_one",
    )
    data = add_analysis_columns(data, threshold=threshold)
    subgroup_metrics = build_subgroup_metrics(
        data,
        min_group_size=min_group_size,
    )
    confident_errors = build_confident_errors(data, top_errors=top_errors)

    report_dir.mkdir(parents=True, exist_ok=True)
    subgroup_path = report_dir / "subgroup_metrics.csv"
    errors_path = report_dir / "most_confident_errors.tsv"
    subgroup_metrics.to_csv(subgroup_path, index=False, float_format="%.10f")
    confident_errors.to_csv(errors_path, sep="\t", index=False, float_format="%.10f")

    plot_paths: dict[str, str] = {}
    for group_feature in ("variant_type", "mutation_class", "chromosome"):
        plot_path = report_dir / f"error_rates_by_{group_feature}.png"
        plot_error_rates(
            subgroup_metrics,
            group_feature,
            plot_path,
            min_group_size=min_group_size,
        )
        if plot_path.is_file():
            plot_paths[group_feature] = str(plot_path.resolve())

    overall = subgroup_metrics.loc[
        subgroup_metrics["group_feature"] == "overall"
    ].iloc[0]
    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "predictions": {
                "path": str(predictions_path.resolve()),
                "sha256": sha256_file(predictions_path),
            },
            "features": {
                "path": str(features_path.resolve()),
                "sha256": sha256_file(features_path),
            },
        },
        "data_policy": {
            "analysis_partition": "validation",
            "test_data_used": False,
            "model_name": model_name,
            "threshold": threshold,
            "minimum_group_size_for_plots": min_group_size,
            "small_group_warning": (
                "Tabloda tüm gruplar bulunur; minimumun altındaki oranlar kararsız "
                "olabileceği için ayrıca işaretlenir ve grafiklere alınmaz."
            ),
        },
        "counts": {
            "validation_rows": len(data),
            "false_negatives": int(overall["false_negatives"]),
            "false_positives": int(overall["false_positives"]),
            "confident_error_examples": len(confident_errors),
        },
        "group_features": list(GROUP_FEATURES),
        "outputs": {
            "subgroup_metrics": str(subgroup_path.resolve()),
            "most_confident_errors": str(errors_path.resolve()),
            "plots": plot_paths,
        },
    }
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    eligible_variant_types = subgroup_metrics.loc[
        (subgroup_metrics["group_feature"] == "variant_type")
        & (~subgroup_metrics["below_minimum_group_size"])
        & subgroup_metrics["false_negative_rate"].notna()
    ]
    worst_false_negative = eligible_variant_types.sort_values(
        "false_negative_rate", ascending=False
    ).iloc[0]
    print("Validation alt grup hata analizi tamamlandı")
    print(f"Model: {model_name}; eşik: {threshold:.4f}")
    print(f"Validation satırı: {len(data):,}")
    print(
        "En yüksek yanlış negatif oranlı yeterli destekli varyant türü: "
        f"{worst_false_negative['group_value']} "
        f"(%{100 * worst_false_negative['false_negative_rate']:.2f}, "
        f"n={int(worst_false_negative['row_count']):,})"
    )
    print("Test verisi kullanılmadı.")
    print(f"Rapor: {manifest_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Ana modelin validation hatalarını alt gruplarda inceler."
    )
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--model-name", default=DEFAULT_MODEL_NAME)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--min-group-size", type=int, default=DEFAULT_MIN_GROUP_SIZE)
    parser.add_argument("--top-errors", type=int, default=DEFAULT_TOP_ERRORS)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    analyze_errors(
        args.predictions,
        args.features,
        args.report_dir,
        model_name=args.model_name,
        threshold=args.threshold,
        min_group_size=args.min_group_size,
        top_errors=args.top_errors,
        chunk_size=args.chunk_size,
    )


if __name__ == "__main__":
    main()
