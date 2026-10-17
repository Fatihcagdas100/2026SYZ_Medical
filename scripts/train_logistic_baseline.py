"""ClinVar temel özellikleriyle sızıntıya kapalı lojistik regresyon eğit."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path
from typing import Sequence

import joblib
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    log_loss,
    matthews_corrcoef,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


DEFAULT_FEATURES = Path("data/processed/clinvar_baseline_features.tsv.gz")
DEFAULT_SPLIT = Path("data/processed/baseline_stratified_split.tsv")
DEFAULT_REPORT_DIR = Path("reports/generated/logistic_baseline")
DEFAULT_ARTIFACT_DIR = Path("artifacts/logistic_baseline")
DEFAULT_SEED = 2027
DEFAULT_CHUNK_SIZE = 100_000

CATEGORICAL_FEATURES = (
    "variant_type",
    "chromosome",
    "substitution",
    "mutation_class",
)
NUMERIC_FEATURES = (
    "reference_length",
    "alternate_length",
    "length_delta",
    "absolute_length_delta",
    "is_canonical_snv",
)
MISSINGNESS_FEATURE = "vcf_alleles_missing"
FEATURE_COLUMNS = (*CATEGORICAL_FEATURES, *NUMERIC_FEATURES, MISSINGNESS_FEATURE)
REQUIRED_FEATURE_COLUMNS = ("VariationID", *FEATURE_COLUMNS, "target")
MODEL_FEATURE_SETS = {
    "without_missingness_flag": (*CATEGORICAL_FEATURES, *NUMERIC_FEATURES),
    "with_missingness_flag": (*CATEGORICAL_FEATURES, *NUMERIC_FEATURES, MISSINGNESS_FEATURE),
}
EXPECTED_SPLITS = {"train", "validation", "test"}


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def load_split_assignments(split_path: Path) -> pd.Series:
    if not split_path.is_file():
        raise FileNotFoundError(f"Ayrım dosyası bulunamadı: {split_path}")

    split_table = pd.read_csv(
        split_path,
        sep="\t",
        usecols=["VariationID", "split"],
        dtype={"VariationID": "string", "split": "string"},
    )
    if split_table["VariationID"].isna().any():
        raise ValueError("Ayrım dosyasında boş VariationID bulundu.")
    if split_table["VariationID"].duplicated().any():
        duplicate = split_table.loc[
            split_table["VariationID"].duplicated(), "VariationID"
        ].iloc[0]
        raise ValueError(f"Ayrım dosyasında tekrarlanan VariationID bulundu: {duplicate}")

    observed_splits = set(split_table["split"].dropna().tolist())
    if observed_splits != EXPECTED_SPLITS:
        raise ValueError(
            "Ayrım değerleri train, validation ve test olmalıdır; bulunan: "
            + ", ".join(sorted(observed_splits))
        )
    return split_table.set_index("VariationID")["split"]


def load_development_data(
    features_path: Path,
    split_path: Path,
    *,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Yalnızca train/validation satırlarını döndür; test satırlarını kilitli tut."""

    if not features_path.is_file():
        raise FileNotFoundError(f"Özellik dosyası bulunamadı: {features_path}")
    if chunk_size <= 0:
        raise ValueError("chunk_size pozitif olmalıdır.")

    assignments = load_split_assignments(split_path)
    retained_chunks: list[pd.DataFrame] = []
    seen_ids: set[str] = set()
    partition_counts: Counter[str] = Counter()
    dtypes = {
        "VariationID": "string",
        "variant_type": "string",
        "chromosome": "string",
        "substitution": "string",
        "mutation_class": "string",
        "vcf_alleles_missing": "Int8",
        "is_canonical_snv": "Int8",
        "target": "Int8",
    }

    for chunk in pd.read_csv(
        features_path,
        sep="\t",
        usecols=list(REQUIRED_FEATURE_COLUMNS),
        dtype=dtypes,
        chunksize=chunk_size,
    ):
        if chunk["VariationID"].isna().any():
            raise ValueError("Özellik dosyasında boş VariationID bulundu.")
        chunk_ids = chunk["VariationID"].tolist()
        duplicate_in_chunk = chunk["VariationID"].duplicated()
        if duplicate_in_chunk.any():
            duplicate = chunk.loc[duplicate_in_chunk, "VariationID"].iloc[0]
            raise ValueError(f"Özellik dosyasında tekrarlanan VariationID: {duplicate}")
        repeated_between_chunks = seen_ids.intersection(chunk_ids)
        if repeated_between_chunks:
            duplicate = min(repeated_between_chunks)
            raise ValueError(f"Özellik dosyasında tekrarlanan VariationID: {duplicate}")
        seen_ids.update(chunk_ids)

        chunk_splits = chunk["VariationID"].map(assignments)
        if chunk_splits.isna().any():
            missing_id = chunk.loc[chunk_splits.isna(), "VariationID"].iloc[0]
            raise ValueError(f"Ayrım ataması olmayan VariationID: {missing_id}")
        partition_counts.update(chunk_splits.tolist())

        development_mask = chunk_splits.isin(["train", "validation"])
        development_chunk = chunk.loc[development_mask].copy()
        development_chunk["split"] = chunk_splits.loc[development_mask].to_numpy()
        retained_chunks.append(development_chunk)

    missing_feature_ids = set(assignments.index).difference(seen_ids)
    if missing_feature_ids:
        missing_id = min(missing_feature_ids)
        raise ValueError(f"Özellik dosyasında bulunmayan ayrım kimliği: {missing_id}")

    development = pd.concat(retained_chunks, ignore_index=True)
    if set(development["split"].unique()) != {"train", "validation"}:
        raise ValueError("Geliştirme verisinde train ve validation birlikte bulunmalıdır.")
    if development["target"].isna().any() or not set(
        development["target"].astype(int).unique()
    ).issubset({0, 1}):
        raise ValueError("Train/validation hedefleri yalnızca 0 veya 1 olmalıdır.")
    for split_name in ("train", "validation"):
        targets = set(
            development.loc[development["split"] == split_name, "target"]
            .astype(int)
            .unique()
        )
        if targets != {0, 1}:
            raise ValueError(f"{split_name} bölümünde iki hedef sınıfı da bulunmalıdır.")

    return development, {
        split_name: partition_counts[split_name]
        for split_name in ("train", "validation", "test")
    }


def build_pipeline(feature_names: Sequence[str], *, seed: int) -> Pipeline:
    categorical = [name for name in CATEGORICAL_FEATURES if name in feature_names]
    numeric = [name for name in feature_names if name not in categorical]
    preprocessing = ColumnTransformer(
        transformers=[
            (
                "categorical",
                Pipeline(
                    steps=[
                        ("imputer", SimpleImputer(strategy="most_frequent")),
                        ("onehot", OneHotEncoder(handle_unknown="ignore")),
                    ]
                ),
                categorical,
            ),
            (
                "numeric",
                Pipeline(
                    steps=[
                        ("imputer", SimpleImputer(strategy="median")),
                        ("scaler", StandardScaler()),
                    ]
                ),
                numeric,
            ),
        ],
        sparse_threshold=1.0,
        verbose_feature_names_out=False,
    )
    classifier = LogisticRegression(
        C=1.0,
        solver="lbfgs",
        max_iter=1_000,
        random_state=seed,
    )
    return Pipeline(
        steps=[("preprocessing", preprocessing), ("classifier", classifier)]
    )


def calculate_metrics(
    targets: pd.Series,
    probabilities: pd.Series,
    *,
    threshold: float,
) -> dict[str, float | int]:
    predictions = (probabilities >= threshold).astype(int)
    true_negatives, false_positives, false_negatives, true_positives = confusion_matrix(
        targets, predictions, labels=[0, 1]
    ).ravel()
    specificity = true_negatives / (true_negatives + false_positives)
    return {
        "threshold": threshold,
        "row_count": len(targets),
        "negative_count": int((targets == 0).sum()),
        "positive_count": int((targets == 1).sum()),
        "roc_auc": roc_auc_score(targets, probabilities),
        "average_precision": average_precision_score(targets, probabilities),
        "log_loss": log_loss(targets, probabilities, labels=[0, 1]),
        "brier_score": brier_score_loss(targets, probabilities),
        "accuracy": accuracy_score(targets, predictions),
        "balanced_accuracy": balanced_accuracy_score(targets, predictions),
        "precision": precision_score(targets, predictions, zero_division=0),
        "recall_sensitivity": recall_score(targets, predictions, zero_division=0),
        "specificity": specificity,
        "f1": f1_score(targets, predictions, zero_division=0),
        "matthews_correlation_coefficient": matthews_corrcoef(targets, predictions),
        "true_negatives": int(true_negatives),
        "false_positives": int(false_positives),
        "false_negatives": int(false_negatives),
        "true_positives": int(true_positives),
    }


def train_baselines(
    features_path: Path,
    split_path: Path,
    report_dir: Path,
    artifact_dir: Path,
    *,
    seed: int = DEFAULT_SEED,
    threshold: float = 0.5,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, object]:
    if not 0 < threshold < 1:
        raise ValueError("Eşik 0 ile 1 arasında olmalıdır.")

    development, partition_counts = load_development_data(
        features_path,
        split_path,
        chunk_size=chunk_size,
    )
    train = development.loc[development["split"] == "train"].copy()
    validation = development.loc[development["split"] == "validation"].copy()
    y_train = train["target"].astype(int)
    y_validation = validation["target"].astype(int)

    report_dir.mkdir(parents=True, exist_ok=True)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    all_metrics: list[dict[str, object]] = []
    prediction_tables: list[pd.DataFrame] = []
    coefficient_tables: list[pd.DataFrame] = []
    model_details: dict[str, object] = {}

    for model_name, feature_names in MODEL_FEATURE_SETS.items():
        pipeline = build_pipeline(feature_names, seed=seed)
        pipeline.fit(train.loc[:, feature_names], y_train)
        probabilities = pd.Series(
            pipeline.predict_proba(validation.loc[:, feature_names])[:, 1],
            index=validation.index,
        )
        metrics = calculate_metrics(
            y_validation,
            probabilities,
            threshold=threshold,
        )
        all_metrics.append({"model_name": model_name, "split": "validation", **metrics})

        predictions = (probabilities >= threshold).astype(int)
        prediction_tables.append(
            pd.DataFrame(
                {
                    "VariationID": validation["VariationID"].to_numpy(),
                    "model_name": model_name,
                    "split": "validation",
                    "target": y_validation.to_numpy(),
                    "probability_pathogenic": probabilities.to_numpy(),
                    "prediction": predictions.to_numpy(),
                }
            )
        )

        preprocessing = pipeline.named_steps["preprocessing"]
        classifier = pipeline.named_steps["classifier"]
        transformed_features = preprocessing.get_feature_names_out()
        coefficients = classifier.coef_[0]
        coefficient_tables.append(
            pd.DataFrame(
                {
                    "model_name": model_name,
                    "transformed_feature": transformed_features,
                    "coefficient": coefficients,
                    "absolute_coefficient": abs(coefficients),
                }
            ).sort_values("absolute_coefficient", ascending=False)
        )

        model_path = artifact_dir / f"{model_name}.joblib"
        joblib.dump(pipeline, model_path)
        model_details[model_name] = {
            "feature_columns": list(feature_names),
            "transformed_feature_count": len(transformed_features),
            "iterations": int(classifier.n_iter_[0]),
            "converged": bool(classifier.n_iter_[0] < classifier.max_iter),
            "artifact_path": str(model_path.resolve()),
            "artifact_sha256": sha256_file(model_path),
        }

    metrics_table = pd.DataFrame(all_metrics)
    predictions_table = pd.concat(prediction_tables, ignore_index=True)
    coefficients_table = pd.concat(coefficient_tables, ignore_index=True)
    metrics_path = report_dir / "validation_metrics.csv"
    predictions_path = report_dir / "validation_predictions.tsv.gz"
    coefficients_path = report_dir / "coefficients.csv"
    metrics_table.to_csv(metrics_path, index=False, float_format="%.10f")
    predictions_table.to_csv(
        predictions_path,
        sep="\t",
        index=False,
        compression="gzip",
        float_format="%.10f",
    )
    coefficients_table.to_csv(coefficients_path, index=False, float_format="%.10f")

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
            "identity_columns_excluded_from_model": ["VariationID"],
            "threshold": threshold,
            "threshold_policy": "fixed_in_advance_not_tuned",
        },
        "partition_counts": partition_counts,
        "models": model_details,
        "validation_metrics": {
            row["model_name"]: {
                key: value
                for key, value in row.items()
                if key not in {"model_name", "split"}
            }
            for row in all_metrics
        },
        "software": {
            package: version(package)
            for package in ("pandas", "scikit-learn", "numpy", "scipy", "joblib")
        },
        "reports": {
            "metrics": str(metrics_path.resolve()),
            "validation_predictions": str(predictions_path.resolve()),
            "coefficients": str(coefficients_path.resolve()),
        },
    }
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("Lojistik regresyon baseline eğitimi tamamlandı")
    print(f"Train satırı: {partition_counts['train']:,}")
    print(f"Validation satırı: {partition_counts['validation']:,}")
    print(f"Kilitli test satırı: {partition_counts['test']:,}")
    for row in all_metrics:
        print(
            f"{row['model_name']}: ROC-AUC={row['roc_auc']:.4f}, "
            f"F1={row['f1']:.4f}, dengeli doğruluk={row['balanced_accuracy']:.4f}"
        )
    print(f"Rapor: {manifest_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Train üzerinde lojistik regresyon eğitir, yalnızca validation üzerinde "
            "değerlendirir ve test bölümünü kilitli tutar."
        )
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
    train_baselines(
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
