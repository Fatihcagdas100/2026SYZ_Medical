"""Tam HGVS özelliklerinde lojistik ve ağaç modelini adil biçimde karşılaştır."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import joblib
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OrdinalEncoder

try:
    from scripts.run_hgvs_ablation import calculate_subgroup_metrics
    from scripts.train_logistic_baseline import (
        CATEGORICAL_FEATURES,
        HGVS_CATEGORICAL_FEATURES,
        HGVS_NUMERIC_FEATURES,
        NUMERIC_FEATURES,
        calculate_metrics,
        load_development_data,
        sha256_file,
        train_baselines,
    )
except ModuleNotFoundError:  # Doğrudan ``python scripts/...`` çalıştırması.
    from run_hgvs_ablation import calculate_subgroup_metrics
    from train_logistic_baseline import (
        CATEGORICAL_FEATURES,
        HGVS_CATEGORICAL_FEATURES,
        HGVS_NUMERIC_FEATURES,
        NUMERIC_FEATURES,
        calculate_metrics,
        load_development_data,
        sha256_file,
        train_baselines,
    )


DEFAULT_FEATURES = Path("data/processed/clinvar_hgvs_features.tsv.gz")
DEFAULT_SPLIT = Path("data/processed/baseline_gene_disjoint_split.tsv")
DEFAULT_REPORT_DIR = Path("reports/generated/tree_model_comparison")
DEFAULT_ARTIFACT_DIR = Path("artifacts/tree_model_comparison")
DEFAULT_SEED = 2027
DEFAULT_THRESHOLD = 0.5
DEFAULT_CHUNK_SIZE = 100_000
PRIMARY_LOGISTIC_MODEL = "without_missingness_flag"
CATEGORICAL_COLUMNS = (*CATEGORICAL_FEATURES, *HGVS_CATEGORICAL_FEATURES)
NUMERIC_COLUMNS = (*NUMERIC_FEATURES, *HGVS_NUMERIC_FEATURES)
FEATURE_COLUMNS = (*CATEGORICAL_COLUMNS, *NUMERIC_COLUMNS)
TREE_PARAMETERS = {
    "learning_rate": 0.08,
    "max_iter": 200,
    "max_leaf_nodes": 31,
    "min_samples_leaf": 50,
    "l2_regularization": 1.0,
    "early_stopping": False,
}
SUBGROUPS = (
    "all_variants",
    "missense",
    "synonymous",
    "frameshift",
    "stop_gained_or_changed",
    "protein_annotation_present",
    "protein_annotation_absent",
)


def build_tree_pipeline(*, seed: int) -> Pipeline:
    categorical_pipeline = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="most_frequent")),
            (
                "ordinal",
                OrdinalEncoder(
                    handle_unknown="use_encoded_value",
                    unknown_value=-1,
                    encoded_missing_value=-1,
                ),
            ),
        ]
    )
    numeric_pipeline = Pipeline(
        steps=[("imputer", SimpleImputer(strategy="median"))]
    )
    preprocessing = ColumnTransformer(
        transformers=[
            ("categorical", categorical_pipeline, list(CATEGORICAL_COLUMNS)),
            ("numeric", numeric_pipeline, list(NUMERIC_COLUMNS)),
        ],
        sparse_threshold=0.0,
        verbose_feature_names_out=False,
    )
    categorical_mask = [True] * len(CATEGORICAL_COLUMNS) + [False] * len(
        NUMERIC_COLUMNS
    )
    classifier = HistGradientBoostingClassifier(
        **TREE_PARAMETERS,
        categorical_features=categorical_mask,
        random_state=seed,
    )
    return Pipeline(
        steps=[("preprocessing", preprocessing), ("classifier", classifier)]
    )


def add_subgroup_rows(
    rows: list[dict[str, object]],
    frame: pd.DataFrame,
    *,
    model_name: str,
    threshold: float,
) -> None:
    masks = (
        ("all_variants", pd.Series(True, index=frame.index)),
        ("missense", frame["protein_effect"].eq("missense")),
        ("synonymous", frame["protein_effect"].eq("synonymous")),
        ("frameshift", frame["protein_effect"].eq("frameshift")),
        (
            "stop_gained_or_changed",
            frame["protein_effect"].eq("stop_gained_or_changed"),
        ),
        (
            "protein_annotation_present",
            frame["has_protein_annotation"].eq(1),
        ),
        (
            "protein_annotation_absent",
            frame["has_protein_annotation"].eq(0),
        ),
    )
    for group_name, mask in masks:
        rows.append(
            {
                "model": model_name,
                **calculate_subgroup_metrics(
                    frame.loc[mask], group_name=group_name, threshold=threshold
                ),
            }
        )


def build_delta_table(subgroups: pd.DataFrame) -> pd.DataFrame:
    metrics = (
        "roc_auc",
        "average_precision",
        "precision",
        "recall_sensitivity",
        "specificity",
        "f1",
    )
    rows: list[dict[str, object]] = []
    for group_name in SUBGROUPS:
        group = subgroups.loc[subgroups["group"] == group_name].set_index("model")
        for metric in metrics:
            rows.append(
                {
                    "group": group_name,
                    "metric": metric,
                    "tree_minus_logistic": float(
                        group.loc["hist_gradient_boosting", metric]
                        - group.loc["logistic_regression", metric]
                    ),
                }
            )
    return pd.DataFrame(rows)


def compare_tree_model(
    features_path: Path,
    split_path: Path,
    report_dir: Path,
    artifact_dir: Path,
    *,
    seed: int = DEFAULT_SEED,
    threshold: float = DEFAULT_THRESHOLD,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, object]:
    if not 0 < threshold < 1:
        raise ValueError("Eşik 0 ile 1 arasında olmalıdır.")
    report_dir.mkdir(parents=True, exist_ok=True)
    artifact_dir.mkdir(parents=True, exist_ok=True)

    logistic_report_dir = report_dir / "logistic_reference"
    logistic_artifact_dir = artifact_dir / "logistic_reference"
    logistic_manifest = train_baselines(
        features_path,
        split_path,
        logistic_report_dir,
        logistic_artifact_dir,
        seed=seed,
        threshold=threshold,
        chunk_size=chunk_size,
        feature_profile="hgvs",
    )
    if logistic_manifest["data_policy"]["test_partition"] != (
        "locked_not_loaded_into_development_dataframe"
    ):
        raise ValueError("Lojistik referans test bölümünü kilitli tutmadı.")

    development, partition_counts = load_development_data(
        features_path,
        split_path,
        feature_columns=FEATURE_COLUMNS,
        categorical_features=CATEGORICAL_COLUMNS,
        chunk_size=chunk_size,
    )
    train = development.loc[development["split"] == "train"].copy()
    validation = development.loc[development["split"] == "validation"].copy()
    y_train = train["target"].astype(int)
    y_validation = validation["target"].astype(int)

    tree_pipeline = build_tree_pipeline(seed=seed)
    tree_pipeline.fit(train.loc[:, FEATURE_COLUMNS], y_train)
    tree_probabilities = pd.Series(
        tree_pipeline.predict_proba(validation.loc[:, FEATURE_COLUMNS])[:, 1],
        index=validation.index,
    )
    tree_metrics = calculate_metrics(
        y_validation, tree_probabilities, threshold=threshold
    )
    tree_predictions = pd.DataFrame(
        {
            "VariationID": validation["VariationID"].astype(str).to_numpy(),
            "model": "hist_gradient_boosting",
            "split": "validation",
            "target": y_validation.to_numpy(),
            "probability_pathogenic": tree_probabilities.to_numpy(),
            "prediction": (tree_probabilities >= threshold).astype(int).to_numpy(),
            "protein_effect": validation["protein_effect"].to_numpy(),
            "has_protein_annotation": validation[
                "has_protein_annotation"
            ].to_numpy(),
        }
    )

    logistic_predictions = pd.read_csv(
        logistic_report_dir / "validation_predictions.tsv.gz",
        sep="\t",
        dtype={"VariationID": "string", "model_name": "string", "split": "string"},
    )
    logistic_predictions = logistic_predictions.loc[
        logistic_predictions["model_name"] == PRIMARY_LOGISTIC_MODEL
    ].copy()
    if set(logistic_predictions["split"].dropna()) != {"validation"}:
        raise ValueError("Lojistik referans yalnızca validation tahmini içermelidir.")
    logistic_predictions = logistic_predictions.merge(
        validation.loc[
            :, ["VariationID", "protein_effect", "has_protein_annotation"]
        ].assign(VariationID=lambda frame: frame["VariationID"].astype(str)),
        on="VariationID",
        how="left",
        validate="one_to_one",
    )
    logistic_predictions = logistic_predictions.rename(
        columns={"model_name": "source_model_name"}
    )
    logistic_predictions["model"] = "logistic_regression"
    if set(logistic_predictions["VariationID"]) != set(tree_predictions["VariationID"]):
        raise ValueError("Ağaç ve lojistik validation kimlikleri eşleşmiyor.")
    if not logistic_predictions["target"].astype(int).equals(
        tree_predictions.set_index("VariationID")
        .loc[logistic_predictions["VariationID"], "target"]
        .reset_index(drop=True)
        .astype(int)
    ):
        raise ValueError("Ağaç ve lojistik validation hedefleri eşleşmiyor.")
    logistic_metrics = calculate_metrics(
        logistic_predictions["target"].astype(int),
        logistic_predictions["probability_pathogenic"].astype(float),
        threshold=threshold,
    )

    comparison = pd.DataFrame(
        [
            {"model": "logistic_regression", **logistic_metrics},
            {"model": "hist_gradient_boosting", **tree_metrics},
        ]
    )
    comparison_path = report_dir / "validation_model_comparison.csv"
    comparison.to_csv(comparison_path, index=False, float_format="%.10f")

    subgroup_rows: list[dict[str, object]] = []
    add_subgroup_rows(
        subgroup_rows,
        logistic_predictions,
        model_name="logistic_regression",
        threshold=threshold,
    )
    add_subgroup_rows(
        subgroup_rows,
        tree_predictions,
        model_name="hist_gradient_boosting",
        threshold=threshold,
    )
    subgroups = pd.DataFrame(subgroup_rows)
    subgroup_path = report_dir / "validation_subgroup_comparison.csv"
    subgroups.to_csv(subgroup_path, index=False, float_format="%.10f")
    deltas = build_delta_table(subgroups)
    delta_path = report_dir / "validation_metric_deltas.csv"
    deltas.to_csv(delta_path, index=False, float_format="%.10f")

    prediction_path = report_dir / "tree_validation_predictions.tsv.gz"
    tree_predictions.drop(
        columns=["protein_effect", "has_protein_annotation"]
    ).to_csv(
        prediction_path,
        sep="\t",
        index=False,
        compression="gzip",
        float_format="%.10f",
    )
    model_path = artifact_dir / "hist_gradient_boosting.joblib"
    joblib.dump(tree_pipeline, model_path)

    plot_groups = ["all_variants", "missense"]
    plot_metrics = ["roc_auc", "average_precision", "recall_sensitivity"]
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for axis, group_name in zip(axes, plot_groups, strict=True):
        subset = subgroups.loc[subgroups["group"] == group_name].set_index("model")
        positions = range(len(plot_metrics))
        width = 0.36
        axis.bar(
            [position - width / 2 for position in positions],
            subset.loc["logistic_regression", plot_metrics],
            width,
            label="Lojistik",
        )
        axis.bar(
            [position + width / 2 for position in positions],
            subset.loc["hist_gradient_boosting", plot_metrics],
            width,
            label="HistGradientBoosting",
        )
        axis.set_xticks(list(positions), ["ROC-AUC", "AP", "Recall"])
        axis.set_title(group_name)
        axis.set_ylim(0, 1.02)
        axis.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Validation skoru")
    axes[1].legend(loc="lower right")
    figure.suptitle("Aynı HGVS özellikleriyle model karşılaştırması")
    figure.tight_layout()
    plot_path = report_dir / "tree_vs_logistic_overall_missense.png"
    figure.savefig(plot_path, dpi=180)
    plt.close(figure)

    classifier = tree_pipeline.named_steps["classifier"]
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
            "feature_profile": "hgvs_without_vcf_missingness_flag",
            "same_features_split_seed_and_threshold": True,
            "hyperparameter_search": False,
            "comparison_role": "single_fixed_tree_baseline_not_model_selection_sweep",
        },
        "partition_counts": partition_counts,
        "features": {
            "categorical": list(CATEGORICAL_COLUMNS),
            "numeric": list(NUMERIC_COLUMNS),
        },
        "tree_model": {
            "estimator": "HistGradientBoostingClassifier",
            "parameters": {**TREE_PARAMETERS, "random_state": seed},
            "completed_iterations": int(classifier.n_iter_),
            "artifact_path": str(model_path.resolve()),
            "artifact_sha256": sha256_file(model_path),
        },
        "validation_metrics": {
            row["model"]: {
                key: value for key, value in row.items() if key != "model"
            }
            for row in comparison.to_dict("records")
        },
        "reports": {
            "model_comparison": str(comparison_path.resolve()),
            "subgroup_comparison": str(subgroup_path.resolve()),
            "metric_deltas": str(delta_path.resolve()),
            "tree_predictions": str(prediction_path.resolve()),
            "plot": str(plot_path.resolve()),
        },
        "software": {
            package: version(package)
            for package in ("pandas", "scikit-learn", "numpy", "joblib")
        },
    }

    report_path = report_dir / "tree_model_comparison.md"
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write("# Ağaç modeli karşılaştırması\n\n")
        handle.write(
            "Modeller aynı tam HGVS özellikleri, gene-disjoint train/validation "
            "ayrımı, seed ve sabit 0.5 eşikle karşılaştırıldı. Test yüklenmedi. "
            "Ağaç ayarları tek koşudan önce sabitlendi; hiperparametre taraması "
            "yapılmadı.\n\n"
        )
        handle.write("| Model | ROC-AUC | AP | Precision | Recall | Özgüllük | F1 |\n")
        handle.write("| --- | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for row in comparison.to_dict("records"):
            handle.write(
                f"| {row['model']} | {row['roc_auc']:.4f} | "
                f"{row['average_precision']:.4f} | {row['precision']:.4f} | "
                f"{row['recall_sensitivity']:.4f} | {row['specificity']:.4f} | "
                f"{row['f1']:.4f} |\n"
            )
        handle.write("\n## Alt gruplar\n\n")
        handle.write("| Model | Grup | n | ROC-AUC | AP | Precision | Recall | Özgüllük |\n")
        handle.write("| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for row in subgroups.to_dict("records"):
            handle.write(
                f"| {row['model']} | {row['group']} | {row['row_count']:,} | "
                f"{row['roc_auc']:.4f} | {row['average_precision']:.4f} | "
                f"{row['precision']:.4f} | {row['recall_sensitivity']:.4f} | "
                f"{row['specificity']:.4f} |\n"
            )
        handle.write(
            "\nSonuçlar yalnızca validation karşılaştırmasıdır. Bu tek ağaç "
            "yapılandırması klinik geçerlilik veya nihai model seçimi anlamına gelmez.\n"
        )
    manifest["reports"]["report"] = str(report_path.resolve())
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("Ağaç modeli karşılaştırması tamamlandı")
    for row in comparison.to_dict("records"):
        print(
            f"{row['model']}: ROC-AUC={row['roc_auc']:.4f}, "
            f"AP={row['average_precision']:.4f}, F1={row['f1']:.4f}"
        )
    missense = subgroups.loc[subgroups["group"] == "missense"]
    for row in missense.to_dict("records"):
        print(
            f"missense/{row['model']}: ROC-AUC={row['roc_auc']:.4f}, "
            f"recall={row['recall_sensitivity']:.4f}"
        )
    print(f"Rapor: {report_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Aynı HGVS özellikleriyle lojistik ve ağaç modelini karşılaştırır."
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
    compare_tree_model(
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
