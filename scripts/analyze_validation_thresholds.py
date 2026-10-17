"""Yalnızca validation tahminleriyle eşik ve ROC/PR analizi üret."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
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
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_recall_curve,
    precision_score,
    recall_score,
    roc_auc_score,
    roc_curve,
)


DEFAULT_PREDICTIONS = Path(
    "reports/generated/logistic_baseline/validation_predictions.tsv.gz"
)
DEFAULT_REPORT_DIR = Path(
    "reports/generated/logistic_baseline/threshold_analysis"
)
DEFAULT_SENSITIVITY_TARGETS = (0.80, 0.90, 0.95)
REQUIRED_COLUMNS = {
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


def load_validation_predictions(predictions_path: Path) -> pd.DataFrame:
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
    missing_columns = sorted(REQUIRED_COLUMNS.difference(predictions.columns))
    if missing_columns:
        raise ValueError("Eksik tahmin sütunları: " + ", ".join(missing_columns))
    if set(predictions["split"].dropna().unique()) != {"validation"}:
        raise ValueError("Eşik analizi yalnızca validation tahminlerini kabul eder.")
    if predictions[list(REQUIRED_COLUMNS)].isna().any().any():
        raise ValueError("Gerekli tahmin alanlarında boş değer bulundu.")
    if not set(predictions["target"].unique()).issubset({0, 1}):
        raise ValueError("Hedef değerleri yalnızca 0 veya 1 olmalıdır.")
    if not np.isfinite(predictions["probability_pathogenic"]).all():
        raise ValueError("Sonlu olmayan olasılık değeri bulundu.")
    if not predictions["probability_pathogenic"].between(0, 1).all():
        raise ValueError("Olasılıklar 0 ile 1 arasında olmalıdır.")

    reference_ids: set[str] | None = None
    reference_targets: pd.Series | None = None
    for model_name, model_rows in predictions.groupby("model_name", sort=True):
        if model_rows["VariationID"].duplicated().any():
            raise ValueError(f"{model_name} için tekrarlanan VariationID bulundu.")
        if set(model_rows["target"].unique()) != {0, 1}:
            raise ValueError(f"{model_name} için iki hedef sınıfı da bulunmalıdır.")

        indexed_targets = model_rows.set_index("VariationID")["target"].sort_index()
        model_ids = set(indexed_targets.index)
        if reference_ids is None:
            reference_ids = model_ids
            reference_targets = indexed_targets
        elif model_ids != reference_ids or not indexed_targets.equals(reference_targets):
            raise ValueError("Modeller aynı validation kimlikleri ve hedeflerini içermiyor.")

    return predictions


def threshold_metrics(
    targets: np.ndarray,
    probabilities: np.ndarray,
    threshold: float,
) -> dict[str, float | int]:
    predictions = (probabilities >= threshold).astype(int)
    true_negatives, false_positives, false_negatives, true_positives = confusion_matrix(
        targets, predictions, labels=[0, 1]
    ).ravel()
    predicted_positives = int(predictions.sum())
    row_count = len(targets)
    specificity = true_negatives / (true_negatives + false_positives)
    negative_predictive_value = (
        true_negatives / (true_negatives + false_negatives)
        if true_negatives + false_negatives
        else 0.0
    )
    return {
        "threshold": float(threshold),
        "row_count": row_count,
        "predicted_positive_count": predicted_positives,
        "predicted_positive_rate": predicted_positives / row_count,
        "sensitivity": recall_score(targets, predictions, zero_division=0),
        "specificity": specificity,
        "precision_positive_predictive_value": precision_score(
            targets, predictions, zero_division=0
        ),
        "negative_predictive_value": negative_predictive_value,
        "f1": f1_score(targets, predictions, zero_division=0),
        "accuracy": accuracy_score(targets, predictions),
        "balanced_accuracy": balanced_accuracy_score(targets, predictions),
        "true_negatives": int(true_negatives),
        "false_positives": int(false_positives),
        "false_negatives": int(false_negatives),
        "true_positives": int(true_positives),
    }


def choose_operating_thresholds(
    targets: np.ndarray,
    probabilities: np.ndarray,
    sensitivity_targets: Sequence[float],
) -> list[tuple[str, float, float | None]]:
    false_positive_rates, true_positive_rates, thresholds = roc_curve(
        targets, probabilities, drop_intermediate=False
    )
    finite = np.isfinite(thresholds)
    finite_indices = np.flatnonzero(finite)
    youden_values = true_positive_rates - false_positive_rates
    best_youden = np.max(youden_values[finite_indices])
    youden_candidates = finite_indices[
        np.isclose(youden_values[finite_indices], best_youden)
    ]
    youden_index = youden_candidates[np.argmax(thresholds[youden_candidates])]

    choices: list[tuple[str, float, float | None]] = [
        ("fixed_0.50", 0.5, None),
        ("maximum_youden_j", float(thresholds[youden_index]), None),
    ]
    for target_sensitivity in sensitivity_targets:
        eligible = finite_indices[true_positive_rates[finite_indices] >= target_sensitivity]
        if len(eligible) == 0:
            raise ValueError(
                f"{target_sensitivity:.0%} duyarlılık hedefini karşılayan eşik yok."
            )
        specificities = 1 - false_positive_rates[eligible]
        best_specificity = np.max(specificities)
        candidates = eligible[np.isclose(specificities, best_specificity)]
        selected = candidates[np.argmax(thresholds[candidates])]
        choices.append(
            (
                f"sensitivity_at_least_{target_sensitivity:.2f}",
                float(thresholds[selected]),
                float(target_sensitivity),
            )
        )
    return choices


def plot_roc_and_precision_recall(
    predictions: pd.DataFrame,
    output_path: Path,
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(12, 5))
    prevalence: float | None = None
    for model_name, rows in predictions.groupby("model_name", sort=True):
        targets = rows["target"].to_numpy()
        probabilities = rows["probability_pathogenic"].to_numpy()
        false_positive_rates, true_positive_rates, _ = roc_curve(targets, probabilities)
        precision, recall, _ = precision_recall_curve(targets, probabilities)
        roc_auc = roc_auc_score(targets, probabilities)
        average_precision = average_precision_score(targets, probabilities)
        axes[0].plot(
            false_positive_rates,
            true_positive_rates,
            label=f"{model_name} (AUC={roc_auc:.4f})",
        )
        axes[1].plot(
            recall,
            precision,
            label=f"{model_name} (AP={average_precision:.4f})",
        )
        prevalence = float(targets.mean())

    axes[0].plot([0, 1], [0, 1], linestyle="--", color="gray", label="Rastgele")
    axes[0].set(title="Validation ROC eğrisi", xlabel="1 - Özgüllük", ylabel="Duyarlılık")
    if prevalence is not None:
        axes[1].axhline(
            prevalence,
            linestyle="--",
            color="gray",
            label=f"Sınıf prevalansı ({prevalence:.3f})",
        )
    axes[1].set(
        title="Validation precision-recall eğrisi",
        xlabel="Duyarlılık (Recall)",
        ylabel="Kesinlik (Precision)",
    )
    for axis in axes:
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    figure.tight_layout()
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)


def plot_threshold_tradeoffs(sweep: pd.DataFrame, output_path: Path) -> None:
    model_names = sorted(sweep["model_name"].unique())
    figure, axes = plt.subplots(
        len(model_names),
        1,
        figsize=(9, 4 * len(model_names)),
        squeeze=False,
        sharex=True,
    )
    for axis, model_name in zip(axes[:, 0], model_names):
        rows = sweep.loc[sweep["model_name"] == model_name]
        for metric in ("sensitivity", "specificity", "precision_positive_predictive_value", "f1"):
            axis.plot(rows["threshold"], rows[metric], label=metric)
        axis.axvline(0.5, linestyle="--", color="black", alpha=0.6, label="eşik=0.5")
        axis.set(title=model_name, ylabel="Metrik değeri")
        axis.set_ylim(0, 1)
        axis.grid(alpha=0.25)
        axis.legend(fontsize=8)
    axes[-1, 0].set_xlabel("Karar eşiği")
    figure.suptitle("Validation karar eşiği dengeleri", y=1.01)
    figure.tight_layout()
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)


def analyze_thresholds(
    predictions_path: Path,
    report_dir: Path,
    *,
    sensitivity_targets: Sequence[float] = DEFAULT_SENSITIVITY_TARGETS,
) -> dict[str, object]:
    if not sensitivity_targets:
        raise ValueError("En az bir duyarlılık hedefi verilmelidir.")
    if any(not 0 < target < 1 for target in sensitivity_targets):
        raise ValueError("Duyarlılık hedefleri 0 ile 1 arasında olmalıdır.")

    predictions = load_validation_predictions(predictions_path)
    report_dir.mkdir(parents=True, exist_ok=True)
    threshold_grid = np.round(np.arange(0.01, 1.00, 0.01), 2)
    sweep_rows: list[dict[str, object]] = []
    operating_rows: list[dict[str, object]] = []
    roc_rows: list[dict[str, object]] = []
    precision_recall_rows: list[dict[str, object]] = []
    model_summary: dict[str, object] = {}

    for model_name, rows in predictions.groupby("model_name", sort=True):
        targets = rows["target"].to_numpy()
        probabilities = rows["probability_pathogenic"].to_numpy()
        for threshold in threshold_grid:
            sweep_rows.append(
                {
                    "model_name": model_name,
                    **threshold_metrics(targets, probabilities, float(threshold)),
                }
            )

        for policy, threshold, sensitivity_target in choose_operating_thresholds(
            targets,
            probabilities,
            sensitivity_targets,
        ):
            operating_rows.append(
                {
                    "model_name": model_name,
                    "policy": policy,
                    "sensitivity_target": sensitivity_target,
                    **threshold_metrics(targets, probabilities, threshold),
                }
            )

        false_positive_rates, true_positive_rates, roc_thresholds = roc_curve(
            targets, probabilities, drop_intermediate=False
        )
        for false_positive_rate, true_positive_rate, threshold in zip(
            false_positive_rates, true_positive_rates, roc_thresholds
        ):
            roc_rows.append(
                {
                    "model_name": model_name,
                    "threshold": threshold if math.isfinite(threshold) else None,
                    "false_positive_rate": false_positive_rate,
                    "true_positive_rate_sensitivity": true_positive_rate,
                    "specificity": 1 - false_positive_rate,
                }
            )

        precision, recall, pr_thresholds = precision_recall_curve(targets, probabilities)
        padded_thresholds = np.append(pr_thresholds, np.nan)
        for precision_value, recall_value, threshold in zip(
            precision, recall, padded_thresholds
        ):
            precision_recall_rows.append(
                {
                    "model_name": model_name,
                    "threshold": threshold if math.isfinite(threshold) else None,
                    "precision": precision_value,
                    "recall_sensitivity": recall_value,
                }
            )

        model_summary[model_name] = {
            "validation_rows": len(rows),
            "negative_rows": int((targets == 0).sum()),
            "positive_rows": int((targets == 1).sum()),
            "roc_auc": roc_auc_score(targets, probabilities),
            "average_precision": average_precision_score(targets, probabilities),
        }

    sweep = pd.DataFrame(sweep_rows)
    operating_points = pd.DataFrame(operating_rows)
    roc_points = pd.DataFrame(roc_rows)
    precision_recall_points = pd.DataFrame(precision_recall_rows)
    sweep_path = report_dir / "threshold_sweep.csv"
    operating_path = report_dir / "operating_points.csv"
    roc_path = report_dir / "roc_curve.csv"
    precision_recall_path = report_dir / "precision_recall_curve.csv"
    curves_plot_path = report_dir / "validation_roc_pr_curves.png"
    tradeoffs_plot_path = report_dir / "validation_threshold_tradeoffs.png"
    sweep.to_csv(sweep_path, index=False, float_format="%.10f")
    operating_points.to_csv(operating_path, index=False, float_format="%.10f")
    roc_points.to_csv(roc_path, index=False, float_format="%.10f")
    precision_recall_points.to_csv(
        precision_recall_path, index=False, float_format="%.10f"
    )
    plot_roc_and_precision_recall(predictions, curves_plot_path)
    plot_threshold_tradeoffs(sweep, tradeoffs_plot_path)

    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "path": str(predictions_path.resolve()),
            "sha256": sha256_file(predictions_path),
            "observed_splits": sorted(predictions["split"].unique().tolist()),
        },
        "data_policy": {
            "analysis_partition": "validation",
            "test_data_used": False,
            "selection_warning": (
                "Validation üzerinde seçilen eşik yalnızca bir kez, nihai model "
                "sabitlendikten sonra test üzerinde değerlendirilmelidir."
            ),
        },
        "sensitivity_targets": list(sensitivity_targets),
        "models": model_summary,
        "outputs": {
            "threshold_sweep": str(sweep_path.resolve()),
            "operating_points": str(operating_path.resolve()),
            "roc_curve": str(roc_path.resolve()),
            "precision_recall_curve": str(precision_recall_path.resolve()),
            "roc_pr_plot": str(curves_plot_path.resolve()),
            "threshold_tradeoffs_plot": str(tradeoffs_plot_path.resolve()),
        },
    }
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("Validation eşik analizi tamamlandı")
    for row in operating_rows:
        if row["model_name"] != "without_missingness_flag":
            continue
        print(
            f"{row['policy']}: eşik={row['threshold']:.4f}, "
            f"duyarlılık={row['sensitivity']:.4f}, "
            f"özgüllük={row['specificity']:.4f}, "
            f"yanlış negatif={row['false_negatives']:,}"
        )
    print("Test verisi kullanılmadı.")
    print(f"Rapor: {manifest_path.resolve()}")
    return manifest


def parse_sensitivity_targets(value: str) -> tuple[float, ...]:
    try:
        targets = tuple(float(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as error:
        raise argparse.ArgumentTypeError("Duyarlılık hedefleri sayı olmalıdır.") from error
    if not targets or any(not 0 < target < 1 for target in targets):
        raise argparse.ArgumentTypeError(
            "Duyarlılık hedefleri virgülle ayrılmış 0-1 arası sayılar olmalıdır."
        )
    return targets


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Yalnızca validation tahminleriyle karar eşiği analizi üretir."
    )
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument(
        "--sensitivity-targets",
        type=parse_sensitivity_targets,
        default=DEFAULT_SENSITIVITY_TARGETS,
        help="Virgülle ayrılmış hedefler (varsayılan: 0.80,0.90,0.95)",
    )
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    analyze_thresholds(
        args.predictions,
        args.report_dir,
        sensitivity_targets=args.sensitivity_targets,
    )


if __name__ == "__main__":
    main()
