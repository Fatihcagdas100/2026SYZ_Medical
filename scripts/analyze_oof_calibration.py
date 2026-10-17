"""Train-only OOF tahminlerinin olasılık kalibrasyonunu analiz et."""

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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, log_loss


DEFAULT_PREDICTIONS = Path(
    "reports/generated/hgvs_gene_disjoint_cv/train_out_of_fold_predictions.tsv.gz"
)
DEFAULT_REPORT_DIR = Path(
    "reports/generated/hgvs_gene_disjoint_cv/calibration"
)
DEFAULT_N_BINS = 10
PROBABILITY_EPSILON = 1e-6
REQUIRED_COLUMNS = {
    "VariationID",
    "gene_group",
    "source_partition",
    "fold",
    "target",
    "probability_pathogenic",
}


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def load_oof_predictions(predictions_path: Path) -> pd.DataFrame:
    if not predictions_path.is_file():
        raise FileNotFoundError(f"OOF tahmin dosyası bulunamadı: {predictions_path}")
    predictions = pd.read_csv(
        predictions_path,
        sep="\t",
        dtype={
            "VariationID": "string",
            "gene_group": "string",
            "source_partition": "string",
            "fold": "Int16",
            "target": "Int8",
            "probability_pathogenic": "float64",
        },
    )
    missing_columns = sorted(REQUIRED_COLUMNS.difference(predictions.columns))
    if missing_columns:
        raise ValueError("Eksik OOF sütunları: " + ", ".join(missing_columns))
    if predictions[list(REQUIRED_COLUMNS)].isna().any().any():
        raise ValueError("Gerekli OOF alanlarında boş değer bulundu.")
    if set(predictions["source_partition"].unique()) != {"train"}:
        raise ValueError("Kalibrasyon analizi yalnızca train OOF tahminlerini kabul eder.")
    if predictions["VariationID"].duplicated().any():
        raise ValueError("OOF tahminlerinde tekrarlanan VariationID bulundu.")
    if set(predictions["target"].astype(int).unique()) != {0, 1}:
        raise ValueError("OOF tahminleri iki hedef sınıfını da içermelidir.")
    probabilities = predictions["probability_pathogenic"]
    if not np.isfinite(probabilities).all() or not probabilities.between(0, 1).all():
        raise ValueError("OOF olasılıkları sonlu ve 0 ile 1 arasında olmalıdır.")
    if predictions.groupby("gene_group")["fold"].nunique().max() != 1:
        raise ValueError("Bir gen grubu birden fazla OOF fold'unda bulundu.")
    if predictions["fold"].nunique() < 2:
        raise ValueError("Fold kararlılığı için en az iki OOF fold'u gereklidir.")
    for fold, rows in predictions.groupby("fold"):
        if set(rows["target"].astype(int).unique()) != {0, 1}:
            raise ValueError(f"Fold {fold} iki hedef sınıfını da içermiyor.")
    return predictions


def reliability_table(
    targets: np.ndarray,
    probabilities: np.ndarray,
    *,
    n_bins: int,
    strategy: str,
) -> pd.DataFrame:
    if n_bins < 2:
        raise ValueError("n_bins en az 2 olmalıdır.")
    if strategy == "quantile":
        assignments = pd.qcut(probabilities, q=n_bins, duplicates="drop")
    elif strategy == "fixed_width":
        assignments = pd.cut(
            probabilities,
            bins=np.linspace(0, 1, n_bins + 1),
            include_lowest=True,
        )
    else:
        raise ValueError(f"Bilinmeyen dilim stratejisi: {strategy}")

    frame = pd.DataFrame(
        {"target": targets, "probability": probabilities, "assignment": assignments}
    )
    grouped = (
        frame.groupby("assignment", observed=True)
        .agg(
            row_count=("target", "size"),
            positive_count=("target", "sum"),
            mean_probability=("probability", "mean"),
            observed_pathogenic_rate=("target", "mean"),
            minimum_probability=("probability", "min"),
            maximum_probability=("probability", "max"),
        )
        .reset_index(drop=True)
    )
    grouped.insert(0, "bin", np.arange(1, len(grouped) + 1))
    grouped.insert(0, "strategy", strategy)
    grouped["calibration_gap"] = (
        grouped["mean_probability"] - grouped["observed_pathogenic_rate"]
    )
    grouped["absolute_calibration_gap"] = grouped["calibration_gap"].abs()
    return grouped


def calibration_line(
    targets: np.ndarray, probabilities: np.ndarray
) -> tuple[float, float]:
    clipped = np.clip(probabilities, PROBABILITY_EPSILON, 1 - PROBABILITY_EPSILON)
    logits = np.log(clipped / (1 - clipped)).reshape(-1, 1)
    model = LogisticRegression(C=np.inf, solver="lbfgs", max_iter=1_000)
    model.fit(logits, targets)
    return float(model.intercept_[0]), float(model.coef_[0, 0])


def calibration_metrics(
    targets: np.ndarray,
    probabilities: np.ndarray,
    *,
    n_bins: int,
) -> dict[str, float | int]:
    quantile_table = reliability_table(
        targets, probabilities, n_bins=n_bins, strategy="quantile"
    )
    weights = quantile_table["row_count"] / len(targets)
    intercept, slope = calibration_line(targets, probabilities)
    return {
        "row_count": len(targets),
        "positive_count": int(targets.sum()),
        "prevalence": float(targets.mean()),
        "mean_predicted_probability": float(probabilities.mean()),
        "brier_score": float(brier_score_loss(targets, probabilities)),
        "log_loss": float(log_loss(targets, probabilities)),
        "expected_calibration_error": float(
            (weights * quantile_table["absolute_calibration_gap"]).sum()
        ),
        "maximum_calibration_error": float(
            quantile_table["absolute_calibration_gap"].max()
        ),
        "calibration_intercept": intercept,
        "calibration_slope": slope,
        "quantile_bins_requested": n_bins,
        "quantile_bins_observed": len(quantile_table),
    }


def plot_calibration(
    predictions: pd.DataFrame,
    overall_quantile: pd.DataFrame,
    fixed_width: pd.DataFrame,
    output_path: Path,
    *,
    n_bins: int,
) -> None:
    figure, axes = plt.subplots(2, 2, figsize=(12, 10))

    axes[0, 0].plot([0, 1], [0, 1], linestyle="--", color="gray", label="İdeal")
    axes[0, 0].plot(
        overall_quantile["mean_probability"],
        overall_quantile["observed_pathogenic_rate"],
        marker="o",
        label="Birleşik OOF",
    )
    axes[0, 0].set_title(f"Eşit-frekans güvenilirlik eğrisi ({n_bins} dilim)")

    axes[0, 1].plot([0, 1], [0, 1], linestyle="--", color="gray", label="İdeal")
    for fold, rows in predictions.groupby("fold", sort=True):
        table = reliability_table(
            rows["target"].to_numpy(dtype=int),
            rows["probability_pathogenic"].to_numpy(dtype=float),
            n_bins=n_bins,
            strategy="quantile",
        )
        axes[0, 1].plot(
            table["mean_probability"],
            table["observed_pathogenic_rate"],
            marker=".",
            linewidth=1,
            label=f"Fold {int(fold)}",
        )
    axes[0, 1].set_title("Fold bazında güvenilirlik")

    for target, label, color in (
        (0, "Benign", "tab:blue"),
        (1, "Pathogenic", "tab:orange"),
    ):
        values = predictions.loc[
            predictions["target"] == target, "probability_pathogenic"
        ]
        axes[1, 0].hist(
            values,
            bins=np.linspace(0, 1, 31),
            alpha=0.6,
            label=label,
            color=color,
        )
    axes[1, 0].set_yscale("log")
    axes[1, 0].set_title("Hedef sınıfa göre OOF olasılık dağılımı")
    axes[1, 0].set_xlabel("P(Pathogenic)")
    axes[1, 0].set_ylabel("Varyant sayısı (log)")
    axes[1, 0].grid(alpha=0.2)
    axes[1, 0].legend()

    colors = np.where(fixed_width["calibration_gap"] >= 0, "tab:red", "tab:blue")
    axes[1, 1].bar(
        fixed_width["mean_probability"],
        fixed_width["calibration_gap"],
        width=0.07,
        color=colors,
        alpha=0.8,
    )
    axes[1, 1].axhline(0, color="black", linewidth=1)
    axes[1, 1].set_title("Sabit genişlikli dilimlerde kalibrasyon farkı")
    axes[1, 1].set_xlabel("Ortalama tahmin olasılığı")
    axes[1, 1].set_ylabel("Tahmin − gözlenen oran")
    axes[1, 1].grid(axis="y", alpha=0.2)

    for axis in (axes[0, 0], axes[0, 1]):
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.set_xlabel("Ortalama tahmin olasılığı")
        axis.set_ylabel("Gözlenen pathogenic oranı")
        axis.grid(alpha=0.2)
        axis.legend(fontsize=8)

    figure.suptitle("Train-only gene-disjoint OOF kalibrasyon analizi", fontsize=14)
    figure.tight_layout()
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)


def analyze_calibration(
    predictions_path: Path,
    report_dir: Path,
    *,
    n_bins: int = DEFAULT_N_BINS,
) -> dict[str, object]:
    predictions = load_oof_predictions(predictions_path)
    targets = predictions["target"].to_numpy(dtype=int)
    probabilities = predictions["probability_pathogenic"].to_numpy(dtype=float)

    overall_quantile = reliability_table(
        targets, probabilities, n_bins=n_bins, strategy="quantile"
    )
    fixed_width = reliability_table(
        targets, probabilities, n_bins=n_bins, strategy="fixed_width"
    )
    overall_metrics = calibration_metrics(targets, probabilities, n_bins=n_bins)
    fold_rows = []
    for fold, rows in predictions.groupby("fold", sort=True):
        fold_rows.append(
            {
                "fold": int(fold),
                **calibration_metrics(
                    rows["target"].to_numpy(dtype=int),
                    rows["probability_pathogenic"].to_numpy(dtype=float),
                    n_bins=n_bins,
                ),
            }
        )
    fold_metrics = pd.DataFrame(fold_rows)

    report_dir.mkdir(parents=True, exist_ok=True)
    quantile_path = report_dir / "quantile_reliability.csv"
    fixed_width_path = report_dir / "fixed_width_reliability.csv"
    fold_metrics_path = report_dir / "fold_calibration_metrics.csv"
    plot_path = report_dir / "calibration_reliability.png"
    overall_quantile.to_csv(quantile_path, index=False, float_format="%.10f")
    fixed_width.to_csv(fixed_width_path, index=False, float_format="%.10f")
    fold_metrics.to_csv(fold_metrics_path, index=False, float_format="%.10f")
    plot_calibration(
        predictions,
        overall_quantile,
        fixed_width,
        plot_path,
        n_bins=n_bins,
    )

    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": {
            "path": str(predictions_path.resolve()),
            "sha256": sha256_file(predictions_path),
        },
        "data_policy": {
            "source_partition": "train_out_of_fold",
            "fixed_validation_used": False,
            "test_used": False,
            "recalibration_applied": False,
        },
        "method": {
            "reliability_bins": n_bins,
            "primary_binning": "equal_frequency_quantile",
            "secondary_binning": "fixed_width",
            "calibration_line": "unregularized_logistic_regression_on_clipped_logit",
            "ideal_calibration_intercept": 0.0,
            "ideal_calibration_slope": 1.0,
        },
        "integrity": {
            "oof_rows": len(predictions),
            "unique_variation_ids": int(predictions["VariationID"].nunique()),
            "gene_groups": int(predictions["gene_group"].nunique()),
            "folds": int(predictions["fold"].nunique()),
            "gene_groups_in_multiple_folds": int(
                predictions.groupby("gene_group")["fold"].nunique().gt(1).sum()
            ),
        },
        "overall_metrics": overall_metrics,
        "outputs": {
            "quantile_reliability": str(quantile_path.resolve()),
            "fixed_width_reliability": str(fixed_width_path.resolve()),
            "fold_calibration_metrics": str(fold_metrics_path.resolve()),
            "calibration_plot": str(plot_path.resolve()),
        },
    }
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("Train-only OOF kalibrasyon analizi tamamlandı")
    print(f"OOF satırı: {len(predictions):,}; fold: {predictions['fold'].nunique()}")
    print(
        f"Brier={overall_metrics['brier_score']:.4f}; "
        f"ECE={overall_metrics['expected_calibration_error']:.4f}; "
        f"MCE={overall_metrics['maximum_calibration_error']:.4f}"
    )
    print(
        f"Kalibrasyon doğrusu: intercept="
        f"{overall_metrics['calibration_intercept']:.4f}, "
        f"slope={overall_metrics['calibration_slope']:.4f}"
    )
    print("Sabit validation ve test verisi kullanılmadı; yeniden kalibrasyon uygulanmadı.")
    print(f"Rapor: {manifest_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train-only gene-disjoint OOF olasılık kalibrasyonunu analiz eder."
    )
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--n-bins", type=int, default=DEFAULT_N_BINS)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    analyze_calibration(args.predictions, args.report_dir, n_bins=args.n_bins)


if __name__ == "__main__":
    main()
