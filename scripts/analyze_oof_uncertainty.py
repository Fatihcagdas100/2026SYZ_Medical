"""Train-only OOF tahminlerinde belirsiz bölge ve alt grup hatalarını incele."""

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
    from scripts.analyze_oof_calibration import load_oof_predictions, sha256_file
except ModuleNotFoundError:  # Doğrudan `python scripts/...` çalıştırması için.
    from analyze_oof_calibration import (  # type: ignore[no-redef]
        load_oof_predictions,
        sha256_file,
    )


DEFAULT_PREDICTIONS = Path(
    "reports/generated/hgvs_gene_disjoint_cv/train_out_of_fold_predictions.tsv.gz"
)
DEFAULT_FEATURES = Path("data/processed/clinvar_hgvs_features.tsv.gz")
DEFAULT_REPORT_DIR = Path("reports/generated/hgvs_gene_disjoint_cv/uncertainty")
DEFAULT_THRESHOLD = 0.5
DEFAULT_LOWER_BOUND = 0.2
DEFAULT_UPPER_BOUND = 0.8
DEFAULT_MIN_GROUP_SIZE = 100
DEFAULT_CHUNK_SIZE = 100_000

FEATURE_COLUMNS = (
    "variant_type",
    "chromosome",
    "mutation_class",
    "hgvs_sequence_type",
    "dna_change_type",
    "coding_region",
    "splice_offset_category",
    "protein_effect",
    "has_protein_annotation",
    "is_canonical_snv",
    "absolute_length_delta",
)
GROUP_FEATURES = (
    "variant_type",
    "mutation_class",
    "hgvs_sequence_type",
    "dna_change_type",
    "coding_region",
    "splice_offset_category",
    "protein_effect",
    "has_protein_annotation",
    "is_canonical_snv",
    "absolute_length_delta_bucket",
)
CONFIDENT_ERROR_GROUP_FEATURES = ("gene_group", *GROUP_FEATURES)
REGION_ORDER = ("low_probability", "uncertain", "high_probability")


def load_oof_features(
    features_path: Path,
    variation_ids: set[str],
    *,
    chunk_size: int,
) -> pd.DataFrame:
    if not features_path.is_file():
        raise FileNotFoundError(f"Özellik dosyası bulunamadı: {features_path}")
    if chunk_size <= 0:
        raise ValueError("chunk_size pozitif olmalıdır.")

    dtypes: dict[str, str] = {
        "VariationID": "string",
        "target": "Int8",
        "has_protein_annotation": "Int8",
        "is_canonical_snv": "Int8",
    }
    dtypes.update(
        {
            feature: "string"
            for feature in FEATURE_COLUMNS
            if feature
            not in {
                "has_protein_annotation",
                "is_canonical_snv",
                "absolute_length_delta",
            }
        }
    )
    retained: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        features_path,
        sep="\t",
        usecols=["VariationID", *FEATURE_COLUMNS, "target"],
        dtype=dtypes,
        chunksize=chunk_size,
    ):
        selected = chunk["VariationID"].isin(variation_ids)
        if selected.any():
            retained.append(chunk.loc[selected].copy())
    if not retained:
        raise ValueError("OOF kimliklerinin hiçbiri özellik dosyasında bulunamadı.")

    features = pd.concat(retained, ignore_index=True)
    if features["VariationID"].duplicated().any():
        duplicate = features.loc[
            features["VariationID"].duplicated(), "VariationID"
        ].iloc[0]
        raise ValueError(f"Özellik dosyasında tekrarlanan VariationID: {duplicate}")
    found_ids = set(features["VariationID"].tolist())
    if found_ids != variation_ids:
        missing = variation_ids.difference(found_ids)
        raise ValueError(f"Özellik dosyasında bulunmayan OOF VariationID: {min(missing)}")
    return features


def prepare_analysis_data(
    predictions: pd.DataFrame,
    features: pd.DataFrame,
    *,
    threshold: float,
    lower_bound: float,
    upper_bound: float,
) -> pd.DataFrame:
    features = features.rename(columns={"target": "feature_target"})
    data = predictions.merge(
        features,
        on="VariationID",
        how="inner",
        validate="one_to_one",
    )
    if not data["target"].astype(int).equals(data["feature_target"].astype(int)):
        raise ValueError("OOF tahminleri ile özellik tablosunun hedefleri eşleşmiyor.")
    data = data.drop(columns="feature_target")
    data["prediction"] = (
        data["probability_pathogenic"] >= threshold
    ).astype("int8")
    data["error_type"] = "correct"
    data.loc[
        (data["target"] == 1) & (data["prediction"] == 0), "error_type"
    ] = "false_negative"
    data.loc[
        (data["target"] == 0) & (data["prediction"] == 1), "error_type"
    ] = "false_positive"
    data["probability_region"] = "uncertain"
    data.loc[
        data["probability_pathogenic"] < lower_bound, "probability_region"
    ] = "low_probability"
    data.loc[
        data["probability_pathogenic"] > upper_bound, "probability_region"
    ] = "high_probability"
    data["absolute_length_delta_bucket"] = pd.cut(
        data["absolute_length_delta"],
        bins=[-np.inf, 0, 1, 5, 10, 50, 100, np.inf],
        labels=["0", "1", "2-5", "6-10", "11-50", "51-100", "101+"],
    ).astype("string").fillna("missing")
    return data


def build_region_summary(data: pd.DataFrame) -> pd.DataFrame:
    total_errors = int((data["error_type"] != "correct").sum())
    total_false_negatives = int((data["error_type"] == "false_negative").sum())
    total_false_positives = int((data["error_type"] == "false_positive").sum())
    rows = []
    for region in REGION_ORDER:
        group = data.loc[data["probability_region"] == region]
        false_negatives = int((group["error_type"] == "false_negative").sum())
        false_positives = int((group["error_type"] == "false_positive").sum())
        error_count = false_negatives + false_positives
        rows.append(
            {
                "probability_region": region,
                "row_count": len(group),
                "row_fraction": len(group) / len(data),
                "positive_count": int((group["target"] == 1).sum()),
                "negative_count": int((group["target"] == 0).sum()),
                "mean_probability": group["probability_pathogenic"].mean(),
                "observed_pathogenic_rate": group["target"].mean(),
                "correct_count": int((group["error_type"] == "correct").sum()),
                "error_count": error_count,
                "error_rate": error_count / len(group) if len(group) else np.nan,
                "false_negatives": false_negatives,
                "false_positives": false_positives,
                "share_of_all_errors": (
                    error_count / total_errors if total_errors else np.nan
                ),
                "share_of_all_false_negatives": (
                    false_negatives / total_false_negatives
                    if total_false_negatives
                    else np.nan
                ),
                "share_of_all_false_positives": (
                    false_positives / total_false_positives
                    if total_false_positives
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def subgroup_metrics(group: pd.DataFrame, *, min_group_size: int) -> dict[str, object]:
    positive_count = int((group["target"] == 1).sum())
    negative_count = int((group["target"] == 0).sum())
    false_negatives = int((group["error_type"] == "false_negative").sum())
    false_positives = int((group["error_type"] == "false_positive").sum())
    error_count = false_negatives + false_positives
    mean_probability = float(group["probability_pathogenic"].mean())
    observed_rate = float(group["target"].mean())
    return {
        "row_count": len(group),
        "below_minimum_group_size": len(group) < min_group_size,
        "positive_count": positive_count,
        "negative_count": negative_count,
        "mean_probability": mean_probability,
        "observed_pathogenic_rate": observed_rate,
        "calibration_gap": mean_probability - observed_rate,
        "error_count": error_count,
        "error_rate": error_count / len(group),
        "false_negatives": false_negatives,
        "false_negative_rate": (
            false_negatives / positive_count if positive_count else np.nan
        ),
        "false_positives": false_positives,
        "false_positive_rate": (
            false_positives / negative_count if negative_count else np.nan
        ),
    }


def build_uncertain_subgroups(
    uncertain: pd.DataFrame, *, min_group_size: int
) -> pd.DataFrame:
    rows: list[dict[str, object]] = [
        {
            "group_feature": "overall",
            "group_value": "all_uncertain_rows",
            **subgroup_metrics(uncertain, min_group_size=min_group_size),
        }
    ]
    for feature in GROUP_FEATURES:
        for value, group in uncertain.groupby(
            feature, dropna=False, observed=True, sort=False
        ):
            rows.append(
                {
                    "group_feature": feature,
                    "group_value": "missing" if pd.isna(value) else str(value),
                    **subgroup_metrics(group, min_group_size=min_group_size),
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["group_feature", "row_count", "group_value"],
        ascending=[True, False, True],
        ignore_index=True,
    )


def build_confident_error_subgroups(errors: pd.DataFrame) -> pd.DataFrame:
    confident = errors.loc[errors["probability_region"] != "uncertain"]
    columns = (
        "error_type",
        "group_feature",
        "group_value",
        "error_count",
        "share_within_error_type",
    )
    rows: list[dict[str, object]] = []
    for error_type, error_rows in confident.groupby("error_type", sort=True):
        error_total = len(error_rows)
        for feature in CONFIDENT_ERROR_GROUP_FEATURES:
            counts = error_rows[feature].fillna("missing").astype(str).value_counts()
            for value, count in counts.items():
                rows.append(
                    {
                        "error_type": error_type,
                        "group_feature": feature,
                        "group_value": value,
                        "error_count": int(count),
                        "share_within_error_type": count / error_total,
                    }
                )
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows, columns=columns).sort_values(
        ["error_type", "group_feature", "error_count", "group_value"],
        ascending=[True, True, False, True],
        ignore_index=True,
    )


def plot_uncertainty_overview(
    region_summary: pd.DataFrame,
    subgroup_table: pd.DataFrame,
    output_path: Path,
    *,
    min_group_size: int,
) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    labels = ["Düşük (<0.20)", "Belirsiz (0.20–0.80)", "Yüksek (>0.80)"]
    positions = np.arange(len(region_summary))
    axes[0].bar(positions, region_summary["correct_count"], label="Doğru")
    axes[0].bar(
        positions,
        region_summary["error_count"],
        bottom=region_summary["correct_count"],
        label="Hatalı",
    )
    axes[0].set_xticks(positions, labels)
    axes[0].set_ylabel("OOF varyant sayısı")
    axes[0].set_title("Olasılık bölgesine göre doğru/hatalı tahminler")
    axes[0].legend()
    axes[0].grid(axis="y", alpha=0.2)

    protein_rows = subgroup_table.loc[
        (subgroup_table["group_feature"] == "protein_effect")
        & (subgroup_table["row_count"] >= min_group_size)
    ].copy()
    protein_rows = protein_rows.sort_values(
        ["error_rate", "row_count"], ascending=[False, False]
    ).head(12).sort_values("error_rate")
    axes[1].barh(
        np.arange(len(protein_rows)),
        protein_rows["error_rate"],
        color="tab:red",
        alpha=0.8,
    )
    axes[1].set_yticks(
        np.arange(len(protein_rows)),
        [
            f"{value} (n={count:,})"
            for value, count in zip(
                protein_rows["group_value"], protein_rows["row_count"]
            )
        ],
    )
    axes[1].set_xlim(0, 1)
    axes[1].set_xlabel("Hata oranı")
    axes[1].set_title("Belirsiz bölgede protein etkisine göre hata")
    axes[1].grid(axis="x", alpha=0.2)

    figure.suptitle("Train-only gene-disjoint OOF belirsizlik analizi", fontsize=14)
    figure.tight_layout()
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)


def analyze_uncertainty(
    predictions_path: Path,
    features_path: Path,
    report_dir: Path,
    *,
    threshold: float = DEFAULT_THRESHOLD,
    lower_bound: float = DEFAULT_LOWER_BOUND,
    upper_bound: float = DEFAULT_UPPER_BOUND,
    min_group_size: int = DEFAULT_MIN_GROUP_SIZE,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, object]:
    if not 0 < lower_bound < threshold < upper_bound < 1:
        raise ValueError("Sınırlar 0 < alt < eşik < üst < 1 koşulunu sağlamalıdır.")
    if min_group_size <= 0:
        raise ValueError("min_group_size pozitif olmalıdır.")

    predictions = load_oof_predictions(predictions_path)
    variation_ids = set(predictions["VariationID"].tolist())
    features = load_oof_features(
        features_path, variation_ids, chunk_size=chunk_size
    )
    data = prepare_analysis_data(
        predictions,
        features,
        threshold=threshold,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
    )
    region_summary = build_region_summary(data)
    uncertain = data.loc[data["probability_region"] == "uncertain"].copy()
    if uncertain.empty:
        raise ValueError("Seçilen sınırlar içinde belirsiz OOF tahmini yok.")
    subgroup_table = build_uncertain_subgroups(
        uncertain, min_group_size=min_group_size
    )
    errors = data.loc[data["error_type"] != "correct"].copy()
    error_columns = [
        "VariationID",
        "gene_group",
        "fold",
        "target",
        "prediction",
        "probability_pathogenic",
        "probability_region",
        "error_type",
        *FEATURE_COLUMNS,
        "absolute_length_delta_bucket",
    ]
    errors = errors.loc[:, error_columns].sort_values(
        ["probability_region", "error_type", "probability_pathogenic"]
    )
    confident_error_subgroups = build_confident_error_subgroups(errors)

    report_dir.mkdir(parents=True, exist_ok=True)
    region_path = report_dir / "probability_region_summary.csv"
    subgroup_path = report_dir / "uncertain_subgroup_metrics.csv"
    errors_path = report_dir / "misclassified_oof_variants.tsv.gz"
    confident_subgroups_path = report_dir / "confident_error_subgroups.csv"
    plot_path = report_dir / "uncertainty_overview.png"
    region_summary.to_csv(region_path, index=False, float_format="%.10f")
    subgroup_table.to_csv(subgroup_path, index=False, float_format="%.10f")
    errors.to_csv(
        errors_path,
        sep="\t",
        index=False,
        compression="gzip",
        float_format="%.10f",
    )
    confident_error_subgroups.to_csv(
        confident_subgroups_path, index=False, float_format="%.10f"
    )
    plot_uncertainty_overview(
        region_summary,
        subgroup_table,
        plot_path,
        min_group_size=min_group_size,
    )

    uncertain_summary = region_summary.loc[
        region_summary["probability_region"] == "uncertain"
    ].iloc[0]
    confident_errors = errors.loc[
        errors["probability_region"] != "uncertain"
    ]
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
            "source_partition": "train_out_of_fold",
            "fixed_validation_used": False,
            "test_used": False,
            "threshold": threshold,
            "uncertain_lower_bound_inclusive": lower_bound,
            "uncertain_upper_bound_inclusive": upper_bound,
            "minimum_group_size": min_group_size,
        },
        "integrity": {
            "oof_rows": len(data),
            "unique_variation_ids": int(data["VariationID"].nunique()),
            "gene_groups": int(data["gene_group"].nunique()),
            "folds": int(data["fold"].nunique()),
            "target_mismatches": 0,
        },
        "counts": {
            "uncertain_rows": int(uncertain_summary["row_count"]),
            "uncertain_row_fraction": float(uncertain_summary["row_fraction"]),
            "uncertain_errors": int(uncertain_summary["error_count"]),
            "uncertain_error_rate": float(uncertain_summary["error_rate"]),
            "errors_captured_by_uncertain_band_fraction": float(
                uncertain_summary["share_of_all_errors"]
            ),
            "all_errors": len(errors),
            "errors_outside_uncertain_band": len(confident_errors),
        },
        "group_features": list(GROUP_FEATURES),
        "outputs": {
            "probability_region_summary": str(region_path.resolve()),
            "uncertain_subgroup_metrics": str(subgroup_path.resolve()),
            "misclassified_oof_variants": str(errors_path.resolve()),
            "confident_error_subgroups": str(confident_subgroups_path.resolve()),
            "uncertainty_plot": str(plot_path.resolve()),
        },
    }
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("Train-only OOF belirsizlik analizi tamamlandı")
    print(
        f"Belirsiz bölge: {len(uncertain):,}/{len(data):,} "
        f"(%{100 * len(uncertain) / len(data):.2f})"
    )
    print(
        f"Belirsiz bölge hata oranı: %{100 * uncertain_summary['error_rate']:.2f}; "
        f"tüm hataların %{100 * uncertain_summary['share_of_all_errors']:.2f}'si"
    )
    print(f"Belirsiz bölge dışındaki hatalar: {len(confident_errors):,}")
    print("Sabit validation ve test verisi kullanılmadı.")
    print(f"Rapor: {manifest_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Train-only OOF tahminlerinde belirsiz bölge hatalarını inceler."
    )
    parser.add_argument("--predictions", type=Path, default=DEFAULT_PREDICTIONS)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--lower-bound", type=float, default=DEFAULT_LOWER_BOUND)
    parser.add_argument("--upper-bound", type=float, default=DEFAULT_UPPER_BOUND)
    parser.add_argument("--min-group-size", type=int, default=DEFAULT_MIN_GROUP_SIZE)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    analyze_uncertainty(
        args.predictions,
        args.features,
        args.report_dir,
        threshold=args.threshold,
        lower_bound=args.lower_bound,
        upper_bound=args.upper_bound,
        min_group_size=args.min_group_size,
        chunk_size=args.chunk_size,
    )


if __name__ == "__main__":
    main()
