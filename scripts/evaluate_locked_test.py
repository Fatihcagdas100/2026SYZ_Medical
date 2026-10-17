"""Dondurulmuş modeli kilitli testte yalnızca bir kez değerlendir."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import joblib
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import ConfusionMatrixDisplay, precision_recall_curve, roc_curve

try:
    from scripts.analyze_oof_calibration import sha256_file
    from scripts.train_logistic_baseline import calculate_metrics
except ModuleNotFoundError:  # Doğrudan `python scripts/...` çalıştırması için.
    from analyze_oof_calibration import sha256_file  # type: ignore[no-redef]
    from train_logistic_baseline import calculate_metrics  # type: ignore[no-redef]


DEFAULT_PROTOCOL_LOCK = Path("reports/generated/final_protocol/protocol_lock.json")
DEFAULT_FEATURES = Path("data/processed/clinvar_hgvs_features.tsv.gz")
DEFAULT_SPLIT = Path("data/processed/baseline_gene_disjoint_split.tsv")
DEFAULT_OUTPUT_DIR = Path("reports/generated/final_test")
DEFAULT_CHUNK_SIZE = 100_000


def load_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(f"Gerekli JSON bulunamadı: {path}")
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"JSON nesne olmalıdır: {path}")
    return value


def nested(mapping: dict[str, object], *keys: str) -> object:
    value: object = mapping
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            raise ValueError("Eksik protokol alanı: " + ".".join(keys))
        value = value[key]
    return value


def verify_pretest_lock(
    lock: dict[str, object],
    *,
    features_path: Path,
    split_path: Path,
) -> tuple[Path, list[str], float]:
    if lock.get("status") != "pre_test_frozen":
        raise ValueError("Protokol pre_test_frozen durumunda değil.")
    if lock.get("test_evaluated") is not False or lock.get("test_labels_read") is not False:
        raise ValueError("Protokol kilidi testin daha önce açılmadığını doğrulamıyor.")
    if nested(lock, "test_policy", "maximum_evaluations") != 1:
        raise ValueError("Protokol yalnızca tek test değerlendirmesine izin vermelidir.")
    if sha256_file(features_path) != nested(
        lock, "data_fingerprints", "features_sha256"
    ):
        raise ValueError("Özellik dosyası hash'i donmuş protokolle eşleşmiyor.")
    if sha256_file(split_path) != nested(
        lock, "data_fingerprints", "gene_disjoint_split_sha256"
    ):
        raise ValueError("Veri ayrımı hash'i donmuş protokolle eşleşmiyor.")

    artifact_path = Path(str(nested(lock, "artifact", "path")))
    if not artifact_path.is_file():
        raise FileNotFoundError(f"Donmuş model artifact'i bulunamadı: {artifact_path}")
    if sha256_file(artifact_path) != nested(lock, "artifact", "sha256"):
        raise ValueError("Model artifact hash'i donmuş protokolle eşleşmiyor.")
    feature_columns = nested(lock, "artifact", "feature_columns")
    if not isinstance(feature_columns, list) or not feature_columns:
        raise ValueError("Donmuş model özellik listesi geçersiz.")
    threshold = float(nested(lock, "frozen_operating_point", "threshold"))
    if not 0 < threshold < 1:
        raise ValueError("Donmuş eşik geçersiz.")
    return artifact_path, [str(column) for column in feature_columns], threshold


def load_test_assignments(split_path: Path) -> pd.DataFrame:
    assignments = pd.read_csv(
        split_path,
        sep="\t",
        usecols=["VariationID", "split", "gene_group"],
        dtype={"VariationID": "string", "split": "string", "gene_group": "string"},
    )
    if assignments.isna().any().any():
        raise ValueError("Ayrım dosyasının gerekli alanlarında boş değer bulundu.")
    if assignments["VariationID"].duplicated().any():
        raise ValueError("Ayrım dosyasında tekrarlanan VariationID bulundu.")
    if set(assignments["split"].unique()) != {"train", "validation", "test"}:
        raise ValueError("Ayrım dosyası train, validation ve test içermelidir.")
    test_assignments = assignments.loc[
        assignments["split"] == "test", ["VariationID", "gene_group"]
    ].copy()
    if test_assignments.empty:
        raise ValueError("Kilitli test bölümünde kayıt yok.")
    return test_assignments


def write_evaluation_started_marker(
    output_dir: Path,
    *,
    protocol_lock_path: Path,
    test_rows: int,
) -> Path:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(
            "Test çıktı klasörü boş değil; tek-sefer kuralı nedeniyle yeniden "
            f"değerlendirme reddedildi: {output_dir}"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    marker_path = output_dir / "EVALUATION_STARTED.json"
    marker = {
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "evaluation_number": 1,
        "protocol_lock_path": str(protocol_lock_path.resolve()),
        "protocol_lock_sha256": sha256_file(protocol_lock_path),
        "expected_test_rows": test_rows,
        "warning": "Test labels are considered opened after this marker is created.",
    }
    with marker_path.open("w", encoding="utf-8") as handle:
        json.dump(marker, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    return marker_path


def load_test_features(
    features_path: Path,
    test_assignments: pd.DataFrame,
    *,
    feature_columns: list[str],
    categorical_columns: list[str] | None = None,
    chunk_size: int,
) -> pd.DataFrame:
    if chunk_size <= 0:
        raise ValueError("chunk_size pozitif olmalıdır.")
    test_ids = set(test_assignments["VariationID"].tolist())
    dtypes = {
        "VariationID": "string",
        "target": "Int8",
        "is_canonical_snv": "Int8",
        "has_protein_annotation": "Int8",
    }
    dtypes.update({column: "string" for column in categorical_columns or []})
    retained: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        features_path,
        sep="\t",
        usecols=["VariationID", *feature_columns, "target"],
        dtype=dtypes,
        chunksize=chunk_size,
    ):
        selected = chunk["VariationID"].isin(test_ids)
        if selected.any():
            retained.append(chunk.loc[selected].copy())
    if not retained:
        raise ValueError("Test kimliklerinin hiçbiri özellik dosyasında bulunamadı.")
    test_data = pd.concat(retained, ignore_index=True)
    if test_data["VariationID"].duplicated().any():
        raise ValueError("Test özelliklerinde tekrarlanan VariationID bulundu.")
    found_ids = set(test_data["VariationID"].tolist())
    if found_ids != test_ids:
        missing = test_ids.difference(found_ids)
        raise ValueError(f"Özellik dosyasında bulunmayan test VariationID: {min(missing)}")
    if test_data["target"].isna().any() or set(
        test_data["target"].astype(int).unique()
    ) != {0, 1}:
        raise ValueError("Kilitli test iki hedef sınıfını da içermelidir.")
    return test_data.merge(
        test_assignments, on="VariationID", how="inner", validate="one_to_one"
    )


def get_model_categorical_columns(model: object) -> list[str]:
    named_steps = getattr(model, "named_steps", None)
    if named_steps is None or "preprocessing" not in named_steps:
        return []
    preprocessing = named_steps["preprocessing"]
    for name, _transformer, columns in preprocessing.transformers_:
        if name == "categorical":
            return [str(column) for column in columns]
    return []


def evaluate_acceptance(
    metrics: dict[str, float | int],
    lock: dict[str, object],
) -> tuple[dict[str, dict[str, object]], float, bool]:
    gates = nested(lock, "test_acceptance")
    if not isinstance(gates, dict):
        raise ValueError("Test başarı kapıları geçersiz.")
    validation_roc = float(nested(lock, "validation_reference", "roc_auc"))
    roc_drop = validation_roc - float(metrics["roc_auc"])
    checks: dict[str, dict[str, object]] = {
        "minimum_roc_auc": {
            "observed": metrics["roc_auc"],
            "required": gates["minimum_roc_auc"],
            "passed": float(metrics["roc_auc"]) >= float(gates["minimum_roc_auc"]),
        },
        "minimum_average_precision": {
            "observed": metrics["average_precision"],
            "required": gates["minimum_average_precision"],
            "passed": float(metrics["average_precision"])
            >= float(gates["minimum_average_precision"]),
        },
        "minimum_sensitivity_at_frozen_threshold": {
            "observed": metrics["recall_sensitivity"],
            "required": gates["minimum_sensitivity_at_frozen_threshold"],
            "passed": float(metrics["recall_sensitivity"])
            >= float(gates["minimum_sensitivity_at_frozen_threshold"]),
        },
        "minimum_specificity_at_frozen_threshold": {
            "observed": metrics["specificity"],
            "required": gates["minimum_specificity_at_frozen_threshold"],
            "passed": float(metrics["specificity"])
            >= float(gates["minimum_specificity_at_frozen_threshold"]),
        },
        "maximum_roc_auc_drop_from_validation": {
            "observed": roc_drop,
            "required_maximum": gates["maximum_roc_auc_drop_from_validation"],
            "passed": roc_drop
            <= float(gates["maximum_roc_auc_drop_from_validation"]),
        },
    }
    return checks, roc_drop, all(bool(check["passed"]) for check in checks.values())


def plot_test_results(
    targets: np.ndarray,
    probabilities: np.ndarray,
    predictions: np.ndarray,
    metrics: dict[str, float | int],
    output_path: Path,
) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    false_positive_rate, true_positive_rate, _ = roc_curve(targets, probabilities)
    axes[0].plot(false_positive_rate, true_positive_rate, color="tab:blue")
    axes[0].plot([0, 1], [0, 1], linestyle="--", color="gray")
    axes[0].set_title(f"Kilitli test ROC (AUC={metrics['roc_auc']:.4f})")
    axes[0].set_xlabel("1 − Özgüllük")
    axes[0].set_ylabel("Duyarlılık")

    precision, recall, _ = precision_recall_curve(targets, probabilities)
    axes[1].plot(recall, precision, color="tab:orange")
    axes[1].axhline(targets.mean(), linestyle="--", color="gray")
    axes[1].set_title(f"Precision-recall (AP={metrics['average_precision']:.4f})")
    axes[1].set_xlabel("Duyarlılık")
    axes[1].set_ylabel("Kesinlik")

    ConfusionMatrixDisplay.from_predictions(
        targets,
        predictions,
        labels=[0, 1],
        display_labels=["Benign", "Pathogenic"],
        colorbar=False,
        ax=axes[2],
    )
    axes[2].set_title(f"Donmuş eşik = {metrics['threshold']:.4f}")
    for axis in axes[:2]:
        axis.set_xlim(0, 1)
        axis.set_ylim(0, 1)
        axis.grid(alpha=0.2)
    figure.suptitle("Tek seferlik kilitli test değerlendirmesi", fontsize=14)
    figure.tight_layout()
    figure.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(figure)


def format_test_report(manifest: dict[str, object]) -> str:
    metrics = manifest["test_metrics"]
    checks = manifest["acceptance_checks"]
    assert isinstance(metrics, dict)
    assert isinstance(checks, dict)
    lines = [
        "# Kilitli test değerlendirmesi",
        "",
        f"Birincil sonuç: **{'GEÇTİ' if manifest['all_acceptance_gates_passed'] else 'GEÇMEDİ'}**",
        "",
        "Model, özellikler ve eşik test açılmadan önce dondurulmuştur. Bu rapor",
        "tek seferlik birincil test sonucudur.",
        "",
        "## Test metrikleri",
        "",
        "| Metrik | Değer |",
        "| --- | ---: |",
        f"| ROC-AUC | {metrics['roc_auc']:.4f} |",
        f"| Average precision | {metrics['average_precision']:.4f} |",
        f"| Duyarlılık | {metrics['recall_sensitivity']:.4f} |",
        f"| Özgüllük | {metrics['specificity']:.4f} |",
        f"| F1 | {metrics['f1']:.4f} |",
        f"| Dengeli doğruluk | {metrics['balanced_accuracy']:.4f} |",
        f"| Brier skoru | {metrics['brier_score']:.4f} |",
        f"| Yanlış negatif | {metrics['false_negatives']:,} |",
        f"| Yanlış pozitif | {metrics['false_positives']:,} |",
        "",
        "## Önceden belirlenmiş başarı kapıları",
        "",
        "| Kapı | Gözlenen | Sonuç |",
        "| --- | ---: | --- |",
    ]
    for name, check in checks.items():
        assert isinstance(check, dict)
        lines.append(
            f"| {name} | {float(check['observed']):.4f} | "
            f"{'GEÇTİ' if check['passed'] else 'GEÇMEDİ'} |"
        )
    lines.extend(
        [
            "",
            "Bu iç test bağımsız dış doğrulamanın veya klinik değerlendirmenin yerine geçmez.",
            "",
        ]
    )
    return "\n".join(lines)


def evaluate_locked_test(
    protocol_lock_path: Path,
    features_path: Path,
    split_path: Path,
    output_dir: Path,
    *,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, object]:
    lock = load_json(protocol_lock_path)
    artifact_path, feature_columns, threshold = verify_pretest_lock(
        lock, features_path=features_path, split_path=split_path
    )
    test_assignments = load_test_assignments(split_path)
    model = joblib.load(artifact_path)

    marker_path = write_evaluation_started_marker(
        output_dir,
        protocol_lock_path=protocol_lock_path,
        test_rows=len(test_assignments),
    )
    test_data = load_test_features(
        features_path,
        test_assignments,
        feature_columns=feature_columns,
        categorical_columns=get_model_categorical_columns(model),
        chunk_size=chunk_size,
    )
    probabilities = model.predict_proba(test_data.loc[:, feature_columns])[:, 1]
    targets = test_data["target"].astype(int).reset_index(drop=True)
    probability_series = pd.Series(probabilities)
    metrics = calculate_metrics(targets, probability_series, threshold=threshold)
    predictions = (probabilities >= threshold).astype(int)
    checks, roc_drop, passed = evaluate_acceptance(metrics, lock)

    predictions_table = pd.DataFrame(
        {
            "VariationID": test_data["VariationID"],
            "gene_group": test_data["gene_group"],
            "source_partition": "test",
            "target": targets,
            "probability_pathogenic": probabilities,
            "prediction": predictions,
        }
    )
    predictions_path = output_dir / "test_predictions.tsv.gz"
    metrics_path = output_dir / "test_metrics.json"
    plot_path = output_dir / "test_roc_pr_confusion.png"
    report_path = output_dir / "test_report.md"
    manifest_path = output_dir / "manifest.json"
    predictions_table.to_csv(
        predictions_path,
        sep="\t",
        index=False,
        compression="gzip",
        float_format="%.10f",
    )
    with metrics_path.open("w", encoding="utf-8") as handle:
        json.dump(metrics, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    plot_test_results(
        targets.to_numpy(), probabilities, predictions, metrics, plot_path
    )

    manifest: dict[str, object] = {
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "evaluation_number": 1,
        "test_evaluated": True,
        "test_labels_read": True,
        "protocol_lock": {
            "path": str(protocol_lock_path.resolve()),
            "sha256": sha256_file(protocol_lock_path),
        },
        "frozen_model": nested(lock, "frozen_operating_point", "model_name"),
        "frozen_threshold": threshold,
        "test_rows": len(test_data),
        "test_gene_groups": int(test_data["gene_group"].nunique()),
        "test_metrics": metrics,
        "validation_roc_auc": nested(lock, "validation_reference", "roc_auc"),
        "roc_auc_drop_from_validation": roc_drop,
        "acceptance_checks": checks,
        "all_acceptance_gates_passed": passed,
        "data_fingerprints": {
            "features_sha256": sha256_file(features_path),
            "gene_disjoint_split_sha256": sha256_file(split_path),
            "artifact_sha256": sha256_file(artifact_path),
        },
        "outputs": {
            "started_marker": str(marker_path.resolve()),
            "test_predictions": str(predictions_path.resolve()),
            "test_metrics": str(metrics_path.resolve()),
            "test_plot": str(plot_path.resolve()),
            "test_report": str(report_path.resolve()),
        },
        "reporting_policy": {
            "model_or_threshold_changed_after_test": False,
            "primary_result_reported_without_retuning": True,
            "external_validation_still_required": True,
        },
    }
    report_path.write_text(format_test_report(manifest), encoding="utf-8")
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("Kilitli test tek seferlik değerlendirildi")
    print(f"Test satırı: {len(test_data):,}; gen grubu: {test_data['gene_group'].nunique():,}")
    print(
        f"ROC-AUC={metrics['roc_auc']:.4f}; AP={metrics['average_precision']:.4f}; "
        f"duyarlılık={metrics['recall_sensitivity']:.4f}; "
        f"özgüllük={metrics['specificity']:.4f}"
    )
    print(f"Önceden belirlenmiş tüm başarı kapıları: {'GEÇTİ' if passed else 'GEÇMEDİ'}")
    print("Model ve eşik değiştirilmeden birincil sonuç kaydedildi.")
    print(f"Rapor: {report_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Dondurulmuş modeli kilitli testte yalnızca bir kez değerlendirir."
    )
    parser.add_argument("--protocol-lock", type=Path, default=DEFAULT_PROTOCOL_LOCK)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--split", type=Path, default=DEFAULT_SPLIT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    evaluate_locked_test(
        args.protocol_lock,
        args.features,
        args.split,
        args.output_dir,
        chunk_size=args.chunk_size,
    )


if __name__ == "__main__":
    main()
