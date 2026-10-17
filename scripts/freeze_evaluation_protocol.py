"""Test açılmadan önce model, eşik ve başarı kapılarını dondur."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

try:
    from scripts.analyze_oof_calibration import sha256_file
except ModuleNotFoundError:  # Doğrudan `python scripts/...` çalıştırması için.
    from analyze_oof_calibration import sha256_file  # type: ignore[no-redef]


DEFAULT_CONFIG = Path("config/final_evaluation_protocol.json")
DEFAULT_MODEL_MANIFEST = Path(
    "reports/generated/logistic_hgvs_gene_disjoint/manifest.json"
)
DEFAULT_THRESHOLD_MANIFEST = Path(
    "reports/generated/logistic_hgvs_gene_disjoint/threshold_analysis/manifest.json"
)
DEFAULT_OPERATING_POINTS = Path(
    "reports/generated/logistic_hgvs_gene_disjoint/threshold_analysis/operating_points.csv"
)
DEFAULT_CV_MANIFEST = Path("reports/generated/hgvs_gene_disjoint_cv/manifest.json")
DEFAULT_CALIBRATION_MANIFEST = Path(
    "reports/generated/hgvs_gene_disjoint_cv/calibration/manifest.json"
)
DEFAULT_UNCERTAINTY_MANIFEST = Path(
    "reports/generated/hgvs_gene_disjoint_cv/uncertainty/manifest.json"
)
DEFAULT_OUTPUT_DIR = Path("reports/generated/final_protocol")


def load_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise FileNotFoundError(f"Gerekli JSON bulunamadı: {path}")
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"JSON nesne olmalıdır: {path}")
    return value


def require_false(value: object, label: str) -> None:
    if value is not False:
        raise ValueError(f"Test kilidi doğrulanamadı ({label} false olmalı).")


def nested(mapping: dict[str, object], *keys: str) -> object:
    value: object = mapping
    for key in keys:
        if not isinstance(value, dict) or key not in value:
            raise ValueError("Eksik manifest alanı: " + ".".join(keys))
        value = value[key]
    return value


def validate_config(config: dict[str, object]) -> None:
    if config.get("status") != "pre_test_frozen":
        raise ValueError("Protokol durumu pre_test_frozen olmalıdır.")
    if config.get("clinical_diagnosis_authorized") is not False:
        raise ValueError("Bu aşamada klinik tanı yetkisi false olmalıdır.")
    acceptance = config.get("test_acceptance")
    if not isinstance(acceptance, dict):
        raise ValueError("test_acceptance nesnesi eksik.")
    expected_gates = {
        "minimum_roc_auc",
        "minimum_average_precision",
        "minimum_sensitivity_at_frozen_threshold",
        "minimum_specificity_at_frozen_threshold",
        "maximum_roc_auc_drop_from_validation",
    }
    if set(acceptance) != expected_gates:
        raise ValueError("Test başarı kapıları eksik veya beklenmeyen alan içeriyor.")
    for name, value in acceptance.items():
        if not isinstance(value, (int, float)) or not 0 <= value <= 1:
            raise ValueError(f"Geçersiz test başarı kapısı: {name}")
    test_policy = config.get("test_policy")
    if not isinstance(test_policy, dict) or test_policy.get("maximum_evaluations") != 1:
        raise ValueError("Kilitli test en fazla bir kez değerlendirilebilir.")


def select_operating_point(
    operating_points_path: Path,
    *,
    model_name: str,
    policy: str,
    sensitivity_target: float,
) -> dict[str, object]:
    if not operating_points_path.is_file():
        raise FileNotFoundError(
            f"Çalışma noktaları dosyası bulunamadı: {operating_points_path}"
        )
    points = pd.read_csv(operating_points_path)
    required = {
        "model_name",
        "policy",
        "sensitivity_target",
        "threshold",
        "sensitivity",
        "specificity",
        "precision_positive_predictive_value",
        "negative_predictive_value",
        "false_positives",
        "false_negatives",
    }
    missing = sorted(required.difference(points.columns))
    if missing:
        raise ValueError("Çalışma noktalarında eksik sütunlar: " + ", ".join(missing))
    selected = points.loc[
        (points["model_name"] == model_name)
        & (points["policy"] == policy)
        & np_isclose(points["sensitivity_target"], sensitivity_target)
    ]
    if len(selected) != 1:
        raise ValueError(
            "Seçilen model/politika/duyarlılık hedefi için tam bir çalışma noktası "
            "bulunmalıdır."
        )
    row = selected.iloc[0]
    if not 0 < float(row["threshold"]) < 1:
        raise ValueError("Donmuş eşik 0 ile 1 arasında olmalıdır.")
    return {
        "model_name": model_name,
        "policy": policy,
        "sensitivity_target": sensitivity_target,
        "threshold": float(row["threshold"]),
        "validation_sensitivity": float(row["sensitivity"]),
        "validation_specificity": float(row["specificity"]),
        "validation_precision": float(row["precision_positive_predictive_value"]),
        "validation_negative_predictive_value": float(row["negative_predictive_value"]),
        "validation_false_positives": int(row["false_positives"]),
        "validation_false_negatives": int(row["false_negatives"]),
    }


def np_isclose(series: pd.Series, value: float) -> pd.Series:
    return (series.astype(float) - value).abs() <= 1e-12


def format_protocol_card(lock: dict[str, object]) -> str:
    operating = lock["frozen_operating_point"]
    acceptance = lock["test_acceptance"]
    validation = lock["validation_reference"]
    assert isinstance(operating, dict)
    assert isinstance(acceptance, dict)
    assert isinstance(validation, dict)
    return f"""# Nihai değerlendirme protokolü

Durum: **TEST AÇILMADI — PROTOKOL DONDURULDU**

## Amaç

Bu model yalnızca araştırma amaçlı varyant önceliklendirme prototipidir. Klinik
tanı koymaz ve uzman değerlendirmesinin yerine geçmez.

## Dondurulan seçimler

- Model: `{operating['model_name']}`
- Özellik profili: `{lock['selected_feature_profile']}`
- Artifact SHA-256: `{lock['artifact']['sha256']}`
- Karar eşiği: `{operating['threshold']:.10f}`
- Eşik politikası: validation üzerinde en az `%{100 * operating['sensitivity_target']:.0f}` duyarlılık
- Validation duyarlılığı: `{operating['validation_sensitivity']:.4f}`
- Validation özgüllüğü: `{operating['validation_specificity']:.4f}`
- Validation ROC-AUC: `{validation['roc_auc']:.4f}`
- Validation average precision: `{validation['average_precision']:.4f}`

## Test başarı kapıları

Tüm koşullar birlikte sağlanmalıdır:

- ROC-AUC ≥ `{acceptance['minimum_roc_auc']:.2f}`
- Average precision ≥ `{acceptance['minimum_average_precision']:.2f}`
- Donmuş eşikte duyarlılık ≥ `{acceptance['minimum_sensitivity_at_frozen_threshold']:.2f}`
- Donmuş eşikte özgüllük ≥ `{acceptance['minimum_specificity_at_frozen_threshold']:.2f}`
- Validation'a göre ROC-AUC düşüşü ≤ `{acceptance['maximum_roc_auc_drop_from_validation']:.2f}`

## Test açma kuralı

Kilitli test yalnızca bir kez değerlendirilir. Test sonucu görüldükten sonra
model, özellikler veya eşik değiştirilmeden birincil sonuç raporlanır. Başarı
kapıları geçilmezse sonraki iddialar için yeni ve bağımsız bir holdout gerekir.
Bu iç test, bağımsız dış doğrulamanın yerine geçmez.

## Bütünlük

- Test etiketi okundu: **hayır**
- Geliştirme raporlarında test kullanımı: **hayır**
- Girdi tutarlılığı: **doğrulandı**
- Protokol sürümü: `{lock['protocol_version']}`
"""


def freeze_protocol(
    config_path: Path,
    model_manifest_path: Path,
    threshold_manifest_path: Path,
    operating_points_path: Path,
    cv_manifest_path: Path,
    calibration_manifest_path: Path,
    uncertainty_manifest_path: Path,
    output_dir: Path,
) -> dict[str, object]:
    config = load_json(config_path)
    validate_config(config)
    model_manifest = load_json(model_manifest_path)
    threshold_manifest = load_json(threshold_manifest_path)
    cv_manifest = load_json(cv_manifest_path)
    calibration_manifest = load_json(calibration_manifest_path)
    uncertainty_manifest = load_json(uncertainty_manifest_path)

    if nested(model_manifest, "data_policy", "test_partition") != (
        "locked_not_loaded_into_development_dataframe"
    ):
        raise ValueError("Model manifestinde test kilidi doğrulanamadı.")
    require_false(
        nested(threshold_manifest, "data_policy", "test_data_used"),
        "threshold_analysis.test_data_used",
    )
    require_false(
        nested(cv_manifest, "data_policy", "test_used"), "cross_validation.test_used"
    )
    require_false(
        nested(calibration_manifest, "data_policy", "test_used"),
        "calibration.test_used",
    )
    require_false(
        nested(uncertainty_manifest, "data_policy", "test_used"),
        "uncertainty.test_used",
    )

    feature_hash = nested(model_manifest, "inputs", "features", "sha256")
    split_hash = nested(model_manifest, "inputs", "split", "sha256")
    if feature_hash != nested(cv_manifest, "inputs", "features", "sha256"):
        raise ValueError("Model ve çapraz doğrulama özellik hash'leri eşleşmiyor.")
    if feature_hash != nested(uncertainty_manifest, "inputs", "features", "sha256"):
        raise ValueError("Model ve belirsizlik analizi özellik hash'leri eşleşmiyor.")
    if split_hash != nested(cv_manifest, "inputs", "gene_disjoint_split", "sha256"):
        raise ValueError("Model ve çapraz doğrulama ayrım hash'leri eşleşmiyor.")

    model_name = str(config["selected_model"])
    feature_profile = str(config["selected_feature_profile"])
    if nested(model_manifest, "data_policy", "feature_profile") != feature_profile:
        raise ValueError("Seçilen özellik profili model manifestiyle eşleşmiyor.")
    models = model_manifest.get("models")
    if not isinstance(models, dict) or model_name not in models:
        raise ValueError(f"Seçilen model model manifestinde yok: {model_name}")
    model_info = models[model_name]
    if not isinstance(model_info, dict):
        raise ValueError("Seçilen model bilgisi geçersiz.")
    artifact_path = Path(str(model_info["artifact_path"]))
    if not artifact_path.is_file():
        raise FileNotFoundError(f"Donmuş model artifact'i bulunamadı: {artifact_path}")
    artifact_hash = sha256_file(artifact_path)
    if artifact_hash != model_info.get("artifact_sha256"):
        raise ValueError("Model artifact hash'i manifest ile eşleşmiyor.")

    threshold_source = Path(str(nested(threshold_manifest, "source", "path")))
    if not threshold_source.is_file():
        raise FileNotFoundError(f"Validation tahmin kaynağı bulunamadı: {threshold_source}")
    if sha256_file(threshold_source) != nested(
        threshold_manifest, "source", "sha256"
    ):
        raise ValueError("Validation tahmin hash'i eşleşmiyor.")

    operating_config = config.get("operating_point")
    if not isinstance(operating_config, dict):
        raise ValueError("operating_point nesnesi eksik.")
    operating = select_operating_point(
        operating_points_path,
        model_name=model_name,
        policy=str(operating_config["policy"]),
        sensitivity_target=float(operating_config["sensitivity_target"]),
    )
    validation_metrics = nested(model_manifest, "validation_metrics", model_name)
    if not isinstance(validation_metrics, dict):
        raise ValueError("Validation metrikleri geçersiz.")

    source_paths = {
        "config": config_path,
        "model_manifest": model_manifest_path,
        "threshold_manifest": threshold_manifest_path,
        "operating_points": operating_points_path,
        "cross_validation_manifest": cv_manifest_path,
        "calibration_manifest": calibration_manifest_path,
        "uncertainty_manifest": uncertainty_manifest_path,
    }
    lock: dict[str, object] = {
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_version": config["protocol_version"],
        "status": "pre_test_frozen",
        "test_evaluated": False,
        "test_labels_read": False,
        "intended_use": config["intended_use"],
        "clinical_diagnosis_authorized": False,
        "selected_feature_profile": feature_profile,
        "frozen_operating_point": operating,
        "validation_reference": {
            "roc_auc": float(validation_metrics["roc_auc"]),
            "average_precision": float(validation_metrics["average_precision"]),
        },
        "cross_validation_reference": nested(cv_manifest, "out_of_fold_metrics"),
        "test_acceptance": config["test_acceptance"],
        "test_policy": config["test_policy"],
        "artifact": {
            "path": str(artifact_path.resolve()),
            "sha256": artifact_hash,
            "feature_columns": model_info["feature_columns"],
        },
        "data_fingerprints": {
            "features_sha256": feature_hash,
            "gene_disjoint_split_sha256": split_hash,
            "validation_predictions_sha256": nested(
                threshold_manifest, "source", "sha256"
            ),
        },
        "source_files": {
            name: {"path": str(path.resolve()), "sha256": sha256_file(path)}
            for name, path in source_paths.items()
        },
        "limitations": config["limitations"],
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    lock_path = output_dir / "protocol_lock.json"
    card_path = output_dir / "protocol_card.md"
    with lock_path.open("w", encoding="utf-8") as handle:
        json.dump(lock, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    card_path.write_text(format_protocol_card(lock), encoding="utf-8")

    print("Nihai değerlendirme protokolü donduruldu")
    print(f"Model: {model_name}")
    print(f"Donmuş eşik: {operating['threshold']:.10f}")
    print("Test etiketi okunmadı; test henüz değerlendirilmedi.")
    print(f"Protokol kilidi: {lock_path.resolve()}")
    return lock


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Test öncesi model, eşik ve başarı ölçütlerini dondurur."
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--model-manifest", type=Path, default=DEFAULT_MODEL_MANIFEST)
    parser.add_argument(
        "--threshold-manifest", type=Path, default=DEFAULT_THRESHOLD_MANIFEST
    )
    parser.add_argument("--operating-points", type=Path, default=DEFAULT_OPERATING_POINTS)
    parser.add_argument("--cv-manifest", type=Path, default=DEFAULT_CV_MANIFEST)
    parser.add_argument(
        "--calibration-manifest", type=Path, default=DEFAULT_CALIBRATION_MANIFEST
    )
    parser.add_argument(
        "--uncertainty-manifest", type=Path, default=DEFAULT_UNCERTAINTY_MANIFEST
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    freeze_protocol(
        args.config,
        args.model_manifest,
        args.threshold_manifest,
        args.operating_points,
        args.cv_manifest,
        args.calibration_manifest,
        args.uncertainty_manifest,
        args.output_dir,
    )


if __name__ == "__main__":
    main()
