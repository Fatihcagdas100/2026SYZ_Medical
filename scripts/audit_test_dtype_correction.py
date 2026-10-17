"""İlk test koşusundaki kategorik dtype uyarısını teknik olarak denetle."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

try:
    from scripts.analyze_oof_calibration import sha256_file
    from scripts.evaluate_locked_test import (
        evaluate_acceptance,
        get_model_categorical_columns,
        load_json,
        load_test_assignments,
        load_test_features,
        nested,
        verify_pretest_lock,
    )
    from scripts.train_logistic_baseline import calculate_metrics
except ModuleNotFoundError:  # Doğrudan `python scripts/...` çalıştırması için.
    from analyze_oof_calibration import sha256_file  # type: ignore[no-redef]
    from evaluate_locked_test import (  # type: ignore[no-redef]
        evaluate_acceptance,
        get_model_categorical_columns,
        load_json,
        load_test_assignments,
        load_test_features,
        nested,
        verify_pretest_lock,
    )
    from train_logistic_baseline import calculate_metrics  # type: ignore[no-redef]


DEFAULT_ORIGINAL_MANIFEST = Path("reports/generated/final_test/manifest.json")
DEFAULT_ORIGINAL_PREDICTIONS = Path(
    "reports/generated/final_test/test_predictions.tsv.gz"
)
DEFAULT_PROTOCOL_LOCK = Path("reports/generated/final_protocol/protocol_lock.json")
DEFAULT_FEATURES = Path("data/processed/clinvar_hgvs_features.tsv.gz")
DEFAULT_SPLIT = Path("data/processed/baseline_gene_disjoint_split.tsv")
DEFAULT_OUTPUT_DIR = Path("reports/generated/final_test/technical_dtype_audit")
DEFAULT_CHUNK_SIZE = 100_000


def format_audit_report(audit: dict[str, object]) -> str:
    original = audit["original_metrics"]
    corrected = audit["corrected_metrics"]
    comparison = audit["prediction_comparison"]
    assert isinstance(original, dict)
    assert isinstance(corrected, dict)
    assert isinstance(comparison, dict)
    return f"""# Kilitli test teknik dtype düzeltme denetimi

İlk koşu korunmuştur ve silinmemiştir. Koşu sırasında `chromosome` alanında
karışık veri tipi uyarısı görülmüş; eğitimde kategorik alanların string olarak
okunduğu doğrulanmıştır. Bu denetim model, özellik seti veya eşik değiştirmeden
yalnızca test okuyucusunu eğitim veri tipi sözleşmesiyle eşleştirir.

| Metrik | İlk koşu | Dtype-düzeltilmiş | Fark |
| --- | ---: | ---: | ---: |
| ROC-AUC | {original['roc_auc']:.6f} | {corrected['roc_auc']:.6f} | {corrected['roc_auc'] - original['roc_auc']:+.6f} |
| Average precision | {original['average_precision']:.6f} | {corrected['average_precision']:.6f} | {corrected['average_precision'] - original['average_precision']:+.6f} |
| Duyarlılık | {original['recall_sensitivity']:.6f} | {corrected['recall_sensitivity']:.6f} | {corrected['recall_sensitivity'] - original['recall_sensitivity']:+.6f} |
| Özgüllük | {original['specificity']:.6f} | {corrected['specificity']:.6f} | {corrected['specificity'] - original['specificity']:+.6f} |
| F1 | {original['f1']:.6f} | {corrected['f1']:.6f} | {corrected['f1'] - original['f1']:+.6f} |

- Olasılığı değişen kayıt: `{comparison['probabilities_changed']:,}`
- Sınıf kararı değişen kayıt: `{comparison['classifications_changed']:,}`
- En büyük mutlak olasılık farkı: `{comparison['maximum_absolute_probability_change']:.10f}`
- Düzeltme sonrası tüm başarı kapıları: **{'GEÇTİ' if audit['all_acceptance_gates_passed_after_correction'] else 'GEÇMEDİ'}**

Bu işlem test sonucuna göre model seçimi veya eşik ayarı yapmamıştır. İlk koşu
ve teknik düzeltme birlikte raporlanmalıdır.
"""


def audit_dtype_correction(
    original_manifest_path: Path,
    original_predictions_path: Path,
    protocol_lock_path: Path,
    features_path: Path,
    split_path: Path,
    output_dir: Path,
    *,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> dict[str, object]:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Teknik denetim çıktı klasörü boş değil: {output_dir}")
    original_manifest = load_json(original_manifest_path)
    if original_manifest.get("evaluation_number") != 1:
        raise ValueError("İlk test değerlendirme manifesti doğrulanamadı.")
    if original_manifest.get("test_evaluated") is not True:
        raise ValueError("İlk manifest testin değerlendirildiğini doğrulamıyor.")
    lock = load_json(protocol_lock_path)
    artifact_path, feature_columns, threshold = verify_pretest_lock(
        lock, features_path=features_path, split_path=split_path
    )
    if nested(original_manifest, "protocol_lock", "sha256") != sha256_file(
        protocol_lock_path
    ):
        raise ValueError("İlk koşu ile protokol kilidi hash'i eşleşmiyor.")
    if not original_predictions_path.is_file():
        raise FileNotFoundError(
            f"İlk test tahminleri bulunamadı: {original_predictions_path}"
        )

    model = joblib.load(artifact_path)
    categorical_columns = get_model_categorical_columns(model)
    if not categorical_columns:
        raise ValueError("Modelde kategorik sütun sözleşmesi bulunamadı.")
    test_assignments = load_test_assignments(split_path)
    test_data = load_test_features(
        features_path,
        test_assignments,
        feature_columns=feature_columns,
        categorical_columns=categorical_columns,
        chunk_size=chunk_size,
    )
    corrected_probabilities = model.predict_proba(
        test_data.loc[:, feature_columns]
    )[:, 1]
    targets = test_data["target"].astype(int).reset_index(drop=True)
    corrected_metrics = calculate_metrics(
        targets, pd.Series(corrected_probabilities), threshold=threshold
    )
    corrected_predictions = (corrected_probabilities >= threshold).astype(int)
    checks, roc_drop, passed = evaluate_acceptance(corrected_metrics, lock)

    original = pd.read_csv(
        original_predictions_path,
        sep="\t",
        dtype={"VariationID": "string"},
    )
    corrected = pd.DataFrame(
        {
            "VariationID": test_data["VariationID"],
            "target": targets,
            "probability_pathogenic": corrected_probabilities,
            "prediction": corrected_predictions,
        }
    )
    comparison = original.merge(
        corrected,
        on="VariationID",
        how="inner",
        validate="one_to_one",
        suffixes=("_original", "_corrected"),
    )
    if len(comparison) != len(test_data):
        raise ValueError("İlk ve düzeltilmiş tahmin kimlikleri eşleşmiyor.")
    if not comparison["target_original"].equals(comparison["target_corrected"]):
        raise ValueError("İlk ve düzeltilmiş test hedefleri eşleşmiyor.")
    probability_change = (
        comparison["probability_pathogenic_corrected"]
        - comparison["probability_pathogenic_original"]
    ).abs()
    prediction_comparison = {
        "rows": len(comparison),
        "probabilities_changed": int((probability_change > 1e-9).sum()),
        "classifications_changed": int(
            (
                comparison["prediction_corrected"]
                != comparison["prediction_original"]
            ).sum()
        ),
        "mean_absolute_probability_change": float(probability_change.mean()),
        "maximum_absolute_probability_change": float(probability_change.max()),
    }

    output_dir.mkdir(parents=True, exist_ok=False)
    corrected_predictions_path = output_dir / "corrected_test_predictions.tsv.gz"
    audit_path = output_dir / "audit.json"
    report_path = output_dir / "audit_report.md"
    corrected.to_csv(
        corrected_predictions_path,
        sep="\t",
        index=False,
        compression="gzip",
        float_format="%.10f",
    )
    audit: dict[str, object] = {
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "audit_reason": (
            "The first evaluator did not force training-time string dtypes for "
            "categorical columns; pandas emitted a mixed-type warning for chromosome."
        ),
        "original_evaluation_preserved": True,
        "model_changed": False,
        "threshold_changed": False,
        "feature_set_changed": False,
        "protocol_lock_sha256": sha256_file(protocol_lock_path),
        "original_manifest_sha256": sha256_file(original_manifest_path),
        "original_predictions_sha256": sha256_file(original_predictions_path),
        "categorical_columns_forced_to_string": categorical_columns,
        "original_metrics": original_manifest["test_metrics"],
        "corrected_metrics": corrected_metrics,
        "roc_auc_drop_from_validation_after_correction": roc_drop,
        "acceptance_checks_after_correction": checks,
        "all_acceptance_gates_passed_after_correction": passed,
        "prediction_comparison": prediction_comparison,
        "outputs": {
            "corrected_predictions": str(corrected_predictions_path.resolve()),
            "audit_report": str(report_path.resolve()),
        },
    }
    with audit_path.open("w", encoding="utf-8") as handle:
        json.dump(audit, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    report_path.write_text(format_audit_report(audit), encoding="utf-8")

    print("Kilitli test dtype teknik düzeltme denetimi tamamlandı")
    print(
        f"Değişen olasılık: {prediction_comparison['probabilities_changed']:,}; "
        f"değişen sınıf: {prediction_comparison['classifications_changed']:,}"
    )
    print(
        f"Düzeltilmiş ROC-AUC={corrected_metrics['roc_auc']:.4f}; "
        f"AP={corrected_metrics['average_precision']:.4f}; "
        f"duyarlılık={corrected_metrics['recall_sensitivity']:.4f}; "
        f"özgüllük={corrected_metrics['specificity']:.4f}"
    )
    print(f"Başarı kapıları: {'GEÇTİ' if passed else 'GEÇMEDİ'}")
    print(f"Denetim raporu: {report_path.resolve()}")
    return audit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="İlk kilitli test koşusundaki kategorik dtype uyarısını denetler."
    )
    parser.add_argument(
        "--original-manifest", type=Path, default=DEFAULT_ORIGINAL_MANIFEST
    )
    parser.add_argument(
        "--original-predictions", type=Path, default=DEFAULT_ORIGINAL_PREDICTIONS
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
    audit_dtype_correction(
        args.original_manifest,
        args.original_predictions,
        args.protocol_lock,
        args.features,
        args.split,
        args.output_dir,
        chunk_size=args.chunk_size,
    )


if __name__ == "__main__":
    main()
