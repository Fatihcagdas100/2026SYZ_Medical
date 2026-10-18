"""Dondurulmuş dış aday modeli kabul edilmiş kohortta tek seferde değerlendir."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

try:
    from scripts.train_logistic_baseline import calculate_metrics, sha256_file
except ModuleNotFoundError:  # Doğrudan `python scripts/...` çalıştırması için.
    from train_logistic_baseline import calculate_metrics, sha256_file  # type: ignore[no-redef]


DEFAULT_CANDIDATE = Path("reports/generated/external_candidate/manifest.json")
DEFAULT_FEATURES = Path("data/external/admitted_missense_features.tsv.gz")
DEFAULT_REPORT_DIR = Path("reports/generated/external_evaluation")
DEFAULT_BOOTSTRAP_REPLICATES = 1_000
DEFAULT_SEED = 2027


def load_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("Manifest JSON nesnesi olmalıdır.")
    return value


def bootstrap_intervals(targets: np.ndarray, probabilities: np.ndarray, *, seed: int, replicates: int) -> dict[str, dict[str, float]]:
    if replicates < 100:
        raise ValueError("Bootstrap tekrar sayısı en az 100 olmalıdır.")
    rng = np.random.default_rng(seed)
    negative = np.flatnonzero(targets == 0)
    positive = np.flatnonzero(targets == 1)
    values: dict[str, list[float]] = {"roc_auc": [], "average_precision": []}
    for _ in range(replicates):
        indices = np.concatenate((rng.choice(negative, len(negative), replace=True), rng.choice(positive, len(positive), replace=True)))
        sample_target = targets[indices]
        sample_probability = probabilities[indices]
        values["roc_auc"].append(float(roc_auc_score(sample_target, sample_probability)))
        values["average_precision"].append(float(average_precision_score(sample_target, sample_probability)))
    return {
        metric: {"lower_95": float(np.quantile(scores, 0.025)), "upper_95": float(np.quantile(scores, 0.975))}
        for metric, scores in values.items()
    }


def evaluate_external_candidate(
    candidate_path: Path,
    features_path: Path,
    report_dir: Path,
    *,
    bootstrap_replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    seed: int = DEFAULT_SEED,
) -> dict[str, object]:
    if report_dir.exists() and any(report_dir.iterdir()):
        raise FileExistsError("Dış değerlendirme rapor klasörü boş olmalıdır; tek-sefer kuralı geçerlidir.")
    candidate = load_json(candidate_path)
    if candidate.get("status") != "frozen_for_external_validation_not_evaluated":
        raise ValueError("Aday model dış değerlendirme için dondurulmuş durumda değil.")
    model_info = candidate.get("model")
    if not isinstance(model_info, dict):
        raise ValueError("Aday model bilgisi eksik.")
    artifact_path = Path(str(model_info["artifact_path"]))
    if not artifact_path.is_file() or sha256_file(artifact_path) != model_info["artifact_sha256"]:
        raise ValueError("Dondurulmuş artifact bulunamadı veya hash'i uyuşmuyor.")
    feature_columns = [str(value) for value in model_info["feature_columns"]]
    categorical_columns = [str(value) for value in model_info["categorical_columns"]]
    required = ["VariationID", *feature_columns, "target", "protein_effect"]
    frame = pd.read_csv(features_path, sep="\t", usecols=required, dtype={"VariationID": "string", "protein_effect": "string", **{column: "string" for column in categorical_columns}})
    if frame["VariationID"].isna().any() or frame["VariationID"].duplicated().any():
        raise ValueError("Dış özellik tablosunda VariationID tekil olmalıdır.")
    targets = frame["target"].astype(int)
    if set(targets.unique()) != {0, 1}:
        raise ValueError("Dış kohort iki hedef sınıfı da içermelidir.")
    model = joblib.load(artifact_path)
    probabilities = model.predict_proba(frame.loc[:, feature_columns])[:, 1]
    threshold = float(candidate["threshold"])
    metrics = calculate_metrics(targets, pd.Series(probabilities), threshold=threshold)
    intervals = bootstrap_intervals(targets.to_numpy(), probabilities, seed=seed, replicates=bootstrap_replicates)
    predictions = (probabilities >= threshold).astype(int)
    report_dir.mkdir(parents=True)
    predictions_path = report_dir / "external_predictions.tsv.gz"
    pd.DataFrame({"VariationID": frame["VariationID"], "target": targets, "probability_pathogenic": probabilities, "prediction": predictions}).to_csv(predictions_path, sep="\t", index=False, compression="gzip", float_format="%.10f")
    manifest: dict[str, object] = {
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "evaluation_number": 1,
        "candidate_manifest": {"path": str(candidate_path.resolve()), "sha256": sha256_file(candidate_path)},
        "external_features": {"path": str(features_path.resolve()), "sha256": sha256_file(features_path)},
        "frozen_artifact": {"path": str(artifact_path.resolve()), "sha256": sha256_file(artifact_path)},
        "data_policy": {"model_refit_after_external_labels": False, "threshold_changed_after_external_labels": False},
        "metrics": metrics,
        "bootstrap_95_ci": intervals,
        "outputs": {"predictions": str(predictions_path.resolve())},
    }
    report_path = report_dir / "external_evaluation.md"
    report_path.write_text(
        "# Dış kohort değerlendirmesi\n\n"
        "Dondurulmuş artifact ve sabit eşik kullanıldı; dış etiketlerden sonra yeniden eğitim veya eşik ayarı yapılmadı.\n\n"
        "| Metrik | Değer | %95 GA |\n| --- | ---: | --- |\n"
        f"| ROC-AUC | {metrics['roc_auc']:.4f} | {intervals['roc_auc']['lower_95']:.4f}–{intervals['roc_auc']['upper_95']:.4f} |\n"
        f"| Average precision | {metrics['average_precision']:.4f} | {intervals['average_precision']['lower_95']:.4f}–{intervals['average_precision']['upper_95']:.4f} |\n"
        f"| Duyarlılık (eşik={threshold:.2f}) | {metrics['recall_sensitivity']:.4f} | — |\n"
        f"| Özgüllük (eşik={threshold:.2f}) | {metrics['specificity']:.4f} | — |\n",
        encoding="utf-8",
    )
    manifest["outputs"]["report"] = str(report_path.resolve())
    with (report_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print("Dış kohort değerlendirmesi tamamlandı")
    print(f"ROC-AUC={metrics['roc_auc']:.4f}; AP={metrics['average_precision']:.4f}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Dondurulmuş adayı dış özellik tablosunda tek seferde ölçer.")
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--bootstrap-replicates", type=int, default=DEFAULT_BOOTSTRAP_REPLICATES)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    evaluate_external_candidate(args.candidate, args.features, args.report_dir, bootstrap_replicates=args.bootstrap_replicates, seed=args.seed)


if __name__ == "__main__":
    main()
