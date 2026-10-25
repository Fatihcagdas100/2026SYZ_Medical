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
DEFAULT_ADMISSION = Path("reports/generated/external_admission/manifest.json")
DEFAULT_FEATURE_MANIFEST = Path("reports/generated/external_features/manifest.json")
DEFAULT_REPORT_DIR = Path("reports/generated/external_evaluation")
DEFAULT_BOOTSTRAP_REPLICATES = 1_000
DEFAULT_SEED = 2027
CI_METRICS = ("roc_auc", "average_precision", "precision", "recall_sensitivity", "specificity", "f1")


def load_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError("Manifest JSON nesnesi olmalıdır.")
    return value


def bootstrap_intervals(targets: np.ndarray, probabilities: np.ndarray, *, seed: int, replicates: int, threshold: float = 0.5) -> dict[str, dict[str, float]]:
    if replicates < 1000:
        raise ValueError("Bootstrap tekrar sayısı en az 1000 olmalıdır.")
    if set(np.unique(targets)) != {0, 1}:
        raise ValueError("Bootstrap için iki sınıf da gereklidir.")
    rng = np.random.default_rng(seed)
    negative = np.flatnonzero(targets == 0)
    positive = np.flatnonzero(targets == 1)
    values: dict[str, list[float]] = {metric: [] for metric in CI_METRICS}
    for _ in range(replicates):
        indices = np.concatenate((rng.choice(negative, len(negative), replace=True), rng.choice(positive, len(positive), replace=True)))
        sample_target = targets[indices]
        sample_probability = probabilities[indices]
        values["roc_auc"].append(float(roc_auc_score(sample_target, sample_probability)))
        values["average_precision"].append(float(average_precision_score(sample_target, sample_probability)))
        prediction = sample_probability >= threshold
        tp = int(np.sum((sample_target == 1) & prediction))
        fp = int(np.sum((sample_target == 0) & prediction))
        fn = len(positive) - tp
        values["precision"].append(tp / (tp + fp) if tp + fp else 0.0)
        values["recall_sensitivity"].append(tp / len(positive))
        values["specificity"].append((len(negative) - fp) / len(negative))
        values["f1"].append(2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0)
    return {
        metric: {"lower_95": float(np.quantile(scores, 0.025)), "upper_95": float(np.quantile(scores, 0.975))}
        for metric, scores in values.items()
    }


def verified_file(info: dict) -> Path:
    digest = info.get("sha256")
    if not isinstance(digest, str) or len(digest) != 64:
        raise ValueError("Dosya SHA-256 alanı eksik veya geçersiz.")
    path = Path(info["path"])
    if not path.is_file() or sha256_file(path) != digest:
        raise ValueError(f"Dosya hash uyuşmazlığı: {path}")
    return path


def start_evaluation(candidate_path: Path, admission: dict, report_dir: Path) -> Path:
    """Aynı ham kohortu farklı çıktı klasörlerinde de kilitle; hata sonrası silme."""
    cohort_hash = admission["inputs"]["external_candidate"]["sha256"]
    registry_dir = candidate_path.parent / "evaluation_registry"
    registry_dir.mkdir(parents=True, exist_ok=True)
    marker_path = registry_dir / f"{cohort_hash}.json"
    marker = {
        "status": "started", "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "candidate_sha256": sha256_file(candidate_path),
        "external_candidate_sha256": cohort_hash,
        "report_dir": str(report_dir.resolve()),
        "policy": "Başarısız veya yarım kalmış koşu da yeniden başlatılmaz.",
    }
    with marker_path.open("x", encoding="utf-8") as handle:
        json.dump(marker, handle, ensure_ascii=False, indent=2)
    report_dir.mkdir(parents=True, exist_ok=True)
    with (report_dir / "EVALUATION_STARTED.json").open("x", encoding="utf-8") as handle:
        json.dump(marker, handle, ensure_ascii=False, indent=2)
    return marker_path


def subgroup_result(frame: pd.DataFrame, *, threshold: float, replicates: int, seed: int) -> dict:
    counts = {"row_count": len(frame), "negative_count": int(frame["target"].eq(0).sum()),
              "positive_count": int(frame["target"].eq(1).sum())}
    if frame["target"].nunique() != 2:
        return {"status": "insufficient_classes", "counts": counts,
                "metrics": None, "bootstrap_95_ci": None}
    return {
        "status": "evaluated", "counts": counts,
        "metrics": calculate_metrics(frame["target"], frame["probability_pathogenic"], threshold=threshold),
        "bootstrap_95_ci": bootstrap_intervals(
            frame["target"].to_numpy(), frame["probability_pathogenic"].to_numpy(),
            threshold=threshold, seed=seed, replicates=replicates,
        ),
    }


def evaluate_external_candidate(
    candidate_path: Path,
    features_path: Path,
    admission_path: Path,
    feature_manifest_path: Path,
    report_dir: Path,
    *,
    bootstrap_replicates: int = DEFAULT_BOOTSTRAP_REPLICATES,
    seed: int = DEFAULT_SEED,
) -> dict[str, object]:
    if bootstrap_replicates < 1000:
        raise ValueError("Bootstrap tekrar sayısı en az 1000 olmalıdır.")
    if report_dir.exists() and any(report_dir.iterdir()):
        raise FileExistsError("Dış değerlendirme rapor klasörü boş olmalıdır; tek-sefer kuralı geçerlidir.")
    candidate = load_json(candidate_path)
    admission = load_json(admission_path)
    feature_manifest = load_json(feature_manifest_path)
    if candidate.get("status") != "frozen_for_external_validation_not_evaluated":
        raise ValueError("Aday model dış değerlendirme için dondurulmuş durumda değil.")
    if admission.get("admission_decision") != "accepted" or admission.get("schema_version") != 2:
        raise ValueError("Referans denetimli v2 kabul manifesti gerekli.")
    if admission.get("reference", {}).get("assembly") != "GRCh38":
        raise ValueError("Kabul referansı GRCh38 değil.")
    verified_file(admission["inputs"]["development"])
    verified_file(admission["inputs"]["external_candidate"])
    admitted_path = verified_file(admission["output"])
    verified_file(feature_manifest["input"])
    verified_file(feature_manifest["output"])
    annotation_path = verified_file(admission["annotations"])
    if admission.get("output", {}).get("sha256") != feature_manifest.get("input", {}).get("sha256"):
        raise ValueError("Kabul edilen kohort ile özellik kaynağı hash zinciri uyuşmuyor.")
    if feature_manifest.get("output", {}).get("sha256") != sha256_file(features_path):
        raise ValueError("Dış özellik dosyası hash'i özellik manifestiyle uyuşmuyor.")
    model_info = candidate.get("model")
    if not isinstance(model_info, dict):
        raise ValueError("Aday model bilgisi eksik.")
    artifact_path = Path(str(model_info["artifact_path"]))
    if not artifact_path.is_file() or sha256_file(artifact_path) != model_info["artifact_sha256"]:
        raise ValueError("Dondurulmuş artifact bulunamadı veya hash'i uyuşmuyor.")
    feature_columns = [str(value) for value in model_info["feature_columns"]]
    categorical_columns = [str(value) for value in model_info["categorical_columns"]]
    threshold = float(candidate["threshold"])
    if not np.isfinite(threshold) or not 0 < threshold < 1:
        raise ValueError("Dondurulmuş eşik geçersiz.")
    marker_path = start_evaluation(candidate_path, admission, report_dir)
    required = list(dict.fromkeys(["VariationID", *feature_columns, "target", "protein_effect", "has_protein_annotation"]))
    frame = pd.read_csv(features_path, sep="\t", usecols=required, dtype={"VariationID": "string", "target": "string", "protein_effect": "string", **{column: "string" for column in categorical_columns}})
    if frame["VariationID"].isna().any() or frame["VariationID"].duplicated().any():
        raise ValueError("Dış özellik tablosunda VariationID tekil olmalıdır.")
    if frame["target"].isna().any() or set(frame["target"].unique()) != {"0", "1"}:
        raise ValueError("Dış hedefler tam olarak 0 ve 1 olmalıdır.")
    frame["target"] = frame["target"].astype(int)
    annotations = pd.read_csv(annotation_path, sep="\t", dtype={"VariationID": "string"})
    frame = frame.merge(annotations, on="VariationID", how="outer", validate="one_to_one", indicator=True)
    if not frame["_merge"].eq("both").all():
        raise ValueError("Kabul anotasyonları ve özellik ID'leri eşleşmiyor.")
    admitted = pd.read_csv(admitted_path, sep="\t", usecols=["VariationID", "target"], dtype={"VariationID": "string"})
    joined = frame[["VariationID", "target"]].merge(admitted, on="VariationID", how="outer", validate="one_to_one", suffixes=("_features", "_admitted"), indicator=True)
    if not joined["_merge"].eq("both").all() or not joined["target_features"].eq(joined["target_admitted"]).all():
        raise ValueError("Kabul kohortu ile özellik hedefleri/ID'leri eşleşmiyor.")
    targets = frame["target"]
    model = joblib.load(artifact_path)
    probabilities = model.predict_proba(frame.loc[:, feature_columns])[:, 1]
    if not np.isfinite(probabilities).all() or not ((probabilities >= 0) & (probabilities <= 1)).all():
        raise ValueError("Model geçersiz olasılık üretti.")
    metrics = calculate_metrics(targets, pd.Series(probabilities), threshold=threshold)
    intervals = bootstrap_intervals(targets.to_numpy(), probabilities, seed=seed, replicates=bootstrap_replicates, threshold=threshold)
    predictions = (probabilities >= threshold).astype(int)
    frame["probability_pathogenic"] = probabilities
    subgroups = {
        name: subgroup_result(group, threshold=threshold, replicates=bootstrap_replicates, seed=seed)
        for name, group in {
            "missense": frame.loc[frame["protein_effect"].eq("missense")],
            "unseen_genes": frame.loc[frame["gene_status"].eq("unseen")],
        }.items()
    }
    predictions_path = report_dir / "external_predictions.tsv.gz"
    pd.DataFrame({"VariationID": frame["VariationID"], "target": targets, "probability_pathogenic": probabilities, "prediction": predictions}).to_csv(predictions_path, sep="\t", index=False, compression="gzip", float_format="%.10f")
    manifest: dict[str, object] = {
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "evaluation_number": 1,
        "cohort_kind": admission["cohort_kind"],
        "frozen_threshold": threshold,
        "candidate_manifest": {"path": str(candidate_path.resolve()), "sha256": sha256_file(candidate_path)},
        "admission_manifest": {"path": str(admission_path.resolve()), "sha256": sha256_file(admission_path)},
        "feature_manifest": {"path": str(feature_manifest_path.resolve()), "sha256": sha256_file(feature_manifest_path)},
        "external_features": {"path": str(features_path.resolve()), "sha256": sha256_file(features_path)},
        "frozen_artifact": {"path": str(artifact_path.resolve()), "sha256": sha256_file(artifact_path)},
        "data_policy": {"model_refit_after_external_labels": False, "threshold_changed_after_external_labels": False},
        "metrics": metrics,
        "bootstrap_95_ci": intervals,
        "subgroups": subgroups,
        "bootstrap": {"replicates": bootstrap_replicates, "seed": seed, "method": "class_stratified_row_percentile", "limitations": "Sınıf sayıları sabittir; gen içi bağımlılık bu aralıklarda modellenmez."},
        "annotation_coverage": {
            "protein_annotation_missing_rows": int(frame["has_protein_annotation"].eq(0).sum()),
            "gene_annotation_missing_rows": int(frame["gene_status"].eq("missing").sum()),
        },
        "outputs": {"predictions": str(predictions_path.resolve()), "started_marker": str(marker_path.resolve())},
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
    with report_path.open("a", encoding="utf-8") as handle:
        handle.write(f"\nKohort türü: `{admission['cohort_kind']}`\n\n")
        for group_name, result in {"all_variants": {"metrics": metrics, "bootstrap_95_ci": intervals}, **subgroups}.items():
            handle.write(f"## {group_name}\n\n")
            if result["metrics"] is None:
                handle.write(f"İki sınıf birlikte yok; metrik hesaplanmadı. Sayımlar: {result['counts']}\n\n")
                continue
            handle.write("| Metrik | Değer | %95 GA |\n| --- | ---: | --- |\n")
            for metric in CI_METRICS:
                interval = result["bootstrap_95_ci"][metric]
                handle.write(f"| {metric} | {result['metrics'][metric]:.4f} | {interval['lower_95']:.4f}–{interval['upper_95']:.4f} |\n")
            m = result["metrics"]
            handle.write(f"\nKarışıklık sayıları: TN={m['true_negatives']}, FP={m['false_positives']}, FN={m['false_negatives']}, TP={m['true_positives']}. n={m['row_count']}\n\n")
        handle.write(f"Anotasyon kapsamı: {manifest['annotation_coverage']}\n\n{manifest['bootstrap']['limitations']}\n")
    manifest["outputs"]["report"] = str(report_path.resolve())
    with (report_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2, allow_nan=False)
        handle.write("\n")
    print("Dış kohort değerlendirmesi tamamlandı")
    print(f"ROC-AUC={metrics['roc_auc']:.4f}; AP={metrics['average_precision']:.4f}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Dondurulmuş adayı dış özellik tablosunda tek seferde ölçer.")
    parser.add_argument("--candidate", type=Path, default=DEFAULT_CANDIDATE)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--admission-manifest", type=Path, default=DEFAULT_ADMISSION)
    parser.add_argument("--feature-manifest", type=Path, default=DEFAULT_FEATURE_MANIFEST)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--bootstrap-replicates", type=int, default=DEFAULT_BOOTSTRAP_REPLICATES)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    evaluate_external_candidate(args.candidate, args.features, args.admission_manifest, args.feature_manifest, args.report_dir, bootstrap_replicates=args.bootstrap_replicates, seed=args.seed)


if __name__ == "__main__":
    main()
