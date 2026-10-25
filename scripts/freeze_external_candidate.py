"""Seçilmiş profili yalnızca geliştirme verisiyle dış doğrulama için dondur."""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import joblib
import pandas as pd

try:
    from scripts.train_logistic_baseline import (
        FEATURE_PROFILES,
        build_pipeline,
        load_development_data,
        sha256_file,
    )
except ModuleNotFoundError:  # Doğrudan `python scripts/...` çalıştırması için.
    from train_logistic_baseline import (  # type: ignore[no-redef]
        FEATURE_PROFILES,
        build_pipeline,
        load_development_data,
        sha256_file,
    )


DEFAULT_FEATURES = Path("data/processed/clinvar_missense_features.tsv.gz")
DEFAULT_SPLIT = Path("data/processed/baseline_gene_disjoint_split.tsv")
DEFAULT_REPORT_DIR = Path("reports/generated/external_candidate")
DEFAULT_ARTIFACT_DIR = Path("artifacts/external_candidate")
DEFAULT_SEED = 2027
DEFAULT_THRESHOLD = 0.5
DEFAULT_CHUNK_SIZE = 100_000
DEFAULT_FEATURE_PROFILE = "hgvs_missense"
MODEL_NAME = "hgvs_missense_development_refit"


def format_candidate_card(manifest: dict[str, object]) -> str:
    model = manifest["model"]
    counts = manifest["partition_counts"]
    assert isinstance(model, dict)
    assert isinstance(counts, dict)
    return f"""# Dış doğrulama adayı model kartı

Durum: **DIŞ KOHORT GELMEDEN ÖNCE DONDURULDU — DEĞERLENDİRİLMEDİ**

Bu model yalnızca araştırma amaçlı varyant önceliklendirme prototipidir; klinik
tanı koymaz ve uzman değerlendirmesinin yerine geçmez.

## Dondurulan tarif

- Model: `{MODEL_NAME}`
- Özellik profili: `{manifest['feature_profile']}`
- Eğitim satırları: `{counts['train']:,}` train + `{counts['validation']:,}` validation
- Dışlanan kilitli test satırları: `{counts['test']:,}`
- Seed: `{manifest['seed']}`
- Karar eşiği: `{manifest['threshold']:.4f}` (önceden sabit; dış kohorta göre ayarlanmayacak)
- Artifact SHA-256: `{model['artifact_sha256']}`

## Veri politikası

Profil seçimi, sabit validation ve train-only gene-disjoint OOF sonuçlarıyla
tamamlanmıştır. Bu artifact daha sonra yalnızca önceden kabul edilmiş, geliştirme
verisiyle varyant/allel/koordinat çakışması olmayan dış kohortta ölçülecektir.
Kilitli iç test bu eğitimde veya seçimde kullanılmamıştır.

## Yasak işlem

Dış kohortun etiketleri görüldükten sonra bu artifact, özellik profili veya karar
eşiği değiştirilmez. Yeni bir deneme gerekirse ayrı bir dış kohort ve ayrı bir
protokol gerekir.
"""


def freeze_external_candidate(
    features_path: Path,
    split_path: Path,
    report_dir: Path,
    artifact_dir: Path,
    *,
    seed: int = DEFAULT_SEED,
    threshold: float = DEFAULT_THRESHOLD,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    feature_profile: str = DEFAULT_FEATURE_PROFILE,
) -> dict[str, object]:
    """Train ve validation üzerinde tek adımda yeniden eğit; test satırı alma."""
    for directory in (report_dir, artifact_dir):
        if directory.exists() and any(directory.iterdir()):
            raise FileExistsError(f"Dondurulmuş adayın üzerine yazılamaz: {directory}")
    if not 0 < threshold < 1:
        raise ValueError("Eşik 0 ile 1 arasında olmalıdır.")
    if feature_profile not in FEATURE_PROFILES:
        raise ValueError(
            f"Bilinmeyen özellik profili: {feature_profile}. "
            f"Seçenekler: {', '.join(FEATURE_PROFILES)}"
        )
    profile = FEATURE_PROFILES[feature_profile]
    categorical_features = tuple(profile["categorical"])
    feature_columns = (*categorical_features, *tuple(profile["numeric"]))
    development, partition_counts = load_development_data(
        features_path,
        split_path,
        feature_columns=feature_columns,
        categorical_features=categorical_features,
        chunk_size=chunk_size,
    )
    if set(development["split"].unique()) != {"train", "validation"}:
        raise RuntimeError("Aday eğitim çerçevesi yalnızca train ve validation içermelidir.")

    pipeline = build_pipeline(
        feature_columns, seed=seed, categorical_features=categorical_features
    )
    pipeline.fit(development.loc[:, feature_columns], development["target"].astype(int))
    classifier = pipeline.named_steps["classifier"]
    preprocessing = pipeline.named_steps["preprocessing"]

    report_dir.mkdir(parents=True, exist_ok=True)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    artifact_path = artifact_dir / f"{MODEL_NAME}.joblib"
    joblib.dump(pipeline, artifact_path)
    coefficients = pd.DataFrame(
        {
            "transformed_feature": preprocessing.get_feature_names_out(),
            "coefficient": classifier.coef_[0],
        }
    )
    coefficients["absolute_coefficient"] = coefficients["coefficient"].abs()
    coefficients = coefficients.sort_values("absolute_coefficient", ascending=False)
    coefficients_path = report_dir / "coefficients.csv"
    coefficients.to_csv(coefficients_path, index=False, float_format="%.10f")

    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "frozen_for_external_validation_not_evaluated",
        "model_name": MODEL_NAME,
        "feature_profile": feature_profile,
        "seed": seed,
        "threshold": threshold,
        "inputs": {
            "features": {
                "path": str(features_path.resolve()),
                "sha256": sha256_file(features_path),
            },
            "gene_disjoint_split": {
                "path": str(split_path.resolve()),
                "sha256": sha256_file(split_path),
            },
        },
        "data_policy": {
            "training_partitions": ["train", "validation"],
            "locked_test_used": False,
            "locked_test_loaded_into_training_dataframe": False,
            "external_labels_used": False,
            "threshold_policy": "fixed_before_external_cohort_evaluation",
            "hyperparameter_search": False,
        },
        "partition_counts": partition_counts,
        "model": {
            "artifact_path": str(artifact_path.resolve()),
            "artifact_sha256": sha256_file(artifact_path),
            "feature_columns": list(feature_columns),
            "categorical_columns": list(categorical_features),
            "transformed_feature_count": len(preprocessing.get_feature_names_out()),
            "iterations": int(classifier.n_iter_[0]),
            "converged": bool(classifier.n_iter_[0] < classifier.max_iter),
        },
        "outputs": {"coefficients": str(coefficients_path.resolve())},
        "limitations": [
            "Bu artifact dış kohort üzerinde henüz değerlendirilmemiştir.",
            "ClinVar tabanlı geliştirme sonuçları bağımsız dış doğrulamanın yerine geçmez.",
            "Çıktı klinik tanı veya bireysel klinik risk olarak yorumlanamaz.",
        ],
    }
    card_path = report_dir / "candidate_card.md"
    card_path.write_text(format_candidate_card(manifest), encoding="utf-8")
    manifest["outputs"]["candidate_card"] = str(card_path.resolve())
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("Dış doğrulama adayı model donduruldu")
    print(f"Geliştirme satırı: {len(development):,}; kilitli test kullanılmadı.")
    print(f"Özellik profili: {feature_profile}; eşik: {threshold:.4f}")
    print(f"Model kartı: {card_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Seçilmiş modeli train+validation ile dış doğrulama için dondurur."
    )
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--split", type=Path, default=DEFAULT_SPLIT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--artifact-dir", type=Path, default=DEFAULT_ARTIFACT_DIR)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--chunk-size", type=int, default=DEFAULT_CHUNK_SIZE)
    parser.add_argument(
        "--feature-profile",
        choices=tuple(FEATURE_PROFILES),
        default=DEFAULT_FEATURE_PROFILE,
    )
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    freeze_external_candidate(
        args.features,
        args.split,
        args.report_dir,
        args.artifact_dir,
        seed=args.seed,
        threshold=args.threshold,
        chunk_size=args.chunk_size,
        feature_profile=args.feature_profile,
    )


if __name__ == "__main__":
    main()
