"""HGVS ayrıştırıcısını ham Name örnekleri ve yeniden üretim denetimiyle incele."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import itertools
import json
import re
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

try:
    from scripts.build_hgvs_features import HGVS_FEATURE_COLUMNS, derive_hgvs_features
except ModuleNotFoundError:  # Doğrudan ``python scripts/...`` çalıştırması.
    from build_hgvs_features import HGVS_FEATURE_COLUMNS, derive_hgvs_features


DEFAULT_SOURCE = Path("data/processed/clinvar_grch38_germline_exact_binary.tsv.gz")
DEFAULT_FEATURES = Path("data/processed/clinvar_hgvs_features.tsv.gz")
DEFAULT_SPLIT = Path("data/processed/baseline_gene_disjoint_split.tsv")
DEFAULT_REPORT_DIR = Path("reports/generated/hgvs_parser_audit")
DEFAULT_SEED = 2027
DEFAULT_SAMPLES_PER_CLASS = 10
AUDIT_CLASSES = (
    "missense",
    "synonymous",
    "frameshift",
    "stop_gained_or_changed",
    "no_protein_annotation",
    "uncertain_or_absent",
)


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def load_split_assignments(path: Path) -> dict[str, str]:
    assignments: dict[str, str] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Ayrım dosyasında başlık bulunamadı.")
        missing = {"VariationID", "split"}.difference(reader.fieldnames)
        if missing:
            raise ValueError("Ayrım dosyasında eksik sütunlar: " + ", ".join(sorted(missing)))
        for row in reader:
            variation_id = row["VariationID"]
            if variation_id in assignments:
                raise ValueError(f"Tekrarlanan split kimliği: {variation_id}")
            if row["split"] not in {"train", "validation", "test"}:
                raise ValueError(f"Beklenmeyen split: {row['split']}")
            assignments[variation_id] = row["split"]
    return assignments


def extract_protein_token(name: str) -> str | None:
    """Örnek denetimi için üretim ayrıştırıcısından bağımsız basit token çıkarımı."""

    start = name.casefold().find("(p.")
    if start < 0:
        return None
    end = name.find(")", start + 3)
    if end < 0:
        return None
    return name[start + 3 : end].strip()


def marker_consistent(expected_class: str, name: str) -> bool:
    token = extract_protein_token(name)
    if expected_class == "no_protein_annotation":
        return token is None
    if token is None:
        return False
    if expected_class == "missense":
        return re.fullmatch(r"[A-Z][a-z]{2}\d+[A-Z][a-z]{2}", token) is not None
    if expected_class == "synonymous":
        return token.endswith("=")
    if expected_class == "frameshift":
        return "fs" in token.casefold()
    if expected_class == "stop_gained_or_changed":
        return "ter" in token.casefold() or "*" in token
    if expected_class == "uncertain_or_absent":
        return "?" in token or token == "0"
    raise ValueError(f"Denetlenmeyen protein sınıfı: {expected_class}")


def stable_rank(seed: int, variation_id: str) -> int:
    payload = f"{seed}:{variation_id}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest(), "big")


def audit_hgvs_parser(
    source_path: Path,
    features_path: Path,
    split_path: Path,
    report_dir: Path,
    *,
    seed: int = DEFAULT_SEED,
    samples_per_class: int = DEFAULT_SAMPLES_PER_CLASS,
) -> dict[str, object]:
    if samples_per_class <= 0:
        raise ValueError("samples_per_class pozitif olmalıdır.")
    for path in (source_path, features_path, split_path):
        if not path.is_file():
            raise FileNotFoundError(f"Girdi bulunamadı: {path}")

    assignments = load_split_assignments(split_path)
    report_dir.mkdir(parents=True, exist_ok=True)
    counters: Counter[str] = Counter()
    mismatch_counts: Counter[str] = Counter()
    class_counts: Counter[str] = Counter()
    # Her sınıfta en küçük sabit hash değerlerine sahip kayıtları tutar.
    sample_heaps: dict[str, list[tuple[int, str, dict[str, str]]]] = {
        effect: [] for effect in AUDIT_CLASSES
    }
    seen_ids: set[str] = set()

    with (
        gzip.open(source_path, "rt", encoding="utf-8", newline="") as source_handle,
        gzip.open(features_path, "rt", encoding="utf-8", newline="") as feature_handle,
    ):
        source_reader = csv.DictReader(source_handle, delimiter="\t")
        feature_reader = csv.DictReader(feature_handle, delimiter="\t")
        if source_reader.fieldnames is None or feature_reader.fieldnames is None:
            raise ValueError("Girdilerden birinde başlık bulunamadı.")
        source_missing = {"VariationID", "Name", "target"}.difference(
            source_reader.fieldnames
        )
        feature_missing = {"VariationID", "target", *HGVS_FEATURE_COLUMNS}.difference(
            feature_reader.fieldnames
        )
        if source_missing or feature_missing:
            raise ValueError(
                "Eksik sütunlar: "
                + ", ".join(sorted(source_missing.union(feature_missing)))
            )

        for source_row, feature_row in itertools.zip_longest(source_reader, feature_reader):
            if source_row is None or feature_row is None:
                raise ValueError("Kaynak ve özellik satır sayıları eşleşmiyor.")
            variation_id = source_row["VariationID"]
            if variation_id != feature_row["VariationID"]:
                raise ValueError(
                    f"Satır sırası eşleşmiyor: {variation_id} != "
                    f"{feature_row['VariationID']}"
                )
            if variation_id in seen_ids:
                raise ValueError(f"Tekrarlanan VariationID: {variation_id}")
            seen_ids.add(variation_id)
            split_name = assignments.get(variation_id)
            if split_name is None:
                raise ValueError(f"Split ataması olmayan VariationID: {variation_id}")
            counters[f"{split_name}_rows_seen"] += 1

            # Test satırları yalnızca dosya hizasını doğrulamak için geçilir; ayrıştırılmaz,
            # sayılmaz, örneklenmez ve hedefleri analize alınmaz.
            if split_name == "test":
                continue
            counters["development_rows_audited"] += 1
            derived = derive_hgvs_features(source_row["Name"])
            for column in HGVS_FEATURE_COLUMNS:
                if str(feature_row[column]) != str(derived[column]):
                    mismatch_counts[column] += 1

            effect = derived["protein_effect"]
            class_counts[effect] += 1
            if effect not in sample_heaps:
                continue
            rank = stable_rank(seed, variation_id)
            sample = {
                "VariationID": variation_id,
                "split": split_name,
                "target": source_row["target"],
                "Name": source_row["Name"],
                **derived,
                "independent_protein_token": extract_protein_token(source_row["Name"])
                or "",
                "automatic_marker_consistent": str(
                    marker_consistent(effect, source_row["Name"])
                ).lower(),
            }
            heap = sample_heaps[effect]
            entry = (-rank, variation_id, sample)
            if len(heap) < samples_per_class:
                heapq.heappush(heap, entry)
            elif rank < -heap[0][0]:
                heapq.heapreplace(heap, entry)

    missing_ids = set(assignments).difference(seen_ids)
    if missing_ids:
        raise ValueError(f"Özellik girdisinde bulunmayan split kimliği: {min(missing_ids)}")

    samples = [
        entry[2]
        for effect in AUDIT_CLASSES
        for entry in sorted(sample_heaps[effect], key=lambda item: -item[0])
    ]
    samples.sort(key=lambda row: (AUDIT_CLASSES.index(row["protein_effect"]), row["VariationID"]))
    sample_path = report_dir / "sampled_raw_name_examples.tsv"
    sample_fields = [
        "VariationID",
        "split",
        "target",
        "Name",
        *HGVS_FEATURE_COLUMNS,
        "independent_protein_token",
        "automatic_marker_consistent",
    ]
    with sample_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sample_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(samples)

    distribution_path = report_dir / "development_protein_effect_counts.csv"
    with distribution_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["protein_effect", "row_count"])
        for effect, count in sorted(class_counts.items(), key=lambda item: (-item[1], item[0])):
            writer.writerow([effect, count])

    marker_mismatches = sum(
        row["automatic_marker_consistent"] != "true" for row in samples
    )
    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "source": {"path": str(source_path.resolve()), "sha256": sha256_file(source_path)},
            "features": {"path": str(features_path.resolve()), "sha256": sha256_file(features_path)},
            "split": {"path": str(split_path.resolve()), "sha256": sha256_file(split_path)},
        },
        "data_policy": {
            "audited_partitions": ["train", "validation"],
            "test_partition": "skipped_after_row_alignment_check",
            "test_target_or_hgvs_features_used": False,
            "sampling": "smallest_sha256(seed:VariationID)_within_each_class",
            "seed": seed,
            "samples_per_class": samples_per_class,
        },
        "counts": dict(counters),
        "development_protein_effect_counts": dict(class_counts),
        "saved_feature_reproduction_mismatches": {
            column: mismatch_counts[column] for column in HGVS_FEATURE_COLUMNS
        },
        "sample_audit": {
            "requested_classes": list(AUDIT_CLASSES),
            "sample_count": len(samples),
            "automatic_marker_mismatch_count": marker_mismatches,
            "limitation": (
                "The marker check is an independent syntactic sanity check, not a "
                "clinical or transcript-aware expert annotation review."
            ),
        },
        "reports": {
            "sampled_raw_name_examples": str(sample_path.resolve()),
            "development_protein_effect_counts": str(distribution_path.resolve()),
        },
    }
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    report_path = report_dir / "audit_report.md"
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write("# HGVS ayrıştırıcı denetimi\n\n")
        handle.write(
            f"Train/validation üzerinde {counters['development_rows_audited']:,} kayıt "
            "yeniden ayrıştırıldı. Test satırları analize alınmadı.\n\n"
        )
        handle.write(
            f"Kaydedilmiş özelliklerle toplam uyuşmazlık: {sum(mismatch_counts.values()):,}. "
            f"Seçilmiş {len(samples)} ham Name örneğinde basit bağımsız işaret denetimi "
            f"uyuşmazlığı: {marker_mismatches}.\n\n"
        )
        handle.write(
            "Bu otomatik kontrol sözdizimsel bir sağlamlık testidir; klinik veya "
            "transkript-duyarlı uzman doğrulamasının yerine geçmez. Ham örnekler "
            "`sampled_raw_name_examples.tsv` dosyasındadır.\n"
        )
    manifest["reports"]["audit_report"] = str(report_path.resolve())
    # Rapor yolu da manifestte yer alsın.
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("HGVS ayrıştırıcı denetimi tamamlandı")
    print(f"Geliştirme kaydı: {counters['development_rows_audited']:,}")
    print(f"Kaydedilmiş özellik uyuşmazlığı: {sum(mismatch_counts.values()):,}")
    print(f"Örnek işaret uyuşmazlığı: {marker_mismatches:,}/{len(samples):,}")
    print(f"Rapor: {report_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="HGVS özellik ayrıştırıcısını denetler.")
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--split", type=Path, default=DEFAULT_SPLIT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--samples-per-class", type=int, default=DEFAULT_SAMPLES_PER_CLASS)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    audit_hgvs_parser(
        args.source,
        args.features,
        args.split,
        args.report_dir,
        seed=args.seed,
        samples_per_class=args.samples_per_class,
    )


if __name__ == "__main__":
    main()
