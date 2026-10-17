"""Temel ClinVar tablosu için deterministik, sınıf-korumalı ayrım üret."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence


DEFAULT_INPUT = Path("data/processed/clinvar_baseline_features.tsv.gz")
DEFAULT_OUTPUT = Path("data/processed/baseline_stratified_split.tsv")
DEFAULT_REPORT_DIR = Path("reports/generated/baseline_split")
DEFAULT_SEED = 2027
SPLIT_RATIOS = {"train": 0.70, "validation": 0.15, "test": 0.15}
TARGET_LABELS = {"0": "Benign", "1": "Pathogenic"}


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def allocation_counts(total: int, ratios: dict[str, float]) -> dict[str, int]:
    """En büyük kalan yöntemiyle oranları toplamı tam koruyan sayılara çevir."""

    if total < 0:
        raise ValueError("Toplam negatif olamaz.")
    if not math.isclose(sum(ratios.values()), 1.0, rel_tol=0, abs_tol=1e-12):
        raise ValueError("Ayrım oranlarının toplamı 1 olmalıdır.")
    raw = {name: total * ratio for name, ratio in ratios.items()}
    counts = {name: math.floor(value) for name, value in raw.items()}
    remaining = total - sum(counts.values())
    order = sorted(
        ratios,
        key=lambda name: (-(raw[name] - counts[name]), list(ratios).index(name)),
    )
    for name in order[:remaining]:
        counts[name] += 1
    return counts


def deterministic_rank(variation_id: str, seed: int) -> bytes:
    return hashlib.sha256(f"{seed}\x1f{variation_id}".encode("utf-8")).digest()


def assign_splits(
    records: Sequence[tuple[str, str]],
    *,
    seed: int,
    ratios: dict[str, float] = SPLIT_RATIOS,
) -> tuple[dict[str, str], dict[str, dict[str, int]]]:
    by_target: defaultdict[str, list[str]] = defaultdict(list)
    for variation_id, target in records:
        by_target[target].append(variation_id)

    assignments: dict[str, str] = {}
    counts_by_target: dict[str, dict[str, int]] = {}
    for target in sorted(by_target):
        ranked_ids = sorted(
            by_target[target],
            key=lambda variation_id: (deterministic_rank(variation_id, seed), variation_id),
        )
        counts = allocation_counts(len(ranked_ids), ratios)
        counts_by_target[target] = counts
        start = 0
        for split_name in ratios:
            stop = start + counts[split_name]
            for variation_id in ranked_ids[start:stop]:
                assignments[variation_id] = split_name
            start = stop
        if start != len(ranked_ids):
            raise RuntimeError("Ayrım sırasında bazı kimlikler atanmadı.")
    return assignments, counts_by_target


def build_split(
    input_path: Path,
    output_path: Path,
    report_dir: Path,
    *,
    seed: int = DEFAULT_SEED,
) -> dict[str, object]:
    if not input_path.is_file():
        raise FileNotFoundError(f"Girdi dosyası bulunamadı: {input_path}")

    records: list[tuple[str, str]] = []
    seen_ids: set[str] = set()
    with gzip.open(input_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Girdi dosyasında başlık bulunamadı.")
        missing = sorted({"VariationID", "target"}.difference(reader.fieldnames))
        if missing:
            raise ValueError("Gerekli sütunlar bulunamadı: " + ", ".join(missing))
        for row in reader:
            variation_id = row["VariationID"]
            target = row["target"]
            if not variation_id:
                raise ValueError("Boş VariationID bulundu.")
            if variation_id in seen_ids:
                raise ValueError(f"Tekrarlanan VariationID bulundu: {variation_id}")
            if target not in TARGET_LABELS:
                raise ValueError(f"Beklenmeyen target değeri: {target}")
            seen_ids.add(variation_id)
            records.append((variation_id, target))

    assignments, counts_by_target = assign_splits(records, seed=seed)
    if len(assignments) != len(records):
        raise RuntimeError("Atanan kimlik sayısı girdiyle uyuşmuyor.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    temporary_handle = tempfile.NamedTemporaryFile(
        prefix="split_",
        suffix=".tsv.tmp",
        dir=output_path.parent,
        delete=False,
    )
    temporary_path = Path(temporary_handle.name)
    temporary_handle.close()

    try:
        with temporary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["VariationID", "target", "split"])
            for variation_id, target in records:
                writer.writerow([variation_id, target, assignments[variation_id]])
        temporary_path.replace(output_path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise

    split_totals: Counter[str] = Counter(assignments.values())
    target_totals: Counter[str] = Counter(target for _, target in records)
    with (report_dir / "split_counts.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["split", "target", "target_label", "count", "percentage_within_target"]
        )
        for target in sorted(counts_by_target):
            for split_name in SPLIT_RATIOS:
                count = counts_by_target[target][split_name]
                percentage = 100 * count / target_totals[target]
                writer.writerow(
                    [split_name, target, TARGET_LABELS[target], count, f"{percentage:.6f}"]
                )

    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "path": str(input_path.resolve()),
            "size_bytes": input_path.stat().st_size,
            "sha256": sha256_file(input_path),
            "row_count": len(records),
        },
        "output": {
            "path": str(output_path.resolve()),
            "size_bytes": output_path.stat().st_size,
            "sha256": sha256_file(output_path),
            "columns": ["VariationID", "target", "split"],
        },
        "method": {
            "name": "deterministic_stratified_hash_ranking",
            "seed": seed,
            "hash": "SHA-256",
            "hash_input": "seed + unit-separator + VariationID",
            "allocation": "largest remainder within each target class",
            "ratios": SPLIT_RATIOS,
        },
        "counts": {
            "total": len(records),
            "by_split": dict(sorted(split_totals.items())),
            "by_target": dict(sorted(target_totals.items())),
            "by_target_and_split": counts_by_target,
        },
        "integrity": {
            "unique_variation_ids": len(seen_ids),
            "assigned_variation_ids": len(assignments),
            "unassigned_variation_ids": len(seen_ids.difference(assignments)),
            "split_overlap_count": 0,
        },
        "evaluation_scope": (
            "Varyant düzeyinde rastgele-benzeri ayrım; görülmemiş genlere "
            "genelleme iddiası taşımaz."
        ),
    }
    with (report_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("Sabit veri ayrımı oluşturuldu")
    print(f"Toplam: {len(records):,}")
    for split_name in SPLIT_RATIOS:
        print(f"{split_name}: {split_totals[split_name]:,}")
    print(f"Seed: {seed}")
    print(f"Çıktı: {output_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Temel özellik tablosu için sabit sınıf-korumalı ayrım üretir."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    build_split(
        args.input,
        args.output,
        args.report_dir,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()

