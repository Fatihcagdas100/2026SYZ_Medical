"""Kabul edilmiş dış ikili kohort için model-uyumlu kimliksiz özellikler üret."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

try:
    from scripts.build_baseline_features import FEATURE_COLUMNS, derive_features, sha256_file
    from scripts.build_hgvs_features import HGVS_FEATURE_COLUMNS, derive_hgvs_features
    from scripts.build_missense_features import (
        MISSENSE_FEATURE_COLUMNS,
        derive_missense_features,
    )
except ModuleNotFoundError:  # Doğrudan `python scripts/...` çalıştırması için.
    from build_baseline_features import FEATURE_COLUMNS, derive_features, sha256_file  # type: ignore[no-redef]
    from build_hgvs_features import HGVS_FEATURE_COLUMNS, derive_hgvs_features  # type: ignore[no-redef]
    from build_missense_features import (  # type: ignore[no-redef]
        MISSENSE_FEATURE_COLUMNS,
        derive_missense_features,
    )


DEFAULT_INPUT = Path("data/external/admitted_binary.tsv.gz")
DEFAULT_OUTPUT = Path("data/external/admitted_missense_features.tsv.gz")
DEFAULT_REPORT_DIR = Path("reports/generated/external_features")
INPUT_COLUMNS = {
    "VariationID", "Type", "Name", "Chromosome", "ReferenceAlleleVCF",
    "AlternateAlleleVCF", "target_label", "target",
}
OUTPUT_COLUMNS = (
    "VariationID", *FEATURE_COLUMNS, *HGVS_FEATURE_COLUMNS,
    *MISSENSE_FEATURE_COLUMNS, "target",
)


def build_external_features(
    input_path: Path,
    output_path: Path,
    report_dir: Path,
    *,
    progress_every: int = 100_000,
) -> dict[str, object]:
    if not input_path.is_file():
        raise FileNotFoundError(f"Dış kohort bulunamadı: {input_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    temp_handle = tempfile.NamedTemporaryFile(
        prefix="external_features_", suffix=".tsv.gz.tmp", dir=output_path.parent,
        delete=False,
    )
    temp_path = Path(temp_handle.name)
    temp_handle.close()
    counts: Counter[str] = Counter()
    targets: Counter[str] = Counter()
    seen_ids: set[str] = set()
    try:
        with (
            gzip.open(input_path, "rt", encoding="utf-8", newline="") as source,
            gzip.open(temp_path, "wt", encoding="utf-8", newline="") as target,
        ):
            reader = csv.DictReader(source, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError("Dış kohortta başlık bulunamadı.")
            missing = sorted(INPUT_COLUMNS.difference(reader.fieldnames))
            if missing:
                raise ValueError("Dış kohortta eksik sütunlar: " + ", ".join(missing))
            writer = csv.DictWriter(target, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                variation_id = row["VariationID"]
                if not variation_id or variation_id in seen_ids:
                    raise ValueError(f"Dış kohort VariationID tekil değil: {variation_id}")
                seen_ids.add(variation_id)
                if row["target_label"] not in {"Benign", "Pathogenic"}:
                    raise ValueError("Dış kohortta yalnızca kesin ikili etiket olmalıdır.")
                expected_target = "1" if row["target_label"] == "Pathogenic" else "0"
                if row["target"] != expected_target:
                    raise ValueError(f"Dış kohort hedef eşleşmiyor: {variation_id}")
                baseline = derive_features(row)
                hgvs = derive_hgvs_features(row["Name"])
                missense = derive_missense_features(row["Name"], hgvs["protein_effect"])
                output = {**baseline, **hgvs, **missense}
                writer.writerow({column: output[column] for column in OUTPUT_COLUMNS})
                counts["rows"] += 1
                targets[row["target"]] += 1
                if hgvs["protein_effect"] == "missense":
                    counts["missense_rows"] += 1
                if hgvs["has_protein_annotation"] == "0":
                    counts["no_protein_annotation_rows"] += 1
                if progress_every and counts["rows"] % progress_every == 0:
                    print(f"İlerleme: {counts['rows']:,} dış satır işlendi.", flush=True)
        if targets.keys() != {"0", "1"}:
            raise ValueError("Dış kohortta iki hedef sınıf da bulunmalıdır.")
        temp_path.replace(output_path)
    except BaseException:
        temp_path.unlink(missing_ok=True)
        raise

    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "input": {"path": str(input_path.resolve()), "sha256": sha256_file(input_path)},
        "output": {"path": str(output_path.resolve()), "sha256": sha256_file(output_path)},
        "counts": {**dict(counts), "benign": targets["0"], "pathogenic": targets["1"]},
        "feature_policy": {
            "raw_name_in_output": False,
            "gene_or_transcript_identity_in_output": False,
            "coordinates_in_output": False,
            "target_used_to_derive_features": False,
            "columns": list(OUTPUT_COLUMNS),
        },
    }
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print("Dış kohort özellik tablosu oluşturuldu")
    print(f"Satır: {counts['rows']:,}; missense: {counts['missense_rows']:,}")
    print(f"Çıktı: {output_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Dış ikili kohorttan model özellikleri üretir.")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--progress-every", type=int, default=100_000)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    build_external_features(args.input, args.output, args.report_dir, progress_every=args.progress_every)


if __name__ == "__main__":
    main()
