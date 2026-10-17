"""Temiz ikili ClinVar alt kümesinden sızıntısı düşük temel özellikler üret."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence


DEFAULT_INPUT = Path(
    "data/processed/clinvar_grch38_germline_exact_binary.tsv.gz"
)
DEFAULT_OUTPUT = Path("data/processed/clinvar_baseline_features.tsv.gz")
DEFAULT_REPORT_DIR = Path("reports/generated/baseline_features")

INPUT_COLUMNS = (
    "VariationID",
    "Type",
    "Chromosome",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
    "target_label",
    "target",
)
FEATURE_COLUMNS = (
    "variant_type",
    "chromosome",
    "reference_length",
    "alternate_length",
    "length_delta",
    "absolute_length_delta",
    "vcf_alleles_missing",
    "is_canonical_snv",
    "substitution",
    "mutation_class",
)
OUTPUT_COLUMNS = ("VariationID", *FEATURE_COLUMNS, "target")
MISSING_MARKERS = {"", "-", "na"}
CANONICAL_BASES = {"A", "C", "G", "T"}
TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def allele_is_missing(value: str) -> bool:
    return value.casefold() in MISSING_MARKERS


def derive_features(row: dict[str, str]) -> dict[str, str]:
    reference = row["ReferenceAlleleVCF"].upper()
    alternate = row["AlternateAlleleVCF"].upper()
    alleles_missing = allele_is_missing(reference) or allele_is_missing(alternate)

    if alleles_missing:
        reference_length = ""
        alternate_length = ""
        length_delta = ""
        absolute_length_delta = ""
    else:
        ref_length_value = len(reference)
        alt_length_value = len(alternate)
        delta_value = alt_length_value - ref_length_value
        reference_length = str(ref_length_value)
        alternate_length = str(alt_length_value)
        length_delta = str(delta_value)
        absolute_length_delta = str(abs(delta_value))

    canonical_snv = (
        not alleles_missing
        and len(reference) == 1
        and len(alternate) == 1
        and reference in CANONICAL_BASES
        and alternate in CANONICAL_BASES
        and reference != alternate
    )
    if canonical_snv:
        substitution = f"{reference}>{alternate}"
        mutation_class = (
            "transition" if (reference, alternate) in TRANSITIONS else "transversion"
        )
    else:
        substitution = "NOT_CANONICAL_SNV"
        mutation_class = "not_canonical_snv"

    return {
        "VariationID": row["VariationID"],
        "variant_type": row["Type"],
        "chromosome": row["Chromosome"],
        "reference_length": reference_length,
        "alternate_length": alternate_length,
        "length_delta": length_delta,
        "absolute_length_delta": absolute_length_delta,
        "vcf_alleles_missing": "1" if alleles_missing else "0",
        "is_canonical_snv": "1" if canonical_snv else "0",
        "substitution": substitution,
        "mutation_class": mutation_class,
        "target": row["target"],
    }


def write_distribution_report(
    report_path: Path,
    distributions: Counter[tuple[str, str, str]],
) -> None:
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["feature", "target", "value", "row_count"])
        for (feature, target, value), count in sorted(
            distributions.items(),
            key=lambda item: (item[0][0], item[0][1], -item[1], item[0][2]),
        ):
            writer.writerow([feature, target, value, count])


def build_features(
    input_path: Path,
    output_path: Path,
    report_dir: Path,
    *,
    progress_every: int = 100_000,
) -> dict[str, object]:
    if not input_path.is_file():
        raise FileNotFoundError(f"Girdi dosyası bulunamadı: {input_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    temporary_output_handle = tempfile.NamedTemporaryFile(
        prefix="baseline_features_",
        suffix=".tsv.gz.tmp",
        dir=output_path.parent,
        delete=False,
    )
    temporary_output_path = Path(temporary_output_handle.name)
    temporary_output_handle.close()

    counters = Counter()
    target_counts: Counter[str] = Counter()
    distributions: Counter[tuple[str, str, str]] = Counter()
    numeric_min: dict[str, int] = {}
    numeric_max: dict[str, int] = {}

    try:
        with (
            gzip.open(input_path, "rt", encoding="utf-8", newline="") as source,
            gzip.open(temporary_output_path, "wt", encoding="utf-8", newline="") as target,
        ):
            reader = csv.DictReader(source, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError("Girdi dosyasında başlık bulunamadı.")
            missing = sorted(set(INPUT_COLUMNS).difference(reader.fieldnames))
            if missing:
                raise ValueError("Gerekli sütunlar bulunamadı: " + ", ".join(missing))

            writer = csv.DictWriter(
                target,
                fieldnames=OUTPUT_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()

            for row in reader:
                counters["input_rows"] += 1
                if row["target_label"] not in {"Benign", "Pathogenic"}:
                    raise ValueError("Beklenmeyen target_label değeri bulundu.")
                expected_target = "1" if row["target_label"] == "Pathogenic" else "0"
                if row["target"] != expected_target:
                    raise ValueError("target ile target_label eşleşmiyor.")

                features = derive_features(row)
                writer.writerow(features)
                counters["output_rows"] += 1
                target_counts[features["target"]] += 1
                if features["vcf_alleles_missing"] == "1":
                    counters["rows_with_missing_vcf_alleles"] += 1
                if features["is_canonical_snv"] == "1":
                    counters["canonical_snv_rows"] += 1

                for feature in (
                    "variant_type",
                    "chromosome",
                    "vcf_alleles_missing",
                    "is_canonical_snv",
                    "substitution",
                    "mutation_class",
                ):
                    distributions[(feature, features["target"], features[feature])] += 1

                for feature in (
                    "reference_length",
                    "alternate_length",
                    "length_delta",
                    "absolute_length_delta",
                ):
                    value = features[feature]
                    if value == "":
                        continue
                    integer_value = int(value)
                    numeric_min[feature] = min(numeric_min.get(feature, math.inf), integer_value)
                    numeric_max[feature] = max(numeric_max.get(feature, -math.inf), integer_value)

                if progress_every > 0 and counters["input_rows"] % progress_every == 0:
                    print(f"İlerleme: {counters['input_rows']:,} satır işlendi.", flush=True)

        temporary_output_path.replace(output_path)
    except BaseException:
        temporary_output_path.unlink(missing_ok=True)
        raise

    write_distribution_report(
        report_dir / "categorical_distributions_by_target.csv",
        distributions,
    )
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "path": str(input_path.resolve()),
            "size_bytes": input_path.stat().st_size,
            "sha256": sha256_file(input_path),
        },
        "output": {
            "path": str(output_path.resolve()),
            "size_bytes": output_path.stat().st_size,
            "sha256": sha256_file(output_path),
            "columns": list(OUTPUT_COLUMNS),
        },
        "counts": {
            **dict(counters),
            "target_0_benign": target_counts["0"],
            "target_1_pathogenic": target_counts["1"],
        },
        "numeric_ranges": {
            feature: {"min": numeric_min[feature], "max": numeric_max[feature]}
            for feature in numeric_min
        },
        "feature_policy": {
            "model_feature_candidates": list(FEATURE_COLUMNS),
            "identity_only": ["VariationID"],
            "target": ["target"],
            "excluded_from_baseline": [
                "#AlleleID",
                "Name",
                "GeneID",
                "GeneSymbol",
                "HGNC_ID",
                "ChromosomeAccession",
                "Start",
                "Stop",
                "PositionVCF",
            ],
            "missing_length_representation": "empty TSV value",
            "non_snv_substitution_value": "NOT_CANONICAL_SNV",
        },
    }
    with (report_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("\nTemel özellik tablosu oluşturuldu")
    print(f"Satır: {counters['output_rows']:,}")
    print(f"Kanonik SNV: {counters.get('canonical_snv_rows', 0):,}")
    print(
        "VCF alleli eksik satır: "
        f"{counters.get('rows_with_missing_vcf_alleles', 0):,}"
    )
    print(f"Çıktı: {output_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Temiz ClinVar alt kümesinden temel özellik tablosu üretir."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--progress-every", type=int, default=100_000)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    build_features(
        args.input,
        args.output,
        args.report_dir,
        progress_every=args.progress_every,
    )


if __name__ == "__main__":
    main()

