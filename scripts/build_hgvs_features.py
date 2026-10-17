"""ClinVar Name alanından kimliksiz, kategorik HGVS özellikleri üret."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import re
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


DEFAULT_SOURCE = Path(
    "data/processed/clinvar_grch38_germline_exact_binary.tsv.gz"
)
DEFAULT_BASELINE = Path("data/processed/clinvar_baseline_features.tsv.gz")
DEFAULT_SPLIT = Path("data/processed/baseline_gene_disjoint_split.tsv")
DEFAULT_OUTPUT = Path("data/processed/clinvar_hgvs_features.tsv.gz")
DEFAULT_REPORT_DIR = Path("reports/generated/hgvs_features")

BASELINE_COLUMNS = (
    "VariationID",
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
HGVS_FEATURE_COLUMNS = (
    "hgvs_sequence_type",
    "dna_change_type",
    "coding_region",
    "splice_offset_category",
    "protein_effect",
    "has_protein_annotation",
)
OUTPUT_COLUMNS = (*BASELINE_COLUMNS, *HGVS_FEATURE_COLUMNS, "target")
DNA_NOTATION = re.compile(r":([cgnmr])\.([^\s]+)", re.IGNORECASE)
PROTEIN_NOTATION = re.compile(r"\(p\.([^)]+)\)", re.IGNORECASE)
INTRON_OFFSET = re.compile(r"(?<=\d)([+-])(\d+)")
STANDARD_AMINO_ACID = (
    r"(?:Ala|Arg|Asn|Asp|Cys|Gln|Glu|Gly|His|Ile|Leu|Lys|Met|Phe|Pro|"
    r"Ser|Thr|Trp|Tyr|Val)"
)
MISSENSE_PROTEIN = re.compile(
    rf"^{STANDARD_AMINO_ACID}\d+{STANDARD_AMINO_ACID}$"
)
SEQUENCE_TYPE_NAMES = {
    "c": "coding_transcript",
    "g": "genomic",
    "n": "noncoding_transcript",
    "m": "mitochondrial",
    "r": "rna",
}


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def classify_dna_change(dna_change: str | None) -> str:
    if dna_change is None:
        return "no_hgvs_dna_notation"
    lowered = dna_change.casefold()
    if "delins" in lowered:
        return "delins"
    if "dup" in lowered:
        return "duplication"
    if "del" in lowered:
        return "deletion"
    if "ins" in lowered:
        return "insertion"
    if "inv" in lowered:
        return "inversion"
    if ">" in dna_change:
        return "substitution"
    if "=" in dna_change:
        return "no_sequence_change"
    return "other"


def classify_coding_region(prefix: str | None, dna_change: str | None) -> str:
    if prefix != "c" or dna_change is None:
        return "not_coding_transcript"
    if dna_change.startswith("-"):
        return "five_prime_utr"
    if dna_change.startswith("*"):
        return "three_prime_utr"
    if INTRON_OFFSET.search(dna_change):
        return "intronic"
    return "coding_exonic"


def classify_splice_offset(prefix: str | None, dna_change: str | None) -> str:
    if prefix != "c" or dna_change is None:
        return "not_coding_transcript"
    offsets = [(sign, int(distance)) for sign, distance in INTRON_OFFSET.findall(dna_change)]
    if not offsets:
        return "not_intronic"
    if any(sign == "+" and distance in {1, 2} for sign, distance in offsets):
        return "donor_core_1_2"
    if any(sign == "-" and distance in {1, 2} for sign, distance in offsets):
        return "acceptor_core_1_2"
    if min(distance for _, distance in offsets) <= 8:
        return "splice_region_3_8"
    return "deep_intronic"


def classify_protein_effect(protein_change: str | None) -> str:
    if protein_change is None:
        return "no_protein_annotation"
    normalized = protein_change.strip()
    lowered = normalized.casefold()
    if normalized.endswith("="):
        return "synonymous"
    if "fs" in lowered:
        return "frameshift"
    if normalized.startswith("Met1"):
        return "start_lost_or_uncertain"
    if (lowered.startswith("ter") or normalized.startswith("*")) and "ext" in lowered:
        return "stop_lost"
    if "ter" in lowered or "*" in normalized:
        return "stop_gained_or_changed"
    if "delins" in lowered:
        return "inframe_delins"
    if "dup" in lowered:
        return "inframe_duplication"
    if "del" in lowered:
        return "inframe_deletion"
    if "ins" in lowered:
        return "inframe_insertion"
    if "?" in normalized or normalized == "0" or "xaa" in lowered:
        return "uncertain_or_absent"
    if MISSENSE_PROTEIN.fullmatch(normalized):
        return "missense"
    return "other_protein_effect"


def derive_hgvs_features(name: str) -> dict[str, str]:
    dna_match = DNA_NOTATION.search(name)
    if dna_match:
        prefix = dna_match.group(1).casefold()
        dna_change = dna_match.group(2)
        sequence_type = SEQUENCE_TYPE_NAMES[prefix]
    else:
        prefix = None
        dna_change = None
        sequence_type = "unknown"

    protein_match = PROTEIN_NOTATION.search(name)
    protein_change = protein_match.group(1) if protein_match else None
    return {
        "hgvs_sequence_type": sequence_type,
        "dna_change_type": classify_dna_change(dna_change),
        "coding_region": classify_coding_region(prefix, dna_change),
        "splice_offset_category": classify_splice_offset(prefix, dna_change),
        "protein_effect": classify_protein_effect(protein_change),
        "has_protein_annotation": "1" if protein_change is not None else "0",
    }


def write_distribution_report(
    path: Path,
    distributions: Counter[tuple[str, str, str, str]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["feature", "split", "target", "value", "row_count"])
        for (feature, split_name, target, value), count in sorted(
            distributions.items(),
            key=lambda item: (
                item[0][0],
                item[0][1],
                item[0][2],
                -item[1],
                item[0][3],
            ),
        ):
            writer.writerow([feature, split_name, target, value, count])


def load_split_assignments(split_path: Path) -> dict[str, str]:
    if not split_path.is_file():
        raise FileNotFoundError(f"Ayrım dosyası bulunamadı: {split_path}")
    assignments: dict[str, str] = {}
    with split_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Ayrım dosyasında başlık bulunamadı.")
        missing = sorted({"VariationID", "split"}.difference(reader.fieldnames))
        if missing:
            raise ValueError("Ayrım dosyasında eksik sütunlar: " + ", ".join(missing))
        for row in reader:
            variation_id = row["VariationID"]
            split_name = row["split"]
            if variation_id in assignments:
                raise ValueError(f"Ayrım dosyasında tekrarlanan VariationID: {variation_id}")
            if split_name not in {"train", "validation", "test"}:
                raise ValueError(f"Beklenmeyen ayrım değeri: {split_name}")
            assignments[variation_id] = split_name
    return assignments


def build_hgvs_features(
    source_path: Path,
    baseline_path: Path,
    split_path: Path,
    output_path: Path,
    report_dir: Path,
    *,
    progress_every: int = 100_000,
) -> dict[str, object]:
    for path in (source_path, baseline_path, split_path):
        if not path.is_file():
            raise FileNotFoundError(f"Girdi dosyası bulunamadı: {path}")

    split_assignments = load_split_assignments(split_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    temporary_handle = tempfile.NamedTemporaryFile(
        prefix="hgvs_features_",
        suffix=".tsv.gz.tmp",
        dir=output_path.parent,
        delete=False,
    )
    temporary_path = Path(temporary_handle.name)
    temporary_handle.close()
    distributions: Counter[tuple[str, str, str, str]] = Counter()
    counters: Counter[str] = Counter()
    seen_ids: set[str] = set()

    try:
        with (
            gzip.open(source_path, "rt", encoding="utf-8", newline="") as source_handle,
            gzip.open(baseline_path, "rt", encoding="utf-8", newline="") as baseline_handle,
            gzip.open(temporary_path, "wt", encoding="utf-8", newline="") as output_handle,
        ):
            source_reader = csv.DictReader(source_handle, delimiter="\t")
            baseline_reader = csv.DictReader(baseline_handle, delimiter="\t")
            if source_reader.fieldnames is None or baseline_reader.fieldnames is None:
                raise ValueError("Girdi dosyalarından birinde başlık bulunamadı.")
            source_missing = sorted(
                {"VariationID", "Name", "target"}.difference(source_reader.fieldnames)
            )
            baseline_missing = sorted(
                set((*BASELINE_COLUMNS, "target")).difference(baseline_reader.fieldnames)
            )
            if source_missing:
                raise ValueError("Kaynakta eksik sütunlar: " + ", ".join(source_missing))
            if baseline_missing:
                raise ValueError(
                    "Temel özelliklerde eksik sütunlar: " + ", ".join(baseline_missing)
                )

            writer = csv.DictWriter(
                output_handle,
                fieldnames=OUTPUT_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for source_row, baseline_row in itertools.zip_longest(
                source_reader,
                baseline_reader,
            ):
                if source_row is None or baseline_row is None:
                    raise ValueError("Kaynak ve temel özellik satır sayıları eşleşmiyor.")
                if source_row["VariationID"] != baseline_row["VariationID"]:
                    raise ValueError(
                        "Kaynak ve temel özellik sırası eşleşmiyor: "
                        f"{source_row['VariationID']} != {baseline_row['VariationID']}"
                    )
                if source_row["target"] != baseline_row["target"]:
                    raise ValueError(
                        f"Hedef uyuşmazlığı: VariationID={source_row['VariationID']}"
                    )
                variation_id = source_row["VariationID"]
                split_name = split_assignments.get(variation_id)
                if split_name is None:
                    raise ValueError(f"Ayrım ataması olmayan VariationID: {variation_id}")
                seen_ids.add(variation_id)

                hgvs_features = derive_hgvs_features(source_row["Name"])
                output_row = {
                    column: baseline_row[column] for column in BASELINE_COLUMNS
                }
                output_row.update(hgvs_features)
                output_row["target"] = baseline_row["target"]
                writer.writerow(output_row)
                counters["rows"] += 1
                counters[f"{split_name}_rows"] += 1
                if hgvs_features["has_protein_annotation"] == "1":
                    counters["rows_with_protein_annotation"] += 1
                if hgvs_features["hgvs_sequence_type"] == "unknown":
                    counters["rows_without_recognized_hgvs_dna"] += 1
                if split_name in {"train", "validation"}:
                    for feature in HGVS_FEATURE_COLUMNS:
                        distributions[
                            (
                                feature,
                                split_name,
                                baseline_row["target"],
                                hgvs_features[feature],
                            )
                        ] += 1

                if progress_every > 0 and counters["rows"] % progress_every == 0:
                    print(f"İlerleme: {counters['rows']:,} satır işlendi.", flush=True)
        missing_feature_ids = set(split_assignments).difference(seen_ids)
        if missing_feature_ids:
            raise ValueError(
                "Özellik girdisinde bulunmayan ayrım kimliği: "
                f"{min(missing_feature_ids)}"
            )
        temporary_path.replace(output_path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise

    distribution_path = report_dir / "categorical_distributions_by_development_partition.csv"
    write_distribution_report(distribution_path, distributions)
    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "clean_dataset": {
                "path": str(source_path.resolve()),
                "sha256": sha256_file(source_path),
            },
            "baseline_features": {
                "path": str(baseline_path.resolve()),
                "sha256": sha256_file(baseline_path),
            },
            "gene_disjoint_split": {
                "path": str(split_path.resolve()),
                "sha256": sha256_file(split_path),
            },
        },
        "output": {
            "path": str(output_path.resolve()),
            "size_bytes": output_path.stat().st_size,
            "sha256": sha256_file(output_path),
            "columns": list(OUTPUT_COLUMNS),
        },
        "counts": dict(counters),
        "feature_policy": {
            "hgvs_features": list(HGVS_FEATURE_COLUMNS),
            "raw_name_in_output": False,
            "gene_symbol_in_output": False,
            "transcript_identifier_in_output": False,
            "parsing_scope": (
                "Only HGVS DNA/protein consequence patterns are retained; raw identity "
                "text is discarded."
            ),
            "distribution_report_partitions": ["train", "validation"],
            "test_target_distribution_reported": False,
        },
        "reports": {
            "categorical_distributions_by_development_partition": str(
                distribution_path.resolve()
            )
        },
    }
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("HGVS özellik tablosu oluşturuldu")
    print(f"Satır: {counters['rows']:,}")
    print(
        "Protein açıklaması bulunan satır: "
        f"{counters.get('rows_with_protein_annotation', 0):,}"
    )
    print(
        "Tanınan HGVS DNA gösterimi olmayan satır: "
        f"{counters.get('rows_without_recognized_hgvs_dna', 0):,}"
    )
    print(f"Çıktı: {output_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Temel özellik tablosuna kimliksiz HGVS sonuç özellikleri ekler."
    )
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--split", type=Path, default=DEFAULT_SPLIT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--progress-every", type=int, default=100_000)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    build_hgvs_features(
        args.source,
        args.baseline,
        args.split,
        args.output,
        args.report_dir,
        progress_every=args.progress_every,
    )


if __name__ == "__main__":
    main()
