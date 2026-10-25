"""Dış kohortu geliştirme verisiyle çakışma denetiminden geçirerek kabul et."""

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
    from scripts.build_baseline_features import MISSING_MARKERS, sha256_file
    from scripts.variant_normalization import InvalidVariant, ReferenceGenome
except ModuleNotFoundError:  # Doğrudan `python scripts/...` çalıştırması için.
    from build_baseline_features import MISSING_MARKERS, sha256_file  # type: ignore[no-redef]
    from variant_normalization import InvalidVariant, ReferenceGenome


DEFAULT_DEVELOPMENT = Path("data/processed/clinvar_grch38_germline_exact_binary.tsv.gz")
DEFAULT_EXTERNAL = Path("data/external/candidate_binary.tsv.gz")
DEFAULT_OUTPUT = Path("data/external/admitted_binary.tsv.gz")
DEFAULT_REPORT_DIR = Path("reports/generated/external_admission")
REQUIRED_COLUMNS = {
    "VariationID", "Type", "Name", "Assembly", "Chromosome", "PositionVCF",
    "ReferenceAlleleVCF", "AlternateAlleleVCF", "target_label", "target",
}
OPTIONAL_ALLELE_ID = "#AlleleID"
OPTIONAL_GENE = "GeneSymbol"


def gene_tokens(value: str) -> set[str]:
    return {token.strip() for token in value.replace(",", ";").split(";")
            if token.strip().casefold() not in MISSING_MARKERS | {"."}}


def validate_label(row: dict[str, str], *, source_name: str) -> None:
    label = row["target_label"]
    target = row["target"]
    if label not in {"Benign", "Pathogenic"}:
        raise ValueError(f"{source_name} kesin ikili etiket içermiyor: {label}")
    expected = "1" if label == "Pathogenic" else "0"
    if target != expected:
        raise ValueError(f"{source_name} hedef eşleşmiyor: {row['VariationID']}")


def load_development_index(path: Path, reference: ReferenceGenome):
    identifiers: set[str] = set()
    allele_ids: set[str] = set()
    coordinates: set[tuple[str, str, str, str]] = set()
    genes: set[str] = set()
    unnormalizable = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Geliştirme kaynağında başlık bulunamadı.")
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames)
        if missing:
            raise ValueError("Geliştirme kaynağında eksik sütunlar: " + ", ".join(sorted(missing)))
        for row in reader:
            variation_id = row["VariationID"]
            if not variation_id or variation_id in identifiers:
                raise ValueError(f"Geliştirme VariationID tekil değil: {variation_id}")
            identifiers.add(variation_id)
            if row.get(OPTIONAL_ALLELE_ID, "").casefold() not in MISSING_MARKERS | {"."}:
                allele_ids.add(row[OPTIONAL_ALLELE_ID])
            genes.update(gene_tokens(row.get(OPTIONAL_GENE, "")))
            try:
                coordinates.add(reference.normalize(row))
            except InvalidVariant as error:
                if str(error) not in {"unsupported_or_missing_alleles", "invalid_position", "coordinate_out_of_range"}:
                    raise ValueError(f"Geliştirme referans kontrolü başarısız: {variation_id}: {error}") from error
                unnormalizable += 1
    return identifiers, allele_ids, coordinates, genes, unnormalizable


def admit_external_cohort(
    development_path: Path,
    external_path: Path,
    output_path: Path,
    report_dir: Path,
    *,
    cohort_kind: str,
    source_version: str,
    source_url: str,
    reference_config: Path,
) -> dict[str, object]:
    with ReferenceGenome(reference_config) as reference:
        return _admit_external_cohort(
            development_path, external_path, output_path, report_dir,
            cohort_kind=cohort_kind, source_version=source_version,
            source_url=source_url, reference=reference,
        )


def _admit_external_cohort(
    development_path: Path, external_path: Path, output_path: Path,
    report_dir: Path, *, cohort_kind: str, source_version: str,
    source_url: str, reference: ReferenceGenome,
) -> dict[str, object]:
    if cohort_kind not in {"independent", "temporal_clinvar"}:
        raise ValueError("cohort_kind independent veya temporal_clinvar olmalıdır.")
    if not source_version.strip() or not source_url.strip():
        raise ValueError("Dış kohort sürümü ve kaynak adresi zorunludur.")
    for path in (development_path, external_path):
        if not path.is_file():
            raise FileNotFoundError(f"Girdi bulunamadı: {path}")
    if output_path.resolve() in {development_path.resolve(), external_path.resolve()}:
        raise ValueError("Kabul çıktısı girdi dosyasının üzerine yazılamaz.")
    development_ids, development_alleles, development_coordinates, development_genes, unnormalizable = load_development_index(development_path, reference)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    temporary_handle = tempfile.NamedTemporaryFile(
        prefix="admitted_external_", suffix=".tsv.gz.tmp", dir=output_path.parent, delete=False
    )
    temporary_path = Path(temporary_handle.name)
    temporary_handle.close()
    counters: Counter[str] = Counter()
    targets: Counter[str] = Counter()
    seen_external_ids: set[str] = set()
    overlap_examples: list[dict[str, str]] = []
    external_genes: set[str] = set()
    annotations: list[dict[str, str]] = []
    seen_coordinates: set[tuple[str, str, str, str]] = set()

    try:
        with (
            gzip.open(external_path, "rt", encoding="utf-8", newline="") as source,
            gzip.open(temporary_path, "wt", encoding="utf-8", newline="") as target,
        ):
            reader = csv.DictReader(source, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError("Dış kohortta başlık bulunamadı.")
            missing = REQUIRED_COLUMNS.difference(reader.fieldnames)
            if missing:
                raise ValueError("Dış kohortta eksik sütunlar: " + ", ".join(sorted(missing)))
            writer = csv.DictWriter(target, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                counters["external_rows"] += 1
                validate_label(row, source_name="Dış kohort")
                variation_id = row["VariationID"]
                if not variation_id or variation_id in seen_external_ids:
                    raise ValueError(f"Dış kohort VariationID tekil değil: {variation_id}")
                seen_external_ids.add(variation_id)
                reasons: list[str] = []
                if variation_id in development_ids:
                    reasons.append("variation_id")
                if (
                    OPTIONAL_ALLELE_ID in reader.fieldnames
                    and row[OPTIONAL_ALLELE_ID].casefold() not in MISSING_MARKERS | {"."}
                    and row[OPTIONAL_ALLELE_ID] in development_alleles
                ):
                    reasons.append("allele_id")
                try:
                    coordinate = reference.normalize(row)
                    if coordinate in development_coordinates:
                        reasons.append("coordinate_allele")
                    if coordinate in seen_coordinates:
                        raise ValueError(f"Dış kohortta aynı normalleştirilmiş varyant tekrarlandı: {variation_id}")
                except InvalidVariant as error:
                    reasons.append(str(error))
                if reasons:
                    counters["excluded_rows"] += 1
                    for reason in reasons:
                        counters[f"excluded_{reason}"] += 1
                    if len(overlap_examples) < 100:
                        overlap_examples.append({"VariationID": variation_id, "exclusion_reasons": ";".join(reasons)})
                    continue
                # Model kromozom kategorileri ClinVar biçimindedir (1, X, MT).
                # Allel uzunluk tarifini koru; normalize anahtar yalnızca eşleştirme içindir.
                row["Chromosome"] = coordinate[0]
                writer.writerow(row)
                seen_coordinates.add(coordinate)
                counters["admitted_rows"] += 1
                targets[row["target"]] += 1
                genes = gene_tokens(row.get(OPTIONAL_GENE, ""))
                external_genes.update(genes)
                annotations.append({
                    "VariationID": variation_id,
                    "unseen_gene": "1" if genes and not genes.intersection(development_genes) else "0",
                    "gene_status": "missing" if not genes else ("unseen" if not genes.intersection(development_genes) else "seen"),
                })
        if targets.keys() != {"0", "1"}:
            raise ValueError("Kabul edilen dış kohort iki hedef sınıfı da içermelidir.")
        temporary_path.replace(output_path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise

    examples_path = report_dir / "excluded_overlap_examples.csv"
    with examples_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["VariationID", "exclusion_reasons"])
        writer.writeheader()
        writer.writerows(overlap_examples)
    annotation_path = report_dir / "admitted_annotations.tsv"
    with annotation_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["VariationID", "unseen_gene", "gene_status"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(annotations)
    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "admission_decision": "accepted",
        "schema_version": 2,
        "reference": reference.metadata,
        "cohort_kind": cohort_kind,
        "source_metadata": {"version": source_version, "url": source_url},
        "inputs": {
            "development": {"path": str(development_path.resolve()), "sha256": sha256_file(development_path)},
            "external_candidate": {"path": str(external_path.resolve()), "sha256": sha256_file(external_path)},
        },
        "output": {"path": str(output_path.resolve()), "sha256": sha256_file(output_path)},
        "counts": {**dict(counters), "benign": targets["0"], "pathogenic": targets["1"], "external_gene_symbols_reported": len(external_genes), "development_rows_without_normalizable_coordinates": unnormalizable},
        "overlap_policy": {
            "keys": ["VariationID", "#AlleleID", "GRCh38_coordinate_reference_alternate"],
            "missing_coordinate_allele_rows": "excluded",
            "development_overlap_rows": "excluded",
            "normalization": "reference_checked_biallelic_left_alignment_and_trimming_v1",
            "unsupported_complex_or_symbolic_alleles": "excluded",
            "normalized_duplicate_external_variants": "abort",
            "development_scope": "entire_previously_used_dataset_including_internal_test",
        },
        "annotations": {"path": str(annotation_path.resolve()), "sha256": sha256_file(annotation_path)},
        "outputs": {"excluded_overlap_examples": str(examples_path.resolve())},
    }
    with (report_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print("Dış kohort çakışma denetiminden geçti")
    print(f"Kabul edilen: {counters['admitted_rows']:,}; dışlanan: {counters['excluded_rows']:,}")
    print(f"Rapor: {(report_dir / 'manifest.json').resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Dış kohortu çakışmasız değerlendirme için kabul eder.")
    parser.add_argument("--development", type=Path, default=DEFAULT_DEVELOPMENT)
    parser.add_argument("--external", type=Path, default=DEFAULT_EXTERNAL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--cohort-kind", choices=("independent", "temporal_clinvar"), required=True)
    parser.add_argument("--source-version", required=True)
    parser.add_argument("--source-url", required=True)
    parser.add_argument("--reference-config", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    admit_external_cohort(args.development, args.external, args.output, args.report_dir, cohort_kind=args.cohort_kind, source_version=args.source_version, source_url=args.source_url, reference_config=args.reference_config)


if __name__ == "__main__":
    main()
