"""eRepo adaylarını sabitlenmiş ClinVar GRCh38 allellerine eşle; etiket kökenini koru."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

try:
    from scripts.build_baseline_features import sha256_file
    from scripts.build_initial_binary_dataset import OUTPUT_SOURCE_COLUMNS
except ModuleNotFoundError:
    from build_baseline_features import sha256_file
    from build_initial_binary_dataset import OUTPUT_SOURCE_COLUMNS


DEFAULT_SCREEN = Path("reports/generated/clingen_erepo_screen/manifest.json")
DEFAULT_RAW = Path("data/raw/variant_summary.txt.gz")
DEFAULT_OUTPUT = Path("data/external/clingen_mapped_candidates.tsv.gz")
DEFAULT_REPORT = Path("reports/generated/clingen_candidate_mapping")
PROVENANCE_COLUMNS = (
    "annotation_source", "clingen_curation_uuids_json", "clingen_conditions_json",
    "clingen_panels_json", "clingen_doc_versions_json", "original_clinvar_labels_json",
)
OUTPUT_COLUMNS = (*OUTPUT_SOURCE_COLUMNS, "target_label", "target", *PROVENANCE_COLUMNS)


def map_candidates(screen_path: Path, raw_path: Path, output_path: Path, report_dir: Path) -> dict:
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    if screen.get("decision") != "screening_only_not_admitted_for_external_evaluation":
        raise ValueError("Beklenen eRepo tarama manifesti bulunamadı.")
    candidate_info = screen["outputs"]["curation_candidates"]
    candidates_path = Path(candidate_info["path"])
    if sha256_file(candidates_path) != candidate_info["sha256"]:
        raise ValueError("eRepo aday tablosu hash'i uyuşmuyor.")
    groups: dict[str, list[dict]] = defaultdict(list)
    with candidates_path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            groups[row["ClinVarVariationID"]].append(row)
    for variation_id, curations in groups.items():
        labels = {row["classification"] for row in curations}
        if len(labels) != 1 or not labels.issubset({"Pathogenic", "Benign"}):
            raise ValueError(f"Adayın tek kesin etiketi yok: {variation_id}")
    if len(groups) != screen["counts"]["id_unmatched_candidate_unique_variants"]:
        raise ValueError("Aday varyant sayısı manifestle uyuşmuyor.")
    if output_path.resolve() in {raw_path.resolve(), candidates_path.resolve()}:
        raise ValueError("Çıktı girdinin üzerine yazılamaz.")
    if output_path.exists():
        raise FileExistsError(f"Mevcut eşleme çıktısı üzerine yazılmayacak: {output_path}")
    representations: dict[str, dict[tuple, dict]] = defaultdict(dict)
    old_labels: dict[str, set[str]] = defaultdict(set)
    counts: Counter[str] = Counter()
    required = set(OUTPUT_SOURCE_COLUMNS) | {"ClinicalSignificance"}
    csv.field_size_limit(10_000_000)
    with gzip.open(raw_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if required - set(header):
            raise ValueError("ClinVar ham eşleme kaynağında gerekli sütunlar eksik.")
        ix = {name: index for index, name in enumerate(header)}
        for values in reader:
            counts["raw_rows_scanned"] += 1
            if len(values) != len(header):
                raise ValueError("ClinVar kaynağında bozuk sütun sayısı.")
            variation_id = values[ix["VariationID"]]
            if variation_id in groups and values[ix["Assembly"]] == "GRCh38" and values[ix["OriginSimple"]] == "germline":
                mapped = {name: values[ix[name]] for name in OUTPUT_SOURCE_COLUMNS}
                key = tuple(mapped[name] for name in ("#AlleleID", "ChromosomeAccession", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF", "Name", "GeneSymbol", "Type"))
                representations[variation_id][key] = mapped
                old_labels[variation_id].add(values[ix["ClinicalSignificance"]])
                counts["candidate_grch38_germline_source_rows"] += 1
            if counts["raw_rows_scanned"] % 1_000_000 == 0:
                print(f"ClinVar eşleme taraması: {counts['raw_rows_scanned']:,} satır", flush=True)
    output_rows = []
    audit_rows = []
    for variation_id, curations in sorted(groups.items(), key=lambda item: int(item[0])):
        matches = representations.get(variation_id, {})
        label = curations[0]["classification"]  # Etiket tekilliği yukarıda doğrulandı.
        if not matches:
            status = "not_found_in_pinned_grch38_germline_source"
        elif len(matches) != 1:
            status = "ambiguous_multiple_source_representations"
        else:
            row = dict(next(iter(matches.values())))
            row.update({
                "target_label": label, "target": "1" if label == "Pathogenic" else "0",
                "annotation_source": "pinned_clinvar_raw_Name_and_genomic_alleles",
                "clingen_curation_uuids_json": json.dumps(sorted({curation["uuid"] for curation in curations})),
                "clingen_conditions_json": json.dumps(sorted({curation["condition"] for curation in curations})),
                "clingen_panels_json": json.dumps(sorted({curation["expert_panel"] for curation in curations})),
                "clingen_doc_versions_json": json.dumps(sorted({curation["doc_version"] for curation in curations})),
                "original_clinvar_labels_json": json.dumps(sorted(old_labels[variation_id])),
            })
            output_rows.append(row)
            status = "mapped_pending_reference_and_overlap_checks"
            counts[f"mapped_{label.lower()}"] += 1
        counts[status] += 1
        audit_rows.append({
            "VariationID": variation_id, "status": status, "distinct_source_representations": len(matches),
            "clingen_label": label, "original_clinvar_labels_json": json.dumps(sorted(old_labels[variation_id])),
        })
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    audit_path = report_dir / "mapping_audit.tsv"
    with audit_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["VariationID", "status", "distinct_source_representations", "clingen_label", "original_clinvar_labels_json"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(audit_rows)
    counts["input_candidate_variants"] = len(groups)
    counts["output_mapped_variants"] = len(output_rows)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "status": "mapped_not_yet_admitted",
        "inputs": {"screen_manifest": {"path": str(screen_path.resolve()), "sha256": sha256_file(screen_path)},
                   "raw_clinvar": {"path": str(raw_path.resolve()), "sha256": sha256_file(raw_path)}},
        "output": {"path": str(output_path.resolve()), "sha256": sha256_file(output_path)},
        "counts": dict(counts),
        "provenance": {"target_source": "ClinGen_eRepo_expert_curation", "annotation_and_allele_source": "pinned_local_ClinVar_raw_snapshot", "source_independence_verified": False,
                       "cohort_kind": "expert_curated_supplementary", "model_loaded": False},
        "outputs": {"mapping_audit": str(audit_path.resolve())},
    }
    (report_dir / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest["counts"], indent=2))
    return manifest


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description="eRepo adaylarının ClinVar GRCh38 allel köprüsünü kurar.")
    parser.add_argument("--screen-manifest", type=Path, default=DEFAULT_SCREEN)
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    map_candidates(args.screen_manifest, args.raw, args.output, args.report_dir)


if __name__ == "__main__":
    main()
