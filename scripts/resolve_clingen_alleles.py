"""eRepo CAID'lerini ClinGen Registry GRCh38 ana kromozom allellerine çöz."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

try:
    from scripts.build_baseline_features import sha256_file
    from scripts.build_initial_binary_dataset import OUTPUT_SOURCE_COLUMNS
    from scripts.map_clingen_candidates import OUTPUT_COLUMNS
    from scripts.variant_normalization import ReferenceGenome
except ModuleNotFoundError:
    from build_baseline_features import sha256_file
    from build_initial_binary_dataset import OUTPUT_SOURCE_COLUMNS
    from map_clingen_candidates import OUTPUT_COLUMNS
    from variant_normalization import ReferenceGenome


def registry_vcf(document: dict, variation_id: str, reference: ReferenceGenome) -> dict:
    """Registry 0-based half-open allelini referans denetimli VCF gösterimine çevir."""
    linked_ids = {str(item["variationId"]) for item in document.get("externalRecords", {}).get("ClinVarVariations", [])}
    if variation_id not in linked_ids:
        raise ValueError("registry_clinvar_identity_not_confirmed")
    options = {}
    for allele in document.get("genomicAlleles", []):
        if allele.get("referenceGenome") != "GRCh38" or len(allele.get("coordinates", [])) != 1:
            continue
        accessions = {hgvs.split(":", 1)[0] for hgvs in allele.get("hgvs", [])}
        accessions = accessions.intersection(reference.primary_accessions)
        if len(accessions) != 1:
            continue
        accession = next(iter(accessions))
        chromosome = reference.primary_accessions[accession]
        coordinate = allele["coordinates"][0]
        start, end = coordinate["start"], coordinate["end"]
        ref, alt = coordinate["referenceAllele"].upper(), coordinate["allele"].upper()
        if not isinstance(start, int) or not isinstance(end, int) or start < 0 or end - start != len(ref):
            raise ValueError("registry_coordinate_length_mismatch")
        if set(ref + alt) - set("ACGT") or ref == alt:
            raise ValueError("registry_unsupported_allele")
        variant_type = "Deletion" if not alt else "Insertion" if not ref else "single nucleotide variant" if len(ref) == len(alt) == 1 else "Indel"
        position = start + 1
        if not ref or not alt:
            if start:
                anchor = reference.fetch(chromosome, start, 1)
                position = start
                ref, alt = anchor + ref, anchor + alt
            else:
                anchor = reference.fetch(chromosome, end + 1, 1)
                position = 1
                ref, alt = ref + anchor, alt + anchor
        row = {"Assembly": "GRCh38", "Chromosome": chromosome, "ChromosomeAccession": accession,
               "PositionVCF": str(position), "ReferenceAlleleVCF": ref, "AlternateAlleleVCF": alt,
               "Start": str(start + 1), "Stop": str(end), "Type": variant_type}
        key = reference.normalize(row)
        options[key] = row
    if len(options) != 1:
        raise ValueError("registry_missing_or_ambiguous_primary_grch38_mapping")
    return next(iter(options.values()))


def resolve_candidates(screen_path: Path, local_mapped_path: Path, reference_config: Path,
                       cache_dir: Path, output_path: Path, report_dir: Path) -> dict:
    if output_path.exists():
        raise FileExistsError(f"Mevcut çıktı korunuyor: {output_path}")
    screen = json.loads(screen_path.read_text(encoding="utf-8"))
    info = screen["outputs"]["curation_candidates"]
    if sha256_file(Path(info["path"])) != info["sha256"]:
        raise ValueError("Aday yorum tablosu hash uyuşmazlığı.")
    groups = defaultdict(list)
    with Path(info["path"]).open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            groups[row["ClinVarVariationID"]].append(row)
    local = {}
    with gzip.open(local_mapped_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            local[row["VariationID"]] = row
    cache_dir.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    caids = sorted({row["CAID"] for rows in groups.values() for row in rows if row["CAID"].startswith("CA") and row["CAID"][2:].isdigit()})
    def download(caid):
        if not caid.startswith("CA") or not caid[2:].isdigit():
            return caid, "invalid_caid"
        path = cache_dir / f"{caid}.json"
        metadata = cache_dir / f"{caid}.receipt.json"
        if path.exists():
            if not metadata.exists() or json.loads(metadata.read_text())["sha256"] != sha256_file(path):
                raise ValueError(f"Registry cache hash uyuşmazlığı: {caid}")
            return caid, None
        url = f"https://reg.genome.network/allele/{caid}"
        try:
            with urlopen(url, timeout=30) as response:
                data = response.read()
            document = json.loads(data)
            if document.get("@id", "").rsplit("/", 1)[-1] != caid:
                raise ValueError("registry_caid_identity_mismatch")
            path.write_bytes(data)
            metadata.write_text(json.dumps({"url": url, "retrieved_at_utc": datetime.now(timezone.utc).isoformat(), "sha256": sha256_file(path)}, indent=2), encoding="utf-8")
            return caid, None
        except Exception as error:
            return caid, f"{type(error).__name__}: {error}"
    failures = {}
    with ThreadPoolExecutor(max_workers=4) as executor:
        for number, (caid, error) in enumerate(executor.map(download, caids), 1):
            if error:
                failures[caid] = error
            if number % 25 == 0 or number == len(caids):
                print(f"ClinGen Registry: {number}/{len(caids)} CAID", flush=True)
    counts = Counter()
    audit = []
    resolved = []
    with ReferenceGenome(reference_config) as reference:
        for variation_id, curations in sorted(groups.items(), key=lambda item: int(item[0])):
            counts["input_variants"] += 1
            caids_for_variant = {row["CAID"] for row in curations}
            try:
                if len(caids_for_variant) != 1:
                    raise ValueError("ambiguous_caid")
                caid = next(iter(caids_for_variant))
                if not caid.startswith("CA") or not caid[2:].isdigit():
                    raise ValueError("missing_or_invalid_caid")
                if caid in failures:
                    raise ValueError("registry_download_failed")
                document = json.loads((cache_dir / f"{caid}.json").read_text())
                mapping = registry_vcf(document, variation_id, reference)
                labels = {row["classification"] for row in curations}
                if len(labels) != 1 or not labels.issubset({"Benign", "Pathogenic"}):
                    raise ValueError("ambiguous_target")
                label = next(iter(labels))
                names = {row["preferred_variant_title"] for row in curations if row["preferred_variant_title"]}
                old = local.get(variation_id)
                if old:
                    name, name_source = old["Name"], "pinned_clinvar_raw_Name"
                elif len(names) == 1:
                    name, name_source = next(iter(names)), "clingen_erepo_preferred_title"
                else:
                    titles = document.get("communityStandardTitle", [])
                    if len(titles) != 1:
                        raise ValueError("missing_or_ambiguous_name_annotation")
                    name, name_source = titles[0], "clingen_registry_community_title"
                allele_ids = {str(item["alleleId"]) for item in document.get("externalRecords", {}).get("ClinVarAlleles", [])}
                if len(allele_ids) != 1:
                    raise ValueError("ambiguous_clinvar_allele_id")
                row = {column: (old or {}).get(column, "-") for column in OUTPUT_SOURCE_COLUMNS}
                row.update(mapping)
                row.update({"VariationID": variation_id, "#AlleleID": next(iter(allele_ids)), "Name": name,
                    "GeneSymbol": ";".join(sorted({curation["gene"] for curation in curations})), "OriginSimple": "germline",
                    "target_label": label, "target": "1" if label == "Pathogenic" else "0",
                    "annotation_source": name_source + ";GRCh38_alleles_from_ClinGen_Registry",
                    "clingen_curation_uuids_json": json.dumps(sorted({curation["uuid"] for curation in curations})),
                    "clingen_conditions_json": json.dumps(sorted({curation["condition"] for curation in curations})),
                    "clingen_panels_json": json.dumps(sorted({curation["expert_panel"] for curation in curations})),
                    "clingen_doc_versions_json": json.dumps(sorted({curation["doc_version"] for curation in curations})),
                    "original_clinvar_labels_json": (old or {}).get("original_clinvar_labels_json", "[]")})
                resolved.append(row)
                status = "resolved_reference_checked"
                counts[f"resolved_{label.lower()}"] += 1
                counts[name_source] += 1
            except ValueError as error:
                status = str(error)
            counts[status] += 1
            audit.append({"VariationID": variation_id, "status": status})
        reference_metadata = reference.metadata
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(resolved)
    with (report_dir / "resolution_audit.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["VariationID", "status"], delimiter="\t")
        writer.writeheader()
        writer.writerows(audit)
    result = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "status": "reference_checked_pending_development_overlap",
        "reference": reference_metadata, "counts": dict(counts), "network_failures": failures,
        "inputs": {"screen_manifest": {"path": str(screen_path.resolve()), "sha256": sha256_file(screen_path)},
                   "local_mapping": {"path": str(local_mapped_path.resolve()), "sha256": sha256_file(local_mapped_path)}},
        "registry_cache": [{"caid": caid, "sha256": sha256_file(cache_dir / f"{caid}.json")} for caid in caids if caid not in failures],
        "output": {"path": str(output_path.resolve()), "sha256": sha256_file(output_path)},
        "provenance": {"target_source": "ClinGen_eRepo_expert_curation", "cohort_kind": "expert_curated_supplementary", "source_independence_verified": False, "model_loaded": False},
    }
    (report_dir / "manifest.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(counts), indent=2))
    return result


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description="ClinGen CAID → GRCh38 ana kromozom allel çözümlemesi.")
    parser.add_argument("--screen-manifest", type=Path, default=Path("reports/generated/clingen_erepo_screen/manifest.json"))
    parser.add_argument("--local-mapping", type=Path, default=Path("data/external/clingen_mapped_candidates.tsv.gz"))
    parser.add_argument("--reference-config", type=Path, default=Path("data/reference/ucsc_hg38_initial/grch38_with_aliases.json"))
    parser.add_argument("--cache-dir", type=Path, default=Path("data/external/clingen_registry_2026-10-26"))
    parser.add_argument("--output", type=Path, default=Path("data/external/clingen_resolved_candidates.tsv.gz"))
    parser.add_argument("--report-dir", type=Path, default=Path("reports/generated/clingen_allele_resolution"))
    args = parser.parse_args()
    resolve_candidates(args.screen_manifest, args.local_mapping, args.reference_config, args.cache_dir, args.output, args.report_dir)


if __name__ == "__main__":
    main()
