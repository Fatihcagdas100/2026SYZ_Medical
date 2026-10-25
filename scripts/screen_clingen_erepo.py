"""eRepo anlık görüntüsünde kesin ikili hedef ve eski ClinVar kimlik örtüşmesini tara.

Bu bir uygunluk taramasıdır; normalleştirilmiş koordinat/allel kabulü veya model
değerlendirmesi değildir. Eşleşmeyen kimlikler yalnızca üst sınır adayıdır.
"""

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
    from scripts.build_hgvs_features import derive_hgvs_features
except ModuleNotFoundError:
    from build_baseline_features import sha256_file
    from build_hgvs_features import derive_hgvs_features


DEFAULT_SNAPSHOT = Path("data/external/clingen_erepo_2026-10-26")
DEFAULT_DEVELOPMENT = Path("data/processed/clinvar_grch38_germline_exact_binary.tsv.gz")
DEFAULT_REPORT = Path("reports/generated/clingen_erepo_screen")
EXACT_LABELS = {"Pathogenic", "Benign"}


def load_snapshot(snapshot_dir: Path) -> tuple[dict, list[dict]]:
    manifest = json.loads((snapshot_dir / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("status") != "complete":
        raise ValueError("eRepo anlık görüntüsü tamamlanmış değil.")
    rows = []
    seen: set[str] = set()
    for item in [manifest["count_before"], manifest["count_after"], manifest["terms_of_use"], *manifest["pages"]]:
        path = (snapshot_dir / item["path"]).resolve()
        if not path.is_relative_to(snapshot_dir.resolve()) or sha256_file(path) != item["sha256"]:
            raise ValueError("Snapshot dosya hash'i/yolu uyuşmuyor.")
    for page in manifest["pages"]:
        with gzip.open(snapshot_dir / page["path"], "rt", encoding="utf-8") as handle:
            payload = json.load(handle)
        if payload.get("status", {}).get("code") != 200 or len(payload["data"]) != page["rows"]:
            raise ValueError("Snapshot sayfa satır sayısı/durumu uyuşmuyor.")
        for row in payload["data"]:
            if row["uuid"] in seen:
                raise ValueError("Snapshot UUID tekrarlandı.")
            seen.add(row["uuid"])
            rows.append(row)
    if len(rows) != manifest["downloaded_unique_classifications"] or len(rows) != manifest["expected_total"]:
        raise ValueError("Snapshot toplam satır sayısı uyuşmuyor.")
    return manifest, rows


def screen_snapshot(snapshot_dir: Path, development_path: Path, report_dir: Path) -> dict:
    snapshot, rows = load_snapshot(snapshot_dir)
    with gzip.open(development_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        previous_ids = set()
        previous_genes = set()
        for row in reader:
            previous_ids.add(row["VariationID"])
            previous_genes.update(gene.strip() for gene in row["GeneSymbol"].replace(",", ";").split(";") if gene.strip() not in {"", "-", "na"})
    report_dir.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    distribution: Counter[str] = Counter()
    by_cvid: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        counts["total_classification_records"] += 1
        label = row.get("classification", "")
        distribution[label] += 1
        if row.get("retracted") is not False:
            counts["retracted_or_unknown_status_records"] += 1
            continue
        counts["active_records"] += 1
        if label in EXACT_LABELS:
            counts["exact_binary_records"] += 1
        cvid = str(row.get("cvId") or "").strip()
        if not cvid.isdigit():
            counts["records_without_single_clinvar_id"] += 1
            if label in EXACT_LABELS:
                counts["exact_records_without_single_clinvar_id"] += 1
            continue
        by_cvid[cvid].append(row)
    candidate_records: list[dict] = []
    summary_rows: list[dict] = []
    for cvid, curations in sorted(by_cvid.items(), key=lambda item: int(item[0])):
        exact = [row for row in curations if row["classification"] in EXACT_LABELS]
        if not exact:
            continue
        counts["exact_binary_unique_clinvar_ids"] += 1
        labels = {row["classification"] for row in curations}
        overlaps = cvid in previous_ids
        if len(labels) > 1:
            eligibility = "excluded_discordant_classifications_across_curations"
            counts["discordant_unique_clinvar_ids"] += 1
        elif overlaps:
            eligibility = "excluded_previous_variation_id"
            counts["previous_id_overlap_unique_variants"] += 1
        else:
            eligibility = "id_unmatched_pending_allele_coordinate_and_provenance"
            counts["id_unmatched_candidate_unique_variants"] += 1
            label = next(iter(labels))
            counts[f"id_unmatched_{label.lower()}_variants"] += 1
            effects = {derive_hgvs_features(row.get("preferredVarTitle", ""))["protein_effect"] for row in exact}
            if effects == {"missense"}:
                counts["id_unmatched_unambiguous_missense_variants"] += 1
                counts[f"id_unmatched_missense_{label.lower()}_variants"] += 1
            genes = {row.get("gene", "") for row in exact} - {""}
            if genes and not genes.intersection(previous_genes):
                counts["id_unmatched_unseen_gene_variants"] += 1
            for row in exact:
                candidate_records.append({
                    "uuid": row["uuid"], "ClinVarVariationID": cvid,
                    "CAID": row.get("caId", ""), "classification": row["classification"],
                    "gene": row.get("gene", ""), "condition": row.get("condition", ""),
                    "expert_panel": row.get("ep", ""), "doc_version": row.get("docVersion", ""),
                    "approved_date": row.get("approvedDate", ""), "published_date": row.get("publishedDate", ""),
                    "preferred_variant_title": row.get("preferredVarTitle", ""),
                    "hgvs_json": json.dumps(row.get("hgvs", []), ensure_ascii=False),
                    "met_codes_json": json.dumps(row.get("metCodes", []), ensure_ascii=False),
                })
        summary_rows.append({
            "ClinVarVariationID": cvid, "curation_count": len(curations),
            "classifications": ";".join(sorted(labels)), "previous_variation_id_overlap": int(overlaps),
            "status": eligibility,
        })
    for key in ("id_unmatched_candidate_unique_variants", "id_unmatched_pathogenic_variants", "id_unmatched_benign_variants", "id_unmatched_unambiguous_missense_variants", "id_unmatched_missense_pathogenic_variants", "id_unmatched_missense_benign_variants", "previous_id_overlap_unique_variants", "discordant_unique_clinvar_ids"):
        counts.setdefault(key, 0)
    candidates_path = report_dir / "id_unmatched_curation_candidates.tsv"
    candidate_columns = ["uuid", "ClinVarVariationID", "CAID", "classification", "gene", "condition", "expert_panel", "doc_version", "approved_date", "published_date", "preferred_variant_title", "hgvs_json", "met_codes_json"]
    for path, records, fields in (
        (candidates_path, candidate_records, candidate_columns),
        (report_dir / "variant_screen.tsv", summary_rows, ["ClinVarVariationID", "curation_count", "classifications", "previous_variation_id_overlap", "status"]),
    ):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(records)
    result = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "decision": "screening_only_not_admitted_for_external_evaluation",
        "snapshot_manifest": {"path": str((snapshot_dir / "manifest.json").resolve()), "sha256": sha256_file(snapshot_dir / "manifest.json")},
        "snapshot_kind": snapshot["snapshot_kind"],
        "development_source": {"path": str(development_path.resolve()), "sha256": sha256_file(development_path), "unique_variation_ids": len(previous_ids)},
        "counts": dict(counts), "classification_record_distribution": dict(distribution),
        "policies": {
            "exact_labels": sorted(EXACT_LABELS),
            "multiple_curation_policy": "exclude_variant_if_any_active_classification_differs",
            "old_id_scope": "all_previously_used_variants_including_internal_test",
            "candidate_status": "upper_bound_only_pending_allele_coordinate_and_label_provenance",
            "coordinate_or_allele_overlap_checked": False,
            "source_independence_verified": False,
            "model_loaded": False, "model_or_threshold_tuned": False,
        },
        "outputs": {"curation_candidates": {"path": str(candidates_path.resolve()), "sha256": sha256_file(candidates_path)}},
    }
    with (report_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
    lines = ["# ClinGen eRepo uygunluk ve kimlik çakışması taraması", "",
        "Durum: yalnızca tarama; dış kohort kabulü ve model değerlendirmesi yapılmadı.", "",
        "| Ölçüm | Sayı |", "| --- | ---: |"]
    for name, count in counts.items():
        lines.append(f"| {name} | {count:,} |")
    lines.extend(["", "Kimliği eşleşmeyen aday sayısı üst sınırdır. Allel/koordinat eşdeğerliği, GRCh38 referansı ve etiket kökeni denetlenmeden bu kayıtlar bağımsız dış veri sayılmaz.", "",
        "İndirme canlı API'den yerel hash'li anlık görüntüdür; resmi sabit sürüm arşivi değildir. Her yorumun sürümü, hastalığı ve uzman paneli korunmuştur."])
    (report_dir / "clingen_erepo_screen.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(result["counts"], ensure_ascii=False, indent=2))
    print(f"Tarama raporu: {(report_dir / 'clingen_erepo_screen.md').resolve()}")
    return result


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description="eRepo snapshot'ını önceki varyant kimliklerine karşı tarar.")
    parser.add_argument("--snapshot-dir", type=Path, default=DEFAULT_SNAPSHOT)
    parser.add_argument("--development", type=Path, default=DEFAULT_DEVELOPMENT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    screen_snapshot(args.snapshot_dir, args.development, args.report_dir)


if __name__ == "__main__":
    main()
