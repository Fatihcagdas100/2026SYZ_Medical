"""Small, label-free GRCh38 SNV annotation pilot; never loads or scores a model.

Only --fetch permits network calls. Cached replies have URL, UTC and SHA256
receipts. A cached replay does not claim that the live service is unchanged.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import time
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

try:
    from scripts.build_baseline_features import sha256_file
    from scripts.variant_normalization import ReferenceGenome
except ModuleNotFoundError:
    from build_baseline_features import sha256_file
    from variant_normalization import ReferenceGenome

ROOT = Path(__file__).resolve().parents[1]
BASE = "https://rest.ensembl.org"
RAW_COLUMNS = ["VariationID", "Name", "GeneSymbol", "Assembly", "ChromosomeAccession",
               "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF"]
FEATURE_COLUMNS = ["biotype", "consequence_terms", "exon_rank", "exon_total",
                   "location", "donor_distance_bp", "acceptor_distance_bp",
                   "nearest_splice_boundary_bp", "protein_relative_position",
                   "pfam_annotation_count", "interpro_annotation_count"]
PARAMETERS = dict(domains=1, canonical=1, mane=1, numbers=1, protein=1, transcript_version=1)


class AnnotationRequestFailure(RuntimeError):
    """Public service failed after retries; not an annotated variant."""


def read_columns(path, columns):
    # Do not load target/clinical significance columns into the pilot dataframe.
    return pd.read_csv(path, sep="\t", usecols=columns, dtype=str, keep_default_na=False)


def verify(item):
    path = ROOT / item["path"]
    if sha256_file(path) != item["sha256"]:
        raise ValueError(f"Input/protected file hash mismatch: {path}")
    return path


def select_training_rows(raw, features, split, strata, per_stratum, seed):
    if split.VariationID.duplicated().any():
        raise ValueError("Duplicate split ID")
    train = split.loc[split.split.eq("train"), ["VariationID"]]
    joined = raw.merge(train, on="VariationID", validate="one_to_one").merge(
        features, on="VariationID", validate="one_to_one")
    # First pilot is deliberately restricted to primary-chromosome biallelic SNVs.
    mask = (joined.ReferenceAlleleVCF.str.fullmatch("[ACGT]") &
            joined.AlternateAlleleVCF.str.fullmatch("[ACGT]") &
            joined.ReferenceAlleleVCF.ne(joined.AlternateAlleleVCF) &
            joined.Chromosome.str.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y)") &
            joined.ChromosomeAccession.str.startswith("NC_"))
    joined = joined.loc[mask].copy()
    joined["stratum"] = joined.splice_offset_category.where(
        joined.splice_offset_category.isin(strata), joined.protein_effect)
    joined["rank"] = joined.VariationID.map(
        lambda v: hashlib.sha256(f"{seed}:{v}".encode()).hexdigest())
    chosen = []
    for stratum in strata:
        rows = joined.loc[joined.stratum.eq(stratum)].sort_values(["rank", "VariationID"]).head(per_stratum)
        if len(rows) != per_stratum:
            raise ValueError(f"Not enough training SNVs for {stratum}")
        chosen.extend(rows.to_dict("records"))
    for row in chosen:
        row["role"] = "train_feasibility"
        row.pop("rank")
    return chosen


class CachedRest:
    def __init__(self, directory, fetch=False):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)
        self.fetch = fetch
        self.used = []

    def get(self, endpoint, params=None, tag=""):
        params = {"content-type": "application/json", **(params or {})}
        url = BASE + endpoint + "?" + urlencode(sorted(params.items()))
        key = hashlib.sha256((url + "#" + tag).encode()).hexdigest()
        body = self.directory / f"{key}.json"
        receipt = self.directory / f"{key}.receipt.json"
        failure = self.directory / f"{key}.failure.json"
        if body.exists() != receipt.exists():
            raise ValueError(f"Incomplete cache pair: {key}")
        if not body.exists():
            if not self.fetch:
                if failure.exists():
                    metadata = json.loads(failure.read_text(encoding="utf-8"))
                    if metadata["url"] != url:
                        raise ValueError("Failure receipt URL mismatch")
                    self.used.append({"failure_receipt": failure.name, "sha256": sha256_file(failure), **metadata})
                    raise AnnotationRequestFailure(metadata["attempts"][-1]["error"])
                raise FileNotFoundError(f"No verified cache; use --fetch: {url}")
            for attempt in range(3):
                try:
                    with urlopen(Request(url, headers={"Accept": "application/json",
                                                         "User-Agent": "2026SYZ-context-pilot/1"}), timeout=25) as response:
                        payload = response.read()
                    json.loads(payload)  # Reject non-JSON errors before storing.
                    break
                except (OSError, ValueError) as error:
                    if attempt == 2:
                        metadata = json.loads(failure.read_text(encoding="utf-8")) if failure.exists() else {"url": url, "attempts": []}
                        metadata["attempts"].append({"retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                                                    "error_type": type(error).__name__, "error": str(error),
                                                    "http_status": getattr(error, "code", None), "retries": 3})
                        failure.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
                        self.used.append({"failure_receipt": failure.name, "sha256": sha256_file(failure), **metadata})
                        raise AnnotationRequestFailure(str(error)) from error
                    time.sleep(1 + attempt)
            body.write_bytes(payload)
            receipt.write_text(json.dumps({"url": url, "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
                                           "sha256": sha256_file(body)}, indent=2), encoding="utf-8")
            time.sleep(0.2)
        metadata = json.loads(receipt.read_text(encoding="utf-8"))
        if metadata["url"] != url or metadata["sha256"] != sha256_file(body):
            raise ValueError(f"Corrupt cache: {key}")
        self.used.append({"body": body.name, "receipt": receipt.name,
                          "receipt_sha256": sha256_file(receipt), **metadata})
        return json.loads(body.read_bytes())


def select_transcript(consequences, gene_symbol):
    symbols = set(gene_symbol.replace(";", ",").split(","))
    valid = [t for t in consequences if t.get("transcript_id")]
    if not valid:
        raise ValueError("No transcript consequence")
    # Target-free, fixed order. Keep every other transcript in the raw cache.
    return min(valid, key=lambda t: (t.get("gene_symbol") not in symbols,
                                    not bool(t.get("mane_select")),
                                    not bool(t.get("canonical")), t["transcript_id"]))


def splice_geometry(position, transcript):
    exons = sorted(transcript.get("Exon", []), key=lambda e: e["start"])
    strand = int(transcript["strand"])
    if strand not in {-1, 1} or not exons:
        raise ValueError("Invalid transcript exon geometry")
    location, donor, acceptor = "outside_transcript", None, None
    for i, exon in enumerate(exons):
        start, end = int(exon["start"]), int(exon["end"])
        if start <= position <= end:
            location = "exonic"
            # Transcript ends are NOT splice junctions; missing is not zero.
            donor = end - position if i < len(exons) - 1 else None
            acceptor = position - start if i > 0 else None
            if strand == -1:
                donor, acceptor = acceptor, donor
            break
        if i < len(exons) - 1 and end < position < int(exons[i + 1]["start"]):
            location = "intronic"
            donor, acceptor = position - end, int(exons[i + 1]["start"]) - position
            if strand == -1:
                donor, acceptor = acceptor, donor
            break
    distances = [x for x in (donor, acceptor) if x is not None]
    return dict(location=location, donor_distance_bp=donor, acceptor_distance_bp=acceptor,
                nearest_splice_boundary_bp=min(distances) if distances else None)


def extract_features(consequence, transcript, position):
    geometry = splice_geometry(position, transcript)
    exon = str(consequence.get("exon", "")).split("/")
    protein_start = consequence.get("protein_start")
    length = transcript.get("Translation", {}).get("length")
    ratio = float(protein_start) / float(length) if protein_start and length else None
    domains = consequence.get("domains", [])
    result = dict(biotype=consequence.get("biotype", "unknown"),
                  consequence_terms=";".join(sorted(consequence.get("consequence_terms", []))),
                  exon_rank=int(exon[0]) if len(exon) == 2 else None,
                  exon_total=int(exon[1]) if len(exon) == 2 else None,
                  protein_relative_position=ratio,
                  pfam_annotation_count=len({d["name"] for d in domains if d.get("db") == "Pfam"}),
                  interpro_annotation_count=len({d["name"] for d in domains if d.get("db") == "InterPro"}),
                  **geometry)
    assert set(result) == set(FEATURE_COLUMNS)
    return result


def validate_vep(result, key):
    chromosome, position, reference, alternate = key
    if (result.get("assembly_name") != "GRCh38" or str(result.get("seq_region_name")) != chromosome
            or int(result.get("start", -1)) != int(position) or int(result.get("end", -1)) != int(position)
            or result.get("allele_string") != reference + "/" + alternate):
        raise ValueError("VEP coordinate/reference/alternate mismatch")


def annotate(row, reference, client):
    key = reference.normalize(row)
    chrom, pos, ref, alt = key
    if len(ref) != 1 or len(alt) != 1:
        raise ValueError("Pilot only supports SNVs")
    replies = client.get(f"/vep/human/region/{chrom}:{pos}-{pos}:1/{alt}", PARAMETERS)
    if not isinstance(replies, list) or len(replies) != 1:
        raise ValueError("Unexpected VEP response")
    result = replies[0]
    validate_vep(result, key)
    consequences = [t for t in result.get("transcript_consequences", []) if t.get("variant_allele") == alt]
    selected = select_transcript(consequences, row["GeneSymbol"])
    transcript_id = selected["transcript_id"]
    transcript = client.get("/lookup/id/" + transcript_id.split(".")[0], {"expand": 1})
    actual_id = f"{transcript['id']}.{transcript['version']}"
    if (actual_id != transcript_id or transcript.get("assembly_name") != "GRCh38"
            or str(transcript.get("seq_region_name")) != chrom
            or int(transcript["strand"]) != int(selected["strand"])):
        raise ValueError("Transcript version/assembly/strand mismatch")
    features = extract_features(selected, transcript, int(pos))
    audit = {k: row[k] for k in ("VariationID", "Name", "GeneSymbol", "role", "stratum")}
    audit.update(selected_transcript=transcript_id, mane_select=selected.get("mane_select", ""),
                 selected_gene=selected.get("gene_symbol", ""),
                 gene_matches=selected.get("gene_symbol") in row["GeneSymbol"].replace(";", ",").split(","),
                 transcript_count=len(consequences),
                 consequence_sets="|".join(sorted({";".join(sorted(t.get("consequence_terms", []))) for t in consequences})),
                 **features)
    return audit


def run(config_path, fetch=False):
    config = json.loads(config_path.read_text(encoding="utf-8"))
    paths = {name: verify(item) for name, item in config["inputs"].items()}
    for item in config["protected"]:
        verify(item)
    raw = read_columns(paths["variants"], RAW_COLUMNS)
    features = read_columns(paths["features"], ["VariationID", "protein_effect", "splice_offset_category"])
    split = read_columns(paths["split"], ["VariationID", "split"])
    rows = select_training_rows(raw, features, split, config["strata"], config["per_stratum"], config["seed"])
    diagnostics = read_columns(paths["diagnostics"], RAW_COLUMNS).set_index("VariationID", drop=False)
    for variant_id in config["diagnostic_ids"]:
        row = diagnostics.loc[variant_id].to_dict()
        row.update(role="seen_clingen_diagnostic", stratum="posthoc_diagnostic")
        rows.append(row)
    client = CachedRest(ROOT / config["cache"], fetch)
    # A resumed --fetch may add replies to an older cache. Check the live release
    # as well as the original receipt, so a service upgrade cannot silently mix
    # annotation releases. Preserve every earlier response and receipt.
    live_tag = datetime.now(timezone.utc).isoformat()
    live_before = client.get("/info/data", tag="live_before_" + live_tag) if fetch else None
    before = client.get("/info/data", tag="release_before")
    if before.get("releases") != [config["ensembl_release"]]:
        raise ValueError("Unexpected Ensembl release")
    if fetch and live_before != before:
        raise ValueError("Live release differs from the original cache release")
    rest_info = client.get("/info/rest")
    software_info = client.get("/info/software")
    audit_rows, failures = [], []
    with ReferenceGenome(ROOT / config["reference_config"]) as reference:
        reference_metadata = reference.metadata
        for row in rows:
            print(f"Annotating {row['VariationID']} ({row['role']}, {row['stratum']})", flush=True)
            try:
                audit_rows.append(annotate(row, reference, client))
            except AnnotationRequestFailure as error:
                failures.append({k: row[k] for k in ("VariationID", "role", "stratum")})
                failures[-1]["error"] = str(error)
                print(f"Annotation unavailable for {row['VariationID']}: {error}", flush=True)
    after = client.get("/info/data", tag="release_after")
    if before != after:
        raise ValueError("Ensembl release changed during pilot")
    if fetch and client.get("/info/data", tag="live_after_" + live_tag) != before:
        raise ValueError("Live Ensembl release changed during pilot")
    for item in config["protected"]:
        verify(item)
    output = ROOT / config["output"]
    output.mkdir(parents=True, exist_ok=True)
    audit_columns = ["VariationID", "Name", "GeneSymbol", "role", "stratum", "selected_transcript",
                     "mane_select", "selected_gene", "gene_matches", "transcript_count", "consequence_sets"] + FEATURE_COLUMNS
    audit = pd.DataFrame(audit_rows, columns=audit_columns)
    audit.to_csv(output / "annotation_audit.tsv", sep="\t", index=False)
    pd.DataFrame(failures, columns=["VariationID", "role", "stratum", "error"]).to_csv(
        output / "annotation_failures.tsv", sep="\t", index=False)
    for role, name in (("train_feasibility", "train_candidate_features.tsv"),
                       ("seen_clingen_diagnostic", "seen_diagnostic_features.tsv")):
        audit.loc[audit.role.eq(role), ["VariationID"] + FEATURE_COLUMNS].to_csv(output / name, sep="\t", index=False)
    coverage = {}
    for role, group in audit.groupby("role"):
        coverage[role] = dict(rows=len(group), mane_select=int(group.mane_select.ne("").sum()),
                             gene_matches=int(group.gene_matches.sum()),
                             splice_geometry=int(group.nearest_splice_boundary_bp.notna().sum()),
                             protein_position=int(group.protein_relative_position.notna().sum()),
                             pfam_annotated=int(group.pfam_annotation_count.gt(0).sum()),
                             interpro_annotated=int(group.interpro_annotation_count.gt(0).sum()),
                             multiple_consequence_sets=int(group.consequence_sets.str.contains(r"\|").sum()))
    manifest = dict(created_at_utc=datetime.now(timezone.utc).isoformat(), status="annotation_feasibility_only",
                    config_sha256=sha256_file(config_path), ensembl_release=before, rest_info=rest_info,
                    script_sha256=sha256_file(Path(__file__)), live_release_checked=fetch,
                    software_info=software_info, reference=reference_metadata, inputs=config["inputs"],
                    protected=config["protected"], protected_hashes_unchanged=True,
                    policy=dict(labels_loaded=False, model_loaded=False, model_refit=False, predictions_recomputed=False,
                                diagnostic_cases_previously_seen=True, gene_and_transcript_ids_are_audit_only=True),
                    transcript_selection="gene symbol match, MANE Select, canonical, lexicographic transcript ID",
                    attempted_counts=dict(pd.Series([r["role"] for r in rows]).value_counts().items()),
                    failures=failures, coverage=coverage, cache=client.used,
                    outputs={p.name: sha256_file(p) for p in sorted(output.glob("*.tsv"))})
    manifest_text = json.dumps(manifest, indent=2, ensure_ascii=False) + "\n"
    (output / "manifest.json").write_text(manifest_text, encoding="utf-8")
    # Keep the live acquisition record when an offline replay is performed.
    (output / ("fetch_manifest.json" if fetch else "replay_manifest.json")).write_text(manifest_text, encoding="utf-8")
    print(json.dumps(coverage, indent=2), flush=True)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/context_annotation_pilot.json")
    parser.add_argument("--fetch", action="store_true", help="Allow retrieval of missing public annotations")
    args = parser.parse_args()
    run(args.config, args.fetch)


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
