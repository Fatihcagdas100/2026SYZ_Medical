"""Compare local GTF annotations with the saved REST pilot, never model scores."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sqlite3
import sys

import pandas as pd

try:
    from scripts import local_vep_pipeline as pipeline
    from scripts.pilot_context_annotations import RAW_COLUMNS, read_columns, select_transcript, extract_features
    from scripts.variant_normalization import ReferenceGenome
except ModuleNotFoundError:
    import local_vep_pipeline as pipeline
    from pilot_context_annotations import RAW_COLUMNS, read_columns, select_transcript, extract_features
    from variant_normalization import ReferenceGenome


def enrich_consequence(consequence, connection):
    result = dict(consequence)
    raw_identifier = result["transcript_id"]
    base = raw_identifier.split(".")[0]
    record = connection.execute("SELECT model FROM transcripts WHERE id = ?", (base,)).fetchone()
    if record is None:
        raise ValueError(f"Transcript not in locked GTF: {raw_identifier}")
    model = json.loads(record[0])
    if {"cds_start_NF", "cds_end_NF"}.intersection(model["tags"]):
        model["Translation"] = {}
    identifier = f"{base}.{model['version']}"
    if "." in raw_identifier and raw_identifier != identifier:
        raise ValueError("VEP/GTF transcript version mismatch")
    if int(result["strand"]) != model["strand"]:
        raise ValueError("VEP/GTF transcript strand mismatch")
    result["transcript_id"] = identifier
    result["gene_symbol"] = model["gene_symbol"]
    if result.get("biotype") and result["biotype"] != model["biotype"]:
        raise ValueError("VEP/GTF transcript biotype mismatch")
    result["biotype"] = model["biotype"]
    result["canonical"] = int("Ensembl_canonical" in model["tags"])
    result["mane_select"] = "gtf:MANE_Select" if "MANE_Select" in model["tags"] else ""
    return result, model


def context_features(record, gene_symbol, position, connection):
    enriched, models = [], {}
    alternate = record["allele_string"].split("/")[1]
    for consequence in record.get("transcript_consequences", []):
        if consequence.get("variant_allele") != alternate:
            continue
        result, model = enrich_consequence(consequence, connection)
        enriched.append(result)
        models[result["transcript_id"]] = model
    selected = select_transcript(enriched, gene_symbol)
    model = models[selected["transcript_id"]]
    if str(record["seq_region_name"]) != model["seq_region_name"]:
        raise ValueError("VEP/GTF transcript chromosome mismatch")
    features = extract_features(selected, model, position)
    # GTF has no domain annotation: missing is not a zero-domain measurement.
    features["pfam_annotation_count"] = None
    features["interpro_annotation_count"] = None
    return dict(selected_transcript=selected["transcript_id"], mane_select=selected["mane_select"], **features)


def equal_value(left, right):
    if pd.isna(left) or pd.isna(right) or left == "" or right == "":
        return (pd.isna(left) or left == "") and (pd.isna(right) or right == "")
    try:
        return math.isclose(float(left), float(right), rel_tol=1e-10, abs_tol=1e-10)
    except (ValueError, TypeError):
        return str(left) == str(right)


def compare(config, config_path):
    pipeline.protect(config)
    lock, lock_hash = pipeline.verify_resources(config, config_path)
    previous = pipeline.ROOT / "reports/generated/context_annotation_pilot"
    fetch_manifest = json.loads((previous / "fetch_manifest.json").read_text(encoding="utf-8"))
    audit_path = pipeline.checked(previous / "annotation_audit.tsv", fetch_manifest["outputs"]["annotation_audit.tsv"])
    old = pd.read_csv(audit_path, sep="\t", dtype=str, keep_default_na=False)
    pilot_config = json.loads((pipeline.ROOT / "config/context_annotation_pilot.json").read_text(encoding="utf-8"))
    raw = read_columns(pipeline.checked(pipeline.ROOT / config["inputs"]["variants"]["path"], config["inputs"]["variants"]["sha256"]), RAW_COLUMNS)
    item = pilot_config["inputs"]["diagnostics"]
    diagnostics = read_columns(pipeline.checked(pipeline.ROOT / item["path"], item["sha256"]), RAW_COLUMNS)
    tables = {"train_feasibility": raw.set_index("VariationID", drop=False),
              "seen_clingen_diagnostic": diagnostics.set_index("VariationID", drop=False)}
    directory = pipeline.ROOT / "reports/generated/gtf_context_pilot"
    inputs = directory / "input"
    inputs.mkdir(parents=True, exist_ok=True)
    lines = []
    with ReferenceGenome(pipeline.ROOT / config["inputs"]["reference_config"]["path"]) as reference:
        for row in old.to_dict("records"):
            source = tables[row["role"]].loc[row["VariationID"]].to_dict()
            chrom, pos, ref, alt = reference.normalize(source)
            lines.append(f"{chrom}\t{pos}\t{row['VariationID']}\t{ref}\t{alt}\t.\t.\t.\n")
    lines.sort(key=lambda line: (int(line.split("\t")[0]) if line.split("\t")[0].isdigit()
                                else {"X": 23, "Y": 24}[line.split("\t")[0]], int(line.split("\t")[1])))
    vcf = inputs / "pilot.vcf"
    pipeline.write_once(vcf, (pipeline.VCF_HEADER + "".join(lines)).encode("utf-8"))
    synonyms = pipeline.ROOT / config["prepared_dir"] / "chromosome_synonyms.tsv"
    pipeline.write_once(inputs / synonyms.name, synonyms.read_bytes())
    output = directory / "annotations"
    output.mkdir(exist_ok=True)
    raw_output = output / "annotations.jsonl"
    receipt_path = directory / "annotation.receipt.json"
    command = pipeline.build_command(config, inputs, {"path": vcf.name}, output, lock)
    if receipt_path.exists():
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        if receipt["resource_lock_sha256"] != lock_hash or receipt["input_sha256"] != pipeline.sha256_file(vcf):
            raise ValueError("Pilot input/resource lock changed")
        if json.loads((directory / "command.json").read_text(encoding="utf-8")) != ["docker"] + command:
            raise ValueError("Pilot annotation command changed")
        pipeline.checked(raw_output, receipt["output_sha256"])
    else:
        pipeline.write_json(directory / "command.json", ["docker"] + command)
        result = pipeline.docker(command, timeout=600)
        (directory / "stdout.txt").write_text(result.stdout, encoding="utf-8")
        (directory / "stderr.txt").write_text(result.stderr, encoding="utf-8")
        if result.returncode:
            raise RuntimeError("Local GTF pilot failed: " + result.stderr)
        pipeline.validate_output(vcf, raw_output)
        pipeline.write_json(receipt_path, dict(resource_lock_sha256=lock_hash, input_sha256=pipeline.sha256_file(vcf),
                                               output_sha256=pipeline.sha256_file(raw_output)))
    records = {str(r["id"]): r for r in (json.loads(line) for line in raw_output.read_text(encoding="utf-8").splitlines())}
    rows = []
    connection = sqlite3.connect(Path(lock["gtf_dir"]) / "transcripts.sqlite")
    try:
        for row in old.to_dict("records"):
            record = records[row["VariationID"]]
            features = context_features(record, row["GeneSymbol"], int(record["start"]), connection)
            rows.append({k: row[k] for k in ["VariationID", "GeneSymbol", "role"]} | features)
    finally:
        connection.close()
    audit = pd.DataFrame(rows)
    audit.to_csv(directory / "local_annotation_audit.tsv", sep="\t", index=False)
    feature_columns = [c for c in rows[0] if c not in {"GeneSymbol", "role", "selected_transcript", "mane_select",
                                                     "pfam_annotation_count", "interpro_annotation_count"}]
    for role, name in [("train_feasibility", "train_context_features.tsv"), ("seen_clingen_diagnostic", "seen_diagnostic_context_features.tsv")]:
        audit.loc[audit.role.eq(role), feature_columns].to_csv(directory / name, sep="\t", index=False)
    comparisons = []
    for local, rest in zip(rows, old.to_dict("records")):
        for column in ["selected_transcript"] + feature_columns[1:]:
            comparisons.append(dict(VariationID=local["VariationID"], column=column, rest_value=rest[column],
                                    local_value=local[column], equal=equal_value(rest[column], local[column])))
    frame = pd.DataFrame(comparisons)
    frame.to_csv(directory / "comparison.tsv", sep="\t", index=False)
    manifest = dict(created_at_utc=pipeline.now(), status="gtf_feasibility_not_model_validation", rows=len(rows),
                    script_sha256=pipeline.sha256_file(Path(__file__)), config_sha256=pipeline.sha256_file(config_path),
                    resource_lock_sha256=lock_hash, saved_rest_audit_sha256=pipeline.sha256_file(audit_path),
                    command_sha256=pipeline.sha256_file(directory / "command.json"),
                    agreement_by_column=frame.groupby("column").equal.agg(["sum", "count"]).to_dict("index"),
                    differences=frame.loc[~frame.equal].to_dict("records"), excluded_features=["pfam_annotation_count", "interpro_annotation_count"],
                    policy=dict(model_loaded=False, labels_loaded=False, predictions_recomputed=False, previously_seen_cases=True),
                    outputs={p.name: pipeline.sha256_file(p) for p in directory.glob("*.tsv")})
    pipeline.protect(config)
    pipeline.write_json(directory / "manifest.json", manifest)
    print(json.dumps({k: manifest[k] for k in ["rows", "agreement_by_column", "differences"]}, indent=2))
    return manifest


def export_completed(config, config_path):
    pipeline.protect(config)
    prepared, bundle = pipeline.verify_bundle(config, config_path)
    lock, lock_hash = pipeline.verify_resources(config, config_path)
    mapping = pd.read_csv(prepared / "variant_map.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("VariationID")
    rows, receipts = [], []
    directory = pipeline.ROOT / config["output_dir"]
    connection = sqlite3.connect(Path(lock["gtf_dir"]) / "transcripts.sqlite")
    try:
        for chunk in bundle["chunks"]:
            chunk_dir = directory / Path(chunk["path"]).stem
            receipt_path = chunk_dir / "complete.json"
            if not receipt_path.exists():
                continue
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            if receipt["input_sha256"] != chunk["sha256"] or receipt["resource_sha256"] != lock_hash:
                raise ValueError("Completed chunk resource/input identity differs")
            output = pipeline.checked(chunk_dir / receipt["output"], receipt["output_sha256"])
            pipeline.validate_output(prepared / chunk["path"], output)
            receipts.append(dict(path=str(receipt_path), sha256=pipeline.sha256_file(receipt_path)))
            with output.open(encoding="utf-8") as handle:
                for line in handle:
                    record = json.loads(line)
                    identifier = str(record["id"])
                    if identifier not in mapping.index:
                        raise ValueError("Annotated ID is outside prepared train mapping")
                    gene = mapping.loc[identifier, "GeneSymbol"]
                    if record.get("transcript_consequences"):
                        features = context_features(record, gene, int(record["start"]), connection)
                    else:
                        features = dict(selected_transcript="", mane_select="", biotype="unknown",
                                        consequence_terms=record["most_severe_consequence"], exon_rank=None, exon_total=None,
                                        location="no_transcript", donor_distance_bp=None, acceptor_distance_bp=None,
                                        nearest_splice_boundary_bp=None, protein_relative_position=None,
                                        pfam_annotation_count=None, interpro_annotation_count=None)
                    rows.append(dict(VariationID=identifier, GeneSymbol=gene, **features))
    finally:
        connection.close()
    if not rows:
        raise ValueError("No completed GTF chunks to export")
    audit = pd.DataFrame(rows)
    if audit.VariationID.duplicated().any():
        raise ValueError("Duplicate annotated train IDs")
    out = pipeline.ROOT / "reports/generated/local_vep_gtf_features"
    out.mkdir(parents=True, exist_ok=True)
    audit.to_csv(out / "train_annotation_audit.tsv", sep="\t", index=False)
    exclude = {"GeneSymbol", "selected_transcript", "mane_select", "pfam_annotation_count", "interpro_annotation_count"}
    audit[[c for c in audit if c not in exclude]].to_csv(out / "train_context_features.tsv", sep="\t", index=False)
    manifest = dict(created_at_utc=pipeline.now(), rows=len(rows), chunks=len(receipts), expected_rows=bundle["counts"]["accepted"],
                    script_sha256=pipeline.sha256_file(Path(__file__)), config_sha256=pipeline.sha256_file(config_path),
                    complete=len(rows) == bundle["counts"]["accepted"], resource_lock_sha256=lock_hash, source_receipts=receipts,
                    available_geometry=int(audit.nearest_splice_boundary_bp.notna().sum()),
                    available_protein_position=int(audit.protein_relative_position.notna().sum()),
                    policy=dict(labels_loaded=False, model_loaded=False, partitions=["train"], domains_unavailable=True),
                    outputs={p.name: pipeline.sha256_file(p) for p in out.glob("*.tsv")})
    pipeline.protect(config)
    pipeline.write_json(out / "manifest.json", manifest)
    print(json.dumps({k: manifest[k] for k in ["rows", "chunks", "complete", "available_geometry", "available_protein_position"]}, indent=2))
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=pipeline.ROOT / "config/local_vep_gtf.json")
    parser.add_argument("--completed-batches", action="store_true", help="Export train features from verified completed GTF chunks")
    args = parser.parse_args()
    action = export_completed if args.completed_batches else compare
    action(pipeline.load_config(args.config), args.config)


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
