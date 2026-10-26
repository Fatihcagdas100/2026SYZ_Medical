"""Prepare label-free train VCF chunks and run digest-pinned, offline VEP.

Commands: prepare, check, lock-resources, run. No model or evaluation is loaded.
Resource installation is separate; this script never silently downloads a cache.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.request import Request, urlopen
from uuid import uuid4

import pandas as pd

try:
    from scripts.pilot_context_annotations import RAW_COLUMNS, read_columns, validate_vep
    from scripts.variant_normalization import InvalidVariant, ReferenceGenome
    from scripts.build_baseline_features import sha256_file
except ModuleNotFoundError:
    from pilot_context_annotations import RAW_COLUMNS, read_columns, validate_vep
    from variant_normalization import InvalidVariant, ReferenceGenome
    from build_baseline_features import sha256_file

ROOT = Path(__file__).resolve().parents[1]
VCF_HEADER = "##fileformat=VCFv4.2\n##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
VEP_FLAGS = ["--offline", "--cache", "--format", "vcf", "--json", "--no_stats", "--no_progress",
             "--check_ref", "--domains", "--canonical", "--mane", "--numbers",
             "--protein", "--transcript_version", "--total_length", "--symbol"]


def now():
    return datetime.now(timezone.utc).isoformat()


def canonical(value):
    return (json.dumps(value, sort_keys=True, ensure_ascii=False, indent=2) + "\n").encode("utf-8")


def checked(path, expected):
    path = Path(path)
    if not path.is_file() or sha256_file(path) != expected:
        raise ValueError(f"Missing file or SHA256 mismatch: {path}")
    return path


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(canonical(value))
    temporary.replace(path)


def write_once(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != payload:
            raise ValueError(f"Existing prepared file differs; use another prepared_dir: {path}")
    else:
        temporary = path.with_name(path.name + ".tmp")
        temporary.write_bytes(payload)
        temporary.replace(path)
    return {"path": path.name, "sha256": sha256_file(path)}


def load_config(path):
    config = json.loads(Path(path).read_text(encoding="utf-8"))
    if (config["scope"] != "train_primary_biallelic_snv_only" or config["release"] != 116
            or not re.fullmatch(r"ensemblorg/ensembl-vep@sha256:[a-f0-9]{64}", config["image"])
            or config["platform"] != "linux/amd64" or config["chunk_size"] < 1
            or config.get("annotation_mode", "cache") not in {"cache", "gtf"}):
        raise ValueError("Unsupported scope/release/platform, unpinned image or invalid chunk size")
    return config


def protect(config):
    for item in config["protected"]:
        checked(ROOT / item["path"], item["sha256"])


def preparation_hash(config):
    # Cache/output locations can change without rebuilding identical train VCFs.
    recipe = {key: config[key] for key in ("schema_version", "scope", "chunk_size", "inputs")}
    return hashlib.sha256(canonical(recipe)).hexdigest()


def verify_bundle(config, config_path):
    directory = ROOT / config["prepared_dir"]
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    if "preparation_sha256" in manifest:
        if manifest["preparation_sha256"] != preparation_hash(config):
            raise ValueError("Prepared input recipe differs")
    elif manifest["config_sha256"] != sha256_file(config_path):
        raise ValueError("Legacy prepared input config differs; verify with original config first")
    for item in manifest["chunks"] + manifest["tables"]:
        checked(directory / item["path"], item["sha256"])
    return directory, manifest


def prepare(config, config_path):
    protect(config)
    paths = {k: checked(ROOT / v["path"], v["sha256"]) for k, v in config["inputs"].items()}
    directory = ROOT / config["prepared_dir"]
    if (directory / "manifest.json").exists():
        _, manifest = verify_bundle(config, config_path)
        if "preparation_sha256" not in manifest:
            manifest["preparation_sha256"] = preparation_hash(config)
            write_json(directory / "manifest.json", manifest)
        print("Verified existing prepared bundle:", manifest["counts"])
        return manifest
    split = read_columns(paths["split"], ["VariationID", "split"])
    if split.VariationID.duplicated().any() or not set(split.split) <= {"train", "validation", "test"}:
        raise ValueError("Invalid split IDs or partitions")
    raw = read_columns(paths["variants"], RAW_COLUMNS)
    if raw.VariationID.duplicated().any() or set(raw.VariationID) != set(split.VariationID):
        raise ValueError("Variant/split IDs must match one-to-one")
    train = raw.merge(split.loc[split.split.eq("train"), ["VariationID"]], on="VariationID", validate="one_to_one")
    accepted, excluded = [], []
    with ReferenceGenome(paths["reference_config"]) as reference:
        metadata = reference.metadata
        for row in train.to_dict("records"):
            ref, alt = row["ReferenceAlleleVCF"], row["AlternateAlleleVCF"]
            reason = None
            if not re.fullmatch("[ACGT]", ref) or not re.fullmatch("[ACGT]", alt) or ref == alt:
                reason = "outside_biallelic_snv_scope"
            elif not re.fullmatch(r"(?:[1-9]|1[0-9]|2[0-2]|X|Y)", row["Chromosome"]):
                reason = "outside_primary_chromosome_scope"
            else:
                try:
                    chrom, pos, normalized_ref, normalized_alt = reference.normalize(row)
                except InvalidVariant as error:
                    if str(error) != "nonprimary_or_unmapped_reference_accession":
                        raise  # REF/assembly/accession mismatch must not be silently excluded.
                    reason = str(error)
            if reason:
                excluded.append(dict(VariationID=row["VariationID"], reason=reason))
            else:
                accepted.append(dict(VariationID=row["VariationID"], GeneSymbol=row["GeneSymbol"],
                                     chromosome=chrom, position=int(pos), reference=normalized_ref, alternate=normalized_alt))
    accepted.sort(key=lambda r: (int(r["chromosome"]) if r["chromosome"].isdigit() else {"X": 23, "Y": 24}[r["chromosome"]],
                                 r["position"], int(r["VariationID"])))
    chunks = []
    for start in range(0, len(accepted), config["chunk_size"]):
        rows = accepted[start:start + config["chunk_size"]]
        name = f"chunk_{len(chunks) + 1:06d}.vcf"
        payload = VCF_HEADER + "".join(f"{r['chromosome']}\t{r['position']}\t{r['VariationID']}\t{r['reference']}\t{r['alternate']}\t.\t.\t.\n" for r in rows)
        item = write_once(directory / name, payload.encode("utf-8"))
        item["rows"] = len(rows)
        chunks.append(item)
    tables = []
    for name, rows, columns in (("variant_map.tsv", accepted, ["VariationID", "GeneSymbol", "chromosome", "position", "reference", "alternate"]),
                                ("excluded_train_rows.tsv", excluded, ["VariationID", "reason"])):
        tables.append(write_once(directory / name, pd.DataFrame(rows, columns=columns).to_csv(sep="\t", index=False, lineterminator="\n").encode("utf-8")))
    synonyms = "".join(f"{c}\tchr{c}\n" for c in [str(i) for i in range(1, 23)] + ["X", "Y"])
    tables.append(write_once(directory / "chromosome_synonyms.tsv", synonyms.encode("ascii")))
    manifest = dict(created_at_utc=now(), config_sha256=sha256_file(config_path), preparation_sha256=preparation_hash(config),
                    script_sha256=sha256_file(Path(__file__)),
                    scope=config["scope"], inputs=config["inputs"], reference=metadata, chunks=chunks, tables=tables,
                    counts=dict(train_rows=len(train), accepted=len(accepted), excluded=len(excluded), chunks=len(chunks)),
                    exclusion_reasons=dict(Counter(r["reason"] for r in excluded)),
                    policy=dict(labels_loaded=False, partitions=["train"], model_loaded=False, diagnostics_included=False))
    assert len(accepted) + len(excluded) == len(train)
    protect(config)
    write_json(directory / "manifest.json", manifest)
    print(json.dumps(manifest["counts"], indent=2))
    return manifest


def docker(arguments, timeout=60):
    return subprocess.run(["docker"] + arguments, capture_output=True, text=True, encoding="utf-8", errors="replace", timeout=timeout)


def existing_parent(path):
    path = Path(path).resolve()
    while not path.exists():
        path = path.parent
    return path


def preflight(config, inspect_remote=False):
    gtf_mode = config.get("annotation_mode") == "gtf"
    archive = dict(config["gtf_source"] if gtf_mode else config["archive"])
    if inspect_remote:
        with urlopen(Request(archive["url"], method="HEAD"), timeout=25) as response:
            archive["live_headers"] = dict(response.headers)
            archive["checked_at_utc"] = now()
            if int(response.headers["Content-Length"]) != archive.get("content_length_bytes", archive.get("observed_content_length_bytes")):
                raise ValueError("Provider archive size changed; review source lock")
    cache = ROOT / config["gtf_dir" if gtf_mode else "cache_dir"]
    cache_info = cache / "manifest.json" if gtf_mode else cache / "homo_sapiens" / f"{config['release']}_GRCh38" / "info.txt"
    free = shutil.disk_usage(existing_parent(cache)).free
    required = (config["working_reserve_bytes"] if cache_info.exists() else config["setup_planning_bytes"]) if gtf_mode else (
        config["working_reserve_bytes"] + (0 if cache_info.exists() else 2 * archive["observed_content_length_bytes"]))
    engine, image, version = False, False, ""
    error = ""
    try:
        result = docker(["version", "--format", "{{json .Server.Version}}"])
        engine = result.returncode == 0
        if engine:
            result = docker(["image", "inspect", config["image"], "--format", "{{.Id}}"])
            image = result.returncode == 0
            if image:
                result = docker(["run", "--rm", "--pull=never", "--network=none", "--platform", config["platform"], config["image"], "vep", "--help"])
                version = (result.stdout + result.stderr).strip()
                image = result.returncode == 0 and bool(re.search(r"\b116(?:\.\d+)?\b", version))
        error = result.stderr.strip() if result.returncode else ""
    except (OSError, subprocess.TimeoutExpired) as exception:
        error = str(exception)
    blockers = []
    for condition, reason in ((engine, "docker_engine_unavailable"), (image, "pinned_vep_image_unavailable_or_wrong_version"),
                              (cache_info.exists(), "gtf_resources_missing" if gtf_mode else "release_116_grch38_cache_missing"), (free >= required, "disk_below_planning_budget"),
                              ((ROOT / config["resource_lock"]).exists(), "resources_not_locked")):
        if not condition:
            blockers.append(reason)
    report = dict(created_at_utc=now(), ready=not blockers, blockers=blockers, docker_engine=engine, pinned_image=image,
                  annotation_mode="gtf" if gtf_mode else "cache", image=config["image"], vep_help=version, docker_error=error, cache_dir=str(cache.resolve()),
                  cache_info_present=cache_info.exists(), free_bytes=free, planning_required_bytes=required,
                  planning_policy=config["disk_planning_policy"], archive=archive,
                  resource_lock_present=(ROOT / config["resource_lock"]).exists(),
                  note="Presence checks only; run verifies every locked resource before annotation.")
    write_json(ROOT / config["readiness_report"], report)
    print(json.dumps({k: report[k] for k in ["ready", "blockers", "free_bytes", "planning_required_bytes"]}, indent=2))
    return report


def cache_inventory(directory):
    return [{"path": p.relative_to(directory).as_posix(), "size": p.stat().st_size, "sha256": sha256_file(p)}
            for p in sorted(Path(directory).rglob("*")) if p.is_file()]


def reference_paths(config):
    item = config["inputs"]["reference_config"]
    path = checked(ROOT / item["path"], item["sha256"])
    reference = json.loads(path.read_text(encoding="utf-8"))
    fasta = (path.parent / reference["fasta"]).resolve()
    checked(fasta, reference["fasta_sha256"])
    checked(Path(str(fasta) + ".fai"), reference["fai_sha256"])
    return fasta


def lock_resources(config, config_path):
    if config.get("annotation_mode") == "gtf":
        directory, manifest = verified_gtf(config)
        fasta = reference_paths(config)
        lock = dict(config_sha256=sha256_file(config_path), image=config["image"], platform=config["platform"],
                    annotation_mode="gtf", gtf_dir=str(directory), fasta=str(fasta),
                    gtf_manifest_sha256=sha256_file(directory / "manifest.json"), files=manifest["files"])
        write_once(ROOT / config["resource_lock"], canonical(lock))
        print("GTF resources locked:", len(manifest["files"]), "files")
        return lock
    cache = (ROOT / config["cache_dir"]).resolve()
    info = cache / "homo_sapiens" / f"{config['release']}_GRCh38" / "info.txt"
    if not info.is_file():
        raise FileNotFoundError(f"Install official release 116 GRCh38 cache first: {info}")
    fasta = reference_paths(config)
    result = docker(["run", "--rm", "--pull=never", "--network=none", "--platform", config["platform"],
                     "--mount", f"type=bind,source={cache},target=/cache,readonly", config["image"], "vep",
                     "--offline", "--cache", "--dir_cache", "/cache", "--cache_version", str(config["release"]),
                     "--species", "homo_sapiens", "--assembly", "GRCh38", "--show_cache_info"])
    if result.returncode:
        raise RuntimeError(result.stderr or result.stdout)
    inventory = cache_inventory(cache / "homo_sapiens" / "116_GRCh38")
    lock = dict(config_sha256=sha256_file(config_path), image=config["image"], platform=config["platform"],
                cache_dir=str(cache), fasta=str(fasta), cache_info=result.stdout, inventory=inventory,
                inventory_sha256=hashlib.sha256(canonical(inventory)).hexdigest())
    write_once(ROOT / config["resource_lock"], canonical(lock))
    print("Resources locked:", len(inventory), "cache files")
    return lock


def verify_resources(config, config_path):
    lock_path = ROOT / config["resource_lock"]
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    if lock["config_sha256"] != sha256_file(config_path) or lock["image"] != config["image"]:
        raise ValueError("Resource lock config/image mismatch")
    fasta = reference_paths(config)
    if config.get("annotation_mode") == "gtf":
        directory, manifest = verified_gtf(config)
        if (lock.get("annotation_mode") != "gtf" or lock["gtf_dir"] != str(directory)
                or lock["fasta"] != str(fasta) or lock["files"] != manifest["files"]
                or lock["gtf_manifest_sha256"] != sha256_file(directory / "manifest.json")):
            raise ValueError("GTF resources changed")
        return lock, sha256_file(lock_path)
    cache = (ROOT / config["cache_dir"]).resolve()
    if lock["cache_dir"] != str(cache) or lock["fasta"] != str(fasta):
        raise ValueError("Resource paths changed")
    inventory = cache_inventory(cache / "homo_sapiens" / "116_GRCh38")
    if inventory != lock["inventory"] or hashlib.sha256(canonical(inventory)).hexdigest() != lock["inventory_sha256"]:
        raise ValueError("Cache inventory changed")
    return lock, sha256_file(lock_path)


def verified_gtf(config):
    directory = (ROOT / config["gtf_dir"]).resolve()
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    if (manifest["release"] != config["release"] or manifest["assembly"] != "GRCh38"
            or manifest["source_url"] != config["gtf_source"]["url"] or manifest["image"] != config["image"]):
        raise ValueError("GTF source/release/image mismatch")
    for item in manifest["files"]:
        checked(directory / item["path"], item["sha256"])
    return directory, manifest


def annotation_flags(config):
    if config.get("annotation_mode") == "gtf":
        # VEP's --offline shortcut forces cache=1, even with --cache 0. With a
        # local GTF custom source and no --database, VEP uses local files only;
        # Docker --network=none enforces the network boundary independently.
        return [flag for flag in VEP_FLAGS if flag not in {"--offline", "--cache", "--domains", "--canonical", "--mane", "--total_length"}]
    return VEP_FLAGS


def build_command(config, prepared, chunk, output, lock):
    fasta = Path(lock["fasta"])
    gtf_mode = config.get("annotation_mode") == "gtf"
    source_mount = (f"type=bind,source={lock['gtf_dir']},target=/gtf,readonly" if gtf_mode
                    else f"type=bind,source={Path(lock['cache_dir'])},target=/cache,readonly")
    # GTF-only source, explicitly disabling cache; the container has no network.
    source_args = ["--cache", "0", "--gtf", "/gtf/annotation.gtf.gz"] if gtf_mode else ["--cache_version", str(config["release"]), "--dir_cache", "/cache"]
    return ["run", "--rm", "--pull=never", "--network=none", "--platform", config["platform"],
            "--mount", f"type=bind,source={prepared.resolve()},target=/input,readonly",
            "--mount", source_mount,
            "--mount", f"type=bind,source={fasta.parent},target=/reference,readonly",
            "--mount", f"type=bind,source={output.resolve()},target=/output", config["image"], "vep",
            *annotation_flags(config), *source_args, "--assembly", "GRCh38", "--species", "homo_sapiens",
            "--fasta", "/reference/" + fasta.name, "--synonyms", "/input/chromosome_synonyms.tsv",
            "--input_file", "/input/" + chunk["path"], "--output_file", "/output/annotations.jsonl"]


def validate_output(vcf, output):
    expected = {}
    for line in Path(vcf).read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        chrom, pos, identifier, ref, alt, *_ = line.split("\t")
        if identifier in expected:
            raise ValueError("Duplicate input ID")
        expected[identifier] = [chrom, pos, identifier, ref, alt]
    seen = set()
    with Path(output).open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            fields = record.get("input", "").split()
            if len(fields) < 5 or fields[2] not in expected or fields[:5] != expected[fields[2]]:
                raise ValueError("Output ID/coordinate/allele does not match VCF")
            if record.get("assembly_name") != "GRCh38" or fields[2] in seen:
                raise ValueError("Output assembly mismatch or duplicate ID")
            validate_vep(record, (fields[0], fields[1], fields[3], fields[4]))
            if str(record.get("id")) != fields[2] or record.get("colocated_variants"):
                raise ValueError("Output variant ID mismatch or unrequested colocated variant data")
            if not any(record.get(name) for name in ("transcript_consequences", "intergenic_consequences",
                                                     "regulatory_feature_consequences", "motif_feature_consequences")):
                raise ValueError("Output has no consequence annotation")
            seen.add(fields[2])
    if seen != set(expected):
        raise ValueError("VEP skipped or omitted input variants")
    return len(seen)


def run_batches(config, config_path, max_chunks=1):
    if max_chunks < 1:
        raise ValueError("max_chunks must be positive")
    protect(config)
    prepared, manifest = verify_bundle(config, config_path)
    lock, resource_hash = verify_resources(config, config_path)
    execution_hash = hashlib.sha256(canonical(dict(script_sha256=sha256_file(Path(__file__)), flags=annotation_flags(config),
                                                    prepared_manifest_sha256=sha256_file(prepared / "manifest.json")))).hexdigest()
    output = ROOT / config["output_dir"]
    processed, reused = 0, 0
    for chunk in manifest["chunks"]:
        chunk_dir = output / Path(chunk["path"]).stem
        receipt_path = chunk_dir / "complete.json"
        if receipt_path.exists():
            receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
            if (receipt["input_sha256"] != chunk["sha256"] or receipt["resource_sha256"] != resource_hash
                    or receipt.get("execution_sha256") != execution_hash):
                raise ValueError("Completed chunk input/resource lock differs")
            checked(chunk_dir / receipt["output"], receipt["output_sha256"])
            reused += 1
            continue
        if processed >= max_chunks:
            break
        if shutil.disk_usage(existing_parent(output)).free < config["working_reserve_bytes"]:
            raise ValueError("Insufficient working disk reserve before next chunk")
        attempt = chunk_dir / ("attempt_" + uuid4().hex)
        attempt.mkdir(parents=True)
        command = build_command(config, prepared, chunk, attempt, lock)
        write_json(attempt / "command.json", ["docker"] + command)
        print("Running", chunk["path"], flush=True)
        result = docker(command, timeout=3600)
        (attempt / "stdout.txt").write_text(result.stdout, encoding="utf-8")
        (attempt / "stderr.txt").write_text(result.stderr, encoding="utf-8")
        if result.returncode:
            raise RuntimeError(f"VEP failed; incomplete attempt retained: {attempt}")
        annotation = attempt / "annotations.jsonl"
        rows = validate_output(prepared / chunk["path"], annotation)
        protect(config)
        write_json(receipt_path, dict(created_at_utc=now(), input_sha256=chunk["sha256"], resource_sha256=resource_hash,
                                      execution_sha256=execution_hash,
                                      output=annotation.relative_to(chunk_dir).as_posix(), output_sha256=sha256_file(annotation), rows=rows))
        processed += 1
    protect(config)
    print(json.dumps(dict(new_chunks=processed, reused_chunks=reused, total_chunks=len(manifest["chunks"])), indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["prepare", "check", "lock-resources", "run"])
    parser.add_argument("--config", type=Path, default=ROOT / "config/local_vep_pipeline.json")
    parser.add_argument("--inspect-remote", action="store_true", help="HEAD request for official archive metadata only")
    parser.add_argument("--max-chunks", type=int, default=1)
    args = parser.parse_args()
    config = load_config(args.config)
    if args.command == "prepare":
        prepare(config, args.config)
    elif args.command == "check":
        if not preflight(config, args.inspect_remote)["ready"]:
            raise SystemExit(2)
    elif args.command == "lock-resources":
        lock_resources(config, args.config)
    else:
        run_batches(config, args.config, args.max_chunks)


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
