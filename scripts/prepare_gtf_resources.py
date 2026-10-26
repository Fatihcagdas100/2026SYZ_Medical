"""Small offline VEP source: Ensembl 116 GTF, tabix and transcript metadata.

The GTF contains no Pfam annotations; unavailable domains must remain missing.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import json
from pathlib import Path
import re
import shutil
import sqlite3
import sys
from urllib.request import Request, urlopen

try:
    from scripts import local_vep_pipeline as pipeline
except ModuleNotFoundError:
    import local_vep_pipeline as pipeline

PRIMARY = {str(i) for i in range(1, 23)} | {"X", "Y"}
ATTRIBUTE = re.compile(r'(\w+) "([^"]*)"')


def parse_attributes(text):
    attributes = {}
    for key, value in ATTRIBUTE.findall(text):
        attributes.setdefault(key, []).append(value)
    return attributes


def download_ranges(url, destination, total, etag, block_size=4 * 1024 * 1024):
    """Bounded, resumable HTTP ranges; reject ignored or mixed-version ranges."""
    ranges = destination.parent / "download_ranges"
    ranges.mkdir(exist_ok=True)

    def fetch(start):
        end = min(start + block_size, total) - 1
        path = ranges / f"{start:012d}.part"
        receipt = path.with_suffix(".json")
        if path.exists() and receipt.exists():
            saved = json.loads(receipt.read_text(encoding="utf-8"))
            if (saved["url"] != url or saved["etag"] != etag or saved["range"] != [start, end]
                    or path.stat().st_size != end - start + 1):
                raise ValueError("Saved GTF range identity mismatch")
            pipeline.checked(path, saved["sha256"])
            return start, path
        with urlopen(Request(url, headers={"Range": f"bytes={start}-{end}"}), timeout=40) as response:
            if (response.status != 206 or response.headers.get("Content-Range") != f"bytes {start}-{end}/{total}"
                    or response.headers.get("ETag") != etag):
                raise ValueError("Provider ignored range or changed GTF version")
            payload = response.read(end - start + 2)
        if len(payload) != end - start + 1:
            raise ValueError("Truncated GTF range")
        path.write_bytes(payload)
        pipeline.write_json(receipt, dict(url=url, etag=etag, range=[start, end], sha256=pipeline.sha256_file(path)))
        return start, path

    completed = {}
    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = [executor.submit(fetch, start) for start in range(0, total, block_size)]
        for future in as_completed(futures):
            start, path = future.result()
            completed[start] = path
            print(f"GTF ranges: {len(completed)}/{len(futures)}", flush=True)
    with destination.open("wb") as output:
        for start in sorted(completed):
            with completed[start].open("rb") as source:
                shutil.copyfileobj(source, output)


def transcript_database(source, destination):
    """Ensembl source is grouped by transcript; reject repeated groups explicitly."""
    temporary = Path(str(destination) + ".tmp")
    # Rebuild an interrupted disposable database without touching the final one.
    if temporary.exists():
        temporary.unlink()
    connection = sqlite3.connect(temporary)
    connection.execute("CREATE TABLE transcripts (id TEXT PRIMARY KEY, model TEXT NOT NULL)")
    current, model, cds_length, count = None, None, 0, 0

    def flush():
        nonlocal count
        if model is None:
            return
        if not model["Exon"]:
            raise ValueError("GTF transcript has no exons")
        complete = not {"cds_start_NF", "cds_end_NF"}.intersection(model["tags"])
        model["Translation"] = {"length": cds_length // 3} if complete and cds_length and cds_length % 3 == 0 else {}
        connection.execute("INSERT INTO transcripts VALUES (?, ?)", (current, json.dumps(model)))
        count += 1

    try:
        with gzip.open(source, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[0] not in PRIMARY or fields[2] not in {"transcript", "exon", "CDS"}:
                    continue
                attributes = parse_attributes(fields[8])
                identifier = attributes.get("transcript_id", [None])[0]
                if identifier is None:
                    raise ValueError("GTF transcript/exon/CDS lacks transcript_id")
                if identifier != current:
                    flush()
                    current, cds_length = identifier, 0
                    model = dict(id=identifier, version=int(attributes.get("transcript_version", ["0"])[0]),
                                 assembly_name="GRCh38", seq_region_name=fields[0], strand=1 if fields[6] == "+" else -1,
                                 gene_symbol=attributes.get("gene_name", [""])[0],
                                 biotype=attributes.get("transcript_biotype", ["unknown"])[0],
                                 tags=attributes.get("tag", []), Exon=[])
                if fields[2] == "exon":
                    model["Exon"].append(dict(start=int(fields[3]), end=int(fields[4])))
                elif fields[2] == "CDS":
                    cds_length += int(fields[4]) - int(fields[3]) + 1
        flush()
        connection.commit()
    finally:
        connection.close()
    temporary.replace(destination)
    return count


def setup(config, fetch=False):
    pipeline.protect(config)
    root = pipeline.ROOT / config["gtf_dir"]
    root.mkdir(parents=True, exist_ok=True)
    manifest_path = root / "manifest.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest["source_url"] != config["gtf_source"]["url"]:
            raise ValueError("GTF source URL changed")
        for item in manifest["files"]:
            pipeline.checked(root / item["path"], item["sha256"])
        print("Verified existing GTF resources:", manifest["transcripts"], "transcripts")
        return manifest
    source = root / "Homo_sapiens.GRCh38.116.gtf.gz"
    receipt_path = root / "download.receipt.json"
    if source.exists() != receipt_path.exists():
        raise ValueError("Incomplete GTF source/receipt pair")
    if not source.exists():
        if not fetch:
            raise FileNotFoundError("Use --fetch for the first public GTF download")
        if shutil.disk_usage(root).free < config["setup_planning_bytes"]:
            raise ValueError("Insufficient disk for GTF sort/build planning budget")
        temporary = source.with_name(source.name + ".part")
        with urlopen(Request(config["gtf_source"]["url"], method="HEAD"), timeout=40) as response:
            headers = dict(response.headers)
            if int(response.headers["Content-Length"]) != config["gtf_source"]["content_length_bytes"]:
                raise ValueError("GTF provider size changed")
            etag = response.headers["ETag"]
        download_ranges(config["gtf_source"]["url"], temporary, config["gtf_source"]["content_length_bytes"], etag)
        if temporary.stat().st_size != config["gtf_source"]["content_length_bytes"]:
            raise ValueError("GTF download truncated")
        temporary.replace(source)
        pipeline.write_json(receipt_path, dict(url=config["gtf_source"]["url"], retrieved_at_utc=pipeline.now(),
                                                headers=headers, sha256=pipeline.sha256_file(source)))
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    pipeline.checked(source, receipt["sha256"])
    if receipt["url"] != config["gtf_source"]["url"]:
        raise ValueError("GTF download receipt URL changed")
    # Only primary chromosomes; sorted, bgzipped, then indexed. Entire command is
    # fixed text executed inside the pinned container, with pipefail enabled.
    command = ("gzip -dc /gtf/Homo_sapiens.GRCh38.116.gtf.gz | "
               "awk 'BEGIN {FS=OFS=\"\\t\"} $1 ~ /^([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print}' | "
               "LC_ALL=C sort -k1,1 -k4,4n -k5,5n -T /gtf -S 256M | "
               "bgzip -c > /gtf/annotation.gtf.gz.tmp && "
               "tabix -f -p gff /gtf/annotation.gtf.gz.tmp")
    result = pipeline.docker(["run", "--rm", "--pull=never", "--network=none", "--platform", config["platform"],
                              "--mount", f"type=bind,source={root.resolve()},target=/gtf", config["image"],
                              "bash", "-o", "pipefail", "-c", command], timeout=1800)
    (root / "build_stdout.txt").write_text(result.stdout, encoding="utf-8")
    (root / "build_stderr.txt").write_text(result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError("GTF sort/index failed: " + result.stderr)
    (root / "annotation.gtf.gz.tmp").replace(root / "annotation.gtf.gz")
    (root / "annotation.gtf.gz.tmp.tbi").replace(root / "annotation.gtf.gz.tbi")
    count = transcript_database(source, root / "transcripts.sqlite")
    names = [source.name, receipt_path.name, "annotation.gtf.gz", "annotation.gtf.gz.tbi", "transcripts.sqlite"]
    files = [dict(path=name, size=(root / name).stat().st_size, sha256=pipeline.sha256_file(root / name)) for name in names]
    manifest = dict(created_at_utc=pipeline.now(), source_url=config["gtf_source"]["url"], release=116,
                    assembly="GRCh38", image=config["image"], transcripts=count, files=files,
                    unavailable=["pfam_domains", "interpro_domains", "pathogenicity_scores"],
                    total_bytes=sum(f["size"] for f in files), build_command=command)
    pipeline.protect(config)
    pipeline.write_json(manifest_path, manifest)
    print(json.dumps({k: manifest[k] for k in ("transcripts", "total_bytes")}, indent=2))
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=pipeline.ROOT / "config/local_vep_gtf.json")
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    setup(pipeline.load_config(args.config), args.fetch)


if __name__ == "__main__":
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    main()
