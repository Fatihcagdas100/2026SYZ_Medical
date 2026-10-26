"""UCSC initial hg38 FASTA'yı doğrula, FAI indeksle ve hash config'i oluştur."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import mmap
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

try:
    from scripts.build_baseline_features import sha256_file
except ModuleNotFoundError:
    from build_baseline_features import sha256_file


BASE = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/"


def add_accession_config(output_dir: Path, config: dict, md5s: dict) -> dict:
    """Çekirdek referans config'ini koru; accession eşlemeli yeni config üret."""
    path = output_dir / "hg38.chromAlias.txt"
    if not path.exists():
        with urlopen(BASE + path.name, timeout=45) as response:
            path.write_bytes(response.read())
    if hashlib.md5(path.read_bytes()).hexdigest() != md5s[path.name]:
        raise ValueError("chromAlias MD5 uyuşmazlığı.")
    aliases = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        contig, *names = line.split("\t")
        if contig not in {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}:
            continue
        for name in names:
            if name.startswith("NC_"):
                aliases[name] = "MT" if contig == "chrM" else contig[3:]
    updated = {**config, "primary_accession_to_chromosome": aliases,
               "chrom_alias_sha256": sha256_file(path), "chrom_alias_file": path.name}
    target = output_dir / "grch38_with_aliases.json"
    if target.exists():
        if json.loads(target.read_text(encoding="utf-8")) != updated:
            raise ValueError("Mevcut accession config'i uyuşmuyor.")
    else:
        target.write_text(json.dumps(updated, indent=2) + "\n", encoding="utf-8")
    print(f"Accession denetimli referans config'i: {target.resolve()}")
    return updated


def build_fai(fasta: Path, expected_sizes: dict[str, int]) -> Path:
    """Sağlayıcının düzenli satırlı FASTA'sını indeksle; uzunlukları ayrı dosyayla doğrula."""
    records = []
    with fasta.open("rb") as handle, mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ) as data:
        header = 0
        while header < len(data):
            if data[header:header + 1] != b">":
                raise ValueError("FASTA başlığı bulunamadı.")
            first = data.find(b"\n", header) + 1
            name = data[header + 1:first].strip().split()[0].decode("ascii")
            next_marker = data.find(b"\n>", first)
            end = next_marker + 1 if next_marker >= 0 else len(data)
            line_end = data.find(b"\n", first, end)
            width = (line_end + 1 - first) if line_end >= 0 else end - first
            bases = len(data[first:first + width].rstrip(b"\r\n"))
            if not bases or width < bases:
                raise ValueError("FASTA satır genişliği geçersiz.")
            full, remainder = divmod(end - first, width)
            length = full * bases
            if remainder:
                length += len(data[end - remainder:end].rstrip(b"\r\n"))
            if expected_sizes.get(name) != length:
                raise ValueError(f"FASTA uzunluğu sağlayıcı chrom.sizes ile eşleşmiyor: {name}")
            records.append((name, length, first, bases, width))
            header = end
    if {record[0] for record in records} != set(expected_sizes):
        raise ValueError("Referans kromozom listesi chrom.sizes ile eşleşmiyor.")
    fai = Path(str(fasta) + ".fai")
    fai.write_text("".join("\t".join(map(str, record)) + "\n" for record in records), encoding="utf-8")
    return fai


def prepare_reference(output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    config_path = output_dir / "grch38.json"
    if config_path.exists():
        config = json.loads(config_path.read_text(encoding="utf-8"))
        fasta = output_dir / config["fasta"]
        if sha256_file(fasta) != config["fasta_sha256"] or sha256_file(Path(str(fasta) + ".fai")) != config["fai_sha256"]:
            raise ValueError("Mevcut referans hash'leri uyuşmuyor.")
        print(f"Sabitlenmiş referans doğrulandı: {config_path.resolve()}")
        md5s = {line.split()[1]: line.split()[0] for line in (output_dir / "md5sum.txt").read_text().splitlines() if line.strip()}
        return add_accession_config(output_dir, config, md5s)
    with urlopen(BASE + "md5sum.txt", timeout=45) as response:
        md5_bytes = response.read()
    (output_dir / "md5sum.txt").write_bytes(md5_bytes)
    md5s = {line.split()[1]: line.split()[0] for line in md5_bytes.decode().splitlines() if line.strip()}
    compressed = output_dir / "hg38.fa.gz"
    expected_md5 = md5s["hg38.fa.gz"]
    if not compressed.exists():
        temporary = output_dir / "hg38.fa.gz.part"
        if temporary.exists():
            raise FileExistsError(f"Tamamlanmamış indirme mevcut; inceleyin: {temporary}")
        with urlopen(BASE + "hg38.fa.gz", timeout=45) as response, temporary.open("xb") as handle:
            total = 0
            next_progress = 64 * 1024 * 1024
            while chunk := response.read(8 * 1024 * 1024):
                handle.write(chunk)
                total += len(chunk)
                if total >= next_progress:
                    print(f"GRCh38 indirme: {total / 1024**2:.0f} MiB", flush=True)
                    next_progress += 64 * 1024 * 1024
        temporary.replace(compressed)
    digest = hashlib.md5()
    with compressed.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    if digest.hexdigest() != expected_md5:
        raise ValueError("UCSC FASTA MD5 doğrulaması başarısız.")
    print("UCSC MD5 doğrulandı; FASTA açılıyor.", flush=True)
    fasta = output_dir / "hg38.fa"
    if fasta.exists():
        raise FileExistsError(f"Config olmadan mevcut FASTA üzerine yazılmayacak: {fasta}")
    with gzip.open(compressed, "rb") as source, fasta.open("xb") as target:
        shutil.copyfileobj(source, target, length=8 * 1024 * 1024)
    with urlopen(BASE + "hg38.chrom.sizes", timeout=45) as response:
        sizes_bytes = response.read()
    (output_dir / "hg38.chrom.sizes").write_bytes(sizes_bytes)
    if hashlib.md5(sizes_bytes).hexdigest() != md5s["hg38.chrom.sizes"]:
        raise ValueError("chrom.sizes MD5 uyuşmazlığı.")
    sizes = {name: int(size) for name, size in (line.split() for line in sizes_bytes.decode().splitlines())}
    fai = build_fai(fasta, sizes)
    config = {
        "assembly": "GRCh38", "fasta": fasta.name,
        "fasta_sha256": sha256_file(fasta), "fai_sha256": sha256_file(fai),
        "source_url": BASE + "hg38.fa.gz", "provider_release": "UCSC hg38 initial / GCA_000001405.15",
        "provider_md5": expected_md5, "md5_list_sha256": sha256_file(output_dir / "md5sum.txt"),
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    print(f"GRCh38 referansı hazır: {config_path.resolve()}")
    return add_accession_config(output_dir, config, md5s)


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description="Sağlayıcı MD5 denetimli GRCh38 FASTA/FAI hazırlar.")
    parser.add_argument("--output-dir", type=Path, default=Path("data/reference/ucsc_hg38_initial"))
    args = parser.parse_args()
    prepare_reference(args.output_dir)


if __name__ == "__main__":
    main()
