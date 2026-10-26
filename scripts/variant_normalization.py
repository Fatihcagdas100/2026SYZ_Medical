"""Sabitlenmiş, FAI indeksli FASTA ile iki allelli küçük varyant normalizasyonu."""

from __future__ import annotations

import json
from pathlib import Path

try:
    from scripts.build_baseline_features import sha256_file
except ModuleNotFoundError:
    from build_baseline_features import sha256_file


class InvalidVariant(ValueError):
    """Kabul edilemeyen genom/allel veya referans uyuşmazlığı."""


def chromosome_alias(value: str) -> str:
    value = value.strip()
    if value.lower().startswith("chr"):
        value = value[3:]
    return "MT" if value.upper() in {"M", "MT"} else value


class ReferenceGenome:
    """Uncompressed FASTA'dan FAI byte-offset hesabıyla küçük aralıklar oku.

    Config assembly, fasta, fasta_sha256 ve fai_sha256 alanlarını içermelidir.
    FASTA yolu config konumuna göredir; dosyaların hash'leri başta doğrulanır.
    """

    def __init__(self, config_path: Path):
        config = json.loads(config_path.read_text(encoding="utf-8"))
        if config.get("assembly") != "GRCh38":
            raise ValueError("Referans config assembly GRCh38 olmalıdır.")
        fasta = Path(config["fasta"])
        if not fasta.is_absolute():
            fasta = config_path.parent / fasta
        fai = Path(str(fasta) + ".fai")
        if fasta.suffix == ".gz":
            raise ValueError("Referans FASTA sıkıştırılmamış olmalıdır.")
        for path, key in ((fasta, "fasta_sha256"), (fai, "fai_sha256")):
            if sha256_file(path) != config.get(key):
                raise ValueError(f"Referans hash uyuşmazlığı: {path}")
        self.contigs: dict[str, tuple[int, int, int, int]] = {}
        for line in fai.read_text(encoding="utf-8").splitlines():
            name, length, offset, line_bases, line_bytes, *_ = line.split("\t")
            key = chromosome_alias(name)
            record = tuple(map(int, (length, offset, line_bases, line_bytes)))
            if key in self.contigs or record[0] <= 0 or record[1] < 0 or record[2] <= 0 or record[3] < record[2]:
                raise ValueError("FAI indeksi geçersiz veya kromozom alias'ları belirsiz.")
            self.contigs[key] = record
        self.metadata = {
            "assembly": "GRCh38", "config_sha256": sha256_file(config_path),
            "fasta_sha256": config["fasta_sha256"], "fai_sha256": config["fai_sha256"],
        }
        self.primary_accessions = config.get("primary_accession_to_chromosome", {})
        if self.primary_accessions:
            alias_path = config_path.parent / config["chrom_alias_file"]
            if sha256_file(alias_path) != config.get("chrom_alias_sha256"):
                raise ValueError("Referans chromAlias hash uyuşmazlığı.")
            self.metadata["chrom_alias_sha256"] = config["chrom_alias_sha256"]
        self.handle = fasta.open("rb")

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.handle.close()

    def fetch(self, chromosome: str, start: int, length: int) -> str:
        """1-based başlangıç; aralık sınırı ve baz sayısını doğrula."""
        if chromosome not in self.contigs:
            raise InvalidVariant("unknown_contig")
        total, offset, bases, width = self.contigs[chromosome]
        if start < 1 or length < 1 or start + length - 1 > total:
            raise InvalidVariant("coordinate_out_of_range")
        first = offset + ((start - 1) // bases) * width + (start - 1) % bases
        last = offset + ((start + length - 2) // bases) * width + (start + length - 2) % bases
        self.handle.seek(first)
        sequence = self.handle.read(last - first + 1).replace(b"\n", b"").replace(b"\r", b"").decode("ascii").upper()
        if len(sequence) != length:
            raise InvalidVariant("invalid_fasta_index")
        return sequence

    def normalize(self, row: dict[str, str]) -> tuple[str, str, str, str]:
        if row.get("Assembly", "").strip() != "GRCh38":
            raise InvalidVariant("assembly_not_grch38")
        chromosome = chromosome_alias(row.get("Chromosome", ""))
        try:
            position = int(row.get("PositionVCF", ""))
        except ValueError:
            raise InvalidVariant("invalid_position") from None
        reference = row.get("ReferenceAlleleVCF", "").strip().upper()
        alternate = row.get("AlternateAlleleVCF", "").strip().upper()
        if not reference or not alternate or set(reference + alternate) - set("ACGT") or reference == alternate:
            raise InvalidVariant("unsupported_or_missing_alleles")
        if position < 1:
            raise InvalidVariant("invalid_position")
        accession = row.get("ChromosomeAccession", "").strip()
        if accession and accession not in {"-", "na"}:
            expected_chromosome = self.primary_accessions.get(accession)
            if expected_chromosome is None:
                raise InvalidVariant("nonprimary_or_unmapped_reference_accession")
            if expected_chromosome != chromosome:
                raise InvalidVariant("accession_chromosome_mismatch")
        if self.fetch(chromosome, position, len(reference)) != reference:
            raise InvalidVariant("reference_allele_mismatch")
        # Ortak son eki kaldır; bir allel boşalacaksa önce soldan referansla uzat.
        # Bu döngü tekrar bölgelerindeki eşdeğer indelleri en sola taşır.
        while reference[-1] == alternate[-1]:
            if min(len(reference), len(alternate)) == 1:
                if position == 1:
                    break
                preceding = self.fetch(chromosome, position - 1, 1)
                if preceding not in "ACGT":
                    raise InvalidVariant("ambiguous_reference_base")
                reference = preceding + reference
                alternate = preceding + alternate
                position -= 1
            reference, alternate = reference[:-1], alternate[:-1]
        while min(len(reference), len(alternate)) > 1 and reference[0] == alternate[0]:
            reference, alternate = reference[1:], alternate[1:]
            position += 1
        return chromosome, str(position), reference, alternate
