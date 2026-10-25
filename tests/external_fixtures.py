"""Gerçek veri/artifact kullanmadan dış hat entegrasyon testleri için küçük girdiler."""

import csv
import gzip
import json
from pathlib import Path

from scripts.build_baseline_features import sha256_file


def write_reference(root: Path, sequence: str = "CAAAAACGTACGTACGTACGTACGT" * 20) -> Path:
    fasta = root / "reference.fa"
    width = 7
    header = b">chr1\r\n"
    body = b"\r\n".join(sequence[i:i + width].encode("ascii") for i in range(0, len(sequence), width)) + b"\r\n"
    fasta.write_bytes(header + body)
    fai = Path(str(fasta) + ".fai")
    fai.write_text(f"chr1\t{len(sequence)}\t{len(header)}\t{width}\t{width + 2}\n", encoding="utf-8")
    config_path = root / "reference.json"
    config_path.write_text(json.dumps({
        "assembly": "GRCh38", "fasta": "reference.fa",
        "fasta_sha256": sha256_file(fasta), "fai_sha256": sha256_file(fai),
    }), encoding="utf-8")
    return config_path


def write_rows(path: Path, rows: list[dict]) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_row(variation_id, position, target, *, reference="A", alternate="G", gene="DEV", assembly="GRCh38") -> dict:
    return {
        "VariationID": str(variation_id), "#AlleleID": str(10000 + int(variation_id)),
        "Assembly": assembly, "GeneSymbol": gene,
        "Type": "single nucleotide variant",
        "Name": "NM_1(G):c.1A>G (p.Arg1Gly)" if target else "NM_1(G):c.1A>G (p.Arg1Cys)",
        "Chromosome": "1", "PositionVCF": str(position),
        "ReferenceAlleleVCF": reference, "AlternateAlleleVCF": alternate,
        "target_label": "Pathogenic" if target else "Benign", "target": str(target),
    }
