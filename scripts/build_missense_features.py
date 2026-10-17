"""Ham HGVS protein değişiminden kimliksiz missense özellikleri üret."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import re
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

try:
    from scripts.build_hgvs_features import PROTEIN_NOTATION
except ModuleNotFoundError:  # Doğrudan ``python scripts/...`` çalıştırması.
    from build_hgvs_features import PROTEIN_NOTATION


DEFAULT_SOURCE = Path("data/processed/clinvar_grch38_germline_exact_binary.tsv.gz")
DEFAULT_FEATURES = Path("data/processed/clinvar_hgvs_features.tsv.gz")
DEFAULT_OUTPUT = Path("data/processed/clinvar_missense_features.tsv.gz")
DEFAULT_REPORT_DIR = Path("reports/generated/missense_features")

MISSENSE_CATEGORICAL_FEATURES = (
    "missense_reference_residue",
    "missense_alternate_residue",
    "missense_substitution",
    "missense_residue_class_transition",
    "missense_charge_transition",
)
MISSENSE_NUMERIC_FEATURES = (
    "missense_features_available",
    "missense_grantham_distance",
    "missense_blosum62_score",
    "missense_hydropathy_delta",
    "missense_absolute_hydropathy_delta",
)
MISSENSE_FEATURE_COLUMNS = (*MISSENSE_CATEGORICAL_FEATURES, *MISSENSE_NUMERIC_FEATURES)
MISSENSE_TOKEN = re.compile(r"^([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")

AA_THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
}
AA_ORDER = "ARNDCQEGHILKMFPSTWYV"

# Grantham (1974), Table 2; order: ARNDCQEGHILKMFPSTWYV.
GRANTHAM_ROWS = (
    (0,112,111,126,195,91,107,60,86,94,96,106,84,113,27,99,58,148,112,64),
    (112,0,86,96,180,43,54,125,29,97,102,26,91,97,103,110,71,101,77,96),
    (111,86,0,23,139,46,42,80,68,149,153,94,142,158,91,46,65,174,143,133),
    (126,96,23,0,154,61,45,94,81,168,172,101,160,177,108,65,85,181,160,152),
    (195,180,139,154,0,154,170,159,174,198,198,202,196,205,169,112,149,215,194,192),
    (91,43,46,61,154,0,29,87,24,109,113,53,101,116,76,68,42,130,99,96),
    (107,54,42,45,170,29,0,98,40,134,138,56,126,140,93,80,65,152,122,121),
    (60,125,80,94,159,87,98,0,98,135,138,127,127,153,42,56,59,184,147,109),
    (86,29,68,81,174,24,40,98,0,94,99,32,87,100,77,89,47,115,83,84),
    (94,97,149,168,198,109,134,135,94,0,5,102,10,21,95,142,89,61,33,29),
    (96,102,153,172,198,113,138,138,99,5,0,107,15,22,98,145,92,61,36,32),
    (106,26,94,101,202,53,56,127,32,102,107,0,95,102,103,121,78,110,85,97),
    (84,91,142,160,196,101,126,127,87,10,15,95,0,28,87,135,81,67,36,21),
    (113,97,158,177,205,116,140,153,100,21,22,102,28,0,114,155,103,40,22,50),
    (27,103,91,108,169,76,93,42,77,95,98,103,87,114,0,74,38,147,110,68),
    (99,110,46,65,112,68,80,56,89,142,145,121,135,155,74,0,58,177,144,124),
    (58,71,65,85,149,42,65,59,47,89,92,78,81,103,38,58,0,128,92,69),
    (148,101,174,181,215,130,152,184,115,61,61,110,67,40,147,177,128,0,37,88),
    (112,77,143,160,194,99,122,147,83,33,36,85,36,22,110,144,92,37,0,55),
    (64,96,133,152,192,96,121,109,84,29,32,97,21,50,68,124,69,88,55,0),
)

# Henikoff & Henikoff (1992), BLOSUM62; order: ARNDCQEGHILKMFPSTWYV.
BLOSUM62_ROWS = (
    (4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0),
    (-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3),
    (-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3),
    (-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3),
    (0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1),
    (-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2),
    (-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2),
    (0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3),
    (-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3),
    (-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3),
    (-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1),
    (-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2),
    (-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1),
    (-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1),
    (-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2),
    (1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2),
    (0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0),
    (-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3),
    (-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1),
    (0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4),
)

# Kyte & Doolittle (1982) amino-acid hydropathy indices.
HYDROPATHY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9,
    "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9,
    "Y": -1.3, "P": -1.6, "H": -3.2, "E": -3.5, "Q": -3.5,
    "D": -3.5, "N": -3.5, "K": -3.9, "R": -4.5,
}
RESIDUE_CLASS = {
    **{aa: "positive" for aa in "RKH"},
    **{aa: "negative" for aa in "DE"},
    **{aa: "polar_uncharged" for aa in "STNQCY"},
    **{aa: "hydrophobic" for aa in "AVILMFW"},
    **{aa: "special" for aa in "GP"},
}
CHARGE_CLASS = {
    **{aa: "positive" for aa in "RKH"},
    **{aa: "negative" for aa in "DE"},
    **{aa: "neutral" for aa in "ACFGILMNPQSTVWY"},
}


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def matrix_value(rows: tuple[tuple[int, ...], ...], reference: str, alternate: str) -> int:
    return rows[AA_ORDER.index(reference)][AA_ORDER.index(alternate)]


def derive_missense_features(name: str, protein_effect: str) -> dict[str, str | int | float]:
    sentinel: dict[str, str | int | float] = {
        "missense_reference_residue": "not_missense",
        "missense_alternate_residue": "not_missense",
        "missense_substitution": "not_missense",
        "missense_residue_class_transition": "not_missense",
        "missense_charge_transition": "not_missense",
        "missense_features_available": 0,
        "missense_grantham_distance": 0,
        "missense_blosum62_score": 0,
        "missense_hydropathy_delta": 0.0,
        "missense_absolute_hydropathy_delta": 0.0,
    }
    if protein_effect != "missense":
        return sentinel
    protein_match = PROTEIN_NOTATION.search(name)
    token_match = MISSENSE_TOKEN.fullmatch(protein_match.group(1)) if protein_match else None
    if token_match is None:
        sentinel.update(
            {
                "missense_reference_residue": "unparsed_missense",
                "missense_alternate_residue": "unparsed_missense",
                "missense_substitution": "unparsed_missense",
                "missense_residue_class_transition": "unparsed_missense",
                "missense_charge_transition": "unparsed_missense",
            }
        )
        return sentinel
    reference = AA_THREE_TO_ONE.get(token_match.group(1))
    alternate = AA_THREE_TO_ONE.get(token_match.group(3))
    if reference is None or alternate is None:
        sentinel.update(
            {
                "missense_reference_residue": "unparsed_missense",
                "missense_alternate_residue": "unparsed_missense",
                "missense_substitution": "unparsed_missense",
                "missense_residue_class_transition": "unparsed_missense",
                "missense_charge_transition": "unparsed_missense",
            }
        )
        return sentinel
    hydropathy_delta = HYDROPATHY[alternate] - HYDROPATHY[reference]
    return {
        "missense_reference_residue": reference,
        "missense_alternate_residue": alternate,
        "missense_substitution": f"{reference}>{alternate}",
        "missense_residue_class_transition": (
            f"{RESIDUE_CLASS[reference]}>{RESIDUE_CLASS[alternate]}"
        ),
        "missense_charge_transition": (
            f"{CHARGE_CLASS[reference]}>{CHARGE_CLASS[alternate]}"
        ),
        "missense_features_available": 1,
        "missense_grantham_distance": matrix_value(
            GRANTHAM_ROWS, reference, alternate
        ),
        "missense_blosum62_score": matrix_value(
            BLOSUM62_ROWS, reference, alternate
        ),
        "missense_hydropathy_delta": hydropathy_delta,
        "missense_absolute_hydropathy_delta": abs(hydropathy_delta),
    }


def build_missense_features(
    source_path: Path,
    features_path: Path,
    output_path: Path,
    report_dir: Path,
    *,
    progress_every: int = 100_000,
) -> dict[str, object]:
    for path in (source_path, features_path):
        if not path.is_file():
            raise FileNotFoundError(f"Girdi bulunamadı: {path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    temporary_handle = tempfile.NamedTemporaryFile(
        prefix="missense_features_",
        suffix=".tsv.gz.tmp",
        dir=output_path.parent,
        delete=False,
    )
    temporary_path = Path(temporary_handle.name)
    temporary_handle.close()
    counters: Counter[str] = Counter()

    try:
        with (
            gzip.open(source_path, "rt", encoding="utf-8", newline="") as source_handle,
            gzip.open(features_path, "rt", encoding="utf-8", newline="") as feature_handle,
            gzip.open(temporary_path, "wt", encoding="utf-8", newline="") as output_handle,
        ):
            source_reader = csv.DictReader(source_handle, delimiter="\t")
            feature_reader = csv.DictReader(feature_handle, delimiter="\t")
            if source_reader.fieldnames is None or feature_reader.fieldnames is None:
                raise ValueError("Girdilerden birinde başlık bulunamadı.")
            if not {"VariationID", "Name", "target"}.issubset(source_reader.fieldnames):
                raise ValueError("Ham kaynakta VariationID, Name veya target eksik.")
            if not {"VariationID", "protein_effect", "target"}.issubset(
                feature_reader.fieldnames
            ):
                raise ValueError("HGVS özellik girdisinde gerekli sütunlar eksik.")
            output_columns = [
                column for column in feature_reader.fieldnames if column != "target"
            ] + [*MISSENSE_FEATURE_COLUMNS, "target"]
            writer = csv.DictWriter(
                output_handle,
                fieldnames=output_columns,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for source_row, feature_row in itertools.zip_longest(
                source_reader, feature_reader
            ):
                if source_row is None or feature_row is None:
                    raise ValueError("Kaynak ve HGVS özellik satır sayıları eşleşmiyor.")
                if source_row["VariationID"] != feature_row["VariationID"]:
                    raise ValueError(
                        f"Satır sırası eşleşmiyor: {source_row['VariationID']} != "
                        f"{feature_row['VariationID']}"
                    )
                if source_row["target"] != feature_row["target"]:
                    raise ValueError(
                        f"Hedef uyuşmazlığı: VariationID={source_row['VariationID']}"
                    )
                derived = derive_missense_features(
                    source_row["Name"], feature_row["protein_effect"]
                )
                output_row = {
                    column: feature_row[column]
                    for column in feature_reader.fieldnames
                    if column != "target"
                }
                output_row.update(derived)
                output_row["target"] = feature_row["target"]
                writer.writerow(output_row)
                counters["rows"] += 1
                if feature_row["protein_effect"] == "missense":
                    counters["missense_rows"] += 1
                    if derived["missense_features_available"] == 1:
                        counters["parsed_missense_rows"] += 1
                    else:
                        counters["unparsed_missense_rows"] += 1
                if progress_every > 0 and counters["rows"] % progress_every == 0:
                    print(f"İlerleme: {counters['rows']:,} satır işlendi.", flush=True)
        temporary_path.replace(output_path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise

    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "source": {"path": str(source_path.resolve()), "sha256": sha256_file(source_path)},
            "hgvs_features": {"path": str(features_path.resolve()), "sha256": sha256_file(features_path)},
        },
        "output": {
            "path": str(output_path.resolve()),
            "sha256": sha256_file(output_path),
            "size_bytes": output_path.stat().st_size,
        },
        "counts": dict(counters),
        "feature_policy": {
            "features": list(MISSENSE_FEATURE_COLUMNS),
            "protein_position_in_output": False,
            "raw_name_in_output": False,
            "gene_or_transcript_identity_in_output": False,
            "target_used_to_derive_features": False,
            "non_missense_numeric_sentinel": 0,
            "non_missense_categorical_sentinel": "not_missense",
        },
        "sources": {
            "grantham": {
                "citation": "Grantham R. Science. 1974;185(4154):862-864.",
                "doi": "10.1126/science.185.4154.862",
            },
            "blosum62": {
                "citation": "Henikoff S, Henikoff JG. PNAS. 1992;89(22):10915-10919.",
                "doi": "10.1073/pnas.89.22.10915",
            },
            "hydropathy": {
                "citation": "Kyte J, Doolittle RF. J Mol Biol. 1982;157(1):105-132.",
                "doi": "10.1016/0022-2836(82)90515-0",
            },
        },
    }
    manifest_path = report_dir / "manifest.json"
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print("Missense biyolojik özellik tablosu oluşturuldu")
    print(f"Toplam satır: {counters['rows']:,}")
    print(
        f"Ayrıştırılan missense: {counters['parsed_missense_rows']:,}/"
        f"{counters['missense_rows']:,}"
    )
    print(f"Çıktı: {output_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="HGVS missense değişimlerinden kimliksiz biyolojik özellikler üretir."
    )
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--features", type=Path, default=DEFAULT_FEATURES)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--progress-every", type=int, default=100_000)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    build_missense_features(
        args.source,
        args.features,
        args.output,
        args.report_dir,
        progress_every=args.progress_every,
    )


if __name__ == "__main__":
    main()
