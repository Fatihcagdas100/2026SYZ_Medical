"""İlk temiz ClinVar ikili alt kümesini iki geçişte oluştur."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sqlite3
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from itertools import islice
from pathlib import Path
from typing import Iterable, Iterator, Sequence


DEFAULT_INPUT = Path("data/raw/variant_summary.txt.gz")
DEFAULT_OUTPUT = Path(
    "data/processed/clinvar_grch38_germline_exact_binary.tsv.gz"
)
DEFAULT_REPORT_DIR = Path("reports/generated/initial_binary_dataset")
TARGET_ASSEMBLY = "GRCh38"
TARGET_ORIGIN = "germline"
LABEL_TO_TARGET = {"Benign": "0", "Pathogenic": "1"}

OUTPUT_SOURCE_COLUMNS = (
    "VariationID",
    "#AlleleID",
    "Type",
    "Name",
    "GeneID",
    "GeneSymbol",
    "HGNC_ID",
    "Assembly",
    "ChromosomeAccession",
    "Chromosome",
    "Start",
    "Stop",
    "PositionVCF",
    "ReferenceAlleleVCF",
    "AlternateAlleleVCF",
    "OriginSimple",
)
OUTPUT_COLUMNS = (*OUTPUT_SOURCE_COLUMNS, "target_label", "target")


def batched(rows: Iterable[list[str]], size: int) -> Iterator[list[list[str]]]:
    iterator = iter(rows)
    while batch := list(islice(iterator, size)):
        yield batch


def is_candidate(row: Sequence[str], indices: dict[str, int]) -> bool:
    return (
        row[indices["Assembly"]] == TARGET_ASSEMBLY
        and row[indices["OriginSimple"]] == TARGET_ORIGIN
        and row[indices["ClinicalSignificance"]] in LABEL_TO_TARGET
    )


def create_index(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous = OFF;
        PRAGMA temp_store = MEMORY;
        PRAGMA cache_size = -65536;

        CREATE TABLE candidates (
            variation_id TEXT NOT NULL,
            label TEXT NOT NULL,
            output_signature TEXT NOT NULL,
            row_count INTEGER NOT NULL,
            PRIMARY KEY (variation_id, label, output_signature)
        ) WITHOUT ROWID;
        """
    )


def update_index(
    connection: sqlite3.Connection,
    occurrences: Counter[tuple[str, str, str]],
) -> None:
    connection.executemany(
        """
        INSERT INTO candidates (variation_id, label, output_signature, row_count)
        VALUES (?, ?, ?, ?)
        ON CONFLICT(variation_id, label, output_signature) DO UPDATE
        SET row_count = row_count + excluded.row_count
        """,
        (
            (variation_id, label, signature, count)
            for (variation_id, label, signature), count in occurrences.items()
        ),
    )
    connection.commit()


def read_index_statistics(connection: sqlite3.Connection) -> dict[str, int]:
    unique_ids, repeated_ids, conflicting_ids, ambiguous_ids, usable_ids = connection.execute(
        """
        SELECT
            COUNT(*),
            COALESCE(SUM(CASE WHEN total_rows > 1 THEN 1 ELSE 0 END), 0),
            COALESCE(SUM(CASE WHEN label_count > 1 THEN 1 ELSE 0 END), 0)
            ,COALESCE(SUM(CASE WHEN signature_count > 1 THEN 1 ELSE 0 END), 0)
            ,COALESCE(SUM(CASE WHEN label_count = 1 AND signature_count = 1 THEN 1 ELSE 0 END), 0)
        FROM (
            SELECT
                variation_id,
                SUM(row_count) AS total_rows,
                COUNT(DISTINCT label) AS label_count,
                COUNT(DISTINCT output_signature) AS signature_count
            FROM candidates
            GROUP BY variation_id
        )
        """
    ).fetchone()
    return {
        "unique_candidate_ids": unique_ids,
        "ids_with_multiple_candidate_rows": repeated_ids,
        "ids_with_conflicting_exact_labels": conflicting_ids,
        "ids_with_multiple_distinct_output_rows": ambiguous_ids,
        "expected_output_rows": usable_ids,
    }


def query_id_set(connection: sqlite3.Connection, distinct_column: str) -> set[str]:
    if distinct_column not in {"label", "output_signature"}:
        raise ValueError("Desteklenmeyen ayrım sütunu.")
    rows = connection.execute(
        f"""
        SELECT variation_id
        FROM candidates
        GROUP BY variation_id
        HAVING COUNT(DISTINCT {distinct_column}) > 1
        """
    )
    return {row[0] for row in rows}


def query_duplicate_ids(connection: sqlite3.Connection) -> set[str]:
    rows = connection.execute(
        """
        SELECT variation_id
        FROM candidates
        GROUP BY variation_id
        HAVING SUM(row_count) > 1
        """
    )
    return {row[0] for row in rows}


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def first_pass(
    input_path: Path,
    connection: sqlite3.Connection,
    *,
    chunk_size: int,
    progress_every: int,
) -> tuple[list[str], dict[str, int], dict[str, int]]:
    counters = Counter()
    next_progress = progress_every
    with gzip.open(input_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValueError("Girdi dosyası boş.") from exc

        required = {
            *OUTPUT_SOURCE_COLUMNS,
            "ClinicalSignificance",
        }
        missing = sorted(required.difference(header))
        if missing:
            raise ValueError("Gerekli sütunlar bulunamadı: " + ", ".join(missing))
        indices = {name: header.index(name) for name in required}

        for chunk in batched(reader, chunk_size):
            occurrences: Counter[tuple[str, str, str]] = Counter()
            for row in chunk:
                counters["rows_scanned_pass_1"] += 1
                if len(row) != len(header):
                    counters["malformed_rows_pass_1"] += 1
                    continue
                if not is_candidate(row, indices):
                    continue
                counters["candidate_rows"] += 1
                variation_id = row[indices["VariationID"]]
                if not variation_id or variation_id == "-":
                    counters["candidate_rows_without_variation_id"] += 1
                    continue
                label = row[indices["ClinicalSignificance"]]
                output_values = tuple(row[indices[name]] for name in OUTPUT_SOURCE_COLUMNS)
                signature = hashlib.sha256(
                    "\x1f".join(output_values).encode("utf-8")
                ).hexdigest()
                occurrences[(variation_id, label, signature)] += 1
            if occurrences:
                update_index(connection, occurrences)
            if progress_every > 0 and counters["rows_scanned_pass_1"] >= next_progress:
                print(
                    f"1. geçiş: {counters['rows_scanned_pass_1']:,} satır tarandı.",
                    flush=True,
                )
                while next_progress <= counters["rows_scanned_pass_1"]:
                    next_progress += progress_every
    return header, indices, dict(counters)


def second_pass(
    input_path: Path,
    output_path: Path,
    report_dir: Path,
    header: Sequence[str],
    indices: dict[str, int],
    *,
    conflicting_ids: set[str],
    ambiguous_ids: set[str],
    duplicate_ids: set[str],
    chunk_size: int,
    progress_every: int,
) -> tuple[dict[str, object], Path]:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    temporary_output = Path(
        tempfile.NamedTemporaryFile(
            prefix="binary_subset_",
            suffix=".tsv.gz.tmp",
            dir=output_path.parent,
            delete=False,
        ).name
    )

    counters = Counter()
    output_by_label: Counter[str] = Counter()
    seen_ids: set[str] = set()
    first_duplicate_signature: dict[str, tuple[str, ...]] = {}
    duplicate_differences: defaultdict[str, set[str]] = defaultdict(set)
    next_progress = progress_every
    duplicate_report_path = report_dir / "duplicate_candidate_rows.tsv"

    try:
        with (
            gzip.open(input_path, "rt", encoding="utf-8", newline="") as source,
            gzip.open(temporary_output, "wt", encoding="utf-8", newline="") as target,
            duplicate_report_path.open("w", encoding="utf-8", newline="") as duplicate_handle,
        ):
            reader = csv.reader(source, delimiter="\t")
            read_header = next(reader)
            if list(header) != read_header:
                raise ValueError("İki geçiş arasında ham dosya başlığı değişti.")
            writer = csv.writer(target, delimiter="\t", lineterminator="\n")
            writer.writerow(OUTPUT_COLUMNS)
            duplicate_writer = csv.writer(
                duplicate_handle, delimiter="\t", lineterminator="\n"
            )
            duplicate_writer.writerow(["selected_as_representative", *header])

            for chunk in batched(reader, chunk_size):
                for row in chunk:
                    counters["rows_scanned_pass_2"] += 1
                    if len(row) != len(header):
                        counters["malformed_rows_pass_2"] += 1
                        continue
                    if not is_candidate(row, indices):
                        continue

                    variation_id = row[indices["VariationID"]]
                    if not variation_id or variation_id == "-":
                        continue
                    label = row[indices["ClinicalSignificance"]]
                    output_values = tuple(row[indices[name]] for name in OUTPUT_SOURCE_COLUMNS)
                    excluded = (
                        variation_id in conflicting_ids or variation_id in ambiguous_ids
                    )
                    selected = not excluded and variation_id not in seen_ids

                    if variation_id in duplicate_ids:
                        duplicate_writer.writerow([str(selected).lower(), *row])
                        if variation_id not in first_duplicate_signature:
                            first_duplicate_signature[variation_id] = output_values
                        else:
                            first_values = first_duplicate_signature[variation_id]
                            for name, first_value, current_value in zip(
                                OUTPUT_SOURCE_COLUMNS, first_values, output_values
                            ):
                                if first_value != current_value:
                                    duplicate_differences[variation_id].add(name)

                    if variation_id in conflicting_ids:
                        counters["rows_excluded_for_label_conflict"] += 1
                        continue
                    if variation_id in ambiguous_ids:
                        counters["rows_excluded_for_multiple_output_rows"] += 1
                        continue

                    if not selected:
                        counters["duplicate_rows_skipped"] += 1
                        continue
                    seen_ids.add(variation_id)
                    writer.writerow([*output_values, label, LABEL_TO_TARGET[label]])
                    counters["output_rows"] += 1
                    output_by_label[label] += 1

                if progress_every > 0 and counters["rows_scanned_pass_2"] >= next_progress:
                    print(
                        f"2. geçiş: {counters['rows_scanned_pass_2']:,} satır tarandı.",
                        flush=True,
                    )
                    while next_progress <= counters["rows_scanned_pass_2"]:
                        next_progress += progress_every

        temporary_output.replace(output_path)
    except BaseException:
        temporary_output.unlink(missing_ok=True)
        raise

    with (report_dir / "duplicate_core_differences.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        difference_writer = csv.writer(handle)
        difference_writer.writerow(["VariationID", "differing_output_columns"])
        for variation_id in sorted(duplicate_differences, key=lambda value: int(value)):
            difference_writer.writerow(
                [variation_id, "|".join(sorted(duplicate_differences[variation_id]))]
            )

    return (
        {
            **dict(counters),
            "output_rows_by_label": dict(sorted(output_by_label.items())),
            "duplicate_ids_with_output_column_differences": len(duplicate_differences),
        },
        duplicate_report_path,
    )


def build_dataset(
    input_path: Path,
    output_path: Path,
    report_dir: Path,
    *,
    chunk_size: int = 5_000,
    progress_every: int = 500_000,
) -> dict[str, object]:
    if not input_path.is_file():
        raise FileNotFoundError(f"Girdi dosyası bulunamadı: {input_path}")
    if chunk_size <= 0:
        raise ValueError("chunk_size sıfırdan büyük olmalıdır.")

    report_dir.mkdir(parents=True, exist_ok=True)
    temporary_directory: tempfile.TemporaryDirectory[str] | None = None
    connection: sqlite3.Connection | None = None
    try:
        temporary_directory = tempfile.TemporaryDirectory(
            prefix="candidate_index_", dir=report_dir
        )
        connection = sqlite3.connect(
            Path(temporary_directory.name) / "candidate_ids.sqlite3"
        )
        create_index(connection)
        header, indices, first_pass_counts = first_pass(
            input_path,
            connection,
            chunk_size=chunk_size,
            progress_every=progress_every,
        )
        index_stats = read_index_statistics(connection)
        conflicting_ids = query_id_set(connection, "label")
        ambiguous_ids = query_id_set(connection, "output_signature")
        duplicate_ids = query_duplicate_ids(connection)

        second_pass_counts, duplicate_report_path = second_pass(
            input_path,
            output_path,
            report_dir,
            header,
            indices,
            conflicting_ids=conflicting_ids,
            ambiguous_ids=ambiguous_ids,
            duplicate_ids=duplicate_ids,
            chunk_size=chunk_size,
            progress_every=progress_every,
        )

        if second_pass_counts["output_rows"] != index_stats["expected_output_rows"]:
            raise RuntimeError(
                "Üretilen satır sayısı beklenen benzersiz kimlik sayısıyla uyuşmuyor."
            )

        manifest = {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "source": {
                "path": str(input_path.resolve()),
                "size_bytes": input_path.stat().st_size,
                "sha256": sha256_file(input_path),
            },
            "output": {
                "path": str(output_path.resolve()),
                "size_bytes": output_path.stat().st_size,
                "sha256": sha256_file(output_path),
                "delimiter": "tab",
                "compression": "gzip",
                "columns": list(OUTPUT_COLUMNS),
            },
            "selection_policy": {
                "assembly_exact": TARGET_ASSEMBLY,
                "origin_simple_exact": TARGET_ORIGIN,
                "label_mapping": LABEL_TO_TARGET,
                "likely_and_combined_labels_included": False,
                "conflicting_ids_included": False,
                "multiple_distinct_output_rows_included": False,
                "identical_duplicate_representative_rule": "first row in the pinned source file",
                "output_signature_algorithm": "SHA-256 over unit-separator-joined output source columns",
            },
            "column_use_policy": {
                "target_columns": ["target_label", "target"],
                "never_model_features": [
                    "VariationID",
                    "#AlleleID",
                    "Assembly",
                    "OriginSimple",
                    "target_label",
                    "target",
                ],
                "requires_explicit_experiment": [
                    "Name",
                    "GeneID",
                    "GeneSymbol",
                    "HGNC_ID",
                ],
                "note": "Bu çıktı henüz doğrudan X matrisi değildir.",
            },
            "first_pass": first_pass_counts,
            "candidate_index": index_stats,
            "second_pass": second_pass_counts,
            "duplicate_report": str(duplicate_report_path.resolve()),
        }
        with (report_dir / "manifest.json").open("w", encoding="utf-8") as handle:
            json.dump(manifest, handle, ensure_ascii=False, indent=2)
            handle.write("\n")

        print("\nİlk ikili veri seti oluşturuldu")
        print(f"Çıktı satırı: {second_pass_counts['output_rows']:,}")
        print(f"Benign: {second_pass_counts['output_rows_by_label'].get('Benign', 0):,}")
        print(f"Pathogenic: {second_pass_counts['output_rows_by_label'].get('Pathogenic', 0):,}")
        print(
            "Atlanan özdeş tekrar satırı: "
            f"{second_pass_counts.get('duplicate_rows_skipped', 0):,}"
        )
        print(
            "Çok-konumlu/farklı çıktı temsili nedeniyle dışlanan ID: "
            f"{index_stats['ids_with_multiple_distinct_output_rows']:,}"
        )
        print(
            "Çıktı sütunlarında farklılık taşıyan tekrar ID: "
            f"{second_pass_counts['duplicate_ids_with_output_column_differences']:,}"
        )
        print(f"Veri dosyası: {output_path.resolve()}")
        print(f"Manifest: {(report_dir / 'manifest.json').resolve()}")
        return manifest
    finally:
        if connection is not None:
            connection.close()
        if temporary_directory is not None:
            temporary_directory.cleanup()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="İlk GRCh38 germline tam Benign/Pathogenic veri setini üretir."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--chunk-size", type=int, default=5_000)
    parser.add_argument("--progress-every", type=int, default=500_000)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    build_dataset(
        args.input,
        args.output,
        args.report_dir,
        chunk_size=args.chunk_size,
        progress_every=args.progress_every,
    )


if __name__ == "__main__":
    main()

