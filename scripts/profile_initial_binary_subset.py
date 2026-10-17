"""İlk ClinVar ikili alt kümesi için filtre hunisini profille.

Bu adım veri seti üretmez. Ham dosyayı parçalar halinde okuyarak seçilecek
satırları, etiket dağılımını ve VariationID tekrarlarını sayar.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sqlite3
import sys
import tempfile
from collections import Counter
from datetime import datetime, timezone
from itertools import islice
from pathlib import Path
from typing import Iterable, Iterator


DEFAULT_INPUT = Path("data/raw/variant_summary.txt.gz")
DEFAULT_OUTPUT_DIR = Path("reports/generated/initial_binary_subset")
TARGET_ASSEMBLY = "GRCh38"
TARGET_ORIGIN = "germline"
TARGET_LABELS = ("Benign", "Pathogenic")
MISSINGNESS_FIELDS = ("PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF")


def batched(rows: Iterable[list[str]], size: int) -> Iterator[list[list[str]]]:
    iterator = iter(rows)
    while batch := list(islice(iterator, size)):
        yield batch


def create_candidate_index(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous = OFF;
        PRAGMA temp_store = MEMORY;
        PRAGMA cache_size = -65536;

        CREATE TABLE candidates (
            variation_id TEXT NOT NULL,
            label TEXT NOT NULL,
            row_count INTEGER NOT NULL,
            PRIMARY KEY (variation_id, label)
        ) WITHOUT ROWID;
        """
    )


def update_candidate_index(
    connection: sqlite3.Connection,
    occurrences: Counter[tuple[str, str]],
) -> None:
    connection.executemany(
        """
        INSERT INTO candidates (variation_id, label, row_count)
        VALUES (?, ?, ?)
        ON CONFLICT(variation_id, label) DO UPDATE
        SET row_count = row_count + excluded.row_count
        """,
        (
            (variation_id, label, count)
            for (variation_id, label), count in occurrences.items()
        ),
    )
    connection.commit()


def candidate_statistics(connection: sqlite3.Connection) -> dict[str, object]:
    unique_ids, repeated_ids, conflicting_ids = connection.execute(
        """
        SELECT
            COUNT(*),
            COALESCE(SUM(CASE WHEN total_rows > 1 THEN 1 ELSE 0 END), 0),
            COALESCE(SUM(CASE WHEN label_count > 1 THEN 1 ELSE 0 END), 0)
        FROM (
            SELECT
                variation_id,
                SUM(row_count) AS total_rows,
                COUNT(*) AS label_count
            FROM candidates
            GROUP BY variation_id
        )
        """
    ).fetchone()

    usable_by_label = dict(
        connection.execute(
            """
            SELECT label, COUNT(*)
            FROM (
                SELECT variation_id, MIN(label) AS label
                FROM candidates
                GROUP BY variation_id
                HAVING COUNT(*) = 1
            )
            GROUP BY label
            ORDER BY label
            """
        )
    )
    return {
        "unique_candidate_variation_ids": unique_ids,
        "candidate_ids_with_multiple_rows": repeated_ids,
        "candidate_ids_with_conflicting_exact_labels": conflicting_ids,
        "usable_nonconflicting_unique_ids": unique_ids - conflicting_ids,
        "usable_unique_ids_by_label": usable_by_label,
    }


def write_counter_csv(
    output_path: Path,
    header: list[str],
    rows: Iterable[tuple[object, ...]],
) -> None:
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def profile_subset(
    input_path: Path,
    output_dir: Path,
    *,
    chunk_size: int = 5_000,
    max_rows: int | None = None,
    progress_every: int = 250_000,
) -> dict[str, object]:
    if not input_path.is_file():
        raise FileNotFoundError(f"Girdi dosyası bulunamadı: {input_path}")
    if chunk_size <= 0:
        raise ValueError("chunk_size sıfırdan büyük olmalıdır.")
    if max_rows is not None and max_rows <= 0:
        raise ValueError("max_rows sıfırdan büyük olmalıdır.")

    output_dir.mkdir(parents=True, exist_ok=True)
    temporary_directory: tempfile.TemporaryDirectory[str] | None = None
    connection: sqlite3.Connection | None = None

    try:
        temporary_directory = tempfile.TemporaryDirectory(
            prefix="candidate_index_", dir=output_dir
        )
        connection = sqlite3.connect(
            Path(temporary_directory.name) / "candidate_ids.sqlite3"
        )
        create_candidate_index(connection)

        counts = Counter()
        germline_labels: Counter[str] = Counter()
        candidate_rows_by_label: Counter[str] = Counter()
        review_status_by_label: Counter[tuple[str, str]] = Counter()
        variant_type_by_label: Counter[tuple[str, str]] = Counter()
        missingness: Counter[tuple[str, str, str]] = Counter()
        next_progress = progress_every
        limit_reached = False

        with gzip.open(input_path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValueError("Girdi dosyası boş.") from exc

            required = {
                "VariationID",
                "Assembly",
                "OriginSimple",
                "ClinicalSignificance",
                "ReviewStatus",
                "Type",
                *MISSINGNESS_FIELDS,
            }
            missing_columns = sorted(required.difference(header))
            if missing_columns:
                raise ValueError(
                    "Gerekli sütunlar bulunamadı: " + ", ".join(missing_columns)
                )
            indices = {name: header.index(name) for name in required}

            for chunk in batched(reader, chunk_size):
                chunk_occurrences: Counter[tuple[str, str]] = Counter()
                for row in chunk:
                    if max_rows is not None and counts["total_rows"] >= max_rows:
                        limit_reached = True
                        break

                    counts["total_rows"] += 1
                    if len(row) != len(header):
                        counts["malformed_rows"] += 1
                        continue

                    if row[indices["Assembly"]] != TARGET_ASSEMBLY:
                        continue
                    counts["grch38_rows"] += 1

                    if row[indices["OriginSimple"]] != TARGET_ORIGIN:
                        continue
                    counts["grch38_germline_rows"] += 1

                    label = row[indices["ClinicalSignificance"]]
                    germline_labels[label] += 1
                    if label not in TARGET_LABELS:
                        continue
                    counts["exact_binary_candidate_rows"] += 1
                    candidate_rows_by_label[label] += 1

                    variation_id = row[indices["VariationID"]]
                    if not variation_id or variation_id == "-":
                        counts["candidate_rows_without_variation_id"] += 1
                        continue
                    counts["candidate_rows_with_variation_id"] += 1
                    chunk_occurrences[(variation_id, label)] += 1

                    review_status_by_label[
                        (label, row[indices["ReviewStatus"]])
                    ] += 1
                    variant_type_by_label[(label, row[indices["Type"]])] += 1
                    for field in MISSINGNESS_FIELDS:
                        value = row[indices[field]]
                        if value == "":
                            missingness[(label, field, "empty")] += 1
                        elif value == "-":
                            missingness[(label, field, "dash")] += 1
                        elif len(value) == 2 and value.lower() == "na":
                            missingness[(label, field, "na_literal")] += 1

                if chunk_occurrences:
                    update_candidate_index(connection, chunk_occurrences)

                if progress_every > 0 and counts["total_rows"] >= next_progress:
                    print(
                        f"İlerleme: {counts['total_rows']:,} satır tarandı.",
                        flush=True,
                    )
                    while next_progress <= counts["total_rows"]:
                        next_progress += progress_every
                if limit_reached:
                    break

        stats = candidate_statistics(connection)
        scan_complete = not limit_reached

        funnel = [
            ("all_rows", "row", counts["total_rows"]),
            ("assembly_grch38", "row", counts["grch38_rows"]),
            ("grch38_and_origin_germline", "row", counts["grch38_germline_rows"]),
            ("exact_pathogenic_or_benign", "row", counts["exact_binary_candidate_rows"]),
            ("candidate_rows_with_variation_id", "row", counts["candidate_rows_with_variation_id"]),
            ("unique_candidate_variation_ids", "id", stats["unique_candidate_variation_ids"]),
            ("usable_nonconflicting_unique_ids", "id", stats["usable_nonconflicting_unique_ids"]),
        ]
        previous_by_unit: dict[str, int] = {}
        funnel_rows = []
        for step, unit, count in funnel:
            previous = previous_by_unit.get(unit)
            percent_previous = 100 * count / previous if previous else 100.0
            percent_total = 100 * count / counts["total_rows"] if counts["total_rows"] else 0
            funnel_rows.append((step, unit, count, f"{percent_previous:.6f}", f"{percent_total:.6f}"))
            previous_by_unit[unit] = count

        write_counter_csv(
            output_dir / "filter_funnel.csv",
            ["step", "unit", "count", "percentage_of_previous_same_unit", "percentage_of_all_rows"],
            funnel_rows,
        )
        write_counter_csv(
            output_dir / "germline_label_distribution.csv",
            ["ClinicalSignificance", "row_count"],
            sorted(germline_labels.items(), key=lambda item: (-item[1], item[0])),
        )
        write_counter_csv(
            output_dir / "candidate_rows_by_label.csv",
            ["label", "row_count"],
            sorted(candidate_rows_by_label.items()),
        )
        write_counter_csv(
            output_dir / "review_status_by_label.csv",
            ["label", "ReviewStatus", "row_count"],
            (
                (label, status, count)
                for (label, status), count in sorted(
                    review_status_by_label.items(),
                    key=lambda item: (item[0][0], -item[1], item[0][1]),
                )
            ),
        )
        write_counter_csv(
            output_dir / "variant_type_by_label.csv",
            ["label", "Type", "row_count"],
            (
                (label, variant_type, count)
                for (label, variant_type), count in sorted(
                    variant_type_by_label.items(),
                    key=lambda item: (item[0][0], -item[1], item[0][1]),
                )
            ),
        )
        write_counter_csv(
            output_dir / "candidate_missingness.csv",
            ["label", "column", "missing_marker", "row_count"],
            (
                (label, field, marker, count)
                for (label, field, marker), count in sorted(missingness.items())
            ),
        )

        result = {
            "profiled_at_utc": datetime.now(timezone.utc).isoformat(),
            "input_file": str(input_path.resolve()),
            "scan_complete": scan_complete,
            "chunk_size": chunk_size,
            "policy": {
                "assembly": TARGET_ASSEMBLY,
                "origin_simple": TARGET_ORIGIN,
                "accepted_clinical_significance_exact": list(TARGET_LABELS),
                "likely_and_combined_labels_included": False,
                "conflicting_exact_labels_included": False,
            },
            "counts": dict(counts),
            "candidate_statistics": stats,
        }
        with (output_dir / "summary.json").open("w", encoding="utf-8") as handle:
            json.dump(result, handle, ensure_ascii=False, indent=2)
            handle.write("\n")

        print("\nİlk ikili alt küme filtre profili")
        print(f"Toplam taranan satır: {counts['total_rows']:,}")
        print(f"GRCh38 satırı: {counts['grch38_rows']:,}")
        print(f"GRCh38 + germline satırı: {counts['grch38_germline_rows']:,}")
        print(f"Tam Pathogenic/Benign aday satırı: {counts['exact_binary_candidate_rows']:,}")
        print(f"Benzersiz aday VariationID: {stats['unique_candidate_variation_ids']:,}")
        print(f"Çelişkili tam etiketli ID: {stats['candidate_ids_with_conflicting_exact_labels']:,}")
        print(f"Kullanılabilir benzersiz ID: {stats['usable_nonconflicting_unique_ids']:,}")
        print(f"Rapor klasörü: {output_dir.resolve()}")
        return result
    finally:
        if connection is not None:
            connection.close()
        if temporary_directory is not None:
            temporary_directory.cleanup()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="İlk GRCh38 germline Pathogenic/Benign alt kümesini profiller."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--chunk-size", type=int, default=5_000)
    parser.add_argument("--max-rows", type=int, default=None)
    parser.add_argument("--progress-every", type=int, default=250_000)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    profile_subset(
        args.input,
        args.output_dir,
        chunk_size=args.chunk_size,
        max_rows=args.max_rows,
        progress_every=args.progress_every,
    )


if __name__ == "__main__":
    main()

