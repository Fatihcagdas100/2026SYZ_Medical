"""ClinVar variant_summary dosyasını belleğe yüklemeden denetle.

Betik yalnızca Python standart kütüphanesini kullanır. Gzip içindeki TSV
satırlarını parçalar halinde okur, sütun profilleri ve seçili kategorik
dağılımlar üretir. VariationID tekrarları geçici bir SQLite indeksiyle izlenir;
bu nedenle bütün kimlikleri RAM'de tutmak gerekmez.
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
from typing import Iterable, Iterator, Sequence


DEFAULT_INPUT = Path("data/raw/variant_summary.txt.gz")
DEFAULT_OUTPUT_DIR = Path("reports/generated/clinvar_audit")
DEFAULT_DISTRIBUTION_COLUMNS = (
    "Assembly",
    "Chromosome",
    "Type",
    "ClinicalSignificance",
    "ClinSigSimple",
    "OriginSimple",
    "ReviewStatus",
    "SomaticClinicalImpact",
    "Oncogenicity",
)
DEFAULT_PREVIEW_COLUMNS = (
    "VariationID",
    "#AlleleID",
    "Assembly",
    "Chromosome",
    "PositionVCF",
    "Type",
    "GeneSymbol",
    "ClinicalSignificance",
    "ReviewStatus",
)


def batched(rows: Iterable[list[str]], size: int) -> Iterator[list[list[str]]]:
    """Bir iterator'ı en fazla ``size`` satırlık listelere böl."""

    iterator = iter(rows)
    while batch := list(islice(iterator, size)):
        yield batch


def normalize_row(values: list[str], column_count: int) -> tuple[list[str], bool]:
    """Eksik/fazla alanlı satırı sütun sayısına getir ve bozukluk bilgisini dön."""

    malformed = len(values) != column_count
    if len(values) < column_count:
        values = values + [""] * (column_count - len(values))
    elif len(values) > column_count:
        values = values[:column_count]
    return values, malformed


def display_value(value: str, width: int = 42) -> str:
    value = value.replace("\r", " ").replace("\n", " ")
    return value if len(value) <= width else value[: width - 1] + "…"


def print_preview(
    header: Sequence[str], sample_rows: Sequence[Sequence[str]], requested: Sequence[str]
) -> None:
    selected = [name for name in requested if name in header]
    if not selected or not sample_rows:
        print("\nÖrnek satır bulunamadı.")
        return

    indices = [header.index(name) for name in selected]
    table = [[display_value(row[index]) for index in indices] for row in sample_rows]
    widths = [
        min(42, max(len(name), *(len(row[i]) for row in table)))
        for i, name in enumerate(selected)
    ]

    def render(values: Sequence[str]) -> str:
        return " | ".join(value.ljust(widths[i]) for i, value in enumerate(values))

    print("\nÖrnek kayıtlar (seçili sütunlar):")
    print(render(selected))
    print("-+-".join("-" * width for width in widths))
    for row in table:
        print(render(row))


def create_variation_index(connection: sqlite3.Connection) -> None:
    connection.executescript(
        """
        PRAGMA journal_mode = OFF;
        PRAGMA synchronous = OFF;
        PRAGMA temp_store = MEMORY;
        PRAGMA cache_size = -65536;

        CREATE TABLE variation_occurrences (
            variation_id TEXT NOT NULL,
            assembly TEXT NOT NULL,
            row_count INTEGER NOT NULL,
            PRIMARY KEY (variation_id, assembly)
        ) WITHOUT ROWID;
        """
    )


def update_variation_index(
    connection: sqlite3.Connection,
    variation_occurrences: Counter[tuple[str, str]],
) -> None:
    connection.executemany(
        """
        INSERT INTO variation_occurrences (variation_id, assembly, row_count)
        VALUES (?, ?, ?)
        ON CONFLICT(variation_id, assembly) DO UPDATE
        SET row_count = row_count + excluded.row_count
        """,
        (
            (variation_id, assembly, count)
            for (variation_id, assembly), count in variation_occurrences.items()
        ),
    )
    connection.commit()


def variation_statistics(connection: sqlite3.Connection) -> dict[str, int]:
    unique_ids, repeated_ids, extra_rows, multiple_assemblies = connection.execute(
        """
        SELECT
            COUNT(*),
            COALESCE(SUM(CASE WHEN total_rows > 1 THEN 1 ELSE 0 END), 0),
            COALESCE(SUM(CASE WHEN total_rows > 1 THEN total_rows - 1 ELSE 0 END), 0),
            COALESCE(SUM(CASE WHEN assembly_count > 1 THEN 1 ELSE 0 END), 0)
        FROM (
            SELECT
                variation_id,
                SUM(row_count) AS total_rows,
                COUNT(*) AS assembly_count
            FROM variation_occurrences
            GROUP BY variation_id
        )
        """
    ).fetchone()
    return {
        "unique_variation_ids": unique_ids,
        "repeated_variation_ids": repeated_ids,
        "rows_beyond_first_per_variation_id": extra_rows,
        "variation_ids_with_multiple_assemblies": multiple_assemblies,
    }


def write_repeat_examples(
    connection: sqlite3.Connection, output_path: Path, limit: int = 100
) -> None:
    query = """
        SELECT
            variation_id,
            SUM(row_count) AS total_rows,
            COUNT(*) AS assembly_count,
            GROUP_CONCAT(assembly, '|') AS assemblies
        FROM variation_occurrences
        GROUP BY variation_id
        HAVING SUM(row_count) > 1
        ORDER BY total_rows DESC, variation_id
        LIMIT ?
    """
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["VariationID", "row_count", "assembly_count", "assemblies"]
        )
        writer.writerows(connection.execute(query, (limit,)))


def write_reports(
    *,
    output_dir: Path,
    input_path: Path,
    header: Sequence[str],
    sample_rows: Sequence[Sequence[str]],
    profiles: Sequence[Sequence[int]],
    distributions: dict[str, Counter[str]],
    row_count: int,
    malformed_rows: int,
    scan_complete: bool,
    chunk_size: int,
    variation_stats: dict[str, int] | None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    with (output_dir / "columns.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["position", "column"])
        writer.writerows(enumerate(header, start=1))

    with (output_dir / "column_profile.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "position",
                "column",
                "populated_count",
                "empty_string_count",
                "dash_count",
                "na_literal_count",
            ]
        )
        for position, (name, profile) in enumerate(zip(header, profiles), start=1):
            writer.writerow(
                [
                    position,
                    name,
                    profile[0],
                    profile[1],
                    profile[2],
                    profile[3],
                ]
            )

    with (output_dir / "categorical_distributions.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["column", "value", "count", "percentage_of_rows"])
        for name, counter in distributions.items():
            for value, count in sorted(counter.items(), key=lambda item: (-item[1], item[0])):
                percentage = (100 * count / row_count) if row_count else 0
                writer.writerow([name, value, count, f"{percentage:.6f}"])

    with (output_dir / "sample_rows.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(sample_rows)

    stat = input_path.stat()
    summary = {
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_file": str(input_path.resolve()),
        "input_file_size_bytes": stat.st_size,
        "input_file_last_modified_utc": datetime.fromtimestamp(
            stat.st_mtime, timezone.utc
        ).isoformat(),
        "scan_complete": scan_complete,
        "rows_scanned": row_count,
        "malformed_rows": malformed_rows,
        "column_count": len(header),
        "chunk_size": chunk_size,
        "sample_row_count": len(sample_rows),
        "distribution_columns": list(distributions),
        "variation_id_statistics": variation_stats,
        "missing_value_note": (
            "Boş metin, '-' ve 'na' ayrı sayılmıştır; biyolojik anlamları "
            "doğrulanmadan tek bir eksik-değer sınıfında birleştirilmemiştir."
        ),
    }
    with (output_dir / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def audit_file(
    input_path: Path,
    output_dir: Path,
    *,
    chunk_size: int = 5_000,
    sample_size: int = 10,
    max_rows: int | None = None,
    progress_every: int = 250_000,
    distribution_columns: Sequence[str] = DEFAULT_DISTRIBUTION_COLUMNS,
    track_variation_ids: bool = True,
) -> dict[str, object]:
    if chunk_size <= 0:
        raise ValueError("chunk_size sıfırdan büyük olmalıdır.")
    if sample_size < 0:
        raise ValueError("sample_size negatif olamaz.")
    if max_rows is not None and max_rows <= 0:
        raise ValueError("max_rows sıfırdan büyük olmalıdır.")
    if not input_path.is_file():
        raise FileNotFoundError(f"Girdi dosyası bulunamadı: {input_path}")

    output_dir.mkdir(parents=True, exist_ok=True)
    temporary_directory: tempfile.TemporaryDirectory[str] | None = None
    connection: sqlite3.Connection | None = None

    try:
        if track_variation_ids:
            temporary_directory = tempfile.TemporaryDirectory(
                prefix="variation_index_", dir=output_dir
            )
            index_path = Path(temporary_directory.name) / "variation_ids.sqlite3"
            connection = sqlite3.connect(index_path)
            create_variation_index(connection)

        with gzip.open(input_path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValueError("Girdi dosyası boş.") from exc

            if len(header) != len(set(header)):
                raise ValueError("Sütun adları benzersiz değil.")

            column_count = len(header)
            indices = {name: header.index(name) for name in distribution_columns if name in header}
            missing_distribution_columns = [
                name for name in distribution_columns if name not in header
            ]
            if missing_distribution_columns:
                print(
                    "Uyarı: bulunamayan dağılım sütunları: "
                    + ", ".join(missing_distribution_columns)
                )

            variation_id_index = header.index("VariationID") if "VariationID" in header else None
            assembly_index = header.index("Assembly") if "Assembly" in header else None
            # Her sütun için sırasıyla: dolu, boş metin, '-', 'na'. Liste
            # kullanımı milyonlarca hücrede Counter'a göre belirgin biçimde hızlıdır.
            profiles = [[0, 0, 0, 0] for _ in header]
            distributions = {name: Counter() for name in indices}
            samples: list[list[str]] = []
            row_count = 0
            malformed_rows = 0
            next_progress = progress_every
            limit_reached = False

            for chunk in batched(reader, chunk_size):
                chunk_variation_occurrences: Counter[tuple[str, str]] = Counter()

                for raw_values in chunk:
                    if max_rows is not None and row_count >= max_rows:
                        limit_reached = True
                        break

                    values, malformed = normalize_row(raw_values, column_count)
                    malformed_rows += int(malformed)
                    row_count += 1

                    if len(samples) < sample_size:
                        samples.append(values.copy())

                    for index, value in enumerate(values):
                        if value == "":
                            profiles[index][1] += 1
                        elif value == "-":
                            profiles[index][2] += 1
                        elif len(value) == 2 and value.lower() == "na":
                            profiles[index][3] += 1
                        else:
                            profiles[index][0] += 1

                    for name, index in indices.items():
                        value = values[index] if values[index] != "" else "<EMPTY>"
                        distributions[name][value] += 1

                    if connection is not None and variation_id_index is not None:
                        variation_id = values[variation_id_index]
                        if variation_id and variation_id != "-":
                            assembly = (
                                values[assembly_index]
                                if assembly_index is not None
                                else "<NO_ASSEMBLY_COLUMN>"
                            )
                            chunk_variation_occurrences[
                                (variation_id, assembly or "<EMPTY>")
                            ] += 1

                if connection is not None and chunk_variation_occurrences:
                    update_variation_index(
                        connection,
                        chunk_variation_occurrences,
                    )

                if progress_every > 0 and row_count >= next_progress:
                    print(f"İlerleme: {row_count:,} satır tarandı.", flush=True)
                    while next_progress <= row_count:
                        next_progress += progress_every

                if limit_reached:
                    break

        scan_complete = not limit_reached
        variation_stats = variation_statistics(connection) if connection is not None else None
        if connection is not None:
            write_repeat_examples(
                connection,
                output_dir / "variation_id_repeat_examples.csv",
            )

        write_reports(
            output_dir=output_dir,
            input_path=input_path,
            header=header,
            sample_rows=samples,
            profiles=profiles,
            distributions=distributions,
            row_count=row_count,
            malformed_rows=malformed_rows,
            scan_complete=scan_complete,
            chunk_size=chunk_size,
            variation_stats=variation_stats,
        )

        print(f"\nSütunlar ({len(header)}):")
        for position, name in enumerate(header, start=1):
            print(f"{position:>2}. {name}")
        print_preview(header, samples, DEFAULT_PREVIEW_COLUMNS)
        print(f"\nTaranan satır: {row_count:,}")
        print(f"Bozuk sütun sayılı satır: {malformed_rows:,}")
        print(f"Tarama tamamlandı mı: {'evet' if scan_complete else 'hayır (sınırlı deneme)'}")
        print(f"Rapor klasörü: {output_dir.resolve()}")

        return {
            "row_count": row_count,
            "column_count": len(header),
            "malformed_rows": malformed_rows,
            "scan_complete": scan_complete,
            "variation_id_statistics": variation_stats,
        }
    finally:
        if connection is not None:
            connection.close()
        if temporary_directory is not None:
            temporary_directory.cleanup()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="ClinVar variant_summary.txt.gz dosyasını parçalı olarak denetler."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--chunk-size", type=int, default=5_000)
    parser.add_argument("--sample-rows", type=int, default=10)
    parser.add_argument(
        "--max-rows",
        type=int,
        default=None,
        help="Yalnızca hızlı deneme için en fazla okunacak veri satırı.",
    )
    parser.add_argument("--progress-every", type=int, default=250_000)
    parser.add_argument(
        "--distribution-column",
        action="append",
        dest="distribution_columns",
        help="Dağılımı sayılacak sütun. Birden çok kez verilebilir.",
    )
    parser.add_argument(
        "--skip-variation-index",
        action="store_true",
        help="VariationID tekrar indeksini atlayarak daha hızlı tarama yap.",
    )
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    audit_file(
        args.input,
        args.output_dir,
        chunk_size=args.chunk_size,
        sample_size=args.sample_rows,
        max_rows=args.max_rows,
        progress_every=args.progress_every,
        distribution_columns=(
            args.distribution_columns
            if args.distribution_columns is not None
            else DEFAULT_DISTRIBUTION_COLUMNS
        ),
        track_variation_ids=not args.skip_variation_index,
    )


if __name__ == "__main__":
    main()

