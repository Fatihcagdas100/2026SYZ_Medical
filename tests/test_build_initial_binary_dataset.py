import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.build_initial_binary_dataset import OUTPUT_COLUMNS, build_dataset


class BuildInitialBinaryDatasetTest(unittest.TestCase):
    def test_deduplication_target_mapping_and_conflict_exclusion(self) -> None:
        header = [
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
            "ClinicalSignificance",
        ]
        rows = [
            ["1", "10", "SNV", "first", "1", "A", "HGNC:1", "GRCh38", "NC_1", "1", "10", "10", "10", "A", "T", "germline", "Pathogenic"],
            ["1", "10", "SNV", "different duplicate", "1", "A", "HGNC:1", "GRCh38", "NC_1", "1", "10", "10", "10", "A", "T", "germline", "Pathogenic"],
            ["2", "20", "SNV", "benign", "2", "B", "HGNC:2", "GRCh38", "NC_2", "2", "20", "20", "20", "G", "C", "germline", "Benign"],
            ["6", "60", "SNV", "identical", "6", "F", "HGNC:6", "GRCh38", "NC_6", "6", "60", "60", "60", "A", "G", "germline", "Pathogenic"],
            ["6", "60", "SNV", "identical", "6", "F", "HGNC:6", "GRCh38", "NC_6", "6", "60", "60", "60", "A", "G", "germline", "Pathogenic"],
            ["3", "30", "SNV", "conflict", "3", "C", "HGNC:3", "GRCh38", "NC_3", "3", "30", "30", "30", "C", "A", "germline", "Pathogenic"],
            ["3", "30", "SNV", "conflict", "3", "C", "HGNC:3", "GRCh38", "NC_3", "3", "30", "30", "30", "C", "A", "germline", "Benign"],
            ["4", "40", "SNV", "old assembly", "4", "D", "HGNC:4", "GRCh37", "NC_4", "4", "40", "40", "40", "T", "G", "germline", "Pathogenic"],
        ]

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "input.tsv.gz"
            output_path = root / "processed" / "output.tsv.gz"
            report_dir = root / "reports"
            with gzip.open(input_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                writer.writerows(rows)

            manifest = build_dataset(
                input_path,
                output_path,
                report_dir,
                chunk_size=2,
                progress_every=0,
            )

            with gzip.open(output_path, "rt", encoding="utf-8", newline="") as handle:
                output_rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(list(output_rows[0]), list(OUTPUT_COLUMNS))
            self.assertEqual(len(output_rows), 2)
            self.assertEqual(output_rows[0]["VariationID"], "2")
            self.assertEqual(output_rows[0]["target"], "0")
            self.assertEqual(output_rows[1]["VariationID"], "6")
            self.assertEqual(output_rows[1]["target"], "1")
            self.assertEqual(manifest["candidate_index"]["ids_with_conflicting_exact_labels"], 1)
            self.assertEqual(
                manifest["candidate_index"]["ids_with_multiple_distinct_output_rows"],
                1,
            )
            self.assertEqual(manifest["second_pass"]["duplicate_rows_skipped"], 1)
            self.assertEqual(
                manifest["second_pass"]["duplicate_ids_with_output_column_differences"],
                1,
            )


if __name__ == "__main__":
    unittest.main()

