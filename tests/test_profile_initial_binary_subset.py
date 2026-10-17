import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.profile_initial_binary_subset import profile_subset


class ProfileInitialBinarySubsetTest(unittest.TestCase):
    def test_filter_funnel_duplicates_and_conflicts(self) -> None:
        header = [
            "VariationID",
            "Assembly",
            "OriginSimple",
            "ClinicalSignificance",
            "ReviewStatus",
            "Type",
            "PositionVCF",
            "ReferenceAlleleVCF",
            "AlternateAlleleVCF",
        ]
        rows = [
            ["1", "GRCh37", "germline", "Pathogenic", "reviewed", "SNV", "1", "A", "T"],
            ["2", "GRCh38", "somatic", "Pathogenic", "reviewed", "SNV", "2", "A", "T"],
            ["3", "GRCh38", "germline", "Pathogenic", "reviewed", "SNV", "3", "A", "T"],
            ["3", "GRCh38", "germline", "Pathogenic", "reviewed", "SNV", "3", "A", "T"],
            ["4", "GRCh38", "germline", "Benign", "reviewed", "SNV", "4", "G", "C"],
            ["5", "GRCh38", "germline", "Pathogenic", "reviewed", "SNV", "5", "G", "C"],
            ["5", "GRCh38", "germline", "Benign", "reviewed", "SNV", "5", "G", "C"],
            ["6", "GRCh38", "germline", "Likely pathogenic", "reviewed", "SNV", "6", "G", "C"],
        ]

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "input.txt.gz"
            output_dir = root / "report"
            with gzip.open(input_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                writer.writerows(rows)

            result = profile_subset(
                input_path,
                output_dir,
                chunk_size=2,
                progress_every=0,
            )

            self.assertEqual(result["counts"]["total_rows"], 8)
            self.assertEqual(result["counts"]["grch38_rows"], 7)
            self.assertEqual(result["counts"]["grch38_germline_rows"], 6)
            self.assertEqual(result["counts"]["exact_binary_candidate_rows"], 5)
            stats = result["candidate_statistics"]
            self.assertEqual(stats["unique_candidate_variation_ids"], 3)
            self.assertEqual(stats["candidate_ids_with_multiple_rows"], 2)
            self.assertEqual(stats["candidate_ids_with_conflicting_exact_labels"], 1)
            self.assertEqual(stats["usable_nonconflicting_unique_ids"], 2)
            self.assertEqual(stats["usable_unique_ids_by_label"], {"Benign": 1, "Pathogenic": 1})


if __name__ == "__main__":
    unittest.main()

