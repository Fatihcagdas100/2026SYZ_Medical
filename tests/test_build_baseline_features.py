import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.build_baseline_features import build_features


class BuildBaselineFeaturesTest(unittest.TestCase):
    def test_lengths_substitutions_and_missing_alleles(self) -> None:
        header = [
            "VariationID",
            "Type",
            "Chromosome",
            "ReferenceAlleleVCF",
            "AlternateAlleleVCF",
            "target_label",
            "target",
        ]
        rows = [
            ["1", "single nucleotide variant", "1", "A", "G", "Pathogenic", "1"],
            ["2", "Deletion", "2", "AT", "A", "Benign", "0"],
            ["3", "Deletion", "3", "na", "na", "Pathogenic", "1"],
        ]

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "input.tsv.gz"
            output_path = root / "output.tsv.gz"
            report_dir = root / "report"
            with gzip.open(input_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                writer.writerows(rows)

            manifest = build_features(
                input_path,
                output_path,
                report_dir,
                progress_every=0,
            )

            with gzip.open(output_path, "rt", encoding="utf-8", newline="") as handle:
                output_rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(output_rows[0]["reference_length"], "1")
            self.assertEqual(output_rows[0]["substitution"], "A>G")
            self.assertEqual(output_rows[0]["mutation_class"], "transition")
            self.assertEqual(output_rows[1]["length_delta"], "-1")
            self.assertEqual(output_rows[1]["absolute_length_delta"], "1")
            self.assertEqual(output_rows[1]["is_canonical_snv"], "0")
            self.assertEqual(output_rows[2]["reference_length"], "")
            self.assertEqual(output_rows[2]["vcf_alleles_missing"], "1")
            self.assertEqual(manifest["counts"]["rows_with_missing_vcf_alleles"], 1)


if __name__ == "__main__":
    unittest.main()

