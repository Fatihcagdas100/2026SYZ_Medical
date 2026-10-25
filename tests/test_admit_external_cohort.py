import csv
import gzip
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.admit_external_cohort import admit_external_cohort


class AdmitExternalCohortTest(unittest.TestCase):
    def test_removes_identifier_and_coordinate_overlap_and_keeps_binary_cohort(self) -> None:
        columns = ["VariationID", "#AlleleID", "Type", "Name", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF", "target_label", "target"]
        development_rows = [
            ["1", "11", "SNV", "NM_1(G):c.1A>G (p.Arg1Gly)", "1", "100", "A", "G", "Benign", "0"],
        ]
        external_rows = [
            ["1", "21", "SNV", "NM_1(G):c.1A>G (p.Arg1Gly)", "1", "100", "A", "G", "Benign", "0"],
            ["2", "22", "SNV", "NM_2(G):c.2C>T (p.Arg1Cys)", "1", "100", "A", "G", "Pathogenic", "1"],
            ["3", "23", "SNV", "NM_3(G):c.3G>A (p.Gly1Arg)", "2", "200", "G", "A", "Benign", "0"],
            ["4", "24", "SNV", "NM_4(G):c.4T>C (p.Val2Ala)", "3", "300", "T", "C", "Pathogenic", "1"],
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            development = root / "development.tsv.gz"
            external = root / "external.tsv.gz"
            for path, rows in ((development, development_rows), (external, external_rows)):
                with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t")
                    writer.writerow(columns)
                    writer.writerows(rows)
            output = root / "admitted.tsv.gz"
            manifest = admit_external_cohort(development, external, output, root / "report", cohort_kind="independent", source_version="v1", source_url="https://example.org")
            admitted = pd.read_csv(output, sep="\t")
            self.assertEqual(set(admitted["VariationID"].astype(str)), {"3", "4"})
            self.assertEqual(manifest["counts"]["excluded_variation_id"], 1)
            self.assertEqual(manifest["counts"]["excluded_coordinate_allele"], 2)
            self.assertEqual(manifest["admission_decision"], "accepted")


if __name__ == "__main__":
    unittest.main()
