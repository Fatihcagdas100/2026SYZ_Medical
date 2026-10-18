import csv
import gzip
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.build_external_features import build_external_features


class BuildExternalFeaturesTest(unittest.TestCase):
    def test_output_matches_model_schema_without_identity_text(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "external.tsv.gz"
            rows = [
                {"VariationID": "1", "Type": "SNV", "Name": "NM_1(G):c.1A>G (p.Arg1Gly)", "Chromosome": "1", "ReferenceAlleleVCF": "A", "AlternateAlleleVCF": "G", "target_label": "Benign", "target": "0"},
                {"VariationID": "2", "Type": "Deletion", "Name": "NM_2(G):c.2del (p.Gly1fs)", "Chromosome": "2", "ReferenceAlleleVCF": "AT", "AlternateAlleleVCF": "A", "target_label": "Pathogenic", "target": "1"},
            ]
            with gzip.open(source, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            output = root / "features.tsv.gz"
            manifest = build_external_features(source, output, root / "report", progress_every=0)
            features = pd.read_csv(output, sep="\t")
            self.assertEqual(len(features), 2)
            self.assertNotIn("Name", features.columns)
            self.assertEqual(features.loc[0, "missense_grantham_distance"], 125)
            self.assertEqual(manifest["counts"]["missense_rows"], 1)


if __name__ == "__main__":
    unittest.main()
