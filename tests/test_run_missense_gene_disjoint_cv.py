import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.run_missense_gene_disjoint_cv import build_oof_subgroup_table


class MissenseGeneDisjointCvTest(unittest.TestCase):
    def test_oof_subgroups_accept_only_train_predictions(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            features_path = root / "features.tsv.gz"
            predictions_path = root / "predictions.tsv.gz"
            pd.DataFrame(
                {
                    "VariationID": ["1", "2", "3", "4"],
                    "protein_effect": ["missense", "missense", "synonymous", "synonymous"],
                    "has_protein_annotation": [1, 1, 1, 1],
                }
            ).to_csv(features_path, sep="\t", index=False, compression="gzip")
            pd.DataFrame(
                {
                    "VariationID": ["1", "2", "3", "4"],
                    "source_partition": ["train"] * 4,
                    "target": [0, 1, 0, 1],
                    "probability_pathogenic": [0.1, 0.9, 0.2, 0.8],
                }
            ).to_csv(predictions_path, sep="\t", index=False, compression="gzip")

            table = build_oof_subgroup_table(
                predictions_path,
                features_path,
                threshold=0.5,
                chunk_size=2,
            )

            missense = table.loc[table["group"] == "missense"].iloc[0]
            self.assertEqual(missense["row_count"], 2)
            self.assertEqual(missense["roc_auc"], 1.0)


if __name__ == "__main__":
    unittest.main()
