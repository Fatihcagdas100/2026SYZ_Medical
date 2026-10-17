import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.cross_validate_gene_disjoint import cross_validate


class CrossValidateGeneDisjointTest(unittest.TestCase):
    def test_cross_validation_uses_only_train_and_keeps_groups_together(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            features_path = root / "features.tsv.gz"
            split_path = root / "split.tsv"
            report_dir = root / "report"

            feature_rows = []
            split_rows = []
            variation_id = 0
            for partition, group_count in (("train", 30), ("validation", 3), ("test", 3)):
                for group_index in range(group_count):
                    gene_group = f"{partition}_gene_{group_index}"
                    for target in (0, 1):
                        variation_id += 1
                        pathogenic = target == 1
                        feature_rows.append(
                            {
                                "VariationID": str(variation_id),
                                "variant_type": "Deletion" if pathogenic else "SNV",
                                "chromosome": "1",
                                "reference_length": 2 if pathogenic else 1,
                                "alternate_length": 1,
                                "length_delta": -1 if pathogenic else 0,
                                "absolute_length_delta": 1 if pathogenic else 0,
                                "vcf_alleles_missing": 0,
                                "is_canonical_snv": 0 if pathogenic else 1,
                                "substitution": "NOT_CANONICAL_SNV" if pathogenic else "A>G",
                                "mutation_class": (
                                    "not_canonical_snv" if pathogenic else "transition"
                                ),
                                "hgvs_sequence_type": "coding_transcript",
                                "dna_change_type": "substitution",
                                "coding_region": "coding_exonic",
                                "splice_offset_category": "not_intronic",
                                "protein_effect": (
                                    "stop_gained_or_changed" if pathogenic else "synonymous"
                                ),
                                "has_protein_annotation": 1,
                                "target": target,
                            }
                        )
                        split_rows.append(
                            {
                                "VariationID": str(variation_id),
                                "target": target,
                                "split": partition,
                                "gene_group": gene_group,
                                "fold": 0,
                            }
                        )

            pd.DataFrame(feature_rows).to_csv(
                features_path, sep="\t", index=False, compression="gzip"
            )
            pd.DataFrame(split_rows).to_csv(split_path, sep="\t", index=False)

            manifest = cross_validate(
                features_path,
                split_path,
                report_dir,
                n_splits=5,
                chunk_size=11,
            )

            self.assertFalse(manifest["data_policy"]["fixed_validation_used"])
            self.assertFalse(manifest["data_policy"]["test_used"])
            self.assertEqual(manifest["counts"]["train_rows"], 60)
            self.assertEqual(manifest["counts"]["out_of_fold_prediction_rows"], 60)
            fold_metrics = pd.read_csv(report_dir / "fold_metrics.csv")
            self.assertEqual(len(fold_metrics), 5)
            predictions = pd.read_csv(
                report_dir / "train_out_of_fold_predictions.tsv.gz", sep="\t"
            )
            train_ids = {
                int(row["VariationID"])
                for row in split_rows
                if row["split"] == "train"
            }
            excluded_ids = {
                int(row["VariationID"])
                for row in split_rows
                if row["split"] != "train"
            }
            self.assertEqual(set(predictions["VariationID"]), train_ids)
            self.assertFalse(set(predictions["VariationID"]) & excluded_ids)
            self.assertEqual(predictions.groupby("gene_group")["fold"].nunique().max(), 1)


if __name__ == "__main__":
    unittest.main()
