import csv
import gzip
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.train_logistic_baseline import train_baselines


class TrainLogisticBaselineTest(unittest.TestCase):
    def test_training_uses_validation_and_keeps_test_locked(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            features_path = root / "features.tsv.gz"
            split_path = root / "split.tsv"
            report_dir = root / "reports"
            artifact_dir = root / "artifacts"

            fieldnames = [
                "VariationID",
                "variant_type",
                "chromosome",
                "reference_length",
                "alternate_length",
                "length_delta",
                "absolute_length_delta",
                "vcf_alleles_missing",
                "is_canonical_snv",
                "substitution",
                "mutation_class",
                "hgvs_sequence_type",
                "dna_change_type",
                "coding_region",
                "splice_offset_category",
                "protein_effect",
                "has_protein_annotation",
                "missense_reference_residue",
                "missense_alternate_residue",
                "missense_substitution",
                "missense_residue_class_transition",
                "missense_charge_transition",
                "missense_features_available",
                "missense_grantham_distance",
                "missense_blosum62_score",
                "missense_hydropathy_delta",
                "missense_absolute_hydropathy_delta",
                "target",
            ]
            rows = []
            split_rows = []
            allocation = {
                "train": (4, 4),
                "validation": (2, 2),
                "test": (1, 1),
            }
            variation_id = 0
            for split_name, class_counts in allocation.items():
                for target, count in enumerate(class_counts):
                    for _ in range(count):
                        variation_id += 1
                        is_pathogenic = target == 1
                        rows.append(
                            {
                                "VariationID": str(variation_id),
                                "variant_type": (
                                    "Deletion" if is_pathogenic else "single nucleotide variant"
                                ),
                                "chromosome": "1" if is_pathogenic else "2",
                                "reference_length": 2 if is_pathogenic else 1,
                                "alternate_length": 1,
                                "length_delta": -1 if is_pathogenic else 0,
                                "absolute_length_delta": 1 if is_pathogenic else 0,
                                "vcf_alleles_missing": 0,
                                "is_canonical_snv": 0 if is_pathogenic else 1,
                                "substitution": (
                                    "NOT_CANONICAL_SNV" if is_pathogenic else "A>G"
                                ),
                                "mutation_class": (
                                    "not_canonical_snv" if is_pathogenic else "transition"
                                ),
                                "hgvs_sequence_type": "coding_transcript",
                                "dna_change_type": "substitution",
                                "coding_region": "coding_exonic",
                                "splice_offset_category": "not_intronic",
                                "protein_effect": (
                                    "stop_gained_or_changed"
                                    if is_pathogenic
                                    else "synonymous"
                                ),
                                "has_protein_annotation": 1,
                                "missense_reference_residue": "R",
                                "missense_alternate_residue": "G",
                                "missense_substitution": "R>G",
                                "missense_residue_class_transition": "positive>special",
                                "missense_charge_transition": "positive>neutral",
                                "missense_features_available": 1,
                                "missense_grantham_distance": 125,
                                "missense_blosum62_score": -2,
                                "missense_hydropathy_delta": 4.1,
                                "missense_absolute_hydropathy_delta": 4.1,
                                "target": target,
                            }
                        )
                        split_rows.append(
                            {
                                "VariationID": str(variation_id),
                                "target": target,
                                "split": split_name,
                            }
                        )

            with gzip.open(features_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerows(rows)
            with split_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["VariationID", "target", "split"],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(split_rows)

            manifest = train_baselines(
                features_path,
                split_path,
                report_dir,
                artifact_dir,
                seed=2027,
                chunk_size=5,
            )

            self.assertEqual(
                manifest["partition_counts"],
                {"train": 8, "validation": 4, "test": 2},
            )
            self.assertEqual(
                manifest["data_policy"]["test_partition"],
                "locked_not_loaded_into_development_dataframe",
            )
            self.assertEqual(
                set(manifest["validation_metrics"]),
                {"without_missingness_flag", "with_missingness_flag"},
            )

            predictions = pd.read_csv(
                report_dir / "validation_predictions.tsv.gz", sep="\t"
            )
            validation_ids = {
                int(row["VariationID"])
                for row in split_rows
                if row["split"] == "validation"
            }
            test_ids = {
                int(row["VariationID"])
                for row in split_rows
                if row["split"] == "test"
            }
            self.assertEqual(set(predictions["split"]), {"validation"})
            self.assertEqual(set(predictions["VariationID"]), validation_ids)
            self.assertFalse(set(predictions["VariationID"]) & test_ids)
            self.assertEqual(len(predictions), 8)
            self.assertEqual(len(list(artifact_dir.glob("*.joblib"))), 2)

            hgvs_manifest = train_baselines(
                features_path,
                split_path,
                report_dir / "hgvs",
                artifact_dir / "hgvs",
                seed=2027,
                chunk_size=5,
                feature_profile="hgvs",
            )
            self.assertEqual(hgvs_manifest["data_policy"]["feature_profile"], "hgvs")
            self.assertIn(
                "protein_effect",
                hgvs_manifest["models"]["without_missingness_flag"]["feature_columns"],
            )

            biological_manifest = train_baselines(
                features_path,
                split_path,
                report_dir / "hgvs_biological",
                artifact_dir / "hgvs_biological",
                seed=2027,
                chunk_size=5,
                feature_profile="hgvs_biological",
            )
            biological_columns = biological_manifest["models"][
                "without_missingness_flag"
            ]["feature_columns"]
            self.assertIn("protein_effect", biological_columns)
            self.assertNotIn("hgvs_sequence_type", biological_columns)
            self.assertNotIn("has_protein_annotation", biological_columns)

            missense_manifest = train_baselines(
                features_path,
                split_path,
                report_dir / "hgvs_missense",
                artifact_dir / "hgvs_missense",
                seed=2027,
                chunk_size=5,
                feature_profile="hgvs_missense",
            )
            missense_columns = missense_manifest["models"][
                "without_missingness_flag"
            ]["feature_columns"]
            self.assertIn("missense_grantham_distance", missense_columns)
            self.assertIn("missense_substitution", missense_columns)


if __name__ == "__main__":
    unittest.main()
