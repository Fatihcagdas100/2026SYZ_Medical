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


if __name__ == "__main__":
    unittest.main()
