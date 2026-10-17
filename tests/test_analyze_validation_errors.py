import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.analyze_validation_errors import analyze_errors


class AnalyzeValidationErrorsTest(unittest.TestCase):
    def test_analysis_joins_only_validation_ids_and_reports_groups(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            predictions_path = root / "predictions.tsv.gz"
            features_path = root / "features.tsv.gz"
            report_dir = root / "reports"

            predictions = pd.DataFrame(
                {
                    "VariationID": [str(index) for index in range(1, 9)],
                    "model_name": ["without_missingness_flag"] * 8,
                    "split": ["validation"] * 8,
                    "target": [0, 0, 0, 0, 1, 1, 1, 1],
                    "probability_pathogenic": [
                        0.05,
                        0.20,
                        0.60,
                        0.90,
                        0.10,
                        0.40,
                        0.70,
                        0.95,
                    ],
                }
            )
            predictions.to_csv(
                predictions_path,
                sep="\t",
                index=False,
                compression="gzip",
            )

            feature_rows = []
            for index in range(1, 10):
                canonical = index % 2
                feature_rows.append(
                    {
                        "VariationID": str(index),
                        "variant_type": "SNV" if canonical else "Deletion",
                        "chromosome": "1" if index <= 5 else "2",
                        "mutation_class": "transition" if canonical else "not_canonical_snv",
                        "substitution": "A>G" if canonical else "NOT_CANONICAL_SNV",
                        "is_canonical_snv": canonical,
                        "vcf_alleles_missing": 0,
                        "reference_length": 1 if canonical else 2,
                        "alternate_length": 1,
                        "length_delta": 0 if canonical else -1,
                        "absolute_length_delta": 0 if canonical else 1,
                        "target": int(index >= 5),
                    }
                )
            pd.DataFrame(feature_rows).to_csv(
                features_path,
                sep="\t",
                index=False,
                compression="gzip",
            )

            manifest = analyze_errors(
                predictions_path,
                features_path,
                report_dir,
                threshold=0.5,
                min_group_size=2,
                top_errors=2,
                chunk_size=3,
            )

            self.assertFalse(manifest["data_policy"]["test_data_used"])
            self.assertEqual(manifest["counts"]["validation_rows"], 8)
            self.assertEqual(manifest["counts"]["false_negatives"], 2)
            self.assertEqual(manifest["counts"]["false_positives"], 2)

            metrics = pd.read_csv(report_dir / "subgroup_metrics.csv")
            self.assertIn("variant_type", set(metrics["group_feature"]))
            self.assertIn("absolute_length_delta_bucket", set(metrics["group_feature"]))
            errors = pd.read_csv(report_dir / "most_confident_errors.tsv", sep="\t")
            self.assertEqual(set(errors["error_type"]), {"false_negative", "false_positive"})
            self.assertNotIn(9, set(errors["VariationID"]))
            self.assertTrue((report_dir / "error_rates_by_variant_type.png").is_file())


if __name__ == "__main__":
    unittest.main()
