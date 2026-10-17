import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.analyze_validation_thresholds import analyze_thresholds


class AnalyzeValidationThresholdsTest(unittest.TestCase):
    def make_predictions(self) -> pd.DataFrame:
        rows = []
        targets = [0, 0, 0, 0, 1, 1, 1, 1]
        probability_sets = {
            "model_a": [0.05, 0.20, 0.45, 0.70, 0.30, 0.55, 0.80, 0.95],
            "model_b": [0.10, 0.25, 0.40, 0.65, 0.35, 0.60, 0.75, 0.90],
        }
        for model_name, probabilities in probability_sets.items():
            for index, (target, probability) in enumerate(
                zip(targets, probabilities), start=1
            ):
                rows.append(
                    {
                        "VariationID": str(index),
                        "model_name": model_name,
                        "split": "validation",
                        "target": target,
                        "probability_pathogenic": probability,
                        "prediction": int(probability >= 0.5),
                    }
                )
        return pd.DataFrame(rows)

    def test_analysis_uses_only_validation_and_writes_operating_points(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            predictions_path = root / "predictions.tsv.gz"
            report_dir = root / "reports"
            self.make_predictions().to_csv(
                predictions_path,
                sep="\t",
                index=False,
                compression="gzip",
            )

            manifest = analyze_thresholds(
                predictions_path,
                report_dir,
                sensitivity_targets=(0.75, 1.0 - 1e-6),
            )

            self.assertFalse(manifest["data_policy"]["test_data_used"])
            self.assertEqual(manifest["source"]["observed_splits"], ["validation"])
            operating_points = pd.read_csv(report_dir / "operating_points.csv")
            self.assertEqual(len(operating_points), 8)
            targeted = operating_points[operating_points["sensitivity_target"].notna()]
            self.assertTrue(
                (targeted["sensitivity"] >= targeted["sensitivity_target"]).all()
            )
            for filename in (
                "threshold_sweep.csv",
                "roc_curve.csv",
                "precision_recall_curve.csv",
                "validation_roc_pr_curves.png",
                "validation_threshold_tradeoffs.png",
                "manifest.json",
            ):
                self.assertTrue((report_dir / filename).is_file(), filename)

    def test_analysis_rejects_test_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            predictions_path = root / "predictions.tsv.gz"
            predictions = self.make_predictions()
            predictions.loc[0, "split"] = "test"
            predictions.to_csv(
                predictions_path,
                sep="\t",
                index=False,
                compression="gzip",
            )

            with self.assertRaisesRegex(ValueError, "yalnızca validation"):
                analyze_thresholds(predictions_path, root / "reports")


if __name__ == "__main__":
    unittest.main()
