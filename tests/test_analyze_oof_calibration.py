import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.analyze_oof_calibration import analyze_calibration


class AnalyzeOofCalibrationTest(unittest.TestCase):
    def make_predictions(self) -> pd.DataFrame:
        rows = []
        variation_id = 0
        for fold in range(1, 5):
            for group_index in range(5):
                for target, probability in ((0, 0.05 + fold * 0.01), (1, 0.94 - fold * 0.01)):
                    variation_id += 1
                    rows.append(
                        {
                            "VariationID": str(variation_id),
                            "gene_group": f"gene_{fold}_{group_index}",
                            "source_partition": "train",
                            "fold": fold,
                            "target": target,
                            "probability_pathogenic": probability,
                            "prediction": target,
                        }
                    )
        return pd.DataFrame(rows)

    def test_analysis_uses_train_oof_and_writes_calibration_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            predictions_path = root / "oof.tsv.gz"
            report_dir = root / "calibration"
            self.make_predictions().to_csv(
                predictions_path, sep="\t", index=False, compression="gzip"
            )

            manifest = analyze_calibration(
                predictions_path, report_dir, n_bins=4
            )

            self.assertFalse(manifest["data_policy"]["fixed_validation_used"])
            self.assertFalse(manifest["data_policy"]["test_used"])
            self.assertFalse(manifest["data_policy"]["recalibration_applied"])
            self.assertEqual(manifest["integrity"]["oof_rows"], 40)
            self.assertEqual(
                manifest["integrity"]["gene_groups_in_multiple_folds"], 0
            )
            self.assertLess(manifest["overall_metrics"]["brier_score"], 0.02)
            self.assertTrue((report_dir / "quantile_reliability.csv").is_file())
            self.assertTrue((report_dir / "fixed_width_reliability.csv").is_file())
            self.assertTrue((report_dir / "fold_calibration_metrics.csv").is_file())
            self.assertTrue((report_dir / "calibration_reliability.png").is_file())

    def test_analysis_rejects_non_train_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            predictions_path = root / "oof.tsv.gz"
            predictions = self.make_predictions()
            predictions.loc[0, "source_partition"] = "validation"
            predictions.to_csv(
                predictions_path, sep="\t", index=False, compression="gzip"
            )

            with self.assertRaisesRegex(ValueError, "yalnızca train OOF"):
                analyze_calibration(predictions_path, root / "report", n_bins=4)


if __name__ == "__main__":
    unittest.main()
