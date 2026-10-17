import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.audit_split_transition import audit_split_transition


class AuditSplitTransitionTest(unittest.TestCase):
    def test_audit_reports_previous_roles_of_new_test(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            old_path = root / "old.tsv"
            new_path = root / "new.tsv"
            report_dir = root / "report"
            ids = [str(index) for index in range(1, 13)]
            targets = [0, 1] * 6
            pd.DataFrame(
                {
                    "VariationID": ids,
                    "target": targets,
                    "split": ["train"] * 6 + ["validation"] * 3 + ["test"] * 3,
                }
            ).to_csv(old_path, sep="\t", index=False)
            pd.DataFrame(
                {
                    "VariationID": ids,
                    "target": targets,
                    "split": [
                        "train",
                        "validation",
                        "train",
                        "validation",
                        "test",
                        "test",
                        "train",
                        "validation",
                        "test",
                        "train",
                        "validation",
                        "test",
                    ],
                }
            ).to_csv(new_path, sep="\t", index=False)

            manifest = audit_split_transition(old_path, new_path, report_dir)

            exposure = manifest["final_test_prior_exposure"]
            self.assertEqual(exposure["final_test_rows"], 4)
            self.assertEqual(exposure["previous_train_rows"], 2)
            self.assertEqual(exposure["previous_validation_rows"], 1)
            self.assertEqual(exposure["previous_test_rows"], 1)
            self.assertAlmostEqual(exposure["previous_development_percentage"], 75.0)
            self.assertFalse(
                manifest["interpretation"]["fully_process_unseen_final_test"]
            )
            self.assertTrue((report_dir / "split_transition_heatmap.png").is_file())
            self.assertTrue((report_dir / "split_transition_report.md").is_file())

    def test_audit_rejects_target_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            old_path = root / "old.tsv"
            new_path = root / "new.tsv"
            base = pd.DataFrame(
                {
                    "VariationID": ["1", "2", "3", "4", "5", "6"],
                    "target": [0, 1, 0, 1, 0, 1],
                    "split": ["train", "train", "validation", "validation", "test", "test"],
                }
            )
            base.to_csv(old_path, sep="\t", index=False)
            changed = base.copy()
            changed.loc[0, "target"] = 1
            changed.to_csv(new_path, sep="\t", index=False)

            with self.assertRaisesRegex(ValueError, "hedef uyuşmazlığı"):
                audit_split_transition(old_path, new_path, root / "report")


if __name__ == "__main__":
    unittest.main()
