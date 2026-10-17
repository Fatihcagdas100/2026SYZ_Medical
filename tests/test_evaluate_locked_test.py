import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import joblib
import pandas as pd

from scripts.audit_test_dtype_correction import audit_dtype_correction
from scripts.evaluate_locked_test import evaluate_locked_test
from scripts.train_logistic_baseline import build_pipeline


class EvaluateLockedTestTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.features_path = self.root / "features.tsv.gz"
        self.split_path = self.root / "split.tsv"
        self.artifact_path = self.root / "model.joblib"
        self.lock_path = self.root / "protocol_lock.json"
        self.output_dir = self.root / "final_test"

        features = pd.DataFrame(
            {
                "VariationID": [str(index) for index in range(1, 13)],
                "category": ["low", "high"] * 6,
                "feature_a": [-3, 3, -2, 2, -1, 1, -2.5, -1.5, 1.5, 2.5, -1, 1],
                "target": [0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1],
            }
        )
        features.to_csv(
            self.features_path, sep="\t", index=False, compression="gzip"
        )
        splits = pd.DataFrame(
            {
                "VariationID": [str(index) for index in range(1, 13)],
                "split": ["train"] * 4 + ["validation"] * 2 + ["test"] * 6,
                "gene_group": [f"gene_{index}" for index in range(1, 13)],
            }
        )
        splits.to_csv(self.split_path, sep="\t", index=False)
        model = build_pipeline(
            ["category", "feature_a"], seed=2027, categorical_features=["category"]
        ).fit(
            features.loc[:5, ["category", "feature_a"]],
            features.loc[:5, "target"],
        )
        joblib.dump(model, self.artifact_path)
        lock = {
            "status": "pre_test_frozen",
            "test_evaluated": False,
            "test_labels_read": False,
            "test_policy": {"maximum_evaluations": 1},
            "frozen_operating_point": {
                "model_name": "synthetic_model",
                "threshold": 0.5,
            },
            "validation_reference": {"roc_auc": 1.0},
            "test_acceptance": {
                "minimum_roc_auc": 0.9,
                "minimum_average_precision": 0.9,
                "minimum_sensitivity_at_frozen_threshold": 0.8,
                "minimum_specificity_at_frozen_threshold": 0.8,
                "maximum_roc_auc_drop_from_validation": 0.1,
            },
            "artifact": {
                "path": str(self.artifact_path),
                "sha256": self.sha(self.artifact_path),
                "feature_columns": ["category", "feature_a"],
            },
            "data_fingerprints": {
                "features_sha256": self.sha(self.features_path),
                "gene_disjoint_split_sha256": self.sha(self.split_path),
            },
        }
        self.lock_path.write_text(json.dumps(lock), encoding="utf-8")

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    @staticmethod
    def sha(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def test_locked_test_runs_once_and_records_primary_result(self) -> None:
        manifest = evaluate_locked_test(
            self.lock_path,
            self.features_path,
            self.split_path,
            self.output_dir,
            chunk_size=3,
        )
        self.assertTrue(manifest["test_evaluated"])
        self.assertTrue(manifest["test_labels_read"])
        self.assertEqual(manifest["test_rows"], 6)
        self.assertTrue(manifest["all_acceptance_gates_passed"])
        self.assertTrue((self.output_dir / "EVALUATION_STARTED.json").is_file())
        self.assertTrue((self.output_dir / "test_report.md").is_file())
        with self.assertRaisesRegex(FileExistsError, "tek-sefer"):
            evaluate_locked_test(
                self.lock_path,
                self.features_path,
                self.split_path,
                self.output_dir,
            )

    def test_locked_test_rejects_changed_feature_file_before_opening(self) -> None:
        changed = pd.read_csv(self.features_path, sep="\t")
        changed.loc[0, "feature_a"] = 99
        changed.to_csv(
            self.features_path, sep="\t", index=False, compression="gzip"
        )
        with self.assertRaisesRegex(ValueError, "Özellik dosyası hash"):
            evaluate_locked_test(
                self.lock_path,
                self.features_path,
                self.split_path,
                self.output_dir,
            )
        self.assertFalse(self.output_dir.exists())

    def test_dtype_audit_preserves_original_and_compares_predictions(self) -> None:
        evaluate_locked_test(
            self.lock_path,
            self.features_path,
            self.split_path,
            self.output_dir,
        )
        audit = audit_dtype_correction(
            self.output_dir / "manifest.json",
            self.output_dir / "test_predictions.tsv.gz",
            self.lock_path,
            self.features_path,
            self.split_path,
            self.root / "dtype_audit",
        )
        self.assertTrue(audit["original_evaluation_preserved"])
        self.assertFalse(audit["model_changed"])
        self.assertEqual(audit["prediction_comparison"]["classifications_changed"], 0)
        self.assertTrue((self.root / "dtype_audit" / "audit_report.md").is_file())


if __name__ == "__main__":
    unittest.main()
