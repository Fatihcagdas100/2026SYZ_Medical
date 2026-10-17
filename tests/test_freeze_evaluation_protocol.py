import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.freeze_evaluation_protocol import freeze_protocol


class FreezeEvaluationProtocolTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.artifact = self.root / "model.joblib"
        self.artifact.write_bytes(b"frozen-model")
        self.validation_predictions = self.root / "validation.tsv.gz"
        self.validation_predictions.write_bytes(b"validation-only")
        self.feature_hash = "feature-hash"
        self.split_hash = "split-hash"

        self.config = self.write_json(
            "config.json",
            {
                "protocol_version": "1.0",
                "status": "pre_test_frozen",
                "intended_use": "research_variant_prioritization",
                "clinical_diagnosis_authorized": False,
                "selected_model": "without_missingness_flag",
                "selected_feature_profile": "hgvs",
                "operating_point": {
                    "policy": "sensitivity_at_least_0.95",
                    "sensitivity_target": 0.95,
                },
                "test_acceptance": {
                    "minimum_roc_auc": 0.97,
                    "minimum_average_precision": 0.97,
                    "minimum_sensitivity_at_frozen_threshold": 0.93,
                    "minimum_specificity_at_frozen_threshold": 0.85,
                    "maximum_roc_auc_drop_from_validation": 0.02,
                },
                "test_policy": {
                    "maximum_evaluations": 1,
                    "model_or_threshold_changes_after_test_are_forbidden_before_primary_reporting": True,
                    "failed_acceptance_requires_a_new_external_holdout_for_future_claims": True,
                },
                "limitations": ["research only"],
            },
        )
        artifact_hash = hashlib.sha256(self.artifact.read_bytes()).hexdigest()
        validation_hash = hashlib.sha256(
            self.validation_predictions.read_bytes()
        ).hexdigest()
        self.model_manifest = self.write_json(
            "model.json",
            {
                "inputs": {
                    "features": {"sha256": self.feature_hash},
                    "split": {"sha256": self.split_hash},
                },
                "data_policy": {
                    "test_partition": "locked_not_loaded_into_development_dataframe",
                    "feature_profile": "hgvs",
                },
                "models": {
                    "without_missingness_flag": {
                        "artifact_path": str(self.artifact),
                        "artifact_sha256": artifact_hash,
                        "feature_columns": ["feature_a"],
                    }
                },
                "validation_metrics": {
                    "without_missingness_flag": {
                        "roc_auc": 0.982,
                        "average_precision": 0.981,
                    }
                },
            },
        )
        self.threshold_manifest_data = {
            "source": {
                "path": str(self.validation_predictions),
                "sha256": validation_hash,
            },
            "data_policy": {"test_data_used": False},
        }
        self.threshold_manifest = self.write_json(
            "threshold.json", self.threshold_manifest_data
        )
        self.cv_manifest = self.write_json(
            "cv.json",
            {
                "inputs": {
                    "features": {"sha256": self.feature_hash},
                    "gene_disjoint_split": {"sha256": self.split_hash},
                },
                "data_policy": {"test_used": False},
                "out_of_fold_metrics": {"roc_auc": 0.983},
            },
        )
        self.calibration_manifest = self.write_json(
            "calibration.json", {"data_policy": {"test_used": False}}
        )
        self.uncertainty_manifest = self.write_json(
            "uncertainty.json",
            {
                "inputs": {"features": {"sha256": self.feature_hash}},
                "data_policy": {"test_used": False},
            },
        )
        self.operating_points = self.root / "operating_points.csv"
        pd.DataFrame(
            [
                {
                    "model_name": "without_missingness_flag",
                    "policy": "sensitivity_at_least_0.95",
                    "sensitivity_target": 0.95,
                    "threshold": 0.36,
                    "sensitivity": 0.951,
                    "specificity": 0.888,
                    "precision_positive_predictive_value": 0.887,
                    "negative_predictive_value": 0.950,
                    "false_positives": 100,
                    "false_negatives": 50,
                }
            ]
        ).to_csv(self.operating_points, index=False)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def write_json(self, name: str, value: dict[str, object]) -> Path:
        path = self.root / name
        path.write_text(json.dumps(value), encoding="utf-8")
        return path

    def run_freeze(self) -> dict[str, object]:
        return freeze_protocol(
            self.config,
            self.model_manifest,
            self.threshold_manifest,
            self.operating_points,
            self.cv_manifest,
            self.calibration_manifest,
            self.uncertainty_manifest,
            self.root / "output",
        )

    def test_freeze_writes_lock_without_test_evaluation(self) -> None:
        lock = self.run_freeze()
        self.assertFalse(lock["test_evaluated"])
        self.assertFalse(lock["test_labels_read"])
        self.assertAlmostEqual(lock["frozen_operating_point"]["threshold"], 0.36)
        self.assertTrue((self.root / "output" / "protocol_lock.json").is_file())
        card = (self.root / "output" / "protocol_card.md").read_text(encoding="utf-8")
        self.assertIn("TEST AÇILMADI", card)

    def test_freeze_rejects_development_report_that_used_test(self) -> None:
        self.threshold_manifest_data["data_policy"]["test_data_used"] = True
        self.threshold_manifest = self.write_json(
            "threshold.json", self.threshold_manifest_data
        )
        with self.assertRaisesRegex(ValueError, "Test kilidi doğrulanamadı"):
            self.run_freeze()


if __name__ == "__main__":
    unittest.main()
