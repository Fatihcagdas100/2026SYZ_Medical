"""Kabul → özellik → dondurulmuş model → rapor hattı; gerçek veri kullanılmaz."""

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from scripts.admit_external_cohort import admit_external_cohort
from scripts.build_external_features import build_external_features
from scripts.evaluate_external_candidate import evaluate_external_candidate, subgroup_result
from scripts.freeze_external_candidate import freeze_external_candidate
from tests.external_fixtures import make_row, write_reference, write_rows


class ExternalPipelineTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        sequence = "CAAAAACGTACGTACGTACGTACGT" * 20
        reference_config = write_reference(self.root, sequence)
        self.development = self.root / "development.tsv.gz"
        development_rows = []
        for index in range(16):
            position = 11 + index
            ref = sequence[position - 1]
            development_rows.append(make_row(index + 1, position, index % 2, reference=ref, alternate="T" if ref != "T" else "G"))
        development_rows[0].update({"PositionVCF": "2", "ReferenceAlleleVCF": "AA", "AlternateAlleleVCF": "A"})
        write_rows(self.development, development_rows)
        development_features = self.root / "development_features.tsv.gz"
        split = self.root / "split.tsv"
        pd.DataFrame({"VariationID": [str(i) for i in range(1, 17)], "split": ["train"] * 8 + ["validation"] * 4 + ["test"] * 4}).to_csv(split, sep="\t", index=False)
        self.raw_external = self.root / "external.tsv.gz"
        external_rows = []
        for index in range(8):
            position = 201 + index
            ref = sequence[position - 1]
            row = make_row(1001 + index, position, index % 2, reference=ref, alternate="T" if ref != "T" else "G", gene="EXTERNAL" if index < 6 else "DEV")
            if index >= 6:
                row["Name"] = "-"
            external_rows.append(row)
        external_rows.append(make_row(2000, 4, 1, reference="AA", alternate="A"))
        write_rows(self.raw_external, external_rows)
        self.admitted = self.root / "admitted.tsv.gz"
        self.features = self.root / "features.tsv.gz"
        self.candidate = self.root / "candidate" / "manifest.json"
        self.admission = self.root / "admission" / "manifest.json"
        self.feature_manifest = self.root / "feature_report" / "manifest.json"
        with contextlib.redirect_stdout(io.StringIO()):
            build_external_features(self.development, development_features, self.root / "dev_feature_report")
            freeze_external_candidate(development_features, split, self.candidate.parent, self.root / "artifacts")
            admit_external_cohort(self.development, self.raw_external, self.admitted, self.admission.parent, cohort_kind="independent", source_version="synthetic-v1", source_url="https://example.org/synthetic", reference_config=reference_config)
            build_external_features(self.admitted, self.features, self.feature_manifest.parent)

    def evaluate(self, report_name="evaluation"):
        with contextlib.redirect_stdout(io.StringIO()):
            return evaluate_external_candidate(self.candidate, self.features, self.admission, self.feature_manifest, self.root / report_name)

    def test_end_to_end_reports_missense_confusion_intervals_and_coverage(self):
        result = self.evaluate()
        self.assertEqual(result["metrics"]["row_count"], 8)
        self.assertEqual(result["subgroups"]["missense"]["metrics"]["row_count"], 6)
        self.assertEqual(result["subgroups"]["unseen_genes"]["metrics"]["row_count"], 6)
        self.assertEqual(result["annotation_coverage"]["protein_annotation_missing_rows"], 2)
        self.assertEqual(result["bootstrap"]["replicates"], 1000)
        self.assertEqual(result["frozen_threshold"], 0.5)
        self.assertIn("specificity", result["subgroups"]["missense"]["bootstrap_95_ci"])
        report = (self.root / "evaluation" / "external_evaluation.md").read_text(encoding="utf-8")
        self.assertIn("missense", report)
        self.assertIn("TN=", report)
        with patch("scripts.evaluate_external_candidate.pd.read_csv") as read:
            with self.assertRaises(FileExistsError):
                self.evaluate("different_report")
            read.assert_not_called()

    def test_interrupted_evaluation_keeps_lock_before_model_prediction(self):
        with patch("scripts.evaluate_external_candidate.joblib.load", side_effect=RuntimeError("simulated interruption")):
            with self.assertRaisesRegex(RuntimeError, "simulated"):
                self.evaluate()
        self.assertTrue((self.root / "evaluation" / "EVALUATION_STARTED.json").is_file())
        self.assertFalse((self.root / "evaluation" / "manifest.json").exists())
        with patch("scripts.evaluate_external_candidate.pd.read_csv") as read:
            with self.assertRaises(FileExistsError):
                self.evaluate("another_report")
            read.assert_not_called()

    def test_tampered_feature_file_rejected_before_opening_labels(self):
        with self.features.open("ab") as handle:
            handle.write(b"tampered")
        with patch("scripts.evaluate_external_candidate.pd.read_csv") as read:
            with self.assertRaisesRegex(ValueError, "hash"):
                self.evaluate()
            read.assert_not_called()
        self.assertFalse((self.candidate.parent / "evaluation_registry").exists())

    def test_missing_admission_hash_rejected(self):
        admission = json.loads(self.admission.read_text(encoding="utf-8"))
        del admission["output"]["sha256"]
        self.admission.write_text(json.dumps(admission), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            self.evaluate()

    def test_tampered_artifact_rejected(self):
        candidate = json.loads(self.candidate.read_text(encoding="utf-8"))
        Path(candidate["model"]["artifact_path"]).write_bytes(b"changed artifact")
        with self.assertRaisesRegex(ValueError, "hash"):
            self.evaluate()
        self.assertFalse((self.candidate.parent / "evaluation_registry").exists())

    def test_single_class_subgroup_explicitly_skips_metrics(self):
        frame = pd.DataFrame({"target": [1, 1], "probability_pathogenic": [0.8, 0.9]})
        result = subgroup_result(frame, threshold=0.5, replicates=1000, seed=2027)
        self.assertEqual(result["status"], "insufficient_classes")
        self.assertIsNone(result["metrics"])

    def test_frozen_candidate_cannot_be_overwritten(self):
        with self.assertRaises(FileExistsError):
            freeze_external_candidate(
                self.root / "development_features.tsv.gz", self.root / "split.tsv",
                self.candidate.parent, self.root / "artifacts",
            )
