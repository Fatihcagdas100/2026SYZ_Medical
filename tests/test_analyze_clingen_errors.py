import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from scripts.analyze_clingen_errors import (
    analyze_errors, feature_collisions, join_saved_tables, subgroup_tables, training_coverage,
    verified_file, wilson_interval,
)
from scripts.train_logistic_baseline import sha256_file
from scripts.download_clingen_erepo import download_snapshot
from tests.test_clingen_erepo_screen import curation, source_stub


class AnalyzeClinGenErrorsTest(unittest.TestCase):
    def tables(self):
        predictions = pd.DataFrame({
            "VariationID": ["1", "2", "3", "4"], "target": [1, 0, 1, 0],
            "probability_pathogenic": [0.1, 0.1, 0.8, 0.8], "prediction": [0, 0, 1, 1],
        })
        features = pd.DataFrame({
            "VariationID": ["3", "1", "4", "2"], "target": [1, 1, 0, 0],
            "protein_effect": ["missense"] * 4, "has_protein_annotation": [1] * 4,
            "missense_substitution": ["A>T", "R>H", "A>T", "R>H"],
            "missense_residue_class_transition": ["other"] * 4,
            "missense_charge_transition": ["other"] * 4,
            "variant_type": ["SNV"] * 4, "coding_region": ["coding"] * 4,
            "hgvs_sequence_type": ["coding_transcript"] * 4,
        })
        cohort = pd.DataFrame({
            "VariationID": ["2", "4", "1", "3"], "target": [0, 0, 1, 1],
            "Name": ["name"] * 4, "GeneSymbol": ["G2", "G4", "G1", "G3"],
            "annotation_source": ["source;alleles"] * 4,
            "clingen_panels_json": ["[]"] * 4, "original_clinvar_labels_json": ["[]"] * 4,
        })
        return predictions, features, cohort

    def test_id_join_and_denominators_and_conflicting_feature_vectors(self):
        data = join_saved_tables(*self.tables(), threshold=0.5)
        self.assertEqual(int(data["FN"].sum()), 1)
        self.assertEqual(int(data["FP"].sum()), 1)
        self.assertEqual(data.loc[data["VariationID"].eq("1"), "GeneSymbol"].iloc[0], "G1")
        data["fit_pair_support"] = "100+"
        data["has_unknown_fit_category"] = 0
        groups = subgroup_tables(data, min_class_count=2)
        rh = groups.loc[groups["scope"].eq("missense") & groups["feature"].eq("missense_substitution") & groups["value"].eq("R>H")].iloc[0]
        self.assertEqual(rh["positives"], 1)
        self.assertEqual(rh["negatives"], 1)
        self.assertEqual(rh["false_negative_rate"], 1.0)
        self.assertFalse(rh["fnr_support_at_least_minimum"])
        collisions = feature_collisions(data, ["protein_effect", "missense_substitution"])
        self.assertEqual(len(collisions), 2)
        self.assertEqual(int(collisions["rows"].sum()), 4)

    def test_targets_ids_and_frozen_decisions_are_validated(self):
        predictions, features, cohort = self.tables()
        features.loc[0, "target"] = 0
        with self.assertRaisesRegex(ValueError, "hedefleri"):
            join_saved_tables(predictions, features, cohort, threshold=0.5)
        predictions, features, cohort = self.tables()
        predictions.loc[0, "prediction"] = 1
        with self.assertRaisesRegex(ValueError, "eşikle"):
            join_saved_tables(predictions, features, cohort, threshold=0.5)
        predictions, features, cohort = self.tables()
        with self.assertRaisesRegex(ValueError, "ID kümesi"):
            join_saved_tables(predictions, features.iloc[:3], cohort, threshold=0.5)

    def test_training_coverage_excludes_test_only_categories_without_target_column(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fit = root / "features.tsv.gz"
            split = root / "split.tsv"
            pd.DataFrame({"VariationID": ["10", "11", "12"], "protein_effect": ["missense"] * 3,
                          "missense_substitution": ["R>H", "A>T", "W>C"]}).to_csv(fit, sep="\t", index=False, compression="gzip")
            pd.DataFrame({"VariationID": ["10", "11", "12"], "split": ["train", "validation", "test"]}).to_csv(split, sep="\t", index=False)
            candidate = {"inputs": {"features": {"path": str(fit), "sha256": sha256_file(fit)},
                                    "gene_disjoint_split": {"path": str(split), "sha256": sha256_file(split)}},
                         "model": {"categorical_columns": ["protein_effect", "missense_substitution"]}}
            data = pd.DataFrame({"protein_effect": ["missense", "missense"], "missense_substitution": ["R>H", "W>C"]})
            covered, counts = training_coverage(data, candidate)
            self.assertEqual(covered["fit_pair_count"].tolist(), [1, 0])
            self.assertEqual(covered["has_unknown_fit_category"].tolist(), [0, 1])
            self.assertEqual(counts["fit_rows"], 2)
            self.assertFalse(counts["test_targets_read"])
            with fit.open("ab") as handle:
                handle.write(b"tampered")
            with self.assertRaisesRegex(ValueError, "hash"):
                verified_file(candidate["inputs"]["features"])

    def test_wilson_handles_absent_class_and_extreme_rates(self):
        self.assertEqual(wilson_interval(0, 0), (None, None))
        lower, upper = wilson_interval(0, 5)
        self.assertAlmostEqual(lower, 0.0)
        self.assertGreater(upper, 0.4)
        lower, upper = wilson_interval(5, 5)
        self.assertLess(lower, 0.6)
        self.assertAlmostEqual(upper, 1.0)

    def test_complete_saved_analysis_never_loads_model_and_preserves_primary_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            predictions, features, cohort = self.tables()
            def info(path):
                return {"path": str(path), "sha256": sha256_file(path)}
            files = {}
            for name, frame in (("predictions", predictions), ("features", features), ("cohort", cohort)):
                path = root / f"{name}.tsv.gz"
                frame.to_csv(path, sep="\t", index=False, compression="gzip")
                files[name] = info(path)
            fit_path = root / "fit.tsv.gz"
            features[["VariationID", "protein_effect", "missense_substitution"]].to_csv(fit_path, sep="\t", index=False, compression="gzip")
            split_path = root / "split.tsv"
            pd.DataFrame({"VariationID": ["1", "2", "3", "4"], "split": ["train", "train", "validation", "test"]}).to_csv(split_path, sep="\t", index=False)
            candidate = {"threshold": 0.5, "inputs": {"features": info(fit_path), "gene_disjoint_split": info(split_path)},
                         "model": {"categorical_columns": ["protein_effect", "missense_substitution"], "feature_columns": ["protein_effect", "missense_substitution"]}}
            candidate_path, admission_path, evaluation_path = (root / name for name in ("candidate.json", "admission.json", "evaluation.json"))
            candidate_path.write_text(json.dumps(candidate))
            admission_path.write_text(json.dumps({"output": files["cohort"]}))
            evaluation_path.write_text(json.dumps({"evaluation_number": 1, "frozen_threshold": 0.5,
                "external_features": files["features"], "admission_manifest": info(admission_path),
                "candidate_manifest": info(candidate_path), "outputs": {"predictions": files["predictions"]["path"]},
                "metrics": {"row_count": 4, "false_negatives": 1, "false_positives": 1}}))
            rows = [curation(i, str(i), "Pathogenic" if i in (1, 3) else "Benign") for i in (1, 2, 3, 4)]
            snapshot = root / "snapshot"
            with patch("scripts.download_clingen_erepo.fetch_bytes", side_effect=source_stub(rows)), contextlib.redirect_stdout(io.StringIO()):
                download_snapshot(snapshot, page_size=2)
            lock = {"analysis_type": "post_evaluation_descriptive_error_analysis", "evaluation_manifest": info(evaluation_path),
                    "predictions": files["predictions"], "curation_snapshot_manifest": info(snapshot / "manifest.json"),
                    "min_class_count": 2, "probability_regions": [0.2, 0.8]}
            lock_path = root / "lock.json"
            lock_path.write_text(json.dumps(lock))
            with patch("joblib.load", side_effect=AssertionError("Model must not load")) as load, contextlib.redirect_stdout(io.StringIO()):
                result = analyze_errors(lock_path, root / "report")
                load.assert_not_called()
            self.assertEqual(result["counts"]["false_negatives"], 1)
            self.assertFalse(result["data_policy"]["threshold_search"])
            self.assertEqual(info(evaluation_path), lock["evaluation_manifest"])
            self.assertEqual(info(Path(files["predictions"]["path"])), files["predictions"])
            self.assertTrue((root / "report" / "error_curation_evidence.tsv").is_file())
            with self.assertRaisesRegex(ValueError, "birincil girdi"):
                analyze_errors(lock_path, root)


if __name__ == "__main__":
    unittest.main()
