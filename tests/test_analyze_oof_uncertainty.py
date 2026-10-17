import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.analyze_oof_uncertainty import analyze_uncertainty


class AnalyzeOofUncertaintyTest(unittest.TestCase):
    def test_analysis_reports_uncertain_errors_and_ignores_extra_features(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            predictions_path = root / "oof.tsv.gz"
            features_path = root / "features.tsv.gz"
            report_dir = root / "report"

            patterns = [
                (0, 0.05),
                (0, 0.30),
                (0, 0.60),
                (1, 0.20),
                (1, 0.70),
                (1, 0.95),
            ]
            prediction_rows = []
            feature_rows = []
            variation_id = 0
            for fold in (1, 2):
                for target, probability in patterns:
                    variation_id += 1
                    prediction_rows.append(
                        {
                            "VariationID": str(variation_id),
                            "gene_group": f"gene_{fold}_{variation_id}",
                            "source_partition": "train",
                            "fold": fold,
                            "target": target,
                            "probability_pathogenic": probability,
                            "prediction": int(probability >= 0.5),
                        }
                    )
                    feature_rows.append(self.feature_row(variation_id, target))
            feature_rows.append(self.feature_row(999, 0))
            pd.DataFrame(prediction_rows).to_csv(
                predictions_path, sep="\t", index=False, compression="gzip"
            )
            pd.DataFrame(feature_rows).to_csv(
                features_path, sep="\t", index=False, compression="gzip"
            )

            manifest = analyze_uncertainty(
                predictions_path,
                features_path,
                report_dir,
                min_group_size=2,
                chunk_size=5,
            )

            self.assertFalse(manifest["data_policy"]["fixed_validation_used"])
            self.assertFalse(manifest["data_policy"]["test_used"])
            self.assertEqual(manifest["integrity"]["oof_rows"], 12)
            self.assertEqual(manifest["counts"]["uncertain_rows"], 8)
            self.assertEqual(manifest["counts"]["uncertain_errors"], 4)
            self.assertEqual(manifest["counts"]["all_errors"], 4)
            self.assertEqual(manifest["counts"]["errors_outside_uncertain_band"], 0)
            errors = pd.read_csv(
                report_dir / "misclassified_oof_variants.tsv.gz", sep="\t"
            )
            self.assertEqual(set(errors["VariationID"]), {3, 4, 9, 10})
            self.assertNotIn(999, set(errors["VariationID"]))
            self.assertTrue((report_dir / "confident_error_subgroups.csv").is_file())
            self.assertTrue((report_dir / "uncertainty_overview.png").is_file())

    def test_analysis_rejects_target_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            predictions = []
            features = []
            for index, (target, probability, fold) in enumerate(
                ((0, 0.2, 1), (1, 0.8, 1), (0, 0.3, 2), (1, 0.7, 2)),
                start=1,
            ):
                predictions.append(
                    {
                        "VariationID": str(index),
                        "gene_group": f"gene_{index}",
                        "source_partition": "train",
                        "fold": fold,
                        "target": target,
                        "probability_pathogenic": probability,
                    }
                )
                features.append(self.feature_row(index, 1 - target if index == 1 else target))
            predictions_path = root / "oof.tsv.gz"
            features_path = root / "features.tsv.gz"
            pd.DataFrame(predictions).to_csv(
                predictions_path, sep="\t", index=False, compression="gzip"
            )
            pd.DataFrame(features).to_csv(
                features_path, sep="\t", index=False, compression="gzip"
            )

            with self.assertRaisesRegex(ValueError, "hedefleri eşleşmiyor"):
                analyze_uncertainty(
                    predictions_path,
                    features_path,
                    root / "report",
                    min_group_size=1,
                )

    @staticmethod
    def feature_row(variation_id: int, target: int) -> dict[str, object]:
        pathogenic = target == 1
        return {
            "VariationID": str(variation_id),
            "variant_type": "Deletion" if pathogenic else "SNV",
            "chromosome": "1",
            "mutation_class": "not_canonical_snv" if pathogenic else "transition",
            "hgvs_sequence_type": "coding_transcript",
            "dna_change_type": "deletion" if pathogenic else "substitution",
            "coding_region": "coding_exonic",
            "splice_offset_category": "not_intronic",
            "protein_effect": "frameshift" if pathogenic else "synonymous",
            "has_protein_annotation": 1,
            "is_canonical_snv": 0 if pathogenic else 1,
            "absolute_length_delta": 1 if pathogenic else 0,
            "target": target,
        }


if __name__ == "__main__":
    unittest.main()
