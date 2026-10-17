import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.build_missense_features import (
    AA_ORDER,
    BLOSUM62_ROWS,
    GRANTHAM_ROWS,
    build_missense_features,
    derive_missense_features,
    matrix_value,
)


class BuildMissenseFeaturesTest(unittest.TestCase):
    def test_known_substitutions_and_matrices(self) -> None:
        arg_to_gly = derive_missense_features(
            "NM_1.1(GENE):c.121A>G (p.Arg41Gly)", "missense"
        )
        self.assertEqual(arg_to_gly["missense_substitution"], "R>G")
        self.assertEqual(arg_to_gly["missense_grantham_distance"], 125)
        self.assertEqual(arg_to_gly["missense_blosum62_score"], -2)
        self.assertAlmostEqual(arg_to_gly["missense_hydropathy_delta"], 4.1)
        self.assertEqual(
            arg_to_gly["missense_residue_class_transition"], "positive>special"
        )
        self.assertEqual(
            arg_to_gly["missense_charge_transition"], "positive>neutral"
        )

        non_missense = derive_missense_features(
            "NM_2.1(GENE):c.6A>G (p.Gly2=)", "synonymous"
        )
        self.assertEqual(non_missense["missense_features_available"], 0)
        self.assertEqual(non_missense["missense_substitution"], "not_missense")

        self.assertEqual(len(GRANTHAM_ROWS), len(AA_ORDER))
        self.assertEqual(len(BLOSUM62_ROWS), len(AA_ORDER))
        for row in (*GRANTHAM_ROWS, *BLOSUM62_ROWS):
            self.assertEqual(len(row), len(AA_ORDER))
        for first in AA_ORDER:
            for second in AA_ORDER:
                self.assertEqual(
                    matrix_value(GRANTHAM_ROWS, first, second),
                    matrix_value(GRANTHAM_ROWS, second, first),
                )
                self.assertEqual(
                    matrix_value(BLOSUM62_ROWS, first, second),
                    matrix_value(BLOSUM62_ROWS, second, first),
                )

    def test_builder_preserves_rows_and_excludes_identity_text(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_path = root / "source.tsv.gz"
            features_path = root / "features.tsv.gz"
            output_path = root / "output.tsv.gz"
            report_dir = root / "report"
            with gzip.open(source_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["VariationID", "Name", "target"])
                writer.writerow(["1", "NM_1.1(G):c.121A>G (p.Arg41Gly)", "1"])
                writer.writerow(["2", "NM_2.1(G):c.6A>G (p.Gly2=)", "0"])
            with gzip.open(features_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["VariationID", "protein_effect", "target"])
                writer.writerow(["1", "missense", "1"])
                writer.writerow(["2", "synonymous", "0"])

            manifest = build_missense_features(
                source_path,
                features_path,
                output_path,
                report_dir,
                progress_every=0,
            )
            with gzip.open(output_path, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                rows = list(reader)
                self.assertNotIn("Name", reader.fieldnames or [])
            self.assertEqual(len(rows), 2)
            self.assertEqual(rows[0]["missense_substitution"], "R>G")
            self.assertEqual(rows[1]["missense_substitution"], "not_missense")
            self.assertEqual(manifest["counts"]["parsed_missense_rows"], 1)
            self.assertFalse(manifest["feature_policy"]["protein_position_in_output"])


if __name__ == "__main__":
    unittest.main()
