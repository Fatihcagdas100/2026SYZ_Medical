import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.build_hgvs_features import (
    BASELINE_COLUMNS,
    build_hgvs_features,
    derive_hgvs_features,
)


class BuildHgvsFeaturesTest(unittest.TestCase):
    def test_hgvs_consequence_classification_ignores_gene_name_substrings(self) -> None:
        cases = {
            "NM_000001.1(WFS1):c.123A>G (p.Arg41Gly)": {
                "hgvs_sequence_type": "coding_transcript",
                "dna_change_type": "substitution",
                "coding_region": "coding_exonic",
                "splice_offset_category": "not_intronic",
                "protein_effect": "missense",
            },
            "NM_000002.1(GENE):c.456C>T (p.Arg152Ter)": {
                "protein_effect": "stop_gained_or_changed"
            },
            "NM_000003.1(GENE):c.100+1G>A": {
                "coding_region": "intronic",
                "splice_offset_category": "donor_core_1_2",
                "protein_effect": "no_protein_annotation",
            },
            "NM_000004.1(GENE):c.90-7A>G": {
                "splice_offset_category": "splice_region_3_8"
            },
            "NM_000005.1(GENE):c.-12A>G": {"coding_region": "five_prime_utr"},
            "NM_000006.1(GENE):c.*12A>G": {"coding_region": "three_prime_utr"},
            "NM_000007.1(GENE):c.300A>G (p.Gly100=)": {
                "protein_effect": "synonymous"
            },
            "NM_000008.1(GENE):c.123del (p.Val42fs)": {
                "dna_change_type": "deletion",
                "protein_effect": "frameshift",
            },
            "NC_000001.1:g.123A>G": {
                "hgvs_sequence_type": "genomic",
                "coding_region": "not_coding_transcript",
            },
        }
        for name, expected in cases.items():
            with self.subTest(name=name):
                observed = derive_hgvs_features(name)
                for feature, value in expected.items():
                    self.assertEqual(observed[feature], value)

    def test_builder_preserves_alignment_and_excludes_raw_name(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_path = root / "source.tsv.gz"
            baseline_path = root / "baseline.tsv.gz"
            split_path = root / "split.tsv"
            output_path = root / "output.tsv.gz"
            report_dir = root / "report"

            with gzip.open(source_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["VariationID", "Name", "target"])
                writer.writerow(["1", "NM_1.1(GENE):c.1A>G (p.Met1Val)", "1"])
                writer.writerow(["2", "NM_2.1(GENE):c.6A>G (p.Gly2=)", "0"])
                writer.writerow(["3", "NM_3.1(GENE):c.7C>T (p.Gln3Ter)", "1"])

            with gzip.open(baseline_path, "wt", encoding="utf-8", newline="") as handle:
                fieldnames = [*BASELINE_COLUMNS, "target"]
                writer = csv.DictWriter(
                    handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                for variation_id, target in (("1", "1"), ("2", "0"), ("3", "1")):
                    row = {column: "value" for column in BASELINE_COLUMNS}
                    row["VariationID"] = variation_id
                    row["target"] = target
                    writer.writerow(row)

            with split_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["VariationID", "split"])
                writer.writerow(["1", "train"])
                writer.writerow(["2", "validation"])
                writer.writerow(["3", "test"])

            manifest = build_hgvs_features(
                source_path,
                baseline_path,
                split_path,
                output_path,
                report_dir,
                progress_every=0,
            )
            with gzip.open(output_path, "rt", encoding="utf-8", newline="") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                rows = list(reader)
                self.assertNotIn("Name", reader.fieldnames or [])
                self.assertNotIn("GeneSymbol", reader.fieldnames or [])
            self.assertEqual(rows[0]["protein_effect"], "start_lost_or_uncertain")
            self.assertEqual(rows[1]["protein_effect"], "synonymous")
            self.assertFalse(manifest["feature_policy"]["raw_name_in_output"])
            self.assertFalse(
                manifest["feature_policy"]["test_target_distribution_reported"]
            )
            distribution_path = (
                report_dir / "categorical_distributions_by_development_partition.csv"
            )
            with distribution_path.open(encoding="utf-8", newline="") as handle:
                distributions = list(csv.DictReader(handle))
            self.assertEqual({row["split"] for row in distributions}, {"train", "validation"})


if __name__ == "__main__":
    unittest.main()
