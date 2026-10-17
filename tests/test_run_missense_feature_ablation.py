import csv
import gzip
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.run_missense_feature_ablation import run_missense_feature_ablation


class RunMissenseFeatureAblationTest(unittest.TestCase):
    def test_ablation_keeps_test_locked_and_compares_missense(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            features_path = root / "features.tsv.gz"
            split_path = root / "split.tsv"
            report_dir = root / "reports"
            artifact_dir = root / "artifacts"
            fieldnames = [
                "VariationID", "variant_type", "chromosome", "reference_length",
                "alternate_length", "length_delta", "absolute_length_delta",
                "vcf_alleles_missing", "is_canonical_snv", "substitution",
                "mutation_class", "hgvs_sequence_type", "dna_change_type",
                "coding_region", "splice_offset_category", "protein_effect",
                "has_protein_annotation", "missense_reference_residue",
                "missense_alternate_residue", "missense_substitution",
                "missense_residue_class_transition", "missense_charge_transition",
                "missense_features_available", "missense_grantham_distance",
                "missense_blosum62_score", "missense_hydropathy_delta",
                "missense_absolute_hydropathy_delta", "target",
            ]
            allocation = {"train": 24, "validation": 8, "test": 4}
            rows = []
            split_rows = []
            variation_id = 0
            for split_name, count in allocation.items():
                for index in range(count):
                    variation_id += 1
                    target = index % 2
                    missense = index % 4 in {0, 1}
                    rows.append(
                        {
                            "VariationID": str(variation_id),
                            "variant_type": "single nucleotide variant",
                            "chromosome": "1" if target else "2",
                            "reference_length": 1, "alternate_length": 1,
                            "length_delta": 0, "absolute_length_delta": 0,
                            "vcf_alleles_missing": 0, "is_canonical_snv": 1,
                            "substitution": "A>G", "mutation_class": "transition",
                            "hgvs_sequence_type": "coding_transcript",
                            "dna_change_type": "substitution",
                            "coding_region": "coding_exonic",
                            "splice_offset_category": "not_intronic",
                            "protein_effect": "missense" if missense else "synonymous",
                            "has_protein_annotation": 1,
                            "missense_reference_residue": "R" if missense else "not_missense",
                            "missense_alternate_residue": "G" if missense else "not_missense",
                            "missense_substitution": "R>G" if missense else "not_missense",
                            "missense_residue_class_transition": "positive>special" if missense else "not_missense",
                            "missense_charge_transition": "positive>neutral" if missense else "not_missense",
                            "missense_features_available": int(missense),
                            "missense_grantham_distance": 125 if missense else 0,
                            "missense_blosum62_score": -2 if missense else 0,
                            "missense_hydropathy_delta": 4.1 if missense else 0,
                            "missense_absolute_hydropathy_delta": 4.1 if missense else 0,
                            "target": target,
                        }
                    )
                    split_rows.append(
                        {"VariationID": str(variation_id), "target": target, "split": split_name}
                    )
            with gzip.open(features_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            with split_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["VariationID", "target", "split"], delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(split_rows)

            manifest = run_missense_feature_ablation(
                features_path, split_path, report_dir, artifact_dir, chunk_size=7
            )
            self.assertEqual(
                manifest["data_policy"]["test_partition"],
                "locked_not_loaded_into_development_dataframe",
            )
            comparison = pd.read_csv(report_dir / "validation_model_comparison.csv")
            self.assertEqual(len(comparison), 2)
            subgroups = pd.read_csv(report_dir / "validation_subgroup_comparison.csv")
            missense = subgroups.loc[subgroups["group"] == "missense"]
            self.assertEqual(len(missense), 2)
            self.assertEqual(set(missense["row_count"]), {4})


if __name__ == "__main__":
    unittest.main()
