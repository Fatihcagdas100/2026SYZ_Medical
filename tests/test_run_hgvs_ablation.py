import csv
import gzip
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.run_hgvs_ablation import run_hgvs_ablation


class RunHgvsAblationTest(unittest.TestCase):
    def test_ablation_uses_validation_and_writes_expected_subgroups(self) -> None:
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
                "has_protein_annotation", "target",
            ]
            effects = [
                ("missense", 1), ("synonymous", 1), ("frameshift", 1),
                ("stop_gained_or_changed", 1), ("no_protein_annotation", 0),
                ("missense", 1), ("synonymous", 1), ("frameshift", 1),
            ]
            allocation = {"train": 16, "validation": 8, "test": 4}
            rows = []
            splits = []
            variation_id = 0
            for split_name, count in allocation.items():
                for index in range(count):
                    variation_id += 1
                    target = index % 2
                    effect, annotated = effects[index % len(effects)]
                    rows.append(
                        {
                            "VariationID": str(variation_id),
                            "variant_type": "Deletion" if target else "single nucleotide variant",
                            "chromosome": "1" if target else "2",
                            "reference_length": 2 if target else 1,
                            "alternate_length": 1,
                            "length_delta": -1 if target else 0,
                            "absolute_length_delta": 1 if target else 0,
                            "vcf_alleles_missing": 0,
                            "is_canonical_snv": 0 if target else 1,
                            "substitution": "NOT_CANONICAL_SNV" if target else "A>G",
                            "mutation_class": "not_canonical_snv" if target else "transition",
                            "hgvs_sequence_type": "coding_transcript" if annotated else "genomic",
                            "dna_change_type": "deletion" if target else "substitution",
                            "coding_region": "coding_exonic" if annotated else "not_coding_transcript",
                            "splice_offset_category": "not_intronic" if annotated else "not_coding_transcript",
                            "protein_effect": effect,
                            "has_protein_annotation": annotated,
                            "target": target,
                        }
                    )
                    splits.append({"VariationID": str(variation_id), "target": target, "split": split_name})
            with gzip.open(features_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            with split_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["VariationID", "target", "split"], delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(splits)

            manifest = run_hgvs_ablation(
                features_path, split_path, report_dir, artifact_dir, chunk_size=5
            )
            self.assertEqual(
                manifest["data_policy"]["test_partition"],
                "locked_not_loaded_into_development_dataframe",
            )
            metrics = pd.read_csv(report_dir / "ablation_metrics.csv")
            self.assertEqual(len(metrics), 3)
            subgroups = pd.read_csv(report_dir / "full_model_validation_subgroups.csv")
            self.assertEqual(
                set(subgroups["group"]),
                {
                    "all_variants", "missense", "synonymous", "frameshift",
                    "stop_gained_or_changed", "protein_annotation_present",
                    "protein_annotation_absent",
                },
            )
            self.assertEqual(
                int(subgroups.loc[subgroups["group"] == "all_variants", "row_count"].iloc[0]),
                8,
            )


if __name__ == "__main__":
    unittest.main()
