import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.freeze_external_candidate import freeze_external_candidate


class FreezeExternalCandidateTest(unittest.TestCase):
    def test_refit_uses_development_only_and_writes_auditable_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            features_path = root / "features.tsv.gz"
            split_path = root / "split.tsv"
            rows = []
            splits = []
            variation_id = 0
            for partition, count_per_class in (("train", 4), ("validation", 2), ("test", 2)):
                for target in (0, 1):
                    for _ in range(count_per_class):
                        variation_id += 1
                        rows.append(
                            {
                                "VariationID": str(variation_id),
                                "variant_type": "Deletion" if target else "SNV",
                                "chromosome": "1",
                                "substitution": "NOT_CANONICAL_SNV" if target else "A>G",
                                "mutation_class": "not_canonical_snv" if target else "transition",
                                "hgvs_sequence_type": "coding_transcript",
                                "dna_change_type": "substitution",
                                "coding_region": "coding_exonic",
                                "splice_offset_category": "not_intronic",
                                "protein_effect": "missense",
                                "reference_length": 2 if target else 1,
                                "alternate_length": 1,
                                "length_delta": -1 if target else 0,
                                "absolute_length_delta": 1 if target else 0,
                                "is_canonical_snv": 0 if target else 1,
                                "has_protein_annotation": 1,
                                "missense_reference_residue": "Arg",
                                "missense_alternate_residue": "Gly",
                                "missense_substitution": "Arg>Gly",
                                "missense_residue_class_transition": "basic>special",
                                "missense_charge_transition": "positive>neutral",
                                "missense_features_available": 1,
                                "missense_grantham_distance": 125,
                                "missense_blosum62_score": -2,
                                "missense_hydropathy_delta": 4.1,
                                "missense_absolute_hydropathy_delta": 4.1,
                                "target": target,
                            }
                        )
                        splits.append(
                            {"VariationID": str(variation_id), "split": partition}
                        )
            pd.DataFrame(rows).to_csv(features_path, sep="\t", index=False, compression="gzip")
            pd.DataFrame(splits).to_csv(split_path, sep="\t", index=False)

            manifest = freeze_external_candidate(
                features_path,
                split_path,
                root / "report",
                root / "artifacts",
                chunk_size=5,
            )

            self.assertFalse(manifest["data_policy"]["locked_test_used"])
            self.assertEqual(manifest["partition_counts"], {"train": 8, "validation": 4, "test": 4})
            self.assertTrue(Path(manifest["model"]["artifact_path"]).is_file())
            self.assertTrue((root / "report" / "candidate_card.md").is_file())


if __name__ == "__main__":
    unittest.main()
