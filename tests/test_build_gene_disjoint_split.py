import csv
import gzip
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

from scripts.build_gene_disjoint_split import build_gene_disjoint_split


class BuildGeneDisjointSplitTest(unittest.TestCase):
    def test_split_is_deterministic_stratified_and_gene_disjoint(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "input.tsv.gz"
            output_1 = root / "split_1.tsv"
            output_2 = root / "split_2.tsv"

            with gzip.open(input_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["VariationID", "GeneSymbol", "target"])
                variation_id = 0
                for group_index in range(40):
                    symbols = (
                        ("GENE_A;GENE_B", "GENE_B;GENE_C")
                        if group_index == 0
                        else (f"GENE_{group_index}", f"GENE_{group_index}")
                    )
                    for target, symbol in enumerate(symbols):
                        variation_id += 1
                        writer.writerow([variation_id, symbol, target])

            manifest_1 = build_gene_disjoint_split(
                input_path,
                output_1,
                root / "report_1",
                seed=2027,
            )
            manifest_2 = build_gene_disjoint_split(
                input_path,
                output_2,
                root / "report_2",
                seed=2027,
            )

            self.assertEqual(output_1.read_bytes(), output_2.read_bytes())
            self.assertEqual(
                manifest_1["output"]["sha256"], manifest_2["output"]["sha256"]
            )
            self.assertEqual(manifest_1["integrity"]["gene_group_overlap_count"], 0)

            rows = []
            with output_1.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            splits_by_group: defaultdict[str, set[str]] = defaultdict(set)
            targets_by_split: defaultdict[str, set[str]] = defaultdict(set)
            for row in rows:
                splits_by_group[row["gene_group"]].add(row["split"])
                targets_by_split[row["split"]].add(row["target"])
            self.assertTrue(all(len(splits) == 1 for splits in splits_by_group.values()))
            self.assertEqual({row["gene_group"] for row in rows[:2]}, {"GENE_A"})
            self.assertEqual(set(targets_by_split), {"train", "validation", "test"})
            self.assertTrue(all(targets == {"0", "1"} for targets in targets_by_split.values()))


if __name__ == "__main__":
    unittest.main()
