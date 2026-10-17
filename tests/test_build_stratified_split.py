import csv
import gzip
import tempfile
import unittest
from collections import Counter, defaultdict
from pathlib import Path

from scripts.build_stratified_split import build_split


class BuildStratifiedSplitTest(unittest.TestCase):
    def test_split_is_stratified_disjoint_and_deterministic(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "features.tsv.gz"
            with gzip.open(input_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["VariationID", "target"])
                for target in ("0", "1"):
                    for index in range(6):
                        writer.writerow([f"{target}-{index}", target])

            output_1 = root / "split_1.tsv"
            output_2 = root / "split_2.tsv"
            manifest_1 = build_split(input_path, output_1, root / "report_1", seed=2027)
            manifest_2 = build_split(input_path, output_2, root / "report_2", seed=2027)

            self.assertEqual(output_1.read_bytes(), output_2.read_bytes())
            self.assertEqual(
                manifest_1["output"]["sha256"], manifest_2["output"]["sha256"]
            )

            by_target: defaultdict[str, Counter[str]] = defaultdict(Counter)
            ids_by_split: defaultdict[str, set[str]] = defaultdict(set)
            with output_1.open(encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    by_target[row["target"]][row["split"]] += 1
                    ids_by_split[row["split"]].add(row["VariationID"])

            for target in ("0", "1"):
                self.assertEqual(
                    by_target[target],
                    Counter({"train": 4, "validation": 1, "test": 1}),
                )
            self.assertFalse(ids_by_split["train"] & ids_by_split["validation"])
            self.assertFalse(ids_by_split["train"] & ids_by_split["test"])
            self.assertFalse(ids_by_split["validation"] & ids_by_split["test"])


if __name__ == "__main__":
    unittest.main()

