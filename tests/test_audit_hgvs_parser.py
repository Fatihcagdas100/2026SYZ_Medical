import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from scripts.audit_hgvs_parser import audit_hgvs_parser
from scripts.build_hgvs_features import HGVS_FEATURE_COLUMNS, derive_hgvs_features


class AuditHgvsParserTest(unittest.TestCase):
    def test_audit_reproduces_features_and_excludes_test_examples(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source_path = root / "source.tsv.gz"
            features_path = root / "features.tsv.gz"
            split_path = root / "split.tsv"
            report_dir = root / "report"
            names = [
                ("1", "NM_1.1(G):c.3A>G (p.Lys1Arg)", "1", "train"),
                ("2", "NM_2.1(G):c.6A>G (p.Gly2=)", "0", "validation"),
                ("3", "NM_3.1(G):c.7del (p.Val3fs)", "1", "train"),
                ("4", "NM_4.1(G):c.9C>T (p.Gln3Ter)", "1", "validation"),
                ("5", "NM_5.1(G):c.10+1G>A", "0", "train"),
                ("6", "NM_6.1(G):c.12A>G (p.?)", "0", "validation"),
                ("8", "NM_8.1(G):c.12A>N (p.Lys4Xaa)", "0", "train"),
                ("7", "NM_TEST.1(G):c.3A>T (p.Lys1Ter)", "1", "test"),
            ]
            with gzip.open(source_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["VariationID", "Name", "target"])
                writer.writerows([row[:3] for row in names])
            with gzip.open(features_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["VariationID", *HGVS_FEATURE_COLUMNS, "target"],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for variation_id, name, target, _ in names:
                    writer.writerow(
                        {
                            "VariationID": variation_id,
                            **derive_hgvs_features(name),
                            "target": target,
                        }
                    )
            with split_path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["VariationID", "split"])
                writer.writerows([(row[0], row[3]) for row in names])

            manifest = audit_hgvs_parser(
                source_path,
                features_path,
                split_path,
                report_dir,
                samples_per_class=2,
            )
            self.assertFalse(manifest["data_policy"]["test_target_or_hgvs_features_used"])
            self.assertEqual(
                sum(manifest["saved_feature_reproduction_mismatches"].values()), 0
            )
            self.assertEqual(
                manifest["sample_audit"]["automatic_marker_mismatch_count"], 0
            )
            with (report_dir / "sampled_raw_name_examples.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                samples = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["split"] for row in samples}, {"train", "validation"})
            self.assertNotIn("NM_TEST", "\n".join(row["Name"] for row in samples))


if __name__ == "__main__":
    unittest.main()
