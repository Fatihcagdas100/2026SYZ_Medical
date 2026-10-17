import csv
import gzip
import json
import tempfile
import unittest
from pathlib import Path

from scripts.audit_variant_summary import audit_file


class AuditVariantSummaryTest(unittest.TestCase):
    def test_streaming_audit_outputs_and_repeat_counts(self) -> None:
        header = [
            "#AlleleID",
            "Type",
            "GeneSymbol",
            "ClinicalSignificance",
            "ClinSigSimple",
            "OriginSimple",
            "Assembly",
            "Chromosome",
            "ReviewStatus",
            "VariationID",
            "SomaticClinicalImpact",
            "Oncogenicity",
        ]
        rows = [
            ["10", "SNV", "GENE1", "Pathogenic", "1", "germline", "GRCh37", "1", "reviewed", "100", "-", "-"],
            ["10", "SNV", "GENE1", "Pathogenic", "1", "germline", "GRCh38", "1", "reviewed", "100", "-", "-"],
            ["11", "Deletion", "GENE2", "Benign", "0", "na", "GRCh38", "2", "", "101", "-", "-"],
        ]

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "input.txt.gz"
            output_dir = root / "report"
            with gzip.open(input_path, "wt", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(header)
                writer.writerows(rows)

            result = audit_file(
                input_path,
                output_dir,
                chunk_size=2,
                sample_size=2,
                progress_every=0,
            )

            self.assertEqual(result["row_count"], 3)
            self.assertEqual(result["column_count"], len(header))
            self.assertTrue(result["scan_complete"])
            self.assertEqual(
                result["variation_id_statistics"]["repeated_variation_ids"], 1
            )
            self.assertEqual(
                result["variation_id_statistics"][
                    "variation_ids_with_multiple_assemblies"
                ],
                1,
            )

            summary = json.loads((output_dir / "summary.json").read_text("utf-8"))
            self.assertEqual(summary["rows_scanned"], 3)
            self.assertEqual(summary["malformed_rows"], 0)
            self.assertEqual(
                len((output_dir / "sample_rows.tsv").read_text("utf-8").splitlines()),
                3,
            )


if __name__ == "__main__":
    unittest.main()

