import tempfile
import unittest
from pathlib import Path

import pandas as pd

from scripts.admit_external_cohort import admit_external_cohort
from tests.external_fixtures import make_row, write_reference, write_rows


class AdmitExternalCohortTest(unittest.TestCase):
    def test_normalized_indel_overlap_bad_build_and_reference_are_excluded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            reference = write_reference(root)
            development, external = root / "dev.tsv.gz", root / "external.tsv.gz"
            write_rows(development, [make_row(1, 2, 0, reference="AA", alternate="A")])
            write_rows(external, [
                make_row(2, 4, 1, reference="AA", alternate="A"),
                make_row(3, 7, 1, reference="C", alternate="T", assembly="GRCh37"),
                make_row(4, 8, 0, reference="T", alternate="G"),
                make_row(5, 7, 0, reference="C", alternate="T"),
                make_row(6, 8, 1, reference="G", alternate="A"),
            ])
            output = root / "admitted.tsv.gz"
            manifest = admit_external_cohort(
                development, external, output, root / "report", cohort_kind="independent",
                source_version="synthetic", source_url="https://example.org/test",
                reference_config=reference,
            )
            self.assertEqual(set(pd.read_csv(output, sep="\t")["VariationID"]), {5, 6})
            self.assertEqual(manifest["counts"]["excluded_coordinate_allele"], 1)
            self.assertEqual(manifest["counts"]["excluded_assembly_not_grch38"], 1)
            self.assertEqual(manifest["counts"]["excluded_reference_allele_mismatch"], 1)
            self.assertEqual(manifest["reference"]["assembly"], "GRCh38")

    def test_duplicate_equivalent_external_indels_abort(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = write_reference(root)
            development, external = root / "dev.tsv.gz", root / "external.tsv.gz"
            write_rows(development, [make_row(1, 7, 0, reference="C", alternate="T")])
            write_rows(external, [make_row(2, 2, 0, reference="AA", alternate="A"), make_row(3, 4, 1, reference="AA", alternate="A")])
            with self.assertRaisesRegex(ValueError, "normalleştirilmiş"):
                admit_external_cohort(development, external, root / "output.tsv.gz", root / "report", cohort_kind="independent", source_version="test", source_url="https://example.org", reference_config=config)
            self.assertFalse((root / "output.tsv.gz").exists())
