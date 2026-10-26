import tempfile
import unittest
from pathlib import Path

from scripts.prepare_grch38_reference import build_fai
from scripts.variant_normalization import ReferenceGenome
from tests.external_fixtures import write_reference


class PrepareGrch38ReferenceTest(unittest.TestCase):
    def test_fai_geometry_matches_crlf_wrapped_fasta_and_partial_last_line(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sequence = "CAAAAACGTACGTACGTACGTA"
            config = write_reference(root, sequence)
            fasta = root / "reference.fa"
            old_fai = Path(str(fasta) + ".fai").read_bytes()
            generated = build_fai(fasta, {"chr1": len(sequence)})
            self.assertEqual(generated.read_bytes(), old_fai)
            with ReferenceGenome(config) as reference:
                self.assertEqual(reference.fetch("1", 1, len(sequence)), sequence)
            with self.assertRaisesRegex(ValueError, "uzunluğu"):
                build_fai(fasta, {"chr1": len(sequence) - 1})
