import tempfile
import unittest
from pathlib import Path

from scripts.variant_normalization import ReferenceGenome
from tests.external_fixtures import make_row, write_reference


class VariantNormalizationTest(unittest.TestCase):
    def test_alias_padding_and_tandem_repeat_normalize_equally(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with ReferenceGenome(write_reference(root, "CATATATACGT")) as reference:
                key = reference.normalize(make_row(1, 1, 0, reference="CAT", alternate="C"))
                shifted = make_row(2, 4, 0, reference="ATA", alternate="A")
                shifted["Chromosome"] = "chr1"
                self.assertEqual(key, reference.normalize(shifted))
                self.assertEqual(reference.fetch("1", 6, 5), "ATACG")

    def test_normalization_preserves_haplotype_and_is_idempotent(self):
        sequence = "CAAAAACGT"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with ReferenceGenome(write_reference(root, sequence)) as reference:
                canonical_by_haplotype = {}
                for position in range(1, len(sequence) + 1):
                    for length in range(1, min(3, len(sequence) - position + 1) + 1):
                        ref = sequence[position - 1:position - 1 + length]
                        for alt in ("A", "C", "G", "T", "AA", "AT", "CA", "GT"):
                            if ref == alt:
                                continue
                            key = reference.normalize(make_row(1, position, 0, reference=ref, alternate=alt))
                            _, normalized_position, normalized_ref, normalized_alt = key
                            index = int(normalized_position) - 1
                            haplotype = sequence[:position - 1] + alt + sequence[position - 1 + length:]
                            self.assertEqual(haplotype, sequence[:index] + normalized_alt + sequence[index + len(normalized_ref):])
                            self.assertEqual(key, reference.normalize(make_row(1, normalized_position, 0, reference=normalized_ref, alternate=normalized_alt)))
                            if haplotype in canonical_by_haplotype:
                                self.assertEqual(key, canonical_by_haplotype[haplotype])
                            canonical_by_haplotype[haplotype] = key

    def test_reference_tampering_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = write_reference(root)
            (root / "reference.fa").write_bytes(b">chr1\nG\n")
            with self.assertRaisesRegex(ValueError, "hash"):
                ReferenceGenome(config)
