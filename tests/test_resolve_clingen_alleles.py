import json
import tempfile
import unittest
from pathlib import Path

from scripts.build_baseline_features import sha256_file
from scripts.resolve_clingen_alleles import registry_vcf
from scripts.variant_normalization import InvalidVariant, ReferenceGenome
from tests.external_fixtures import make_row, write_reference


class ResolveClinGenAllelesTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        root = Path(self.temporary.name)
        config_path = write_reference(root, "CAAAAACGTACGTACGT")
        config = json.loads(config_path.read_text())
        alias_file = root / "chromAlias.txt"
        alias_file.write_text("chr1\t1\tNC_000001.11\n")
        config.update({"primary_accession_to_chromosome": {"NC_000001.11": "1"},
                       "chrom_alias_file": alias_file.name, "chrom_alias_sha256": sha256_file(alias_file)})
        config_path.write_text(json.dumps(config))
        self.reference = ReferenceGenome(config_path)
        self.addCleanup(self.reference.handle.close)

    def document(self, start, ref, alt):
        return {"externalRecords": {"ClinVarVariations": [{"variationId": 123}]},
                "genomicAlleles": [{"referenceGenome": "GRCh38", "chromosome": "1",
                    "hgvs": ["NC_000001.11:g.example"],
                    "coordinates": [{"start": start, "end": start + len(ref), "referenceAllele": ref, "allele": alt}]}]}

    def test_registry_snv_is_converted_from_zero_based_coordinates(self):
        row = registry_vcf(self.document(6, "C", "T"), "123", self.reference)
        self.assertEqual((row["PositionVCF"], row["ReferenceAlleleVCF"], row["AlternateAlleleVCF"]), ("7", "C", "T"))

    def test_deletion_and_insertion_are_anchored_and_reference_checked(self):
        deletion = registry_vcf(self.document(3, "AA", ""), "123", self.reference)
        self.assertEqual((deletion["PositionVCF"], deletion["ReferenceAlleleVCF"], deletion["AlternateAlleleVCF"]), ("3", "AAA", "A"))
        self.assertEqual(self.reference.normalize(deletion), ("1", "1", "CAA", "C"))
        insertion = registry_vcf(self.document(3, "", "A"), "123", self.reference)
        self.assertEqual(self.reference.normalize(insertion), ("1", "1", "C", "CA"))
        boundary = registry_vcf(self.document(0, "C", ""), "123", self.reference)
        self.assertEqual((boundary["PositionVCF"], boundary["ReferenceAlleleVCF"], boundary["AlternateAlleleVCF"]), ("1", "CA", "A"))

    def test_wrong_identity_and_reference_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "identity"):
            registry_vcf(self.document(6, "C", "T"), "456", self.reference)
        with self.assertRaisesRegex(ValueError, "reference_allele_mismatch"):
            registry_vcf(self.document(6, "A", "T"), "123", self.reference)

    def test_patch_accession_is_not_treated_as_primary_chromosome(self):
        row = make_row(1, 7, 0, reference="C", alternate="T")
        row["ChromosomeAccession"] = "NW_009646201.1"
        with self.assertRaisesRegex(InvalidVariant, "nonprimary"):
            self.reference.normalize(row)
        row["ChromosomeAccession"] = "NC_000001.11"
        row["Chromosome"] = "2"
        with self.assertRaisesRegex(InvalidVariant, "accession_chromosome_mismatch"):
            self.reference.normalize(row)
