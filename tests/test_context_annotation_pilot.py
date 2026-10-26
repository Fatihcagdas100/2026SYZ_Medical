import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from scripts.pilot_context_annotations import (
    CachedRest, FEATURE_COLUMNS, RAW_COLUMNS, annotate, extract_features, run, select_training_rows,
    select_transcript, splice_geometry, validate_vep,
)


def transcript(strand=1):
    return {"id": "ENST1", "version": 2, "assembly_name": "GRCh38", "seq_region_name": "1",
            "strand": strand, "Exon": [{"start": 100, "end": 109}, {"start": 200, "end": 209}],
            "Translation": {"length": 100}}


class ContextAnnotationPilotTest(unittest.TestCase):
    def test_sampling_is_train_only_label_free_and_order_independent(self):
        raw = pd.DataFrame([dict(VariationID=str(i), Chromosome="1", ChromosomeAccession="NC_000001.11",
                                 ReferenceAlleleVCF="A", AlternateAlleleVCF="G") for i in range(8)])
        features = pd.DataFrame([dict(VariationID=str(i), protein_effect="missense",
                                      splice_offset_category="not_intronic") for i in range(8)])
        split = pd.DataFrame({"VariationID": [str(i) for i in range(8)],
                              "split": ["train"] * 5 + ["validation", "test", "test"]})
        first = select_training_rows(raw, features, split, ["missense"], 2, 2027)
        second = select_training_rows(raw.iloc[::-1], features.iloc[::-1], split, ["missense"], 2, 2027)
        self.assertEqual(first, second)
        self.assertTrue(all(int(r["VariationID"]) < 5 for r in first))
        with self.assertRaisesRegex(ValueError, "Not enough"):
            select_training_rows(raw, features, split, ["missense"], 6, 2027)

    def test_transcript_choice_does_not_follow_severity_or_clinical_fields(self):
        candidates = [dict(transcript_id="Z", gene_symbol="G", canonical=1, consequence_terms=["stop_gained"]),
                      dict(transcript_id="B", gene_symbol="G", mane_select="NM_1", clinical_significance="benign"),
                      dict(transcript_id="A", gene_symbol="OTHER", mane_select="NM_2")]
        self.assertEqual(select_transcript(candidates, "G")["transcript_id"], "B")
        self.assertEqual(select_transcript(candidates[::-1], "G")["transcript_id"], "B")

    def test_junction_distances_are_strand_aware_and_exclude_transcript_ends(self):
        plus = splice_geometry(109, transcript())
        self.assertEqual(plus["donor_distance_bp"], 0)
        self.assertIsNone(plus["acceptor_distance_bp"])
        minus = splice_geometry(109, transcript(-1))
        self.assertIsNone(minus["donor_distance_bp"])
        self.assertEqual(minus["acceptor_distance_bp"], 0)
        self.assertEqual(splice_geometry(110, transcript())["donor_distance_bp"], 1)
        self.assertEqual(splice_geometry(110, transcript(-1))["acceptor_distance_bp"], 1)
        single = transcript()
        single["Exon"] = single["Exon"][:1]
        self.assertIsNone(splice_geometry(105, single)["nearest_splice_boundary_bp"])
        self.assertEqual(splice_geometry(99, transcript())["location"], "outside_transcript")

    def test_feature_allowlist_excludes_ids_clinical_labels_and_scores(self):
        consequence = dict(transcript_id="ENST1.2", gene_symbol="G", exon="1/2", biotype="protein_coding",
                           consequence_terms=["synonymous_variant", "splice_region_variant"], protein_start=10,
                           clin_sig=["pathogenic"], sift_score=0, polyphen_score=1,
                           domains=[{"db": "Pfam", "name": "PF1"}] * 2 + [{"db": "ENSP_mappings", "name": "structure"}])
        features = extract_features(consequence, transcript(), 109)
        self.assertEqual(set(features), set(FEATURE_COLUMNS))
        self.assertEqual(features["pfam_annotation_count"], 1)
        self.assertEqual(features["interpro_annotation_count"], 0)
        self.assertEqual(features["protein_relative_position"], 0.1)
        self.assertFalse(set(features) & {"target", "clin_sig", "gene_symbol", "transcript_id", "sift_score"})

    def test_coordinate_and_allele_mismatch_fail_closed(self):
        reply = dict(assembly_name="GRCh38", seq_region_name="1", start=109, end=109, allele_string="A/G")
        validate_vep(reply, ("1", "109", "A", "G"))
        for field, wrong in (("assembly_name", "GRCh37"), ("seq_region_name", "2"),
                             ("start", 110), ("end", 110), ("allele_string", "C/G")):
            with self.assertRaises(ValueError):
                validate_vep({**reply, field: wrong}, ("1", "109", "A", "G"))

    def test_offline_cache_replays_and_detects_tampering_without_network(self):
        with tempfile.TemporaryDirectory() as directory, patch("scripts.pilot_context_annotations.urlopen") as network:
            client = CachedRest(directory)
            with self.assertRaises(FileNotFoundError):
                client.get("/info/data")
            url = "https://rest.ensembl.org/info/data?content-type=application%2Fjson"
            key = hashlib.sha256((url + "#").encode()).hexdigest()
            body = Path(directory) / f"{key}.json"
            body.write_text('{"releases": [116]}', encoding="utf-8")
            receipt = Path(directory) / f"{key}.receipt.json"
            receipt.write_text(json.dumps({"url": url, "sha256": hashlib.sha256(body.read_bytes()).hexdigest(),
                                            "retrieved_at_utc": "2026-10-26T00:00:00Z"}), encoding="utf-8")
            self.assertEqual(client.get("/info/data"), {"releases": [116]})
            body.write_text('{"releases": [117]}', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Corrupt cache"):
                client.get("/info/data")
            network.assert_not_called()

    def test_annotation_checks_transcript_version_and_keeps_selection_in_audit(self):
        class Reference:
            def normalize(self, row):
                return "1", "109", "A", "G"

        class Client:
            version = 2

            def get(self, endpoint, params):
                if endpoint.startswith("/vep"):
                    return [dict(assembly_name="GRCh38", seq_region_name="1", start=109, end=109,
                                 allele_string="A/G", transcript_consequences=[dict(
                                     transcript_id="ENST1.2", gene_symbol="G", variant_allele="G", strand=1,
                                     biotype="protein_coding", consequence_terms=["synonymous_variant"])])]
                return {**transcript(), "version": self.version}

        row = dict(VariationID="1", Name="name", GeneSymbol="G", role="train_feasibility", stratum="synonymous")
        client = Client()
        audit = annotate(row, Reference(), client)
        self.assertEqual(audit["selected_transcript"], "ENST1.2")
        self.assertEqual(audit["nearest_splice_boundary_bp"], 0)
        client.version = 3
        with self.assertRaisesRegex(ValueError, "Transcript version"):
            annotate(row, Reference(), client)

    def test_pipeline_keeps_diagnostics_separate_and_preserves_protected_files(self):
        class Reference:
            metadata = {"fixture": True}

            def __init__(self, path):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *args):
                pass

            def normalize(self, row):
                return "1", "109", "A", "G"

        class Client:
            used = []

            def __init__(self, directory, fetch):
                pass

            def get(self, endpoint, params=None, tag=""):
                if endpoint == "/info/data":
                    return {"releases": [116]}
                if endpoint.startswith("/info/"):
                    return {"release": "fixture"}
                if endpoint.startswith("/vep"):
                    return [dict(assembly_name="GRCh38", seq_region_name="1", start=109, end=109,
                                 allele_string="A/G", transcript_consequences=[dict(
                                     transcript_id="ENST1.2", gene_symbol="G", variant_allele="G", strand=1,
                                     mane_select="NM_1", biotype="protein_coding",
                                     consequence_terms=["synonymous_variant"])])]
                return transcript()

        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            template = dict(zip(RAW_COLUMNS, ["1", "NM_1:c.1A>G", "G", "GRCh38", "NC_000001.11", "1", "109", "A", "G"]))
            pd.DataFrame([template, {**template, "VariationID": "2"}]).assign(target=[0, 1]).to_csv(root / "variants.tsv", sep="\t", index=False)
            pd.DataFrame([dict(VariationID=str(i), protein_effect="missense", splice_offset_category="not_intronic")
                          for i in [1, 2]]).to_csv(root / "features.tsv", sep="\t", index=False)
            pd.DataFrame({"VariationID": ["1", "2"], "split": ["train", "test"], "target": [0, 1]}).to_csv(root / "split.tsv", sep="\t", index=False)
            pd.DataFrame([{**template, "VariationID": "3", "target": 1}]).to_csv(root / "diagnostics.tsv", sep="\t", index=False)
            protected = root / "protected.bin"
            protected.write_bytes(b"frozen artifact fixture")

            def pin(name):
                return dict(path=name, sha256=hashlib.sha256((root / name).read_bytes()).hexdigest())

            config = dict(inputs={k: pin(k + ".tsv") for k in ["variants", "features", "split", "diagnostics"]},
                          protected=[pin("protected.bin")], diagnostic_ids=["3"], strata=["missense"],
                          per_stratum=1, seed=2027, ensembl_release=116, reference_config="unused.json",
                          cache="cache", output="output")
            path = root / "config.json"
            path.write_text(json.dumps(config), encoding="utf-8")
            with patch("scripts.pilot_context_annotations.ROOT", root), \
                    patch("scripts.pilot_context_annotations.ReferenceGenome", Reference), \
                    patch("scripts.pilot_context_annotations.CachedRest", Client), \
                    patch("scripts.pilot_context_annotations.urlopen") as network, \
                    patch("builtins.print"):
                manifest = run(path)
            network.assert_not_called()
            self.assertEqual(protected.read_bytes(), b"frozen artifact fixture")
            train = pd.read_csv(root / "output/train_candidate_features.tsv", sep="\t", dtype=str)
            diagnostics = pd.read_csv(root / "output/seen_diagnostic_features.tsv", sep="\t", dtype=str)
            self.assertEqual(train.VariationID.tolist(), ["1"])
            self.assertEqual(diagnostics.VariationID.tolist(), ["3"])
            self.assertEqual(train.columns.tolist(), ["VariationID"] + FEATURE_COLUMNS)
            self.assertFalse(manifest["policy"]["labels_loaded"])
            original_get = Client.get

            def upgraded_service(self, endpoint, params=None, tag=""):
                if tag.startswith("live_before_"):
                    return {"releases": [117]}
                return original_get(self, endpoint, params, tag)

            with patch("scripts.pilot_context_annotations.ROOT", root), \
                    patch("scripts.pilot_context_annotations.CachedRest", Client), \
                    patch.object(Client, "get", upgraded_service), patch("builtins.print"):
                with self.assertRaisesRegex(ValueError, "Live release differs"):
                    run(path, fetch=True)


if __name__ == "__main__":
    unittest.main()
