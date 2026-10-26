import gzip
from contextlib import closing
import io
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch

from scripts import local_vep_pipeline as pipeline
from scripts.compare_gtf_context_pilot import context_features, enrich_consequence, export_completed
from scripts.prepare_gtf_resources import download_ranges, parse_attributes, transcript_database


class GtfContextTest(unittest.TestCase):
    def source(self, path, partial=False):
        attrs = 'gene_id "G1"; gene_name "G"; transcript_id "ENST1"; transcript_version "2"; transcript_biotype "protein_coding"; tag "MANE_Select"; tag "Ensembl_canonical";'
        if partial:
            attrs += ' tag "cds_end_NF";'
        lines = [("transcript", 100, 209), ("exon", 100, 109), ("CDS", 104, 109), ("exon", 200, 209), ("CDS", 200, 205)]
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write("#!genome-version GRCh38\n")
            for feature, start, end in lines:
                handle.write(f"1\thavana\t{feature}\t{start}\t{end}\t.\t+\t0\t{attrs}\n")

    def test_metadata_preserves_tags_version_strand_exons_and_cds_length(self):
        with tempfile.TemporaryDirectory() as directory:
            source, db = Path(directory) / "source.gz", Path(directory) / "db.sqlite"
            self.source(source)
            self.assertEqual(transcript_database(source, db), 1)
            with closing(sqlite3.connect(db)) as connection:
                model = json.loads(connection.execute("SELECT model FROM transcripts").fetchone()[0])
            self.assertEqual(model["tags"], ["MANE_Select", "Ensembl_canonical"])
            self.assertEqual(model["version"], 2)
            self.assertEqual(model["Translation"]["length"], 4)
            self.assertEqual(len(model["Exon"]), 2)

    def test_partial_cds_does_not_invent_a_complete_protein_length(self):
        with tempfile.TemporaryDirectory() as directory:
            source, db = Path(directory) / "source.gz", Path(directory) / "db.sqlite"
            self.source(source, partial=True)
            transcript_database(source, db)
            with closing(sqlite3.connect(db)) as connection:
                model = json.loads(connection.execute("SELECT model FROM transcripts").fetchone()[0])
            self.assertEqual(model["Translation"], {})

    def test_gtf_features_keep_domains_missing_and_check_transcript_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            source, db = Path(directory) / "source.gz", Path(directory) / "db.sqlite"
            self.source(source)
            transcript_database(source, db)
            consequence = dict(transcript_id="ENST1", strand=1, variant_allele="G", protein_start=1,
                               exon="1/2", consequence_terms=["missense_variant"])
            record = dict(allele_string="A/G", seq_region_name="1", transcript_consequences=[consequence])
            with closing(sqlite3.connect(db)) as connection:
                features = context_features(record, "G", 109, connection)
                self.assertEqual(features["selected_transcript"], "ENST1.2")
                self.assertEqual(features["biotype"], "protein_coding")
                self.assertEqual(features["nearest_splice_boundary_bp"], 0)
                self.assertEqual(features["protein_relative_position"], 0.25)
                self.assertIsNone(features["pfam_annotation_count"])
                self.assertIsNone(features["interpro_annotation_count"])
                with self.assertRaisesRegex(ValueError, "version mismatch"):
                    enrich_consequence({**consequence, "transcript_id": "ENST1.3"}, connection)
                with self.assertRaisesRegex(ValueError, "strand mismatch"):
                    enrich_consequence({**consequence, "strand": -1}, connection)

    def test_range_download_is_bounded_verified_and_reusable_offline(self):
        payload = b"abcdefghijklmnopqrstuvwxyz"

        class Response(io.BytesIO):
            status = 206

            def __init__(self, start, end):
                super().__init__(payload[start:end + 1])
                self.headers = {"Content-Range": f"bytes {start}-{end}/{len(payload)}", "ETag": "version1"}

        def open_range(request, timeout):
            start, end = map(int, request.headers["Range"].split("=")[1].split("-"))
            return Response(start, end)

        with tempfile.TemporaryDirectory() as directory, patch("builtins.print"):
            path = Path(directory) / "source"
            with patch("scripts.prepare_gtf_resources.urlopen", side_effect=open_range):
                download_ranges("https://provider/source", path, len(payload), "version1", block_size=7)
            self.assertEqual(path.read_bytes(), payload)
            with patch("scripts.prepare_gtf_resources.urlopen") as network:
                download_ranges("https://provider/source", path, len(payload), "version1", block_size=7)
                network.assert_not_called()
            next((path.parent / "download_ranges").glob("*.part")).write_bytes(b"x" * 7)
            with self.assertRaises(ValueError), patch("scripts.prepare_gtf_resources.urlopen"):
                download_ranges("https://provider/source", path, len(payload), "version1", block_size=7)

    def test_range_download_rejects_a_full_response_or_different_etag(self):
        class Response(io.BytesIO):
            status = 200
            headers = {"ETag": "other", "Content-Range": "bytes 0-6/7"}

        with tempfile.TemporaryDirectory() as directory, patch("scripts.prepare_gtf_resources.urlopen", return_value=Response(b"ignored")):
            with self.assertRaisesRegex(ValueError, "ignored range or changed"):
                download_ranges("https://provider/source", Path(directory) / "source", 7, "version1", block_size=7)

    def test_gtf_command_omits_cache_and_unavailable_fields_and_reuses_inputs(self):
        root = Path(__file__).resolve().parents[1]
        config = pipeline.load_config(root / "config/local_vep_gtf.json")
        cache_config = pipeline.load_config(root / "config/local_vep_pipeline.json")
        self.assertEqual(pipeline.preparation_hash(config), pipeline.preparation_hash(cache_config))
        command = pipeline.build_command(config, root / "input", {"path": "chunk.vcf"}, root / "output",
                                         dict(fasta=str(root / "ref/fa"), gtf_dir=str(root / "gtf")))
        self.assertIn("--gtf", command)
        self.assertNotIn("--offline", command)
        self.assertNotIn("--database", command)
        self.assertIn("--network=none", command)
        self.assertEqual(command[command.index("--cache") + 1], "0")
        self.assertNotIn("--domains", command)
        self.assertNotIn("--mane", command)

    def test_export_keeps_partial_training_coverage_explicit_and_excludes_domains(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.source(root / "source.gz")
            transcript_database(root / "source.gz", root / "transcripts.sqlite")
            prepared = root / "prepared"
            prepared.mkdir()
            line = "1\t109\t1\tA\tG\t.\t.\t.\n"
            (prepared / "chunk.vcf").write_bytes((pipeline.VCF_HEADER + line).encode())
            (prepared / "variant_map.tsv").write_bytes(b"VariationID\tGeneSymbol\n1\tG\n2\tG\n")
            chunk = dict(path="chunk.vcf", sha256=pipeline.sha256_file(prepared / "chunk.vcf"))
            batch = root / "batch/chunk"
            batch.mkdir(parents=True)
            record = dict(id="1", input=line.strip(), seq_region_name="1", start=109, end=109,
                          assembly_name="GRCh38", allele_string="A/G", transcript_consequences=[dict(
                              transcript_id="ENST1.2", strand=1, variant_allele="G", protein_start=1,
                              exon="1/2", biotype="protein_coding", consequence_terms=["missense_variant"])])
            (batch / "output.jsonl").write_text(json.dumps(record) + "\n")
            pipeline.write_json(batch / "complete.json", dict(input_sha256=chunk["sha256"], resource_sha256="lock",
                                                              output="output.jsonl", output_sha256=pipeline.sha256_file(batch / "output.jsonl")))
            config = dict(output_dir="batch", protected=[])
            config_path = root / "config.json"
            pipeline.write_json(config_path, config)
            bundle = dict(chunks=[chunk], counts=dict(accepted=2))
            with patch.object(pipeline, "ROOT", root), patch.object(pipeline, "verify_bundle", return_value=(prepared, bundle)), \
                    patch.object(pipeline, "verify_resources", return_value=({"gtf_dir": str(root)}, "lock")), patch("builtins.print"):
                manifest = export_completed(config, config_path)
            self.assertEqual(manifest["rows"], 1)
            self.assertFalse(manifest["complete"])
            header = (root / "reports/generated/local_vep_gtf_features/train_context_features.tsv").read_text().splitlines()[0]
            self.assertNotIn("pfam", header)
            self.assertNotIn("target", header)
            self.assertNotIn("selected_transcript", header)


if __name__ == "__main__":
    unittest.main()
