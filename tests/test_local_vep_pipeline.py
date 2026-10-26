import json
from pathlib import Path
from subprocess import CompletedProcess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import pandas as pd

from scripts import local_vep_pipeline as pipeline


class LocalVepPipelineTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.addCleanup(patch.stopall)
        patch.object(pipeline, "ROOT", self.root).start()
        patch("builtins.print").start()
        self.config = json.loads((Path(__file__).resolve().parents[1] / "config/local_vep_pipeline.json").read_text(encoding="utf-8"))
        fasta = self.root / "tiny.fa"
        fasta.write_bytes((">chr1\n" + "A" * 200 + "\n").encode("ascii"))
        Path(str(fasta) + ".fai").write_bytes(b"chr1\t200\t6\t200\t201\n")
        (self.root / "aliases.tsv").write_text("chr1\tNC_000001.11\n", encoding="ascii")
        reference = dict(assembly="GRCh38", fasta="tiny.fa", fasta_sha256=pipeline.sha256_file(fasta),
                         fai_sha256=pipeline.sha256_file(Path(str(fasta) + ".fai")),
                         primary_accession_to_chromosome={"NC_000001.11": "1"}, chrom_alias_file="aliases.tsv",
                         chrom_alias_sha256=pipeline.sha256_file(self.root / "aliases.tsv"))
        pipeline.write_json(self.root / "reference.json", reference)
        template = dict(VariationID="1", Name="NM_1:c.1A>G", GeneSymbol="G", Assembly="GRCh38",
                        ChromosomeAccession="NC_000001.11", Chromosome="1", PositionVCF="1",
                        ReferenceAlleleVCF="A", AlternateAlleleVCF="G", target="1")
        rows = [template, {**template, "VariationID": "2", "PositionVCF": "2"},
                {**template, "VariationID": "3", "PositionVCF": "3"},
                {**template, "VariationID": "4", "AlternateAlleleVCF": "AA"},
                {**template, "VariationID": "5", "Chromosome": "MT"},
                {**template, "VariationID": "6", "ChromosomeAccession": "NW_1"},
                {**template, "VariationID": "7", "PositionVCF": "99"}]
        pd.DataFrame(rows).to_csv(self.root / "variants.tsv", sep="\t", index=False)
        pd.DataFrame({"VariationID": [str(i) for i in range(1, 8)],
                      "split": ["train", "test", "validation", "train", "train", "train", "train"],
                      "target": [1] * 7}).to_csv(self.root / "split.tsv", sep="\t", index=False)
        (self.root / "frozen.bin").write_bytes(b"frozen model fixture")
        self.config.update(chunk_size=1, prepared_dir="prepared", output_dir="outputs", cache_dir="cache",
                           resource_lock="resources.lock.json", readiness_report="preflight.json")
        self.config["inputs"] = {"variants": self.pin("variants.tsv"), "split": self.pin("split.tsv"),
                                 "reference_config": self.pin("reference.json")}
        self.config["protected"] = [self.pin("frozen.bin")]
        cache = self.root / "cache/homo_sapiens/116_GRCh38"
        cache.mkdir(parents=True)
        (cache / "info.txt").write_text("version\t116\nassembly\tGRCh38\n", encoding="ascii")
        self.config_path = self.root / "config.json"
        pipeline.write_json(self.config_path, self.config)

    def pin(self, path):
        return dict(path=path, sha256=pipeline.sha256_file(self.root / path))

    def prepare(self):
        return pipeline.prepare(self.config, self.config_path)

    def fake_docker(self, arguments, timeout=60):
        if "--show_cache_info" in arguments:
            return CompletedProcess(arguments, 0, "cache version 116 GRCh38\n", "")
        output_mount = next(a for a in arguments if a.startswith("type=bind,source=") and ",target=/output" in a)
        output = Path(output_mount.removeprefix("type=bind,source=").split(",target=")[0])
        vcf_name = Path(arguments[arguments.index("--input_file") + 1]).name
        vcf = self.root / "prepared" / vcf_name
        records = []
        for line in vcf.read_text().splitlines():
            if line.startswith("#"):
                continue
            chrom, pos, identifier, ref, alt, *_ = line.split("\t")
            records.append(dict(input=line, id=identifier, assembly_name="GRCh38", seq_region_name=chrom,
                                start=int(pos), end=int(pos), allele_string=ref + "/" + alt,
                                transcript_consequences=[dict(transcript_id="ENST1.1", consequence_terms=["missense_variant"])]))
        (output / "annotations.jsonl").write_text("".join(json.dumps(r) + "\n" for r in records), encoding="utf-8")
        return CompletedProcess(arguments, 0, "", "")

    def lock(self):
        with patch.object(pipeline, "docker", side_effect=self.fake_docker):
            return pipeline.lock_resources(self.config, self.config_path)

    def test_prepare_accounts_for_all_train_rows_and_excludes_labels_and_other_roles(self):
        with patch.object(pd, "read_csv", wraps=pd.read_csv) as reader:
            manifest = self.prepare()
        self.assertEqual(manifest["counts"], dict(train_rows=5, accepted=2, excluded=3, chunks=2))
        self.assertTrue(all("target" not in call.kwargs["usecols"] for call in reader.call_args_list))
        mapping = pd.read_csv(self.root / "prepared/variant_map.tsv", sep="\t", dtype=str)
        self.assertEqual(mapping.VariationID.tolist(), ["1", "7"])
        self.assertNotIn("target", mapping)
        self.assertEqual(sum(manifest["exclusion_reasons"].values()), 3)
        self.assertEqual((self.root / "frozen.bin").read_bytes(), b"frozen model fixture")

    def test_prepare_is_idempotent_and_refuses_tampered_chunks(self):
        first = self.prepare()
        self.assertEqual(self.prepare(), first)
        chunk = self.root / "prepared" / first["chunks"][0]["path"]
        chunk.write_text(chunk.read_text() + "corruption\n")
        with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
            self.prepare()

    def test_runtime_cache_path_can_change_without_rebuilding_train_inputs(self):
        first = self.prepare()
        self.config["cache_dir"] = "larger_disk_cache"
        pipeline.write_json(self.config_path, self.config)
        self.assertEqual(self.prepare(), first)
        self.config["chunk_size"] = 2
        pipeline.write_json(self.config_path, self.config)
        with self.assertRaisesRegex(ValueError, "Prepared input recipe differs"):
            self.prepare()

    def test_reference_mismatch_aborts_instead_of_becoming_an_exclusion(self):
        frame = pd.read_csv(self.root / "variants.tsv", sep="\t", dtype=str)
        frame.loc[frame.VariationID.eq("1"), "ReferenceAlleleVCF"] = "C"
        frame.to_csv(self.root / "variants.tsv", sep="\t", index=False)
        self.config["inputs"]["variants"] = self.pin("variants.tsv")
        pipeline.write_json(self.config_path, self.config)
        with self.assertRaisesRegex(ValueError, "reference_allele_mismatch"):
            self.prepare()
        self.assertFalse((self.root / "prepared/manifest.json").exists())

    def test_docker_command_pins_image_and_disables_network_with_readonly_inputs(self):
        command = pipeline.build_command(self.config, self.root / "prepared", {"path": "chunk_1.vcf"},
                                         self.root / "out with spaces", dict(fasta=str(self.root / "tiny.fa"), cache_dir=str(self.root / "cache")))
        self.assertIn("--network=none", command)
        self.assertIn("--pull=never", command)
        self.assertIn(self.config["image"], command)
        self.assertIn("--offline", command)
        self.assertIn("--symbol", command)
        self.assertNotIn("--check_existing", command)
        mounts = [arg for arg in command if arg.startswith("type=bind")]
        self.assertEqual(sum(m.endswith(",readonly") for m in mounts), 3)
        self.assertTrue(any("out with spaces" in m for m in mounts))

    def test_output_requires_exact_ids_coordinates_alleles_assembly_and_annotations(self):
        self.prepare()
        vcf = self.root / "prepared/chunk_000001.vcf"
        line = [line for line in vcf.read_text().splitlines() if not line.startswith("#")][0]
        record = dict(input=line, id="1", assembly_name="GRCh38", seq_region_name="1", start=1, end=1,
                      allele_string="A/G", intergenic_consequences=[{"consequence_terms": ["intergenic_variant"]}])
        output = self.root / "out.jsonl"
        output.write_text(json.dumps(record) + "\n")
        self.assertEqual(pipeline.validate_output(vcf, output), 1)
        wrong_records = [{**record, "start": 2}, {**record, "allele_string": "C/G"},
                         {**record, "assembly_name": "GRCh37"}, {**record, "id": "999"},
                         {**record, "intergenic_consequences": []}, {**record, "colocated_variants": [{"clin_sig": ["pathogenic"]}]}]
        for bad in wrong_records:
            output.write_text(json.dumps(bad) + "\n")
            with self.assertRaises(ValueError):
                pipeline.validate_output(vcf, output)
        output.write_text(2 * (json.dumps(record) + "\n"))
        with self.assertRaisesRegex(ValueError, "duplicate"):
            pipeline.validate_output(vcf, output)
        output.write_text("")
        with self.assertRaisesRegex(ValueError, "skipped or omitted"):
            pipeline.validate_output(vcf, output)

    def test_cache_inventory_detects_additions_and_modifications(self):
        self.lock()
        pipeline.verify_resources(self.config, self.config_path)
        (self.root / "cache/homo_sapiens/116_GRCh38/info.txt").write_text("changed")
        with self.assertRaisesRegex(ValueError, "Cache inventory changed"):
            pipeline.verify_resources(self.config, self.config_path)

    def test_interrupted_chunk_retries_then_completed_chunk_is_reused(self):
        self.prepare()
        self.lock()
        with patch.object(pipeline, "docker", return_value=CompletedProcess([], 1, "partial", "failed")):
            with self.assertRaisesRegex(RuntimeError, "incomplete attempt retained"):
                pipeline.run_batches(self.config, self.config_path)
        self.assertFalse((self.root / "outputs/chunk_000001/complete.json").exists())
        with patch.object(pipeline, "docker", side_effect=self.fake_docker) as engine:
            pipeline.run_batches(self.config, self.config_path)
            self.assertEqual(engine.call_count, 1)
            pipeline.run_batches(self.config, self.config_path)
            self.assertEqual(engine.call_count, 2)  # first chunk reused, second newly run
            pipeline.run_batches(self.config, self.config_path)
            self.assertEqual(engine.call_count, 2)  # all chunks reused
        self.assertEqual(len(list((self.root / "outputs/chunk_000001").glob("attempt_*"))), 2)
        receipt = json.loads((self.root / "outputs/chunk_000001/complete.json").read_text())
        (self.root / "outputs/chunk_000001" / receipt["output"]).write_text("tampered")
        with self.assertRaisesRegex(ValueError, "SHA256 mismatch"):
            pipeline.run_batches(self.config, self.config_path)
        self.assertEqual((self.root / "frozen.bin").read_bytes(), b"frozen model fixture")

    def test_changed_execution_flags_cannot_reuse_completed_chunks(self):
        self.prepare()
        self.lock()
        with patch.object(pipeline, "docker", side_effect=self.fake_docker):
            pipeline.run_batches(self.config, self.config_path)
        with patch.object(pipeline, "VEP_FLAGS", pipeline.VEP_FLAGS + ["--everything"]):
            with self.assertRaisesRegex(ValueError, "Completed chunk"):
                pipeline.run_batches(self.config, self.config_path)

    def test_preflight_reports_missing_cache_and_disk_budget_without_downloading(self):
        self.config["cache_dir"] = "absent_cache"
        with patch.object(pipeline.shutil, "disk_usage", return_value=SimpleNamespace(free=100)), \
                patch.object(pipeline, "docker", return_value=CompletedProcess([], 0, "ensembl-vep : 116.1", "")), \
                patch.object(pipeline, "urlopen") as network:
            report = pipeline.preflight(self.config)
        network.assert_not_called()
        self.assertFalse(report["ready"])
        self.assertIn("release_116_grch38_cache_missing", report["blockers"])
        self.assertIn("disk_below_planning_budget", report["blockers"])
        self.assertIn("resources_not_locked", report["blockers"])


if __name__ == "__main__":
    unittest.main()
