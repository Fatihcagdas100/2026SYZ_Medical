import contextlib
import csv
import gzip
import io
import json
import tempfile
import unittest
from pathlib import Path

from scripts.build_baseline_features import sha256_file
from scripts.build_initial_binary_dataset import OUTPUT_SOURCE_COLUMNS
from scripts.map_clingen_candidates import map_candidates
from tests.external_fixtures import make_row, write_rows


class MapClinGenCandidatesTest(unittest.TestCase):
    def test_label_provenance_and_ambiguous_or_missing_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            candidates = root / "candidates.tsv"
            rows = [{"ClinVarVariationID": str(variation_id), "classification": label,
                     "uuid": str(variation_id), "condition": "condition", "expert_panel": "panel", "doc_version": "2.0"}
                    for variation_id, label in ((1, "Pathogenic"), (2, "Benign"), (3, "Pathogenic"))]
            with candidates.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            screen = root / "screen.json"
            screen.write_text(json.dumps({"decision": "screening_only_not_admitted_for_external_evaluation",
                "outputs": {"curation_candidates": {"path": str(candidates), "sha256": sha256_file(candidates)}},
                "counts": {"id_unmatched_candidate_unique_variants": 3}}), encoding="utf-8")
            raw = root / "raw.tsv.gz"
            source_rows = []
            for variation_id, position in ((1, 5), (2, 6), (2, 7)):
                values = {key: "-" for key in OUTPUT_SOURCE_COLUMNS}
                values.update({key: value for key, value in make_row(variation_id, position, 0).items() if key in values})
                values.update({"OriginSimple": "germline", "ClinicalSignificance": "Uncertain significance"})
                source_rows.append(values)
            write_rows(raw, source_rows)
            output = root / "output.tsv.gz"
            with contextlib.redirect_stdout(io.StringIO()):
                result = map_candidates(screen, raw, output, root / "report")
            self.assertEqual(result["counts"]["output_mapped_variants"], 1)
            self.assertEqual(result["counts"]["ambiguous_multiple_source_representations"], 1)
            self.assertEqual(result["counts"]["not_found_in_pinned_grch38_germline_source"], 1)
            with gzip.open(output, "rt", encoding="utf-8") as handle:
                mapped = list(csv.DictReader(handle, delimiter="\t"))[0]
            self.assertEqual(mapped["target_label"], "Pathogenic")
            self.assertEqual(mapped["target"], "1")
            self.assertIn("Uncertain significance", mapped["original_clinvar_labels_json"])
            self.assertFalse(result["provenance"]["source_independence_verified"])

