import contextlib
import gzip
import io
import json
import tempfile
import unittest
from pathlib import Path
from urllib.parse import parse_qs, urlparse
from unittest.mock import patch

from scripts.download_clingen_erepo import TERMS_URL, download_snapshot
from scripts.screen_clingen_erepo import screen_snapshot
from tests.external_fixtures import make_row, write_rows


def curation(uuid, cvid, classification, **extra):
    return {
        "uuid": str(uuid), "cvId": cvid, "classification": classification,
        "retracted": False, "gene": "GENE", "caId": f"CA{cvid}",
        "ep": "Test panel", "condition": "Test condition", "docVersion": "1.0.0",
        "preferredVarTitle": "NM_1(GENE):c.1A>G (p.Arg1Gly)", **extra,
    }


def source_stub(rows, *, final_count=None):
    count_calls = 0
    def fetch(url):
        nonlocal count_calls
        http = {"url": url, "http_status": 200}
        if url == TERMS_URL:
            return b"<html>Synthetic terms fixture</html>", http
        if url.endswith("/count"):
            count_calls += 1
            total = len(rows) if count_calls == 1 or final_count is None else final_count
            return json.dumps({"status": {"code": 200}, "data": {"total": total}}).encode(), http
        query = parse_qs(urlparse(url).query)
        page, size = int(query["pg"][0]), int(query["pgSize"][0])
        return json.dumps({"status": {"code": 200}, "data": rows[(page - 1) * size:page * size]}).encode(), http
    return fetch


class ClinGenErepoScreenTest(unittest.TestCase):
    def test_pagination_conflicts_missing_ids_and_multiple_curations_preserved(self):
        rows = [
            curation(1, "10", "Pathogenic"), curation(2, "20", "Benign"),
            curation(3, "30", "Pathogenic"), curation(4, "30", "Benign"),
            curation(5, None, "Benign"), curation(6, "40", "Likely Pathogenic"),
            curation(7, "50", "Pathogenic"), curation(8, "50", "Pathogenic", condition="Other condition"),
            curation(9, "60", "Benign", retracted=True),
            curation(10, "70", "Pathogenic"), curation(11, "70", "Likely Pathogenic"),
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            snapshot = root / "snapshot"
            development = root / "development.tsv.gz"
            write_rows(development, [make_row(10, 1, 0)])
            with patch("scripts.download_clingen_erepo.fetch_bytes", side_effect=source_stub(rows)), contextlib.redirect_stdout(io.StringIO()):
                manifest = download_snapshot(snapshot, page_size=2)
            self.assertEqual(manifest["status"], "complete")
            self.assertEqual(len(manifest["pages"]), 6)
            with contextlib.redirect_stdout(io.StringIO()):
                result = screen_snapshot(snapshot, development, root / "report")
            counts = result["counts"]
            self.assertEqual(counts["id_unmatched_candidate_unique_variants"], 2)
            self.assertEqual(counts["id_unmatched_benign_variants"], 1)
            self.assertEqual(counts["id_unmatched_pathogenic_variants"], 1)
            self.assertEqual(counts["previous_id_overlap_unique_variants"], 1)
            self.assertEqual(counts["discordant_unique_clinvar_ids"], 2)
            self.assertEqual(counts["exact_records_without_single_clinvar_id"], 1)
            candidates = (root / "report" / "id_unmatched_curation_candidates.tsv").read_text(encoding="utf-8")
            self.assertIn("Other condition", candidates)
            self.assertEqual(len(candidates.splitlines()), 4)  # header + 3 curations for 2 variants
            self.assertFalse(result["policies"]["model_loaded"])
            with self.assertRaises(FileExistsError):
                download_snapshot(snapshot)
            with (snapshot / "page_0001.json.gz").open("ab") as handle:
                handle.write(b"corruption")
            with self.assertRaisesRegex(ValueError, "hash"):
                screen_snapshot(snapshot, development, root / "changed_report")

    def test_changing_source_count_marks_snapshot_failed(self):
        rows = [curation(1, "10", "Benign"), curation(2, "20", "Pathogenic")]
        with tempfile.TemporaryDirectory() as directory:
            snapshot = Path(directory) / "snapshot"
            with patch("scripts.download_clingen_erepo.fetch_bytes", side_effect=source_stub(rows, final_count=3)), contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError, "sayımları"):
                    download_snapshot(snapshot, page_size=2)
            manifest = json.loads((snapshot / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "failed")
            with self.assertRaisesRegex(ValueError, "tamamlanmış"):
                screen_snapshot(snapshot, Path(directory) / "unused.tsv.gz", Path(directory) / "report")


if __name__ == "__main__":
    unittest.main()
