"""ClinGen eRepo canlı API'sini sayfa hash'leriyle yerel tarama anlık görüntüsüne al."""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

try:
    from scripts.build_baseline_features import sha256_file
except ModuleNotFoundError:
    from build_baseline_features import sha256_file


API_BASE = "https://erepo.clinicalgenome.org/evrepo/api/summary/classifications"
TERMS_URL = "https://www.clinicalgenome.org/docs/terms-of-use/"
DEFAULT_OUTPUT = Path("data/external/clingen_erepo_snapshot")


def fetch_bytes(url: str) -> tuple[bytes, dict]:
    request = Request(url, headers={"User-Agent": "2026SYZ-research-cohort-screen/1.0", "Accept": "application/json,text/html"})
    with urlopen(request, timeout=45) as response:
        return response.read(), {
            "url": response.url, "http_status": response.status,
            "content_type": response.headers.get("Content-Type"),
            "last_modified": response.headers.get("Last-Modified"),
            "etag": response.headers.get("ETag"),
        }


def read_api(data: bytes) -> dict:
    document = json.loads(data)
    if document.get("status", {}).get("code") != 200 or "data" not in document:
        raise ValueError("eRepo API başarılı yanıt vermedi.")
    return document


def download_snapshot(output_dir: Path, *, page_size: int = 500) -> dict:
    if not 1 <= page_size <= 500:
        raise ValueError("Sayfa boyutu 1–500 arasında olmalıdır.")
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"Anlık görüntü klasörü boş olmalıdır: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema_version": 1, "status": "downloading",
        "started_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": API_BASE, "snapshot_kind": "retrieval_snapshot_not_official_release",
        "page_size": page_size, "sort": "uuid", "order": "asc", "pages": [],
        "limitation": "Canlı API'nin işlemsel snapshot garantisi yoktur; sayfa hash'leri, benzersizlik ve önce/sonra sayımları doğrulanır.",
    }
    manifest_path = output_dir / "manifest.json"
    def save():
        manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    save()
    try:
        raw_count, count_http = fetch_bytes(API_BASE + "/count")
        expected_total = int(read_api(raw_count)["data"]["total"])
        if expected_total <= 0:
            raise ValueError("eRepo toplam kayıt sayısı pozitif olmalıdır.")
        count_path = output_dir / "count_before.json"
        count_path.write_bytes(raw_count)
        manifest["count_before"] = {"path": count_path.name, "sha256": sha256_file(count_path), **count_http}
        manifest["expected_total"] = expected_total
        terms_data, terms_http = fetch_bytes(TERMS_URL)
        terms_path = output_dir / "terms_of_use.html"
        terms_path.write_bytes(terms_data)
        manifest["terms_of_use"] = {"path": terms_path.name, "sha256": sha256_file(terms_path), **terms_http}
        seen: set[str] = set()
        for page_number in range(1, (expected_total + page_size - 1) // page_size + 1):
            url = API_BASE + "?" + urlencode({"pgSize": page_size, "pg": page_number, "sort": "uuid", "order": "asc"})
            raw, http = fetch_bytes(url)
            payload = read_api(raw)
            rows = payload["data"]
            expected_rows = min(page_size, expected_total - len(seen))
            if not isinstance(rows, list) or len(rows) != expected_rows:
                raise ValueError(f"Sayfa {page_number} kayıt sayısı beklenenden farklı.")
            for row in rows:
                identifier = row.get("uuid")
                if not isinstance(identifier, str) or not identifier or identifier in seen:
                    raise ValueError("eRepo snapshot UUID eksik/tekrarlı; canlı kaynak değişmiş olabilir.")
                seen.add(identifier)
            page_path = output_dir / f"page_{page_number:04d}.json.gz"
            with gzip.open(page_path, "wb") as handle:
                handle.write(raw)
            manifest["pages"].append({
                "page": page_number, "path": page_path.name, "sha256": sha256_file(page_path),
                "rows": len(rows), "server_metadata": payload.get("metadata"), **http,
            })
            save()
            print(f"eRepo indirme: {len(seen):,}/{expected_total:,} yorum ({page_number}. sayfa)", flush=True)
        raw_after, after_http = fetch_bytes(API_BASE + "/count")
        after_total = int(read_api(raw_after)["data"]["total"])
        count_after_path = output_dir / "count_after.json"
        count_after_path.write_bytes(raw_after)
        manifest["count_after"] = {"path": count_after_path.name, "sha256": sha256_file(count_after_path), **after_http}
        if after_total != expected_total or len(seen) != expected_total:
            raise ValueError("İndirme öncesi/sonrası eRepo sayımları uyuşmuyor.")
        manifest.update({"status": "complete", "downloaded_unique_classifications": len(seen), "finished_at_utc": datetime.now(timezone.utc).isoformat()})
        save()
        return manifest
    except BaseException as error:
        manifest.update({"status": "failed", "error": f"{type(error).__name__}: {error}"})
        save()
        raise


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description="eRepo yorumlarını hash'li yerel snapshot olarak indirir.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--page-size", type=int, default=500)
    args = parser.parse_args()
    download_snapshot(args.output_dir, page_size=args.page_size)


if __name__ == "__main__":
    main()
