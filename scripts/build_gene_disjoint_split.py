"""Gen bileşenleri çakışmayan deterministik train/validation/test ayrımı üret."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

import numpy as np
from sklearn.model_selection import StratifiedGroupKFold


DEFAULT_INPUT = Path(
    "data/processed/clinvar_grch38_germline_exact_binary.tsv.gz"
)
DEFAULT_OUTPUT = Path("data/processed/baseline_gene_disjoint_split.tsv")
DEFAULT_REPORT_DIR = Path("reports/generated/gene_disjoint_split")
DEFAULT_SEED = 2027
N_FOLDS = 20
FOLD_TO_SPLIT = {
    **{fold: "train" for fold in range(14)},
    **{fold: "validation" for fold in range(14, 17)},
    **{fold: "test" for fold in range(17, 20)},
}
EXPECTED_RATIOS = {"train": 0.70, "validation": 0.15, "test": 0.15}
MISSING_GENE_MARKERS = {"", "-", "na"}
GENE_SEPARATOR = re.compile(r"[;,]")


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def find(self, item: str) -> str:
        self.parent.setdefault(item, item)
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            parent = self.parent[item]
            self.parent[item] = root
            item = parent
        return root

    def union(self, first: str, second: str) -> str:
        first_root = self.find(first)
        second_root = self.find(second)
        if first_root == second_root:
            return first_root
        canonical, other = sorted((first_root, second_root))
        self.parent[other] = canonical
        return canonical


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_size):
            digest.update(block)
    return digest.hexdigest()


def parse_gene_tokens(raw_value: str) -> tuple[str, ...]:
    normalized = raw_value.strip()
    if normalized.casefold() in MISSING_GENE_MARKERS:
        return ("__UNKNOWN_GENE__",)
    if normalized.casefold().startswith("covers "):
        return (f"__DESCRIPTIVE__:{normalized.casefold()}",)
    tokens = tuple(
        sorted(
            {
                token.strip().upper()
                for token in GENE_SEPARATOR.split(normalized)
                if token.strip()
            }
        )
    )
    return tokens or ("__UNKNOWN_GENE__",)


def read_records(
    input_path: Path,
) -> tuple[list[tuple[str, str, tuple[str, ...]]], UnionFind, Counter[str]]:
    if not input_path.is_file():
        raise FileNotFoundError(f"Girdi dosyası bulunamadı: {input_path}")

    union_find = UnionFind()
    counters: Counter[str] = Counter()
    records: list[tuple[str, str, tuple[str, ...]]] = []
    seen_ids: set[str] = set()
    with gzip.open(input_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Girdi dosyasında başlık bulunamadı.")
        missing = sorted({"VariationID", "GeneSymbol", "target"}.difference(reader.fieldnames))
        if missing:
            raise ValueError("Gerekli sütunlar bulunamadı: " + ", ".join(missing))

        for row in reader:
            variation_id = row["VariationID"]
            target = row["target"]
            raw_gene = row["GeneSymbol"]
            if not variation_id:
                raise ValueError("Boş VariationID bulundu.")
            if variation_id in seen_ids:
                raise ValueError(f"Tekrarlanan VariationID bulundu: {variation_id}")
            if target not in {"0", "1"}:
                raise ValueError(f"Beklenmeyen target değeri: {target}")
            seen_ids.add(variation_id)

            tokens = parse_gene_tokens(raw_gene)
            for token in tokens:
                union_find.find(token)
            for token in tokens[1:]:
                union_find.union(tokens[0], token)
            if len(tokens) > 1:
                counters["multi_gene_rows"] += 1
            if tokens[0] == "__UNKNOWN_GENE__":
                counters["unknown_gene_rows"] += 1
            if tokens[0].startswith("__DESCRIPTIVE__:"):
                counters["descriptive_gene_rows"] += 1
            counters["rows"] += 1
            records.append((variation_id, target, tokens))
    return records, union_find, counters


def assign_gene_disjoint_splits(
    records: Sequence[tuple[str, str, tuple[str, ...]]],
    union_find: UnionFind,
    *,
    seed: int,
) -> tuple[list[str], list[str], list[int]]:
    if not records:
        raise ValueError("Ayrım için en az bir kayıt gereklidir.")
    targets = np.array([int(target) for _, target, _ in records], dtype=np.int8)
    gene_groups = [union_find.find(tokens[0]) for _, _, tokens in records]
    unique_groups = set(gene_groups)
    if len(unique_groups) < N_FOLDS:
        raise ValueError(
            f"En az {N_FOLDS} gen bileşeni gereklidir; bulunan: {len(unique_groups)}"
        )

    fold_assignments = np.full(len(records), -1, dtype=np.int16)
    splitter = StratifiedGroupKFold(
        n_splits=N_FOLDS,
        shuffle=True,
        random_state=seed,
    )
    dummy_features = np.zeros((len(records), 1), dtype=np.int8)
    for fold, (_, fold_indices) in enumerate(
        splitter.split(dummy_features, targets, groups=gene_groups)
    ):
        fold_assignments[fold_indices] = fold
    if (fold_assignments < 0).any():
        raise RuntimeError("Bazı kayıtlar hiçbir kata atanmadı.")

    split_assignments = [FOLD_TO_SPLIT[int(fold)] for fold in fold_assignments]
    group_splits: defaultdict[str, set[str]] = defaultdict(set)
    for gene_group, split_name in zip(gene_groups, split_assignments):
        group_splits[gene_group].add(split_name)
    overlapping = [group for group, splits in group_splits.items() if len(splits) > 1]
    if overlapping:
        raise RuntimeError(f"Bir gen bileşeni birden çok ayrıma düştü: {overlapping[0]}")
    return split_assignments, gene_groups, fold_assignments.tolist()


def build_gene_disjoint_split(
    input_path: Path,
    output_path: Path,
    report_dir: Path,
    *,
    seed: int = DEFAULT_SEED,
) -> dict[str, object]:
    records, union_find, counters = read_records(input_path)
    split_assignments, gene_groups, fold_assignments = assign_gene_disjoint_splits(
        records,
        union_find,
        seed=seed,
    )

    split_counts: Counter[str] = Counter(split_assignments)
    target_counts: Counter[tuple[str, str]] = Counter()
    group_counts: Counter[str] = Counter(gene_groups)
    groups_by_split: defaultdict[str, set[str]] = defaultdict(set)
    for (_, target, _), split_name, gene_group in zip(
        records, split_assignments, gene_groups
    ):
        target_counts[(split_name, target)] += 1
        groups_by_split[split_name].add(gene_group)

    split_names = ("train", "validation", "test")
    for index, first in enumerate(split_names):
        for second in split_names[index + 1 :]:
            if groups_by_split[first].intersection(groups_by_split[second]):
                raise RuntimeError(f"{first} ve {second} arasında gen grubu çakışması var.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    temporary_handle = tempfile.NamedTemporaryFile(
        prefix="gene_disjoint_split_",
        suffix=".tsv.tmp",
        dir=output_path.parent,
        delete=False,
    )
    temporary_path = Path(temporary_handle.name)
    temporary_handle.close()
    try:
        with temporary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["VariationID", "target", "split", "gene_group", "fold"])
            for (variation_id, target, _), split_name, gene_group, fold in zip(
                records,
                split_assignments,
                gene_groups,
                fold_assignments,
            ):
                writer.writerow([variation_id, target, split_name, gene_group, fold])
        temporary_path.replace(output_path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise

    with (report_dir / "split_counts.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "split",
                "row_count",
                "percentage",
                "benign_count",
                "pathogenic_count",
                "pathogenic_percentage",
                "gene_group_count",
            ]
        )
        for split_name in split_names:
            rows = split_counts[split_name]
            benign = target_counts[(split_name, "0")]
            pathogenic = target_counts[(split_name, "1")]
            writer.writerow(
                [
                    split_name,
                    rows,
                    f"{100 * rows / len(records):.6f}",
                    benign,
                    pathogenic,
                    f"{100 * pathogenic / rows:.6f}",
                    len(groups_by_split[split_name]),
                ]
            )

    with (report_dir / "largest_gene_groups.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["gene_group", "row_count", "split"])
        group_to_split = {
            group: split_name
            for split_name, groups in groups_by_split.items()
            for group in groups
        }
        for gene_group, row_count in group_counts.most_common(100):
            writer.writerow([gene_group, row_count, group_to_split[gene_group]])

    total_target_counts = Counter(target for _, target, _ in records)
    manifest: dict[str, object] = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "source": {
            "path": str(input_path.resolve()),
            "size_bytes": input_path.stat().st_size,
            "sha256": sha256_file(input_path),
            "row_count": len(records),
        },
        "output": {
            "path": str(output_path.resolve()),
            "size_bytes": output_path.stat().st_size,
            "sha256": sha256_file(output_path),
            "columns": ["VariationID", "target", "split", "gene_group", "fold"],
        },
        "method": {
            "name": "stratified_group_k_fold_with_connected_gene_components",
            "seed": seed,
            "n_folds": N_FOLDS,
            "fold_allocation": {"train": 14, "validation": 3, "test": 3},
            "expected_ratios": EXPECTED_RATIOS,
            "multi_gene_policy": (
                "Semicolon/comma separated gene symbols are unioned; any shared "
                "symbol keeps the full connected component in one split."
            ),
            "unknown_gene_policy": "All missing-like gene values share one conservative group.",
            "descriptive_gene_policy": (
                "Values beginning with 'covers ' are treated as whole pseudo-gene groups."
            ),
        },
        "counts": {
            **dict(counters),
            "unique_gene_tokens": len(union_find.parent),
            "unique_gene_components": len(group_counts),
            "target_0_benign": total_target_counts["0"],
            "target_1_pathogenic": total_target_counts["1"],
            "by_split": dict(split_counts),
            "by_split_and_target": {
                split_name: {
                    "0": target_counts[(split_name, "0")],
                    "1": target_counts[(split_name, "1")],
                }
                for split_name in split_names
            },
            "gene_groups_by_split": {
                split_name: len(groups_by_split[split_name]) for split_name in split_names
            },
        },
        "integrity": {
            "gene_group_overlap_count": 0,
            "unassigned_rows": 0,
            "all_rows_assigned_once": True,
        },
        "evaluation_scope": (
            "Measures generalization to gene components not observed during training; "
            "it is stricter than the variant-level split."
        ),
    }
    with (report_dir / "manifest.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2)
        handle.write("\n")

    print("Gen-grubu çakışmasız veri ayrımı oluşturuldu")
    print(f"Toplam satır: {len(records):,}")
    print(f"Gen bileşeni: {len(group_counts):,}")
    for split_name in split_names:
        print(
            f"{split_name}: {split_counts[split_name]:,} satır, "
            f"{len(groups_by_split[split_name]):,} gen grubu"
        )
    print("Gen grubu çakışması: 0")
    print(f"Çıktı: {output_path.resolve()}")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Bağlı gen bileşenleri train/validation/test arasında çakışmayan ayrım üretir."
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return parser.parse_args()


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    args = parse_args()
    build_gene_disjoint_split(
        args.input,
        args.output,
        args.report_dir,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
