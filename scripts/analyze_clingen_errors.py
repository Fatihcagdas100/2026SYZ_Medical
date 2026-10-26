"""Kaydedilmiş ClinGen tahminlerinde açıklayıcı hata analizi; model yüklemez."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from scripts.train_logistic_baseline import sha256_file
    from scripts.screen_clingen_erepo import load_snapshot
except ModuleNotFoundError:
    from train_logistic_baseline import sha256_file
    from screen_clingen_erepo import load_snapshot


DEFAULT_LOCK = Path("config/clingen_error_analysis_inputs.json")
DEFAULT_REPORT = Path("reports/generated/clingen_error_analysis")
GROUP_FEATURES = (
    "protein_effect", "has_protein_annotation", "annotation_source", "GeneSymbol",
    "missense_substitution", "missense_residue_class_transition", "missense_charge_transition",
    "variant_type", "coding_region", "hgvs_sequence_type", "fit_pair_support",
    "has_unknown_fit_category",
)


def verified_file(info: dict) -> Path:
    path = Path(info["path"])
    digest = info.get("sha256")
    if not isinstance(digest, str) or len(digest) != 64 or sha256_file(path) != digest:
        raise ValueError(f"Hata analizi girdi hash uyuşmazlığı: {path}")
    return path


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def wilson_interval(events: int, denominator: int) -> tuple[float | None, float | None]:
    if not denominator:
        return None, None
    z = 1.959963984540054
    fraction = events / denominator
    divisor = 1 + z * z / denominator
    center = (fraction + z * z / (2 * denominator)) / divisor
    radius = z * math.sqrt(fraction * (1 - fraction) / denominator + z * z / (4 * denominator**2)) / divisor
    return max(0.0, center - radius), min(1.0, center + radius)


def join_saved_tables(predictions: pd.DataFrame, features: pd.DataFrame, cohort: pd.DataFrame, *, threshold: float) -> pd.DataFrame:
    for name, frame in (("tahmin", predictions), ("özellik", features), ("kohort", cohort)):
        if frame["VariationID"].isna().any() or frame["VariationID"].duplicated().any():
            raise ValueError(f"{name} kimlikleri eksik veya tekrarlı.")
        if set(frame["VariationID"]) != set(predictions["VariationID"]):
            raise ValueError(f"{name} ile tahmin ID kümesi uyuşmuyor.")
    data = predictions.merge(features.rename(columns={"target": "feature_target"}), on="VariationID", validate="one_to_one")
    metadata = cohort[["VariationID", "Name", "GeneSymbol", "annotation_source", "target", "clingen_panels_json", "original_clinvar_labels_json"]].rename(columns={"target": "cohort_target"})
    data = data.merge(metadata, on="VariationID", validate="one_to_one")
    if not data["target"].eq(data["feature_target"]).all() or not data["target"].eq(data["cohort_target"]).all():
        raise ValueError("Kaydedilmiş tahmin, özellik ve kohort hedefleri uyuşmuyor.")
    probability = data["probability_pathogenic"]
    if not np.isfinite(probability).all() or not probability.between(0, 1).all():
        raise ValueError("Kaydedilmiş olasılık geçersiz.")
    if not (probability.ge(threshold).astype(int) == data["prediction"]).all():
        raise ValueError("Kaydedilmiş kararlar donmuş eşikle uyuşmuyor.")
    data["FN"] = (data["target"].eq(1) & data["prediction"].eq(0)).astype(int)
    data["FP"] = (data["target"].eq(0) & data["prediction"].eq(1)).astype(int)
    data["error_type"] = np.select([data["FN"].eq(1), data["FP"].eq(1)], ["false_negative", "false_positive"], default="correct")
    data["annotation_source"] = data["annotation_source"].str.split(";").str[0]
    return data


def training_coverage(data: pd.DataFrame, candidate: dict) -> tuple[pd.DataFrame, dict]:
    """Yalnızca fit=train+validation özellik desteğini say; hedef sütunu okuma."""
    fit_path = verified_file(candidate["inputs"]["features"])
    split_path = verified_file(candidate["inputs"]["gene_disjoint_split"])
    split = pd.read_csv(split_path, sep="\t", usecols=["VariationID", "split"], dtype="string")
    if split["VariationID"].duplicated().any() or split.isna().any().any():
        raise ValueError("Fit ayrım kimlikleri geçersiz.")
    fit_ids = set(split.loc[split["split"].isin(["train", "validation"]), "VariationID"])
    categorical = candidate["model"]["categorical_columns"]
    columns = list(dict.fromkeys(["VariationID", "protein_effect", "missense_substitution", *categorical]))
    categories = {column: set() for column in categorical}
    pair_counts: Counter[str] = Counter()
    found = set()
    for chunk in pd.read_csv(fit_path, sep="\t", usecols=columns, dtype="string", chunksize=100_000):
        selected = chunk.loc[chunk["VariationID"].isin(fit_ids)]
        ids = set(selected["VariationID"])
        if found.intersection(ids) or selected["VariationID"].duplicated().any():
            raise ValueError("Fit özellik ID'si tekrarlandı.")
        found.update(ids)
        for column in categorical:
            categories[column].update(selected[column].dropna())
        pair_counts.update(selected.loc[selected["protein_effect"].eq("missense"), "missense_substitution"].dropna())
    if found != fit_ids:
        raise ValueError("Fit özellik ID'leri ayrımla uyuşmuyor.")
    data = data.copy()
    missense = data["protein_effect"].eq("missense")
    data["fit_pair_count"] = data["missense_substitution"].map(pair_counts).fillna(0).astype(int)
    data["fit_pair_support"] = pd.cut(data["fit_pair_count"], bins=[-1, 0, 19, 99, float("inf")], labels=["unseen", "1-19", "20-99", "100+"]).astype("string")
    data.loc[~missense, "fit_pair_support"] = "not_missense"
    unknowns = pd.DataFrame({column: ~data[column].astype("string").isin(values) & data[column].notna() for column, values in categories.items()})
    data["has_unknown_fit_category"] = unknowns.any(axis=1).astype(int)
    data["unknown_fit_columns_json"] = [json.dumps(list(unknowns.columns[mask])) for mask in unknowns.to_numpy()]
    return data, {
        "fit_rows": len(found), "test_targets_read": False,
        "missense_unseen_pair_rows": int((missense & data["fit_pair_count"].eq(0)).sum()),
        "missense_seen_pair_rows": int((missense & data["fit_pair_count"].gt(0)).sum()),
        "rows_with_unknown_categorical_value": int(data["has_unknown_fit_category"].sum()),
        "unknown_category_counts": {column: int(unknowns[column].sum()) for column in categorical},
    }


def subgroup_tables(data: pd.DataFrame, min_class_count: int) -> pd.DataFrame:
    rows = []
    for scope, subset in (("all_variants", data), ("missense", data.loc[data["protein_effect"].eq("missense")])):
        total_fn, total_fp = int(subset["FN"].sum()), int(subset["FP"].sum())
        for feature in GROUP_FEATURES:
            for value, group in subset.groupby(feature, dropna=False, sort=False):
                positives, negatives = int(group["target"].eq(1).sum()), int(group["target"].eq(0).sum())
                fn, fp = int(group["FN"].sum()), int(group["FP"].sum())
                fn_low, fn_high = wilson_interval(fn, positives)
                fp_low, fp_high = wilson_interval(fp, negatives)
                rows.append({
                    "scope": scope, "feature": feature, "value": str(value), "rows": len(group),
                    "positives": positives, "negatives": negatives, "false_negatives": fn, "false_positives": fp,
                    "false_negative_rate": fn / positives if positives else None,
                    "false_positive_rate": fp / negatives if negatives else None,
                    "fnr_wilson_lower_95": fn_low, "fnr_wilson_upper_95": fn_high,
                    "fpr_wilson_lower_95": fp_low, "fpr_wilson_upper_95": fp_high,
                    "fnr_support_at_least_minimum": positives >= min_class_count,
                    "fpr_support_at_least_minimum": negatives >= min_class_count,
                    "share_of_scope_false_negatives": fn / total_fn if total_fn else None,
                    "share_of_scope_false_positives": fp / total_fp if total_fp else None,
                })
    return pd.DataFrame(rows)


def feature_collisions(data: pd.DataFrame, feature_columns: list[str]) -> pd.DataFrame:
    """Model girdileri birebir aynı, etiketleri farklı missense gruplarını bul."""
    rows = []
    missense = data.loc[data["protein_effect"].eq("missense")]
    for _, group in missense.groupby(feature_columns, dropna=False, sort=False):
        if group["target"].nunique() < 2:
            continue
        if group["probability_pathogenic"].max() - group["probability_pathogenic"].min() > 1e-9:
            raise ValueError("Aynı model girdilerinde kaydedilmiş olasılıklar farklı; girdi/tahmin ilişkisini inceleyin.")
        rows.append({
            "rows": len(group), "positives": int(group["target"].eq(1).sum()),
            "negatives": int(group["target"].eq(0).sum()),
            "false_negatives": int(group["FN"].sum()), "false_positives": int(group["FP"].sum()),
            "probability_min": float(group["probability_pathogenic"].min()),
            "probability_max": float(group["probability_pathogenic"].max()),
            "missense_substitution": group["missense_substitution"].iloc[0],
            "genes_json": json.dumps(sorted(set(group["GeneSymbol"]))),
            "variation_ids_json": json.dumps(sorted(set(group["VariationID"]))),
        })
    return pd.DataFrame(rows, columns=["rows", "positives", "negatives", "false_negatives", "false_positives", "probability_min", "probability_max", "missense_substitution", "genes_json", "variation_ids_json"])


def saved_curation_evidence(snapshot_manifest: Path, errors: pd.DataFrame) -> pd.DataFrame:
    """Hata örneklerinin uzman özetlerini aynı tarihli snapshot'tan ilişkilendir."""
    _, rows = load_snapshot(snapshot_manifest.parent)
    targets = errors.set_index("VariationID")["target"].to_dict()
    evidence = []
    for row in rows:
        variation_id = str(row.get("cvId"))
        if variation_id not in targets or row.get("retracted") is not False:
            continue
        expected_label = "Pathogenic" if targets[variation_id] else "Benign"
        if row["classification"] != expected_label:
            raise ValueError("Hata örneğinin snapshot uzman etiketi farklı.")
        evidence.append({
            "VariationID": variation_id, "uuid": row["uuid"], "classification": row["classification"],
            "gene": row.get("gene", ""), "condition": row.get("condition", ""),
            "expert_panel": row.get("ep", ""), "doc_version": row.get("docVersion", ""),
            "summary_desc": row.get("summaryDesc", ""),
        })
    if {row["VariationID"] for row in evidence} != set(targets):
        raise ValueError("Bütün hata kimlikleri için aynı snapshot'ta uzman yorumu bulunamadı.")
    return pd.DataFrame(evidence, columns=["VariationID", "uuid", "classification", "gene", "condition", "expert_panel", "doc_version", "summary_desc"])


def analyze_errors(input_lock: Path, report_dir: Path) -> dict:
    lock = load_json(input_lock)
    if lock.get("analysis_type") != "post_evaluation_descriptive_error_analysis":
        raise ValueError("Hata analizi girdi kilidi geçersiz.")
    evaluation_path = verified_file(lock["evaluation_manifest"])
    predictions_path = verified_file(lock["predictions"])
    evaluation = load_json(evaluation_path)
    if evaluation.get("evaluation_number") != 1:
        raise ValueError("Beklenen birincil değerlendirme kaydı bulunamadı.")
    if predictions_path.resolve() != Path(evaluation["outputs"]["predictions"]).resolve():
        raise ValueError("Tahmin yolu birincil değerlendirmeyle uyuşmuyor.")
    features_path = verified_file(evaluation["external_features"])
    admission = load_json(verified_file(evaluation["admission_manifest"]))
    cohort_path = verified_file(admission["output"])
    candidate_path = verified_file(evaluation["candidate_manifest"])
    candidate = load_json(candidate_path)
    protected_parents = {
        path.resolve().parent for path in (
            evaluation_path, predictions_path, features_path, cohort_path,
            candidate_path, Path(evaluation["admission_manifest"]["path"]), input_lock,
        )
    }
    if report_dir.resolve() in protected_parents:
        raise ValueError("Hata analizi çıktısı birincil girdi klasörüne yazılamaz.")
    threshold = evaluation["frozen_threshold"]
    if threshold != candidate["threshold"]:
        raise ValueError("Kaydedilmiş eşik adayla uyuşmuyor.")
    read_args = {"sep": "\t", "dtype": {"VariationID": "string"}}
    data = join_saved_tables(pd.read_csv(predictions_path, **read_args), pd.read_csv(features_path, **read_args), pd.read_csv(cohort_path, **read_args), threshold=threshold)
    metrics = evaluation["metrics"]
    if len(data) != metrics["row_count"] or int(data["FN"].sum()) != metrics["false_negatives"] or int(data["FP"].sum()) != metrics["false_positives"]:
        raise ValueError("Kayıtlı hata sayıları birincil sonuçla uyuşmuyor.")
    data, coverage = training_coverage(data, candidate)
    min_class_count = int(lock["min_class_count"])
    if min_class_count < 1:
        raise ValueError("Asgari sınıf sayısı pozitif olmalıdır.")
    lower, upper = map(float, lock["probability_regions"])
    if not 0 <= lower < threshold < upper <= 1:
        raise ValueError("Tanımlı olasılık bölgeleri donmuş eşikle uyumsuz.")
    groups = subgroup_tables(data, min_class_count)
    collisions = feature_collisions(data, candidate["model"]["feature_columns"])
    errors = data.loc[data["FN"].eq(1) | data["FP"].eq(1)].copy()
    errors["wrong_class_confidence"] = np.where(errors["FN"].eq(1), 1 - errors["probability_pathogenic"], errors["probability_pathogenic"])
    errors = errors.sort_values(["wrong_class_confidence", "VariationID"], ascending=[False, True])
    evidence = saved_curation_evidence(verified_file(lock["curation_snapshot_manifest"]), errors)
    report_dir.mkdir(parents=True, exist_ok=True)
    groups_path, collisions_path, errors_path = (report_dir / name for name in ("subgroup_error_rates.csv", "missense_feature_collisions.tsv", "saved_error_examples.tsv"))
    groups.to_csv(groups_path, index=False, float_format="%.10f")
    collisions.to_csv(collisions_path, sep="\t", index=False, float_format="%.10f")
    example_columns = ["VariationID", "target", "prediction", "probability_pathogenic", "error_type", "GeneSymbol", "Name", "protein_effect", "missense_substitution", "annotation_source", "fit_pair_count", "unknown_fit_columns_json", "clingen_panels_json", "original_clinvar_labels_json"]
    errors[example_columns].to_csv(errors_path, sep="\t", index=False, float_format="%.10f")
    evidence_path = report_dir / "error_curation_evidence.tsv"
    evidence.to_csv(evidence_path, sep="\t", index=False)
    counts = {
        "cohort_rows": len(data), "false_negatives": int(data["FN"].sum()), "false_positives": int(data["FP"].sum()),
        "missense_false_negatives": int(data.loc[data["protein_effect"].eq("missense"), "FN"].sum()),
        "low_probability_false_negatives": int((data["FN"].eq(1) & data["probability_pathogenic"].lt(lower)).sum()),
        "high_probability_false_positives": int((data["FP"].eq(1) & data["probability_pathogenic"].gt(upper)).sum()),
        "conflicting_label_missense_feature_groups": len(collisions),
        "rows_in_conflicting_label_feature_groups": int(collisions["rows"].sum()),
        "false_negatives_in_conflicting_label_feature_groups": int(collisions["false_negatives"].sum()),
        "false_positives_in_conflicting_label_feature_groups": int(collisions["false_positives"].sum()),
    }
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "analysis_kind": lock["analysis_type"],
        "inputs": {"input_lock": {"path": str(input_lock.resolve()), "sha256": sha256_file(input_lock)},
                   "evaluation_manifest": lock["evaluation_manifest"], "predictions": lock["predictions"],
                   "features": evaluation["external_features"], "cohort": admission["output"],
                   "curation_snapshot_manifest": lock["curation_snapshot_manifest"],
                   "fit_features": candidate["inputs"]["features"], "fit_split": candidate["inputs"]["gene_disjoint_split"]},
        "frozen_threshold": threshold, "counts": counts, "training_coverage": coverage,
        "data_policy": {"saved_predictions_only": True, "model_loaded": False, "model_refit": False, "threshold_search": False, "internal_test_targets_read": False, "interpretation": "posthoc_descriptive_not_causal"},
        "subgroup_policy": {"min_class_count": min_class_count, "gene_groups": "exact_GeneSymbol_field_without_explosion", "interval_method": "Wilson_row_binomial_not_gene_clustered", "multiple_comparison_correction": False},
        "outputs": {"subgroups": str(groups_path.resolve()), "collisions": str(collisions_path.resolve()), "error_examples": str(errors_path.resolve()), "curation_evidence": str(evidence_path.resolve())},
    }
    report_path = report_dir / "clingen_error_analysis.md"
    lines = ["# ClinGen kaydedilmiş tahminlerinde hata analizi", "", "Bu, sonuç sonrası açıklayıcı analizdir. Model yeniden yüklenmedi; eşik taranmadı.", "", f"Sayımlar: {counts}", "", f"Fit desteği: {coverage}", ""]
    for scope, feature in (("all_variants", "protein_effect"), ("all_variants", "annotation_source"), ("all_variants", "GeneSymbol"), ("missense", "missense_substitution"), ("missense", "annotation_source")):
        subset = groups.loc[groups["scope"].eq(scope) & groups["feature"].eq(feature)]
        subset = subset.sort_values(["false_negatives", "positives"], ascending=[False, False])
        if feature in {"GeneSymbol", "missense_substitution"}:
            subset = subset.loc[subset["fnr_support_at_least_minimum"]].head(12)
        lines.extend([f"## {scope} / {feature}", "", "| Grup | n | Pathogenic | Benign | FN | FP | FN/Pathogenic |", "| --- | ---: | ---: | ---: | ---: | ---: | ---: |"])
        for row in subset.to_dict("records"):
            rate = f"{row['false_negative_rate']:.4f}" if row["positives"] else "—"
            lines.append(f"| {row['value']} | {row['rows']} | {row['positives']} | {row['negatives']} | {row['false_negatives']} | {row['false_positives']} | {rate} |")
        lines.append("")
    lines.extend(["Gruplar gen, etiket ve anotasyon dağılımları açısından farklıdır; oranlar nedensel etki değildir. Küçük gruplar ve çoklu inceleme keşifsel yorum gerektirir.", "", "Birebir aynı model girdileri farklı etiketler taşıyabilir; bu bulgu mevcut gösterimin bağlam ayrımını sınırladığını gösterir. Yeni bir özellik veya eşik, bu görülmüş kohortta tekrar doğrulanmış sayılamaz."])
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest["outputs"]["report"] = str(report_path.resolve())
    # Salt okunur birincil girdilerin analiz boyunca değişmediğini doğrula.
    verified_file(lock["evaluation_manifest"])
    verified_file(lock["predictions"])
    (report_dir / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"counts": counts, "training_coverage": coverage}, ensure_ascii=False, indent=2))
    print(f"Hata analizi: {report_path.resolve()}")
    return manifest


def main():
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    parser = argparse.ArgumentParser(description="Dondurulmuş ClinGen sonucunun kayıtlı hatalarını çözümler.")
    parser.add_argument("--input-lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT)
    args = parser.parse_args()
    analyze_errors(args.input_lock, args.report_dir)


if __name__ == "__main__":
    main()
