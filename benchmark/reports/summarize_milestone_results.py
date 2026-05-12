"""summarize_milestone_results.py
================================
Step 11: Summarize milestone benchmark results by provider_id, label_mode,
and cell_state_key, without mixing consensus / sensitivity / legacy scGPT outputs.

Writes four CSVs to benchmark/reports/ (or --output-dir):

  core_summary.csv
      One row per (result_dir x provider_id x label_mode x metric).
      Columns: dataset_id, method, scenario, run_id, result_dir,
               metric_family, metric_name, metric_value,
               provider_id, label_mode, cell_state_key, reference_graph_path,
               result_class, formal_benchmark, analysis_role,
               is_primary, is_sensitivity, is_legacy, source_json

  lineage_summary.csv
      One row per run_dir that has a milestone-tagged lineage_metrics.json.
      Columns: dataset_id, method, scenario, run_id, result_dir,
               provider_id, label_mode, cell_state_key, reference_graph_path,
               n_reference_edges, auroc, auprc, jaccard_topk, status,
               result_class, formal_benchmark, analysis_role,
               is_primary, is_sensitivity, is_legacy, source_json

  embedding_summary.csv
      One row per (result_dir x label_mode) from embedding_milestone_eval/.
      Columns: dataset_id, method, scenario, run_id, result_dir,
               provider_id, label_mode, state_key, embedding_key, cluster_source,
               adjusted_rand_index, mean_normalized_entropy,
               weighted_mean_normalized_entropy, n_cells_evaluated, status,
               result_class, formal_benchmark, analysis_role,
               is_primary, is_sensitivity, is_legacy, source_json

  provider_agreement_summary.csv
      One row per (result_dir x provider pair).
      Columns: dataset_id, result_dir, pair,
               pairwise_adjusted_rand_index, pairwise_exact_match_fraction,
               n_cells_compared, source_json

Classification rules:
  primary     label_mode in ('official_silver', 'consensus') (not legacy)
  sensitivity label_mode in ('embedding_based', 'classifier_based')
  legacy      provider_id contains 'scgpt' OR
              cell_state_key == 'scgpt_pseudostate_provisional' OR
              label_mode == 'legacy_scgpt_pseudostate'
  unknown     anything else -- never classified as primary

Smoke directory absence:
  Does NOT require benchmark/results/smoke to exist. Missing dirs are skipped.

Usage
-----
  python -m benchmark.reports.summarize_milestone_results
  python -m benchmark.reports.summarize_milestone_results \\
      --results-dir benchmark/results \\
      --output-dir  benchmark/reports
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

EMBEDDING_MODES = (
    "official_silver",
    "consensus",
    "embedding_based",
    "classifier_based",
)

VALID_RESULT_CLASSES = {
    "official", "formal", "pilot", "pilot_backup",
    "reduced_validation", "smoke", "hpc_validation",
    "diagnostic", "external_validation", "archive",
    "marker_fm_transition_silver_formal",
}
FORMAL_RESULT_CLASSES = {"official", "formal", "marker_fm_transition_silver_formal"}

LINEAGE_METRIC_NAMES_MAP = {
    "auroc": "auroc",
    "auprc": "auprc",
    "jaccard_similarity_topk": "jaccard_topk",
    "single_step_recovery": "single_step_recovery",
    "multi_step_recovery": "multi_step_recovery",
}

EMBEDDING_METRIC_NAMES_MAP = {
    "adjusted_rand_index": "adjusted_rand_index",
    "mean_prediction_entropy": "mean_prediction_entropy",
    "weighted_prediction_entropy": "weighted_prediction_entropy",
    "mean_normalized_entropy": "mean_normalized_entropy",
    "weighted_mean_normalized_entropy": "weighted_mean_normalized_entropy",
    "n_cells_evaluated": "n_cells_evaluated",
}

CORE_FIELDS = [
    "dataset_id", "method", "scenario", "run_id", "result_dir",
    "metric_family", "metric_name", "metric_value",
    "provider_id", "label_mode", "cell_state_key", "reference_graph_path",
    "result_class", "formal_benchmark", "analysis_role",
    "is_primary", "is_sensitivity", "is_legacy", "source_json",
]

LINEAGE_FIELDS = [
    "dataset_id", "method", "scenario", "run_id", "result_dir",
    "provider_id", "label_mode", "cell_state_key", "reference_graph_path",
    "n_reference_edges", "auroc", "auprc", "jaccard_topk", "status",
    "result_class", "formal_benchmark", "analysis_role",
    "is_primary", "is_sensitivity", "is_legacy", "source_json",
]

EMBEDDING_FIELDS = [
    "dataset_id", "method", "scenario", "run_id", "result_dir",
    "provider_id", "label_mode", "state_key", "embedding_key", "cluster_source",
    "adjusted_rand_index", "mean_prediction_entropy", "weighted_prediction_entropy",
    "entropy_basis", "n_probability_classes", "prediction_entropy_note",
    "mean_normalized_entropy", "weighted_mean_normalized_entropy",
    "hard_label_entropy_note", "n_cells_evaluated", "status",
    "result_class", "formal_benchmark", "analysis_role",
    "is_primary", "is_sensitivity", "is_legacy", "source_json",
]

AGREEMENT_FIELDS = [
    "dataset_id", "result_dir", "pair",
    "pairwise_adjusted_rand_index", "pairwise_exact_match_fraction",
    "n_cells_compared", "source_json",
]


# ---------------------------------------------------------------------------
# Project root
# ---------------------------------------------------------------------------

def _project_root() -> Path:
    here = Path(__file__).resolve().parent
    for c in [here, *here.parents]:
        if (c / "benchmark").exists() and (c / "benchmark" / "evaluation").exists():
            return c
    return here.parent.parent


# ---------------------------------------------------------------------------
# JSON loading
# ---------------------------------------------------------------------------

def _load_json_safe(path: Path) -> Optional[Dict]:
    if not path.exists():
        return None
    try:
        raw = path.read_text(encoding="utf-8").rstrip("\x00﻿\n\r\t ")
        return json.loads(raw) if raw else None
    except Exception as exc:
        print(f"[summarize_milestone] Warning: cannot read {path}: {exc}", file=sys.stderr)
        return None


# ---------------------------------------------------------------------------
# Classification
# ---------------------------------------------------------------------------

def _make_analysis_flags(
    label_mode: Optional[str],
    provider_id: Optional[str],
    cell_state_key: Optional[str] = None,
) -> Dict[str, Any]:
    lm = (label_mode or "").strip().lower()
    pid = (provider_id or "").lower()
    csk = (cell_state_key or "").lower()

    # Legacy first
    if (lm == "legacy_scgpt_pseudostate"
            or "scgpt" in pid
            or csk == "scgpt_pseudostate_provisional"):
        return {"analysis_role": "legacy",
                "is_primary": False, "is_sensitivity": False, "is_legacy": True}

    # Primary
    if lm in ("official_silver", "consensus"):
        return {"analysis_role": "primary",
                "is_primary": True, "is_sensitivity": False, "is_legacy": False}

    # Sensitivity
    if lm in ("embedding_based", "classifier_based"):
        return {"analysis_role": "sensitivity",
                "is_primary": False, "is_sensitivity": True, "is_legacy": False}

    # Fallback via provider_id
    if "_marker_fm_transition_silver_" in pid or "_milestone_consensus_" in pid:
        return {"analysis_role": "primary",
                "is_primary": True, "is_sensitivity": False, "is_legacy": False}
    if "_milestone_embedding_" in pid or "_milestone_classifier_" in pid:
        return {"analysis_role": "sensitivity",
                "is_primary": False, "is_sensitivity": True, "is_legacy": False}

    return {"analysis_role": "unknown",
            "is_primary": False, "is_sensitivity": False, "is_legacy": False}


def _run_is_legacy_scgpt(run_meta: Optional[Dict]) -> bool:
    """Return True when the run itself used the legacy scGPT pseudostate system.

    Some older run directories may contain later diagnostic milestone files under
    embedding_milestone_eval/. Those files should not promote the whole run back
    into primary milestone summaries.
    """
    if not run_meta:
        return False
    provider_id = (
        run_meta.get("provider_id")
        or (run_meta.get("ground_truth") or {}).get("provider_id")
        or ""
    )
    label_mode = run_meta.get("label_mode") or ""
    cell_state_key = (
        run_meta.get("cell_state_key")
        or (run_meta.get("ground_truth") or {}).get("state_key")
        or ""
    )
    flags = _make_analysis_flags(str(label_mode), str(provider_id), str(cell_state_key))
    return bool(flags["is_legacy"])


def _infer_result_class(meta: Optional[Dict]) -> str:
    declared = (meta or {}).get("result_class", "")
    if declared in VALID_RESULT_CLASSES:
        return str(declared)
    sampling = (meta or {}).get("sampling") or {}
    if sampling.get("subsample_per_timepoint") not in (None, "", False):
        return "pilot"
    return "official"


def _is_formal_benchmark(result_class: str) -> bool:
    return result_class in FORMAL_RESULT_CLASSES


def _safe(v: Any) -> Any:
    return "" if v is None else v


def _norm_path_text(v: Any) -> str:
    return str(v or "").replace("\\", "/").strip()


def _run_source_h5ad(run_meta: Optional[Dict]) -> str:
    if not run_meta:
        return ""
    gt = run_meta.get("ground_truth") or {}
    return _norm_path_text(
        gt.get("source_h5ad")
        or run_meta.get("input_h5ad")
        or run_meta.get("h5ad_path")
        or (run_meta.get("dataset") or {}).get("h5ad_path")
    )


def _embedding_matches_run_source(emb: Dict, run_meta: Optional[Dict]) -> bool:
    """Skip stale embedding diagnostics generated against an older h5ad."""
    run_source = _run_source_h5ad(run_meta)
    emb_source = _norm_path_text(
        emb.get("input_h5ad")
        or emb.get("reference_h5ad")
        or (emb.get("metadata") or {}).get("input_h5ad")
    )
    if not run_source or not emb_source:
        return True
    return run_source == emb_source


# ---------------------------------------------------------------------------
# Dataset ID / scenario helpers
# ---------------------------------------------------------------------------

def _infer_dataset_id(run_meta: Optional[Dict], metric_doc: Optional[Dict], run_dir: Optional[Path] = None) -> str:
    for d in (metric_doc, run_meta):
        if not d:
            continue
        for key in ("dataset_id", "dataset", "benchmark_dataset"):
            v = d.get(key)
            if v:
                return str(v)
        provider_id = d.get("provider_id") or (d.get("ground_truth") or {}).get("provider_id")
        if provider_id:
            pid = str(provider_id).lower()
            if "gse230659" in pid:
                return "GSE230659"
            if "gse178325" in pid:
                return "GSE178325"
        source_h5ad = (d.get("ground_truth") or {}).get("source_h5ad") or d.get("input_h5ad")
        if source_h5ad:
            source = str(source_h5ad).lower()
            if "gse230659" in source:
                return "GSE230659"
            if "gse178325" in source:
                return "GSE178325"
    if run_dir:
        name = str(run_dir).lower()
        if "gse230659" in name:
            return "GSE230659"
        if "gse178325" in name:
            return "GSE178325"
    return ""


def _parse_scenario(run_dir_name: str, run_meta: Optional[Dict] = None) -> str:
    if run_meta:
        sc = run_meta.get("scenario")
        if sc:
            return str(sc)
    if run_dir_name.startswith("scenario_") and len(run_dir_name) > len("scenario_"):
        return run_dir_name[len("scenario_"):]
    return run_dir_name


# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------

def _is_milestone_run_dir(run_dir: Path) -> bool:
    if (run_dir / "projected_milestone_annotation_metadata.json").exists():
        return True
    emb_dir = run_dir / "embedding_milestone_eval"
    if emb_dir.is_dir():
        for mode in EMBEDDING_MODES:
            if (emb_dir / f"embedding_metrics_{mode}.json").exists():
                return True
    lm_path = run_dir / "lineage_metrics.json"
    if lm_path.exists():
        d = _load_json_safe(lm_path)
        if d and d.get("label_mode"):
            return True
    return False


def _discover_milestone_runs(results_root: Path) -> List[Tuple[str, Path]]:
    runs: List[Tuple[str, Path]] = []
    if not results_root.exists():
        print(f"[summarize_milestone] Note: results root does not exist: {results_root}",
              file=sys.stderr)
        return runs

    # Smoke dir may be absent after cleanup -- not an error
    smoke_dir = results_root / "smoke"
    if not smoke_dir.exists():
        print("[summarize_milestone] Note: benchmark/results/smoke is absent "
              "(outputs may have been cleaned). Proceeding without smoke results.",
              file=sys.stderr)

    for method_dir in sorted(p for p in results_root.iterdir() if p.is_dir()):
        method = method_dir.name
        try:
            children = sorted(p for p in method_dir.iterdir() if p.is_dir())
        except PermissionError as exc:
            print(f"[summarize_milestone] Warning: cannot list {method_dir}: {exc}",
                  file=sys.stderr)
            continue
        for run_dir in children:
            try:
                if _is_milestone_run_dir(run_dir):
                    runs.append((method, run_dir))
            except Exception as exc:
                print(f"[summarize_milestone] Warning: error checking {run_dir}: {exc}",
                      file=sys.stderr)

    return runs


# ---------------------------------------------------------------------------
# Base row builder
# ---------------------------------------------------------------------------

def _base_row(method: str, run_dir: Path, run_meta: Optional[Dict],
              metric_doc: Optional[Dict] = None) -> Dict:
    result_class = _infer_result_class(run_meta)
    scenario = _parse_scenario(run_dir.name, run_meta)
    dataset_id = _infer_dataset_id(run_meta, metric_doc, run_dir)
    return {
        "dataset_id": dataset_id,
        "method": method,
        "scenario": scenario,
        "run_id": run_dir.name,
        "result_dir": str(run_dir),
        "result_class": result_class,
        "formal_benchmark": _is_formal_benchmark(result_class),
    }


# ---------------------------------------------------------------------------
# core_summary builders
# ---------------------------------------------------------------------------

def _lineage_core_rows(base: Dict, lm: Dict, source_json: str) -> List[Dict]:
    label_mode = str(lm.get("label_mode") or "")
    provider_id = str(lm.get("provider_id") or "")
    cell_state_key = str(lm.get("cell_state_key") or "")
    gt = lm.get("ground_truth") or {}
    reference_graph_path = str(gt.get("reference_graph_path") or "")
    flags = _make_analysis_flags(label_mode, provider_id, cell_state_key)
    rows: List[Dict] = []
    for json_key, metric_name in LINEAGE_METRIC_NAMES_MAP.items():
        val = lm.get(json_key)
        if val is None:
            continue
        row = dict(base)
        row.update({
            "metric_family": "lineage_fidelity",
            "metric_name": metric_name,
            "metric_value": val,
            "provider_id": provider_id,
            "label_mode": label_mode,
            "cell_state_key": cell_state_key,
            "reference_graph_path": reference_graph_path,
            "source_json": source_json,
        })
        row.update(flags)
        rows.append(row)
    return rows


def _embedding_core_rows(base: Dict, emb: Dict, mode: str, source_json: str) -> List[Dict]:
    label_mode = str(emb.get("label_mode") or mode)
    provider_id = str(emb.get("provider_id") or "")
    state_key = str(emb.get("state_key") or "")
    flags = _make_analysis_flags(label_mode, provider_id, state_key)
    rows: List[Dict] = []
    for json_key, metric_name in EMBEDDING_METRIC_NAMES_MAP.items():
        val = emb.get(json_key)
        if val is None:
            continue
        row = dict(base)
        row.update({
            "metric_family": "embedding_coherence",
            "metric_name": metric_name,
            "metric_value": val,
            "provider_id": provider_id,
            "label_mode": label_mode,
            "cell_state_key": state_key,
            "reference_graph_path": "",
            "source_json": source_json,
        })
        row.update(flags)
        rows.append(row)
    return rows


def _collect_core_rows(method: str, run_dir: Path, run_meta: Optional[Dict]) -> List[Dict]:
    base = _base_row(method, run_dir, run_meta)
    rows: List[Dict] = []

    # Generic lineage_metrics.json (only if label_mode is set)
    lm_path = run_dir / "lineage_metrics.json"
    lm = _load_json_safe(lm_path)
    if lm and lm.get("label_mode"):
        rows.extend(_lineage_core_rows(_base_row(method, run_dir, run_meta, lm), lm, str(lm_path)))

    # Mode-specific lineage metrics
    for mode in EMBEDDING_MODES:
        lm_mode_path = run_dir / f"lineage_metrics_{mode}.json"
        lm_mode = _load_json_safe(lm_mode_path)
        if lm_mode and lm_mode.get("label_mode"):
            rows.extend(_lineage_core_rows(_base_row(method, run_dir, run_meta, lm_mode), lm_mode, str(lm_mode_path)))

    # Embedding milestone eval. If the run itself is legacy scGPT, ignore any
    # later diagnostic milestone embedding files left in that historical folder.
    emb_dir = run_dir / "embedding_milestone_eval"
    if emb_dir.is_dir() and not _run_is_legacy_scgpt(run_meta):
        for mode in EMBEDDING_MODES:
            emb_path = emb_dir / f"embedding_metrics_{mode}.json"
            emb = _load_json_safe(emb_path)
            if emb is not None and _embedding_matches_run_source(emb, run_meta):
                rows.extend(_embedding_core_rows(_base_row(method, run_dir, run_meta, emb), emb, mode, str(emb_path)))

    return rows


# ---------------------------------------------------------------------------
# lineage_summary builders
# ---------------------------------------------------------------------------

def _collect_lineage_rows(method: str, run_dir: Path, run_meta: Optional[Dict]) -> List[Dict]:
    rows: List[Dict] = []
    lineage_files: List[Path] = [run_dir / "lineage_metrics.json"]
    for mode in EMBEDDING_MODES:
        lineage_files.append(run_dir / f"lineage_metrics_{mode}.json")

    for lm_path in lineage_files:
        lm = _load_json_safe(lm_path)
        if not lm:
            continue
        label_mode = str(lm.get("label_mode") or "")
        if not label_mode:
            continue  # Not milestone-tagged

        provider_id = str(lm.get("provider_id") or "")
        cell_state_key = str(lm.get("cell_state_key") or "")
        gt = lm.get("ground_truth") or {}
        reference_graph_path = str(gt.get("reference_graph_path") or "")
        flags = _make_analysis_flags(label_mode, provider_id, cell_state_key)
        base = _base_row(method, run_dir, run_meta, lm)
        row = dict(base)
        row.update({
            "provider_id": provider_id,
            "label_mode": label_mode,
            "cell_state_key": cell_state_key,
            "reference_graph_path": reference_graph_path,
            "n_reference_edges": _safe(lm.get("n_reference_edges")),
            "auroc": _safe(lm.get("auroc")),
            "auprc": _safe(lm.get("auprc")),
            "jaccard_topk": _safe(lm.get("jaccard_similarity_topk")),
            "status": _safe(lm.get("status")),
            "source_json": str(lm_path),
        })
        row.update(flags)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# embedding_summary builders
# ---------------------------------------------------------------------------

def _collect_embedding_rows(method: str, run_dir: Path, run_meta: Optional[Dict]) -> List[Dict]:
    if _run_is_legacy_scgpt(run_meta):
        return []
    emb_dir = run_dir / "embedding_milestone_eval"
    if not emb_dir.is_dir():
        return []
    rows: List[Dict] = []
    for mode in EMBEDDING_MODES:
        emb_path = emb_dir / f"embedding_metrics_{mode}.json"
        emb = _load_json_safe(emb_path)
        if emb is None:
            continue
        if not _embedding_matches_run_source(emb, run_meta):
            continue
        label_mode = str(emb.get("label_mode") or mode)
        provider_id = str(emb.get("provider_id") or "")
        state_key = str(emb.get("state_key") or "")
        embedding_key = str(emb.get("embedding_key") or "")
        cluster_source = str(emb.get("cluster_source") or "")
        flags = _make_analysis_flags(label_mode, provider_id, state_key)
        base = _base_row(method, run_dir, run_meta, emb)
        row = dict(base)
        row.update({
            "provider_id": provider_id,
            "label_mode": label_mode,
            "state_key": state_key,
            "embedding_key": embedding_key,
            "cluster_source": cluster_source,
            "adjusted_rand_index": _safe(emb.get("adjusted_rand_index")),
            "mean_prediction_entropy": _safe(emb.get("mean_prediction_entropy")),
            "weighted_prediction_entropy": _safe(emb.get("weighted_prediction_entropy")),
            "entropy_basis": _safe(emb.get("entropy_basis")),
            "n_probability_classes": _safe(emb.get("n_probability_classes")),
            "prediction_entropy_note": _safe(emb.get("prediction_entropy_note")),
            "mean_normalized_entropy": _safe(emb.get("mean_normalized_entropy")),
            "weighted_mean_normalized_entropy": _safe(
                emb.get("weighted_mean_normalized_entropy")),
            "hard_label_entropy_note": _safe(emb.get("entropy_note") or emb.get("note")),
            "n_cells_evaluated": _safe(emb.get("n_cells_evaluated")),
            "status": _safe(emb.get("ari_status")),
            "source_json": str(emb_path),
        })
        row.update(flags)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# provider_agreement_summary builders
# ---------------------------------------------------------------------------

def _collect_agreement_rows(method: str, run_dir: Path, run_meta: Optional[Dict]) -> List[Dict]:
    if _run_is_legacy_scgpt(run_meta):
        return []
    pa_path = run_dir / "embedding_milestone_eval" / "provider_agreement_metrics.json"
    pa = _load_json_safe(pa_path)
    if not pa:
        return []
    if not _embedding_matches_run_source(pa, run_meta):
        return []
    dataset_id = _infer_dataset_id(run_meta, pa, run_dir)
    n_cells_compared = _safe(pa.get("n_cells_compared"))
    ari_map: Dict = pa.get("pairwise_adjusted_rand_index") or {}
    emf_map: Dict = pa.get("pairwise_exact_match_fraction") or {}
    all_pairs = sorted(set(list(ari_map.keys()) + list(emf_map.keys())))
    rows: List[Dict] = []
    for pair in all_pairs:
        rows.append({
            "dataset_id": dataset_id,
            "result_dir": str(run_dir),
            "pair": pair,
            "pairwise_adjusted_rand_index": _safe(ari_map.get(pair)),
            "pairwise_exact_match_fraction": _safe(emf_map.get(pair)),
            "n_cells_compared": n_cells_compared,
            "source_json": str(pa_path),
        })
    return rows


# ---------------------------------------------------------------------------
# CSV writer
# ---------------------------------------------------------------------------

def _write_csv(path: Path, rows: List[Dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


# ---------------------------------------------------------------------------
# Sort helpers
# ---------------------------------------------------------------------------

def _row_sort_key(r: Dict) -> Tuple:
    return (r.get("method", ""), r.get("scenario", ""),
            r.get("label_mode", ""), r.get("run_id", ""))


def _core_sort_key(r: Dict) -> Tuple:
    return (r.get("method", ""), r.get("scenario", ""),
            r.get("label_mode", ""), r.get("metric_family", ""),
            r.get("metric_name", ""), r.get("run_id", ""))


def _count_data_rows(path: Path) -> int:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
        return max(0, len(lines) - 1)
    except Exception:
        return 0


# ---------------------------------------------------------------------------
# Main runner
# ---------------------------------------------------------------------------

def run_summary(results_dir: str, output_dir: str) -> Dict[str, Any]:
    results_root = Path(results_dir)
    out_root = Path(output_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    print(f"[summarize_milestone] Scanning {results_root} ...")
    runs = _discover_milestone_runs(results_root)
    print(f"[summarize_milestone] Found {len(runs)} milestone run directories.")

    core_rows: List[Dict] = []
    lineage_rows: List[Dict] = []
    embedding_rows: List[Dict] = []
    agreement_rows: List[Dict] = []

    for method, run_dir in runs:
        run_meta = _load_json_safe(run_dir / "run_metadata.json") or {}
        core_rows.extend(_collect_core_rows(method, run_dir, run_meta))
        lineage_rows.extend(_collect_lineage_rows(method, run_dir, run_meta))
        embedding_rows.extend(_collect_embedding_rows(method, run_dir, run_meta))
        agreement_rows.extend(_collect_agreement_rows(method, run_dir, run_meta))

    core_rows.sort(key=_core_sort_key)
    lineage_rows.sort(key=_row_sort_key)
    embedding_rows.sort(key=_row_sort_key)
    agreement_rows.sort(key=lambda r: (r.get("result_dir", ""), r.get("pair", "")))

    paths = {
        "core_summary": out_root / "core_summary.csv",
        "lineage_summary": out_root / "lineage_summary.csv",
        "embedding_summary": out_root / "embedding_summary.csv",
        "provider_agreement_summary": out_root / "provider_agreement_summary.csv",
    }
    _write_csv(paths["core_summary"], core_rows, CORE_FIELDS)
    _write_csv(paths["lineage_summary"], lineage_rows, LINEAGE_FIELDS)
    _write_csv(paths["embedding_summary"], embedding_rows, EMBEDDING_FIELDS)
    _write_csv(paths["provider_agreement_summary"], agreement_rows, AGREEMENT_FIELDS)

    for name, path in paths.items():
        n = _count_data_rows(path)
        print(f"[summarize_milestone] Wrote {name}: {n} data row(s) -> {path}")

    return {
        "n_runs": len(runs),
        "n_core_rows": len(core_rows),
        "n_lineage_rows": len(lineage_rows),
        "n_embedding_rows": len(embedding_rows),
        "n_agreement_rows": len(agreement_rows),
        "paths": {k: str(v) for k, v in paths.items()},
    }


def _print_summary_table(result: Dict) -> None:
    print(f"\n{'=' * 60}")
    print(f"  Milestone Summary -- {result['n_runs']} milestone run(s) found")
    print(f"{'=' * 60}")
    print(f"  core_summary rows          : {result['n_core_rows']}")
    print(f"  lineage_summary rows       : {result['n_lineage_rows']}")
    print(f"  embedding_summary rows     : {result['n_embedding_rows']}")
    print(f"  provider_agreement rows    : {result['n_agreement_rows']}")
    print(f"{'=' * 60}")
    for name, path in result["paths"].items():
        print(f"  {name:35s}: {path}")
    print()


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--results-dir", default=None,
        help="Root directory to scan. Defaults to <project_root>/benchmark/results.",
    )
    parser.add_argument(
        "--output-dir", default=None,
        help="Output directory for CSVs. Defaults to <project_root>/benchmark/reports.",
    )
    args = parser.parse_args(argv)

    root = _project_root()
    results_dir = args.results_dir or str(root / "benchmark" / "results")
    output_dir = args.output_dir or str(root / "benchmark" / "reports")

    try:
        result = run_summary(results_dir=results_dir, output_dir=output_dir)
    except Exception as exc:
        print(f"[summarize_milestone] ERROR: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

    _print_summary_table(result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
