#!/usr/bin/env python
"""Compute the scTimeBench-style Spearman correlation baseline for Lineage Fidelity.

scTimeBench reports a Correlation baseline (Spearman, maximum-vote) that a method
must beat to demonstrate it has learned more than adjacent-timepoint expression
similarity. The machinery already lives in benchmark/evaluation/eval_lineage.py
(`compute_correlation_baseline`), but the per-run evaluator was invoked without the
AnnData, so the baseline was skipped. This standalone runner supplies the AnnData
and the frozen reference graph for each marker-silver dataset and writes the baseline
graph-sim metrics. CPU-only (no GPU / no method retraining).

Output: benchmark/reports/official_silver/spearman_baseline.csv (+ .json).
Run on Shirokane (or locally) where the input h5ads + providers are present.
"""
from __future__ import annotations

import csv
import glob
import json
import os
import sys
from pathlib import Path

ROOT = Path(os.environ.get("TRAJ_PROJECT_ROOT", os.getcwd()))
sys.path.insert(0, str(ROOT))

import anndata as ad  # noqa: E402
from benchmark.evaluation.eval_lineage import (  # noqa: E402
    load_reference_graph,
    compute_correlation_baseline,
)

# dataset_id -> (input-h5ad glob, provider id). Observed-time (abs_day) axis.
DATASETS = {
    "GSE178325": ("benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/*HVG2000_benchmark_input.h5ad",
                  "gse178325_marker_fm_transition_silver_v1"),
    "GSE230659": ("benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/*HVG2000_benchmark_input.h5ad",
                  "gse230659_marker_fm_transition_silver_v1"),
    "GSE218855": ("benchmark/inputs/gse218855_marker_fm_transition_silver_hvg2000/*HVG2000_benchmark_input.h5ad",
                  "gse218855_marker_fm_transition_silver_v1"),
    "GSE298212": ("benchmark/inputs/gse298212_marker_fm_transition_silver_hvg2000/*HVG2000_benchmark_input.h5ad",
                  "gse298212_marker_fm_transition_silver_v1"),
}
CELL_STATE_KEY = "final_milestone_label_coarse"
TIME_KEY = "abs_day"


def flat(label, d):
    """Pull single/multi-step graph-sim metrics out of the baseline dict."""
    gm = (d or {}).get("graph_metrics") or {}
    ss = gm.get("single_step", {}) or {}
    ms = gm.get("multi_step", {}) or {}
    return {
        "dataset": label,
        "baseline": "spearman_max_vote",
        "single_step_auroc": ss.get("auc_roc", d.get("auroc")),
        "single_step_auprc": ss.get("auc_prc", d.get("auprc")),
        "single_step_jaccard": ss.get("jaccard", d.get("jaccard_similarity")),
        "multi_step_auroc": ms.get("auc_roc"),
        "multi_step_auprc": ms.get("auc_prc"),
        "note": d.get("note", ""),
    }


def main() -> None:
    rows = []
    full = {}
    for ds, (pattern, provider) in DATASETS.items():
        hits = glob.glob(str(ROOT / pattern))
        graph = ROOT / "benchmark" / "ground_truth" / "providers" / provider / "reference_graph.json"
        if not hits or not graph.exists():
            print(f"[spearman_baseline] SKIP {ds}: input or provider missing "
                  f"(h5ad={bool(hits)}, graph={graph.exists()})")
            continue
        print(f"[spearman_baseline] {ds}: {os.path.basename(hits[0])}")
        adata = ad.read_h5ad(hits[0])
        ref_matrix, ref_edges, ref_nodes = load_reference_graph(str(graph), edge_confidence_mode="all")
        res = compute_correlation_baseline(
            adata=adata, reference_matrix=ref_matrix, reference_edges=ref_edges,
            ref_node_ids=ref_nodes, cell_state_key=CELL_STATE_KEY, time_key=TIME_KEY,
        )
        full[ds] = res
        rows.append(flat(ds, res))
        print(f"   single-step baseline AUROC = {rows[-1]['single_step_auroc']}")

    repdir = ROOT / "benchmark" / "reports" / "official_silver"
    repdir.mkdir(parents=True, exist_ok=True)
    keys = ["dataset", "baseline", "single_step_auroc", "single_step_auprc",
            "single_step_jaccard", "multi_step_auroc", "multi_step_auprc", "note"]
    with open(repdir / "spearman_baseline.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in keys})
    (repdir / "spearman_baseline.json").write_text(json.dumps(full, indent=2, default=str))
    print(f"[spearman_baseline] wrote {len(rows)} rows -> {repdir/'spearman_baseline.csv'}")


if __name__ == "__main__":
    main()
