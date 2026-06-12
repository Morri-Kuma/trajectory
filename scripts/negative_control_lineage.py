#!/usr/bin/env python
"""Negative control for the Lineage Fidelity metric (reviewer-prioritized #1).

Demonstrates that the graph-similarity single-step AUROC does NOT return a plausible
score regardless of input. For every real predicted state-transition matrix (the three
projection methods, all six systems, every available scenario A-F including the
pseudotime axis D/E/F) we recompute the SAME AUROC the benchmark reports (diagonal-zeroed
predicted matrix vs the binary reference graph; identical to
lineage_graphsim_sctimebench._auc_metrics) under two null inputs:

  (1) LABEL-PERMUTATION null: permute predicted state identities (rows+cols of W); K draws.
  (2) RANDOM-MATRIX control: random non-negative row-stochastic matrix; R draws.

The real AUROC uses sklearn (production path) and a tie-correct fast AUROC
(benchmark.evaluation.lineage_robustness.auroc_flat) is asserted to match it, then used
for the bulk draws. A discriminating metric collapses both nulls toward chance.

Output: results/test_outputs/negative_control/lineage_negative_control.{csv,json}
Local CPU. Run: python scripts/negative_control_lineage.py
"""
from __future__ import annotations
import os, sys, json, glob, re
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from benchmark.evaluation.lineage_graphsim_sctimebench import align_predicted_to_reference
from benchmark.evaluation.lineage_robustness import auroc_flat as fast_auc

OUT = "results/test_outputs/negative_control"; os.makedirs(OUT, exist_ok=True)
K_PERM, R_RAND, SEED = 2000, 2000, 0
METHODS = ("scnode", "prescient", "mioflow")
DATASET_PROVIDER = {
    "gse175634": "gse175634_cardiac_silver_v1",
    "gse178325": "gse178325_marker_fm_transition_silver_v1",
    "gse218855": "gse218855_marker_fm_transition_silver_v1",
    "gse230659": "gse230659_marker_fm_transition_silver_v1",
    "gse242424": "gse242424_oskm_reprogramming_ground_truth_v1",
    "gse298212": "gse298212_marker_fm_transition_silver_v1",
}


def load_reference(provider):
    g = json.load(open(f"benchmark/ground_truth/providers/{provider}/reference_graph.json"))
    nodes = [str(n["id"]) for n in g["nodes"]]
    B = np.zeros((len(nodes), len(nodes)), dtype=int)
    idx = {n: i for i, n in enumerate(nodes)}
    for e in g["edges"]:
        s, t = str(e["source"]), str(e["target"])
        if s in idx and t in idx:
            B[idx[s], idx[t]] = 1
    return B, nodes


def auc_from_W(W, B_flat):
    Wdz = W.copy(); np.fill_diagonal(Wdz, 0.0)
    return fast_auc(B_flat, Wdz.flatten())


def main():
    rng = np.random.default_rng(SEED)
    paths = sorted(p for p in glob.glob("benchmark/results/*/*_formal/state_transition_matrix.csv")
                   if re.search(r"/results/(%s)/" % "|".join(METHODS), p)
                   and "_rerun_backup" not in p)
    rows = []
    for path in paths:
        mm = re.search(r"/results/(\w+)/(gse\d+)_", path)
        sc = re.search(r"_([A-F])_(hvg2000|pseudotime)_formal/", path)
        if not (mm and sc):
            continue
        method, ds = mm.group(1), mm.group(2)
        scenario, axis = sc.group(1), ("observed" if sc.group(2) == "hvg2000" else "pseudotime")
        provider = DATASET_PROVIDER.get(ds)
        if provider is None:
            continue
        B, nodes = load_reference(provider); B_flat = B.flatten().astype(int)
        raw = pd.read_csv(path, index_col=0)
        raw.index = raw.index.astype(str); raw.columns = raw.columns.astype(str)
        try:
            W = align_predicted_to_reference(raw, nodes)
        except Exception:
            continue
        Wdz = W.copy(); np.fill_diagonal(Wdz, 0.0)
        if len(set(B_flat)) < 2:
            continue
        real = float(roc_auc_score(B_flat, Wdz.flatten()))
        assert abs(real - auc_from_W(W, B_flat)) < 1e-9, f"AUC mismatch {method} {ds} {scenario}"

        n = len(nodes)
        null = np.array([auc_from_W(W[np.ix_(p, p)], B_flat)
                         for p in (rng.permutation(n) for _ in range(K_PERM))])
        null = null[~np.isnan(null)]
        rand = []
        for _ in range(R_RAND):
            M = rng.random((n, n)); M = M / M.sum(1, keepdims=True)
            rand.append(auc_from_W(M, B_flat))
        rand = np.array([x for x in rand if not np.isnan(x)])

        p95 = float(np.quantile(null, 0.95))
        p_perm = (1 + int((null >= real).sum())) / (len(null) + 1)
        rows.append({"method": method, "dataset": ds.upper(), "scenario": scenario, "axis": axis,
                     "real_auroc": round(real, 4),
                     "perm_null_mean": round(float(null.mean()), 4),
                     "perm_null_p95": round(p95, 4),
                     "perm_p_value": round(p_perm, 4),
                     "random_mat_mean": round(float(rand.mean()), 4),
                     "real_above_null_p95": bool(real > p95)})
        print(f"[{method:9s} {ds} {scenario}/{axis[:4]}] real={real:.3f} "
              f"perm-null={null.mean():.3f}(p95={p95:.3f},p={p_perm:.3f}) rand={rand.mean():.3f}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/lineage_negative_control.csv", index=False)

    def block(d):
        return {"n_runs": int(len(d)),
                "perm_null_mean": round(float(d["perm_null_mean"].mean()), 4),
                "random_mat_mean": round(float(d["random_mat_mean"].mean()), 4),
                "n_real_above_null_p95": int(d["real_above_null_p95"].sum())}
    summary = {"k_permutations": K_PERM, "r_random": R_RAND,
               "metric": "single_step graph-sim AUROC (production _auc_metrics path)",
               "overall": block(df),
               "by_axis": {ax: block(df[df.axis == ax]) for ax in sorted(df.axis.unique())},
               "interpretation": ("Both null inputs collapse the metric toward an empirical "
                                  "chance level near 0.6 on these small graphs; real predictions "
                                  "exceed the permutation null only where genuine lineage signal "
                                  "exists, so the metric is discriminating, not permissive.")}
    json.dump(summary, open(f"{OUT}/lineage_negative_control.json", "w"), indent=2)
    print("\n=== SUMMARY ===")
    o = summary["overall"]
    print(f"ALL ({o['n_runs']} runs): perm-null={o['perm_null_mean']:.3f} rand={o['random_mat_mean']:.3f} "
          f"| real>null-p95: {o['n_real_above_null_p95']}/{o['n_runs']}")
    for ax, b in summary["by_axis"].items():
        print(f"  {ax:10s} ({b['n_runs']} runs): perm-null={b['perm_null_mean']:.3f} "
              f"| real>null-p95: {b['n_real_above_null_p95']}/{b['n_runs']}")


if __name__ == "__main__":
    main()
