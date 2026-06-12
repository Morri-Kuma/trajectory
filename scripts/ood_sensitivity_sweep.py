#!/usr/bin/env python
"""k / quantile sensitivity sweep for the geometric OOD reference gate.

Reviewer request (R1, P1.2): show the accept/reject decision of the geometric OOD gate
is not an artifact of the k and quantile choices. Sweeps k in {5,15,30,50} and quantile
in {0.90,0.95,0.99} on the SAME real prepped inputs and the SAME PCA-latent proxy used
by scripts/demo_ood_recalibration.py (a local CPU stand-in for the server scANVI latent).

Because the reference self-flags exactly (1-q) of its own cells by construction, the
decision-relevant, quantile-aware statistic is ENRICHMENT = frac_ood / (1-q): how many
times more often a query cell is flagged than a reference cell.

PCA-latent proxy, NOT the server scANVI latent; absolute fractions are lower than the
scANVI-latent values in the main text (0.73-0.76 at k=15,q=0.95), but enrichment is
large and k-stable across q in {0.90, 0.95}.

Usage: python scripts/ood_sensitivity_sweep.py
Outputs: results/test_outputs/ood_recalibration/ood_sensitivity_sweep.{csv,json}
"""
from __future__ import annotations
import os, json, sys
import numpy as np
import anndata as ad
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.annotation import ood_score_and_flag

OUT = "results/test_outputs/ood_recalibration"
os.makedirs(OUT, exist_ok=True)
N = 10000
REF = "results/annotation_branch/inputs/gse242424_ref_fullgene_labelled.h5ad"
QUERIES = {
    "gse178325": "results/annotation_branch/inputs/gse178325_query_fullgene_labelled.h5ad",
    "gse230659": "results/annotation_branch/inputs/gse230659_query_fullgene_labelled.h5ad",
}
K_GRID = [5, 15, 30, 50]
Q_GRID = [0.90, 0.95, 0.99]
ENRICH_REJECT = 3.0


def lognorm(A):
    import scipy.sparse as sp
    X = A.layers["counts"] if "counts" in A.layers else A.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X, float)
    lib = X.sum(1, keepdims=True)
    lib[lib == 0] = 1
    return np.log1p(X / lib * 1e4).astype(np.float32)


def sub(A, n, seed=0):
    if A.n_obs <= n:
        return A
    idx = np.random.default_rng(seed).choice(A.n_obs, n, replace=False)
    idx.sort()
    return A[idx].copy()


def main():
    ref = sub(ad.read_h5ad(REF), N)
    Xr = lognorm(ref)
    scaler = StandardScaler(with_mean=True).fit(Xr)
    pca = PCA(n_components=50, random_state=0).fit(scaler.transform(Xr))
    ref_lat = pca.transform(scaler.transform(Xr)).astype(np.float32)

    rows = []
    for name, path in QUERIES.items():
        q_lat = pca.transform(scaler.transform(lognorm(sub(ad.read_h5ad(path), N)))).astype(np.float32)
        for k in K_GRID:
            for qn in Q_GRID:
                res = ood_score_and_flag(ref_lat, q_lat, k=k, quantile=qn)
                enrich = res["frac_ood"] / (1.0 - qn)
                rows.append({
                    "dataset": name, "k": k, "quantile": qn,
                    "frac_ood": round(res["frac_ood"], 4),
                    "enrichment_vs_ref": round(enrich, 2),
                    "threshold": round(res["threshold"], 4),
                    "decision": "reject" if enrich > ENRICH_REJECT else "accept",
                })
                print("[%s] k=%2d q=%.2f -> frac_ood=%5.1f%%  enrich=%5.1fx  %s"
                      % (name, k, qn, res["frac_ood"] * 100, enrich, rows[-1]["decision"]), flush=True)

    cols = ["dataset", "k", "quantile", "frac_ood", "enrichment_vs_ref", "threshold", "decision"]
    with open(os.path.join(OUT, "ood_sensitivity_sweep.csv"), "w") as f:
        f.write(",".join(cols) + "\n")
        for r in rows:
            f.write(",".join(str(r[c]) for c in cols) + "\n")

    summ = {"latent": "PCA-latent proxy (local CPU stand-in for server scANVI latent)",
            "reject_rule": "enrichment_vs_ref > 3x", "by_dataset": {}}
    for name in QUERIES:
        ds = [r for r in rows if r["dataset"] == name]
        decisions = {r["decision"] for r in ds}
        enr = [r["enrichment_vs_ref"] for r in ds]
        summ["by_dataset"][name] = {
            "n_settings": len(ds),
            "unanimous_decision": (list(decisions)[0] if len(decisions) == 1 else "MIXED"),
            "enrichment_min": min(enr), "enrichment_max": max(enr),
        }
    with open(os.path.join(OUT, "ood_sensitivity_sweep.json"), "w") as f:
        json.dump(summ, f, indent=2)

    print("\n=== decision stability (PCA-latent proxy, enrichment rule) ===")
    for name, s in summ["by_dataset"].items():
        print("[%s] %d (k,q) settings -> %s (enrichment %.1f-%.1fx over ref baseline)"
              % (name, s["n_settings"], s["unanimous_decision"], s["enrichment_min"], s["enrichment_max"]))


if __name__ == "__main__":
    main()
