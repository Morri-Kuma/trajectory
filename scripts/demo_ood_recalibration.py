#!/usr/bin/env python
"""Demonstrate geometric OOD recalibration on the REAL prepped inputs (local, CPU).

Fits a PCA latent on the GSE242424 (OSKM) reference, projects each chemical query
into it, and scores how far each query cell is from the reference manifold. This
is a local PCA-latent proxy for the server's scANVI latent; it shows the chemical
query cells ARE out-of-distribution — which scANVI's posterior confidence missed
(it flagged ~0.04% OOD while collapsing 76-96% of cells onto `hOSK`).
"""
from __future__ import annotations
import os, json, sys
import numpy as np
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.annotation import ood_score_and_flag, knn_reference_distance

OUT = "results/test_outputs/ood_recalibration"; os.makedirs(OUT, exist_ok=True)
N = 10000  # subsample per dataset for a fast local demo
REF = "results/annotation_branch/inputs/gse242424_ref_fullgene_labelled.h5ad"
QUERIES = {
    "gse178325": "results/annotation_branch/inputs/gse178325_query_fullgene_labelled.h5ad",
    "gse230659": "results/annotation_branch/inputs/gse230659_query_fullgene_labelled.h5ad",
}

def lognorm(A):
    import scipy.sparse as sp
    X = A.layers["counts"] if "counts" in A.layers else A.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X, float)
    lib = X.sum(1, keepdims=True); lib[lib == 0] = 1
    return np.log1p(X / lib * 1e4)

def sub(A, n, seed=0):
    if A.n_obs <= n: return A
    idx = np.random.default_rng(seed).choice(A.n_obs, n, replace=False); idx.sort()
    return A[idx].copy()

ref = sub(ad.read_h5ad(REF), N)
Xr = lognorm(ref)
scaler = StandardScaler(with_mean=True).fit(Xr)
pca = PCA(n_components=50, random_state=0).fit(scaler.transform(Xr))
ref_lat = pca.transform(scaler.transform(Xr))

summary = {}
for name, path in QUERIES.items():
    q = sub(ad.read_h5ad(path), N)
    q_lat = pca.transform(scaler.transform(lognorm(q)))
    res = ood_score_and_flag(ref_lat, q_lat, k=15, quantile=0.95)
    dq, dr = knn_reference_distance(ref_lat, q_lat, k=15)
    summary[name] = {"n_query": int(q.n_obs), "frac_ood_geometric": res["frac_ood"],
                     "query_median_knn_dist": res["query_median_knn_dist"],
                     "ref_median_knn_dist": res["ref_median_knn_dist"],
                     "ood_threshold_p95": res["threshold"],
                     "scanvi_frac_ood_from_run": 0.0004 if name == "gse178325" else 0.0030}
    fig, axA = plt.subplots(figsize=(6, 4))
    axA.hist(dr, bins=60, alpha=0.6, density=True, label="reference (GSE242424, OSKM)")
    axA.hist(dq, bins=60, alpha=0.6, density=True, label=f"query ({name}, chemical)")
    axA.axvline(res["threshold"], color="k", ls="--", lw=1, label="ref p95 OOD threshold")
    axA.set_xlabel("mean distance to 15 nearest reference cells (PCA latent)")
    axA.set_ylabel("density"); axA.set_title(f"Geometric OOD: {name} vs OSKM reference")
    axA.legend(fontsize=7)
    fig.savefig(f"{OUT}/ood_hist_{name}.png", dpi=120, bbox_inches="tight"); plt.close(fig)

json.dump(summary, open(f"{OUT}/ood_recalibration_summary.json", "w"), indent=2)
print("=== Geometric OOD recalibration on REAL data (PCA-latent proxy) ===")
for k, v in summary.items():
    print(f"[{k}] geometric OOD flagged = {v['frac_ood_geometric']*100:.1f}%  "
          f"(scANVI confidence flagged only {v['scanvi_frac_ood_from_run']*100:.2f}%)  "
          f"| query med dist {v['query_median_knn_dist']:.2f} vs ref {v['ref_median_knn_dist']:.2f}")
print(f"outputs -> {OUT}/")
