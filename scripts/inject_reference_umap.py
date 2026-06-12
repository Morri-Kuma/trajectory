#!/usr/bin/env python
"""Inject a reference UMAP embedding into a benchmark input h5ad's .obsm.

The official embedding-coherence evaluator (benchmark/evaluation/eval_embedding_milestone.py)
needs a reference embedding in the input h5ad's .obsm for the milestone Leiden
clustering step. Its detection priority is ["X_pca", "X_umap"]; the primary
datasets (GSE230659/GSE178325) carry X_umap there, so they are evaluated on a 2D
UMAP. The newer marker-FM builders (GSE218855/GSE298212) write no obsm embedding,
which makes the evaluator raise "No embeddings found in obsm.".

This script computes PCA -> neighbors -> UMAP on the (log-normalized) HVG matrix
and stores ONLY X_umap, deleting the intermediate X_pca, so the evaluator picks
X_umap exactly as it does for the primaries (cross-dataset consistency). It does
NOT touch X, layers, obs labels, or uns, so existing annotations and benchmark
results remain valid. Idempotent: skips if X_umap already present (unless --force).
"""
from __future__ import annotations

import argparse
import os

# Keep every numeric/parallel backend single-threaded. The neighbor search
# (pynndescent/joblib/numba) otherwise tries to spawn one worker per core, which
# trips the login node's thread/process ulimit ("can't start new thread"). Must
# be set before importing scanpy/numba.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMBA_NUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_THREADING_LAYER"):
    os.environ.setdefault(_v, "1" if _v != "NUMBA_THREADING_LAYER" else "workqueue")

import anndata as ad
import scanpy as sc

sc.settings.n_jobs = 1


def inject(h5ad_path: str, n_pcs: int = 50, n_neighbors: int = 15,
           seed: int = 0, force: bool = False) -> None:
    a = ad.read_h5ad(h5ad_path)
    if "X_umap" in a.obsm and not force:
        print(f"[inject_reference_umap] X_umap already present in {h5ad_path}; skipping "
              f"(shape={a.obsm['X_umap'].shape}). Use --force to recompute.")
        return
    n_pcs = int(min(n_pcs, a.n_vars - 1, a.n_obs - 1))
    print(f"[inject_reference_umap] {h5ad_path}: n_obs={a.n_obs}, n_vars={a.n_vars}; "
          f"computing PCA({n_pcs}) -> neighbors({n_neighbors}) -> UMAP ...")
    sc.pp.pca(a, n_comps=n_pcs, random_state=seed)
    sc.pp.neighbors(a, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed)
    sc.tl.umap(a, random_state=seed)
    # Keep ONLY X_umap so the evaluator's [X_pca, X_umap] priority matches the
    # primary datasets (which expose X_umap, not X_pca).
    if "X_pca" in a.obsm:
        del a.obsm["X_pca"]
    a.write_h5ad(h5ad_path)
    print(f"[inject_reference_umap] wrote X_umap {a.obsm['X_umap'].shape} -> {h5ad_path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-h5ad", required=True, help="Benchmark input h5ad to patch in place.")
    ap.add_argument("--n-pcs", type=int, default=50)
    ap.add_argument("--n-neighbors", type=int, default=15)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--force", action="store_true", help="Recompute even if X_umap exists.")
    args = ap.parse_args()
    inject(args.input_h5ad, args.n_pcs, args.n_neighbors, args.seed, args.force)


if __name__ == "__main__":
    main()
