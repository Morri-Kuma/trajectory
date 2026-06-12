"""Preprocessing: normalize, HVG, PCA, kNN graph.

Production path uses scanpy when available; otherwise a dependency-light
scikit-learn fallback runs the same conceptual steps so the smoke test works
on a CPU-only laptop without scanpy. Both paths populate:
    adata.layers['counts']   raw counts
    adata.X                  log1p-normalized
    adata.obsm['X_pca']      PCA embedding
    adata.obsp['distances']  kNN distance graph (CSR)
    adata.uns['knn_indices'] kNN index matrix (custom key; surrogate trajectory)
"""
from __future__ import annotations

from typing import Optional

import numpy as np
from scipy import sparse


def _to_dense(X):
    return X.toarray() if sparse.issparse(X) else np.asarray(X)


def select_hvg(adata, n_top: Optional[int]):
    """Keep the ``n_top`` most variable genes (variance on log-normalized data).
    No-op if n_top is None or >= n_genes."""
    if n_top is None or n_top >= adata.n_vars:
        return adata
    X = _to_dense(adata.X)
    var = X.var(axis=0)
    keep = np.argsort(var)[::-1][:n_top]
    keep.sort()
    return adata[:, keep].copy()


def basic_preprocess(
    adata,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    max_genes: Optional[int] = None,
    target_sum: float = 1e4,
    seed: int = 0,
    prefer_scanpy: bool = True,
):
    """Normalize -> log1p -> (HVG) -> scale -> PCA -> kNN graph."""
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    used_scanpy = False
    if prefer_scanpy:
        try:
            import scanpy as sc

            sc.pp.normalize_total(adata, target_sum=target_sum)
            sc.pp.log1p(adata)
            if max_genes is not None and max_genes < adata.n_vars:
                sc.pp.highly_variable_genes(adata, n_top_genes=max_genes)
                adata = adata[:, adata.var["highly_variable"]].copy()
            sc.pp.scale(adata, max_value=10)
            sc.tl.pca(adata, n_comps=min(n_pcs, adata.n_vars - 1), random_state=seed)
            sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=seed)
            used_scanpy = True
        except Exception:
            used_scanpy = False

    if not used_scanpy:
        adata = _sklearn_preprocess(adata, n_pcs, n_neighbors, max_genes, target_sum, seed)

    adata.uns.setdefault("preprocess", {})["backend"] = "scanpy" if used_scanpy else "sklearn"
    return adata


def _sklearn_preprocess(adata, n_pcs, n_neighbors, max_genes, target_sum, seed):
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors

    X = _to_dense(adata.layers["counts"]).astype(float)
    libsize = X.sum(axis=1, keepdims=True)
    libsize[libsize == 0] = 1.0
    Xn = np.log1p(X / libsize * target_sum)

    if max_genes is not None and max_genes < Xn.shape[1]:
        var = Xn.var(axis=0)
        keep = np.argsort(var)[::-1][:max_genes]
        keep.sort()
        Xn = Xn[:, keep]
        adata = adata[:, keep].copy()
    adata.X = Xn

    mu = Xn.mean(axis=0)
    sd = Xn.std(axis=0)
    sd[sd == 0] = 1.0
    Xs = np.clip((Xn - mu) / sd, -10, 10)

    n_comps = int(min(n_pcs, Xs.shape[1] - 1, Xs.shape[0] - 1))
    pca = PCA(n_components=n_comps, random_state=seed)
    adata.obsm["X_pca"] = pca.fit_transform(Xs)

    k = int(min(n_neighbors, adata.n_obs - 1))
    nn = NearestNeighbors(n_neighbors=k).fit(adata.obsm["X_pca"])
    dist, idx = nn.kneighbors(adata.obsm["X_pca"])
    adata.obsp["distances"] = _knn_to_sparse(dist, idx, adata.n_obs)
    # custom (non-reserved) keys: anndata will not relocate these into obsp.
    adata.uns["knn_indices"] = idx
    adata.uns["knn_n_neighbors"] = k
    return adata


def _knn_to_sparse(dist, idx, n):
    rows = np.repeat(np.arange(n), idx.shape[1])
    cols = idx.ravel()
    vals = dist.ravel()
    return sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
