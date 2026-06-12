"""Pseudotime + state-graph (PAGA-like) inference.

Production path uses scanpy (DPT, PAGA). A dependency-light surrogate runs the
same idea (geodesic distance from a root on the kNN graph; inter-group
connectivity) so robustness comparisons execute without scanpy.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse.csgraph import dijkstra


def _root_indices(adata, root_state: str, label_key: str) -> np.ndarray:
    labels = adata.obs[label_key].astype(str).values
    idx = np.where(labels == root_state)[0]
    if len(idx) == 0:  # fall back to first cell
        idx = np.array([0])
    return idx


def compute_pseudotime(
    adata, root_state: str, label_key: str, n_dcs: int = 15, prefer_scanpy: bool = True
) -> np.ndarray:
    """Return a per-cell pseudotime in [0, 1]."""
    if prefer_scanpy:
        try:
            import scanpy as sc

            roots = _root_indices(adata, root_state, label_key)
            adata.uns["iroot"] = int(roots[0])
            if "neighbors" not in adata.uns:
                sc.pp.neighbors(adata, use_rep="X_pca")
            sc.tl.diffmap(adata, n_comps=min(n_dcs, adata.obsm["X_pca"].shape[1]))
            sc.tl.dpt(adata)
            return _normalize01(np.asarray(adata.obs["dpt_pseudotime"].values, dtype=float))
        except Exception:
            pass
    return _geodesic_pseudotime(adata, root_state, label_key)


def _build_knn_graph(adata):
    from sklearn.neighbors import NearestNeighbors

    X = adata.obsm["X_pca"]
    k = int(min(15, adata.n_obs - 1))
    d, idx = NearestNeighbors(n_neighbors=k).fit(X).kneighbors(X)
    rows = np.repeat(np.arange(adata.n_obs), k)
    return sparse.csr_matrix((d.ravel(), (rows, idx.ravel())),
                             shape=(adata.n_obs, adata.n_obs)), idx


def _geodesic_pseudotime(adata, root_state, label_key) -> np.ndarray:
    """Surrogate DPT: shortest-path distance from root cells on the kNN graph."""
    g = adata.obsp["distances"] if "distances" in adata.obsp else None
    if g is None or g.shape[0] != adata.n_obs:
        g, _ = _build_knn_graph(adata)
    g = g.maximum(g.T)  # symmetrize
    roots = _root_indices(adata, root_state, label_key)
    dist = dijkstra(g, directed=False, indices=roots, min_only=True)
    finite = np.isfinite(dist)
    if finite.any():
        dist[~finite] = np.nanmax(dist[finite])
    else:
        dist[:] = 0.0
    return _normalize01(dist)


def _normalize01(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    finite = np.isfinite(x)
    if not finite.any():
        return np.zeros_like(x)
    lo, hi = np.nanmin(x[finite]), np.nanmax(x[finite])
    if hi - lo < 1e-12:
        return np.zeros_like(x)
    return (x - lo) / (hi - lo)


def compute_state_graph(adata, label_key: str, prefer_scanpy: bool = True) -> Dict:
    """Return {'connectivity': DataFrame, 'edges': [(s1,s2,w)], 'backend': str}.

    Surrogate PAGA: connectivity(s1,s2) = kNN edges between groups s1 and s2,
    normalized by group sizes then by the max.
    """
    if prefer_scanpy:
        try:
            import scanpy as sc

            if "neighbors" not in adata.uns:
                sc.pp.neighbors(adata, use_rep="X_pca")
            sc.tl.paga(adata, groups=label_key)
            conn = np.asarray(adata.uns["paga"]["connectivities"].todense())
            cats = list(adata.obs[label_key].astype("category").cat.categories)
            df = pd.DataFrame(conn, index=cats, columns=cats)
            return {"connectivity": df, "edges": _edges_from_conn(df, 0.05), "backend": "scanpy"}
        except Exception:
            pass

    labels = adata.obs[label_key].astype(str).values
    cats = sorted(np.unique(labels))
    cat_idx = {c: i for i, c in enumerate(cats)}
    n = len(cats)
    M = np.zeros((n, n))
    idx = adata.uns.get("knn_indices")
    if idx is None or getattr(idx, "shape", [0])[0] != adata.n_obs:
        _, idx = _build_knn_graph(adata)
    for i in range(adata.n_obs):
        ci = cat_idx[labels[i]]
        for j in idx[i]:
            M[ci, cat_idx[labels[j]]] += 1
    M = (M + M.T) / 2.0
    sizes = np.array([(labels == c).sum() for c in cats]).reshape(-1, 1)
    conn = M / np.maximum(sizes, 1)
    np.fill_diagonal(conn, 0.0)
    if conn.max() > 0:
        conn = conn / conn.max()
    df = pd.DataFrame(conn, index=cats, columns=cats)
    return {"connectivity": df, "edges": _edges_from_conn(df, 0.1), "backend": "surrogate"}


def _edges_from_conn(df, thresh: float) -> List[Tuple[str, str, float]]:
    cats = list(df.index)
    edges = []
    for i in range(len(cats)):
        for j in range(i + 1, len(cats)):
            w = float(df.iloc[i, j])
            if w >= thresh:
                edges.append((cats[i], cats[j], w))
    return edges


def graph_similarity(edges_a, edges_b) -> Dict:
    """Jaccard similarity of undirected edge sets (ignoring weights)."""
    sa = {frozenset((s1, s2)) for s1, s2, _ in edges_a}
    sb = {frozenset((s1, s2)) for s1, s2, _ in edges_b}
    inter = len(sa & sb)
    union = len(sa | sb)
    return {"edge_jaccard": float(inter / union) if union else float("nan"),
            "n_edges_a": len(sa), "n_edges_b": len(sb), "shared_edges": inter}
