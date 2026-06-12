import numpy as np
import anndata as ad
from scipy.stats import spearmanr
from src.trajectory import graph_similarity, compute_pseudotime, compute_state_graph

def test_graph_similarity_jaccard():
    a = [("x", "y", 0.5), ("y", "z", 0.5)]
    b = [("x", "y", 0.1)]
    assert abs(graph_similarity(a, b)["edge_jaccard"] - 0.5) < 1e-9

def test_pseudotime_monotonic_on_line():
    n = 30
    x = np.linspace(0, 10, n)
    A = ad.AnnData(X=np.zeros((n, 2)))
    A.obsm["X_pca"] = np.c_[x, np.zeros(n)]
    A.obs["lab"] = ["root"] + ["mid"] * (n - 2) + ["end"]
    pt = compute_pseudotime(A, "root", "lab", prefer_scanpy=False)
    assert spearmanr(pt, x).statistic > 0.95

def test_state_graph_surrogate_runs():
    rng = np.random.default_rng(0)
    X = np.vstack([rng.normal(0, 0.1, (20, 2)), rng.normal(6, 0.1, (20, 2))])
    A = ad.AnnData(X=np.zeros((40, 2)))
    A.obsm["X_pca"] = X
    A.obs["lab"] = ["a"] * 20 + ["b"] * 20
    g = compute_state_graph(A, "lab", prefer_scanpy=False)
    assert g["backend"] == "surrogate"
    assert g["connectivity"].shape == (2, 2)
