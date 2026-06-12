import numpy as np
import anndata as ad
from src.evaluation import (composition_compare, agreement_matrix,
                            marker_validation, pseudotime_compare, disagreement_cells)

def _adata(a, b, t, X=None, genes=None):
    n = len(a)
    X = np.zeros((n, 2)) if X is None else np.asarray(X, dtype=float)
    A = ad.AnnData(X=X)
    A.obs["A"] = list(a); A.obs["B"] = list(b); A.obs["time"] = list(t)
    if genes is not None:
        A.var_names = genes
    return A

def test_pseudotime_identical():
    v = np.linspace(0, 1, 25)
    assert abs(pseudotime_compare(v, v)["spearman"] - 1.0) < 1e-9

def test_agreement_perfect_correspondence():
    A = _adata(["hADSCs", "hCiPS", "hADSCs"], ["Fibroblast", "iPSC", "Fibroblast"], ["D0"] * 3)
    res = agreement_matrix(A, "A", "B", {"Fibroblast": "hADSCs", "iPSC": "hCiPS"})
    assert abs(res["metrics"]["correspondence_overall_agreement"] - 1.0) < 1e-9
    assert res["metrics"]["ari"] > 0.9

def test_composition_zero_tv_when_identical():
    A = _adata(["hADSCs", "hADSCs"], ["hADSCs", "hADSCs"], ["D0", "D0"])
    assert composition_compare(A, "A", "B", "time")["overall_tv"] < 1e-9

def test_marker_validation_auroc_one():
    X = [[10, 0], [10, 0], [0, 0], [0, 0]]
    A = _adata(["S", "S", "T", "T"], ["x"] * 4, ["t"] * 4, X=X, genes=["G0", "G1"])
    df = marker_validation(A, {"prog": ["G0"]}, "A")
    auc = df[(df.state == "S") & (df.program == "prog")]["auroc"].iloc[0]
    assert auc > 0.99

def test_disagreement_count():
    A = _adata(["hADSCs", "hCiPS"], ["Fibroblast", "Fibroblast"], ["D0", "D8"])
    d = disagreement_cells(A, "A", "B", {"Fibroblast": "hADSCs", "iPSC": "hCiPS"})
    assert d.shape[0] == 1
