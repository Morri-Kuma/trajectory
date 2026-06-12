import numpy as np
import anndata as ad
from src.preprocessing import basic_preprocess

def test_sklearn_preprocess_outputs():
    rng = np.random.default_rng(0)
    A = ad.AnnData(X=rng.poisson(2, size=(50, 30)).astype(float))
    A = basic_preprocess(A, n_pcs=10, n_neighbors=8, prefer_scanpy=False)
    assert A.obsm["X_pca"].shape == (50, 10)
    assert "counts" in A.layers
    assert A.uns["preprocess"]["backend"] == "sklearn"
