import numpy as np
from scipy import sparse as sp
from src.data import reconstruct_counts

def _make(seed, shape, lam):
    rng = np.random.default_rng(seed)
    C = rng.poisson(lam, size=shape).astype(float)
    tot = C.sum(1); tot[tot == 0] = 1.0
    norm = np.log1p(C / tot[:, None] * 1e4)
    return C, tot, norm

def test_roundtrip_dense_exact():
    C, tot, norm = _make(0, (40, 12), 3)
    rec = reconstruct_counts(norm, tot)
    assert np.allclose(rec, C)

def test_roundtrip_sparse_exact_and_integer():
    C, tot, norm = _make(1, (60, 20), 2)
    rec = reconstruct_counts(sp.csr_matrix(norm), tot)
    assert sp.issparse(rec)
    assert np.allclose(rec.toarray(), C)
    assert np.allclose(rec.data, np.round(rec.data))
