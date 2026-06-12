"""
Anti-leakage tests for build_representation_inputs.

Core guarantee (Work-Plan Step 2 / Risk 2): with reducer_fit_scope='train_only',
the scaler/PCA must be fit on TRAINING cells only — held-out cells must not
influence the fitted transform.
"""

import numpy as np
import pytest

from benchmark.representations.build_representation_inputs import (
    ALL_CELLS_EXPLORATORY,
    TRAIN_ONLY,
    build_hvg_pca,
    build_scfm_pca,
)


def _toy_embedding(n, d, seed):
    rng = np.random.default_rng(seed)
    return rng.normal(size=(n, d)).astype(np.float64)


def test_scfm_train_only_ignores_heldout_cells():
    """
    Fitting on train cells only must give an identical transform regardless of
    what the held-out cells contain. We build two datasets that share the same
    training rows but have wildly different held-out rows; the X_rep values on
    the training rows must match to numerical precision.
    """
    d = 8
    train = _toy_embedding(40, d, seed=1)

    held_a = _toy_embedding(20, d, seed=2)
    held_b = _toy_embedding(20, d, seed=3) * 100.0 + 50.0  # very different scale

    emb_a = np.vstack([train, held_a])
    emb_b = np.vstack([train, held_b])
    mask = np.array([True] * 40 + [False] * 20)

    res_a = build_scfm_pca(emb_a, mask, "rep_scgpt_cls_pca50", n_components=5, seed=0)
    res_b = build_scfm_pca(emb_b, mask, "rep_scgpt_cls_pca50", n_components=5, seed=0)

    assert res_a.reducer_fit_scope == TRAIN_ONLY
    # Training-row embeddings must be identical (PCA sign can differ, compare abs)
    np.testing.assert_allclose(
        np.abs(res_a.X_rep[:40]), np.abs(res_b.X_rep[:40]), rtol=1e-5, atol=1e-5
    )


def test_scfm_all_cells_mode_does_depend_on_heldout():
    """Sanity check: exploratory all-cell fitting DOES depend on held-out cells,
    confirming the train-only test above is actually exercising leakage."""
    d = 8
    train = _toy_embedding(40, d, seed=1)
    held_a = _toy_embedding(20, d, seed=2)
    held_b = _toy_embedding(20, d, seed=3) * 100.0 + 50.0
    emb_a = np.vstack([train, held_a])
    emb_b = np.vstack([train, held_b])

    res_a = build_scfm_pca(emb_a, None, "rep_scgpt_cls_pca50", n_components=5)
    res_b = build_scfm_pca(emb_b, None, "rep_scgpt_cls_pca50", n_components=5)
    assert res_a.reducer_fit_scope == ALL_CELLS_EXPLORATORY
    # Training-row embeddings should now differ because the fit saw held-out data
    assert not np.allclose(np.abs(res_a.X_rep[:40]), np.abs(res_b.X_rep[:40]),
                           rtol=1e-3, atol=1e-3)


def test_scfm_output_shape_and_metadata():
    emb = _toy_embedding(30, 10, seed=5)
    mask = np.array([True] * 20 + [False] * 10)
    res = build_scfm_pca(emb, mask, "rep_scfoundation_pca50", n_components=4, seed=0)
    assert res.X_rep.shape == (30, 4)
    assert res.metadata["reducer_fit_scope"] == TRAIN_ONLY
    assert res.metadata["n_train_cells_fit"] == 20
    assert res.metadata["n_total_cells"] == 30
    assert res.metadata["input_space"] == "scfm_embedding"


def test_hvg_pca_train_only_ignores_heldout_cells():
    rng = np.random.default_rng(0)
    counts_train = rng.poisson(2.0, size=(40, 50)).astype(float)
    held_a = rng.poisson(2.0, size=(15, 50)).astype(float)
    held_b = rng.poisson(20.0, size=(15, 50)).astype(float)  # different depth
    counts_a = np.vstack([counts_train, held_a])
    counts_b = np.vstack([counts_train, held_b])
    mask = np.array([True] * 40 + [False] * 15)

    # Provide fixed HVG indices so the comparison isolates scaler/PCA leakage.
    hvg = np.arange(30)
    res_a = build_hvg_pca(counts_a, mask, n_components=5, hvg_indices=hvg, seed=0)
    res_b = build_hvg_pca(counts_b, mask, n_components=5, hvg_indices=hvg, seed=0)
    assert res_a.reducer_fit_scope == TRAIN_ONLY
    np.testing.assert_allclose(
        np.abs(res_a.X_rep[:40]), np.abs(res_b.X_rep[:40]), rtol=1e-5, atol=1e-5
    )


def test_empty_train_mask_raises():
    emb = _toy_embedding(10, 4, seed=1)
    with pytest.raises(ValueError, match="zero cells"):
        build_scfm_pca(emb, np.zeros(10, dtype=bool), "rep_scgpt_cls_pca50",
                       n_components=2)
