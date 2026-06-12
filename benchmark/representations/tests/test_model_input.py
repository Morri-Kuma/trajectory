"""
Tests for representation config parsing and get_model_input_matrix, including
the guarantee that the existing expression-space path is unchanged.
"""

import numpy as np
import pytest

from benchmark.representations.model_input import (
    RepresentationConfig,
    get_model_input_matrix,
    parse_representation_config,
    representation_is_active,
)


class _FakeAdata:
    def __init__(self, X, obsm=None):
        self.X = X
        self.obsm = obsm or {}


def test_no_representation_block_returns_X_unchanged():
    X = np.arange(12, dtype=np.float32).reshape(4, 3)
    adata = _FakeAdata(X, {"X_rep": np.zeros((4, 2))})
    out = get_model_input_matrix(adata, {})
    assert out is X  # exact same object: byte-identical legacy behaviour
    assert not representation_is_active({})


def test_disabled_representation_returns_X():
    X = np.zeros((4, 3), dtype=np.float32)
    adata = _FakeAdata(X, {"X_rep": np.ones((4, 2))})
    cfg = {"representation": {"enabled": False, "input_mode": "obsm",
                             "obsm_key": "X_rep"}}
    assert get_model_input_matrix(adata, cfg) is X


def test_enabled_obsm_returns_representation():
    X = np.zeros((4, 3), dtype=np.float32)
    rep = np.arange(8, dtype=np.float32).reshape(4, 2)
    adata = _FakeAdata(X, {"X_rep": rep})
    cfg = {"representation": {"enabled": True, "input_mode": "obsm",
                             "obsm_key": "X_rep", "final_dim": 2}}
    out = get_model_input_matrix(adata, cfg)
    assert out.shape == (4, 2)
    np.testing.assert_allclose(out, rep)
    assert representation_is_active(cfg)


def test_missing_obsm_key_raises():
    adata = _FakeAdata(np.zeros((4, 3)), {})
    cfg = {"representation": {"enabled": True, "input_mode": "obsm",
                             "obsm_key": "X_rep"}}
    with pytest.raises(KeyError, match="X_rep"):
        get_model_input_matrix(adata, cfg)


def test_final_dim_mismatch_raises():
    adata = _FakeAdata(np.zeros((4, 3)), {"X_rep": np.zeros((4, 5))})
    cfg = {"representation": {"enabled": True, "input_mode": "obsm",
                             "obsm_key": "X_rep", "final_dim": 50}}
    with pytest.raises(ValueError, match="final_dim=50"):
        get_model_input_matrix(adata, cfg)


def test_parse_defaults():
    rep = parse_representation_config(None)
    assert isinstance(rep, RepresentationConfig)
    assert rep.enabled is False
    assert rep.input_mode == "X"
    assert rep.uses_obsm is False


def test_parse_full_block():
    cfg = {"representation": {
        "enabled": True, "input_mode": "obsm", "obsm_key": "X_rep",
        "representation_id": "rep_scgpt_cls_pca50", "input_space": "scfm_embedding",
        "output_space": "representation", "final_dim": 50,
        "reducer_fit_scope": "train_only",
    }}
    rep = parse_representation_config(cfg)
    assert rep.uses_obsm
    assert rep.representation_id == "rep_scgpt_cls_pca50"
    assert rep.final_dim == 50
    assert rep.reducer_fit_scope == "train_only"


def test_parse_invalid_input_mode_raises():
    with pytest.raises(ValueError, match="input_mode"):
        parse_representation_config({"representation": {"input_mode": "bogus"}})


def test_parse_invalid_fit_scope_raises():
    with pytest.raises(ValueError, match="reducer_fit_scope"):
        parse_representation_config(
            {"representation": {"reducer_fit_scope": "leaky"}}
        )
