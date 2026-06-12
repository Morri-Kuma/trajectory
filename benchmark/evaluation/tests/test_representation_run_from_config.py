"""
Unit tests for the representation metric backfill helper
(benchmark.evaluation.eval_representation_run_from_config).

These tests are self-contained: they build a tiny synthetic AnnData + config +
result directory and never import the Geneformer/scGPT packages. The numpy
backend is forced so geomloss/torch are not required.
"""

import json

import numpy as np
import pytest

ad = pytest.importorskip("anndata")
yaml = pytest.importorskip("yaml")

from benchmark.evaluation.eval_representation_run_from_config import (
    resolve_eval_timepoints,
    resolve_representation_metadata,
    resolve_state_key,
    resolve_time_key,
    run_from_config,
)


def _make_adata(seed=0, n_per_tp=10, d=4):
    rng = np.random.default_rng(seed)
    times = np.repeat([0.0, 1.0, 2.0, 3.0], n_per_tp)
    n = times.shape[0]
    # representation drifts with time so diagnostics are well-defined
    X_rep = rng.normal(size=(n, d)) + times[:, None]
    states = np.where(times < 2.0, "early", "late")
    adata = ad.AnnData(X=np.zeros((n, 3), dtype=np.float32))
    adata.obs["abs_day"] = times
    adata.obs["final_milestone_label_coarse"] = states
    adata.obsm["X_rep"] = X_rep.astype(np.float32)
    adata.uns["active_representation_id"] = "rep_geneformer_cls_pca50"
    adata.uns["representation_inputs_metadata"] = {
        "rep_geneformer_cls_pca50": {
            "representation_id": "rep_geneformer_cls_pca50",
            "reducer_fit_scope": "train_only",
            "n_components": d,
        }
    }
    adata.uns["scfm_embedding_metadata"] = {
        "geneformer": {
            "model_name": "Geneformer",
            "embedding_dim": 512,
            "emb_mode": "cls",
            "model_version": "V2",
        }
    }
    return adata


# ---------------------------------------------------------------------------
# Metadata / key resolution
# ---------------------------------------------------------------------------

def test_resolve_metadata_from_anndata_object():
    adata = _make_adata()
    meta = resolve_representation_metadata(adata, "rep_geneformer_cls_pca50")
    assert meta["representation_id"] == "rep_geneformer_cls_pca50"
    assert meta["reducer_fit_scope"] == "train_only"
    assert "scfm" in meta
    assert meta["scfm"]["model_key"] == "geneformer"
    assert meta["scfm"]["model_version"] == "V2"


def test_resolve_metadata_from_temp_h5ad(tmp_path):
    adata = _make_adata()
    p = tmp_path / "tiny.h5ad"
    adata.write_h5ad(p)
    back = ad.read_h5ad(p)
    meta = resolve_representation_metadata(back, "rep_geneformer_cls_pca50")
    assert meta["representation_id"] == "rep_geneformer_cls_pca50"
    assert meta["scfm"]["model_name"] == "Geneformer"


def test_resolve_keys_defaults():
    cfg = {
        "dataset": {"time_key": "abs_day"},
        "lineage": {"cell_state_key": "final_milestone_label_coarse"},
    }
    assert resolve_time_key(cfg, None) == "abs_day"
    assert resolve_state_key(cfg, None) == "final_milestone_label_coarse"
    assert resolve_time_key(cfg, "override_t") == "override_t"


def test_resolve_eval_timepoints_priority():
    observed = np.array([0.0, 1.0, 2.0, 3.0])
    # heldout wins
    cfg = {"scenario_params": {"train_times": [0.0, 1.0], "heldout_times": [2.0, 3.0]}}
    assert resolve_eval_timepoints(cfg, observed) == [2.0, 3.0]
    # fall back to observed-minus-train
    cfg2 = {"scenario_params": {"train_times": [0.0, 1.0]}}
    assert resolve_eval_timepoints(cfg2, observed) == [2.0, 3.0]
    # fall back to all observed
    cfg3 = {"scenario_params": {}}
    assert resolve_eval_timepoints(cfg3, observed) == [0.0, 1.0, 2.0, 3.0]


# ---------------------------------------------------------------------------
# End-to-end on a tiny synthetic config + result dir
# ---------------------------------------------------------------------------

def test_run_from_config_writes_both_outputs(tmp_path):
    adata = _make_adata()
    h5ad_path = tmp_path / "rep_input.h5ad"
    adata.write_h5ad(h5ad_path)

    config = {
        "dataset": {"id": "TINY", "h5ad_path": str(h5ad_path), "time_key": "abs_day"},
        "scenario_params": {"train_times": [0.0, 1.0], "heldout_times": [2.0, 3.0]},
        "cell_state_key": "final_milestone_label_coarse",
        "representation": {"enabled": True, "representation_id": "rep_geneformer_cls_pca50"},
        "evaluation": {"mode": "representation"},
    }
    config_path = tmp_path / "cfg.yaml"
    config_path.write_text(yaml.safe_dump(config))

    result_dir = tmp_path / "result"
    result_dir.mkdir()
    # 3D projected representation: (n_eval_tp, n_pred_cells, d)
    rng = np.random.default_rng(1)
    projected = rng.normal(size=(2, 8, 4)) + np.array([2.0, 3.0])[:, None, None]
    np.save(result_dir / "projected_expression.npy", projected)

    out = run_from_config(config_path, result_dir, backend="numpy")

    fc_path = result_dir / "representation_forecast_metrics.json"
    csv_path = result_dir / "representation_forecast_per_timepoint.csv"
    ts_path = result_dir / "representation_temporal_signal.json"
    assert fc_path.exists() and csv_path.exists() and ts_path.exists()

    fc = json.loads(fc_path.read_text())
    assert fc["metric_space"] == "representation"
    assert fc["status"] == "completed"
    for k in ("rep_wasserstein_distance", "rep_gaussian_mmd",
              "rep_energy_distance_mmd", "rep_hausdorff_loss"):
        assert k in fc
    # representation metadata embedded, including flat scfm block for reports
    assert fc["representation_metadata"]["scfm"]["model_version"] == "V2"

    ts = json.loads(ts_path.read_text())
    assert ts["metric_space"] == "representation"
    assert ts["representation_id"] == "rep_geneformer_cls_pca50"
    assert "temporal_signal_ratio" in ts
    assert out["forecast"]["status"] == "completed"


def test_run_from_config_missing_x_rep_raises(tmp_path):
    adata = _make_adata()
    del adata.obsm["X_rep"]
    h5ad_path = tmp_path / "no_xrep.h5ad"
    adata.write_h5ad(h5ad_path)
    config = {"dataset": {"h5ad_path": str(h5ad_path), "time_key": "abs_day"}}
    config_path = tmp_path / "cfg.yaml"
    config_path.write_text(yaml.safe_dump(config))
    result_dir = tmp_path / "result"
    result_dir.mkdir()
    np.save(result_dir / "projected_expression.npy", np.zeros((2, 4, 4)))
    with pytest.raises(KeyError):
        run_from_config(config_path, result_dir, backend="numpy")
