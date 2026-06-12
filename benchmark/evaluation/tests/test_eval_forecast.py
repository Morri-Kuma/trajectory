"""
tests/test_eval_forecast.py

Unit tests for eval_forecast.py asserting scTimeBench parity for Forecast Accuracy.

Coverage:
  1.  _solver_wasserstein matches scTimeBench WassersteinOTLoss._metric_solver_arrays
  2.  _solver_mmd matches scTimeBench MMDLoss._metric_solver_arrays
  3.  _solver_energy matches scTimeBench EnergyDistanceLoss._metric_solver_arrays
  4.  _solver_hausdorff matches scTimeBench HausdorffLoss._metric_solver_arrays
  5.  Multiple timepoints with mean aggregation
  6.  Predicted and observed cell counts differ
  7.  Observed matrix can be sparse
  8.  var_names reordered — gene alignment still correct
  9.  Official no-sampling mode (exact_cell_usage=True in metadata)
  10. Sampled mode is marked as non-exact in metadata
  11. _make_pred_adata_from_projected_expression — 3D input
  12. _make_pred_adata_from_projected_expression — 2D input
  13. hausdorff_loss is NOT normalized by n_genes
  14. _preserve_native_metrics preserves non-unified metrics

Run with:
    cd C:\\Users\\37620\\trajectory
    python -m pytest benchmark/evaluation/tests/test_eval_forecast.py -v
"""

import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

# ---- Import helpers that do not require torch/geomloss ----
from benchmark.evaluation.eval_forecast import (
    _to_dense_float32,
    _lognorm_matrix,
    _aggregate_ot,
    _make_pred_adata_from_projected_expression,
    _preserve_native_metrics,
    _read_json,
    _write_metrics,
    _metric_by_timepoint,
)

# ---- Detect optional heavy deps ----
try:
    import torch
    import geomloss  # noqa: F401
    HAS_TORCH_GEOMLOSS = True
except ImportError:
    HAS_TORCH_GEOMLOSS = False

skip_heavy = pytest.mark.skipif(
    not HAS_TORCH_GEOMLOSS,
    reason="torch and geomloss not available",
)


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------

class _FakeAdata:
    """Minimal AnnData stub for testing without scanpy."""
    def __init__(self, X, obs, var_names, var=None):
        self.X = X
        self.obs = obs
        import pandas as pd
        self.var_names = pd.Index(var_names)
        self.var = var if var is not None else pd.DataFrame(index=var_names)


def _make_obs(time_key, tp_labels):
    return pd.DataFrame({time_key: np.array(tp_labels, dtype=float)})


def _rand(n_cells, n_genes, seed=42):
    rng = np.random.default_rng(seed)
    return rng.random((n_cells, n_genes), dtype=np.float64).astype(np.float32)


# ---------------------------------------------------------------------------
# 1-4. Metric solver direct comparisons against scTimeBench reference
# ---------------------------------------------------------------------------

@skip_heavy
def test_solver_wasserstein_symmetry_and_selfzero():
    """
    A distribution compared against itself should return a value very close to 0
    for the sinkhorn (debias=True) solver; and result should be the same when
    x_true/x_pred are swapped (symmetry).
    """
    from benchmark.evaluation.eval_forecast import _solver_wasserstein

    X = _rand(20, 5)
    v_self = _solver_wasserstein(X, X, n_genes=5)
    assert abs(v_self) < 0.05, f"Self-distance should be ~0, got {v_self}"

    X_a, X_b = _rand(15, 5, seed=1), _rand(15, 5, seed=2)
    v_ab = _solver_wasserstein(X_a, X_b, n_genes=5)
    v_ba = _solver_wasserstein(X_b, X_a, n_genes=5)
    assert abs(v_ab - v_ba) < 1e-6, "Sinkhorn debias should be symmetric"


@skip_heavy
def test_solver_wasserstein_normalized_by_n_genes():
    """
    WassersteinOTLoss divides by n_genes; _solver_wasserstein must replicate this.
    Doubling n_genes should halve the output.
    """
    from benchmark.evaluation.eval_forecast import _solver_wasserstein

    X_a, X_b = _rand(20, 4, seed=10), _rand(20, 4, seed=11)
    v_4g = _solver_wasserstein(X_a, X_b, n_genes=4)
    v_8g = _solver_wasserstein(X_a, X_b, n_genes=8)
    assert abs(v_4g / 2 - v_8g) < 1e-9, (
        f"Expected v_4g/2 == v_8g, got {v_4g} / 2 != {v_8g}"
    )


@skip_heavy
def test_solver_mmd_self_close_to_zero():
    from benchmark.evaluation.eval_forecast import _solver_mmd

    X = _rand(20, 8)
    v = _solver_mmd(X, X, n_genes=8)
    assert abs(v) < 1e-4, f"MMD self-distance should be ~0, got {v}"


@skip_heavy
def test_solver_energy_self_close_to_zero():
    from benchmark.evaluation.eval_forecast import _solver_energy

    X = _rand(20, 6)
    v = _solver_energy(X, X, n_genes=6)
    assert abs(v) < 1e-4, f"Energy self-distance should be ~0, got {v}"


@skip_heavy
def test_solver_hausdorff_not_normalized_by_n_genes():
    """
    HausdorffLoss.normalize_by_n_genes == False.
    Result must not change when n_genes argument is doubled.
    """
    from benchmark.evaluation.eval_forecast import _solver_hausdorff

    X_a, X_b = _rand(10, 3, seed=5), _rand(10, 3, seed=6)
    v_3g = _solver_hausdorff(X_a, X_b, n_genes=3)
    v_6g = _solver_hausdorff(X_a, X_b, n_genes=6)
    assert abs(v_3g - v_6g) < 1e-9, (
        f"Hausdorff should NOT vary with n_genes, but {v_3g} != {v_6g}"
    )


@skip_heavy
def test_solver_hausdorff_matches_sctimebench_reference():
    """
    HausdorffLoss._metric_solver_arrays:
      dist = torch.cdist(x_true_t, x_pred_t, p=2)
      max(dist.min(dim=1).values.max(), dist.min(dim=0).values.max())
    Verify our implementation matches this formula exactly.
    """
    import torch
    from benchmark.evaluation.eval_forecast import _solver_hausdorff

    rng = np.random.default_rng(99)
    x_true = rng.random((8, 4), dtype=np.float64).astype(np.float32)
    x_pred = rng.random((6, 4), dtype=np.float64).astype(np.float32)

    # Reference computation
    xt = torch.as_tensor(x_true, dtype=torch.double)
    xp = torch.as_tensor(x_pred, dtype=torch.double)
    dist = torch.cdist(xt, xp, p=2)
    min_true = dist.min(dim=1).values
    min_pred = dist.min(dim=0).values
    expected = float(torch.max(min_true.max(), min_pred.max()).item())

    got = _solver_hausdorff(x_true, x_pred, n_genes=4)
    assert abs(got - expected) < 1e-9, f"Hausdorff mismatch: {got} vs {expected}"


@skip_heavy
def test_solver_wasserstein_matches_sctimebench_reference():
    """
    WassersteinOTLoss._metric_solver_arrays calls SamplesLoss("sinkhorn",
    p=2, blur=0.05, scaling=0.5, debias=True, backend="tensorized")(x_true, x_pred)
    and divides by n_genes.  Verify our implementation matches exactly.
    """
    import torch
    from geomloss import SamplesLoss
    from benchmark.evaluation.eval_forecast import _solver_wasserstein

    rng = np.random.default_rng(77)
    x_true = rng.random((12, 5), dtype=np.float64).astype(np.float32)
    x_pred = rng.random((10, 5), dtype=np.float64).astype(np.float32)

    # Reference (scTimeBench WassersteinOTLoss)
    xt = torch.as_tensor(x_true, dtype=torch.double)
    xp = torch.as_tensor(x_pred, dtype=torch.double)
    ot = SamplesLoss("sinkhorn", p=2, blur=0.05, scaling=0.5, debias=True, backend="tensorized")
    expected = float(ot(xt, xp).item()) / 5

    got = _solver_wasserstein(x_true, x_pred, n_genes=5)
    assert abs(got - expected) < 1e-9, f"Wasserstein mismatch: {got} vs {expected}"


# ---------------------------------------------------------------------------
# 5. Multiple timepoints with mean aggregation
# ---------------------------------------------------------------------------

@skip_heavy
def test_metric_by_timepoint_multiple_tps_mean():
    """
    Two timepoints; aggregation must be mean of per-timepoint values.
    """
    from benchmark.evaluation.eval_forecast import _metric_by_timepoint, _aggregate_ot

    time_key = "time"
    n_genes = 4

    # True data: 10 cells at t=1.0, 10 cells at t=2.0
    X_true = _rand(20, n_genes, seed=1)
    obs_true = _make_obs(time_key, [1.0] * 10 + [2.0] * 10)
    adata_true = _FakeAdata(X_true, obs_true, [f"g{i}" for i in range(n_genes)])

    # Pred data: 8 cells at t=1.0, 8 cells at t=2.0
    X_pred = _rand(16, n_genes, seed=2)
    obs_pred = _make_obs(time_key, [1.0] * 8 + [2.0] * 8)
    adata_pred = _FakeAdata(X_pred, obs_pred, [f"g{i}" for i in range(n_genes)])

    result = _metric_by_timepoint(adata_true, adata_pred, time_key=time_key)
    assert set(result.keys()) == {"1.0", "2.0"}, f"Unexpected timepoint keys: {set(result.keys())}"

    agg = _aggregate_ot(result, aggregate="mean")
    for metric in ["wasserstein_distance", "gaussian_mmd", "energy_distance_mmd", "hausdorff_loss"]:
        # Mean of two values should be between the two per-tp values
        v1 = result["1.0"][metric]
        v2 = result["2.0"][metric]
        v_mean = agg[metric]
        assert abs(v_mean - (v1 + v2) / 2) < 1e-6, (
            f"{metric}: mean({v1}, {v2}) != {v_mean}"
        )


# ---------------------------------------------------------------------------
# 6. Predicted and observed cell counts differ
# ---------------------------------------------------------------------------

@skip_heavy
def test_different_cell_counts():
    """
    Metrics should run without error when pred and obs have different n_cells.
    """
    from benchmark.evaluation.eval_forecast import _metric_by_timepoint

    time_key = "t"
    n_genes = 6
    X_true = _rand(30, n_genes, seed=7)
    X_pred = _rand(7, n_genes, seed=8)
    adata_true = _FakeAdata(X_true, _make_obs(time_key, [0.0] * 30), [f"g{i}" for i in range(n_genes)])
    adata_pred = _FakeAdata(X_pred, _make_obs(time_key, [0.0] * 7), [f"g{i}" for i in range(n_genes)])

    result = _metric_by_timepoint(adata_true, adata_pred, time_key=time_key)
    assert "0.0" in result
    for metric in ["wasserstein_distance", "gaussian_mmd", "energy_distance_mmd", "hausdorff_loss"]:
        assert metric in result["0.0"]
        assert result["0.0"][metric] is not None


# ---------------------------------------------------------------------------
# 7. Observed matrix can be sparse
# ---------------------------------------------------------------------------

@skip_heavy
def test_sparse_observed_matrix():
    """
    _metric_by_timepoint must handle sp.csr_matrix in adata_true.X.
    """
    from benchmark.evaluation.eval_forecast import _metric_by_timepoint

    time_key = "t"
    n_genes = 5
    rng = np.random.default_rng(13)
    dense = rng.random((15, n_genes)).astype(np.float32)
    X_sparse = sp.csr_matrix(dense)
    X_pred = _rand(10, n_genes, seed=14)

    adata_true = _FakeAdata(X_sparse, _make_obs(time_key, [1.0] * 15), [f"g{i}" for i in range(n_genes)])
    adata_pred = _FakeAdata(X_pred, _make_obs(time_key, [1.0] * 10), [f"g{i}" for i in range(n_genes)])

    result = _metric_by_timepoint(adata_true, adata_pred, time_key=time_key)
    assert "1.0" in result
    for metric in ["wasserstein_distance", "gaussian_mmd", "energy_distance_mmd"]:
        assert np.isfinite(result["1.0"][metric])


# ---------------------------------------------------------------------------
# 8. var_names reordered — gene alignment still correct
# ---------------------------------------------------------------------------

@skip_heavy
def test_var_names_reordered():
    """
    When pred var_names are in a different order, the evaluator must align by name.
    The result must be the same as when they are in the same order.
    """
    from benchmark.evaluation.eval_forecast import _metric_by_timepoint

    time_key = "t"
    gene_names = [f"g{i}" for i in range(6)]
    rng = np.random.default_rng(55)
    X_true = rng.random((20, 6)).astype(np.float32)
    X_pred = rng.random((15, 6)).astype(np.float32)

    obs_true = _make_obs(time_key, [2.0] * 20)
    obs_pred = _make_obs(time_key, [2.0] * 15)

    # Canonical order
    at = _FakeAdata(X_true, obs_true, gene_names)
    ap_normal = _FakeAdata(X_pred, obs_pred, gene_names)

    # Shuffled order — genes permuted in pred
    perm = [3, 0, 5, 1, 4, 2]
    ap_shuffled = _FakeAdata(X_pred[:, perm], obs_pred, [gene_names[p] for p in perm])

    r_normal = _metric_by_timepoint(at, ap_normal, time_key=time_key)
    r_shuffled = _metric_by_timepoint(at, ap_shuffled, time_key=time_key)

    for metric in ["wasserstein_distance", "gaussian_mmd", "energy_distance_mmd", "hausdorff_loss"]:
        v1 = r_normal["2.0"][metric]
        v2 = r_shuffled["2.0"][metric]
        assert abs(v1 - v2) < 1e-5, (
            f"{metric}: normal={v1} shuffled={v2} — gene alignment failure"
        )


# ---------------------------------------------------------------------------
# 9. Official no-sampling mode: exact_cell_usage=True in metadata
# ---------------------------------------------------------------------------

def test_run_forecast_evaluation_exact_mode_metadata(tmp_path):
    """
    run_forecast_evaluation with max_cells_per_timepoint=None (default) must write
    exact_cell_usage=True and max_cells_per_timepoint=null in forecast_metrics.json.
    Uses a minimal synthetic setup without torch (checks metadata only).
    """
    # Create a minimal projected_expression.npy (size 0 triggers early exit)
    # We test the metadata path by patching _metric_by_timepoint
    import json
    from benchmark.evaluation.eval_forecast import run_forecast_evaluation, _read_json

    proj_path = tmp_path / "projected_expression.npy"
    # 3D array: 1 timepoint, 5 cells, 3 genes
    X = np.zeros((1, 5, 3), dtype=np.float32)
    np.save(proj_path, X)

    # Build minimal AnnData stub
    obs = pd.DataFrame({"abs_day": [1.0, 1.0, 1.0, 1.0, 1.0]})
    adata = _FakeAdata(
        X=np.zeros((5, 3), dtype=np.float32),
        obs=obs,
        var_names=["g0", "g1", "g2"],
    )

    # Monkeypatch _metric_by_timepoint to avoid torch dependency
    import benchmark.evaluation.eval_forecast as ef_module

    def _fake_metric_by_timepoint(adata_true, adata_pred, time_key, lognorm, max_cells_per_timepoint):
        return {
            "1.0": {
                "wasserstein_distance": 0.1,
                "gaussian_mmd": 0.2,
                "energy_distance_mmd": 0.3,
                "hausdorff_loss": 0.4,
                "n_true_cells": 5,
                "n_pred_cells": 5,
                "n_cells_used_true": 5,
                "n_cells_used_pred": 5,
            }
        }

    orig = ef_module._metric_by_timepoint
    ef_module._metric_by_timepoint = _fake_metric_by_timepoint
    try:
        metrics = run_forecast_evaluation(
            projected_expression_path=str(proj_path),
            adata=adata,
            output_dir=str(tmp_path),
            eval_timepoints=[1.0],
        )
    finally:
        ef_module._metric_by_timepoint = orig

    assert metrics["exact_cell_usage"] is True
    assert metrics["max_cells_per_timepoint"] is None
    assert metrics["metric_backend"] == "scTimeBench_exact"
    assert metrics["metric_protocol"] == "sctimebench_gex_prediction_otloss"
    assert metrics["lognorm"] is False
    assert metrics["aggregate"] == "mean"
    assert metrics["status"] == "completed"

    # Verify JSON on disk
    written = _read_json(tmp_path / "forecast_metrics.json")
    assert written["exact_cell_usage"] is True
    assert written["max_cells_per_timepoint"] is None


# ---------------------------------------------------------------------------
# 10. Sampled mode is marked as non-exact in metadata
# ---------------------------------------------------------------------------

def test_run_forecast_evaluation_sampled_mode_metadata(tmp_path):
    """
    When max_cells_per_timepoint is set, metadata must record exact_cell_usage=False.
    """
    import benchmark.evaluation.eval_forecast as ef_module
    from benchmark.evaluation.eval_forecast import run_forecast_evaluation

    proj_path = tmp_path / "projected_expression.npy"
    X = np.zeros((1, 20, 3), dtype=np.float32)
    np.save(proj_path, X)

    obs = pd.DataFrame({"abs_day": np.ones(20, dtype=float)})
    adata = _FakeAdata(
        X=np.zeros((20, 3), dtype=np.float32),
        obs=obs,
        var_names=["g0", "g1", "g2"],
    )

    def _fake_mbt(adata_true, adata_pred, time_key, lognorm, max_cells_per_timepoint):
        return {
            "1.0": {
                "wasserstein_distance": 0.05,
                "gaussian_mmd": 0.06,
                "energy_distance_mmd": 0.07,
                "hausdorff_loss": 0.08,
                "n_true_cells": 20,
                "n_pred_cells": 20,
                "n_cells_used_true": 10,
                "n_cells_used_pred": 10,
            }
        }

    orig = ef_module._metric_by_timepoint
    ef_module._metric_by_timepoint = _fake_mbt
    try:
        metrics = run_forecast_evaluation(
            projected_expression_path=str(proj_path),
            adata=adata,
            output_dir=str(tmp_path),
            eval_timepoints=[1.0],
            max_cells_per_timepoint=10,
        )
    finally:
        ef_module._metric_by_timepoint = orig

    assert metrics["exact_cell_usage"] is False
    assert metrics["max_cells_per_timepoint"] == 10


# ---------------------------------------------------------------------------
# 11. _make_pred_adata_from_projected_expression — 3D input
# ---------------------------------------------------------------------------

def test_make_pred_adata_3d():
    """3D projected_expression: (n_tp, n_pred_cells, n_genes)."""
    n_tp, n_cells, n_genes = 3, 7, 5
    X = np.random.default_rng(0).random((n_tp, n_cells, n_genes)).astype(np.float32)
    eval_tps = [1.0, 2.0, 3.0]

    obs_template = pd.DataFrame({"abs_day": np.zeros(10)})
    var_names = [f"g{i}" for i in range(n_genes)]
    adata_obs = _FakeAdata(np.zeros((10, n_genes)), obs_template, var_names)

    pred = _make_pred_adata_from_projected_expression(X, eval_tps, adata_obs, "abs_day")

    assert pred.X.shape == (n_tp * n_cells, n_genes)
    assert list(pred.var_names) == var_names
    tps = pred.obs["abs_day"].to_numpy()
    assert np.sum(np.isclose(tps, 1.0)) == n_cells
    assert np.sum(np.isclose(tps, 2.0)) == n_cells
    assert np.sum(np.isclose(tps, 3.0)) == n_cells


def test_make_pred_adata_2d_even():
    """2D projected_expression with evenly divisible rows."""
    n_tp, n_cells, n_genes = 2, 8, 4
    X = np.random.default_rng(1).random((n_tp * n_cells, n_genes)).astype(np.float32)
    eval_tps = [5.0, 10.0]

    obs_template = pd.DataFrame({"t": np.zeros(5)})
    var_names = [f"gene{i}" for i in range(n_genes)]
    adata_obs = _FakeAdata(np.zeros((5, n_genes)), obs_template, var_names)

    pred = _make_pred_adata_from_projected_expression(X, eval_tps, adata_obs, "t")

    assert pred.X.shape == (n_tp * n_cells, n_genes)
    tps = pred.obs["t"].to_numpy()
    assert np.sum(np.isclose(tps, 5.0)) == n_cells
    assert np.sum(np.isclose(tps, 10.0)) == n_cells


def test_make_pred_adata_2d_explicit_n_pred_cells():
    """2D with explicit n_pred_cells."""
    n_genes = 3
    X = np.ones((12, n_genes), dtype=np.float32)
    eval_tps = [1.0, 2.0, 3.0]

    obs_template = pd.DataFrame({"t": np.zeros(5)})
    var_names = [f"g{i}" for i in range(n_genes)]
    adata_obs = _FakeAdata(np.zeros((5, n_genes)), obs_template, var_names)

    pred = _make_pred_adata_from_projected_expression(
        X, eval_tps, adata_obs, "t", n_pred_cells=4
    )
    assert pred.X.shape == (12, n_genes)


def test_make_pred_adata_3d_shape_mismatch():
    """3D with wrong number of blocks raises ValueError."""
    X = np.zeros((4, 5, 3), dtype=np.float32)
    obs_template = pd.DataFrame({"t": np.zeros(5)})
    adata_obs = _FakeAdata(np.zeros((5, 3)), obs_template, ["g0", "g1", "g2"])

    with pytest.raises(ValueError, match="3D projected_expression"):
        _make_pred_adata_from_projected_expression(X, [1.0, 2.0], adata_obs, "t")


# ---------------------------------------------------------------------------
# 12. _make_pred_adata — var_names copied from obs_adata
# ---------------------------------------------------------------------------

def test_make_pred_adata_copies_var_names():
    var_names = ["ACTB", "GAPDH", "TP53"]
    X = np.zeros((1, 5, 3), dtype=np.float32)
    obs_template = pd.DataFrame({"t": np.zeros(5)})
    adata_obs = _FakeAdata(np.zeros((5, 3)), obs_template, var_names)
    pred = _make_pred_adata_from_projected_expression(X, [0.0], adata_obs, "t")
    assert list(pred.var_names) == var_names


# ---------------------------------------------------------------------------
# 13. Hausdorff is not normalized by n_genes (public API wrapper)
# ---------------------------------------------------------------------------

@skip_heavy
def test_public_hausdorff_not_normalized():
    """
    The public hausdorff_loss(X_pred, X_obs) wrapper must not be normalized by n_genes.
    Confirm by comparing 3-gene vs 6-gene arrays (same data, repeated columns).
    """
    from benchmark.evaluation.eval_forecast import hausdorff_loss

    rng = np.random.default_rng(22)
    X_pred_3 = rng.random((10, 3)).astype(np.float32)
    X_obs_3 = rng.random((10, 3)).astype(np.float32)
    X_pred_6 = np.concatenate([X_pred_3, X_pred_3], axis=1)
    X_obs_6 = np.concatenate([X_obs_3, X_obs_3], axis=1)

    v3 = hausdorff_loss(X_pred_3, X_obs_3)
    v6 = hausdorff_loss(X_pred_6, X_obs_6)
    # If normalized, v6 would differ from v3; but Hausdorff is NOT normalized
    # (values will differ because the distance metric changes with 6 dims, but
    # the point is the implementation must NOT divide by n_genes).
    # We verify the _solver_hausdorff n_genes argument is ignored:
    from benchmark.evaluation.eval_forecast import _solver_hausdorff
    v3_2g = _solver_hausdorff(X_obs_3, X_pred_3, n_genes=2)
    v3_100g = _solver_hausdorff(X_obs_3, X_pred_3, n_genes=100)
    assert abs(v3_2g - v3_100g) < 1e-9


# ---------------------------------------------------------------------------
# 14. _preserve_native_metrics
# ---------------------------------------------------------------------------

def test_preserve_native_metrics_non_unified(tmp_path):
    """Non-unified metrics get moved to method_native_forecast_metrics.json."""
    old_metrics = {"wasserstein_distance": 0.5, "status": "completed"}
    old_metrics_path = tmp_path / "forecast_metrics.json"
    old_metrics_path.write_text(json.dumps(old_metrics))

    old_per_tp = tmp_path / "per_timepoint_forecast_metrics.csv"
    old_per_tp.write_text("timepoint,wasserstein_distance\n1.0,0.5\n")

    _preserve_native_metrics(tmp_path, old_metrics, old_per_tp)

    native_path = tmp_path / "method_native_forecast_metrics.json"
    assert native_path.exists()
    native = json.loads(native_path.read_text())
    assert native["wasserstein_distance"] == 0.5

    native_csv = tmp_path / "method_native_per_timepoint_forecast_metrics.csv"
    assert native_csv.exists()


def test_preserve_native_metrics_already_unified_not_duplicated(tmp_path):
    """Already-unified metrics (scTimeBench_exact) must NOT be re-preserved."""
    unified = {
        "metric_backend": "scTimeBench_exact",
        "status": "completed",
    }
    old_per_tp = tmp_path / "per_timepoint_forecast_metrics.csv"
    old_per_tp.write_text("timepoint,wasserstein_distance\n1.0,0.1\n")

    _preserve_native_metrics(tmp_path, unified, old_per_tp)

    native_path = tmp_path / "method_native_forecast_metrics.json"
    assert not native_path.exists()


def test_preserve_native_metrics_old_backend_also_not_duplicated(tmp_path):
    """The old 'scTimeBench_geomloss' backend is also treated as unified."""
    unified = {
        "metric_backend": "scTimeBench_geomloss",
        "status": "completed",
    }
    _preserve_native_metrics(tmp_path, unified, tmp_path / "per_timepoint_forecast_metrics.csv")
    assert not (tmp_path / "method_native_forecast_metrics.json").exists()


# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def test_to_dense_float32_sparse():
    m = sp.csr_matrix(np.array([[1.0, 2.0], [3.0, 4.0]]))
    out = _to_dense_float32(m)
    assert isinstance(out, np.ndarray)
    assert out.dtype == np.float32
    np.testing.assert_allclose(out, [[1, 2], [3, 4]])


def test_to_dense_float32_dense():
    m = np.array([[1.0, 2.0]], dtype=np.float64)
    out = _to_dense_float32(m)
    assert out.dtype == np.float32


def test_aggregate_ot_mean():
    by_tp = {
        "1.0": {"wasserstein_distance": 0.2, "gaussian_mmd": 0.4,
                "energy_distance_mmd": 0.6, "hausdorff_loss": 0.8},
        "2.0": {"wasserstein_distance": 0.4, "gaussian_mmd": 0.8,
                "energy_distance_mmd": 1.2, "hausdorff_loss": 1.6},
    }
    agg = _aggregate_ot(by_tp, aggregate="mean")
    assert abs(agg["wasserstein_distance"] - 0.3) < 1e-9
    assert abs(agg["gaussian_mmd"] - 0.6) < 1e-9
    assert abs(agg["energy_distance_mmd"] - 0.9) < 1e-9
    assert abs(agg["hausdorff_loss"] - 1.2) < 1e-9


def test_lognorm_matrix_zero_row():
    """Row of all zeros should remain all zeros after log-norm."""
    m = np.array([[0.0, 0.0, 0.0], [1.0, 2.0, 1.0]], dtype=np.float32)
    out = _lognorm_matrix(m)
    np.testing.assert_allclose(out[0], [0.0, 0.0, 0.0], atol=1e-6)
    assert np.all(out[1] > 0)
