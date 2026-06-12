"""
eval_forecast.py  --  Forecast Accuracy evaluator.
Framework reference: docs/framework/experimental_framework_v2.md §7, §14 Step 4.

Implements OT-based gene-expression prediction metrics exactly matching
scTimeBench OTLossMetric / WassersteinOTLoss / MMDLoss / EnergyDistanceLoss /
HausdorffLoss (base.py, wass.py, mmd.py, energy.py, hausdorff.py).

Alignment guarantees (must not drift from scTimeBench):
  - shared timepoints intersected between observed and predicted
  - shared genes by var_names intersected; columns aligned per timepoint
  - per-timepoint metric aggregated by mean (default)
  - lognorm=False by default
  - normalize_by_n_genes: True for Wasserstein/MMD/Energy, False for Hausdorff
  - no cell subsampling by default (exact_cell_usage=True in official mode)
  - solver argument order: (x_true, x_pred, n_genes) -- true first

STATUS: INACTIVE -- WOT and CellRank2 are NOT eligible (no unseen-timepoint
projection). The dispatcher enforces this via capability flags.
"""

import json
import shutil
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
from typing import Optional


# ------------------------------------------------------------------
# Metric solvers -- matching scTimeBench _metric_solver_arrays convention
# ------------------------------------------------------------------

def _solver_wasserstein(x_true: np.ndarray, x_pred: np.ndarray, n_genes: int) -> float:
    """WassersteinOTLoss: sinkhorn p=2 blur=0.05 scaling=0.5 debias / n_genes."""
    import torch
    from geomloss import SamplesLoss
    xt = torch.as_tensor(x_true, dtype=torch.double)
    xp = torch.as_tensor(x_pred, dtype=torch.double)
    ot = SamplesLoss("sinkhorn", p=2, blur=0.05, scaling=0.5,
                     debias=True, backend="tensorized")
    return float(ot(xt, xp).item() / n_genes)


def _solver_mmd(x_true: np.ndarray, x_pred: np.ndarray, n_genes: int) -> float:
    """MMDLoss: gaussian blur=1.0 debias / n_genes."""
    import torch
    from geomloss import SamplesLoss
    xt = torch.as_tensor(x_true, dtype=torch.double)
    xp = torch.as_tensor(x_pred, dtype=torch.double)
    mmd = SamplesLoss("gaussian", blur=1.0, debias=True, backend="tensorized")
    return float(mmd(xt, xp).item() / n_genes)


def _solver_energy(x_true: np.ndarray, x_pred: np.ndarray, n_genes: int) -> float:
    """EnergyDistanceLoss: energy blur=1.0 debias / n_genes."""
    import torch
    from geomloss import SamplesLoss
    xt = torch.as_tensor(x_true, dtype=torch.double)
    xp = torch.as_tensor(x_pred, dtype=torch.double)
    energy = SamplesLoss("energy", blur=1.0, debias=True, backend="tensorized")
    return float(energy(xt, xp).item() / n_genes)


def _solver_hausdorff(x_true: np.ndarray, x_pred: np.ndarray, n_genes: int) -> float:
    """HausdorffLoss: cdist bidirectional, NOT normalized by n_genes."""
    import torch
    if x_true.size == 0 or x_pred.size == 0:
        return float("nan")
    xt = torch.as_tensor(x_true, dtype=torch.double)
    xp = torch.as_tensor(x_pred, dtype=torch.double)
    dist = torch.cdist(xt, xp, p=2)
    if dist.numel() == 0:
        return float("nan")
    min_true = dist.min(dim=1).values  # each true cell: min dist to any pred
    min_pred = dist.min(dim=0).values  # each pred cell: min dist to any true
    return float(torch.max(min_true.max(), min_pred.max()).item())


# ------------------------------------------------------------------
# Public metric wrappers (API: X_pred first, X_obs second)
# Internally reorder to (x_true=X_obs, x_pred=X_pred) per scTimeBench.
# ------------------------------------------------------------------

def wasserstein_distance(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    """scTimeBench WassersteinOTLoss."""
    n_genes = X_pred.shape[1] if X_pred.ndim == 2 else 1
    return _solver_wasserstein(
        np.asarray(X_obs, dtype=np.float32),
        np.asarray(X_pred, dtype=np.float32),
        n_genes,
    )


def gaussian_mmd(X_pred: np.ndarray, X_obs: np.ndarray,
                 sigma: Optional[float] = None) -> float:
    """scTimeBench MMDLoss. sigma is ignored (kept for API compat)."""
    _ = sigma
    n_genes = X_pred.shape[1] if X_pred.ndim == 2 else 1
    return _solver_mmd(
        np.asarray(X_obs, dtype=np.float32),
        np.asarray(X_pred, dtype=np.float32),
        n_genes,
    )


def energy_distance_mmd(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    """scTimeBench EnergyDistanceLoss."""
    n_genes = X_pred.shape[1] if X_pred.ndim == 2 else 1
    return _solver_energy(
        np.asarray(X_obs, dtype=np.float32),
        np.asarray(X_pred, dtype=np.float32),
        n_genes,
    )


def hausdorff_loss(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    """scTimeBench HausdorffLoss. Not normalized by n_genes."""
    return _solver_hausdorff(
        np.asarray(X_obs, dtype=np.float32),
        np.asarray(X_pred, dtype=np.float32),
        n_genes=0,  # ignored: normalize_by_n_genes=False
    )


# ------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------

def _to_dense_float32(X) -> np.ndarray:
    if sp.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=np.float32)


def _sample_rows(X: np.ndarray, max_cells: int, seed: int) -> np.ndarray:
    if X.shape[0] <= max_cells:
        return X
    rng = np.random.default_rng(seed)
    idx = rng.choice(X.shape[0], size=max_cells, replace=False)
    return X[np.sort(idx)]


def _lognorm_matrix(matrix: np.ndarray) -> np.ndarray:
    """Log-normalise dense float32 matrix (matches scTimeBench base.py)."""
    totals = np.sum(matrix, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.where(totals == 0, 0.0, 1e4 / totals)
    normed = matrix * scale[:, None]
    return np.log1p(normed)


def _aggregate_ot(metric_by_tp: dict, aggregate: str = "mean") -> dict:
    """Aggregate per-timepoint metrics. Matches OTLossMetric._aggregate_ot."""
    metric_keys = ["wasserstein_distance", "gaussian_mmd",
                   "energy_distance_mmd", "hausdorff_loss"]
    agg = {}
    for key in metric_keys:
        vals = [v[key] for v in metric_by_tp.values()
                if isinstance(v, dict) and key in v and v[key] is not None]
        arr = np.array(vals, dtype=float)
        if len(arr) == 0:
            agg[key] = None
        elif aggregate == "mean":
            agg[key] = float(np.mean(arr))
        elif aggregate == "median":
            agg[key] = float(np.median(arr))
        elif aggregate == "sum":
            agg[key] = float(np.sum(arr))
        else:
            raise ValueError(f"Invalid aggregate: {aggregate!r}")
    return agg


class _PredAdata:
    """Minimal AnnData-like container for predicted expression."""
    def __init__(self, X, obs, var_names, var):
        self.X = X
        self.obs = obs
        self.var_names = var_names
        self.var = var


def _make_pred_adata_from_projected_expression(
    X_pred: np.ndarray,
    eval_timepoints: list,
    obs_adata,
    time_key: str,
    n_pred_cells: Optional[int] = None,
) -> _PredAdata:
    """
    Convert projected_expression.npy into a _PredAdata object.

    Supports:
      - 3D: (n_eval_timepoints, n_pred_cells, n_genes)
      - 2D: (n_eval_timepoints * n_pred_cells, n_genes) concatenated blocks

    var_names and var are copied from obs_adata.
    obs[time_key] is set to the corresponding eval_timepoint for each row.
    """
    n_tp = len(eval_timepoints)
    if X_pred.ndim == 3:
        if X_pred.shape[0] != n_tp:
            raise ValueError(
                f"3D projected_expression has {X_pred.shape[0]} blocks "
                f"but {n_tp} eval_timepoints were provided."
            )
        blocks = [_to_dense_float32(X_pred[i]) for i in range(n_tp)]
        cells_per_tp = [b.shape[0] for b in blocks]
    elif X_pred.ndim == 2:
        if n_pred_cells is None:
            if n_tp == 0:
                raise ValueError("eval_timepoints must be non-empty for 2D input.")
            if X_pred.shape[0] % n_tp != 0:
                raise ValueError(
                    f"Cannot infer n_pred_cells: {X_pred.shape[0]} rows "
                    f"for {n_tp} timepoints."
                )
            n_pred_cells = X_pred.shape[0] // n_tp
        blocks = [
            _to_dense_float32(X_pred[i * n_pred_cells:(i + 1) * n_pred_cells])
            for i in range(n_tp)
        ]
        cells_per_tp = [n_pred_cells] * n_tp
    else:
        raise ValueError(
            f"Unsupported shape {X_pred.shape}; expected 2D or 3D."
        )
    X_all = np.concatenate(blocks, axis=0)
    tp_labels = np.concatenate(
        [np.full(cells_per_tp[i], float(eval_timepoints[i])) for i in range(n_tp)]
    )
    obs_df = pd.DataFrame({time_key: tp_labels})
    obs_df.index = [str(k) for k in range(len(tp_labels))]
    return _PredAdata(
        X=X_all,
        obs=obs_df,
        var_names=obs_adata.var_names.copy(),
        var=obs_adata.var.copy(),
    )


def _metric_by_timepoint(
    adata_true,
    adata_pred,
    time_key: str,
    lognorm: bool = False,
    max_cells_per_timepoint: Optional[int] = None,
) -> dict:
    """
    Per-timepoint metric computation matching OTLossMetric._metric_by_timepoint.

    Steps:
      1. Intersect timepoints (np.intersect1d on obs[time_key]).
      2. Intersect var_names; compute aligned column indices.
      3. For each shared timepoint: extract dense matrices, optional lognorm,
         call all four solvers with (x_true, x_pred, n_genes).
      4. Return dict keyed str(tp) -> per-metric dict.

    max_cells_per_timepoint=None means exact mode (all cells, no sampling).
    """
    true_tps = np.unique(adata_true.obs[time_key].astype(float))
    pred_tps = np.unique(adata_pred.obs[time_key].astype(float))
    shared_tps = np.intersect1d(true_tps, pred_tps)
    if len(shared_tps) == 0:
        raise ValueError(
            "No overlapping timepoints between observed and predicted data."
        )
    true_var_names = np.asarray(adata_true.var_names)
    pred_var_names = np.asarray(adata_pred.var_names)
    shared_genes = np.intersect1d(true_var_names, pred_var_names)
    if len(shared_genes) == 0:
        raise ValueError("No overlapping genes between observed and predicted.")
    true_gene_lookup = {gene: i for i, gene in enumerate(true_var_names)}
    pred_gene_lookup = {gene: i for i, gene in enumerate(pred_var_names)}
    true_gene_idx = np.array([true_gene_lookup[gene] for gene in shared_genes])
    pred_gene_idx = np.array([pred_gene_lookup[gene] for gene in shared_genes])
    n_genes = len(shared_genes)

    results = {}
    for tp in shared_tps:
        true_mask = np.isclose(adata_true.obs[time_key].astype(float).to_numpy(), tp)
        pred_mask = np.isclose(adata_pred.obs[time_key].astype(float).to_numpy(), tp)
        if not np.any(true_mask) or not np.any(pred_mask):
            continue
        x_true = _to_dense_float32(adata_true.X[true_mask][:, true_gene_idx])
        x_pred = _to_dense_float32(adata_pred.X[pred_mask][:, pred_gene_idx])
        if lognorm:
            x_true = _lognorm_matrix(x_true)
            x_pred = _lognorm_matrix(x_pred)
        n_true_orig = int(x_true.shape[0])
        n_pred_orig = int(x_pred.shape[0])
        if max_cells_per_timepoint is not None:
            tp_idx = int(np.searchsorted(shared_tps, tp))
            x_true = _sample_rows(x_true, max_cells_per_timepoint, 2000 + tp_idx)
            x_pred = _sample_rows(x_pred, max_cells_per_timepoint, 1000 + tp_idx)
        results[str(tp)] = {
            "wasserstein_distance": _solver_wasserstein(x_true, x_pred, n_genes),
            "gaussian_mmd": _solver_mmd(x_true, x_pred, n_genes),
            "energy_distance_mmd": _solver_energy(x_true, x_pred, n_genes),
            "hausdorff_loss": _solver_hausdorff(x_true, x_pred, n_genes),
            "n_true_cells": n_true_orig,
            "n_pred_cells": n_pred_orig,
            "n_cells_used_true": int(x_true.shape[0]),
            "n_cells_used_pred": int(x_pred.shape[0]),
        }
    return results


# ------------------------------------------------------------------
# Main evaluation entry point
# ------------------------------------------------------------------

def run_forecast_evaluation(
    projected_expression_path: str,
    adata,
    output_dir: str,
    held_out_time_labels: Optional[list] = None,
    time_key: str = "abs_day",
    eval_timepoints: Optional[list] = None,
    n_pred_cells: Optional[int] = None,
    max_cells_per_timepoint: Optional[int] = None,
    lognorm: bool = False,
    aggregate: str = "mean",
) -> dict:
    """
    Evaluate Forecast Accuracy for one method x scenario run.

    Reproduces OTLossMetric._gex_eval:
      1. Load projected_expression.npy -> _PredAdata.
      2. Intersect timepoints and genes.
      3. Per-timepoint metrics using all cells (exact mode by default).
      4. Mean aggregation across timepoints.
      5. Write forecast_metrics.json + per_timepoint_forecast_metrics.csv.

    Official default: max_cells_per_timepoint=None (no sampling, exact mode).
    Set max_cells_per_timepoint to a positive int only for debug/low-memory;
    the output will record exact_cell_usage=False in that case.

    IMPORTANT: Only call for methods with supports_unseen_timepoint_projection=True.
    The dispatcher enforces this. Do not call for WOT or CellRank2.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    proj_path = Path(projected_expression_path)
    if not proj_path.exists():
        m = {
            "wasserstein_distance": None, "gaussian_mmd": None,
            "energy_distance_mmd": None, "hausdorff_loss": None,
            "status": "error: projected_expression.npy not found",
        }
        _write_metrics(m, out_dir)
        return m

    X_pred_raw = np.load(proj_path, allow_pickle=True)

    if X_pred_raw.size == 0:
        m = {
            "wasserstein_distance": None, "gaussian_mmd": None,
            "energy_distance_mmd": None, "hausdorff_loss": None,
            "status": "inactive: no projected expression data.",
        }
        _write_metrics(m, out_dir)
        return m

    existing_metrics = _read_json(out_dir / "forecast_metrics.json")
    existing_per_tp = out_dir / "per_timepoint_forecast_metrics.csv"
    _preserve_native_metrics(out_dir, existing_metrics, existing_per_tp)

    if eval_timepoints is None:
        eval_timepoints = held_out_time_labels or existing_metrics.get("eval_timepoints")
    if eval_timepoints is None:
        raise ValueError(
            "eval_timepoints must be provided or present in forecast_metrics.json"
        )
    eval_timepoints = [float(t) for t in eval_timepoints]

    adata_pred = _make_pred_adata_from_projected_expression(
        X_pred=X_pred_raw,
        eval_timepoints=eval_timepoints,
        obs_adata=adata,
        time_key=time_key,
        n_pred_cells=n_pred_cells,
    )

    metric_by_tp = _metric_by_timepoint(
        adata_true=adata,
        adata_pred=adata_pred,
        time_key=time_key,
        lognorm=lognorm,
        max_cells_per_timepoint=max_cells_per_timepoint,
    )

    agg = _aggregate_ot(metric_by_tp, aggregate=aggregate)

    per_tp_rows = [
        {"timepoint": float(tp_str), **tp_m}
        for tp_str, tp_m in metric_by_tp.items()
    ]
    pd.DataFrame(per_tp_rows).to_csv(
        out_dir / "per_timepoint_forecast_metrics.csv", index=False
    )

    exact_mode = (max_cells_per_timepoint is None)
    metrics = {
        "wasserstein_distance": agg["wasserstein_distance"],
        "gaussian_mmd": agg["gaussian_mmd"],
        "energy_distance_mmd": agg["energy_distance_mmd"],
        "hausdorff_loss": agg["hausdorff_loss"],
        "n_eval_timepoints": len(eval_timepoints),
        "eval_timepoints": eval_timepoints,
        "scenario_type": existing_metrics.get("scenario_type"),
        "metric_backend": "scTimeBench_exact",
        "metric_protocol": "sctimebench_gex_prediction_otloss",
        "lognorm": lognorm,
        "aggregate": aggregate,
        "exact_cell_usage": exact_mode,
        "max_cells_per_timepoint": max_cells_per_timepoint,
        "metric_definitions": {
            "wasserstein_distance": {
                "class": "WassersteinOTLoss",
                "loss": "sinkhorn", "p": 2, "blur": 0.05,
                "scaling": 0.5, "debias": True, "backend": "tensorized",
                "normalize_by_n_genes": True,
            },
            "gaussian_mmd": {
                "class": "MMDLoss",
                "loss": "gaussian", "blur": 1.0,
                "debias": True, "backend": "tensorized",
                "normalize_by_n_genes": True,
            },
            "energy_distance_mmd": {
                "class": "EnergyDistanceLoss",
                "loss": "energy", "blur": 1.0,
                "debias": True, "backend": "tensorized",
                "normalize_by_n_genes": True,
            },
            "hausdorff_loss": {
                "class": "HausdorffLoss",
                "implementation": "torch.cdist(x_true,x_pred,p=2)_bidirectional",
                "normalize_by_n_genes": False,
            },
        },
        "status": "completed",
    }
    _write_metrics(metrics, out_dir)
    return metrics


# ------------------------------------------------------------------
# I/O helpers
# ------------------------------------------------------------------

def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        return json.load(f) or {}


def _write_metrics(metrics: dict, out_dir: Path) -> None:
    p = out_dir / "forecast_metrics.json"
    with open(p, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"[eval_forecast] Metrics written to {p}")


def _preserve_native_metrics(
    out_dir: Path,
    existing_metrics: dict,
    existing_per_tp: Path,
) -> None:
    """
    If existing forecast_metrics.json was written by a method adapter (not this
    unified evaluator), copy it to method_native_forecast_metrics.json.
    Runs with metric_backend in {scTimeBench_exact, scTimeBench_geomloss} are
    already unified and are not re-preserved.
    """
    if not existing_metrics:
        return
    if existing_metrics.get("metric_backend") in (
        "scTimeBench_exact", "scTimeBench_geomloss"
    ):
        return
    native = out_dir / "method_native_forecast_metrics.json"
    with open(native, "w", encoding="utf-8") as f:
        json.dump(existing_metrics, f, indent=2)
    if existing_per_tp.exists() and existing_per_tp.stat().st_size > 0:
        native_csv = out_dir / "method_native_per_timepoint_forecast_metrics.csv"
        shutil.copy2(existing_per_tp, native_csv)
