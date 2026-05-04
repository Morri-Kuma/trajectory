"""
eval_forecast.py
Forecast Accuracy evaluator — framework module (currently inactive).
Framework reference: experimental framework v2.md §7, §14 Step 4

STATUS: INACTIVE at the current benchmark stage.

WOT and CellRank2 are NOT eligible for Forecast Accuracy because they do not
generate projected cells at unseen future time points. This evaluator will be
activated only when generative / forecasting models are added to the benchmark.

This module is implemented as a clean framework stub so that the evaluation
logic is ready when eligible models arrive. It must NOT be called for WOT or
CellRank2. The dispatcher (eval_dispatch.py) enforces this via capability flags.

Metrics (per v2 §7.3):
  - Wasserstein Distance
  - Gaussian MMD
  - Energy Distance MMD
  - Hausdorff Loss

Aggregation (per v2 §7.4):
  Rank methods within each metric, then average ranks for the
  Forecast Accuracy rank within each scenario.
"""

import json
import shutil
import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path
from typing import Optional


# ------------------------------------------------------------------
# Metric functions
# ------------------------------------------------------------------

def wasserstein_distance(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    """
    scTimeBench Wasserstein OT loss.

    Matches scTimeBench's OTLossMetric/WassersteinOTLoss definition:
    geomloss SamplesLoss("sinkhorn", p=2, blur=0.05, scaling=0.5,
    debias=True, backend="tensorized"), normalized by n_genes.
    """
    return _geomloss_metric(
        X_pred,
        X_obs,
        loss="sinkhorn",
        normalize_by_n_genes=True,
        p=2,
        blur=0.05,
        scaling=0.5,
        debias=True,
        backend="tensorized",
    )


def gaussian_mmd(X_pred: np.ndarray, X_obs: np.ndarray,
                  sigma: Optional[float] = None) -> float:
    """
    scTimeBench Gaussian MMD loss.

    The sigma argument is retained for API compatibility but ignored; scTimeBench
    uses GeomLoss with blur=1.0.
    """
    _ = sigma
    return _geomloss_metric(
        X_pred,
        X_obs,
        loss="gaussian",
        normalize_by_n_genes=True,
        blur=1.0,
        debias=True,
        backend="tensorized",
    )


def energy_distance_mmd(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    """
    scTimeBench Energy Distance loss.
    """
    return _geomloss_metric(
        X_pred,
        X_obs,
        loss="energy",
        normalize_by_n_genes=True,
        blur=1.0,
        debias=True,
        backend="tensorized",
    )


def hausdorff_loss(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    """
    scTimeBench Hausdorff loss.

    Bidirectional nearest-neighbor Hausdorff distance computed with torch.cdist.
    Not normalized by n_genes by default, matching scTimeBench.
    """
    if X_pred.size == 0 or X_obs.size == 0:
        return float("nan")
    import torch

    pred = torch.as_tensor(np.asarray(X_pred), dtype=torch.double)
    obs = torch.as_tensor(np.asarray(X_obs), dtype=torch.double)
    d_xy = torch.cdist(pred, obs)
    return float(
        torch.maximum(
            d_xy.min(dim=1).values.max(),
            d_xy.min(dim=0).values.max(),
        ).item()
    )


def _geomloss_metric(
    X_pred: np.ndarray,
    X_obs: np.ndarray,
    loss: str,
    normalize_by_n_genes: bool,
    **samples_loss_kwargs,
) -> float:
    if X_pred.size == 0 or X_obs.size == 0:
        return float("nan")
    import torch
    from geomloss import SamplesLoss

    pred = torch.as_tensor(np.asarray(X_pred), dtype=torch.double)
    obs = torch.as_tensor(np.asarray(X_obs), dtype=torch.double)
    metric = SamplesLoss(loss, **samples_loss_kwargs)
    value = metric(pred, obs)
    if normalize_by_n_genes:
        value = value / pred.shape[1]
    return float(value.item())


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


# ------------------------------------------------------------------
# Main evaluation function
# ------------------------------------------------------------------

def run_forecast_evaluation(
    projected_expression_path: str,
    adata,
    output_dir: str,
    held_out_time_labels: Optional[list] = None,
    time_key: str = "abs_day",
    eval_timepoints: Optional[list] = None,
    n_pred_cells: Optional[int] = None,
    max_cells_per_timepoint: int = 1000,
) -> dict:
    """
    Evaluate Forecast Accuracy for one method × scenario run.

    IMPORTANT: This function must only be called for methods where
    supports_unseen_timepoint_projection = True.
    The dispatcher (eval_dispatch.py) enforces this. Do not call this
    function directly for WOT or CellRank2.

    Parameters
    ----------
    projected_expression_path : str
        Path to projected_expression.npy produced by the method adapter.
        Shape: (n_cells_projected, n_genes).
    adata : AnnData
        Benchmark input with observed cells at the evaluation time points.
    output_dir : str
        Directory where forecast_metrics.json and
        per_timepoint_forecast_metrics.csv will be written.
    held_out_time_labels : list, optional
        Time points that were held out during training and used as evaluation
        targets. If None, inferred from the scenario config.

    Returns
    -------
    dict with keys: wasserstein_distance, gaussian_mmd, energy_distance_mmd,
                    hausdorff_loss, forecast_accuracy_rank, status
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    proj_path = Path(projected_expression_path)
    if not proj_path.exists():
        metrics = {
            "wasserstein_distance": None,
            "gaussian_mmd": None,
            "energy_distance_mmd": None,
            "hausdorff_loss": None,
            "status": "error: projected_expression.npy not found",
        }
        _write_metrics(metrics, out_dir)
        return metrics

    X_pred = np.load(proj_path, allow_pickle=True)

    if X_pred.size == 0:
        metrics = {
            "wasserstein_distance": None,
            "gaussian_mmd": None,
            "energy_distance_mmd": None,
            "hausdorff_loss": None,
            "status": "inactive: no projected expression data. "
                       "This evaluator is inactive for the current benchmark stage.",
        }
        _write_metrics(metrics, out_dir)
        return metrics

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

    if X_pred.ndim == 3:
        if X_pred.shape[0] != len(eval_timepoints):
            raise ValueError(
                f"Projected expression has {X_pred.shape[0]} timepoint blocks "
                f"but {len(eval_timepoints)} eval_timepoints were provided."
            )
        n_pred_cells = int(X_pred.shape[1])
    elif X_pred.ndim == 2:
        if n_pred_cells is None:
            if len(eval_timepoints) == 0:
                raise ValueError("No eval_timepoints available.")
            if X_pred.shape[0] % len(eval_timepoints) != 0:
                raise ValueError(
                    f"Cannot infer n_pred_cells: {X_pred.shape[0]} projected rows "
                    f"for {len(eval_timepoints)} timepoints."
                )
            n_pred_cells = X_pred.shape[0] // len(eval_timepoints)
    else:
        raise ValueError(
            f"Unsupported projected expression shape {X_pred.shape}; expected 2D or 3D."
        )

    per_tp_existing = (
        pd.read_csv(existing_per_tp)
        if existing_per_tp.exists() and existing_per_tp.stat().st_size > 0
        else pd.DataFrame()
    )
    per_tp_rows = []
    cell_tps = adata.obs[time_key].astype(float).to_numpy()

    for i, tp in enumerate(eval_timepoints):
        if X_pred.ndim == 3:
            X_pred_tp = _to_dense_float32(X_pred[i])
        else:
            start = i * n_pred_cells
            stop = start + n_pred_cells
            X_pred_tp = _to_dense_float32(X_pred[start:stop])
        obs_mask = np.isclose(cell_tps, tp)
        X_obs_tp = _to_dense_float32(adata[obs_mask].X)
        X_pred_sample = _sample_rows(
            X_pred_tp, max_cells_per_timepoint, seed=1000 + i
        )
        X_obs_sample = _sample_rows(
            X_obs_tp, max_cells_per_timepoint, seed=2000 + i
        )
        row = {"timepoint": tp}
        if not per_tp_existing.empty and "timepoint" in per_tp_existing.columns:
            match = per_tp_existing[np.isclose(per_tp_existing["timepoint"].astype(float), tp)]
            if not match.empty:
                row.update(match.iloc[0].to_dict())
        row.update({
            "wasserstein_distance": wasserstein_distance(X_pred_sample, X_obs_sample),
            "gaussian_mmd": gaussian_mmd(X_pred_sample, X_obs_sample),
            "energy_distance_mmd": energy_distance_mmd(X_pred_sample, X_obs_sample),
            "hausdorff_loss": hausdorff_loss(X_pred_sample, X_obs_sample),
            "metric_sample_pred_cells": int(X_pred_sample.shape[0]),
            "metric_sample_obs_cells": int(X_obs_sample.shape[0]),
        })
        per_tp_rows.append(row)

    per_tp = pd.DataFrame(per_tp_rows)
    per_tp.to_csv(out_dir / "per_timepoint_forecast_metrics.csv", index=False)

    def _mean_col(col):
        return float(np.nanmean(per_tp[col].astype(float))) if col in per_tp else None

    metrics = {
        "wasserstein_distance": _mean_col("wasserstein_distance"),
        "gaussian_mmd": _mean_col("gaussian_mmd"),
        "energy_distance_mmd": _mean_col("energy_distance_mmd"),
        "hausdorff_loss": _mean_col("hausdorff_loss"),
        "n_eval_timepoints": len(eval_timepoints),
        "eval_timepoints": eval_timepoints,
        "scenario_type": existing_metrics.get("scenario_type"),
        "max_cells_per_timepoint": max_cells_per_timepoint,
        "metric_backend": "scTimeBench_geomloss",
        "metric_definitions": {
            "wasserstein_distance": {
                "loss": "sinkhorn",
                "p": 2,
                "blur": 0.05,
                "scaling": 0.5,
                "debias": True,
                "backend": "tensorized",
                "normalize_by_n_genes": True,
            },
            "gaussian_mmd": {
                "loss": "gaussian",
                "blur": 1.0,
                "debias": True,
                "backend": "tensorized",
                "normalize_by_n_genes": True,
            },
            "energy_distance_mmd": {
                "loss": "energy",
                "blur": 1.0,
                "debias": True,
                "backend": "tensorized",
                "normalize_by_n_genes": True,
            },
            "hausdorff_loss": {
                "implementation": "torch.cdist_bidirectional_hausdorff",
                "normalize_by_n_genes": False,
            },
        },
        "aggregation": "mean_over_eval_timepoints",
        "status": "completed",
    }
    _write_metrics(metrics, out_dir)
    return metrics


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        return json.load(f) or {}


def _write_metrics(metrics: dict, out_dir: Path):
    metrics_path = out_dir / "forecast_metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"[eval_forecast] Metrics written to {metrics_path}")


def _preserve_native_metrics(
    out_dir: Path,
    existing_metrics: dict,
    existing_per_tp: Path,
) -> None:
    if not existing_metrics:
        return
    if existing_metrics.get("metric_backend") == "scTimeBench_geomloss":
        return

    native_metrics_path = out_dir / "method_native_forecast_metrics.json"
    with open(native_metrics_path, "w", encoding="utf-8") as f:
        json.dump(existing_metrics, f, indent=2)

    if existing_per_tp.exists() and existing_per_tp.stat().st_size > 0:
        native_per_tp_path = out_dir / "method_native_per_timepoint_forecast_metrics.csv"
        shutil.copy2(existing_per_tp, native_per_tp_path)
