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
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional


# ------------------------------------------------------------------
# Metric functions
# ------------------------------------------------------------------

def wasserstein_distance(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    """
    Compute Wasserstein-1 distance between predicted and observed distributions.
    Uses a sliced approximation for high-dimensional gene expression data.
    """
    # TODO: implement sliced Wasserstein or use scipy.stats.wasserstein_distance
    # for each gene then aggregate.
    raise NotImplementedError("Wasserstein distance not yet implemented.")


def gaussian_mmd(X_pred: np.ndarray, X_obs: np.ndarray,
                  sigma: float = 1.0) -> float:
    """
    Compute Gaussian Maximum Mean Discrepancy between predicted and observed cells.
    """
    # TODO: implement Gaussian kernel MMD.
    raise NotImplementedError("Gaussian MMD not yet implemented.")


def energy_distance_mmd(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    """
    Compute Energy Distance MMD between predicted and observed cells.
    """
    # TODO: implement energy distance.
    raise NotImplementedError("Energy Distance MMD not yet implemented.")


def hausdorff_loss(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    """
    Compute Hausdorff distance between predicted and observed cell sets.
    """
    # TODO: implement Hausdorff distance.
    raise NotImplementedError("Hausdorff loss not yet implemented.")


# ------------------------------------------------------------------
# Main evaluation function
# ------------------------------------------------------------------

def run_forecast_evaluation(
    projected_expression_path: str,
    adata,
    output_dir: str,
    held_out_time_labels: Optional[list] = None,
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

    # TODO: extract observed cells at the held-out time points from adata.
    # X_obs = adata[adata.obs["time_label"].isin(held_out_time_labels)].X

    # TODO: compute all four metrics.
    # wd = wasserstein_distance(X_pred, X_obs)
    # mmd_g = gaussian_mmd(X_pred, X_obs)
    # mmd_e = energy_distance_mmd(X_pred, X_obs)
    # hl = hausdorff_loss(X_pred, X_obs)

    metrics = {
        "wasserstein_distance": None,
        "gaussian_mmd": None,
        "energy_distance_mmd": None,
        "hausdorff_loss": None,
        "status": "not_implemented: metric computation not yet wired.",
    }
    _write_metrics(metrics, out_dir)
    return metrics


def _write_metrics(metrics: dict, out_dir: Path):
    metrics_path = out_dir / "forecast_metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    per_tp = pd.DataFrame(columns=["timepoint", "wasserstein_distance",
                                    "gaussian_mmd", "energy_distance_mmd",
                                    "hausdorff_loss"])
    per_tp.to_csv(out_dir / "per_timepoint_forecast_metrics.csv", index=False)
    print(f"[eval_forecast] Metrics written to {metrics_path}")
