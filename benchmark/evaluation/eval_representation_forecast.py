"""
eval_representation_forecast.py
===============================

Representation-space forecast evaluator (Work-Plan Step 5).

Compares a *predicted* representation distribution against the *observed*
held-out representation distribution at each target timepoint, using the four
``rep_``-prefixed distribution metrics:

    rep_wasserstein_distance   (lower better)
    rep_gaussian_mmd           (lower better)
    rep_energy_distance_mmd    (lower better)
    rep_hausdorff_loss         (lower better)

Naming uses the ``rep_`` prefix so representation-space results are never
confused with the official gene-expression forecast metrics (Risk 1).

Backends
--------
``--backend auto`` (default): use the scTimeBench-aligned geomloss solvers from
``eval_forecast`` when torch+geomloss are importable (parity with the official
metrics), otherwise fall back to the self-contained numpy implementations here.
``--backend numpy`` / ``--backend geomloss`` force a specific backend.

The numpy implementations are exact for MMD / energy / Hausdorff and use a
stable log-domain Sinkhorn for the entropic Wasserstein distance, so the
evaluator (and its tests) run without any heavy dependency.

Inputs (per the work plan)
--------------------------
    projected_representation.npy   (n_pred_cells_total, d) or (n_tp, n_cells, d)
    observed_representation.npy    (n_obs_cells, d)
    --observed-times / metadata    timepoint per observed cell
    --eval-timepoints              target timepoints to score
    representation_metadata.json   (optional provenance, copied into output)
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import numpy as np

METRIC_KEYS = (
    "rep_wasserstein_distance",
    "rep_gaussian_mmd",
    "rep_energy_distance_mmd",
    "rep_hausdorff_loss",
)


# ===========================================================================
# Numpy metric implementations (exact; no torch required)
# ===========================================================================

def _pairwise_sq_dists(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Squared euclidean distances between rows of a and b."""
    a2 = np.sum(a * a, axis=1)[:, None]
    b2 = np.sum(b * b, axis=1)[None, :]
    d = a2 + b2 - 2.0 * a @ b.T
    return np.maximum(d, 0.0)


def numpy_gaussian_mmd(X: np.ndarray, Y: np.ndarray, blur: float = 1.0,
                       normalize_dims: bool = True) -> float:
    """Unbiased squared MMD with a Gaussian kernel exp(-||x-y||^2/(2 blur^2)).

    Divided by n_dims when ``normalize_dims`` (mirrors eval_forecast's
    normalize_by_n_genes for parity of scale).
    """
    X = np.asarray(X, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    d = X.shape[1]
    gamma = 1.0 / (2.0 * blur * blur)
    Kxx = np.exp(-gamma * _pairwise_sq_dists(X, X))
    Kyy = np.exp(-gamma * _pairwise_sq_dists(Y, Y))
    Kxy = np.exp(-gamma * _pairwise_sq_dists(X, Y))
    n, m = X.shape[0], Y.shape[0]
    # unbiased estimators (drop diagonal)
    sum_xx = (Kxx.sum() - np.trace(Kxx)) / (n * (n - 1)) if n > 1 else 0.0
    sum_yy = (Kyy.sum() - np.trace(Kyy)) / (m * (m - 1)) if m > 1 else 0.0
    sum_xy = Kxy.mean()
    mmd2 = sum_xx + sum_yy - 2.0 * sum_xy
    val = float(max(mmd2, 0.0))
    return val / d if normalize_dims else val


def numpy_energy_distance(X: np.ndarray, Y: np.ndarray,
                          normalize_dims: bool = True) -> float:
    """Energy distance E = 2 E||X-Y|| - E||X-X'|| - E||Y-Y'|| (>=0)."""
    X = np.asarray(X, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    d = X.shape[1]
    dxy = np.sqrt(_pairwise_sq_dists(X, Y)).mean()
    dxx = np.sqrt(_pairwise_sq_dists(X, X)).mean()
    dyy = np.sqrt(_pairwise_sq_dists(Y, Y)).mean()
    val = float(max(2.0 * dxy - dxx - dyy, 0.0))
    return val / d if normalize_dims else val


def numpy_hausdorff(X: np.ndarray, Y: np.ndarray) -> float:
    """Bidirectional Hausdorff: max over each set of nearest-neighbor distance.

    Matches eval_forecast._solver_hausdorff (NOT normalized by n_dims).
    """
    X = np.asarray(X, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    if X.size == 0 or Y.size == 0:
        return float("nan")
    dist = np.sqrt(_pairwise_sq_dists(X, Y))
    min_x = dist.min(axis=1)  # each X to nearest Y
    min_y = dist.min(axis=0)  # each Y to nearest X
    return float(max(min_x.max(), min_y.max()))


def numpy_sinkhorn_wasserstein(X: np.ndarray, Y: np.ndarray, blur: float = 0.05,
                               n_iter: int = 200, normalize_dims: bool = True) -> float:
    """
    Debiased entropic 2-Wasserstein (Sinkhorn) cost with squared-euclidean
    ground cost, uniform marginals, stable log-domain iterations.

    Debiasing subtracts the self-transport terms (Sinkhorn divergence), so the
    distance of a sample to itself is ~0 — matching geomloss debias=True used
    in eval_forecast. Divided by n_dims when ``normalize_dims``.
    """
    X = np.asarray(X, dtype=np.float64)
    Y = np.asarray(Y, dtype=np.float64)
    d = X.shape[1]
    eps = blur * blur  # geomloss relates blur to sqrt(eps) for p=2

    def _ot(a_pts, b_pts):
        C = _pairwise_sq_dists(a_pts, b_pts)
        n, m = a_pts.shape[0], b_pts.shape[0]
        log_a = np.full(n, -np.log(n))
        log_b = np.full(m, -np.log(m))
        f = np.zeros(n)
        g = np.zeros(m)
        for _ in range(n_iter):
            # f_i = -eps * logsumexp_j( log_b_j + (g_j - C_ij)/eps )
            f = -eps * _logsumexp(log_b[None, :] + (g[None, :] - C) / eps, axis=1)
            # g_j = -eps * logsumexp_i( log_a_i + (f_i - C_ij)/eps )
            g = -eps * _logsumexp(log_a[:, None] + (f[:, None] - C) / eps, axis=0)
        # transport plan: log P_ij = log_a_i + log_b_j + (f_i + g_j - C_ij)/eps
        log_P = log_a[:, None] + log_b[None, :] + (f[:, None] + g[None, :] - C) / eps
        P = np.exp(log_P)
        return float(np.sum(P * C))

    cost_xy = _ot(X, Y)
    cost_xx = _ot(X, X)
    cost_yy = _ot(Y, Y)
    div = max(cost_xy - 0.5 * (cost_xx + cost_yy), 0.0)
    return div / d if normalize_dims else div


def _logsumexp(a: np.ndarray, axis: int) -> np.ndarray:
    amax = np.max(a, axis=axis, keepdims=True)
    amax = np.where(np.isfinite(amax), amax, 0.0)
    out = np.log(np.sum(np.exp(a - amax), axis=axis, keepdims=True)) + amax
    return np.squeeze(out, axis=axis)


# ===========================================================================
# Backend selection
# ===========================================================================

def _resolve_backend(backend: str):
    if backend == "numpy":
        return "numpy"
    try:
        import geomloss  # noqa: F401
        import torch  # noqa: F401
    except Exception:
        if backend == "geomloss":
            raise RuntimeError(
                "backend='geomloss' requested but torch/geomloss not importable."
            )
        return "numpy"
    return "geomloss"


def compute_rep_metrics(X_pred: np.ndarray, X_obs: np.ndarray,
                        backend: str = "auto") -> dict:
    """Compute the four rep_ metrics for one timepoint's pred vs obs clouds."""
    chosen = _resolve_backend(backend)
    if chosen == "geomloss":
        from benchmark.evaluation.eval_forecast import (
            energy_distance_mmd,
            gaussian_mmd,
            hausdorff_loss,
            wasserstein_distance,
        )
        return {
            "rep_wasserstein_distance": wasserstein_distance(X_pred, X_obs),
            "rep_gaussian_mmd": gaussian_mmd(X_pred, X_obs),
            "rep_energy_distance_mmd": energy_distance_mmd(X_pred, X_obs),
            "rep_hausdorff_loss": hausdorff_loss(X_pred, X_obs),
            "backend": "geomloss_sctimebench",
        }
    return {
        "rep_wasserstein_distance": numpy_sinkhorn_wasserstein(X_pred, X_obs),
        "rep_gaussian_mmd": numpy_gaussian_mmd(X_pred, X_obs),
        "rep_energy_distance_mmd": numpy_energy_distance(X_pred, X_obs),
        "rep_hausdorff_loss": numpy_hausdorff(X_pred, X_obs),
        "backend": "numpy_exact",
    }


# ===========================================================================
# I/O orchestration
# ===========================================================================

def _split_pred(X_pred_raw: np.ndarray, eval_timepoints, n_pred_cells=None):
    """Return dict tp -> array. Supports 3D (n_tp,n_cells,d) or 2D blocks."""
    n_tp = len(eval_timepoints)
    if X_pred_raw.ndim == 3:
        if X_pred_raw.shape[0] != n_tp:
            raise ValueError(
                f"3D projected_representation has {X_pred_raw.shape[0]} blocks "
                f"but {n_tp} eval_timepoints."
            )
        return {float(eval_timepoints[i]): np.asarray(X_pred_raw[i], dtype=np.float64)
                for i in range(n_tp)}
    if X_pred_raw.ndim == 2:
        if n_pred_cells is None:
            if n_tp == 0 or X_pred_raw.shape[0] % n_tp != 0:
                raise ValueError(
                    "Cannot infer n_pred_cells for 2D projected_representation; "
                    "pass --n-pred-cells or provide a 3D array."
                )
            n_pred_cells = X_pred_raw.shape[0] // n_tp
        return {
            float(eval_timepoints[i]):
                np.asarray(X_pred_raw[i * n_pred_cells:(i + 1) * n_pred_cells],
                           dtype=np.float64)
            for i in range(n_tp)
        }
    raise ValueError(f"Unsupported projected_representation shape {X_pred_raw.shape}.")


def run_representation_forecast_evaluation(
    projected_representation_path: str,
    observed_representation_path: str,
    observed_times: np.ndarray,
    eval_timepoints: list,
    output_dir: str,
    *,
    backend: str = "auto",
    n_pred_cells: Optional[int] = None,
    representation_metadata: Optional[dict] = None,
    aggregate: str = "mean",
) -> dict:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    X_pred_raw = np.load(projected_representation_path, allow_pickle=True)
    X_obs_all = np.load(observed_representation_path, allow_pickle=True)
    observed_times = np.asarray(observed_times, dtype=float)
    eval_timepoints = [float(t) for t in eval_timepoints]

    pred_by_tp = _split_pred(X_pred_raw, eval_timepoints, n_pred_cells)

    per_tp = {}
    for tp in eval_timepoints:
        obs_mask = np.isclose(observed_times, tp)
        X_obs = np.asarray(X_obs_all[obs_mask], dtype=np.float64)
        X_pred = pred_by_tp[tp]
        if X_obs.shape[0] == 0 or X_pred.shape[0] == 0:
            continue
        m = compute_rep_metrics(X_pred, X_obs, backend=backend)
        m["n_pred_cells"] = int(X_pred.shape[0])
        m["n_obs_cells"] = int(X_obs.shape[0])
        per_tp[str(tp)] = m

    agg = {}
    for key in METRIC_KEYS:
        vals = [v[key] for v in per_tp.values() if v.get(key) is not None]
        arr = np.array(vals, dtype=float)
        agg[key] = float(np.mean(arr)) if len(arr) and aggregate == "mean" else (
            float(np.median(arr)) if len(arr) else None
        )

    backend_used = next(iter(per_tp.values()))["backend"] if per_tp else _resolve_backend(backend)
    metrics = {
        **agg,
        "n_eval_timepoints": len(eval_timepoints),
        "eval_timepoints": eval_timepoints,
        "aggregate": aggregate,
        "metric_space": "representation",
        "metric_backend": backend_used,
        "metric_definitions": {
            "rep_wasserstein_distance": "debiased entropic 2-Wasserstein / n_dims",
            "rep_gaussian_mmd": "unbiased gaussian MMD^2 (blur=1.0) / n_dims",
            "rep_energy_distance_mmd": "energy distance / n_dims",
            "rep_hausdorff_loss": "bidirectional Hausdorff (not normalized)",
        },
        "representation_metadata": representation_metadata or {},
        "status": "completed" if per_tp else "no_overlapping_timepoints",
    }

    with open(out_dir / "representation_forecast_metrics.json", "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)
    import csv
    with open(out_dir / "representation_forecast_per_timepoint.csv", "w",
              newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["timepoint", *METRIC_KEYS, "n_pred_cells", "n_obs_cells", "backend"])
        for tp_str, v in per_tp.items():
            w.writerow([tp_str, *[v.get(k) for k in METRIC_KEYS],
                        v.get("n_pred_cells"), v.get("n_obs_cells"), v.get("backend")])
    return metrics


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Representation-space forecast evaluator.")
    p.add_argument("--projected-representation", required=True)
    p.add_argument("--observed-representation", required=True)
    p.add_argument("--observed-times", required=True,
                   help="Path to .npy of per-observed-cell timepoints.")
    p.add_argument("--eval-timepoints", required=True, nargs="+", type=float)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--backend", default="auto", choices=["auto", "numpy", "geomloss"])
    p.add_argument("--n-pred-cells", type=int, default=None)
    p.add_argument("--representation-metadata", default=None,
                   help="Optional representation_metadata.json to embed in output.")
    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    rep_meta = None
    if args.representation_metadata and Path(args.representation_metadata).exists():
        with open(args.representation_metadata, encoding="utf-8") as f:
            rep_meta = json.load(f)
    run_representation_forecast_evaluation(
        projected_representation_path=args.projected_representation,
        observed_representation_path=args.observed_representation,
        observed_times=np.load(args.observed_times, allow_pickle=True),
        eval_timepoints=args.eval_timepoints,
        output_dir=args.output_dir,
        backend=args.backend,
        n_pred_cells=args.n_pred_cells,
        representation_metadata=rep_meta,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
