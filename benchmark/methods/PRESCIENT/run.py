"""
PRESCIENT runner for the scTimeBench-aligned benchmark.

This module wraps the upstream PRESCIENT implementation without modifying the
author's repository. PRESCIENT trains and simulates in PCA space, so this
runner owns the benchmark-specific data conversion, inverse transformation to
gene expression, embedding metrics, and state-level lineage aggregation.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import random
import sys
import time
import traceback
import types
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from scipy.spatial.distance import cdist as scipy_cdist
from scipy.stats import wasserstein_distance as scipy_wasserstein
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler


LOCAL_PRESCIENT_REPO = Path(__file__).resolve().parent / "prescient_module"
FALLBACK_PRESCIENT_REPO = Path(r"C:\Users\37620\Documents\GitHub\prescient")


def _find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env and Path(env).exists():
        return Path(env)
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "benchmark").is_dir() and (parent / "data").is_dir():
            return parent
    raise RuntimeError("Cannot locate project root. Set TRAJ_PROJECT_ROOT.")


def _find_prescient_repo(project_root: Path) -> Path:
    candidates = []
    env = os.environ.get("PRESCIENT_REPO")
    if env:
        env_path = Path(env)
        if (env_path / "prescient" / "train" / "model.py").exists():
            return env_path
        raise RuntimeError(
            "PRESCIENT_REPO is set but does not contain "
            f"prescient/train/model.py: {env_path}"
        )
    candidates.extend([
        LOCAL_PRESCIENT_REPO,
        FALLBACK_PRESCIENT_REPO,
        project_root / "external" / "prescient",
    ])
    for candidate in candidates:
        if (candidate / "prescient" / "train" / "model.py").exists():
            return candidate
    raise RuntimeError(
        "Cannot locate the PRESCIENT repository. Set PRESCIENT_REPO to the "
        "directory that contains prescient/train/model.py."
    )


def _ensure_prescient_import(project_root: Path) -> Path:
    repo = _find_prescient_repo(project_root)
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))
    if str(repo) not in sys.path:
        sys.path.insert(0, str(repo))
    _bootstrap_minimal_prescient(repo)
    return repo


def _load_module(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {module_name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _bootstrap_minimal_prescient(repo: Path) -> None:
    """
    Load only the PRESCIENT modules needed for training/simulation.

    The upstream package's top-level ``prescient.__init__`` imports utility
    modules that require optional visualization/ANN dependencies such as
    ``annoy``. The benchmark runner does not use those utilities, so we create
    a minimal package namespace and load train/model/run directly.
    """
    pkg_dir = repo / "prescient"
    train_dir = pkg_dir / "train"
    commands_dir = pkg_dir / "commands"
    if "prescient" not in sys.modules:
        pkg = types.ModuleType("prescient")
        pkg.__path__ = [str(pkg_dir)]
        sys.modules["prescient"] = pkg
    if "prescient.train" not in sys.modules:
        train_pkg = types.ModuleType("prescient.train")
        train_pkg.__path__ = [str(train_dir)]
        sys.modules["prescient.train"] = train_pkg
    _load_module("prescient.train.util", train_dir / "util.py")
    _load_module("prescient.train.model", train_dir / "model.py")
    run_mod = _load_module("prescient.train.run", train_dir / "run.py")
    train_pkg = sys.modules["prescient.train"]
    for mod_name in ("prescient.train.util", "prescient.train.model"):
        mod = sys.modules[mod_name]
        for name, value in vars(mod).items():
            if not name.startswith("_"):
                setattr(train_pkg, name, value)
    train_pkg.run = run_mod.run
    if "prescient.commands" not in sys.modules:
        commands_pkg = types.ModuleType("prescient.commands")
        commands_pkg.__path__ = [str(commands_dir)]
        sys.modules["prescient.commands"] = commands_pkg
    train_model_mod = _load_module(
        "prescient.commands.train_model", commands_dir / "train_model.py"
    )
    setattr(sys.modules["prescient.commands"], "train_model", train_model_mod)


def _to_dense_float32(X) -> np.ndarray:
    if sp.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=np.float32)


def _materialize_training_adata(
    adata,
    time_key: str,
    train_times: list[float] | None,
    subsample_per_timepoint: int | None = None,
    seed: int = 0,
):
    if not train_times:
        base = adata.to_memory() if getattr(adata, "isbacked", False) else adata
    else:
        train_set = set(float(t) for t in train_times)
        mask = adata.obs[time_key].astype(float).isin(train_set).values
        subset = adata[mask]
        base = (
            subset.to_memory()
            if hasattr(subset, "to_memory") and getattr(subset, "isbacked", False)
            else subset.copy()
        )
    if not subsample_per_timepoint or int(subsample_per_timepoint) <= 0:
        return base

    rng = np.random.default_rng(seed)
    obs_times = base.obs[time_key].astype(float).to_numpy()
    keep = []
    for tp in sorted(np.unique(obs_times)):
        idx = np.where(np.isclose(obs_times, float(tp)))[0]
        if idx.size > int(subsample_per_timepoint):
            idx = rng.choice(idx, size=int(subsample_per_timepoint), replace=False)
        keep.append(np.sort(idx))
    keep_idx = np.sort(np.concatenate(keep)) if keep else np.array([], dtype=int)
    print(
        f"[PRESCIENT] subsample_per_timepoint={subsample_per_timepoint}: "
        f"{base.n_obs} -> {keep_idx.size} training cells"
    )
    subset = base[keep_idx]
    return (
        subset.to_memory()
        if hasattr(subset, "to_memory") and getattr(subset, "isbacked", False)
        else subset.copy()
    )


def _materialize_full_adata(adata):
    return adata.to_memory() if getattr(adata, "isbacked", False) else adata


def _fit_prescient_data(
    adata_train,
    time_key: str,
    cell_state_key: str,
    num_pcs: int,
    output_dir: Path,
) -> tuple[dict, Path]:
    """Build and cache a PRESCIENT-style torch object plus scaler metadata."""
    output_dir.mkdir(parents=True, exist_ok=True)
    data_path = output_dir / "data.pt"
    metadata_path = output_dir / "data_metadata.json"

    X = _to_dense_float32(adata_train.X)
    times = adata_train.obs[time_key].astype(float).to_numpy()
    states = adata_train.obs[cell_state_key].astype(str).to_numpy()
    train_times = sorted(np.unique(times).astype(float).tolist())
    time_to_code = {float(t): i for i, t in enumerate(train_times)}
    codes = np.array([time_to_code[float(t)] for t in times], dtype=int)

    n_components = min(int(num_pcs), X.shape[0] - 1, X.shape[1])
    if n_components < 2:
        raise ValueError(f"PRESCIENT requires at least 2 PCs, got {n_components}.")

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    pca = PCA(n_components=n_components, svd_solver="randomized", random_state=0)
    xp_all = pca.fit_transform(X_scaled).astype(np.float32)
    y = list(range(len(train_times)))

    x_split = [
        torch.from_numpy(X_scaled[codes == code].astype(np.float32))
        for code in y
    ]
    xp_split = [
        torch.from_numpy(xp_all[codes == code].astype(np.float32))
        for code in y
    ]
    # PRESCIENT stores UMAPs for visualization. The benchmark does not use
    # them, so keep a cheap two-dimensional PCA slice as a placeholder.
    xu_split = [
        torch.from_numpy(xp_all[codes == code, :2].astype(np.float32))
        for code in y
    ]
    # Uniform log-growth rates. The upstream train_model exponentiates these
    # into per-source-cell sampling weights, so zeros become uniform weights.
    w = [np.zeros(int((codes == code).sum()), dtype=np.float32) for code in y]

    genes = (
        adata_train.var["gene_symbol"].astype(str).to_numpy()
        if "gene_symbol" in adata_train.var.columns
        else adata_train.var_names.astype(str).to_numpy()
    )
    data_pt = {
        "data": X,
        "genes": genes,
        "celltype": states,
        "tps": codes,
        "x": x_split,
        "xp": xp_split,
        "xu": xu_split,
        "y": y,
        "pca": pca,
        "scaler": scaler,
        "w": w,
        "time_mapping": {
            "train_times": train_times,
            "time_to_code": {str(k): v for k, v in time_to_code.items()},
        },
    }
    torch.save(data_pt, data_path)
    metadata_path.write_text(
        json.dumps(
            {
                "n_cells": int(X.shape[0]),
                "n_genes": int(X.shape[1]),
                "num_pcs": int(n_components),
                "train_times": train_times,
                "growth_mode": "uniform_log_growth_zero",
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return data_pt, data_path


def _target_code_from_real_time(target_time: float, train_times: list[float]) -> float:
    """Map real observed time to the compressed PRESCIENT training axis."""
    train_times_arr = np.asarray(train_times, dtype=float)
    train_codes = np.arange(len(train_times_arr), dtype=float)
    target_time = float(target_time)
    if len(train_times_arr) == 1:
        return 0.0
    if target_time <= train_times_arr[0]:
        return float(train_codes[0])
    if target_time >= train_times_arr[-1]:
        slope = (train_codes[-1] - train_codes[-2]) / (
            train_times_arr[-1] - train_times_arr[-2]
        )
        return float(train_codes[-1] + slope * (target_time - train_times_arr[-1]))
    return float(np.interp(target_time, train_times_arr, train_codes))


def _load_model(seed_dir: Path, checkpoint_name: str = "best", device=None):
    from prescient.train.model import AutoGenerator

    device = device or torch.device("cpu")
    config = SimpleNamespace(**_torch_load_trusted(seed_dir / "config.pt", map_location=device))
    model = AutoGenerator(config)
    if checkpoint_name == "best":
        ckpt_path = seed_dir / "train.best.pt"
    else:
        ckpt_path = seed_dir / f"train.epoch_{checkpoint_name}.pt"
    checkpoint = _torch_load_trusted(ckpt_path, map_location=device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.to(device)
    model.eval()
    return model, config


def _torch_load_trusted(path, **kwargs):
    """Load benchmark-generated PRESCIENT files under PyTorch 2.6+."""
    try:
        return torch.load(path, weights_only=False, **kwargs)
    except TypeError:
        return torch.load(path, **kwargs)


def _call_with_trusted_torch_load(func, *args, **kwargs):
    """Temporarily make upstream PRESCIENT torch.load calls PyTorch-2.6 safe."""
    original_load = torch.load

    def trusted_load(*load_args, **load_kwargs):
        load_kwargs.setdefault("weights_only", False)
        return original_load(*load_args, **load_kwargs)

    torch.load = trusted_load
    try:
        return func(*args, **kwargs)
    finally:
        torch.load = original_load


def _simulate_pc(
    model,
    config,
    source_pc: np.ndarray,
    target_code: float,
    source_code: float = 0.0,
    device=None,
    seed: int = 0,
) -> np.ndarray:
    """Simulate source PCs to a target code on the PRESCIENT time axis."""
    device = device or torch.device("cpu")
    rng = torch.Generator(device="cpu")
    rng.manual_seed(int(seed))
    dt = float(config.train_dt)
    delta = max(0.0, float(target_code) - float(source_code))
    n_steps = int(np.round(delta / dt))
    residual = float(delta - n_steps * dt)
    x = torch.from_numpy(source_pc.astype(np.float32)).to(device)
    # AutoGenerator._drift uses autograd.grad even during simulation, so we
    # must keep gradient tracking enabled for each step. Detach after each
    # step to avoid retaining a long inference graph.
    with torch.enable_grad():
        for _ in range(n_steps):
            x = x.detach()
            z = torch.randn(
                x.shape[0], x.shape[1], generator=rng, dtype=x.dtype
            ).to(device) * float(config.train_sd)
            x = model._step(x.float(), dt=dt, z=z).detach()
        if residual > 1e-8:
            x = x.detach()
            z = torch.randn(
                x.shape[0], x.shape[1], generator=rng, dtype=x.dtype
            ).to(device) * float(config.train_sd)
            x = model._step(x.float(), dt=residual, z=z).detach()
    return x.detach().cpu().numpy()


def _inverse_expression(pc: np.ndarray, pca: PCA, scaler: StandardScaler) -> np.ndarray:
    X_scaled = pca.inverse_transform(pc)
    X = scaler.inverse_transform(X_scaled)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    return np.clip(X, 0.0, None).astype(np.float32)


def _sample_rows(X: np.ndarray, n: int, seed: int) -> np.ndarray:
    if X.shape[0] <= n:
        return X
    rng = np.random.default_rng(seed)
    idx = rng.choice(X.shape[0], size=n, replace=False)
    return X[np.sort(idx)]


def _mean_gene_wasserstein(X_pred: np.ndarray, X_obs: np.ndarray) -> float:
    if X_pred.size == 0 or X_obs.size == 0:
        return float("nan")
    vals = [
        scipy_wasserstein(X_pred[:, j], X_obs[:, j])
        for j in range(X_pred.shape[1])
    ]
    return float(np.mean(vals))


def _write_forecast_outputs(
    model,
    config,
    data_pt: dict,
    adata_full,
    time_key: str,
    eval_times: list[float],
    n_sim_cells: int,
    metric_sample_cells: int,
    output_dir: Path,
    seed: int,
) -> None:
    from benchmark.evaluation.eval_forecast import (
        energy_distance_mmd,
        gaussian_mmd,
        hausdorff_loss,
    )

    train_times = data_pt["time_mapping"]["train_times"]
    pca = data_pt["pca"]
    scaler = data_pt["scaler"]
    t0_code = 0
    xp_t0 = data_pt["xp"][t0_code].numpy()
    rng = np.random.default_rng(seed)
    replace = xp_t0.shape[0] < n_sim_cells
    src_idx = rng.choice(xp_t0.shape[0], size=n_sim_cells, replace=replace)
    source_pc = xp_t0[src_idx]

    full_times = adata_full.obs[time_key].astype(float).to_numpy()
    all_pred_expr = []
    per_tp_rows = []
    for i, tp in enumerate(eval_times):
        target_code = _target_code_from_real_time(float(tp), train_times)
        pred_pc = _simulate_pc(
            model, config, source_pc, target_code, seed=seed + 1000 + i
        )
        X_pred = _inverse_expression(pred_pc, pca, scaler)
        obs_mask = np.isclose(full_times, float(tp))
        X_obs = _to_dense_float32(adata_full[obs_mask].X)
        X_pred_s = _sample_rows(X_pred, metric_sample_cells, seed=seed + i)
        X_obs_s = _sample_rows(X_obs, metric_sample_cells, seed=seed + 100 + i)
        all_pred_expr.append(X_pred)
        per_tp_rows.append({
            "timepoint": float(tp),
            "target_code": float(target_code),
            "wasserstein_distance": _mean_gene_wasserstein(X_pred_s, X_obs_s),
            "gaussian_mmd": gaussian_mmd(X_pred_s, X_obs_s),
            "energy_distance_mmd": energy_distance_mmd(X_pred_s, X_obs_s),
            "hausdorff_loss": hausdorff_loss(X_pred_s, X_obs_s),
            "n_obs_cells": int(X_obs.shape[0]),
            "n_pred_cells": int(X_pred.shape[0]),
            "metric_sample_pred_cells": int(X_pred_s.shape[0]),
            "metric_sample_obs_cells": int(X_obs_s.shape[0]),
        })

    proj = np.stack(all_pred_expr, axis=0) if all_pred_expr else np.array([])
    np.save(output_dir / "projected_expression.npy", proj)
    per_tp = pd.DataFrame(per_tp_rows)
    per_tp.to_csv(output_dir / "per_timepoint_forecast_metrics.csv", index=False)

    def _mean(col):
        return float(np.nanmean(per_tp[col].astype(float))) if col in per_tp else None

    metrics = {
        "wasserstein_distance": _mean("wasserstein_distance"),
        "gaussian_mmd": _mean("gaussian_mmd"),
        "energy_distance_mmd": _mean("energy_distance_mmd"),
        "hausdorff_loss": _mean("hausdorff_loss"),
        "n_eval_timepoints": len(eval_times),
        "eval_timepoints": [float(t) for t in eval_times],
        "scenario_type": "forecast" if eval_times != train_times else "reconstruction",
        "metric_sample_cells": int(metric_sample_cells),
        "status": "completed",
    }
    (output_dir / "forecast_metrics.json").write_text(
        json.dumps(metrics, indent=2), encoding="utf-8"
    )


def _write_embedding_outputs(
    model,
    config,
    data_pt: dict,
    adata_full,
    time_key: str,
    cell_state_key: str,
    eval_times: list[float],
    n_sim_cells: int,
    output_dir: Path,
    seed: int,
    max_reference_cells_per_timepoint: int | None = None,
) -> None:
    train_times = data_pt["time_mapping"]["train_times"]
    xp_all = data_pt["pca"].transform(
        data_pt["scaler"].transform(_to_dense_float32(adata_full.X))
    ).astype(np.float32)
    full_times = adata_full.obs[time_key].astype(float).to_numpy()
    full_states = adata_full.obs[cell_state_key].astype(str).to_numpy()
    train_states = np.asarray(data_pt["celltype"]).astype(str)
    xp_train = np.vstack([x.numpy() for x in data_pt["xp"]])

    labels_unique = sorted(set(train_states))
    centroids = {
        s: xp_train[train_states == s].mean(axis=0)
        for s in labels_unique
        if np.any(train_states == s)
    }
    clf = LogisticRegression(max_iter=500, class_weight="balanced", random_state=0)
    clf.fit(xp_train, train_states)

    xp_t0 = data_pt["xp"][0].numpy()
    rng = np.random.default_rng(seed)
    replace = xp_t0.shape[0] < n_sim_cells
    src_idx = rng.choice(xp_t0.shape[0], size=n_sim_cells, replace=replace)
    source_pc = xp_t0[src_idx]

    all_pred_pc = []
    label_rows = []
    ari_rows = []
    entropy_values = []
    centroid_labels = list(centroids)
    centroid_mat = np.stack([centroids[s] for s in centroid_labels])
    for i, tp in enumerate(eval_times):
        target_code = _target_code_from_real_time(float(tp), train_times)
        pred_pc = _simulate_pc(
            model, config, source_pc, target_code, seed=seed + 2000 + i
        )
        all_pred_pc.append(pred_pc)
        dists = scipy_cdist(pred_pc, centroid_mat, metric="euclidean")
        pred_labels = np.asarray([centroid_labels[j] for j in np.argmin(dists, axis=1)])

        probs = clf.predict_proba(pred_pc)
        if probs.shape[1] > 1:
            eps = 1e-12
            entropy = -np.sum(probs * np.log(probs + eps), axis=1)
            entropy_values.extend((entropy / np.log(probs.shape[1])).tolist())

        obs_mask = np.isclose(full_times, float(tp))
        if obs_mask.any():
            obs_idx = np.where(obs_mask)[0]
            if (
                max_reference_cells_per_timepoint
                and obs_idx.size > int(max_reference_cells_per_timepoint)
            ):
                ref_rng = np.random.default_rng(seed + 5000 + i)
                obs_idx = np.sort(
                    ref_rng.choice(
                        obs_idx,
                        size=int(max_reference_cells_per_timepoint),
                        replace=False,
                    )
                )
            obs_pc = xp_all[obs_idx]
            obs_labels = full_states[obs_idx]
            nn_idx = np.argmin(scipy_cdist(pred_pc, obs_pc), axis=1)
            ref_labels = obs_labels[nn_idx]
            ari_rows.append({
                "timepoint": float(tp),
                "adjusted_rand_index": float(adjusted_rand_score(ref_labels, pred_labels)),
                "n_projected_cells": int(pred_pc.shape[0]),
                "n_observed_reference_cells": int(obs_pc.shape[0]),
            })

        for j, label in enumerate(pred_labels):
            label_rows.append({
                "projected_timepoint": float(tp),
                "cell_idx": int(j),
                "projected_cluster_label": str(label),
            })

    emb = np.vstack(all_pred_pc) if all_pred_pc else np.array([])
    np.save(output_dir / "projected_embedding.npy", emb)
    pd.DataFrame(label_rows).to_csv(
        output_dir / "projected_cluster_labels.csv", index=False
    )
    pd.DataFrame(ari_rows).to_csv(
        output_dir / "per_timepoint_embedding_metrics.csv", index=False
    )
    metrics = {
        "adjusted_rand_index": (
            float(np.mean([r["adjusted_rand_index"] for r in ari_rows]))
            if ari_rows else None
        ),
        "avg_normalized_classifier_entropy": (
            float(np.mean(entropy_values)) if entropy_values else None
        ),
        "projected_latent_dim": int(emb.shape[1]) if emb.size else None,
        "n_eval_timepoints": len(eval_times),
        "eval_timepoints": [float(t) for t in eval_times],
        "status": "completed" if ari_rows and entropy_values else "embedding_metrics_incomplete",
        "embedding_space": "PRESCIENT PCA space",
    }
    (output_dir / "embedding_metrics.json").write_text(
        json.dumps(metrics, indent=2), encoding="utf-8"
    )


def _softmax_assign(pred_pc: np.ndarray, centroids: np.ndarray, temperature: float):
    dists = scipy_cdist(pred_pc, centroids, metric="euclidean")
    scores = -dists / float(temperature)
    row_max = np.max(scores, axis=1, keepdims=True)
    weights = np.exp(scores - row_max)
    weights_sum = weights.sum(axis=1, keepdims=True)
    return np.divide(weights, weights_sum, out=np.zeros_like(weights), where=weights_sum > 0)


def _write_lineage_outputs(
    model,
    config,
    data_pt: dict,
    cell_state_key: str,
    output_dir: Path,
    assignment_temperature: float,
    seed: int,
) -> None:
    states = np.asarray(data_pt["celltype"]).astype(str)
    unique_states = sorted(set(states))
    state_to_idx = {s: i for i, s in enumerate(unique_states)}
    xp_by_code = [x.numpy() for x in data_pt["xp"]]
    state_by_code = [states[np.asarray(data_pt["tps"]) == code] for code in data_pt["y"]]
    xp_train = np.vstack(xp_by_code)
    centroids = np.stack([
        xp_train[states == s].mean(axis=0) if np.any(states == s)
        else np.zeros(xp_train.shape[1], dtype=np.float32)
        for s in unique_states
    ])

    stm_raw = np.zeros((len(unique_states), len(unique_states)), dtype=np.float64)
    stm_counts = np.zeros(len(unique_states), dtype=np.float64)
    for code in range(len(xp_by_code) - 1):
        source_pc = xp_by_code[code]
        pred_pc = _simulate_pc(
            model, config, source_pc, target_code=code + 1,
            source_code=code, seed=seed + 3000 + code
        )
        weights = _softmax_assign(pred_pc, centroids, assignment_temperature)
        for i, src_state in enumerate(state_by_code[code]):
            src_idx = state_to_idx[str(src_state)]
            stm_raw[src_idx] += weights[i]
            stm_counts[src_idx] += 1.0

    stm = np.zeros_like(stm_raw)
    nonzero = stm_counts > 0
    stm[nonzero] = stm_raw[nonzero] / stm_counts[nonzero, None]
    stm = np.nan_to_num(stm, nan=0.0, posinf=0.0, neginf=0.0)
    stm_df = pd.DataFrame(
        stm,
        index=pd.Index(unique_states, name="source_state"),
        columns=unique_states,
    )
    stm_df.to_csv(output_dir / "state_transition_matrix.csv")
    edge_rows = [
        {"source_state": src, "target_state": tgt, "weight": float(stm_df.loc[src, tgt])}
        for src in unique_states
        for tgt in unique_states
        if float(stm_df.loc[src, tgt]) > 1e-8
    ]
    pd.DataFrame(edge_rows).to_csv(output_dir / "lineage_graph_edges.csv", index=False)
    diagnostics = {
        "cell_state_key": cell_state_key,
        "assignment_temperature": float(assignment_temperature),
        "source_transition_counts": {
            s: int(stm_counts[state_to_idx[s]]) for s in unique_states
        },
    }
    (output_dir / "lineage_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2), encoding="utf-8"
    )


def train_or_load_prescient(
    data_path: Path,
    output_dir: Path,
    prescient_cfg: dict,
    device,
) -> Path:
    from prescient.commands import train_model
    import prescient.train as prescient_train

    seed = int(prescient_cfg.get("seed", 42))
    weight_name = str(prescient_cfg.get("weight_name", "uniform"))
    activation = str(prescient_cfg.get("activation", "softplus"))
    layers = int(prescient_cfg.get("layers", 1))
    k_dim = int(prescient_cfg.get("k_dim", 200))
    train_tau = float(prescient_cfg.get("train_tau", 1e-6))
    model_parent = output_dir / (
        f"{weight_name}-{activation}_{layers}_{k_dim}-{train_tau}"
    )
    seed_dir = model_parent / f"seed_{seed}"
    if (seed_dir / "train.best.pt").exists() and (seed_dir / "config.pt").exists():
        print(f"[PRESCIENT] Loading cached model from {seed_dir}")
        return seed_dir

    args = SimpleNamespace(
        no_cuda=(device.type == "cpu"),
        gpu=int(prescient_cfg.get("gpu", 0)),
        out_dir=str(output_dir),
        seed=seed,
        data_path=str(data_path),
        weight_name=weight_name,
        weight=weight_name,
        loss=str(prescient_cfg.get("loss", "euclidean")),
        k_dim=k_dim,
        activation=activation,
        layers=layers,
        pretrain_epochs=int(prescient_cfg.get("pretrain_epochs", 50)),
        train_epochs=int(prescient_cfg.get("train_epochs", 200)),
        train_lr=float(prescient_cfg.get("train_lr", 0.01)),
        train_dt=float(prescient_cfg.get("train_dt", 0.1)),
        train_sd=float(prescient_cfg.get("train_sd", 0.5)),
        train_tau=train_tau,
        train_batch=float(prescient_cfg.get("train_batch", 0.1)),
        train_clip=float(prescient_cfg.get("train_clip", 0.25)),
        save=int(prescient_cfg.get("save", 100)),
        pretrain=bool(prescient_cfg.get("pretrain", True)),
        train=bool(prescient_cfg.get("train", True)),
        config=None,
    )
    print(f"[PRESCIENT] Training model in {output_dir}")
    result = _call_with_trusted_torch_load(
        prescient_train.run,
        args,
        train_model.train_init,
    )
    if result is not None and not (seed_dir / "train.best.pt").exists():
        model, best_state, config, _loss_history = result
        seed_dir.mkdir(parents=True, exist_ok=True)
        torch.save(config.__dict__, seed_dir / "config.pt")
        checkpoint = best_state or {"model_state_dict": model.state_dict()}
        if "model_state_dict" not in checkpoint:
            checkpoint = {"model_state_dict": checkpoint}
        torch.save(checkpoint, seed_dir / "train.best.pt")
    return seed_dir


def run_pipeline(
    adata,
    scenario_id: str,
    scenario_config: dict,
    output_dir: str | Path,
) -> dict:
    project_root = _find_project_root()
    prescient_repo = _ensure_prescient_import(project_root)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    time_key = scenario_config.get("time_key", "abs_day")
    cell_state_key = scenario_config.get("cell_state_key", "final_milestone_label_coarse")
    dataset_id = scenario_config.get("dataset_id", "GSE230659")
    prescient_cfg = scenario_config.get("prescient_params", {}) or {}
    scenario_params = scenario_config.get("scenario_params", {}) or {}
    train_times_raw = scenario_params.get("train_times")
    heldout_times = [float(t) for t in scenario_params.get("heldout_times", [])]
    train_times = [float(t) for t in train_times_raw] if train_times_raw else None
    seed = int(prescient_cfg.get("seed", 42))
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    if time_key not in adata.obs.columns or cell_state_key not in adata.obs.columns:
        raise KeyError(
            f"Missing required obs columns for PRESCIENT: {time_key!r}, "
            f"{cell_state_key!r}."
        )

    t0 = time.time()
    status = "failed"
    notes = ""
    try:
        train_subsample = prescient_cfg.get("subsample_per_timepoint")
        train_adata = _materialize_training_adata(
            adata,
            time_key,
            train_times,
            subsample_per_timepoint=(
                int(train_subsample) if train_subsample not in (None, "", False) else None
            ),
            seed=seed,
        )
        if heldout_times:
            full_adata = _materialize_full_adata(adata)
        else:
            full_adata = train_adata
        internal_dir = output_dir / "prescient_internal"
        data_pt, data_path = _fit_prescient_data(
            train_adata,
            time_key=time_key,
            cell_state_key=cell_state_key,
            num_pcs=int(prescient_cfg.get("num_pcs", 50)),
            output_dir=internal_dir,
        )

        use_cuda = bool(prescient_cfg.get("use_cuda", False)) and torch.cuda.is_available()
        gpu = int(prescient_cfg.get("gpu", 0))
        device = torch.device(f"cuda:{gpu}" if use_cuda else "cpu")
        model_seed_dir = train_or_load_prescient(
            data_path=data_path,
            output_dir=internal_dir / "models",
            prescient_cfg=prescient_cfg,
            device=device,
        )
        model, model_config = _load_model(model_seed_dir, device=device)

        eval_times = heldout_times or data_pt["time_mapping"]["train_times"]
        n_sim_cells = int(prescient_cfg.get("n_sim_cells", 2000))
        metric_sample_cells = int(prescient_cfg.get("metric_sample_cells", 1000))

        _write_forecast_outputs(
            model, model_config, data_pt, full_adata, time_key, eval_times,
            n_sim_cells, metric_sample_cells, output_dir, seed
        )
        _write_embedding_outputs(
            model, model_config, data_pt, full_adata, time_key, cell_state_key,
            eval_times, n_sim_cells, output_dir, seed,
            max_reference_cells_per_timepoint=prescient_cfg.get(
                "embedding_reference_cells_per_timepoint"
            ),
        )
        _write_lineage_outputs(
            model, model_config, data_pt, cell_state_key, output_dir,
            assignment_temperature=float(prescient_cfg.get("assignment_temperature", 1.0)),
            seed=seed,
        )
        status = "completed"
    except Exception:
        notes = traceback.format_exc()
        print(f"[PRESCIENT] ERROR:\n{notes}")
        raise
    finally:
        metadata = {
            "method": "prescient",
            "dataset": dataset_id,
            "scenario": scenario_id,
            "result_class": scenario_config.get("result_class", "official"),
            "formal_benchmark": bool(scenario_config.get("formal_benchmark", True)),
            "capability_flags": {
                "supports_unseen_timepoint_projection": True,
                "supports_lineage_inference": True,
            },
            "dimensions_executed": [
                "forecast_accuracy",
                "embedding_coherence",
                "lineage_fidelity",
            ],
            "runtime_seconds": round(time.time() - t0, 2),
            "status": status,
            "time_key": time_key,
            "cell_state_key": cell_state_key,
            "train_times": train_times if train_times else "all",
            "heldout_times": heldout_times,
            "n_sim_cells": prescient_cfg.get("n_sim_cells", 2000),
            "sampling": {
                "subsample_per_timepoint": prescient_cfg.get("subsample_per_timepoint"),
                "embedding_reference_cells_per_timepoint": prescient_cfg.get(
                    "embedding_reference_cells_per_timepoint"
                ),
            },
            "prescient_repo": str(prescient_repo),
            "prescient_params": prescient_cfg,
            "notes": notes,
        }
        (output_dir / "run_metadata.json").write_text(
            json.dumps(metadata, indent=2), encoding="utf-8"
        )

    return {
        "state_transition_matrix": str(output_dir / "state_transition_matrix.csv"),
        "lineage_graph_edges": str(output_dir / "lineage_graph_edges.csv"),
        "projected_expression": str(output_dir / "projected_expression.npy"),
        "forecast_metrics": str(output_dir / "forecast_metrics.json"),
        "per_timepoint_forecast_metrics": str(output_dir / "per_timepoint_forecast_metrics.csv"),
        "projected_embedding": str(output_dir / "projected_embedding.npy"),
        "embedding_metrics": str(output_dir / "embedding_metrics.json"),
        "projected_cluster_labels": str(output_dir / "projected_cluster_labels.csv"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run PRESCIENT benchmark pipeline.")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    import anndata
    import yaml

    project_root = _find_project_root()
    with open(args.config, encoding="utf-8-sig") as f:
        cfg = yaml.safe_load(f) or {}
    dataset_cfg = cfg.get("dataset", {}) or {}
    scenario_params = cfg.get("scenario_params", {}) or {}
    h5ad_path = Path(dataset_cfg["h5ad_path"])
    if not h5ad_path.is_absolute():
        h5ad_path = project_root / h5ad_path
    backed = "r" if scenario_params.get("train_times") else None
    adata = anndata.read_h5ad(h5ad_path, backed=backed)

    scenario_config = {
        "dataset_id": dataset_cfg.get("id", "GSE230659"),
        "time_key": dataset_cfg.get("time_key", "abs_day"),
        "cell_state_key": (cfg.get("ground_truth") or {}).get(
            "state_key",
            (cfg.get("lineage") or {}).get(
                "cell_state_key", "final_milestone_label_coarse"
            ),
        ),
        "scenario_params": scenario_params,
        "prescient_params": cfg.get("prescient_params", {}) or {},
        "result_class": cfg.get("result_class", "official"),
        "formal_benchmark": cfg.get("formal_benchmark", True),
    }
    out_base = Path((cfg.get("output") or {})["base_dir"])
    if not out_base.is_absolute():
        out_base = project_root / out_base
    run_pipeline(
        adata=adata,
        scenario_id=str(cfg.get("scenario", "A")),
        scenario_config=scenario_config,
        output_dir=out_base,
    )


if __name__ == "__main__":
    main()
