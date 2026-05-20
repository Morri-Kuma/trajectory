"""
MIOFlow runner for the scTimeBench-aligned benchmark.

This wrapper follows scTimeBench's practical treatment of MIOFlow: train the
neural ODE in PCA space, then inverse-transform projected points back to gene
expression for Forecast Accuracy. It uses the local vendored MIOFlow package in
this directory and avoids PHATE/GAGA for the first benchmark integration.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
from scipy.stats import wasserstein_distance as scipy_wasserstein
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import pairwise_distances


def _find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env and Path(env).exists():
        return Path(env)
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "benchmark").is_dir():
            return parent
    raise RuntimeError("Cannot locate project root. Set TRAJ_PROJECT_ROOT.")


_MIOFLOW_REPO = Path(__file__).resolve().parent
if str(_MIOFLOW_REPO) not in sys.path:
    sys.path.insert(0, str(_MIOFLOW_REPO))

from mioflow.core.datasets import TimeSeriesDataset  # noqa: E402
from mioflow.core.models.ode_model import ODEFunc  # noqa: E402
from mioflow.mioflow import train_mioflow  # noqa: E402
from torchdiffeq import odeint  # noqa: E402


def _to_dense_float32(X) -> np.ndarray:
    if sp.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=np.float32)


def _materialize(adata):
    return adata.to_memory() if getattr(adata, "isbacked", False) else adata


def _sample_rows(X: np.ndarray, max_cells: int, seed: int, replace: bool = False) -> np.ndarray:
    if X.shape[0] == 0:
        return X
    if X.shape[0] <= max_cells and not replace:
        return X
    rng = np.random.default_rng(seed)
    idx = rng.choice(X.shape[0], size=max_cells, replace=replace or X.shape[0] < max_cells)
    return X[np.sort(idx)]


def _safe_softmax_neg_distance(dists: np.ndarray, temperature: float = 1.0) -> np.ndarray:
    scores = -dists / max(float(temperature), 1e-8)
    finite = np.isfinite(scores)
    row_has_finite = finite.any(axis=1)
    row_max = np.zeros((scores.shape[0], 1), dtype=float)
    if row_has_finite.any():
        row_max[row_has_finite, 0] = np.max(
            np.where(finite[row_has_finite], scores[row_has_finite], -np.inf),
            axis=1,
        )
    shifted = np.where(finite, scores - row_max, -np.inf)
    weights = np.exp(shifted)
    row_sums = weights.sum(axis=1, keepdims=True)
    good = row_has_finite & np.isfinite(row_sums[:, 0]) & (row_sums[:, 0] > 0)
    weights = np.divide(weights, row_sums, out=np.zeros_like(weights), where=good[:, None])
    for i in np.where(~good)[0]:
        finite_d = np.isfinite(dists[i])
        if finite_d.any():
            weights[i, int(np.argmin(np.where(finite_d, dists[i], np.inf)))] = 1.0
    return weights


def _forecast_metrics(X_pred: np.ndarray, X_obs: np.ndarray) -> dict:
    if X_pred.size == 0 or X_obs.size == 0:
        return {
            "wasserstein_distance": None,
            "gaussian_mmd": None,
            "energy_distance_mmd": None,
            "hausdorff_loss": None,
        }

    wd = float(np.mean([
        scipy_wasserstein(X_pred[:, j], X_obs[:, j])
        for j in range(X_pred.shape[1])
    ]))

    X = np.vstack([X_pred, X_obs])
    X_med = _sample_rows(X, 500, seed=123)
    d = pairwise_distances(X_med, metric="euclidean")
    tri = d[np.triu_indices_from(d, k=1)]
    sigma = float(np.median(tri[tri > 0])) if np.any(tri > 0) else 1.0
    sigma = sigma if sigma > 0 else 1.0
    gamma = 1.0 / (2.0 * sigma * sigma)
    d_xx = pairwise_distances(X_pred, metric="sqeuclidean")
    d_yy = pairwise_distances(X_obs, metric="sqeuclidean")
    d_xy = pairwise_distances(X_pred, X_obs, metric="sqeuclidean")
    mmd = float(max(np.exp(-gamma * d_xx).mean()
                    + np.exp(-gamma * d_yy).mean()
                    - 2.0 * np.exp(-gamma * d_xy).mean(), 0.0))

    e_xy = pairwise_distances(X_pred, X_obs, metric="euclidean")
    e_xx = pairwise_distances(X_pred, metric="euclidean")
    e_yy = pairwise_distances(X_obs, metric="euclidean")
    energy = float(2.0 * e_xy.mean() - e_xx.mean() - e_yy.mean())
    hausdorff = float(max(e_xy.min(axis=1).max(), e_xy.min(axis=0).max()))
    return {
        "wasserstein_distance": wd,
        "gaussian_mmd": mmd,
        "energy_distance_mmd": energy,
        "hausdorff_loss": hausdorff,
    }


class MIOFlowBenchmarkModel:
    def __init__(self, model, pca, mean, std, time_to_idx, device: str):
        self.model = model
        self.pca = pca
        self.mean = mean
        self.std = std
        self.time_to_idx = time_to_idx
        self.device = device

    def project_pca(self, X_pca: np.ndarray, source_time: float, target_time: float) -> np.ndarray:
        if np.isclose(float(source_time), float(target_time)):
            return np.asarray(X_pca, dtype=np.float32)
        self.model.eval()
        x0 = ((X_pca - self.mean) / self.std).astype(np.float32)
        t = torch.tensor(
            [self.time_to_idx[float(source_time)], self.time_to_idx[float(target_time)]],
            dtype=torch.float32,
            device=self.device,
        )
        x0_t = torch.tensor(x0, dtype=torch.float32, device=self.device)
        with torch.no_grad():
            self.model.reset_momentum()
            pred = odeint(self.model, x0_t, t)[1].detach().cpu().numpy()
        return pred * self.std + self.mean

    def pca_to_expression(self, X_pca: np.ndarray) -> np.ndarray:
        return self.pca.inverse_transform(X_pca).astype(np.float32)


def train_or_load(adata_train, full_times: list[float], time_key: str, cfg: dict,
                  output_dir: Path) -> MIOFlowBenchmarkModel:
    seed = int(cfg.get("seed", 0))
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    use_cuda = bool(cfg.get("use_cuda", True)) and torch.cuda.is_available()
    device = "cuda" if use_cuda else "cpu"
    cache_path = output_dir / "trained_mioflow_model.pth"

    if cache_path.exists() and bool(cfg.get("use_cache", True)):
        payload = torch.load(cache_path, map_location="cpu", weights_only=False)
        model = ODEFunc(
            input_dim=int(payload["input_dim"]),
            hidden_dim=int(payload["hidden_dim"]),
            momentum_beta=float(payload.get("momentum_beta", 0.0)),
        )
        model.load_state_dict(payload["model_state"])
        model.to(device)
        return MIOFlowBenchmarkModel(
            model=model,
            pca=payload["pca"],
            mean=payload["mean"],
            std=payload["std"],
            time_to_idx={float(k): float(v) for k, v in payload["time_to_idx"].items()},
            device=device,
        )

    X_train = _to_dense_float32(adata_train.X)
    n_components = min(int(cfg.get("pca_dims", 50)), X_train.shape[0] - 1, X_train.shape[1])
    if n_components < 2:
        raise ValueError(f"Not enough data/features for PCA: n_components={n_components}")
    pca = PCA(n_components=n_components, svd_solver="randomized", random_state=seed)
    X_pca = pca.fit_transform(X_train).astype(np.float32)

    train_tps = adata_train.obs[time_key].astype(float).to_numpy()
    time_to_idx = {float(t): float(i) for i, t in enumerate(sorted(float(t) for t in full_times))}

    mean = X_pca.mean(axis=0)
    std = np.where(X_pca.std(axis=0) == 0, 1.0, X_pca.std(axis=0))
    X_norm = (X_pca - mean) / std
    time_series_data = []
    for t in sorted(set(float(t) for t in train_tps)):
        mask = np.isclose(train_tps, t)
        time_series_data.append((X_norm[mask].astype(np.float32), time_to_idx[float(t)]))
    dataset = TimeSeriesDataset(time_series_data)

    hidden_dim = int(cfg.get("hidden_dim", 64))
    momentum_beta = float(cfg.get("momentum_beta", 0.0))
    model = ODEFunc(input_dim=n_components, hidden_dim=hidden_dim, momentum_beta=momentum_beta)

    train_mioflow(
        model=model,
        dataset=dataset,
        num_epochs=int(cfg.get("n_epochs", 40)),
        batch_size=cfg.get("batch_size", 64),
        learning_rate=float(cfg.get("learning_rate", 1e-3)),
        device=device,
        lambda_ot=float(cfg.get("lambda_ot", 1.0)),
        lambda_density=float(cfg.get("lambda_density", 0.0))
        if bool(cfg.get("use_density_loss", False)) else 0.0,
        lambda_energy=float(cfg.get("lambda_energy", 0.01)),
        energy_time_steps=int(cfg.get("energy_time_steps", 10)),
        scheduler_type=cfg.get("scheduler_type"),
    )
    model.to(device)

    torch.save({
        "model_state": model.state_dict(),
        "pca": pca,
        "mean": mean,
        "std": std,
        "time_to_idx": {str(k): float(v) for k, v in time_to_idx.items()},
        "input_dim": n_components,
        "hidden_dim": hidden_dim,
        "momentum_beta": momentum_beta,
    }, cache_path)

    return MIOFlowBenchmarkModel(model, pca, mean, std, time_to_idx, device)


def run_outputs(model: MIOFlowBenchmarkModel, adata_full, adata_train, time_key: str,
                cell_state_key: str, full_times: list[float], train_times: list[float],
                heldout_times: list[float], cfg: dict, output_dir: Path) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    n_sim_cells = int(cfg.get("n_sim_cells", 2000))
    metric_sample_cells = int(cfg.get("metric_sample_cells", 1000))
    rng = np.random.default_rng(int(cfg.get("seed", 0)))

    X_full = _to_dense_float32(adata_full.X)
    X_train = _to_dense_float32(adata_train.X)
    X_full_pca = model.pca.transform(X_full).astype(np.float32)
    X_train_pca = model.pca.transform(X_train).astype(np.float32)
    tps_full = adata_full.obs[time_key].astype(float).to_numpy()
    tps_train = adata_train.obs[time_key].astype(float).to_numpy()
    states_full = adata_full.obs[cell_state_key].astype(str).to_numpy()
    states_train = adata_train.obs[cell_state_key].astype(str).to_numpy()

    eval_tps = heldout_times if heldout_times else full_times
    source_for_eval = {}
    for target in eval_tps:
        if heldout_times:
            prev = [t for t in train_times if float(t) <= float(target)]
            source_for_eval[float(target)] = float(prev[-1] if prev else train_times[0])
        else:
            source_for_eval[float(target)] = float(train_times[0])

    pred_expr_blocks = []
    pred_pca_blocks = []
    per_tp_rows = []
    for i, target in enumerate(eval_tps):
        source = source_for_eval[float(target)]
        source_mask = np.isclose(tps_full, source)
        source_pca = X_full_pca[source_mask]
        source_pca = _sample_rows(source_pca, n_sim_cells, seed=10_000 + i, replace=True)
        pred_pca = model.project_pca(source_pca, source, float(target))
        pred_expr = model.pca_to_expression(pred_pca)
        pred_pca_blocks.append(pred_pca)
        pred_expr_blocks.append(pred_expr)

        obs_expr = X_full[np.isclose(tps_full, float(target))]
        pred_metric = _sample_rows(pred_expr, metric_sample_cells, seed=20_000 + i)
        obs_metric = _sample_rows(obs_expr, metric_sample_cells, seed=30_000 + i)
        metrics = _forecast_metrics(pred_metric, obs_metric)
        per_tp_rows.append({
            "timepoint": float(target),
            "source_timepoint": source,
            **metrics,
            "n_pred_cells": int(pred_expr.shape[0]),
            "n_obs_cells": int(obs_expr.shape[0]),
            "metric_sample_pred_cells": int(pred_metric.shape[0]),
            "metric_sample_obs_cells": int(obs_metric.shape[0]),
        })

    proj_expr = np.stack(pred_expr_blocks, axis=0) if pred_expr_blocks else np.array([])
    proj_emb = np.vstack(pred_pca_blocks) if pred_pca_blocks else np.array([])
    projected_expression_path = output_dir / "projected_expression.npy"
    projected_embedding_path = output_dir / "projected_embedding.npy"
    np.save(projected_expression_path, proj_expr)
    np.save(projected_embedding_path, proj_emb)

    per_tp_df = pd.DataFrame(per_tp_rows)
    per_tp_forecast_path = output_dir / "per_timepoint_forecast_metrics.csv"
    per_tp_df.to_csv(per_tp_forecast_path, index=False)
    forecast_summary = {
        key: float(per_tp_df[key].mean()) if not per_tp_df.empty else None
        for key in [
            "wasserstein_distance",
            "gaussian_mmd",
            "energy_distance_mmd",
            "hausdorff_loss",
        ]
    }
    forecast_summary.update({
        "n_eval_timepoints": len(eval_tps),
        "eval_timepoints": [float(t) for t in eval_tps],
        "scenario_type": "forecast" if heldout_times else "reconstruction",
        "status": "completed",
    })
    forecast_metrics_path = output_dir / "forecast_metrics.json"
    with open(forecast_metrics_path, "w", encoding="utf-8") as f:
        json.dump(forecast_summary, f, indent=2)

    unique_states = sorted(set(states_full))
    centroids = np.vstack([
        X_full_pca[states_full == state].mean(axis=0) for state in unique_states
    ])
    clf = LogisticRegression(max_iter=500, class_weight="balanced", random_state=0)
    clf.fit(X_full_pca, states_full)

    label_rows = []
    ari_rows = []
    entropy_values = []
    offset = 0
    for i, target in enumerate(eval_tps):
        pred_pca = pred_pca_blocks[i]
        dists = pairwise_distances(pred_pca, centroids)
        nearest = np.argmin(dists, axis=1)
        pred_labels = np.array([unique_states[j] for j in nearest])
        probs = clf.predict_proba(pred_pca)
        if probs.shape[1] > 1:
            entropy = -np.sum(probs * np.log(probs + 1e-12), axis=1) / np.log(probs.shape[1])
            entropy_values.extend(entropy.tolist())

        obs_mask = np.isclose(tps_full, float(target))
        obs_pca = X_full_pca[obs_mask]
        obs_labels = states_full[obs_mask]
        if obs_pca.shape[0] > 0:
            nn = np.argmin(pairwise_distances(pred_pca, obs_pca), axis=1)
            ref_labels = obs_labels[nn]
            ari = float(adjusted_rand_score(ref_labels, pred_labels))
            ari_rows.append({
                "timepoint": float(target),
                "adjusted_rand_index": ari,
                "n_projected_cells": int(pred_pca.shape[0]),
                "n_observed_reference_cells": int(obs_pca.shape[0]),
            })
        for j, label in enumerate(pred_labels):
            label_rows.append({
                "projected_timepoint": float(target),
                "cell_idx": int(offset + j),
                "projected_cluster_label": label,
            })
        offset += pred_pca.shape[0]

    labels_path = output_dir / "projected_cluster_labels.csv"
    pd.DataFrame(label_rows).to_csv(labels_path, index=False)
    per_tp_embedding_path = output_dir / "per_timepoint_embedding_metrics.csv"
    pd.DataFrame(ari_rows).to_csv(per_tp_embedding_path, index=False)
    embedding_metrics = {
        "adjusted_rand_index": float(np.mean([r["adjusted_rand_index"] for r in ari_rows]))
        if ari_rows else None,
        "avg_normalized_classifier_entropy": float(np.mean(entropy_values))
        if entropy_values else None,
        "projected_latent_dim": int(proj_emb.shape[1]) if proj_emb.size else None,
        "n_eval_timepoints": len(eval_tps),
        "eval_timepoints": [float(t) for t in eval_tps],
        "status": "completed" if ari_rows and entropy_values else "embedding_metrics_incomplete",
        "latent_space": "PCA",
    }
    embedding_metrics_path = output_dir / "embedding_metrics.json"
    with open(embedding_metrics_path, "w", encoding="utf-8") as f:
        json.dump(embedding_metrics, f, indent=2)

    train_states = sorted(set(states_train))
    state_to_idx = {state: i for i, state in enumerate(train_states)}
    train_centroids = np.vstack([
        X_train_pca[states_train == state].mean(axis=0) for state in train_states
    ])
    stm_raw = np.zeros((len(train_states), len(train_states)), dtype=float)
    stm_counts = np.zeros(len(train_states), dtype=float)
    for source, target in zip(train_times[:-1], train_times[1:]):
        mask = np.isclose(tps_train, float(source))
        if not mask.any():
            continue
        src_pca = X_train_pca[mask]
        src_states = states_train[mask]
        pred_next = model.project_pca(src_pca, float(source), float(target))
        weights = _safe_softmax_neg_distance(pairwise_distances(pred_next, train_centroids))
        for row_i, src_state in enumerate(src_states):
            src_idx = state_to_idx[src_state]
            stm_raw[src_idx] += weights[row_i]
            stm_counts[src_idx] += 1

    stm_norm = np.zeros_like(stm_raw)
    for i in range(len(train_states)):
        if stm_counts[i] > 0:
            stm_norm[i] = stm_raw[i] / stm_counts[i]
    stm_path = output_dir / "state_transition_matrix.csv"
    pd.DataFrame(
        stm_norm,
        index=pd.Index(train_states, name="source_state"),
        columns=train_states,
    ).to_csv(stm_path)
    edge_rows = [
        {"source_state": s, "target_state": t, "weight": float(stm_norm[state_to_idx[s], state_to_idx[t]])}
        for s in train_states
        for t in train_states
        if stm_norm[state_to_idx[s], state_to_idx[t]] > 1e-8
    ]
    edges_path = output_dir / "lineage_graph_edges.csv"
    pd.DataFrame(edge_rows).to_csv(edges_path, index=False)
    with open(output_dir / "lineage_diagnostics.json", "w", encoding="utf-8") as f:
        json.dump({
            "source_transition_counts": {
                state: int(stm_counts[state_to_idx[state]]) for state in train_states
            },
            "latent_space": "PCA",
        }, f, indent=2)

    return {
        "projected_expression": str(projected_expression_path),
        "projected_embedding": str(projected_embedding_path),
        "projected_cluster_labels": str(labels_path),
        "forecast_metrics": str(forecast_metrics_path),
        "per_timepoint_forecast_metrics": str(per_tp_forecast_path),
        "embedding_metrics": str(embedding_metrics_path),
        "per_timepoint_embedding_metrics": str(per_tp_embedding_path),
        "state_transition_matrix": str(stm_path),
        "lineage_graph_edges": str(edges_path),
    }


def run_pipeline(adata, scenario_id: str, scenario_config: dict, output_dir) -> dict:
    t0 = time.time()
    status = "failed"
    notes = ""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    time_key = scenario_config.get("time_key", "abs_day")
    cell_state_key = scenario_config.get("cell_state_key", "final_milestone_label_coarse")
    dataset_id = scenario_config.get("dataset_id", "GSE230659")
    cfg = scenario_config.get("mioflow_params", {}) or {}
    scenario_params = scenario_config.get("scenario_params", {}) or {}

    try:
        adata_full = _materialize(adata)
        if time_key not in adata_full.obs.columns:
            raise KeyError(f"time_key not found in adata.obs: {time_key}")
        if cell_state_key not in adata_full.obs.columns:
            raise KeyError(f"cell_state_key not found in adata.obs: {cell_state_key}")

        all_times = sorted(float(t) for t in adata_full.obs[time_key].astype(float).unique())
        train_times = [float(t) for t in scenario_params.get("train_times", [])]
        if not train_times:
            train_times = all_times
        heldout_times = [float(t) for t in scenario_params.get("heldout_times", [])]
        if not heldout_times and train_times != all_times:
            train_set = set(train_times)
            heldout_times = [t for t in all_times if t not in train_set]

        train_set = set(train_times)
        train_mask = adata_full.obs[time_key].astype(float).isin(train_set).values
        adata_train = adata_full[train_mask].copy()
        subsample = cfg.get("subsample_per_timepoint")
        if subsample:
            keep = []
            rng = np.random.default_rng(int(cfg.get("seed", 0)))
            train_tps = adata_train.obs[time_key].astype(float).to_numpy()
            for tp in sorted(set(train_tps)):
                idx = np.where(np.isclose(train_tps, tp))[0]
                if len(idx) > int(subsample):
                    idx = rng.choice(idx, size=int(subsample), replace=False)
                keep.extend(idx.tolist())
            adata_train = adata_train[np.sort(keep)].copy()

        print(
            f"[MIOFlow] train_cells={adata_train.n_obs}, full_cells={adata_full.n_obs}, "
            f"train_times={train_times}, heldout_times={heldout_times}"
        )
        model = train_or_load(adata_train, all_times, time_key, cfg, output_dir)
        result = run_outputs(
            model=model,
            adata_full=adata_full,
            adata_train=adata_train,
            time_key=time_key,
            cell_state_key=cell_state_key,
            full_times=all_times,
            train_times=train_times,
            heldout_times=heldout_times,
            cfg=cfg,
            output_dir=output_dir,
        )
        status = "completed"
        return result
    except Exception:
        notes = traceback.format_exc()
        print(f"[MIOFlow] ERROR:\n{notes}")
        raise
    finally:
        metadata = {
            "method": "mioflow",
            "dataset": dataset_id,
            "scenario": scenario_id,
            "result_class": scenario_config.get("result_class", "unknown"),
            "formal_benchmark": bool(scenario_config.get("formal_benchmark", False)),
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
            "mioflow_params": cfg,
            "notes": notes,
        }
        with open(output_dir / "run_metadata.json", "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run MIOFlow benchmark pipeline")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    import anndata
    import yaml

    with open(args.config, encoding="utf-8-sig") as f:
        cfg = yaml.safe_load(f)
    root = _find_project_root()
    dataset_cfg = cfg.get("dataset", {})
    adata_path = root / dataset_cfg["h5ad_path"]
    adata = anndata.read_h5ad(adata_path)
    scenario_cfg = {}
    scenario_cfg.update(cfg)
    scenario_cfg.update({
        "dataset_id": dataset_cfg.get("id", "GSE230659"),
        "time_key": dataset_cfg.get("time_key", "abs_day"),
        "cell_state_key": (cfg.get("lineage") or {}).get(
            "cell_state_key", "final_milestone_label_coarse"
        ),
    })
    run_pipeline(
        adata=adata,
        scenario_id=str(cfg.get("scenario", "A")),
        scenario_config=scenario_cfg,
        output_dir=root / cfg["output"]["base_dir"],
    )


if __name__ == "__main__":
    main()
