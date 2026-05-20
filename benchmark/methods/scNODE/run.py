"""
benchmark/methods/scNODE/run.py
scNODE standalone runner for the scTimeBench-aligned benchmark.
Framework reference: experimental framework v2.md 搂10.3, 搂14 Steps 5-7

Capability flags (per v2 搂6):
    scNODE: supports_unseen_timepoint_projection = True
            supports_lineage_inference            = True
    鈫?evaluated on ALL three benchmark dimensions:
        Forecast Accuracy, Embedding Coherence, Lineage Fidelity

Architecture note:
    scNODE is a VAE + Neural ODE generative model (rsinghlab/scNODE).
    All temporal predictions are rooted at t0 (first observed training
    timepoint, day 0.5 for GSE230659). The ODE propagates the latent
    trajectory from t0 forward; held-out cells are NEVER used as model
    input 鈥?they are only used as ground truth when computing forecast
    metrics.

Usage:
    conda activate traj_env
    # From project root (or set TRAJ_PROJECT_ROOT=/path/to/project):
    python benchmark/methods/scNODE/run.py \
        --config benchmark/configs/scnode_gse230659_marker_fm_silver_A_hvg2000_formal.yaml

    # On HPC (SGE/SLURM), export TRAJ_PROJECT_ROOT first, then:
    python benchmark/methods/scNODE/run.py \
        --config benchmark/configs/scnode_gse230659_marker_fm_silver_A_hvg2000_formal.yaml

    Or via train_and_test.sh (called from benchmark/):
        bash benchmark/methods/scNODE/train_and_test.sh \
            benchmark/configs/scnode_gse230659_marker_fm_silver_A_hvg2000_formal.yaml

What this script does:
    1. Reads the YAML config.
    2. Loads adata (backed='r' for Scenarios B/C; materialises after filter).
    3. Applies scenario_params.train_times filter if specified.
    4. Trains scNODE (or loads cached trained_scnode_model.pth).
    5. Forecast Accuracy:  predict GEX at held-out tps; compute metrics.
    6. Embedding Coherence: ODE latent at eval tps; state centroid assignment.
    7. Lineage Fidelity:   per-pair batch next-tp ODE 鈫?STM 鈫?eval_lineage.
    8. Writes run_metadata.json.
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
from scipy.sparse import issparse

# ---------------------------------------------------------------------------
# Project-root discovery (same logic as WOT/run.py)
# ---------------------------------------------------------------------------

def _find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "benchmark").is_dir() and (parent / "data").is_dir():
            return parent
    raise RuntimeError(
        "Cannot locate project root. "
        "Set the TRAJ_PROJECT_ROOT environment variable."
    )

# ---------------------------------------------------------------------------
# scNODE module import (via local submodule path)
# ---------------------------------------------------------------------------

_SCNODE_MODULE_DIR = Path(__file__).parent / "scNODE_module"
if _SCNODE_MODULE_DIR.exists() and str(_SCNODE_MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(_SCNODE_MODULE_DIR))

from optim.running import constructscNODEModel, scNODETrainWithPreTrain  # noqa: E402
from optim.evaluation import globalEvaluation  # noqa: E402

# ---------------------------------------------------------------------------
# Defaults (match scTimeBench's run.py and original benchmark scripts)
# ---------------------------------------------------------------------------

_DEFAULTS = dict(
    latent_dim=50,
    drift_latent_size=[50, 50],
    enc_latent_list=[64, 64],
    dec_latent_list=[64, 64],
    pretrain_iters=200,
    epochs=10,
    iters=100,
    batch_size=32,
    lr=1e-3,
    latent_coeff=1.0,
    kl_coeff=0.0,
    seed=42,
    n_sim_cells=None,          # None 鈫?use actual cell count at t0 (capped at 2000)
    n_sim_cells_cap=2000,      # safety cap for memory-limited machines
    metric_sample_cells=1000,  # cap observed/predicted cells for OT metrics
)

# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

def prepare_data(adata, time_key: str, train_times=None):
    """
    Convert adata to scNODE training format.

    Parameters
    ----------
    adata       : AnnData (already materialised, not backed)
    time_key    : obs column holding the real-valued timepoint (e.g. "abs_day")
    train_times : list of float or None
        If None, all unique timepoints in adata are treated as training tps.

    Returns
    -------
    train_data      : list[torch.FloatTensor]  鈥?one tensor per training tp
    train_tps       : torch.FloatTensor        鈥?timepoint values (float)
    all_unique_tps  : list[float]              鈥?all tps present in adata
    heldout_tps     : list[float]              鈥?tps in adata but not in train_times
    data_np         : np.ndarray               鈥?full expression matrix (n_cells, n_genes)
    cell_tps_np     : np.ndarray               鈥?per-cell timepoint values
    """
    data_np = adata.X.toarray() if issparse(adata.X) else np.array(adata.X)
    cell_tps_np = adata.obs[time_key].values.astype(float)
    all_unique_tps = sorted(np.unique(cell_tps_np).tolist())

    if train_times is None:
        train_times = all_unique_tps
    train_times_f = sorted(float(t) for t in train_times)
    train_set = set(train_times_f)
    heldout_tps = [t for t in all_unique_tps if t not in train_set]

    train_data = []
    for t in train_times_f:
        mask = np.isclose(cell_tps_np, t)
        tp_data = data_np[mask, :]
        train_data.append(torch.FloatTensor(tp_data))
        print(f"  [data] tp={t:.2f}: {tp_data.shape[0]} cells")
    train_tps = torch.FloatTensor(train_times_f)

    return train_data, train_tps, all_unique_tps, heldout_tps, data_np, cell_tps_np

# ---------------------------------------------------------------------------
# Model construction
# ---------------------------------------------------------------------------

def build_model(n_genes: int, scnode_cfg: dict):
    """Instantiate an scNODE model from config (or defaults)."""
    def g(k):
        return scnode_cfg.get(k, _DEFAULTS[k])

    return constructscNODEModel(
        n_genes,
        latent_dim=g("latent_dim"),
        enc_latent_list=g("enc_latent_list"),
        dec_latent_list=g("dec_latent_list"),
        drift_latent_size=g("drift_latent_size"),
        latent_enc_act="none",
        latent_dec_act="relu",
        drift_act="relu",
        ode_method="euler",
    )

# ---------------------------------------------------------------------------
# Training
# ---------------------------------------------------------------------------

def train_or_load(train_data, train_tps, n_genes: int,
                  scnode_cfg: dict, cache_path: Path):
    """
    Train scNODE or load from cache.

    Cache behaviour: if trained_scnode_model.pth exists in output_dir,
    model weights are loaded and training is skipped.  This matches the
    scTimeBench run.py cache pattern and is critical for HPC reruns.
    """
    model = build_model(n_genes, scnode_cfg)

    if cache_path.exists():
        print(f"[scNODE] Loading cached model from {cache_path}")
        state = torch.load(str(cache_path), map_location="cpu")
        model.load_state_dict(state)
        model.eval()
        return model

    def g(k):
        return scnode_cfg.get(k, _DEFAULTS[k])

    seed = g("seed")
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)

    print(f"[scNODE] Training: pretrain_iters={g('pretrain_iters')}, "
          f"epochs={g('epochs')}, seed={seed}")

    model, loss_list, _, _, _ = scNODETrainWithPreTrain(
        train_data, train_tps, model,
        latent_coeff=g("latent_coeff"),
        epochs=g("epochs"),
        iters=g("iters"),
        batch_size=g("batch_size"),
        lr=g("lr"),
        kl_coeff=g("kl_coeff"),
        pretrain_iters=g("pretrain_iters"),
        pretrain_lr=g("lr"),
    )

    torch.save(model.state_dict(), str(cache_path))
    print(f"[scNODE] Model saved 鈫?{cache_path}")
    model.eval()
    return model

# ---------------------------------------------------------------------------
# Forecast Accuracy
# ---------------------------------------------------------------------------

def _sample_rows_for_metrics(X: np.ndarray, max_cells: int | None, seed: int) -> np.ndarray:
    """Deterministically cap rows before pairwise OT-style metric computation."""
    if max_cells is None or int(max_cells) <= 0 or X.shape[0] <= int(max_cells):
        return X
    rng = np.random.default_rng(seed)
    idx = rng.choice(X.shape[0], size=int(max_cells), replace=False)
    return X[idx]

def run_forecast_accuracy(
    model,
    adata_full,           # FULL adata (all tps) for observed-cell access
    time_key: str,
    all_unique_tps: list,
    heldout_tps: list,
    n_sim_cells: int,
    output_dir: Path,
    metric_sample_cells: int = 1000,
    seed: int = 42,
):
    """
    Predict GEX at held-out (or all) timepoints and compute forecast metrics.

    Strategy
    --------
    - Use model.predict(first_tp_cells, all_tps, n_cells) once.
    - Extract recon_obs[:, tp_idx, :] for each eval timepoint.
    - Compare predicted distribution against observed cells using
      scNODE's own globalEvaluation (L2, cosine, correlation, OT/Wasserstein).
    - If heldout_tps is empty (Scenario A with no holdout), evaluate on
      ALL training timepoints as a within-distribution reconstruction check.

    Outputs
    -------
    projected_expression.npy          shape (n_eval_tps, n_sim_cells, n_genes)
    forecast_metrics.json             mean over eval tps
    per_timepoint_forecast_metrics.csv one row per eval tp
    """
    model.eval()
    data_full = (adata_full.X.toarray() if issparse(adata_full.X)
                 else np.array(adata_full.X))
    cell_tps_full = adata_full.obs[time_key].values.astype(float)

    # Cells at first training timepoint (t0 = day 0.5 for GSE230659)
    t0 = all_unique_tps[0]
    t0_mask = np.isclose(cell_tps_full, t0)
    first_tp_data = torch.FloatTensor(data_full[t0_mask, :])

    eval_tps = heldout_tps if heldout_tps else all_unique_tps

    # Build the full ODE integration axis: union of train tps + held-out tps.
    # For Scenario A (no holdout) this equals all_unique_tps.
    # For Scenarios B/C the ODE must also integrate to held-out time values so
    # that recon_obs[:, tp_idx, :] is available for every held-out timepoint.
    all_eval_tps = sorted(set(list(all_unique_tps) + list(heldout_tps)))
    all_eval_tps_tensor = torch.FloatTensor(all_eval_tps)

    print(f"[Forecast] Predicting from t0={t0} across {len(all_eval_tps)} tps "
          f"(train={len(all_unique_tps)}, heldout={len(heldout_tps)}, "
          f"n_sim_cells={n_sim_cells}, metric_sample_cells={metric_sample_cells}) ...")
    with torch.no_grad():
        _, _, recon_obs = model.predict(first_tp_data, all_eval_tps_tensor, n_cells=n_sim_cells)
    recon_np = recon_obs.detach().numpy()   # (n_sim_cells, len(all_eval_tps), n_genes)

    per_tp_rows, all_pred = [], []
    for t_h in eval_tps:
        if t_h not in all_eval_tps:
            print(f"  [Forecast] WARNING: t={t_h} not in all_eval_tps, skipping")
            continue
        tp_idx = all_eval_tps.index(t_h)
        X_pred = recon_np[:, tp_idx, :]            # (n_sim_cells, n_genes)

        obs_mask = np.isclose(cell_tps_full, t_h)
        X_obs = data_full[obs_mask, :]             # (n_obs_cells, n_genes)

        if X_obs.shape[0] == 0:
            print(f"  [Forecast] t={t_h}: no observed cells, skipping")
            continue

        X_obs_metric = _sample_rows_for_metrics(
            X_obs, metric_sample_cells, seed=seed + 10_000 + len(per_tp_rows)
        )
        X_pred_metric = _sample_rows_for_metrics(
            X_pred, metric_sample_cells, seed=seed + 20_000 + len(per_tp_rows)
        )

        metrics_raw = globalEvaluation(X_obs_metric, X_pred_metric)
        per_tp_rows.append({
            "timepoint": t_h,
            "wasserstein_ot": metrics_raw["ot"],
            "l2_dist":        metrics_raw["l2"],
            "cosine_dist":    metrics_raw["cos"],
            "correlation_dist": metrics_raw["corr"],
            "n_obs_cells":    int(X_obs.shape[0]),
            "n_pred_cells":   int(X_pred.shape[0]),
            "metric_sample_obs_cells": int(X_obs_metric.shape[0]),
            "metric_sample_pred_cells": int(X_pred_metric.shape[0]),
        })
        all_pred.append(X_pred)
        print(f"  [Forecast] t={t_h:.2f}: Wasserstein={metrics_raw['ot']:.4f} "
              f"L2={metrics_raw['l2']:.4f} "
              f"(metric_obs={X_obs_metric.shape[0]}, metric_pred={X_pred_metric.shape[0]})")

    # Save projected_expression.npy
    proj_arr = np.stack(all_pred, axis=0) if all_pred else np.array([])
    proj_path = output_dir / "projected_expression.npy"
    np.save(str(proj_path), proj_arr)

    # Save per-timepoint CSV
    per_tp_df = pd.DataFrame(per_tp_rows)
    per_tp_path = output_dir / "per_timepoint_forecast_metrics.csv"
    per_tp_df.to_csv(str(per_tp_path), index=False)

    # Summary JSON  (Wasserstein is the primary framework metric; others deferred)
    mean_ot = float(np.mean([r["wasserstein_ot"] for r in per_tp_rows])) if per_tp_rows else None
    summary = {
        "wasserstein_distance": mean_ot,
        "gaussian_mmd":        None,   # implement in eval_forecast.py
        "energy_distance_mmd": None,   # implement in eval_forecast.py
        "hausdorff_loss":      None,   # implement in eval_forecast.py
        "n_eval_timepoints":   len(per_tp_rows),
        "eval_timepoints":     eval_tps,
        "scenario_type":       "forecast" if heldout_tps else "reconstruction",
        "metric_sample_cells":  int(metric_sample_cells),
        "status":              "completed",
    }
    metrics_path = output_dir / "forecast_metrics.json"
    with open(str(metrics_path), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    print(f"[Forecast] Done. Mean Wasserstein={mean_ot}")
    return str(proj_path), str(metrics_path), str(per_tp_path)

# ---------------------------------------------------------------------------
# Embedding Coherence
# ---------------------------------------------------------------------------

def run_embedding_coherence(
    model,
    adata_full,
    time_key: str,
    all_unique_tps: list,
    heldout_tps: list,
    n_sim_cells: int,
    cell_state_key: str,
    output_dir: Path,
):
    """
    Produce VAE latent embeddings for projected cells at eval timepoints.

    Strategy
    --------
    - Use model.predict(first_tp_cells, all_tps, n_cells) to obtain
      latent_seq[:, tp_idx, :] 鈥?the ODE-propagated latent at each tp.
    - Assign each projected latent to the nearest observed state centroid
      (computed from vaeReconstruct on training cells).
    - ARI and classifier entropy are deferred (eval_embedding.py stubs).

    Outputs
    -------
    projected_embedding.npy       (n_eval_tps 脳 n_sim_cells, latent_dim)
    embedding_metrics.json        ARI=None, entropy=None (deferred)
    projected_cluster_labels.csv  projected_timepoint | cell_idx | cluster_label
    """
    from scipy.spatial.distance import cdist as scipy_cdist
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import adjusted_rand_score

    model.eval()
    data_full = (adata_full.X.toarray() if issparse(adata_full.X)
                 else np.array(adata_full.X))
    cell_tps_full = adata_full.obs[time_key].values.astype(float)

    t0 = all_unique_tps[0]
    t0_mask = np.isclose(cell_tps_full, t0)
    first_tp_data = torch.FloatTensor(data_full[t0_mask, :])

    eval_tps = heldout_tps if heldout_tps else all_unique_tps

    # Compute state centroids from all observed cells in adata_full
    observed_state_labels = []
    observed_centroids = {}
    if cell_state_key and cell_state_key in adata_full.obs.columns:
        state_labels_all = adata_full.obs[cell_state_key].values
        with torch.no_grad():
            latent_list, _ = model.vaeReconstruct([data_full])
        all_latents_np = latent_list[0].detach().numpy()
        observed_state_labels = sorted(set(state_labels_all))
        for state in observed_state_labels:
            mask = state_labels_all == state
            if mask.sum() > 0:
                observed_centroids[state] = all_latents_np[mask].mean(axis=0)
        print(f"[Embedding] Computed centroids for {len(observed_centroids)} states.")
    else:
        state_labels_all = None
        all_latents_np = None

    clf = None
    if state_labels_all is not None and all_latents_np is not None:
        try:
            clf = LogisticRegression(
                max_iter=500,
                class_weight="balanced",
                random_state=0,
            )
            clf.fit(all_latents_np, state_labels_all)
        except Exception as exc:
            print(f"[Embedding] Classifier entropy skipped: {exc}")
            clf = None

    # ODE latent trajectory from t0.
    # Must include held-out tps in the integration axis so latent_seq[:, tp_idx, :]
    # is available for every held-out timepoint (same fix as run_forecast_accuracy).
    all_eval_tps = sorted(set(list(all_unique_tps) + list(heldout_tps)))
    all_eval_tps_tensor = torch.FloatTensor(all_eval_tps)
    with torch.no_grad():
        _, latent_seq, _ = model.predict(first_tp_data, all_eval_tps_tensor, n_cells=n_sim_cells)
    latent_seq_np = latent_seq.detach().numpy()   # (n_sim_cells, len(all_eval_tps), latent_dim)

    all_pred_latent, all_label_rows = [], []
    ari_rows = []
    entropy_values = []
    for t_h in eval_tps:
        if t_h not in all_eval_tps:
            continue
        tp_idx = all_eval_tps.index(t_h)
        proj_latent = latent_seq_np[:, tp_idx, :]   # (n_sim_cells, latent_dim)
        all_pred_latent.append(proj_latent)

        if observed_centroids:
            centroid_mat = np.stack(
                [observed_centroids[s] for s in observed_state_labels])
            dists = scipy_cdist(proj_latent, centroid_mat, metric="euclidean")
            nearest_idx = np.argmin(dists, axis=1)
            labels = [observed_state_labels[i] for i in nearest_idx]
        else:
            labels = ["unknown"] * proj_latent.shape[0]

        if (
            state_labels_all is not None
            and all_latents_np is not None
            and np.isclose(cell_tps_full, float(t_h)).any()
        ):
            obs_mask_t = np.isclose(cell_tps_full, float(t_h))
            obs_latent_t = all_latents_np[obs_mask_t]
            obs_labels_t = state_labels_all[obs_mask_t]
            nn_dists = scipy_cdist(proj_latent, obs_latent_t, metric="euclidean")
            nn_idx = np.argmin(nn_dists, axis=1)
            reference_labels = obs_labels_t[nn_idx]
            ari = float(adjusted_rand_score(reference_labels, labels))
            ari_rows.append({
                "timepoint": float(t_h),
                "adjusted_rand_index": ari,
                "n_projected_cells": int(proj_latent.shape[0]),
                "n_observed_reference_cells": int(obs_latent_t.shape[0]),
            })

        if clf is not None:
            probs = clf.predict_proba(proj_latent)
            if probs.shape[1] > 1:
                eps = 1e-12
                entropy = -np.sum(probs * np.log(probs + eps), axis=1)
                entropy = entropy / np.log(probs.shape[1])
                entropy_values.extend(entropy.tolist())

        for i, lbl in enumerate(labels):
            all_label_rows.append({
                "projected_timepoint": t_h,
                "cell_idx": i,
                "projected_cluster_label": lbl,
            })

    # Save outputs
    emb_arr = np.vstack(all_pred_latent) if all_pred_latent else np.array([])
    emb_path = output_dir / "projected_embedding.npy"
    np.save(str(emb_path), emb_arr)

    labels_df = pd.DataFrame(all_label_rows)
    labels_path = output_dir / "projected_cluster_labels.csv"
    labels_df.to_csv(str(labels_path), index=False)

    per_tp_embedding_path = output_dir / "per_timepoint_embedding_metrics.csv"
    pd.DataFrame(ari_rows).to_csv(str(per_tp_embedding_path), index=False)

    mean_ari = (
        float(np.mean([row["adjusted_rand_index"] for row in ari_rows]))
        if ari_rows else None
    )
    mean_entropy = float(np.mean(entropy_values)) if entropy_values else None

    emb_metrics = {
        "adjusted_rand_index":             mean_ari,
        "avg_normalized_classifier_entropy": mean_entropy,
        "projected_latent_dim":            int(latent_seq_np.shape[2]),
        "n_eval_timepoints":               len(eval_tps),
        "eval_timepoints":                 eval_tps,
        "status":                          (
            "completed" if mean_ari is not None and mean_entropy is not None
            else "embedding_metrics_incomplete"
        ),
        "ari_reference": (
            "nearest observed cell at the same evaluation timepoint in "
            "scNODE VAE latent space"
        ),
        "entropy_reference": (
            "balanced logistic regression trained on observed scNODE VAE "
            "latent embeddings and provider state labels"
        ),
    }
    emb_metrics_path = output_dir / "embedding_metrics.json"
    with open(str(emb_metrics_path), "w", encoding="utf-8") as f:
        json.dump(emb_metrics, f, indent=2)

    print(f"[Embedding] projected_embedding.npy shape: {emb_arr.shape}")
    return str(emb_path), str(emb_metrics_path), str(labels_path)

# ---------------------------------------------------------------------------
# Lineage Fidelity
# ---------------------------------------------------------------------------

def run_lineage_fidelity(
    model,
    adata_train,          # TRAINING cells only (post-filter)
    time_key: str,
    train_unique_tps: list,
    cell_state_key: str,
    output_dir: Path,
):
    """
    Build a state-level transition matrix from batch next-timepoint ODE predictions.

    Algorithm (per test plan 搂4.7)
    --------------------------------
    For each consecutive training pair (t_i, t_{i+1}):
      1. Batch-predict latent embeddings at t_{i+1} for all cells at t_i.
         model.predict(cells_at_t_i, [t_i, t_{i+1}], n_cells)
         鈫?latent_seq[:, 1, :] is the predicted latent at t_{i+1}.
      2. Soft-assign each predicted latent to observed state centroids via
         softmax over inverse Euclidean distances (temperature T=1.0).
      3. Accumulate weighted transitions to the source-state row.
    Row-normalise to get a stochastic STM.

    This is more efficient than per-cell single-step prediction
    (O(n_consecutive_pairs) model calls vs O(n_cells)).

    Outputs
    -------
    state_transition_matrix.csv     n_states 脳 n_states, source_state index
    lineage_graph_edges.csv         source_state | target_state | weight
    """
    from scipy.spatial.distance import cdist as scipy_cdist

    model.eval()
    data_train = (adata_train.X.toarray() if issparse(adata_train.X)
                  else np.array(adata_train.X))
    cell_tps = adata_train.obs[time_key].values.astype(float)

    if cell_state_key not in adata_train.obs.columns:
        raise KeyError(
            f"[Lineage] cell_state_key '{cell_state_key}' not found in adata.obs. "
            f"Available columns: {list(adata_train.obs.columns)}"
        )
    state_labels = adata_train.obs[cell_state_key].values
    unique_states = sorted(set(state_labels))
    n_states = len(unique_states)
    state_to_idx = {s: i for i, s in enumerate(unique_states)}
    print(f"[Lineage] {n_states} states, {len(train_unique_tps)} training tps")

    # Step 1: State centroids from observed cells (VAE latent space)
    print("[Lineage] Computing state centroids ...")
    with torch.no_grad():
        latent_list, _ = model.vaeReconstruct([data_train])
    all_latents_np = latent_list[0].detach().numpy()   # (n_cells, latent_dim)
    state_centroids = np.zeros((n_states, all_latents_np.shape[1]))
    for state, idx in state_to_idx.items():
        mask = state_labels == state
        if mask.sum() > 0:
            state_centroids[idx] = all_latents_np[mask].mean(axis=0)

    # Step 2 & 3: Batch next-tp predictions, accumulate STM
    sorted_tps = sorted(train_unique_tps)
    stm_raw    = np.zeros((n_states, n_states))
    stm_counts = np.zeros(n_states)
    total_fallback_count = 0
    total_skipped_assignment_count = 0

    for i in range(len(sorted_tps) - 1):
        t_cur = float(sorted_tps[i])
        t_nxt = float(sorted_tps[i + 1])

        cur_mask = np.isclose(cell_tps, t_cur)
        if cur_mask.sum() == 0:
            print(f"  [Lineage] t={t_cur}: no cells, skip")
            continue

        cur_data   = torch.FloatTensor(data_train[cur_mask, :])
        cur_states = state_labels[cur_mask]
        n_cur      = cur_data.shape[0]

        print(f"  [Lineage] {t_cur:.2f} 鈫?{t_nxt:.2f}: {n_cur} cells")
        with torch.no_grad():
            _, latent_seq, _ = model.predict(
                cur_data,
                torch.FloatTensor([t_cur, t_nxt]),
                n_cells=n_cur,
            )
        # latent_seq: (n_cur, 2, latent_dim)  [0]=t_cur, [1]=t_nxt
        pred_nxt = latent_seq[:, 1, :].detach().numpy()   # (n_cur, latent_dim)

        # Soft assignment: stable softmax over negative distances.
        # The naive exp(-distance) can underflow to all zeros on full-scale
        # runs, which then creates NaNs during row normalisation.
        T = 1.0
        dists   = scipy_cdist(pred_nxt, state_centroids, metric="euclidean")
        scores = -dists / T
        finite_scores = np.isfinite(scores)
        row_has_finite = finite_scores.any(axis=1)
        row_max = np.zeros((scores.shape[0], 1), dtype=float)
        if row_has_finite.any():
            row_max[row_has_finite, 0] = np.max(
                np.where(finite_scores[row_has_finite], scores[row_has_finite], -np.inf),
                axis=1,
            )
        shifted = np.where(finite_scores, scores - row_max, -np.inf)
        weights = np.exp(shifted)
        row_sums = weights.sum(axis=1, keepdims=True)
        good_rows = row_has_finite & np.isfinite(row_sums[:, 0]) & (row_sums[:, 0] > 0)
        weights = np.divide(
            weights,
            row_sums,
            out=np.zeros_like(weights),
            where=good_rows[:, None],
        )

        # If numerical underflow still leaves a row without mass, fall back to
        # a deterministic nearest-centroid assignment when possible. Rows with
        # entirely non-finite predictions remain all zero and are skipped below.
        bad_rows = np.where(~good_rows)[0]
        fallback_count = 0
        skipped_assignment_count = 0
        for bad_i in bad_rows:
            finite_d = np.isfinite(dists[bad_i])
            if finite_d.any():
                nearest = int(np.argmin(np.where(finite_d, dists[bad_i], np.inf)))
                weights[bad_i, nearest] = 1.0
                fallback_count += 1
            else:
                skipped_assignment_count += 1

        for cell_i, src_state in enumerate(cur_states):
            if weights[cell_i].sum() <= 0:
                continue
            src_idx = state_to_idx[src_state]
            stm_raw[src_idx]    += weights[cell_i]
            stm_counts[src_idx] += 1

        if fallback_count or skipped_assignment_count:
            total_fallback_count += fallback_count
            total_skipped_assignment_count += skipped_assignment_count
            print(
                "    [Lineage] soft-assignment fallback: "
                f"nearest={fallback_count}, skipped={skipped_assignment_count}"
            )

    # Row-normalise
    stm_norm = np.zeros_like(stm_raw)
    for idx in range(n_states):
        if stm_counts[idx] > 0:
            stm_norm[idx] = stm_raw[idx] / stm_counts[idx]
    stm_norm = np.nan_to_num(stm_norm, nan=0.0, posinf=0.0, neginf=0.0)
    unsupported_source_states = [
        state for state, idx in state_to_idx.items()
        if stm_counts[idx] == 0
    ]

    # Write state_transition_matrix.csv
    stm_df = pd.DataFrame(
        stm_norm,
        index=pd.Index(unique_states, name="source_state"),
        columns=unique_states,
    )
    stm_path = output_dir / "state_transition_matrix.csv"
    stm_df.to_csv(str(stm_path))

    # Write lineage_graph_edges.csv
    edge_rows = [
        {"source_state": s, "target_state": t,
         "weight": float(stm_norm[state_to_idx[s], state_to_idx[t]])}
        for s in unique_states
        for t in unique_states
        if stm_norm[state_to_idx[s], state_to_idx[t]] > 1e-8
    ]
    edges_df  = pd.DataFrame(edge_rows)
    edges_path = output_dir / "lineage_graph_edges.csv"
    edges_df.to_csv(str(edges_path), index=False)

    print(f"[Lineage] STM written  ({n_states}脳{n_states}): {stm_path}")
    print(f"[Lineage] Edges written ({len(edge_rows)} edges): {edges_path}")
    diagnostics = {
        "unsupported_source_states": unsupported_source_states,
        "n_unsupported_source_states": len(unsupported_source_states),
        "source_transition_counts": {
            state: int(stm_counts[idx]) for state, idx in state_to_idx.items()
        },
        "soft_assignment_nearest_fallback_count": int(total_fallback_count),
        "soft_assignment_skipped_count": int(total_skipped_assignment_count),
        "note": (
            "States with zero source_transition_counts are written as all-zero "
            "STM rows because no valid source-cell transitions were available "
            "for those states in the training timepoint pairs."
        ),
    }
    diag_path = output_dir / "lineage_diagnostics.json"
    with open(diag_path, "w", encoding="utf-8") as f:
        json.dump(diagnostics, f, indent=2)
    if unsupported_source_states:
        print(
            "[Lineage] Unsupported source states written as zero rows: "
            f"{unsupported_source_states}"
        )
    return str(stm_path), str(edges_path)

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="scNODE benchmark runner (GSE230659)")
    parser.add_argument(
        "--config", "--yaml_config", required=True,
        metavar="YAML",
        help="Path to YAML config (e.g. benchmark/configs/scnode_gse230659_marker_fm_silver_A_hvg2000_formal.yaml)"
    )
    args = parser.parse_args()

    import yaml
    import anndata

    with open(args.config, encoding="utf-8-sig") as f:
        cfg = yaml.safe_load(f)

    project_root = _find_project_root()
    # Add project root to path so eval_lineage can be imported
    if str(project_root) not in sys.path:
        sys.path.insert(0, str(project_root))

    # 鈹€鈹€ Parse config 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
    run_id        = cfg.get("run_id", "scnode_run")
    scenario_id   = str(cfg.get("scenario", "A"))
    dataset_cfg   = cfg.get("dataset", {})
    scnode_cfg    = cfg.get("scnode_params", {})
    lineage_cfg   = cfg.get("lineage", {})
    scenario_params = cfg.get("scenario_params", {})
    output_cfg    = cfg.get("output", {})

    time_key           = dataset_cfg.get("time_key", "abs_day")
    cell_state_key     = lineage_cfg.get("cell_state_key", "final_milestone_label_coarse")
    ref_graph_path     = lineage_cfg.get("reference_graph_path")
    edge_conf_mode     = lineage_cfg.get("edge_confidence_mode", "medium_and_above")
    excl_uncertain     = lineage_cfg.get("exclude_uncertain_states", False)

    train_times_raw = scenario_params.get("train_times")
    train_times     = [float(t) for t in train_times_raw] if train_times_raw else None
    heldout_times   = [float(t) for t in scenario_params.get("heldout_times", [])]

    n_sim_cells_cfg = scnode_cfg.get("n_sim_cells")   # None means auto
    metric_sample_cells = int(scnode_cfg.get("metric_sample_cells", _DEFAULTS["metric_sample_cells"]))
    seed = int(scnode_cfg.get("seed", _DEFAULTS["seed"]))

    h5ad_raw = dataset_cfg.get("h5ad_path", "")
    h5ad_path = (Path(h5ad_raw) if Path(h5ad_raw).is_absolute()
                 else project_root / h5ad_raw)

    out_base   = output_cfg.get("base_dir",
                                f"benchmark/results/scnode/scenario_{scenario_id}")
    output_dir = (Path(out_base) if Path(out_base).is_absolute()
                  else project_root / out_base)
    output_dir.mkdir(parents=True, exist_ok=True)

    model_cache = output_dir / "trained_scnode_model.pth"

    print(f"\n{'='*60}")
    print(f"[scNODE] run_id      : {run_id}")
    print(f"[scNODE] scenario    : {scenario_id}")
    print(f"[scNODE] h5ad        : {h5ad_path}")
    print(f"[scNODE] time_key    : {time_key}")
    print(f"[scNODE] state_key   : {cell_state_key}")
    print(f"[scNODE] train_times : {train_times}")
    print(f"[scNODE] heldout     : {heldout_times}")
    print(f"[scNODE] output_dir  : {output_dir}")
    print(f"{'='*60}\n")

    t_wall   = time.time()
    status   = "failed"
    err_notes = ""

    try:
        # 鈹€鈹€ Load adata 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        # Use backed='r' when a train-time filter will be applied immediately
        # (avoids materialising the full matrix, matching CellRank2 Adapter Step 0).
        has_filter = bool(train_times)
        backed     = "r" if has_filter else None
        print(f"[scNODE] Loading adata (backed={backed}) ...")
        adata = anndata.read_h5ad(str(h5ad_path), backed=backed)

        # 鈹€鈹€ Apply train-time filter 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        if has_filter:
            train_set = set(train_times)
            raw_tps   = adata.obs[time_key].values.astype(float)
            mask      = np.array([t in train_set for t in raw_tps])
            adata     = adata[mask].to_memory()
            print(f"[scNODE] After train filter: {adata.n_obs} cells 脳 {adata.n_vars} genes")
        else:
            # No filter: materialise fully (if backed, materialise now)
            if backed:
                adata = adata.to_memory()

        # 鈹€鈹€ Prepare training data 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        train_data, train_tps, all_unique_tps, _, data_np, cell_tps_np = prepare_data(
            adata, time_key, train_times
        )
        n_genes = train_data[0].shape[1]

        # Resolve n_sim_cells: config override, else cap-bounded actual t0 count
        t0_count = train_data[0].shape[0]
        if n_sim_cells_cfg is not None:
            n_sim_cells = int(n_sim_cells_cfg)
        else:
            n_sim_cells = min(t0_count, int(_DEFAULTS["n_sim_cells_cap"]))
        print(f"[scNODE] n_genes={n_genes}, n_sim_cells={n_sim_cells} "
              f"(t0_count={t0_count})")

        # 鈹€鈹€ Train (or load cache) 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        model = train_or_load(train_data, train_tps, n_genes, scnode_cfg, model_cache)

        # 鈹€鈹€ Load full adata for Forecast + Embedding (need held-out cells) 鈹€鈹€鈹€鈹€
        # self.adata is train-filtered; held-out cells live only in the full file.
        if heldout_times:
            print(f"[scNODE] Re-loading full adata for held-out cell access ...")
            adata_full = anndata.read_h5ad(str(h5ad_path))
        else:
            adata_full = adata   # Scenario A: no holdout 鈫?same object is fine

        # 鈹€鈹€ Forecast Accuracy 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        print(f"\n[scNODE] 鈹€鈹€ Forecast Accuracy 鈹€鈹€")
        run_forecast_accuracy(
            model, adata_full, time_key,
            all_unique_tps, heldout_times, n_sim_cells, output_dir,
            metric_sample_cells=metric_sample_cells,
            seed=seed,
        )

        # 鈹€鈹€ Embedding Coherence 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        print(f"\n[scNODE] 鈹€鈹€ Embedding Coherence 鈹€鈹€")
        run_embedding_coherence(
            model, adata_full, time_key,
            all_unique_tps, heldout_times, n_sim_cells, cell_state_key, output_dir,
        )

        # 鈹€鈹€ Lineage Fidelity 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        print(f"\n[scNODE] 鈹€鈹€ Lineage Fidelity 鈹€鈹€")
        lf_stm, lf_edges = run_lineage_fidelity(
            model, adata, time_key, all_unique_tps, cell_state_key, output_dir,
        )

        # 鈹€鈹€ Evaluate lineage metrics if reference graph available 鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€鈹€
        if ref_graph_path:
            ref_abs = (Path(ref_graph_path) if Path(ref_graph_path).is_absolute()
                       else project_root / ref_graph_path)
            from benchmark.evaluation.eval_lineage import run_lineage_evaluation
            lf_metrics = run_lineage_evaluation(
                state_transition_matrix_path=lf_stm,
                lineage_graph_edges_path=lf_edges,
                output_dir=str(output_dir),
                reference_graph_path=str(ref_abs),
                edge_confidence_mode=edge_conf_mode,
                exclude_uncertain_states=excl_uncertain,
                cell_state_key=cell_state_key,
                adata=adata,
            )
            print(f"[scNODE] Lineage metrics: {lf_metrics}")
        else:
            print("[scNODE] No reference_graph_path in config 鈥?"
                  "lineage metrics not computed.")

        status = "completed"

    except Exception:
        err_notes = traceback.format_exc()
        print(f"[scNODE] ERROR:\n{err_notes}")
        raise

    finally:
        elapsed = time.time() - t_wall
        metadata = {
            "method":          "scnode",
            "run_id":          run_id,
            "dataset":         dataset_cfg.get("id", "unknown"),
            "scenario":        scenario_id,
            "capability_flags": {
                "supports_unseen_timepoint_projection": True,
                "supports_lineage_inference":           True,
            },
            "dimensions_executed": [
                "forecast_accuracy",
                "embedding_coherence",
                "lineage_fidelity",
            ],
            "runtime_seconds": round(elapsed, 2),
            "status":          status,
            "train_times":     train_times if train_times else "all",
            "heldout_times":   heldout_times,
            "n_sim_cells":     n_sim_cells if "n_sim_cells" in dir() else None,
            "notes":           err_notes,
        }
        meta_path = output_dir / "run_metadata.json"
        with open(str(meta_path), "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)
        print(f"\n[scNODE] Finished in {elapsed:.1f}s  鈫? status={status}")
        print(f"[scNODE] Results: {output_dir}")


if __name__ == "__main__":
    main()
