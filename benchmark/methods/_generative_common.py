"""benchmark/methods/_generative_common.py

Model-AGNOSTIC runner shared by the projection-capable generative methods
(scIMF, PI-SDE, Squidiff). The Forecast Accuracy, Embedding Coherence, and
Lineage Fidelity logic here is ported verbatim (in behaviour) from the proven
benchmark/methods/scNODE/run.py, with scNODE's model.predict / model.vaeReconstruct
replaced by a small, uniform ``Engine`` protocol so any model can plug in.

A method's run.py only needs to:
  1. implement an Engine (train + project + embed), and
  2. expose the 5 contract functions by delegating to the helpers below.

Engine protocol (all arrays numpy, gene-expression space unless noted):
  engine.embed(X)                      -> latent  (n_cells, d)
  engine.project(X0, target_tps, n)    -> (recon, latent_seq)
        recon      : (n, len(target_tps), n_genes)   gene-expression predictions
        latent_seq : (n, len(target_tps), d)         latent at each target tp
  where X0 is the (n0, n_genes) cells at the first observed timepoint, target_tps
  is a sorted list of float timepoints, and n is the number of cells to simulate.

NOTE on forecast metric: the official Forecast Accuracy number is recomputed
centrally (metric_backend=scTimeBench_exact) from the standardized
``projected_expression.npy`` this module writes. The wasserstein value written into
forecast_metrics.json here is a self-contained sanity value (POT if available, else a
sliced-Wasserstein fallback) and is overwritten by the central evaluator.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from scipy.sparse import issparse


# --------------------------------------------------------------------------- #
# Data preparation (identical contract to scNODE/run.py prepare_data)
# --------------------------------------------------------------------------- #
def prepare_data(adata, time_key: str, train_times=None):
    data_np = adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)
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
        train_data.append(np.asarray(data_np[mask, :], dtype=np.float32))
        print(f"  [data] tp={t:.2f}: {int(mask.sum())} cells")
    return train_data, train_times_f, all_unique_tps, heldout_tps, data_np, cell_tps_np


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _dense(adata):
    return adata.X.toarray() if issparse(adata.X) else np.asarray(adata.X)


def _sample_rows(X: np.ndarray, max_cells, seed: int) -> np.ndarray:
    if max_cells is None or int(max_cells) <= 0 or X.shape[0] <= int(max_cells):
        return X
    rng = np.random.default_rng(seed)
    return X[rng.choice(X.shape[0], size=int(max_cells), replace=False)]


def _wasserstein(X_obs: np.ndarray, X_pred: np.ndarray) -> float:
    """Self-contained OT-Wasserstein (POT if available, else sliced fallback).
    The central evaluator recomputes the official exact metric from
    projected_expression.npy; this is only a sanity value."""
    try:
        import ot  # POT
        a = np.ones(X_obs.shape[0]) / X_obs.shape[0]
        b = np.ones(X_pred.shape[0]) / X_pred.shape[0]
        M = ot.dist(X_obs, X_pred, metric="euclidean")
        return float(ot.emd2(a, b, M))
    except Exception:
        # sliced-Wasserstein-1 fallback (deterministic projections)
        rng = np.random.default_rng(0)
        d = X_obs.shape[1]
        n_proj = 50
        dirs = rng.standard_normal((n_proj, d))
        dirs /= np.linalg.norm(dirs, axis=1, keepdims=True) + 1e-12
        po = X_obs @ dirs.T
        pp = X_pred @ dirs.T
        po.sort(axis=0); pp.sort(axis=0)
        m = min(po.shape[0], pp.shape[0])
        idx_o = np.linspace(0, po.shape[0] - 1, m).astype(int)
        idx_p = np.linspace(0, pp.shape[0] - 1, m).astype(int)
        return float(np.mean(np.abs(po[idx_o] - pp[idx_p])))


# --------------------------------------------------------------------------- #
# Forecast Accuracy
# --------------------------------------------------------------------------- #
def run_forecast_accuracy(model, adata_full, time_key, all_unique_tps, heldout_tps,
                          n_sim_cells, output_dir, metric_sample_cells=1000, seed=42):
    output_dir = Path(output_dir)
    data_full = _dense(adata_full)
    cell_tps = adata_full.obs[time_key].values.astype(float)
    t0 = all_unique_tps[0]
    X0 = data_full[np.isclose(cell_tps, t0), :].astype(np.float32)
    eval_tps = list(heldout_tps) if heldout_tps else list(all_unique_tps)
    all_eval_tps = sorted(set(list(all_unique_tps) + list(heldout_tps)))

    recon, _ = model.project(X0, all_eval_tps, n_sim_cells)   # (n, T, g)
    recon = np.asarray(recon)

    rows, preds = [], []
    for t_h in eval_tps:
        if t_h not in all_eval_tps:
            continue
        X_pred = recon[:, all_eval_tps.index(t_h), :]
        X_obs = data_full[np.isclose(cell_tps, t_h), :]
        if X_obs.shape[0] == 0:
            continue
        Xo = _sample_rows(X_obs, metric_sample_cells, seed + 10_000 + len(rows))
        Xp = _sample_rows(X_pred, metric_sample_cells, seed + 20_000 + len(rows))
        wd = _wasserstein(Xo, Xp)
        rows.append({"timepoint": t_h, "wasserstein_ot": wd,
                     "n_obs_cells": int(X_obs.shape[0]), "n_pred_cells": int(X_pred.shape[0])})
        preds.append(X_pred)
        print(f"  [Forecast] t={t_h:.2f}: Wasserstein={wd:.4f}")

    np.save(str(output_dir / "projected_expression.npy"),
            np.stack(preds, axis=0) if preds else np.array([]))
    pd.DataFrame(rows).to_csv(str(output_dir / "per_timepoint_forecast_metrics.csv"), index=False)
    mean_wd = float(np.mean([r["wasserstein_ot"] for r in rows])) if rows else None
    summary = {"wasserstein_distance": mean_wd, "gaussian_mmd": None,
               "energy_distance_mmd": None, "hausdorff_loss": None,
               "n_eval_timepoints": len(rows), "eval_timepoints": eval_tps,
               "scenario_type": "forecast" if heldout_tps else "reconstruction",
               "metric_sample_cells": int(metric_sample_cells),
               "metric_backend": "self_contained_sanity (official = central recompute)",
               "status": "completed"}
    with open(output_dir / "forecast_metrics.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    print(f"[Forecast] Done. Mean Wasserstein={mean_wd}")
    return (str(output_dir / "projected_expression.npy"),
            str(output_dir / "forecast_metrics.json"),
            str(output_dir / "per_timepoint_forecast_metrics.csv"))


# --------------------------------------------------------------------------- #
# Embedding Coherence
# --------------------------------------------------------------------------- #
def run_embedding_coherence(model, adata_full, time_key, all_unique_tps, heldout_tps,
                            n_sim_cells, cell_state_key, output_dir):
    from scipy.spatial.distance import cdist
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import adjusted_rand_score

    output_dir = Path(output_dir)
    data_full = _dense(adata_full)
    cell_tps = adata_full.obs[time_key].values.astype(float)
    t0 = all_unique_tps[0]
    X0 = data_full[np.isclose(cell_tps, t0), :].astype(np.float32)
    eval_tps = list(heldout_tps) if heldout_tps else list(all_unique_tps)

    obs_latent = np.asarray(model.embed(data_full))        # (n_cells, d)
    np.save(str(output_dir / "embedding.npy"), obs_latent)

    state_labels = adata_full.obs[cell_state_key].values if cell_state_key in adata_full.obs.columns else None
    centroids, states = {}, []
    if state_labels is not None:
        states = sorted(set(state_labels))
        for s in states:
            m = state_labels == s
            if m.sum() > 0:
                centroids[s] = obs_latent[m].mean(axis=0)
    clf = None
    if state_labels is not None:
        try:
            clf = LogisticRegression(max_iter=500, class_weight="balanced", random_state=0)
            clf.fit(obs_latent, state_labels)
        except Exception as exc:
            print(f"[Embedding] classifier entropy skipped: {exc}")

    all_eval_tps = sorted(set(list(all_unique_tps) + list(heldout_tps)))
    _, latent_seq = model.project(X0, all_eval_tps, n_sim_cells)
    latent_seq = np.asarray(latent_seq)

    proj_latents, label_rows, ari_rows, entropy = [], [], [], []
    for t_h in eval_tps:
        if t_h not in all_eval_tps:
            continue
        pl = latent_seq[:, all_eval_tps.index(t_h), :]
        proj_latents.append(pl)
        if centroids:
            cmat = np.stack([centroids[s] for s in states])
            labels = [states[i] for i in np.argmin(cdist(pl, cmat), axis=1)]
        else:
            labels = ["unknown"] * pl.shape[0]
        if state_labels is not None and np.isclose(cell_tps, float(t_h)).any():
            om = np.isclose(cell_tps, float(t_h))
            ref = state_labels[om][np.argmin(cdist(pl, obs_latent[om]), axis=1)]
            ari_rows.append({"timepoint": float(t_h),
                             "adjusted_rand_index": float(adjusted_rand_score(ref, labels)),
                             "n_projected_cells": int(pl.shape[0])})
        if clf is not None:
            probs = clf.predict_proba(pl)
            if probs.shape[1] > 1:
                e = -np.sum(probs * np.log(probs + 1e-12), axis=1) / np.log(probs.shape[1])
                entropy.extend(e.tolist())
        for i, lbl in enumerate(labels):
            label_rows.append({"projected_timepoint": t_h, "cell_idx": i,
                               "projected_cluster_label": lbl})

    emb = np.vstack(proj_latents) if proj_latents else np.array([])
    np.save(str(output_dir / "projected_embedding.npy"), emb)
    np.save(str(output_dir / "next_timepoint_embedding.npy"), emb)
    pd.DataFrame(label_rows).to_csv(str(output_dir / "projected_cluster_labels.csv"), index=False)
    pd.DataFrame(ari_rows).to_csv(str(output_dir / "per_timepoint_embedding_metrics.csv"), index=False)
    mean_ari = float(np.mean([r["adjusted_rand_index"] for r in ari_rows])) if ari_rows else None
    mean_ent = float(np.mean(entropy)) if entropy else None
    metrics = {"adjusted_rand_index": mean_ari, "avg_normalized_classifier_entropy": mean_ent,
               "projected_latent_dim": int(latent_seq.shape[2]), "n_eval_timepoints": len(eval_tps),
               "eval_timepoints": eval_tps,
               "status": "completed" if mean_ari is not None else "embedding_metrics_incomplete",
               "ari_reference": "nearest observed cell at the same tp in the model latent space",
               "entropy_reference": "balanced logistic regression on observed latent + provider labels"}
    with open(output_dir / "embedding_metrics.json", "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)
    print(f"[Embedding] projected_embedding.npy shape: {emb.shape}")
    return (str(output_dir / "projected_embedding.npy"),
            str(output_dir / "embedding_metrics.json"),
            str(output_dir / "projected_cluster_labels.csv"))


# --------------------------------------------------------------------------- #
# Lineage Fidelity
# --------------------------------------------------------------------------- #
def run_lineage_fidelity(model, adata_train, time_key, train_unique_tps,
                         cell_state_key, output_dir):
    from scipy.spatial.distance import cdist
    output_dir = Path(output_dir)
    data = _dense(adata_train)
    cell_tps = adata_train.obs[time_key].values.astype(float)
    if cell_state_key not in adata_train.obs.columns:
        raise KeyError(f"[Lineage] cell_state_key '{cell_state_key}' not in obs "
                       f"({list(adata_train.obs.columns)})")
    state_labels = adata_train.obs[cell_state_key].values
    states = sorted(set(state_labels))
    n = len(states); s2i = {s: i for i, s in enumerate(states)}

    obs_latent = np.asarray(model.embed(data))
    centroids = np.zeros((n, obs_latent.shape[1]))
    for s, i in s2i.items():
        m = state_labels == s
        if m.sum() > 0:
            centroids[i] = obs_latent[m].mean(axis=0)

    tps = sorted(train_unique_tps)
    stm = np.zeros((n, n)); counts = np.zeros(n)
    for k in range(len(tps) - 1):
        tc, tn = float(tps[k]), float(tps[k + 1])
        cm = np.isclose(cell_tps, tc)
        if cm.sum() == 0:
            continue
        cur = data[cm, :].astype(np.float32)
        _, lat = model.project(cur, [tc, tn], cur.shape[0])
        pred_nxt = np.asarray(lat)[:, 1, :]
        d = cdist(pred_nxt, centroids, metric="euclidean")
        scores = -d
        finite = np.isfinite(scores)
        rmax = np.where(finite.any(1), np.max(np.where(finite, scores, -np.inf), axis=1), 0.0).reshape(-1, 1)
        w = np.where(finite, np.exp(np.where(finite, scores - rmax, -np.inf)), 0.0)
        rs = w.sum(1, keepdims=True)
        good = (rs[:, 0] > 0) & np.isfinite(rs[:, 0])
        w = np.divide(w, rs, out=np.zeros_like(w), where=good[:, None])
        for bad in np.where(~good)[0]:
            fd = np.isfinite(d[bad])
            if fd.any():
                w[bad, int(np.argmin(np.where(fd, d[bad], np.inf)))] = 1.0
        cur_states = state_labels[cm]
        for ci, src in enumerate(cur_states):
            if w[ci].sum() > 0:
                stm[s2i[src]] += w[ci]; counts[s2i[src]] += 1

    norm = np.zeros_like(stm)
    for i in range(n):
        if counts[i] > 0:
            norm[i] = stm[i] / counts[i]
    norm = np.nan_to_num(norm)
    stm_df = pd.DataFrame(norm, index=pd.Index(states, name="source_state"), columns=states)
    stm_df.to_csv(str(output_dir / "state_transition_matrix.csv"))
    edges = [{"source_state": s, "target_state": t, "weight": float(norm[s2i[s], s2i[t]])}
             for s in states for t in states if norm[s2i[s], s2i[t]] > 1e-8]
    pd.DataFrame(edges).to_csv(str(output_dir / "lineage_graph_edges.csv"), index=False)
    print(f"[Lineage] STM ({n}x{n}) + {len(edges)} edges written.")
    return (str(output_dir / "state_transition_matrix.csv"),
            str(output_dir / "lineage_graph_edges.csv"))
