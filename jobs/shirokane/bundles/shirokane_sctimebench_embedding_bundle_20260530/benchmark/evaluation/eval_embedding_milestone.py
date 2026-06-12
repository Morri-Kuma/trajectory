"""
eval_embedding_milestone.py
============================
Step 8 / Step 9: Milestone-aware Embedding Coherence evaluator.

Two operating modes
-------------------
observed-smoke
    Uses milestone labels present in the benchmark h5ad obs columns.
    Cluster labels from an existing obs column or KMeans. (Step 8)

run-output
    Uses method-side embedding.npy as the observed embedding space and
    next_timepoint_embedding.npy/projected_embedding.npy as projected cells.
    Projected labels are assigned by kNN transfer from observed embeddings, and
    entropy is computed from classifier probability vectors.
    (Step 9 path — projected cells via annotate_projected_cells.py)

CLI
---
    python -m benchmark.evaluation.eval_embedding_milestone \\
        --input-h5ad <path> \\
        --output-dir <dir> \\
        --dataset-id GSE230659 \\
        --mode observed-smoke \\
        --exclude-label ambiguous

    python -m benchmark.evaluation.eval_embedding_milestone \\
        --input-h5ad <path> \\
        --output-dir <dir> \\
        --dataset-id GSE230659 \\
        --mode run-output \\
        --run-output-dir <method_run_dir> \\
        --exclude-label ambiguous
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Label-mode catalogue
# ---------------------------------------------------------------------------

LABEL_MODES: Dict[str, str] = {
    "official_silver": "final_milestone_label_coarse",
    "expanded_sensitivity": "final_milestone_label_expanded",
}

LABEL_MODE_PROVIDER_SUFFIX: Dict[str, str] = {
    "official_silver": "marker_fm_transition_silver_v1",
    "expanded_sensitivity": "marker_fm_transition_silver_v1",
}

DEFAULT_LABEL_MODES: List[str] = ["official_silver"]

CLUSTER_CANDIDATES: List[str] = [
    "leiden",
    "louvain",
]


# ---------------------------------------------------------------------------
# h5py-based h5ad loader (avoids null-encoding AnnData issue)
# ---------------------------------------------------------------------------

def _load_h5ad_h5py(h5ad_path: Path):
    """Return (obs_dict, obsm_dict, n_obs) without importing anndata."""
    try:
        import h5py
    except ImportError as exc:
        raise ImportError("h5py is required for h5ad loading: pip install h5py") from exc

    obs_dict: Dict[str, np.ndarray] = {}
    obsm_dict: Dict[str, np.ndarray] = {}

    with h5py.File(h5ad_path, "r") as f:
        obs_grp = f["obs"]
        for key in obs_grp.keys():
            obj = obs_grp[key]
            import h5py as _h
            if isinstance(obj, _h.Dataset):
                raw = obj[:]
                obs_dict[key] = np.array([
                    v.decode() if isinstance(v, bytes) else v for v in raw
                ])
            elif isinstance(obj, _h.Group):
                if "categories" in obj and "codes" in obj:
                    cats = [c.decode() if isinstance(c, bytes) else str(c)
                            for c in obj["categories"][:]]
                    codes = obj["codes"][:]
                    obs_dict[key] = np.array([
                        cats[int(c)] if (c >= 0 and int(c) < len(cats)) else "NA"
                        for c in codes
                    ])
        if "obsm" in f:
            for key in f["obsm"].keys():
                try:
                    obsm_dict[key] = f["obsm"][key][:]
                except Exception:
                    pass

    n_obs = len(next(iter(obs_dict.values()))) if obs_dict else (
        next(iter(obsm_dict.values())).shape[0] if obsm_dict else 0
    )
    return obs_dict, obsm_dict, n_obs


# ---------------------------------------------------------------------------
# Projected milestone labels loader (Step 9)
# ---------------------------------------------------------------------------

def _load_projected_milestone_obs(csv_path: Path) -> Dict[str, np.ndarray]:
    """
    Read projected_milestone_labels.csv and return obs-style dict
    keyed on the milestone label columns. Only milestone columns are extracted.
    """
    with open(csv_path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError(f"projected_milestone_labels.csv is empty: {csv_path}")

    result: Dict[str, np.ndarray] = {}
    for col in LABEL_MODES.values():
        if col in rows[0]:
            result[col] = np.array([r[col] for r in rows])
    for col in rows[0]:
        if col.startswith("prob_"):
            result[col] = np.array([float(r[col]) if r[col] != "" else np.nan for r in rows])
    return result


# ---------------------------------------------------------------------------
# Embedding / clustering helpers
# ---------------------------------------------------------------------------

def _detect_embedding(obsm_dict: Dict[str, np.ndarray],
                       preferred_key: Optional[str] = None) -> Tuple[str, np.ndarray]:
    priority = ["X_pca", "X_umap"]
    if preferred_key:
        if preferred_key not in obsm_dict:
            raise ValueError(
                f"Embedding key {preferred_key!r} not found in obsm. "
                f"Available: {sorted(obsm_dict)}"
            )
        return preferred_key, obsm_dict[preferred_key]
    for k in priority:
        if k in obsm_dict:
            return k, obsm_dict[k]
    if obsm_dict:
        k = next(iter(obsm_dict))
        return k, obsm_dict[k]
    raise ValueError("No embeddings found in obsm.")


def _get_cluster_labels(
    obs_dict: Dict[str, np.ndarray],
    embedding: np.ndarray,
    preferred_key: Optional[str] = None,
    n_obs: int = 0,
) -> Tuple[np.ndarray, str]:
    if preferred_key:
        if preferred_key not in obs_dict:
            raise ValueError(
                f"Cluster key {preferred_key!r} not found in obs. "
                f"Available: {sorted(obs_dict)}"
            )
        return obs_dict[preferred_key], f"obs_column:{preferred_key}"
    for k in CLUSTER_CANDIDATES:
        if k in obs_dict:
            return obs_dict[k], f"obs_column:{k}"
    try:
        from sklearn.cluster import KMeans
    except ImportError as exc:
        raise RuntimeError(
            "No cluster column found and sklearn unavailable for KMeans. "
            "Install scikit-learn or supply --cluster-key."
        ) from exc
    n_clusters = min(10, max(2, int(math.sqrt(n_obs))))
    print(f"[eval_embedding_milestone] KMeans(n_clusters={n_clusters}) from embedding.")
    km = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = km.fit_predict(embedding).astype(str)
    return labels, f"kmeans_n{n_clusters}_from_embedding"


def _leiden_from_embedding_array(
    emb: np.ndarray,
    embedding_name: str,
    source_prefix: str,
    n_neighbors: int = 15,
    resolution: float = 0.5,
) -> Tuple[np.ndarray, str]:
    """Run scTimeBench-style embedding -> kNN graph -> Leiden clustering."""
    if emb.ndim > 2:
        emb = emb.reshape(-1, emb.shape[-1])
    if emb.ndim != 2:
        raise ValueError(f"Expected 2D projected embedding, got shape={emb.shape}")
    if emb.shape[0] < 2:
        raise ValueError(f"Need at least two cells for Leiden, got {emb.shape[0]}")

    try:
        import anndata as ad
        import scanpy as sc
    except ImportError as exc:
        raise ImportError(
            "scanpy and anndata are required for Leiden clustering in run-output mode."
        ) from exc

    k = max(1, min(int(n_neighbors), emb.shape[0] - 1))
    a = ad.AnnData(X=np.zeros((emb.shape[0], 1), dtype=np.float32))
    a.obsm["X_projected"] = emb.astype(np.float32, copy=False)
    sc.pp.neighbors(a, n_neighbors=k, use_rep="X_projected")
    try:
        sc.tl.leiden(
            a,
            resolution=float(resolution),
            key_added="leiden_projected",
            random_state=42,
            flavor="igraph",
            directed=False,
            n_iterations=2,
        )
    except TypeError:
        sc.tl.leiden(
            a,
            resolution=float(resolution),
            key_added="leiden_projected",
            random_state=42,
        )

    labels = a.obs["leiden_projected"].astype(str).to_numpy()
    src = (
        f"{source_prefix}:{embedding_name}:"
        f"n_neighbors={k}:resolution={float(resolution)}"
    )
    return labels, src


def _leiden_from_projected_embedding(
    embedding_path: Path,
    n_neighbors: int = 15,
    resolution: float = 0.5,
) -> Tuple[np.ndarray, str]:
    """Run scTimeBench-style clustering on projected_embedding.npy."""
    if not embedding_path.exists():
        raise FileNotFoundError(f"projected_embedding.npy not found: {embedding_path}")

    emb = np.load(embedding_path, allow_pickle=False)
    return _leiden_from_embedding_array(
        emb,
        embedding_name=embedding_path.name,
        source_prefix="leiden_from_projected_embedding",
        n_neighbors=n_neighbors,
        resolution=resolution,
    )


def _load_embedding_array(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Embedding file not found: {path}")
    emb = np.load(path, allow_pickle=False)
    if emb.ndim > 2:
        emb = emb.reshape(-1, emb.shape[-1])
    if emb.ndim != 2:
        raise ValueError(f"Expected 2D embedding array from {path}, got shape={emb.shape}")
    return emb


# ---------------------------------------------------------------------------
# Core metric functions
# ---------------------------------------------------------------------------

def _adjusted_rand_index(a: np.ndarray, b: np.ndarray) -> float:
    try:
        from sklearn.metrics import adjusted_rand_score
        return float(adjusted_rand_score(a, b))
    except ImportError as exc:
        raise RuntimeError("scikit-learn required for ARI.") from exc


def _normalized_label_entropy_per_cluster(
    cluster_labels: np.ndarray,
    milestone_labels: np.ndarray,
) -> Dict[str, Any]:
    clusters = np.unique(cluster_labels)
    label_set = np.unique(milestone_labels)
    n_labels = len(label_set)
    if n_labels <= 1:
        return {
            "per_cluster_entropy": {},
            "mean_normalized_entropy": None,
            "weighted_mean_normalized_entropy": None,
            "n_clusters": int(len(clusters)),
            "n_labels": int(n_labels),
            "note": "entropy undefined: only one label present",
        }
    log_n = math.log(n_labels)
    per_cluster: Dict[str, float] = {}
    cluster_sizes: Dict[str, int] = {}
    for cl in clusters:
        mask = cluster_labels == cl
        sub = milestone_labels[mask]
        n = len(sub)
        cluster_sizes[str(cl)] = n
        counts = Counter(sub.tolist())
        entropy = 0.0
        for cnt in counts.values():
            p = cnt / n
            if p > 0:
                entropy -= p * math.log(p)
        per_cluster[str(cl)] = round(float(entropy / log_n), 6)
    total = sum(cluster_sizes.values())
    mean_ent = float(np.mean(list(per_cluster.values())))
    weighted_ent = sum(
        per_cluster[k] * cluster_sizes[k] / total for k in per_cluster
    )
    return {
        "per_cluster_entropy": per_cluster,
        "mean_normalized_entropy": round(mean_ent, 6),
        "weighted_mean_normalized_entropy": round(weighted_ent, 6),
        "n_clusters": int(len(clusters)),
        "n_labels": int(n_labels),
    }


def _prediction_entropy_from_probabilities(
    probability_matrix: Optional[np.ndarray],
    cluster_labels: np.ndarray,
) -> Dict[str, Any]:
    """scTimeBench-style per-cell normalized prediction entropy."""
    if probability_matrix is None:
        return {
            "mean_prediction_entropy": None,
            "weighted_prediction_entropy": None,
            "per_cluster_prediction_entropy": {},
            "n_probability_classes": 0,
            "prediction_entropy_note": "probability vector unavailable",
            "entropy_basis": "unavailable",
        }

    probs = np.asarray(probability_matrix, dtype=np.float64)
    if probs.ndim != 2 or probs.shape[0] != len(cluster_labels):
        return {
            "mean_prediction_entropy": None,
            "weighted_prediction_entropy": None,
            "per_cluster_prediction_entropy": {},
            "n_probability_classes": int(probs.shape[1]) if probs.ndim == 2 else 0,
            "prediction_entropy_note": (
                f"probability matrix shape {probs.shape} does not match "
                f"{len(cluster_labels)} cluster labels"
            ),
            "entropy_basis": "milestone_probability_vector",
        }

    k = int(probs.shape[1])
    if k <= 1:
        return {
            "mean_prediction_entropy": None,
            "weighted_prediction_entropy": None,
            "per_cluster_prediction_entropy": {},
            "n_probability_classes": k,
            "prediction_entropy_note": "prediction entropy undefined: probability vector has <=1 class",
            "entropy_basis": "milestone_probability_vector",
        }

    probs = np.clip(probs, 0.0, None)
    row_sums = probs.sum(axis=1, keepdims=True)
    zero_rows = row_sums[:, 0] <= 0
    probs = probs / np.where(row_sums <= 0, 1.0, row_sums)
    if zero_rows.any():
        probs[zero_rows, :] = 1.0 / k

    log_k = math.log(k)
    cell_entropy = -np.sum(
        np.where(probs > 0, probs * np.log(np.where(probs > 0, probs, 1.0)), 0.0),
        axis=1,
    ) / log_k

    clusters = np.unique(cluster_labels)
    per_cluster: Dict[str, float] = {}
    cluster_sizes: Dict[str, int] = {}
    for cl in clusters:
        mask = cluster_labels == cl
        cluster_sizes[str(cl)] = int(np.sum(mask))
        per_cluster[str(cl)] = round(float(np.mean(cell_entropy[mask])), 6)

    mean_ent = float(np.mean(cell_entropy)) if len(cell_entropy) else float("nan")
    total = sum(cluster_sizes.values())
    weighted_ent = (
        sum(per_cluster[k_] * cluster_sizes[k_] / total for k_ in per_cluster)
        if total else float("nan")
    )
    return {
        "mean_prediction_entropy": round(mean_ent, 6),
        "weighted_prediction_entropy": round(float(weighted_ent), 6),
        "per_cluster_prediction_entropy": per_cluster,
        "n_probability_classes": k,
        "prediction_entropy_note": "",
        "entropy_basis": "milestone_probability_vector",
    }


def _entropy_vector_from_probabilities(
    probabilities: np.ndarray,
    eps: float = 1e-10,
) -> Tuple[np.ndarray, np.ndarray]:
    probs = np.asarray(probabilities, dtype=np.float64)
    raw_entropy = -np.sum(probs * np.log(probs + eps), axis=1)
    normalized_entropy = raw_entropy / math.log(probs.shape[1])
    return raw_entropy, normalized_entropy


def _knn_transfer_labels(
    reference_embedding: np.ndarray,
    reference_labels: np.ndarray,
    projected_embedding: np.ndarray,
    n_neighbors: int,
) -> np.ndarray:
    """Transfer labels to projected embeddings exactly like scTimeBench ARI."""
    try:
        from sklearn.neighbors import NearestNeighbors
    except ImportError as exc:
        raise RuntimeError("scikit-learn required for kNN label transfer.") from exc

    k = min(int(n_neighbors), reference_embedding.shape[0])
    if k <= 0:
        raise ValueError("Need at least one reference cell for kNN label transfer.")
    knn_model = NearestNeighbors(n_neighbors=k)
    knn_model.fit(reference_embedding)
    _, neighbor_indices = knn_model.kneighbors(projected_embedding)

    transferred: List[Any] = []
    for neighbors in neighbor_indices:
        neighbor_labels = reference_labels[neighbors]
        labels, counts = np.unique(neighbor_labels, return_counts=True)
        transferred.append(labels[counts.argmax()])
    return np.asarray(transferred).astype(str)


def _sctimebench_classifier_entropy(
    reference_embedding: np.ndarray,
    reference_labels: np.ndarray,
    projected_embedding: np.ndarray,
    cluster_labels: np.ndarray,
    test_size: float = 0.2,
    random_state: int = 42,
) -> Dict[str, Any]:
    """scTimeBench ClassificationEntropy over classifier probability vectors."""
    if reference_embedding.shape[0] != len(reference_labels):
        return {
            "avg_entropy": None,
            "std_entropy": None,
            "avg_normalized_entropy": None,
            "std_normalized_entropy": None,
            "pred_tp_avg_entropy": None,
            "pred_tp_avg_normalized_entropy": None,
            "pred_tp_std_entropy": None,
            "mean_prediction_entropy": None,
            "weighted_prediction_entropy": None,
            "per_cluster_prediction_entropy": {},
            "n_probability_classes": 0,
            "prediction_entropy_note": (
                "reference embedding and label lengths do not match"
            ),
            "entropy_basis": "classifier_probability_vector",
        }
    if projected_embedding.shape[0] != len(cluster_labels):
        return {
            "avg_entropy": None,
            "std_entropy": None,
            "avg_normalized_entropy": None,
            "std_normalized_entropy": None,
            "pred_tp_avg_entropy": None,
            "pred_tp_avg_normalized_entropy": None,
            "pred_tp_std_entropy": None,
            "mean_prediction_entropy": None,
            "weighted_prediction_entropy": None,
            "per_cluster_prediction_entropy": {},
            "n_probability_classes": 0,
            "prediction_entropy_note": (
                "projected embedding and Leiden cluster lengths do not match"
            ),
            "entropy_basis": "classifier_probability_vector",
        }

    labels_used = np.unique(reference_labels)
    if len(labels_used) <= 1:
        return {
            "avg_entropy": None,
            "std_entropy": None,
            "avg_normalized_entropy": None,
            "std_normalized_entropy": None,
            "pred_tp_avg_entropy": None,
            "pred_tp_avg_normalized_entropy": None,
            "pred_tp_std_entropy": None,
            "mean_prediction_entropy": None,
            "weighted_prediction_entropy": None,
            "per_cluster_prediction_entropy": {},
            "n_probability_classes": int(len(labels_used)),
            "prediction_entropy_note": (
                "classification entropy undefined: <=1 reference label"
            ),
            "entropy_basis": "classifier_probability_vector",
        }

    try:
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.model_selection import train_test_split
    except ImportError as exc:
        raise RuntimeError("scikit-learn required for classifier entropy.") from exc

    try:
        x_train, x_test, y_train, _ = train_test_split(
            reference_embedding,
            reference_labels,
            test_size=test_size,
            random_state=random_state,
        )
        classifier = RandomForestClassifier(
            n_estimators=100,
            max_depth=None,
            random_state=random_state,
        )
        classifier.fit(x_train, y_train)
        test_probas = classifier.predict_proba(x_test)
        projected_probas = classifier.predict_proba(projected_embedding)
    except Exception as exc:
        return {
            "avg_entropy": None,
            "std_entropy": None,
            "avg_normalized_entropy": None,
            "std_normalized_entropy": None,
            "pred_tp_avg_entropy": None,
            "pred_tp_avg_normalized_entropy": None,
            "pred_tp_std_entropy": None,
            "mean_prediction_entropy": None,
            "weighted_prediction_entropy": None,
            "per_cluster_prediction_entropy": {},
            "n_probability_classes": int(len(labels_used)),
            "prediction_entropy_note": f"classifier entropy error: {exc}",
            "entropy_basis": "classifier_probability_vector",
        }

    num_classes = int(test_probas.shape[1])
    if num_classes <= 1:
        return {
            "avg_entropy": None,
            "std_entropy": None,
            "avg_normalized_entropy": None,
            "std_normalized_entropy": None,
            "pred_tp_avg_entropy": None,
            "pred_tp_avg_normalized_entropy": None,
            "pred_tp_std_entropy": None,
            "mean_prediction_entropy": None,
            "weighted_prediction_entropy": None,
            "per_cluster_prediction_entropy": {},
            "n_probability_classes": num_classes,
            "prediction_entropy_note": (
                "classification entropy undefined: classifier has <=1 class"
            ),
            "entropy_basis": "classifier_probability_vector",
        }

    test_entropy, test_norm_entropy = _entropy_vector_from_probabilities(test_probas)
    pred_entropy, pred_norm_entropy = _entropy_vector_from_probabilities(projected_probas)

    per_cluster: Dict[str, float] = {}
    cluster_sizes: Dict[str, int] = {}
    for cl in np.unique(cluster_labels):
        mask = cluster_labels == cl
        cluster_sizes[str(cl)] = int(np.sum(mask))
        per_cluster[str(cl)] = round(float(np.mean(pred_norm_entropy[mask])), 6)
    total = sum(cluster_sizes.values())
    weighted_pred_norm = (
        sum(per_cluster[k] * cluster_sizes[k] / total for k in per_cluster)
        if total else float("nan")
    )

    return {
        "avg_entropy": float(np.mean(test_entropy)),
        "std_entropy": float(np.std(test_entropy)),
        "avg_normalized_entropy": float(np.mean(test_norm_entropy)),
        "std_normalized_entropy": float(np.std(test_norm_entropy)),
        "num_classes": num_classes,
        "classifier_classes": [str(c) for c in classifier.classes_.tolist()],
        "pred_tp_avg_entropy": float(np.mean(pred_entropy)),
        "pred_tp_avg_normalized_entropy": float(np.mean(pred_norm_entropy)),
        "pred_tp_std_entropy": float(np.std(pred_entropy)),
        "mean_prediction_entropy": round(float(np.mean(pred_norm_entropy)), 6),
        "weighted_prediction_entropy": round(float(weighted_pred_norm), 6),
        "per_cluster_prediction_entropy": per_cluster,
        "n_probability_classes": num_classes,
        "prediction_entropy_note": "",
        "entropy_basis": "classifier_probability_vector",
        "classifier_model": "RandomForestClassifier",
        "classifier_n_estimators": 100,
        "classifier_random_state": random_state,
        "classifier_test_size": test_size,
    }


def _compute_reference_ari(
    obs_dict: Dict[str, np.ndarray],
    state_key: str,
    embedding: np.ndarray,
    embedding_key: str,
    exclude_labels: List[str],
    n_neighbors: int,
    resolution: float,
) -> Dict[str, Any]:
    """Compute observed/reference ARI using observed embedding and labels."""
    result: Dict[str, Any] = {
        "reference_embedding_key": embedding_key,
        "reference_cluster_source": None,
        "reference_adjusted_rand_index": None,
        "reference_ari_status": "not_computed",
        "reference_n_cells_total": 0,
        "reference_n_cells_evaluated": 0,
        "reference_n_clusters": None,
        "reference_n_labels": None,
        "reference_labels_used": [],
    }

    if state_key not in obs_dict:
        result["reference_ari_status"] = f"error: obs column {state_key!r} not found"
        return result

    labels_full = obs_dict[state_key]
    result["reference_n_cells_total"] = int(len(labels_full))
    if len(labels_full) != embedding.shape[0]:
        result["reference_ari_status"] = (
            f"error: length mismatch: labels={len(labels_full)}, "
            f"embedding_rows={embedding.shape[0]}"
        )
        return result

    if exclude_labels:
        keep = ~np.isin(labels_full, exclude_labels)
    else:
        keep = np.ones(len(labels_full), dtype=bool)

    labels = labels_full[keep]
    emb = embedding[keep]
    labels_used = sorted(set(labels.tolist()))
    result["reference_n_cells_evaluated"] = int(np.sum(keep))
    result["reference_n_labels"] = int(len(labels_used))
    result["reference_labels_used"] = labels_used

    if emb.shape[0] < 2:
        result["reference_ari_status"] = "skipped: fewer than two cells after exclusions"
        return result
    if len(labels_used) < 2:
        result["reference_ari_status"] = "skipped: fewer than two labels after exclusions"
        return result

    try:
        ref_clusters, ref_src = _leiden_from_embedding_array(
            emb,
            embedding_name=embedding_key,
            source_prefix="leiden_from_observed_embedding",
            n_neighbors=n_neighbors,
            resolution=resolution,
        )
        ref_ari = _adjusted_rand_index(ref_clusters, labels)
        result.update({
            "reference_cluster_source": ref_src,
            "reference_adjusted_rand_index": round(float(ref_ari), 6),
            "reference_ari_status": "ok",
            "reference_n_clusters": int(len(np.unique(ref_clusters))),
        })
    except Exception as exc:
        result["reference_ari_status"] = f"error: {exc}"
    return result


# ---------------------------------------------------------------------------
# Per-mode evaluator
# ---------------------------------------------------------------------------

def _evaluate_one_mode(
    label_mode: str,
    state_key: str,
    obs_dict: Dict[str, np.ndarray],
    cluster_labels_full: np.ndarray,
    cluster_source: str,
    embedding_key: str,
    reference_embedding: np.ndarray,
    exclude_labels: List[str],
    dataset_id: str,
    input_h5ad: str,
    output_dir: Path,
    provider_id_prefix: str,
    leiden_n_neighbors: int,
    leiden_resolution: float,
    milestone_labels_override: Optional[np.ndarray] = None,
    probability_matrix_override: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Evaluate one label mode; write embedding_metrics_<label_mode>.json.

    milestone_labels_override
        If provided (Step 9 run-output), use these projected milestone labels
        instead of obs_dict[state_key]. Must have same length as cluster_labels_full.
    """
    provider_suffix = LABEL_MODE_PROVIDER_SUFFIX[label_mode]
    provider_id = f"{provider_id_prefix}_{provider_suffix}"
    reference_ari_result = _compute_reference_ari(
        obs_dict=obs_dict,
        state_key=state_key,
        embedding=reference_embedding,
        embedding_key=embedding_key,
        exclude_labels=exclude_labels,
        n_neighbors=leiden_n_neighbors,
        resolution=leiden_resolution,
    )

    # Resolve milestone labels
    if milestone_labels_override is not None:
        milestone_labels_full = milestone_labels_override
        label_source = "projected_milestone_labels_csv"
    elif state_key in obs_dict:
        milestone_labels_full = obs_dict[state_key]
        label_source = "reference_h5ad_obs"
    else:
        result = {
            "dataset_id": dataset_id,
            "label_mode": label_mode,
            "state_key": state_key,
            "provider_id": provider_id,
            "input_h5ad": input_h5ad,
            "embedding_key": embedding_key,
            "cluster_source": cluster_source,
            "status": f"error: obs column {state_key!r} not found",
            "adjusted_rand_index": None,
            "mean_normalized_entropy": None,
            "weighted_mean_normalized_entropy": None,
            **reference_ari_result,
        }
        _write_mode_json(result, label_mode, output_dir)
        return result

    # Length check
    if len(milestone_labels_full) != len(cluster_labels_full):
        result = {
            "dataset_id": dataset_id,
            "label_mode": label_mode,
            "state_key": state_key,
            "provider_id": provider_id,
            "input_h5ad": input_h5ad,
            "embedding_key": embedding_key,
            "cluster_source": cluster_source,
            "status": (
                f"error: length mismatch — milestone_labels has "
                f"{len(milestone_labels_full)} cells, cluster_labels has "
                f"{len(cluster_labels_full)} cells."
            ),
            "adjusted_rand_index": None,
            "mean_normalized_entropy": None,
            "weighted_mean_normalized_entropy": None,
            **reference_ari_result,
        }
        _write_mode_json(result, label_mode, output_dir)
        return result

    n_cells_total = len(milestone_labels_full)

    # Exclusion mask
    if exclude_labels:
        keep = ~np.isin(milestone_labels_full, exclude_labels)
    else:
        keep = np.ones(n_cells_total, dtype=bool)

    milestone_labels = milestone_labels_full[keep]
    cluster_labels = cluster_labels_full[keep]
    n_cells_evaluated = int(np.sum(keep))

    label_counts = {k: int(v) for k, v in Counter(milestone_labels_full.tolist()).items()}
    cluster_counts = {k: int(v) for k, v in Counter(cluster_labels_full.tolist()).items()}
    ambiguous_count = int(np.sum(milestone_labels_full == "ambiguous"))
    missing_label_count = int(np.sum(milestone_labels_full == "NA"))
    labels_used = sorted(set(milestone_labels.tolist()))
    n_labels = len(labels_used)

    # ARI
    try:
        ari = _adjusted_rand_index(cluster_labels, milestone_labels)
        ari_status = "ok"
    except Exception as exc:
        ari = None
        ari_status = f"error: {exc}"

    reference_ari = reference_ari_result.get("reference_adjusted_rand_index")
    if ari is not None and isinstance(reference_ari, (int, float)) and reference_ari > 0:
        ari_retention_fraction = round(float(ari) / float(reference_ari), 6)
        ari_retention_status = "ok"
    elif ari is None:
        ari_retention_fraction = None
        ari_retention_status = "skipped: projected ARI unavailable"
    elif reference_ari is None:
        ari_retention_fraction = None
        ari_retention_status = "skipped: reference ARI unavailable"
    else:
        ari_retention_fraction = None
        ari_retention_status = "skipped: reference ARI is not positive"

    # Entropy
    if n_cells_evaluated > 0 and n_labels > 0:
        ent_result = _normalized_label_entropy_per_cluster(cluster_labels, milestone_labels)
    else:
        ent_result = {
            "per_cluster_entropy": {},
            "mean_normalized_entropy": None,
            "weighted_mean_normalized_entropy": None,
            "n_clusters": 0,
            "n_labels": 0,
        }
    probability_matrix = (
        probability_matrix_override[keep, :]
        if probability_matrix_override is not None and len(probability_matrix_override) == n_cells_total
        else None
    )
    pred_ent_result = _prediction_entropy_from_probabilities(
        probability_matrix=probability_matrix,
        cluster_labels=cluster_labels,
    )

    result: Dict[str, Any] = {
        "dataset_id": dataset_id,
        "label_mode": label_mode,
        "state_key": state_key,
        "provider_id": provider_id,
        "input_h5ad": input_h5ad,
        "embedding_key": embedding_key,
        "cluster_source": cluster_source,
        "label_source": label_source,
        "n_cells_total": n_cells_total,
        "n_cells_evaluated": n_cells_evaluated,
        "excluded_labels": list(exclude_labels),
        "exclusion_reason": (
            "label not present in marker-defined milestone graph"
            if exclude_labels else "none"
        ),
        "n_labels": n_labels,
        "labels_used": labels_used,
        "label_counts": label_counts,
        "cluster_counts": cluster_counts,
        "n_clusters": ent_result.get("n_clusters", 0),
        "ambiguous_count": ambiguous_count,
        "ambiguous_fraction": round(ambiguous_count / n_cells_total, 4) if n_cells_total > 0 else 0.0,
        "missing_label_count": missing_label_count,
        "adjusted_rand_index": round(ari, 6) if ari is not None else None,
        "ari_status": ari_status,
        **reference_ari_result,
        "ari_retention_fraction": ari_retention_fraction,
        "ari_retention_status": ari_retention_status,
        "mean_prediction_entropy": pred_ent_result.get("mean_prediction_entropy"),
        "weighted_prediction_entropy": pred_ent_result.get("weighted_prediction_entropy"),
        "per_cluster_prediction_entropy": pred_ent_result.get("per_cluster_prediction_entropy", {}),
        "prediction_entropy_note": pred_ent_result.get("prediction_entropy_note", ""),
        "entropy_basis": pred_ent_result.get("entropy_basis", ""),
        "n_probability_classes": pred_ent_result.get("n_probability_classes", 0),
        "mean_normalized_entropy": ent_result.get("mean_normalized_entropy"),
        "weighted_mean_normalized_entropy": ent_result.get("weighted_mean_normalized_entropy"),
        "per_cluster_entropy": ent_result.get("per_cluster_entropy", {}),
        "entropy_note": ent_result.get("note", ""),
        "note": ent_result.get("note", ""),
        "status": "completed",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    _write_mode_json(result, label_mode, output_dir)
    return result


def _evaluate_one_mode_sctimebench(
    label_mode: str,
    state_key: str,
    obs_dict: Dict[str, np.ndarray],
    cluster_labels_full: np.ndarray,
    cluster_source: str,
    embedding_key: str,
    reference_embedding: np.ndarray,
    projected_embedding: np.ndarray,
    projected_embedding_source: str,
    n_projected_total: int,
    exclude_labels: List[str],
    dataset_id: str,
    input_h5ad: str,
    output_dir: Path,
    provider_id_prefix: str,
    leiden_n_neighbors: int,
    leiden_resolution: float,
    milestone_labels_override: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """Evaluate one label mode with the scTimeBench embedding protocol."""
    provider_suffix = LABEL_MODE_PROVIDER_SUFFIX[label_mode]
    provider_id = f"{provider_id_prefix}_{provider_suffix}"
    reference_ari_result = _compute_reference_ari(
        obs_dict=obs_dict,
        state_key=state_key,
        embedding=reference_embedding,
        embedding_key=embedding_key,
        exclude_labels=exclude_labels,
        n_neighbors=leiden_n_neighbors,
        resolution=leiden_resolution,
    )

    def _error_result(status: str) -> Dict[str, Any]:
        result = {
            "dataset_id": dataset_id,
            "label_mode": label_mode,
            "state_key": state_key,
            "provider_id": provider_id,
            "input_h5ad": input_h5ad,
            "embedding_key": embedding_key,
            "embedding_metric_protocol": "sctimebench",
            "projected_embedding_source": projected_embedding_source,
            "cluster_source": cluster_source,
            "status": status,
            "adjusted_rand_index": None,
            "ari_next_timepoint": None,
            "pred_tp_avg_normalized_entropy": None,
            **reference_ari_result,
        }
        _write_mode_json(result, label_mode, output_dir)
        return result

    if state_key not in obs_dict:
        return _error_result(f"error: obs column {state_key!r} not found")

    reference_labels_full = obs_dict[state_key].astype(str)
    if len(reference_labels_full) != reference_embedding.shape[0]:
        return _error_result(
            "error: length mismatch: observed labels have "
            f"{len(reference_labels_full)} cells, observed embedding has "
            f"{reference_embedding.shape[0]} cells."
        )
    if len(cluster_labels_full) != projected_embedding.shape[0]:
        return _error_result(
            "error: length mismatch: Leiden clusters have "
            f"{len(cluster_labels_full)} cells, projected embedding has "
            f"{projected_embedding.shape[0]} cells."
        )

    if exclude_labels:
        reference_keep = ~np.isin(reference_labels_full, exclude_labels)
    else:
        reference_keep = np.ones(len(reference_labels_full), dtype=bool)
    reference_embedding_eval = reference_embedding[reference_keep]
    reference_labels = reference_labels_full[reference_keep]
    n_reference_cells_evaluated = int(np.sum(reference_keep))
    if n_reference_cells_evaluated == 0:
        return _error_result("error: no reference cells remain after label exclusions")

    try:
        transferred_labels = _knn_transfer_labels(
            reference_embedding=reference_embedding_eval,
            reference_labels=reference_labels,
            projected_embedding=projected_embedding,
            n_neighbors=leiden_n_neighbors,
        )
    except Exception as exc:
        return _error_result(f"error: kNN projected-label transfer failed: {exc}")

    cluster_labels = cluster_labels_full
    labels_used = sorted(set(transferred_labels.tolist()))
    n_labels = len(labels_used)
    try:
        ari = _adjusted_rand_index(cluster_labels, transferred_labels)
        ari_status = "ok"
    except Exception as exc:
        ari = None
        ari_status = f"error: {exc}"

    reference_ari = reference_ari_result.get("reference_adjusted_rand_index")
    if ari is not None and isinstance(reference_ari, (int, float)) and reference_ari > 0:
        ari_retention_fraction = round(float(ari) / float(reference_ari), 6)
        ari_retention_status = "ok"
    elif ari is None:
        ari_retention_fraction = None
        ari_retention_status = "skipped: projected ARI unavailable"
    elif reference_ari is None:
        ari_retention_fraction = None
        ari_retention_status = "skipped: reference ARI unavailable"
    else:
        ari_retention_fraction = None
        ari_retention_status = "skipped: reference ARI is not positive"

    pred_ent_result = _sctimebench_classifier_entropy(
        reference_embedding=reference_embedding_eval,
        reference_labels=reference_labels,
        projected_embedding=projected_embedding,
        cluster_labels=cluster_labels,
    )

    pred_tp_avg_normalized_entropy = pred_ent_result.get(
        "pred_tp_avg_normalized_entropy"
    )
    status = (
        "completed"
        if ari is not None and pred_tp_avg_normalized_entropy is not None
        else "error: scTimeBench ARI/classifier entropy incomplete"
    )
    n_cells_total = int(n_projected_total)
    n_cells_evaluated = int(projected_embedding.shape[0])
    label_counts = {k: int(v) for k, v in Counter(transferred_labels.tolist()).items()}
    reference_label_counts = {
        k: int(v) for k, v in Counter(reference_labels_full.tolist()).items()
    }
    cluster_counts = {k: int(v) for k, v in Counter(cluster_labels.tolist()).items()}
    ambiguous_count = int(np.sum(transferred_labels == "ambiguous"))
    missing_label_count = int(np.sum(transferred_labels == "NA"))

    result: Dict[str, Any] = {
        "dataset_id": dataset_id,
        "label_mode": label_mode,
        "state_key": state_key,
        "provider_id": provider_id,
        "input_h5ad": input_h5ad,
        "embedding_key": embedding_key,
        "embedding_metric_protocol": "sctimebench",
        "reference_embedding_source": embedding_key,
        "projected_embedding_source": projected_embedding_source,
        "cluster_source": cluster_source,
        "label_source": "knn_transfer_from_observed_embedding",
        "n_cells_total": n_cells_total,
        "n_cells_evaluated": n_cells_evaluated,
        "n_projected_cells_total": n_cells_total,
        "n_projected_cells_valid": n_cells_evaluated,
        "n_reference_cells_total": int(len(reference_labels_full)),
        "n_reference_cells_evaluated": n_reference_cells_evaluated,
        "excluded_labels": list(exclude_labels),
        "exclusion_reason": (
            "label not present in marker-defined milestone graph"
            if exclude_labels else "none"
        ),
        "n_labels": n_labels,
        "labels_used": labels_used,
        "label_counts": label_counts,
        "reference_label_counts": reference_label_counts,
        "cluster_counts": cluster_counts,
        "n_clusters": int(len(set(cluster_labels.tolist()))),
        "ambiguous_count": ambiguous_count,
        "ambiguous_fraction": (
            round(ambiguous_count / n_cells_evaluated, 4)
            if n_cells_evaluated > 0 else 0.0
        ),
        "missing_label_count": missing_label_count,
        "adjusted_rand_index": round(ari, 6) if ari is not None else None,
        "ari_next_timepoint": round(ari, 6) if ari is not None else None,
        "ari_status": ari_status,
        **reference_ari_result,
        "ari_ground_truth": reference_ari_result.get("reference_adjusted_rand_index"),
        "ari_retention_fraction": ari_retention_fraction,
        "ari_retention_status": ari_retention_status,
        "avg_entropy": pred_ent_result.get("avg_entropy"),
        "std_entropy": pred_ent_result.get("std_entropy"),
        "avg_normalized_entropy": pred_ent_result.get("avg_normalized_entropy"),
        "std_normalized_entropy": pred_ent_result.get("std_normalized_entropy"),
        "num_classes": pred_ent_result.get("num_classes"),
        "classifier_classes": pred_ent_result.get("classifier_classes", []),
        "pred_tp_avg_entropy": pred_ent_result.get("pred_tp_avg_entropy"),
        "pred_tp_avg_normalized_entropy": pred_tp_avg_normalized_entropy,
        "pred_tp_std_entropy": pred_ent_result.get("pred_tp_std_entropy"),
        "entropy_basis": pred_ent_result.get("entropy_basis", ""),
        "classifier_model": pred_ent_result.get("classifier_model"),
        "classifier_n_estimators": pred_ent_result.get("classifier_n_estimators"),
        "classifier_random_state": pred_ent_result.get("classifier_random_state"),
        "classifier_test_size": pred_ent_result.get("classifier_test_size"),
        "status": status,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    _write_mode_json(result, label_mode, output_dir)
    return result


def _write_mode_json(result: Dict[str, Any], label_mode: str, output_dir: Path):
    # Step 10: ensure compatibility_note is present so validators and summary
    # scripts can identify these as mode-specific milestone embedding outputs.
    if "compatibility_note" not in result:
        result["compatibility_note"] = (
            f"Mode-specific embedding coherence metrics for label_mode={label_mode!r}. "
            "See provider_agreement_metrics.json for cross-provider comparison. "
            "Step 10 output naming policy: these files are the preferred outputs "
            "for official silver reporting; embedding_metrics.json is not written."
        )
    path = output_dir / f"embedding_metrics_{label_mode}.json"
    _atomic_write_json(path, result)
    print(
        f"[eval_embedding_milestone] Wrote {path.name}  "
        f"(ARI={result.get('adjusted_rand_index')}, "
        f"reference_ARI={result.get('reference_adjusted_rand_index')}, "
        f"n_cells={result.get('n_cells_evaluated')})"
    )


# ---------------------------------------------------------------------------
# Provider agreement
# ---------------------------------------------------------------------------

def _compute_provider_agreement(
    obs_dict: Dict[str, np.ndarray],
    exclude_labels: List[str],
    dataset_id: str,
    input_h5ad: str,
    output_dir: Path,
    provider_id_prefix: str,
    label_modes: Optional[List[str]] = None,
    projected_milestone_obs: Optional[Dict[str, np.ndarray]] = None,
) -> Dict[str, Any]:
    """Compare milestone label columns pairwise; write provider_agreement_metrics.json.

    projected_milestone_obs: if provided (Step 9), use projected labels instead of obs.
    """
    modes = list(label_modes or LABEL_MODES.keys())
    state_keys = [LABEL_MODES[m] for m in modes]
    provider_ids = [f"{provider_id_prefix}_{LABEL_MODE_PROVIDER_SUFFIX[m]}" for m in modes]

    if len(modes) < 2:
        result = {
            "dataset_id": dataset_id,
            "input_h5ad": input_h5ad,
            "provider_ids": provider_ids,
            "label_modes": modes,
            "state_keys": state_keys,
            "status": "skipped: fewer than two label modes requested",
        }
        _write_agreement_json(result, output_dir)
        return result

    label_source = projected_milestone_obs if projected_milestone_obs is not None else obs_dict
    missing = [sk for sk in state_keys if sk not in label_source]
    if missing:
        result = {
            "dataset_id": dataset_id,
            "input_h5ad": input_h5ad,
            "provider_ids": provider_ids,
            "label_modes": modes,
            "state_keys": state_keys,
            "status": f"error: missing label columns: {missing}",
        }
        _write_agreement_json(result, output_dir)
        return result

    arrays = {sk: label_source[sk] for sk in state_keys}
    n_cells_total = len(next(iter(arrays.values())))
    all_labels_mat = np.stack(list(arrays.values()), axis=1)
    if exclude_labels:
        keep = ~np.any(np.isin(all_labels_mat, exclude_labels), axis=1)
    else:
        keep = np.ones(n_cells_total, dtype=bool)
    n_cells_compared = int(np.sum(keep))

    pairwise_ari: Dict[str, Any] = {}
    pairwise_exact: Dict[str, float] = {}
    pairwise_crosstabs: Dict[str, Any] = {}
    pairs = [
        (modes[i], modes[j])
        for i in range(len(modes))
        for j in range(i + 1, len(modes))
    ]
    for m1, m2 in pairs:
        sk1, sk2 = LABEL_MODES[m1], LABEL_MODES[m2]
        a1 = arrays[sk1][keep]
        a2 = arrays[sk2][keep]
        pair_key = f"{m1}_vs_{m2}"
        try:
            ari_val = round(_adjusted_rand_index(a1, a2), 6)
        except Exception as exc:
            ari_val = None
        exact_frac = round(float(np.mean(a1 == a2)), 6)
        xtab: Dict[str, Dict[str, int]] = {}
        for v1 in sorted(set(a1.tolist())):
            xtab[v1] = {}
            for v2 in sorted(set(a2.tolist())):
                xtab[v1][v2] = int(np.sum((a1 == v1) & (a2 == v2)))
        pairwise_ari[pair_key] = ari_val
        pairwise_exact[pair_key] = exact_frac
        pairwise_crosstabs[pair_key] = xtab

    result = {
        "dataset_id": dataset_id,
        "input_h5ad": input_h5ad,
        "provider_ids": provider_ids,
        "label_modes": modes,
        "state_keys": state_keys,
        "n_cells_total": n_cells_total,
        "n_cells_compared": n_cells_compared,
        "excluded_labels": list(exclude_labels),
        "exclusion_reason": (
            "label not present in marker-defined milestone graph"
            if exclude_labels else "none"
        ),
        "pairwise_adjusted_rand_index": pairwise_ari,
        "pairwise_exact_match_fraction": pairwise_exact,
        "pairwise_crosstabs": pairwise_crosstabs,
        "status": "completed",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    _write_agreement_json(result, output_dir)
    return result


def _write_agreement_json(result: Dict[str, Any], output_dir: Path):
    path = output_dir / "provider_agreement_metrics.json"
    _atomic_write_json(path, result)
    print(f"[eval_embedding_milestone] Wrote {path.name}")


def _atomic_write_json(path: Path, payload: Dict[str, Any]) -> None:
    tmp_path = path.with_name(f"{path.name}.tmp")
    with open(tmp_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
        f.write("\n")
    tmp_path.replace(path)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_milestone_embedding_evaluation(
    input_h5ad: str,
    output_dir: str,
    dataset_id: str,
    mode: str = "observed-smoke",
    embedding_key: Optional[str] = None,
    cluster_key: Optional[str] = None,
    exclude_labels: Optional[List[str]] = None,
    run_output_dir: Optional[str] = None,
    label_modes: Optional[List[str]] = None,
    leiden_n_neighbors: int = 15,
    leiden_resolution: float = 0.5,
) -> Dict[str, Any]:
    """
    Main evaluation function for Step 8/9 milestone embedding coherence.

    Parameters
    ----------
    input_h5ad : Path to milestone-annotated benchmark h5ad.
    output_dir : Directory for output JSONs.
    dataset_id : e.g. 'GSE230659'.
    mode : 'observed-smoke' or 'run-output'.
    embedding_key : obsm key (auto-detected if None).
    cluster_key : obs column (auto-detected if None).
    exclude_labels : Labels to exclude from ARI/entropy.
    run_output_dir : For 'run-output' mode, directory containing method-side
        embedding.npy plus next_timepoint_embedding.npy or projected_embedding.npy.
    """
    h5ad_path = Path(input_h5ad)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    exclude_labels = list(exclude_labels) if exclude_labels else []
    provider_id_prefix = dataset_id.lower()
    selected_label_modes = list(label_modes or DEFAULT_LABEL_MODES)
    unknown_modes = [m for m in selected_label_modes if m not in LABEL_MODES]
    if unknown_modes:
        raise ValueError(
            f"Unknown label mode(s): {unknown_modes}. "
            f"Available: {sorted(LABEL_MODES)}"
        )

    print(f"[eval_embedding_milestone] Loading h5ad: {h5ad_path}")
    obs_dict, obsm_dict, n_obs = _load_h5ad_h5py(h5ad_path)
    print(f"[eval_embedding_milestone] n_obs={n_obs}, obsm keys={sorted(obsm_dict)}")

    # --- Embedding (always from reference h5ad for now) ---
    emb_key, embedding = _detect_embedding(obsm_dict, preferred_key=embedding_key)
    print(f"[eval_embedding_milestone] embedding_key={emb_key!r}, shape={embedding.shape}")

    # --- Cluster labels + projected milestone labels ---
    projected_milestone_obs: Optional[Dict[str, np.ndarray]] = None
    projected_embedding_for_eval: Optional[np.ndarray] = None
    projected_embedding_source = emb_key
    n_projected_total = int(embedding.shape[0])
    use_sctimebench_protocol = False

    if mode == "run-output" and run_output_dir:
        rdir = Path(run_output_dir)

        observed_embedding_path = rdir / "embedding.npy"
        if not observed_embedding_path.exists():
            raise FileNotFoundError(
                "scTimeBench-aligned embedding coherence requires method-side "
                f"observed embeddings at {observed_embedding_path}. Rerun the "
                "method adapter after the embedding.npy output update."
            )
        embedding = _load_embedding_array(observed_embedding_path)
        emb_key = observed_embedding_path.name
        print(
            f"[eval_embedding_milestone] observed method embedding={emb_key!r}, "
            f"shape={embedding.shape}"
        )

        next_embedding_path = rdir / "next_timepoint_embedding.npy"
        projected_embedding_path = (
            next_embedding_path
            if next_embedding_path.exists()
            else rdir / "projected_embedding.npy"
        )
        projected_embedding_full = _load_embedding_array(projected_embedding_path)
        valid_projected = ~np.isnan(projected_embedding_full).any(axis=1)
        projected_embedding_for_eval = projected_embedding_full[valid_projected]
        n_projected_total = int(projected_embedding_full.shape[0])
        if projected_embedding_for_eval.shape[0] == 0:
            raise ValueError(
                f"No valid projected embeddings in {projected_embedding_path}"
            )
        projected_embedding_source = projected_embedding_path.name

        # scTimeBench ARI path: next-timepoint embedding -> kNN graph -> Leiden.
        source_prefix = (
            "leiden_from_next_timepoint_embedding"
            if projected_embedding_path.name == "next_timepoint_embedding.npy"
            else "leiden_from_projected_embedding"
        )
        cluster_labels_full, cluster_src = _leiden_from_embedding_array(
            projected_embedding_for_eval,
            embedding_name=projected_embedding_path.name,
            source_prefix=source_prefix,
            n_neighbors=leiden_n_neighbors,
            resolution=leiden_resolution,
        )
        n_projected = len(cluster_labels_full)
        use_sctimebench_protocol = True
        print(
            f"[eval_embedding_milestone] run-output mode: "
            f"{n_projected}/{n_projected_total} valid projected cells"
        )

        # Load projected milestone labels if available; these are diagnostics only
        # for the scTimeBench-aligned main ARI/entropy computation.
        proj_milestone_path = rdir / "projected_milestone_labels.csv"
        if proj_milestone_path.exists():
            projected_milestone_obs = _load_projected_milestone_obs(proj_milestone_path)
            n_pm = len(next(iter(projected_milestone_obs.values())))
            print(f"[eval_embedding_milestone] projected_milestone_labels.csv loaded: "
                  f"{n_pm} cells (diagnostic)")
            if n_pm not in {n_projected, n_projected_total}:
                raise ValueError(
                    f"projected_milestone_labels.csv has {n_pm} rows but "
                    f"projected embedding has {n_projected_total} rows and "
                    f"Leiden clusters have {n_projected} rows."
                )
        else:
            print(f"[eval_embedding_milestone] No projected_milestone_labels.csv found "
                  f"in {rdir}; continuing with scTimeBench kNN label transfer.")

    else:
        cluster_labels_full, cluster_src = _get_cluster_labels(
            obs_dict, embedding,
            preferred_key=cluster_key,
            n_obs=n_obs,
        )
        projected_embedding_for_eval = embedding

    print(f"[eval_embedding_milestone] cluster_source={cluster_src!r}, "
          f"n_clusters={len(np.unique(cluster_labels_full))}")

    # --- Per-mode evaluation ---
    all_results: Dict[str, Any] = {}
    probability_override: Optional[np.ndarray] = None
    if projected_milestone_obs is not None:
        prob_cols = sorted(k for k in projected_milestone_obs if k.startswith("prob_"))
        if prob_cols:
            probability_override = np.column_stack([projected_milestone_obs[k] for k in prob_cols])
            print(
                f"[eval_embedding_milestone] probability columns loaded: "
                f"{prob_cols} shape={probability_override.shape}"
            )
    for label_mode in selected_label_modes:
        state_key = LABEL_MODES[label_mode]
        print(f"\n[eval_embedding_milestone] --- mode={label_mode}, key={state_key} ---")
        override = (
            projected_milestone_obs.get(state_key)
            if projected_milestone_obs is not None else None
        )
        if use_sctimebench_protocol:
            assert projected_embedding_for_eval is not None
            res = _evaluate_one_mode_sctimebench(
                label_mode=label_mode,
                state_key=state_key,
                obs_dict=obs_dict,
                cluster_labels_full=cluster_labels_full,
                cluster_source=cluster_src,
                embedding_key=emb_key,
                reference_embedding=embedding,
                projected_embedding=projected_embedding_for_eval,
                projected_embedding_source=projected_embedding_source,
                n_projected_total=n_projected_total,
                exclude_labels=exclude_labels,
                dataset_id=dataset_id,
                input_h5ad=str(h5ad_path),
                output_dir=out_dir,
                provider_id_prefix=provider_id_prefix,
                leiden_n_neighbors=leiden_n_neighbors,
                leiden_resolution=leiden_resolution,
                milestone_labels_override=override,
            )
        else:
            res = _evaluate_one_mode(
                label_mode=label_mode,
                state_key=state_key,
                obs_dict=obs_dict,
                cluster_labels_full=cluster_labels_full,
                cluster_source=cluster_src,
                embedding_key=emb_key,
                reference_embedding=embedding,
                exclude_labels=exclude_labels,
                dataset_id=dataset_id,
                input_h5ad=str(h5ad_path),
                output_dir=out_dir,
                provider_id_prefix=provider_id_prefix,
                leiden_n_neighbors=leiden_n_neighbors,
                leiden_resolution=leiden_resolution,
                milestone_labels_override=override,
                probability_matrix_override=probability_override,
            )
        all_results[label_mode] = res

    # --- Provider agreement ---
    print("\n[eval_embedding_milestone] --- provider agreement ---")
    agreement = _compute_provider_agreement(
        obs_dict=obs_dict,
        exclude_labels=exclude_labels,
        dataset_id=dataset_id,
        input_h5ad=str(h5ad_path),
        output_dir=out_dir,
        provider_id_prefix=provider_id_prefix,
        label_modes=selected_label_modes,
        projected_milestone_obs=projected_milestone_obs,
    )
    all_results["provider_agreement"] = agreement

    print(f"\n[eval_embedding_milestone] Done. Outputs in {out_dir}")
    return all_results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Step 8/9: Milestone-aware Embedding Coherence evaluation."
    )
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--mode", default="observed-smoke",
                        choices=["observed-smoke", "run-output"])
    parser.add_argument("--embedding-key", default=None)
    parser.add_argument("--cluster-key", default=None)
    parser.add_argument("--exclude-label", action="append", dest="exclude_labels",
                        default=[], metavar="LABEL")
    parser.add_argument("--run-output-dir", default=None)
    parser.add_argument("--leiden-n-neighbors", type=int, default=15)
    parser.add_argument("--leiden-resolution", type=float, default=0.5)
    parser.add_argument(
        "--label-mode",
        action="append",
        choices=sorted(LABEL_MODES),
        default=None,
        help=(
            "Label mode to evaluate. Repeatable. Defaults to official_silver; "
            "for the current official benchmark use --label-mode official_silver."
        ),
    )
    args = parser.parse_args(argv)

    try:
        run_milestone_embedding_evaluation(
            input_h5ad=args.input_h5ad,
            output_dir=args.output_dir,
            dataset_id=args.dataset_id,
            mode=args.mode,
            embedding_key=args.embedding_key,
            cluster_key=args.cluster_key,
            exclude_labels=args.exclude_labels,
            run_output_dir=args.run_output_dir,
            label_modes=args.label_mode,
            leiden_n_neighbors=args.leiden_n_neighbors,
            leiden_resolution=args.leiden_resolution,
        )
    except Exception as exc:
        print(f"[eval_embedding_milestone] ERROR: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
