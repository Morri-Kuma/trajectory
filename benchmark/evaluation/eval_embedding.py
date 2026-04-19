"""
eval_embedding.py
Embedding Coherence evaluator — framework module (currently inactive).
Framework reference: experimental framework v2.md §8, §14 Step 4

STATUS: INACTIVE at the current benchmark stage.

WOT and CellRank2 are NOT eligible for Embedding Coherence because they do not
generate projected cells at unseen future time points. This evaluator will be
activated only when generative / forecasting models are added to the benchmark.

This module is implemented as a clean framework stub so that the evaluation
logic is ready when eligible models arrive. It must NOT be called for WOT or
CellRank2. The dispatcher (eval_dispatch.py) enforces this via capability flags.

Metrics (per v2 §8.3):
  - Adjusted Rand Index (ARI)
  - Average normalized classifier entropy

Aggregation (per v2 §8.4):
  Rank methods on ARI, rank methods on entropy, average the ranks
  for the Embedding Coherence rank within each scenario.
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional


# ------------------------------------------------------------------
# Metric functions
# ------------------------------------------------------------------

def adjusted_rand_index(cluster_labels_pred: np.ndarray,
                         cluster_labels_obs: np.ndarray) -> float:
    """
    Compute Adjusted Rand Index between predicted and observed cluster labels.
    """
    from sklearn.metrics import adjusted_rand_score
    return float(adjusted_rand_score(cluster_labels_obs, cluster_labels_pred))


def avg_normalized_classifier_entropy(classifier_probs: np.ndarray) -> float:
    """
    Compute average normalized classifier entropy across projected cells.

    classifier_probs : np.ndarray, shape (n_cells, n_clusters)
        Soft class assignment probabilities from a classifier trained on
        observed cells and applied to projected cells.
    """
    # Normalized entropy per cell: H / log(n_clusters)
    n_clusters = classifier_probs.shape[1]
    if n_clusters <= 1:
        return float("nan")
    eps = 1e-12
    entropy_per_cell = -np.sum(
        classifier_probs * np.log(classifier_probs + eps), axis=1
    )
    normalized_entropy = entropy_per_cell / np.log(n_clusters)
    return float(np.mean(normalized_entropy))


# ------------------------------------------------------------------
# Main evaluation function
# ------------------------------------------------------------------

def run_embedding_evaluation(
    projected_embedding_path: str,
    projected_cluster_labels_path: str,
    adata,
    output_dir: str,
) -> dict:
    """
    Evaluate Embedding Coherence for one method × scenario run.

    IMPORTANT: This function must only be called for methods where
    supports_unseen_timepoint_projection = True.
    The dispatcher (eval_dispatch.py) enforces this. Do not call this
    function directly for WOT or CellRank2.

    Parameters
    ----------
    projected_embedding_path : str
        Path to projected_embedding.npy produced by the method adapter.
    projected_cluster_labels_path : str
        Path to projected_cluster_labels.csv produced by the method adapter.
    adata : AnnData
        Benchmark input with observed cells.
    output_dir : str
        Directory where embedding_metrics.json will be written.

    Returns
    -------
    dict with keys: adjusted_rand_index, avg_normalized_classifier_entropy,
                    embedding_coherence_rank, status
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    emb_path = Path(projected_embedding_path)
    labels_path = Path(projected_cluster_labels_path)

    if not emb_path.exists():
        metrics = {
            "adjusted_rand_index": None,
            "avg_normalized_classifier_entropy": None,
            "status": "error: projected_embedding.npy not found",
        }
        _write_metrics(metrics, out_dir)
        return metrics

    X_emb = np.load(emb_path, allow_pickle=True)

    if X_emb.size == 0:
        metrics = {
            "adjusted_rand_index": None,
            "avg_normalized_classifier_entropy": None,
            "status": "inactive: no projected embedding data. "
                       "This evaluator is inactive for the current benchmark stage.",
        }
        _write_metrics(metrics, out_dir)
        return metrics

    # TODO: load cluster labels and compute ARI.
    # TODO: train classifier on observed cells; apply to projected cells; compute entropy.

    metrics = {
        "adjusted_rand_index": None,
        "avg_normalized_classifier_entropy": None,
        "status": "not_implemented: metric computation not yet wired.",
    }
    _write_metrics(metrics, out_dir)
    return metrics


def _write_metrics(metrics: dict, out_dir: Path):
    metrics_path = out_dir / "embedding_metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"[eval_embedding] Metrics written to {metrics_path}")
