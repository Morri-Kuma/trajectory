"""Threshold-free / robustness companions for Lineage Fidelity evaluation.

Reviewer request (R1, P2.3): on the small reference graphs (3-4 edges for the
marker-silver systems, 9 for OSKM) point estimates of AUROC/Jaccard are
threshold- and negative-set-sensitive. The main evaluator (eval_lineage.py)
already reports the threshold-free AUROC and AUPRC; this module adds the missing
companions:

  * edge_confusion        -- the explicit edge-level TP/FP/FN/TN over all ordered
                             off-diagonal node pairs (the negative set is made
                             explicit, not implicit).
  * pr_curve_points       -- precision/recall arrays for a threshold-free PR curve.
  * bootstrap_edge_metric -- bootstrap CI over the candidate edge set for any
                             score->metric, so a 3-4 edge graph reports an
                             interval rather than a bare point estimate.

Dependency-light (numpy + scikit-learn) and unit-tested, so a reviewer can run it
without GPU. It scores the SAME predicted transition matrices the main evaluator
consumes; regenerating the per-run CIs across the benchmark is a downstream sweep.
"""
from __future__ import annotations

from typing import Callable, Dict, Iterable, List, Tuple
import numpy as np


def _ordered_offdiag_pairs(nodes: List[str]) -> List[Tuple[str, str]]:
    return [(a, b) for a in nodes for b in nodes if a != b]


def edge_confusion(pred_edges: Iterable[Tuple[str, str]],
                   ref_edges: Iterable[Tuple[str, str]],
                   nodes: Iterable[str]) -> Dict[str, object]:
    """Edge-level confusion over all ordered off-diagonal node pairs.

    The negative set is every ordered pair that is not a reference edge, made
    explicit so AUROC/Jaccard are not at the mercy of an implicit negative count.
    """
    nodes = sorted(set(nodes))
    pred = {tuple(e) for e in pred_edges}
    ref = {tuple(e) for e in ref_edges}
    universe = set(_ordered_offdiag_pairs(nodes))
    pred &= universe
    ref &= universe
    tp = len(pred & ref)
    fp = len(pred - ref)
    fn = len(ref - pred)
    tn = len(universe - (pred | ref))
    prec = tp / (tp + fp) if (tp + fp) else 0.0
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    jacc = tp / (tp + fp + fn) if (tp + fp + fn) else 0.0
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn, "n_pairs": len(universe),
            "precision": prec, "recall": rec, "f1": f1, "jaccard": jacc}


def pr_curve_points(y_true: np.ndarray, y_score: np.ndarray) -> Dict[str, np.ndarray]:
    """Threshold-free PR curve points + AUPRC for the edge-candidate scores."""
    from sklearn.metrics import precision_recall_curve, average_precision_score
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, float)
    precision, recall, thresholds = precision_recall_curve(y_true, y_score)
    auprc = float(average_precision_score(y_true, y_score)) if y_true.any() else float("nan")
    return {"precision": precision, "recall": recall, "thresholds": thresholds, "auprc": auprc}


def bootstrap_edge_metric(y_true: np.ndarray, y_score: np.ndarray,
                          metric: Callable[[np.ndarray, np.ndarray], float] = None,
                          n_boot: int = 1000, seed: int = 0,
                          ci: float = 0.95) -> Dict[str, float]:
    """Bootstrap CI over the candidate-edge set for a score->metric.

    Resamples the off-diagonal candidate pairs (the unit of a small lineage graph)
    with replacement, so a 3-4 edge reference reports an interval. Defaults to
    AUPRC, the threshold-free metric most informative on sparse graphs.
    """
    from sklearn.metrics import average_precision_score, roc_auc_score
    if metric is None:
        metric = average_precision_score
    y_true = np.asarray(y_true).astype(int)
    y_score = np.asarray(y_score, float)
    n = len(y_true)
    rng = np.random.default_rng(seed)
    point = float(metric(y_true, y_score)) if len(set(y_true.tolist())) > 1 else float("nan")
    boots = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yt, ys = y_true[idx], y_score[idx]
        if len(set(yt.tolist())) < 2:
            continue
        boots.append(float(metric(yt, ys)))
    if not boots:
        return {"point": point, "lo": float("nan"), "hi": float("nan"),
                "mean": float("nan"), "n_boot_valid": 0}
    boots = np.asarray(boots)
    a = (1 - ci) / 2
    return {"point": point, "mean": float(boots.mean()),
            "lo": float(np.quantile(boots, a)), "hi": float(np.quantile(boots, 1 - a)),
            "n_boot_valid": int(boots.size)}
