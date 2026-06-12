"""
lineage_graphsim_sctimebench.py
scTimeBench-style graph similarity metrics for Lineage Fidelity evaluation.

Ports the logic from scTimeBench graph_sim module:
  utils.py, confusion_matrix.py, jaccard_similarity.py, base.py,
  methods/correlation/run.py

Two threshold criteria matching scTimeBench defaults:
  simple    -- direct weighted adjacency vs binary reference
  all_paths -- Floyd-Warshall closure on both sides

scTimeBench parity notes:
  Diagonal: floyd_warshall_closure sets diagonal=1 (dist[i][i]=0 < inf).
    Full-matrix flatten used throughout; diagonal NOT manually excluded.
  AUC: W_pred diagonal zeroed before roc_auc_score/average_precision_score,
    matching get_threshold_roc/prc in scTimeBench confusion_matrix.py.
    Threshold selection uses unmodified scores (no diagonal zeroing).
  all_paths binary pred: threshold original W_pred (not W_reach).
    W_reach is for threshold selection only (modified_floyd_warshall step).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from typing import Any, Dict, List, Optional, Tuple


PREDICTION_LABEL_SOURCE_SCTIMEBENCH_ZERO_FILL = (
    "reference_graph_nodes_sctimebench_zero_filled"
)


def _coerce_axis_to_str(values) -> List[str]:
    return [str(v) for v in values]


def _expanded_label_examples(labels: List[str]) -> List[str]:
    return sorted({label for label in labels if label.startswith("stage_")})[:20]


def _string_labeled_matrix(matrix: pd.DataFrame, matrix_name: str) -> pd.DataFrame:
    """Return a copy with string labels and fail if coercion creates duplicates."""
    labeled = matrix.copy()
    labeled.index = pd.Index(_coerce_axis_to_str(labeled.index), name=matrix.index.name)
    labeled.columns = pd.Index(_coerce_axis_to_str(labeled.columns), name=matrix.columns.name)
    if not labeled.index.is_unique:
        duplicates = sorted(set(labeled.index[labeled.index.duplicated()].tolist()))
        raise ValueError(
            f"{matrix_name} has duplicate source-state labels after string coercion: "
            f"{duplicates[:20]}"
        )
    if not labeled.columns.is_unique:
        duplicates = sorted(set(labeled.columns[labeled.columns.duplicated()].tolist()))
        raise ValueError(
            f"{matrix_name} has duplicate target-state labels after string coercion: "
            f"{duplicates[:20]}"
        )
    return labeled


def validate_predicted_label_space(
    predicted_matrix: pd.DataFrame,
    ref_node_ids: list,
    *,
    matrix_name: str = "state_transition_matrix.csv",
) -> Dict[str, Any]:
    """Validate a method STM against the reference graph's coarse nodes.

    The reference graph defines the canonical state space.  The method matrix may
    contain extra excluded coarse states such as ``ambiguous``, but it must not
    use expanded ``stage_*`` labels.  Missing reference nodes are aligned the
    same way as scTimeBench: the evaluator reindexes predictions onto the full
    reference node set and fills absent rows/columns with zeros.
    """
    ref_nodes = _coerce_axis_to_str(ref_node_ids)
    ref_set = set(ref_nodes)
    labeled = _string_labeled_matrix(predicted_matrix, matrix_name)
    row_labels = _coerce_axis_to_str(labeled.index)
    col_labels = _coerce_axis_to_str(labeled.columns)
    row_set = set(row_labels)
    col_set = set(col_labels)

    expanded = _expanded_label_examples(row_labels + col_labels)
    if expanded:
        raise ValueError(
            f"{matrix_name} uses expanded/stage labels, but official lineage "
            "evaluation requires final_milestone_label_coarse labels from the "
            f"reference graph. Examples: {expanded}"
        )

    missing_rows = [node for node in ref_nodes if node not in row_set]
    missing_cols = [node for node in ref_nodes if node not in col_set]
    missing_set = set(missing_rows) | set(missing_cols)
    missing_union = [node for node in ref_nodes if node in missing_set]

    return {
        "reference_nodes": ref_nodes,
        "extra_source_states": sorted(row_set - ref_set),
        "extra_target_states": sorted(col_set - ref_set),
        "missing_source_reference_nodes": missing_rows,
        "missing_target_reference_nodes": missing_cols,
        "missing_reference_nodes": missing_union,
        "n_missing_source_reference_nodes": len(missing_rows),
        "n_missing_target_reference_nodes": len(missing_cols),
        "n_missing_reference_nodes": len(missing_union),
        "n_reference_nodes": len(ref_nodes),
        "n_source_states_in_matrix": len(row_labels),
        "n_target_states_in_matrix": len(col_labels),
        "label_space": "final_milestone_label_coarse/reference_graph_nodes",
        "alignment_policy": "sctimebench_zero_fill_missing_reference_nodes",
        "zero_filled_missing_reference_nodes": bool(missing_union),
    }


# ---------------------------------------------------------------------------
# Graph algorithms
# ---------------------------------------------------------------------------

def modified_floyd_warshall(W: np.ndarray) -> np.ndarray:
    n = W.shape[0]
    R = W.astype(np.float64, copy=True)
    for k in range(n):
        via_k = np.minimum(R[:, k:k+1], R[k:k+1, :])
        R = np.maximum(R, via_k)
    return R


def floyd_warshall_closure(B: np.ndarray) -> np.ndarray:
    n = B.shape[0]
    C = (B > 0).astype(np.float64)
    for k in range(n):
        C = np.clip(C + np.outer(C[:, k], C[k, :]), 0.0, 1.0)
    C = C.astype(np.int32)
    np.fill_diagonal(C, 1)
    return C


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------

def _select_threshold(
    scores_flat: np.ndarray,
    labels_flat: np.ndarray,
    prc_threshold: bool = True,
) -> Tuple[float, str]:
    from sklearn.metrics import precision_recall_curve, roc_curve
    thresh_type = "prc" if prc_threshold else "roc"
    y_true  = labels_flat.astype(int)
    y_score = scores_flat.astype(np.float64)
    if len(np.unique(y_true)) < 2:
        finite = y_score[np.isfinite(y_score)]
        return (float(finite.max()) if finite.size > 0 else float("nan")), thresh_type
    try:
        if prc_threshold:
            precision, recall, thresholds = precision_recall_curve(y_true, y_score)
            if len(thresholds) == 0:
                return float("nan"), thresh_type
            obj = precision[:-1] + recall[:-1]
            return float(thresholds[int(np.argmax(obj))]), thresh_type
        else:
            fpr, tpr, thresholds = roc_curve(y_true, y_score)
            return float(thresholds[int(np.argmax(tpr - fpr))]), thresh_type
    except Exception:
        finite = y_score[np.isfinite(y_score)]
        return (float(finite.max()) if finite.size > 0 else float("nan")), thresh_type


def _auc_metrics(y_true_flat: np.ndarray, y_score_flat: np.ndarray) -> Dict[str, float]:
    from sklearn.metrics import roc_auc_score, average_precision_score
    if len(np.unique(y_true_flat)) < 2:
        return {"auc_roc": float("nan"), "auc_prc": float("nan")}
    try:
        auc_roc = float(roc_auc_score(y_true_flat.astype(int), y_score_flat.astype(np.float64)))
    except Exception:
        auc_roc = float("nan")
    try:
        auc_prc = float(average_precision_score(y_true_flat.astype(int), y_score_flat.astype(np.float64)))
    except Exception:
        auc_prc = float("nan")
    return {"auc_roc": auc_roc, "auc_prc": auc_prc}


def _confusion_metrics(pred_flat: np.ndarray, ref_flat: np.ndarray) -> Dict[str, float]:
    p = pred_flat.astype(int)
    r = ref_flat.astype(int)
    tp = int(((p == 1) & (r == 1)).sum())
    fp = int(((p == 1) & (r == 0)).sum())
    fn = int(((p == 0) & (r == 1)).sum())
    tn = int(((p == 0) & (r == 0)).sum())
    total = tp + fp + fn + tn
    accuracy  = (tp + tn) / total if total > 0 else float("nan")
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall    = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2.0 * precision * recall / (precision + recall) if (precision + recall) > 0.0 else 0.0
    return {"accuracy": float(accuracy), "precision": float(precision),
            "recall": float(recall), "f1": float(f1)}


def _jaccard_flat(pred_flat: np.ndarray, ref_flat: np.ndarray) -> float:
    p = pred_flat.astype(bool)
    r = ref_flat.astype(bool)
    intersection = int((p & r).sum())
    union        = int((p | r).sum())
    return float(intersection / union) if union > 0 else float("nan")


# ---------------------------------------------------------------------------
# Per-criterion metrics
# ---------------------------------------------------------------------------

def _compute_simple_metrics(
    W_pred: np.ndarray, B_ref: np.ndarray,
    edge_threshold: float, auto_threshold: bool, prc_threshold: bool,
) -> Dict[str, Any]:
    W_pred_dz = W_pred.copy()
    np.fill_diagonal(W_pred_dz, 0.0)
    y_true = B_ref.flatten().astype(int)
    auc_bundle = _auc_metrics(y_true, W_pred_dz.flatten().astype(np.float64))
    y_score_thresh = W_pred.flatten().astype(np.float64)
    if auto_threshold:
        thresh, thresh_type = _select_threshold(y_score_thresh, y_true, prc_threshold=prc_threshold)
    else:
        thresh, thresh_type = float(edge_threshold), ("prc" if prc_threshold else "roc")
    pred_binary = (W_pred >= thresh).astype(np.int32)
    np.fill_diagonal(pred_binary, 0)   # mirrors scTimeBench: source_id != target_id
    cm_bundle = _confusion_metrics(pred_binary.flatten(), B_ref.flatten())
    jaccard   = _jaccard_flat(pred_binary.flatten(), B_ref.flatten())
    return {"threshold": thresh, "threshold_type": thresh_type,
            **cm_bundle, **auc_bundle, "jaccard_similarity": jaccard}


def _compute_all_paths_metrics(
    W_pred_original: np.ndarray, W_reach: np.ndarray, B_reach: np.ndarray,
    edge_threshold: float, auto_threshold: bool, prc_threshold: bool,
) -> Dict[str, Any]:
    y_true = B_reach.flatten().astype(int)
    W_pred_dz = W_pred_original.copy()
    np.fill_diagonal(W_pred_dz, 0.0)
    auc_bundle = _auc_metrics(y_true, W_pred_dz.flatten().astype(np.float64))
    y_score_thresh = W_reach.flatten().astype(np.float64)
    if auto_threshold:
        thresh, thresh_type = _select_threshold(y_score_thresh, y_true, prc_threshold=prc_threshold)
    else:
        thresh, thresh_type = float(edge_threshold), ("prc" if prc_threshold else "roc")
    pred_binary_direct = (W_pred_original >= thresh).astype(np.int32)
    np.fill_diagonal(pred_binary_direct, 0)   # mirrors scTimeBench: source_id != target_id
    pred_binary_closed = floyd_warshall_closure(pred_binary_direct)
    cm_bundle = _confusion_metrics(pred_binary_closed.flatten(), B_reach.flatten())
    jaccard   = _jaccard_flat(pred_binary_closed.flatten(), B_reach.flatten())
    return {"threshold": thresh, "threshold_type": thresh_type,
            **cm_bundle, **auc_bundle, "jaccard_similarity": jaccard}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def align_predicted_to_reference(
    predicted_matrix: pd.DataFrame,
    ref_node_ids: list,
    *,
    matrix_name: str = "state_transition_matrix.csv",
) -> np.ndarray:
    ref_nodes = _coerce_axis_to_str(ref_node_ids)
    labeled = _string_labeled_matrix(predicted_matrix, matrix_name)
    validate_predicted_label_space(
        labeled,
        ref_nodes,
        matrix_name=matrix_name,
    )
    aligned = labeled.reindex(index=ref_nodes, columns=ref_nodes, fill_value=0.0)
    return aligned.values.astype(np.float64)


def compute_graph_sim_metrics(
    predicted_matrix: pd.DataFrame,
    reference_matrix: pd.DataFrame,
    ref_node_ids: list,
    *,
    threshold_criteria: Optional[list] = None,
    auto_threshold: bool = True,
    edge_threshold: float = 0.1,
    prc_threshold: bool = True,
) -> Dict[str, Any]:
    if threshold_criteria is None:
        threshold_criteria = ["simple", "all_paths"]
    W_pred = align_predicted_to_reference(predicted_matrix, ref_node_ids)
    B_ref  = reference_matrix.reindex(
        index=ref_node_ids, columns=ref_node_ids, fill_value=0
    ).values.astype(np.int32)
    result: Dict[str, Any] = {}
    for criterion in threshold_criteria:
        if criterion == "simple":
            result["single_step"] = _compute_simple_metrics(
                W_pred, B_ref, edge_threshold, auto_threshold, prc_threshold
            )
        elif criterion == "all_paths":
            W_reach = modified_floyd_warshall(W_pred)
            B_reach = floyd_warshall_closure(B_ref)
            result["multi_step"] = _compute_all_paths_metrics(
                W_pred, W_reach, B_reach, edge_threshold, auto_threshold, prc_threshold
            )
        else:
            raise ValueError(f"Unknown criterion {criterion!r}. Expected simple or all_paths.")
    return result


def graph_sim_from_baseline_df(baseline_df, reference_matrix, ref_node_ids, **kwargs):
    return compute_graph_sim_metrics(
        predicted_matrix=baseline_df,
        reference_matrix=reference_matrix,
        ref_node_ids=ref_node_ids,
        **kwargs,
    )
