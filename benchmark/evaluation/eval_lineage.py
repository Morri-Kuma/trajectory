"""
eval_lineage.py
Lineage Fidelity evaluator — the first active core evaluator.
Framework reference: experimental framework v2.md §9, §14 Step 2

Lineage Fidelity is the ONLY active benchmark dimension for the current
WOT vs CellRank2 stage. This evaluator is called by eval_dispatch.py
after each method's adapter produces state_transition_matrix.csv
and lineage_graph_edges.csv.

Metrics (per v2 §9.5):
  - AUROC
  - AUPRC
  - Jaccard Similarity
  - Single-step lineage recovery
  - Multi-step lineage recovery

Baseline (per v2 §9.6):
  A correlation-based baseline analogous to scTimeBench's lineage baseline
  is included so that results are judged relative to a non-trivial baseline.

Reference requirement (per v2 §9.3):
  A formal Lineage Fidelity benchmark requires a defined cell-state system
  and a frozen reference lineage graph. These are benchmark prerequisites
  that must be explicitly frozen before this evaluator is fully activated.

Reference graph format (scGPT v1):
  JSON with top-level keys: _meta, nodes, edges.
  Each edge carries: source, target, weight, confidence, source_status,
  target_status. edge_confidence_mode controls which edges are used as
  ground truth (see load_reference_graph() for details).
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, List


# ------------------------------------------------------------------
# Reference graph loader
# ------------------------------------------------------------------

def load_reference_graph(
    reference_graph_path: str,
    edge_confidence_mode: str = "all",
    exclude_uncertain_states: bool = False,
) -> Tuple[pd.DataFrame, set, List[str]]:
    """
    Load and parse a reference lineage graph JSON (scGPT v1 format).

    The scGPT v1 JSON has top-level keys: _meta, nodes, edges.
    Each edge carries: source, target, weight, confidence,
    source_status, target_status.

    Parameters
    ----------
    reference_graph_path : str
        Path to the reference graph JSON (e.g. scgpt_reference_graph_v1.json).
    edge_confidence_mode : str
        Controls which edges are included as ground-truth positives:

        "all"
            Use all edges regardless of confidence (most permissive).
            GSE230659 v1: 42 edges (3 high + 28 medium + 11 low).

        "medium_and_above"
            Exclude edges where confidence == "low".
            Removes edges that involve at least one uncertain-status state.
            GSE230659 v1: 31 edges (3 high + 28 medium).
            RECOMMENDED for initial benchmark runs.

        "high_only"
            Use only edges where confidence == "high".
            Both endpoints must be confirmed-status states.
            GSE230659 v1: 3 edges. Very sparse; use for signal checks only.

    exclude_uncertain_states : bool
        If True, additionally exclude any edge where source or target node
        has status == "uncertain", regardless of confidence_mode.
        Equivalent to removing PS_11 and PS_12 from the reference in the
        current GSE230659 v1 graph.

    Returns
    -------
    reference_matrix : pd.DataFrame
        Binary adjacency matrix.  Rows = source states, columns = target states.
        Value 1 = edge exists (passed the filter), 0 = no edge.
        Index and columns are the full node-id list from the graph.
    reference_edges : set
        Set of (source_id, target_id) tuples for edges that passed the filter.
    node_ids : list[str]
        Ordered list of all state IDs in the graph (from the nodes block).
    """
    _confidence_rank = {"high": 3, "medium": 2, "low": 1}
    _min_rank = {
        "all": 0,
        "medium_and_above": 2,
        "high_only": 3,
    }

    if edge_confidence_mode not in _min_rank:
        raise ValueError(
            f"edge_confidence_mode must be one of {list(_min_rank.keys())}, "
            f"got {edge_confidence_mode!r}"
        )
    min_rank = _min_rank[edge_confidence_mode]

    with open(reference_graph_path, encoding="utf-8") as f:
        graph = json.load(f)

    # Extract nodes
    nodes = graph.get("nodes", [])
    node_ids = [n["id"] for n in nodes]
    node_status = {n["id"]: n.get("status", "unknown") for n in nodes}

    # Build binary reference matrix (all zeros initially)
    reference_matrix = pd.DataFrame(
        0, index=node_ids, columns=node_ids, dtype=int
    )
    reference_edges: set = set()

    for edge in graph.get("edges", []):
        src = edge.get("source", "")
        tgt = edge.get("target", "")
        conf = edge.get("confidence", "low")

        # Confidence filter
        if _confidence_rank.get(conf, 0) < min_rank:
            continue

        # Optional uncertain-state exclusion
        if exclude_uncertain_states:
            if (node_status.get(src) == "uncertain" or
                    node_status.get(tgt) == "uncertain"):
                continue

        # Guard: skip edges whose nodes are not in the node list
        if src not in reference_matrix.index or tgt not in reference_matrix.columns:
            continue

        reference_matrix.loc[src, tgt] = 1
        reference_edges.add((src, tgt))

    return reference_matrix, reference_edges, node_ids


# ------------------------------------------------------------------
# Metric functions
# ------------------------------------------------------------------

def compute_auroc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Compute AUROC for lineage edge recovery."""
    from sklearn.metrics import roc_auc_score
    if len(np.unique(y_true)) < 2:
        return float("nan")
    return float(roc_auc_score(y_true, y_score))


def compute_auprc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    """Compute AUPRC for lineage edge recovery."""
    from sklearn.metrics import average_precision_score
    if len(np.unique(y_true)) < 2:
        return float("nan")
    return float(average_precision_score(y_true, y_score))


def compute_jaccard(predicted_edges: set, reference_edges: set) -> float:
    """Compute Jaccard similarity between predicted and reference lineage edges."""
    if not predicted_edges and not reference_edges:
        return float("nan")
    intersection = len(predicted_edges & reference_edges)
    union = len(predicted_edges | reference_edges)
    return float(intersection / union) if union > 0 else float("nan")


def compute_single_step_recovery(predicted_matrix: pd.DataFrame,
                                  reference_matrix: pd.DataFrame) -> float:
    """
    Single-step lineage recovery: fraction of reference edges in the top-k
    predicted edges (where k = number of reference edges).
    """
    if reference_matrix.empty or predicted_matrix.empty:
        return float("nan")
    ref_edges = set(zip(*np.where(reference_matrix.values > 0)))
    n_ref = len(ref_edges)
    if n_ref == 0:
        return float("nan")
    flat_pred = predicted_matrix.values.flatten()
    threshold = np.sort(flat_pred)[-n_ref] if n_ref <= len(flat_pred) else 0
    pred_edges = set(zip(*np.where(predicted_matrix.values >= threshold)))
    recovery = len(ref_edges & pred_edges) / n_ref
    return float(recovery)


def compute_multi_step_recovery(predicted_matrix: pd.DataFrame,
                                 reference_matrix: pd.DataFrame,
                                 steps: int = 2) -> float:
    """
    Multi-step lineage recovery: same as single-step but evaluated on
    multi-hop paths through the reference lineage.
    """
    if reference_matrix.empty or predicted_matrix.empty:
        return float("nan")
    # Multi-step reference graph via matrix exponentiation.
    ref_vals = reference_matrix.values
    multi_step_ref = np.linalg.matrix_power(ref_vals > 0, steps).astype(float)
    multi_step_ref_df = pd.DataFrame(multi_step_ref,
                                      index=reference_matrix.index,
                                      columns=reference_matrix.columns)
    return compute_single_step_recovery(predicted_matrix, multi_step_ref_df)


# ------------------------------------------------------------------
# Correlation baseline
# ------------------------------------------------------------------

def compute_correlation_baseline(adata,
                                  reference_matrix: pd.DataFrame) -> dict:
    """
    Correlation-based baseline for Lineage Fidelity (per v2 §9.6).

    Analogous to scTimeBench's lineage baseline: uses pairwise gene expression
    correlation between cell states as a naive predictor of lineage connectivity,
    then evaluates against the reference lineage using the same metrics.
    """
    if adata is None or reference_matrix.empty:
        return {metric: None for metric in
                ["auroc", "auprc", "jaccard_similarity",
                 "single_step_recovery", "multi_step_recovery"]}

    # TODO: compute mean gene expression per cell state and build correlation matrix.
    # state_means = ...
    # corr_matrix = np.corrcoef(state_means)
    # Evaluate corr_matrix against reference_matrix using the same metrics.

    return {
        "auroc": None,
        "auprc": None,
        "jaccard_similarity": None,
        "single_step_recovery": None,
        "multi_step_recovery": None,
        "note": "Correlation baseline not yet implemented. Requires frozen cell-state system.",
    }


# ------------------------------------------------------------------
# Main evaluation function
# ------------------------------------------------------------------

def run_lineage_evaluation(
    state_transition_matrix_path: str,
    lineage_graph_edges_path: str,
    output_dir: str,
    reference_graph_path: Optional[str] = None,
    adata=None,
    edge_confidence_mode: str = "all",
    exclude_uncertain_states: bool = False,
) -> dict:
    """
    Evaluate Lineage Fidelity for one method × scenario run.

    Parameters
    ----------
    state_transition_matrix_path : str
        Path to state_transition_matrix.csv produced by the method adapter.
    lineage_graph_edges_path : str
        Path to lineage_graph_edges.csv produced by the method adapter.
    output_dir : str
        Directory where lineage_metrics.json will be written.
    reference_graph_path : str, optional
        Path to the reference lineage graph JSON (scGPT v1 format: _meta,
        nodes, edges with confidence field). If None, metric computation is
        deferred and a status note is recorded (backward-compatible default).
    adata : AnnData, optional
        Required for the correlation baseline computation.
    edge_confidence_mode : str
        Controls which reference graph edges count as ground-truth positives.
        Options: "all", "medium_and_above" (recommended), "high_only".
        See load_reference_graph() for full documentation.
        Defaults to "all" for backward compatibility with configs that do not
        set this field.
    exclude_uncertain_states : bool
        If True, exclude edges involving uncertain-status nodes from the
        reference. See load_reference_graph() for details. Default False.

    Returns
    -------
    dict with keys: auroc, auprc, jaccard_similarity, single_step_recovery,
                    multi_step_recovery, baseline, status,
                    edge_confidence_mode, n_reference_edges
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load predicted outputs
    stm_path = Path(state_transition_matrix_path)
    edges_path = Path(lineage_graph_edges_path)

    predicted_matrix = pd.DataFrame()
    predicted_edges = set()

    if stm_path.exists() and stm_path.stat().st_size > 0:
        predicted_matrix = pd.read_csv(stm_path, index_col=0)

    if edges_path.exists() and edges_path.stat().st_size > 0:
        edges_df = pd.read_csv(edges_path)
        if not edges_df.empty:
            predicted_edges = set(
                zip(edges_df["source_state"], edges_df["target_state"])
            )

    # Check reference availability
    if reference_graph_path is None:
        status = (
            "deferred: reference lineage graph not yet provided. "
            "Per v2 §9.3, a frozen reference graph is required before "
            "Lineage Fidelity metrics can be computed."
        )
        metrics = {
            "auroc": None,
            "auprc": None,
            "jaccard_similarity": None,
            "single_step_recovery": None,
            "multi_step_recovery": None,
            "baseline": None,
            "status": status,
            "edge_confidence_mode": edge_confidence_mode,
            "n_reference_edges": None,
        }
        metrics_path = out_dir / "lineage_metrics.json"
        with open(metrics_path, "w") as f:
            json.dump(metrics, f, indent=2)
        print(f"[eval_lineage] {status}")
        return metrics

    # Load and parse reference graph (scGPT v1 JSON format).
    # load_reference_graph() handles confidence filtering and returns a
    # binary adjacency DataFrame indexed by state IDs.
    reference_matrix, reference_edges, ref_node_ids = load_reference_graph(
        reference_graph_path,
        edge_confidence_mode=edge_confidence_mode,
        exclude_uncertain_states=exclude_uncertain_states,
    )
    print(
        f"[eval_lineage] Reference graph loaded: {len(ref_node_ids)} nodes, "
        f"{len(reference_edges)} edges "
        f"(mode={edge_confidence_mode!r})"
    )

    # Determine whether the method produced real predictions.
    # predicted_matrix is empty when the adapter wrote scaffold/placeholder outputs
    # (e.g. WOT failed and called _write_scaffold_outputs(), or CellRank2 is
    # scaffold-only). In that case the metrics are all None and we must NOT call
    # this run "completed" — that would hide the fact that no real predictions
    # were made.
    has_predictions = not predicted_matrix.empty or bool(predicted_edges)

    # Compute metrics
    metrics = {}

    if has_predictions and not reference_matrix.empty:
        # Both method output and reference are available — compute real metrics.
        all_states = list(
            set(predicted_matrix.index) | set(reference_matrix.index)
        )
        pred_aligned = predicted_matrix.reindex(
            index=all_states, columns=all_states, fill_value=0.0
        )
        ref_aligned = reference_matrix.reindex(
            index=all_states, columns=all_states, fill_value=0
        )
        y_true = ref_aligned.values.flatten().astype(int)
        y_score = pred_aligned.values.flatten().astype(float)

        metrics["auroc"] = compute_auroc(y_true, y_score)
        metrics["auprc"] = compute_auprc(y_true, y_score)
        metrics["jaccard_similarity"] = compute_jaccard(predicted_edges, reference_edges)
        metrics["single_step_recovery"] = compute_single_step_recovery(
            pred_aligned, ref_aligned
        )
        metrics["multi_step_recovery"] = compute_multi_step_recovery(
            pred_aligned, ref_aligned
        )
        completion_status = "completed"
    else:
        # Either the method wrote no predictions (scaffold/failure path) or the
        # reference matrix is empty.  All metrics are None, but we clearly label
        # WHY so the status is not mistaken for a successful evaluation.
        metrics.update({
            "auroc": None, "auprc": None, "jaccard_similarity": None,
            "single_step_recovery": None, "multi_step_recovery": None,
        })
        if not has_predictions:
            completion_status = (
                "completed_with_empty_predictions: "
                "reference graph was loaded successfully but the method adapter "
                "produced no state transitions (scaffold or runtime failure path). "
                "Check run_metadata.json for the method execution status."
            )
        else:
            # has_predictions but reference_matrix is empty — should not happen
            # with a valid reference graph, but guard defensively.
            completion_status = (
                "completed_with_empty_reference: "
                "method predictions exist but reference matrix is empty after "
                "applying the current edge_confidence_mode filter."
            )

    # Correlation baseline
    metrics["baseline"] = compute_correlation_baseline(adata, reference_matrix)
    metrics["status"] = completion_status
    metrics["edge_confidence_mode"] = edge_confidence_mode
    metrics["n_reference_edges"] = len(reference_edges)

    metrics_path = out_dir / "lineage_metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"[eval_lineage] Metrics written to {metrics_path}")
    return metrics
