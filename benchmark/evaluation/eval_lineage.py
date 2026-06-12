"""
eval_lineage.py
Lineage Fidelity evaluator — the first active core evaluator.
Framework reference: docs/framework/experimental_framework_v2.md §9, §14 Step 2

Lineage Fidelity is the ONLY active benchmark dimension for the current
WOT vs CellRank2 stage. This evaluator is called by eval_dispatch.py
after each method's adapter produces state_transition_matrix.csv
and lineage_graph_edges.csv.

Metric protocol (v2 §9.5, updated to scTimeBench graph similarity):
  Metrics are computed by lineage_graphsim_sctimebench.py using two
  threshold criteria that mirror scTimeBench's graph_sim module:
    - simple    (single_step): direct weighted adjacency vs binary reference
    - all_paths (multi_step):  Floyd-Warshall reachability vs transitive closure

  Each criterion produces:
    threshold, threshold_type, accuracy, precision, recall, f1,
    auc_roc, auc_prc, jaccard_similarity

  Defaults match scTimeBench:
    threshold_criteria: ["simple", "all_paths"]
    auto_threshold: True
    edge_threshold: 0.1
    prc_threshold: True

Baseline (per v2 §9.6):
  scTimeBench-style Spearman maximum-vote correlation baseline evaluated
  with the same graph-sim metric functions.

Reference requirement (per v2 §9.3):
  A formal Lineage Fidelity benchmark requires a defined cell-state system
  and a frozen reference lineage graph. These are benchmark prerequisites
  that must be explicitly frozen before this evaluator is fully activated.

Reference graph format:
  JSON with top-level keys: _meta, nodes, edges.
  Each edge carries: source, target, weight, confidence, source_status,
  target_status. edge_confidence_mode controls which edges are used as
  ground truth (see load_reference_graph() for details).

lineage_metrics.json schema (metric_protocol: "sctimebench_graph_sim"):
  graph_metrics.single_step.*   ← simple criterion
  graph_metrics.multi_step.*    ← all_paths criterion
  Backward-compat top-level aliases:
    auroc                 = graph_metrics.single_step.auc_roc
    auprc                 = graph_metrics.single_step.auc_prc
    jaccard_similarity    = graph_metrics.single_step.jaccard_similarity
    single_step_recovery  = graph_metrics.single_step.recall
    multi_step_recovery   = graph_metrics.multi_step.recall
  Legacy (diagnostic only, not used for ranking):
    legacy_jaccard_similarity_topk
"""

import json
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple, List
from scipy.sparse import issparse
from scipy.stats import rankdata

from benchmark.evaluation.lineage_graphsim_sctimebench import (
    PREDICTION_LABEL_SOURCE_SCTIMEBENCH_ZERO_FILL,
    compute_graph_sim_metrics,
    graph_sim_from_baseline_df,
    validate_predicted_label_space,
)


OFFICIAL_CELL_STATE_KEY = "final_milestone_label_coarse"


def _mark_run_metadata_lineage_failed(out_dir: Path, reason: str) -> None:
    """Annotate method metadata when post-hoc lineage evaluation fails."""
    metadata_path = out_dir / "run_metadata.json"
    if not metadata_path.exists():
        return
    try:
        with open(metadata_path, encoding="utf-8") as f:
            metadata = json.load(f)
        metadata["status"] = "failed_lineage_evaluation"
        prior_notes = metadata.get("notes") or ""
        metadata["notes"] = (
            f"{prior_notes}\n{reason}".strip()
            if prior_notes else reason
        )
        with open(metadata_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)
    except Exception as exc:
        print(
            "[eval_lineage] Warning: failed to update run_metadata.json "
            f"after lineage evaluation failure: {exc}"
        )


def _mark_run_metadata_ground_truth(out_dir: Path, ground_truth: dict) -> None:
    """Attach ground-truth provider metadata to run_metadata.json when present."""
    if not ground_truth:
        return
    metadata_path = out_dir / "run_metadata.json"
    if not metadata_path.exists():
        return
    try:
        with open(metadata_path, encoding="utf-8") as f:
            metadata = json.load(f)
        metadata["ground_truth"] = ground_truth
        with open(metadata_path, "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)
    except Exception as exc:
        print(
            "[eval_lineage] Warning: failed to update run_metadata.json "
            f"with ground-truth metadata: {exc}"
        )


def _nonfinite_matrix_message(matrix: pd.DataFrame, matrix_name: str) -> str:
    """Return a compact diagnostic message for NaN/Inf values in a matrix."""
    values = matrix.values.astype(float)
    bad = ~np.isfinite(values)
    rows = list(matrix.index)
    cols = list(matrix.columns)
    bad_positions = np.argwhere(bad)
    examples = [
        f"{rows[i]}->{cols[j]}={float(values[i, j])!r}"
        for i, j in bad_positions[:10]
    ]
    bad_rows = [rows[i] for i in np.where(bad.any(axis=1))[0]]
    return (
        f"{matrix_name} contains {int(bad.sum())} non-finite values "
        f"across {len(bad_rows)} source rows. "
        f"Bad source rows: {bad_rows[:20]}. "
        f"Examples: {examples}. "
        "Fix the method output before computing AUROC/AUPRC; for source "
        "states with no valid source-cell transitions, write an all-zero row "
        "and record the unsupported states in diagnostics."
    )


def _load_reference_graph_meta(reference_graph_path: str) -> dict:
    with open(reference_graph_path, encoding="utf-8") as f:
        graph = json.load(f)
    meta = graph.get("_meta") or {}
    return meta if isinstance(meta, dict) else {}


def _is_official_silver_contract(reference_meta: dict, ground_truth: dict) -> bool:
    return (
        reference_meta.get("label_mode") == "official_silver"
        or ground_truth.get("label_mode") == "official_silver"
        or reference_meta.get("label_type") == "frozen_silver_standard"
        or ground_truth.get("label_type") == "frozen_silver_standard"
    )


def _validate_reference_state_key(
    *,
    reference_meta: dict,
    ground_truth: dict,
    cell_state_key: Optional[str],
) -> str:
    """Resolve and validate the canonical cell-state key for lineage metrics."""
    reference_state_key = (
        reference_meta.get("state_key")
        or ground_truth.get("state_key")
        or cell_state_key
    )
    if not reference_state_key:
        raise ValueError(
            "Lineage evaluation requires a cell_state_key from the reference "
            "graph or ground-truth provider."
        )
    if cell_state_key and cell_state_key != reference_state_key:
        raise ValueError(
            f"cell_state_key={cell_state_key!r} does not match reference graph "
            f"state_key={reference_state_key!r}."
        )
    if (
        _is_official_silver_contract(reference_meta, ground_truth)
        and reference_state_key != OFFICIAL_CELL_STATE_KEY
    ):
        raise ValueError(
            f"official_silver lineage evaluation must use "
            f"{OFFICIAL_CELL_STATE_KEY!r}, got {reference_state_key!r}."
        )
    return str(reference_state_key)


def _expanded_edge_labels(predicted_edges: set) -> list:
    labels = {
        str(label)
        for edge in predicted_edges
        for label in edge
        if str(label).startswith("stage_")
    }
    return sorted(labels)[:20]


# ------------------------------------------------------------------
# Reference graph loader
# ------------------------------------------------------------------

def load_reference_graph(
    reference_graph_path: str,
    edge_confidence_mode: str = "all",
    exclude_uncertain_states: bool = False,
) -> Tuple[pd.DataFrame, set, List[str]]:
    """
    Load and parse a reference lineage graph JSON.

    The provider JSON has top-level keys: _meta, nodes, edges.
    Each edge carries: source, target, weight, confidence,
    source_status, target_status.

    Parameters
    ----------
    reference_graph_path : str
        Path to the reference graph JSON.
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
# Legacy metric helpers (kept for backward compat / diagnostic use)
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


# ------------------------------------------------------------------
# Correlation baseline
# ------------------------------------------------------------------

def _state_mean_expression(adata, cell_state_key: str,
                            states: list,
                            chunk_size: int = 4096) -> pd.DataFrame:
    """
    Compute mean gene expression per cell state in a memory-friendly way.

    Works whether ``adata.X`` is a sparse CSR, dense ndarray, or a backed
    h5ad ``_CSRDataset`` (AnnData backed='r'). Iterates over cells in chunks
    so we never materialise the full n_cells × n_genes matrix in RAM, which
    is essential at the 75k×23k GSE230659 scale.

    Returns
    -------
    pd.DataFrame indexed by state id, columns = gene indices (positional).
    """
    import scipy.sparse as sp

    cell_states = adata.obs[cell_state_key].astype(str).values
    n_cells = adata.n_obs
    n_genes = adata.n_vars
    state_to_idx = {str(s): i for i, s in enumerate(states)}

    sums = np.zeros((len(states), n_genes), dtype=np.float64)
    counts = np.zeros(len(states), dtype=np.int64)

    X = adata.X
    for start in range(0, n_cells, chunk_size):
        end = min(start + chunk_size, n_cells)
        chunk = X[start:end]
        if sp.issparse(chunk):
            chunk_dense = chunk.toarray()
        else:
            chunk_dense = np.asarray(chunk)
        chunk_states = cell_states[start:end]
        for s_str, i in state_to_idx.items():
            mask = chunk_states == s_str
            if not mask.any():
                continue
            sums[i] += chunk_dense[mask].sum(axis=0)
            counts[i] += int(mask.sum())

    means = np.zeros_like(sums)
    nonzero = counts > 0
    means[nonzero] = sums[nonzero] / counts[nonzero, None]

    return pd.DataFrame(means, index=list(states))


def _resolve_time_key(adata, time_key: Optional[str]) -> Optional[str]:
    """Resolve the time column used by the scTimeBench-style baseline."""
    if time_key and time_key in adata.obs.columns:
        return time_key
    for candidate in ("timepoint", "abs_day", "time_label", "day"):
        if candidate in adata.obs.columns:
            return candidate
    return None


def _sorted_timepoints(values) -> list:
    """Sort timepoints numerically when possible, otherwise lexicographically."""
    unique = pd.Series(values).dropna().unique().tolist()
    try:
        return sorted(unique, key=lambda x: float(x))
    except (TypeError, ValueError):
        return sorted(unique, key=lambda x: str(x))


def _as_dense_float32(X) -> np.ndarray:
    """Materialize a matrix slice as dense float32 for correlation math."""
    X = X.toarray() if issparse(X) else np.asarray(X)
    return X.astype(np.float32, copy=False)


def _rank_zscore_cells(X) -> np.ndarray:
    """Rank genes within each cell and z-score rows for Spearman correlation."""
    ranked = rankdata(X, axis=1, method="average").astype(np.float32, copy=False)
    mean = ranked.mean(axis=1, keepdims=True)
    std = ranked.std(axis=1, keepdims=True)
    std[std == 0] = 1.0
    return (ranked - mean) / std


def _build_pr_threshold_edges(score_df: pd.DataFrame,
                              reference_edges: set) -> Tuple[pd.DataFrame, set, float]:
    """
    Build a sparse graph from a weighted matrix using an automatic PR threshold.

    Candidate thresholds are unique finite off-diagonal scores. The selected
    threshold maximizes precision + recall against the reference edge set.
    """
    vals = score_df.values.astype(float, copy=True)
    np.fill_diagonal(vals, -np.inf)
    finite_scores = np.unique(vals[np.isfinite(vals)])
    if finite_scores.size == 0:
        empty = pd.DataFrame(columns=["source_state", "target_state", "weight"])
        return empty, set(), float("nan")

    states = list(score_df.index)
    ref = set(reference_edges)
    best_threshold = float(finite_scores.max())
    best_objective = -np.inf
    best_edges: set = set()

    for threshold in finite_scores:
        pred = {
            (states[i], states[j])
            for i, j in zip(*np.where(vals >= threshold))
        }
        if not pred:
            continue
        tp = len(pred & ref)
        precision = tp / len(pred)
        recall = tp / len(ref) if ref else 0.0
        objective = precision + recall
        if (
            objective > best_objective
            or (
                objective == best_objective
                and (len(pred), -float(threshold)) < (len(best_edges), -best_threshold)
            )
        ):
            best_objective = objective
            best_threshold = float(threshold)
            best_edges = pred

    rows = [
        {
            "source_state": src,
            "target_state": tgt,
            "weight": float(score_df.loc[src, tgt]),
        }
        for src, tgt in sorted(best_edges)
    ]
    edges_df = pd.DataFrame(rows, columns=["source_state", "target_state", "weight"])
    return edges_df, best_edges, best_threshold


def compute_correlation_baseline(adata,
                                  reference_matrix: pd.DataFrame,
                                  reference_edges: set,
                                  ref_node_ids: list,
                                  cell_state_key: Optional[str] = None,
                                  output_dir: Optional[str] = None,
                                  time_key: Optional[str] = None,
                                  source_chunk_size: int = 512) -> dict:
    """
    scTimeBench-style correlation baseline for Lineage Fidelity.

    Reproduces the scTimeBench Correlation method with
    ``correlation_method: spearmanr`` and ``averaging_method: maximum``:
    adjacent timepoint pairs, cell-level Spearman correlation, maximum score
    per target state, one best-target vote per source cell, then row-normalized
    source-state by target-state votes.

    The resulting score matrix is evaluated using the same scTimeBench
    graph-sim metric functions as method outputs.
    """
    none_result = {
        "graph_metrics": None,
        "auroc": None, "auprc": None, "jaccard_similarity": None,
        "single_step_recovery": None, "multi_step_recovery": None,
    }

    if adata is None:
        return {**none_result,
                "note": "Correlation baseline skipped: no AnnData was passed to the evaluator."}
    if reference_matrix is None or reference_matrix.empty:
        return {**none_result,
                "note": "Correlation baseline skipped: reference matrix is empty under the current edge_confidence_mode."}
    if cell_state_key is None:
        return {**none_result,
                "note": "Correlation baseline skipped: cell_state_key was not supplied to the evaluator."}
    if cell_state_key not in adata.obs.columns:
        return {**none_result,
                "note": f"Correlation baseline skipped: adata.obs does not contain {cell_state_key!r}."}

    resolved_time_key = _resolve_time_key(adata, time_key)
    if resolved_time_key is None:
        return {**none_result,
                "note": "Correlation baseline skipped: no timepoint column was available."}
    if source_chunk_size <= 0:
        source_chunk_size = 512

    states = list(reference_matrix.index)
    state_to_idx = {str(s): i for i, s in enumerate(states)}
    obs = adata.obs
    state_arr = obs[cell_state_key].astype(str).to_numpy()
    time_arr = obs[resolved_time_key].to_numpy()
    timepoints = _sorted_timepoints(time_arr)
    if len(timepoints) < 2:
        return {**none_result,
                "note": (
                    "Correlation baseline skipped: fewer than 2 timepoints are "
                    f"represented in adata under {resolved_time_key!r}."
                )}

    votes = np.zeros((len(states), len(states)), dtype=np.float64)
    n_pairs_used = 0
    n_source_cell_votes = 0

    try:
        for t0, t1 in zip(timepoints[:-1], timepoints[1:]):
            idx_t0 = np.where(time_arr == t0)[0]
            idx_t1 = np.where(time_arr == t1)[0]
            if idx_t0.size == 0 or idx_t1.size == 0:
                continue

            ct_t0 = state_arr[idx_t0]
            ct_t1 = state_arr[idx_t1]
            present_dst_states = [
                s for s in states
                if np.any(ct_t1 == str(s))
            ]
            if not present_dst_states:
                continue

            dst_cache = {}
            for dst_state in present_dst_states:
                dst_pos = idx_t1[ct_t1 == str(dst_state)]
                X_dst = _as_dense_float32(adata.X[dst_pos])
                dst_cache[str(dst_state)] = _rank_zscore_cells(X_dst)

            n_features = next(iter(dst_cache.values())).shape[1]
            n_pairs_used += 1

            for start in range(0, idx_t0.size, source_chunk_size):
                end = min(start + source_chunk_size, idx_t0.size)
                src_pos = idx_t0[start:end]
                src_states = ct_t0[start:end]
                X_src = _as_dense_float32(adata.X[src_pos])
                z_src = _rank_zscore_cells(X_src)

                best_scores = np.full(z_src.shape[0], -np.inf, dtype=np.float32)
                best_dst = np.full(z_src.shape[0], "", dtype=object)
                for dst_state, z_dst in dst_cache.items():
                    corr = (z_src @ z_dst.T) / n_features
                    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
                    corr = np.clip(corr, -1.0, 1.0)
                    scores = corr.max(axis=1)
                    update = scores > best_scores
                    best_scores[update] = scores[update]
                    best_dst[update] = dst_state

                for src_state, dst_state in zip(src_states, best_dst):
                    src_i = state_to_idx.get(str(src_state))
                    dst_i = state_to_idx.get(str(dst_state))
                    if src_i is None or dst_i is None:
                        continue
                    votes[src_i, dst_i] += 1.0
                    n_source_cell_votes += 1
    except Exception as exc:  # pragma: no cover - defensive
        return {**none_result,
                "note": f"Correlation baseline failed during Spearman vote computation: {exc!s}"}

    if n_source_cell_votes == 0:
        return {**none_result,
                "note": "Correlation baseline skipped: no valid source-cell votes were accumulated."}

    row_sums = votes.sum(axis=1, keepdims=True)
    nonzero_rows = row_sums[:, 0] > 0
    score = np.zeros_like(votes, dtype=np.float64)
    score[nonzero_rows] = votes[nonzero_rows] / row_sums[nonzero_rows]
    baseline_df = pd.DataFrame(score, index=states, columns=states, dtype=float)

    # Evaluate with scTimeBench graph-sim metrics (same as method outputs)
    graph_metrics = graph_sim_from_baseline_df(
        baseline_df=baseline_df,
        reference_matrix=reference_matrix,
        ref_node_ids=ref_node_ids,
    )
    single = graph_metrics.get("single_step", {}) or {}
    multi  = graph_metrics.get("multi_step",  {}) or {}

    metrics = {
        "graph_metrics": graph_metrics,
        # backward-compat aliases
        "auroc":               single.get("auc_roc"),
        "auprc":               single.get("auc_prc"),
        "jaccard_similarity":  single.get("jaccard_similarity"),
        "single_step_recovery": single.get("recall"),
        "multi_step_recovery":  multi.get("recall"),
        "method": "correlation_baseline",
        "baseline_style": "sctimebench_correlation",
        "correlation_method": "spearmanr",
        "averaging_method": "maximum",
        "cell_state_key": cell_state_key,
        "time_key": resolved_time_key,
        "n_states_used": int(nonzero_rows.sum()),
        "n_timepoint_pairs_used": int(n_pairs_used),
        "n_source_cell_votes": int(n_source_cell_votes),
        "note": (
            "scTimeBench-style cell-level Spearman correlation baseline with "
            "maximum target-state aggregation and row-normalized source-state votes. "
            "Evaluated with sctimebench_graph_sim metric protocol."
        ),
    }

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        baseline_df.to_csv(out / "baseline_state_transition_matrix.csv")
        # Build and save a threshold-based edge list for traceability
        edges_df, _, _ = _build_pr_threshold_edges(baseline_df, reference_edges)
        edges_df.to_csv(out / "baseline_lineage_graph_edges.csv", index=False)

    return metrics


# ------------------------------------------------------------------
# Compatibility fields  (Step 10 output naming)
# ------------------------------------------------------------------

def _enrich_lineage_metrics_compat(
    metrics: dict,
    ground_truth: dict,
    reference_graph_path: Optional[str],
) -> dict:
    """Promote Step 10 compatibility fields to top-level of lineage_metrics.json.

    The fields label_mode, provider_id, reference_graph_path, and
    compatibility_note are required at the top level so that summary scripts
    and the output naming validator can distinguish consensus / embedding_based /
    classifier_based results without inspecting the nested ground_truth block.

    This helper is intentionally non-destructive: it only adds keys that are
    not already present at the top level.
    """
    gt = ground_truth or {}
    if "label_mode" not in metrics and gt.get("label_mode"):
        metrics["label_mode"] = gt["label_mode"]
    if "provider_id" not in metrics and gt.get("provider_id"):
        metrics["provider_id"] = gt["provider_id"]
    if "reference_graph_path" not in metrics:
        graph_path = gt.get("reference_graph_path") or reference_graph_path
        if graph_path:
            metrics["reference_graph_path"] = str(graph_path)
    if "compatibility_note" not in metrics:
        metrics["compatibility_note"] = (
            "label_mode and provider_id are promoted from ground_truth for "
            "Step 10 output naming compatibility. "
            "See smoke_metadata.json for full evaluation context."
        )
    return metrics


def _write_failed_lineage_metrics(
    *,
    out_dir: Path,
    status: str,
    edge_confidence_mode: str,
    n_reference_edges: Optional[int],
    cell_state_key: Optional[str],
    time_key: Optional[str],
    ground_truth: dict,
    reference_graph_path: Optional[str],
    reason_for_metadata: str,
) -> dict:
    metrics = {
        "metric_protocol": "sctimebench_graph_sim",
        "graph_metrics": None,
        "auroc": None,
        "auprc": None,
        "jaccard_similarity": None,
        "single_step_recovery": None,
        "multi_step_recovery": None,
        "legacy_jaccard_similarity_topk": None,
        "baseline": None,
        "status": status,
        "edge_confidence_mode": edge_confidence_mode,
        "n_reference_edges": n_reference_edges,
        "cell_state_key": cell_state_key,
        "prediction_state_key": cell_state_key,
        "prediction_label_source": "invalid_or_unvalidated",
        "time_key": time_key,
        "ground_truth": ground_truth,
    }
    _enrich_lineage_metrics_compat(metrics, ground_truth, reference_graph_path)
    with open(out_dir / "lineage_metrics.json", "w", encoding="utf-8") as f:
        json.dump(metrics, f, indent=2)
    _mark_run_metadata_lineage_failed(out_dir, reason_for_metadata)
    return metrics


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
    cell_state_key: Optional[str] = None,
    time_key: Optional[str] = None,
    ground_truth: Optional[dict] = None,
) -> dict:
    """
    Evaluate Lineage Fidelity for one method × scenario run.

    Parameters
    ----------
    state_transition_matrix_path : str
        Path to state_transition_matrix.csv produced by the method adapter.
    lineage_graph_edges_path : str
        Path to lineage_graph_edges.csv produced by the method adapter.
        Retained as a diagnostic artifact; no longer the primary metric source.
    output_dir : str
        Directory where lineage_metrics.json will be written.
    reference_graph_path : str, optional
        Path to the reference lineage graph JSON (_meta, nodes, edges with
        confidence field). If None, metric computation is deferred.
    adata : AnnData, optional
        Required for the correlation baseline computation.
    edge_confidence_mode : str
        Controls which reference graph edges count as ground-truth positives.
        Options: "all", "medium_and_above" (recommended), "high_only".
    exclude_uncertain_states : bool
        If True, exclude edges involving uncertain-status nodes from the
        reference. Default False.
    cell_state_key : str, optional
        Name of the obs column holding the cell-state label.
    time_key : str, optional
        Name of the obs column holding the temporal coordinate.
    ground_truth : dict, optional
        Provider metadata for the annotation/reference graph used by this run.

    Returns
    -------
    dict — lineage_metrics.json content, with keys:
        metric_protocol, graph_metrics, status, edge_confidence_mode,
        n_reference_edges, baseline, ground_truth, cell_state_key, time_key,
        and top-level backward-compat aliases (auroc, auprc, jaccard_similarity,
        single_step_recovery, multi_step_recovery).
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ground_truth = ground_truth or {}
    _mark_run_metadata_ground_truth(out_dir, ground_truth)

    # Load predicted outputs
    stm_path   = Path(state_transition_matrix_path)
    edges_path = Path(lineage_graph_edges_path)

    predicted_matrix = pd.DataFrame()
    predicted_edges  = set()   # diagnostic only (lineage_graph_edges.csv)

    if stm_path.exists() and stm_path.stat().st_size > 0:
        predicted_matrix = pd.read_csv(stm_path, index_col=0)

    if edges_path.exists() and edges_path.stat().st_size > 0:
        edges_df = pd.read_csv(edges_path)
        if not edges_df.empty:
            predicted_edges = set(
                zip(edges_df["source_state"], edges_df["target_state"])
            )

    # ── No reference graph: defer ─────────────────────────────────────────────
    if reference_graph_path is None:
        status = (
            "deferred: reference lineage graph not yet provided. "
            "Per v2 §9.3, a frozen reference graph is required before "
            "Lineage Fidelity metrics can be computed."
        )
        metrics = {
            "metric_protocol": "sctimebench_graph_sim",
            "graph_metrics": None,
            # backward-compat aliases
            "auroc": None, "auprc": None,
            "jaccard_similarity": None,
            "single_step_recovery": None, "multi_step_recovery": None,
            "legacy_jaccard_similarity_topk": None,
            "baseline": None,
            "status": status,
            "edge_confidence_mode": edge_confidence_mode,
            "n_reference_edges": None,
            "cell_state_key": cell_state_key,
            "time_key": time_key,
            "ground_truth": ground_truth,
        }
        _enrich_lineage_metrics_compat(metrics, ground_truth, reference_graph_path)
        with open(out_dir / "lineage_metrics.json", "w") as f:
            json.dump(metrics, f, indent=2)
        print(f"[eval_lineage] {status}")
        return metrics

    # ── Load reference ────────────────────────────────────────────────────────
    reference_matrix, reference_edges, ref_node_ids = load_reference_graph(
        reference_graph_path,
        edge_confidence_mode=edge_confidence_mode,
        exclude_uncertain_states=exclude_uncertain_states,
    )
    reference_meta = _load_reference_graph_meta(reference_graph_path)
    try:
        cell_state_key = _validate_reference_state_key(
            reference_meta=reference_meta,
            ground_truth=ground_truth,
            cell_state_key=cell_state_key,
        )
    except ValueError as exc:
        reason = str(exc)
        _write_failed_lineage_metrics(
            out_dir=out_dir,
            status=f"failed_reference_state_key: {reason}",
            edge_confidence_mode=edge_confidence_mode,
            n_reference_edges=len(reference_edges),
            cell_state_key=cell_state_key,
            time_key=time_key,
            ground_truth=ground_truth,
            reference_graph_path=reference_graph_path,
            reason_for_metadata=reason,
        )
        raise
    print(
        f"[eval_lineage] Reference graph loaded: {len(ref_node_ids)} nodes, "
        f"{len(reference_edges)} edges "
        f"(mode={edge_confidence_mode!r}, state_key={cell_state_key!r})"
    )

    has_predictions = not predicted_matrix.empty or bool(predicted_edges)
    prediction_label_report: dict = {}

    # ── Compute graph-sim metrics ─────────────────────────────────────────────
    graph_metrics = None
    single_step   = {}
    multi_step    = {}

    if has_predictions and not reference_matrix.empty:
        expanded_edges = _expanded_edge_labels(predicted_edges)
        if expanded_edges:
            reason = (
                "lineage_graph_edges.csv uses expanded/stage labels, but "
                "official lineage evaluation requires "
                f"{OFFICIAL_CELL_STATE_KEY} labels. Examples: {expanded_edges}"
            )
            _write_failed_lineage_metrics(
                out_dir=out_dir,
                status=f"failed_prediction_label_space: {reason}",
                edge_confidence_mode=edge_confidence_mode,
                n_reference_edges=len(reference_edges),
                cell_state_key=cell_state_key,
                time_key=time_key,
                ground_truth=ground_truth,
                reference_graph_path=reference_graph_path,
                reason_for_metadata=reason,
            )
            raise ValueError(f"[eval_lineage] {reason}")

        try:
            prediction_label_report = validate_predicted_label_space(
                predicted_matrix,
                ref_node_ids,
                matrix_name="state_transition_matrix.csv",
            )
        except ValueError as exc:
            reason = str(exc)
            _write_failed_lineage_metrics(
                out_dir=out_dir,
                status=f"failed_prediction_label_space: {reason}",
                edge_confidence_mode=edge_confidence_mode,
                n_reference_edges=len(reference_edges),
                cell_state_key=cell_state_key,
                time_key=time_key,
                ground_truth=ground_truth,
                reference_graph_path=reference_graph_path,
                reason_for_metadata=reason,
            )
            raise

        # Validate finiteness before handing to graph-sim
        if not np.isfinite(predicted_matrix.values.astype(float)).all():
            reason = _nonfinite_matrix_message(
                predicted_matrix, "state_transition_matrix.csv"
            )
            metrics = {
                "metric_protocol": "sctimebench_graph_sim",
                "graph_metrics": None,
                "auroc": None, "auprc": None,
                "jaccard_similarity": None,
                "single_step_recovery": None, "multi_step_recovery": None,
                "legacy_jaccard_similarity_topk": None,
                "baseline": None,
                "status": f"failed_nonfinite_predictions: {reason}",
                "edge_confidence_mode": edge_confidence_mode,
                "n_reference_edges": len(reference_edges),
                "cell_state_key": cell_state_key,
                "prediction_state_key": cell_state_key,
                "prediction_label_source": "invalid_or_unvalidated",
                "time_key": time_key,
                "ground_truth": ground_truth,
            }
            _enrich_lineage_metrics_compat(metrics, ground_truth, reference_graph_path)
            with open(out_dir / "lineage_metrics.json", "w", encoding="utf-8") as f:
                json.dump(metrics, f, indent=2)
            _mark_run_metadata_lineage_failed(out_dir, reason)
            raise ValueError(f"[eval_lineage] {reason}")

        graph_metrics = compute_graph_sim_metrics(
            predicted_matrix=predicted_matrix,
            reference_matrix=reference_matrix,
            ref_node_ids=ref_node_ids,
        )
        single_step = graph_metrics.get("single_step", {}) or {}
        multi_step  = graph_metrics.get("multi_step",  {}) or {}
        completion_status = "completed"

        # Legacy diagnostic: raw-edge Jaccard from lineage_graph_edges.csv
        # (dense adapter outputs make this near-zero; kept for traceability)
        ref_set_for_diag = {
            (src, tgt) for src, tgt in reference_edges
            if src in predicted_matrix.index and tgt in predicted_matrix.columns
        }
        legacy_topk_jaccard = compute_jaccard(predicted_edges, ref_set_for_diag)

    else:
        legacy_topk_jaccard = None
        if not has_predictions:
            completion_status = (
                "completed_with_empty_predictions: "
                "reference graph was loaded successfully but the method adapter "
                "produced no state transitions (scaffold or runtime failure path). "
                "Check run_metadata.json for the method execution status."
            )
        else:
            completion_status = (
                "completed_with_empty_reference: "
                "method predictions exist but reference matrix is empty after "
                "applying the current edge_confidence_mode filter."
            )

    # ── Correlation baseline ──────────────────────────────────────────────────
    baseline = compute_correlation_baseline(
        adata=adata,
        reference_matrix=reference_matrix,
        reference_edges=reference_edges,
        ref_node_ids=ref_node_ids,
        cell_state_key=cell_state_key,
        time_key=time_key,
        output_dir=str(out_dir),
    )

    # ── Assemble lineage_metrics.json ─────────────────────────────────────────
    metrics = {
        "metric_protocol": "sctimebench_graph_sim",
        "graph_metrics": graph_metrics,
        # Top-level backward-compat aliases (required by eval_dispatch, summary scripts)
        "auroc":               single_step.get("auc_roc"),
        "auprc":               single_step.get("auc_prc"),
        "jaccard_similarity":  single_step.get("jaccard_similarity"),
        "single_step_recovery": single_step.get("recall"),
        "multi_step_recovery":  multi_step.get("recall"),
        # Legacy diagnostic: raw-edge Jaccard from lineage_graph_edges.csv.
        # NOT used for official ranking.  Dense adapter outputs produce values
        # near n_ref/n_total rather than method-discriminating scores.
        "legacy_jaccard_similarity_topk": legacy_topk_jaccard,
        "baseline": baseline,
        "status": completion_status,
        "edge_confidence_mode": edge_confidence_mode,
        "n_reference_edges": len(reference_edges),
        "cell_state_key": cell_state_key,
        "prediction_state_key": cell_state_key,
        "prediction_label_source": PREDICTION_LABEL_SOURCE_SCTIMEBENCH_ZERO_FILL,
        "prediction_label_report": prediction_label_report,
        "time_key": time_key,
        "ground_truth": ground_truth,
    }
    _enrich_lineage_metrics_compat(metrics, ground_truth, reference_graph_path)

    metrics_path = out_dir / "lineage_metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"[eval_lineage] Metrics written to {metrics_path}")
    return metrics
