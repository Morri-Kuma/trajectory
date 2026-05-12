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
from scipy.sparse import issparse
from scipy.stats import rankdata


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
        # Group-sum by state within this chunk.
        chunk_states = cell_states[start:end]
        # Precompute per-state masks via unique for speed.
        for s_str, i in state_to_idx.items():
            mask = chunk_states == s_str
            if not mask.any():
                continue
            sums[i] += chunk_dense[mask].sum(axis=0)
            counts[i] += int(mask.sum())

    # Safe per-state mean.
    means = np.zeros_like(sums)
    nonzero = counts > 0
    means[nonzero] = sums[nonzero] / counts[nonzero, None]

    return pd.DataFrame(means, index=list(states))


def _build_baseline_edges(corr_df: pd.DataFrame,
                           n_top: int) -> Tuple[pd.DataFrame, set]:
    """
    Binarize the correlation "transition" matrix into a sparse edge list by
    selecting the top-N off-diagonal entries (matching the single-step rule).

    This gives the baseline a predicted-edges set comparable to what the
    method adapters emit, so the edge-based Jaccard is not trivially 1.0
    for the full dense matrix.
    """
    flat = corr_df.values.copy()
    # Mask diagonal so self-correlation (1.0) is not treated as an edge.
    np.fill_diagonal(flat, -np.inf)
    n = flat.size
    if n_top <= 0 or n_top > n:
        n_top = max(1, min(n_top, n))
    threshold = np.sort(flat.flatten())[-n_top]
    edge_rows: list = []
    edges: set = set()
    states = list(corr_df.index)
    for i, src in enumerate(states):
        for j, tgt in enumerate(corr_df.columns):
            v = flat[i, j]
            if v >= threshold and np.isfinite(v):
                edge_rows.append({"source_state": src,
                                  "target_state": tgt,
                                  "weight": float(corr_df.iloc[i, j])})
                edges.add((src, tgt))
    edges_df = pd.DataFrame(
        edge_rows if edge_rows else [],
        columns=["source_state", "target_state", "weight"],
    )
    return edges_df, edges


def compute_state_mean_pearson_baseline(adata,
                                         reference_matrix: pd.DataFrame,
                                         reference_edges: set,
                                         cell_state_key: Optional[str] = None,
                                         output_dir: Optional[str] = None) -> dict:
    """
    Correlation-based baseline for Lineage Fidelity (per v2 §9.6).

    Legacy diagnostic only: uses pairwise Pearson
    correlation between per-state mean gene expression vectors as a naive
    predictor of lineage connectivity, then evaluates the resulting
    (symmetric) score matrix against the reference lineage using the same
    metrics as the main evaluator.

    Notes on interpretation
    -----------------------
    The baseline is intentionally naïve and **undirected**: correlation is
    symmetric, so the baseline cannot distinguish A→B from B→A.  That is the
    point — a non-trivial method must beat a direction-agnostic similarity
    signal to demonstrate real lineage inference value.  Diagonal entries
    are masked out so self-correlation (=1.0) does not dominate the
    predicted-edges list.

    Traceability
    ------------
    If ``output_dir`` is provided, the baseline writes two side files so the
    baseline is fully reproducible from disk:
      - ``baseline_state_transition_matrix.csv`` (the correlation matrix)
      - ``baseline_lineage_graph_edges.csv`` (top-N correlation edges)
    """
    none_result = {
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

    # Align state list with the reference matrix so predicted and reference
    # matrices have identical index/columns when metrics are computed.
    states = list(reference_matrix.index)

    # Per-state mean expression.
    try:
        means_df = _state_mean_expression(adata, cell_state_key, states)
    except Exception as exc:  # pragma: no cover — defensive
        return {**none_result,
                "note": f"Correlation baseline failed during mean-expression computation: {exc!s}"}

    # States present in adata (non-zero mean vector).
    row_norms = np.linalg.norm(means_df.values, axis=1)
    present_mask = row_norms > 0

    if present_mask.sum() < 2:
        return {**none_result,
                "note": (
                    "Correlation baseline skipped: fewer than 2 reference states "
                    f"are represented in adata under {cell_state_key!r}."
                )}

    # Pearson correlation across states. Clip NaN rows (zero-variance) to 0.
    corr = np.corrcoef(means_df.values)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    # Shift from [-1, 1] → [0, 1] so it behaves like a probability score for
    # AUROC/AUPRC (the sklearn implementations care about ranking, so the
    # shift is metric-preserving, but it keeps the score non-negative for
    # downstream consistency with the WOT/CR2 STM conventions).
    score = (corr + 1.0) / 2.0
    corr_df = pd.DataFrame(score, index=states, columns=states, dtype=float)

    # Metrics — aligned with the main evaluator so numbers are directly
    # comparable across method vs baseline rows in the summary.
    y_true = reference_matrix.values.flatten().astype(int)
    y_score = corr_df.values.flatten().astype(float)
    metrics = {
        "auroc": compute_auroc(y_true, y_score),
        "auprc": compute_auprc(y_true, y_score),
    }

    # Edge-based metrics require a sparse predicted-edges set.
    n_ref = int(reference_matrix.values.sum())
    edges_df, edges_set = _build_baseline_edges(corr_df, n_top=n_ref)
    metrics["jaccard_similarity"] = compute_jaccard(edges_set, reference_edges)
    metrics["single_step_recovery"] = compute_single_step_recovery(
        corr_df, reference_matrix
    )
    metrics["multi_step_recovery"] = compute_multi_step_recovery(
        corr_df, reference_matrix
    )

    # Provenance fields so the baseline row is traceable in the summary.
    metrics["method"] = "state_mean_pearson_baseline"
    metrics["baseline_style"] = "legacy_state_mean_pearson"
    metrics["cell_state_key"] = cell_state_key
    metrics["n_states_used"] = int(present_mask.sum())
    metrics["n_predicted_edges_topk"] = len(edges_set)
    metrics["note"] = (
        "Per-state mean-expression Pearson correlation (symmetric), "
        f"top-{n_ref} off-diagonal entries as predicted lineage edges."
    )

    # Optional on-disk traceability artifacts.
    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        corr_df.to_csv(out / "baseline_state_transition_matrix.csv")
        edges_df.to_csv(out / "baseline_lineage_graph_edges.csv", index=False)

    return metrics


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
    """
    none_result = {
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

    y_true = reference_matrix.values.flatten().astype(int)
    y_score = baseline_df.values.flatten().astype(float)
    metrics = {
        "auroc": compute_auroc(y_true, y_score),
        "auprc": compute_auprc(y_true, y_score),
    }

    edges_df, edges_set, graph_threshold = _build_pr_threshold_edges(
        baseline_df, reference_edges
    )
    metrics["jaccard_similarity"] = compute_jaccard(edges_set, reference_edges)
    metrics["single_step_recovery"] = compute_single_step_recovery(
        baseline_df, reference_matrix
    )
    metrics["multi_step_recovery"] = compute_multi_step_recovery(
        baseline_df, reference_matrix
    )

    metrics["method"] = "correlation_baseline"
    metrics["baseline_style"] = "sctimebench_correlation"
    metrics["correlation_method"] = "spearmanr"
    metrics["averaging_method"] = "maximum"
    metrics["cell_state_key"] = cell_state_key
    metrics["time_key"] = resolved_time_key
    metrics["n_states_used"] = int(nonzero_rows.sum())
    metrics["n_timepoint_pairs_used"] = int(n_pairs_used)
    metrics["n_source_cell_votes"] = int(n_source_cell_votes)
    metrics["n_predicted_edges_thresholded"] = len(edges_set)
    metrics["graph_threshold"] = graph_threshold
    metrics["note"] = (
        "scTimeBench-style cell-level Spearman correlation baseline with "
        "maximum target-state aggregation and row-normalized source-state votes."
    )

    if output_dir is not None:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        baseline_df.to_csv(out / "baseline_state_transition_matrix.csv")
        edges_df.to_csv(out / "baseline_lineage_graph_edges.csv", index=False)

    return metrics


# ------------------------------------------------------------------
# Main evaluation function
# ------------------------------------------------------------------

def _topk_predicted_edges(predicted_matrix: pd.DataFrame, k: int) -> set:
    """
    Build a top-k predicted-edges set from a (possibly dense) state-transition
    matrix by selecting the k off-diagonal entries with the largest weight.

    This mirrors the rule used by ``compute_single_step_recovery`` so that the
    edge-based Jaccard is computed on a *sparse, comparable* edge list rather
    than on whatever the adapter happened to dump into ``lineage_graph_edges.csv``
    (WOT and CellRank2 currently emit fully dense 14×14 = 196 edges, which
    forces the raw-edge Jaccard to collapse to ``n_ref / 196`` for any method
    whose STM is dense — see benchmark/docs/jaccard_singlestep_diagnosis.md).
    Self-loops are excluded.
    """
    if predicted_matrix.empty or k <= 0:
        return set()
    vals = predicted_matrix.values.copy().astype(float)
    np.fill_diagonal(vals, -np.inf)
    flat = vals.flatten()
    finite = flat[np.isfinite(flat)]
    if finite.size == 0:
        return set()
    k = min(k, finite.size)
    threshold = np.sort(finite)[-k]
    edges: set = set()
    states_idx = list(predicted_matrix.index)
    states_col = list(predicted_matrix.columns)
    for i, src in enumerate(states_idx):
        for j, tgt in enumerate(states_col):
            if np.isfinite(vals[i, j]) and vals[i, j] >= threshold:
                edges.add((src, tgt))
    return edges


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
    cell_state_key : str, optional
        Name of the obs column that holds the cell-state label used by the
        method (e.g. ``"scgpt_pseudostate_provisional"``). Required by the
        correlation baseline; if absent the baseline reports a clear "skipped"
        status instead of returning silently empty values.
    time_key : str, optional
        Name of the obs column that holds the temporal coordinate used by the
        scTimeBench-style correlation baseline. If absent, common names such as
        ``timepoint``, ``abs_day``, and ``time_label`` are tried in order.
    ground_truth : dict, optional
        Provider metadata for the annotation/reference graph used by this run.
        Stored in lineage_metrics.json for downstream sensitivity analysis.

    Returns
    -------
    dict with keys: auroc, auprc, jaccard_similarity, jaccard_similarity_topk,
                    single_step_recovery, multi_step_recovery, baseline,
                    status, edge_confidence_mode, n_reference_edges
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ground_truth = ground_truth or {}
    _mark_run_metadata_ground_truth(out_dir, ground_truth)

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
            "jaccard_similarity_topk": None,
            "single_step_recovery": None,
            "multi_step_recovery": None,
            "baseline": None,
            "status": status,
            "edge_confidence_mode": edge_confidence_mode,
            "n_reference_edges": None,
            "cell_state_key": cell_state_key,
            "time_key": time_key,
            "ground_truth": ground_truth,
        }
        _enrich_lineage_metrics_compat(metrics, ground_truth, reference_graph_path)
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
        all_states = sorted(set(predicted_matrix.index) | set(reference_matrix.index))
        pred_aligned = predicted_matrix.reindex(
            index=all_states, columns=all_states, fill_value=0.0
        )
        ref_aligned = reference_matrix.reindex(
            index=all_states, columns=all_states, fill_value=0
        )
        y_true = ref_aligned.values.flatten().astype(int)
        y_score = pred_aligned.values.flatten().astype(float)
        if not np.isfinite(y_score).all():
            reason = _nonfinite_matrix_message(
                pred_aligned, "state_transition_matrix.csv"
            )
            metrics.update({
                "auroc": None,
                "auprc": None,
                "jaccard_similarity": None,
                "jaccard_similarity_topk": None,
                "single_step_recovery": None,
                "multi_step_recovery": None,
                "baseline": None,
                "status": f"failed_nonfinite_predictions: {reason}",
                "edge_confidence_mode": edge_confidence_mode,
                "n_reference_edges": len(reference_edges),
                "cell_state_key": cell_state_key,
                "time_key": time_key,
                "ground_truth": ground_truth,
            })
            _enrich_lineage_metrics_compat(metrics, ground_truth, reference_graph_path)
            metrics_path = out_dir / "lineage_metrics.json"
            with open(metrics_path, "w", encoding="utf-8") as f:
                json.dump(metrics, f, indent=2)
            _mark_run_metadata_lineage_failed(out_dir, reason)
            raise ValueError(f"[eval_lineage] {reason}")

        metrics["auroc"] = compute_auroc(y_true, y_score)
        metrics["auprc"] = compute_auprc(y_true, y_score)
        # Raw edge Jaccard — treats every nonzero entry in
        # lineage_graph_edges.csv as a predicted edge.  This is the
        # scTimeBench-style set-Jaccard.  For current adapters it collapses
        # to n_ref / 196 because WOT/CR2 emit dense 14×14 edge lists.
        metrics["jaccard_similarity"] = compute_jaccard(predicted_edges, reference_edges)
        # Top-k edge Jaccard — restricts the predicted edge set to the
        # top n_ref off-diagonal entries of the STM (same rule as
        # single_step_recovery). This is the method-discriminating
        # Jaccard used for ranking; see jaccard_singlestep_diagnosis.md.
        ref_set_for_topk = {
            (src, tgt) for src, tgt in reference_edges
            if src in pred_aligned.index and tgt in pred_aligned.columns
        }
        topk_edges = _topk_predicted_edges(pred_aligned, len(ref_set_for_topk))
        metrics["jaccard_similarity_topk"] = compute_jaccard(
            topk_edges, ref_set_for_topk
        )
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
            "auroc": None, "auprc": None,
            "jaccard_similarity": None, "jaccard_similarity_topk": None,
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

    # Correlation baseline (per v2 §9.6).
    # Uses the same scTimeBench Correlation policy as spearman_max.yaml:
    # adjacent-time cell-level Spearman, maximum target-state score, and
    # row-normalized source-state votes. It is evaluated against the SAME
    # reference graph with the SAME metric functions as method outputs.
    metrics["baseline"] = compute_correlation_baseline(
        adata=adata,
        reference_matrix=reference_matrix,
        reference_edges=reference_edges,
        cell_state_key=cell_state_key,
        time_key=time_key,
        output_dir=str(out_dir),
    )
    metrics["status"] = completion_status
    metrics["edge_confidence_mode"] = edge_confidence_mode
    metrics["n_reference_edges"] = len(reference_edges)
    metrics["cell_state_key"] = cell_state_key
    metrics["time_key"] = time_key
    metrics["ground_truth"] = ground_truth
    _enrich_lineage_metrics_compat(metrics, ground_truth, reference_graph_path)

    metrics_path = out_dir / "lineage_metrics.json"
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"[eval_lineage] Metrics written to {metrics_path}")
    return metrics
