"""marker_utils.py
==================
Lightweight helper utilities for milestone marker-gene scoring.

Used by build_marker_seed_labels.py (Step 3) and later annotation scripts.
No heavy optional libraries are imported at module level.
Standard library only at import time; numpy/scipy are imported inside
functions where needed (they are always available when anndata is present).

Public API
----------
  load_markers_yaml(path)
      Load and return the parsed milestone_markers.yaml as a dict.

  build_gene_index_map(var_names, var_df=None, gene_symbol_column=None,
                       use_var_index=False)
      Return a dict mapping uppercase gene symbol -> list[int] of column
      indices in the AnnData var axis.

  resolve_milestone_genes(marker_config, dataset_id, milestones,
                          gene_index_map, min_overlap, strict)
      Validate marker overlap and return per-milestone matched-gene indices.

  compute_mean_marker_scores(X, milestone_gene_indices)
      Given an (n_cells, n_genes) array (dense or sparse) and a dict of
      milestone -> list[int] column indices, return a dict of
      milestone -> np.ndarray of shape (n_cells,) with mean expression scores.

  assign_labels(score_matrix, milestone_names, min_top_score, min_score_margin)
      Given a (n_cells, n_milestones) score array, return arrays of
      labels, top scores, and margins.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# YAML loading
# ---------------------------------------------------------------------------

def load_markers_yaml(path: str | Path) -> dict:
    """Load milestone_markers.yaml; raise with a clear message if PyYAML absent."""
    try:
        import yaml
    except ImportError as exc:
        raise ImportError(
            "PyYAML is required for loading milestone_markers.yaml.  "
            "Install it with: pip install pyyaml"
        ) from exc
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Markers YAML not found: {path}")
    with path.open(encoding="ascii") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise ValueError(f"Expected top-level mapping in {path}, got {type(data)}")
    return data


# ---------------------------------------------------------------------------
# Gene-index resolution
# ---------------------------------------------------------------------------

def build_gene_index_map(
    var_names: Any,
    var_df: Any = None,
    gene_symbol_column: str | None = None,
    use_var_index: bool = False,
) -> dict[str, list[int]]:
    """Map uppercase gene symbol -> list of integer column indices.

    Resolution order:
      1. Use var_df[gene_symbol_column] if gene_symbol_column is given and
         the column exists in var_df.
      2. If that fails and use_var_index is True, fall back to var_names.
      3. Raise ValueError if neither source is available.

    Parameters
    ----------
    var_names : array-like
        adata.var_names (pandas Index or list of strings).
    var_df : pandas.DataFrame or None
        adata.var (may be None if only var_names are used).
    gene_symbol_column : str or None
        Column name in var_df to use as gene symbols.
    use_var_index : bool
        Allow fallback to var_names when the column is absent.

    Returns
    -------
    dict[str, list[int]]
        Symbol (uppercase) -> list of matching column indices (0-based).
        Duplicates are included; all matching columns are listed.
    """
    symbol_source: list[str] | None = None

    # Attempt 1: named var column.
    if gene_symbol_column and var_df is not None:
        if gene_symbol_column in var_df.columns:
            symbol_source = list(var_df[gene_symbol_column].astype(str))
        else:
            if not use_var_index:
                raise ValueError(
                    f"Column '{gene_symbol_column}' not found in adata.var "
                    f"(columns: {list(var_df.columns)}).  "
                    "Pass --use-var-index-if-needed to fall back to var_names."
                )
            # Fall through to var_names.

    # Attempt 2: var_names fallback.
    if symbol_source is None:
        if use_var_index or not gene_symbol_column:
            symbol_source = list(str(v) for v in var_names)
        else:
            raise ValueError(
                "Could not resolve gene symbols: no gene_symbol_column provided "
                "and --use-var-index-if-needed was not set."
            )

    # Build the map (support duplicates).
    index_map: dict[str, list[int]] = {}
    for idx, sym in enumerate(symbol_source):
        key = sym.upper().strip()
        if not key:
            continue
        index_map.setdefault(key, []).append(idx)

    return index_map


# ---------------------------------------------------------------------------
# Milestone gene resolution
# ---------------------------------------------------------------------------

def resolve_milestone_genes(
    marker_config: dict,
    dataset_id: str,
    milestones: list[str],
    gene_index_map: dict[str, list[int]],
    min_overlap: int,
    strict: bool,
) -> tuple[dict[str, list[int]], dict[str, dict]]:
    """Resolve marker genes to column indices for each milestone.

    Parameters
    ----------
    marker_config : dict
        Parsed top-level YAML dict.
    dataset_id : str
        Dataset key (e.g. "GSE178325").
    milestones : list[str]
        Milestone names to resolve.
    gene_index_map : dict[str, list[int]]
        From build_gene_index_map.
    min_overlap : int
        Minimum number of matched marker genes required.
    strict : bool
        If True, raise on insufficient overlap.  If False, warn and skip.

    Returns
    -------
    milestone_indices : dict[str, list[int]]
        milestone name -> flat list of matched column indices.
    overlap_summary : dict[str, dict]
        Per-milestone overlap diagnostics (counts, matched genes, etc.).
    """
    ds = marker_config.get(dataset_id, {})
    marker_sets = ds.get("marker_sets", {})

    milestone_indices: dict[str, list[int]] = {}
    overlap_summary: dict[str, dict] = {}

    for ms in milestones:
        if ms not in marker_sets:
            msg = (
                f"Milestone '{ms}' has no entry in marker_sets for {dataset_id}"
            )
            if strict:
                raise ValueError(msg)
            print(f"  WARN: {msg}; skipping.")
            continue

        genes_required: list[str] = [
            g.upper().strip()
            for g in (marker_sets[ms].get("genes") or [])
            if str(g).strip()
        ]

        matched_genes: list[str] = []
        matched_indices: list[int] = []
        missing_genes: list[str] = []

        for gene in genes_required:
            cols = gene_index_map.get(gene)
            if cols:
                matched_genes.append(gene)
                matched_indices.extend(cols)
            else:
                missing_genes.append(gene)

        n_matched = len(matched_genes)
        overlap_summary[ms] = {
            "required_genes": genes_required,
            "matched_genes": matched_genes,
            "missing_genes": missing_genes,
            "n_required": len(genes_required),
            "n_matched": n_matched,
            "n_column_indices": len(matched_indices),
        }

        if n_matched < min_overlap:
            msg = (
                f"Milestone '{ms}' has only {n_matched} matched marker gene(s) "
                f"(need >= {min_overlap}).  "
                f"Matched: {matched_genes}.  Missing: {missing_genes}."
            )
            if strict:
                raise ValueError(msg)
            print(f"  WARN: {msg}  Skipping this milestone.")
            continue

        milestone_indices[ms] = matched_indices

    return milestone_indices, overlap_summary


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

def compute_mean_marker_scores(
    X: Any,
    milestone_gene_indices: dict[str, list[int]],
) -> dict[str, Any]:
    """Compute mean expression of matched marker-gene columns per cell.

    Parameters
    ----------
    X : np.ndarray or scipy.sparse matrix
        Expression matrix of shape (n_cells, n_genes).
    milestone_gene_indices : dict[str, list[int]]
        milestone name -> list of column indices to average.

    Returns
    -------
    dict[str, np.ndarray]
        milestone name -> 1-D array of shape (n_cells,) with mean scores.
    """
    import numpy as np

    try:
        import scipy.sparse as sp
        is_sparse = sp.issparse(X)
    except ImportError:
        is_sparse = False

    scores: dict[str, Any] = {}
    for ms, indices in milestone_gene_indices.items():
        if not indices:
            continue
        # Slice only the relevant columns; convert sparse slice to dense.
        if is_sparse:
            col_slice = X[:, indices].toarray()  # (n_cells, n_matched)
        else:
            col_slice = np.asarray(X[:, indices], dtype=np.float64)
        scores[ms] = np.nanmean(col_slice, axis=1)  # (n_cells,)

    return scores


# ---------------------------------------------------------------------------
# Label assignment
# ---------------------------------------------------------------------------

def assign_labels(
    score_dict: dict[str, Any],
    milestone_names: list[str],
    min_top_score: float,
    min_score_margin: float,
) -> tuple[Any, Any, Any]:
    """Assign per-cell milestone labels from score arrays.

    Parameters
    ----------
    score_dict : dict[str, np.ndarray]
        milestone name -> 1-D score array (n_cells,).
        May not contain all milestone_names if some were skipped.
    milestone_names : list[str]
        Ordered list of milestone names (defines column order).
    min_top_score : float
        Cells with top score < this threshold are labelled 'ambiguous'.
    min_score_margin : float
        Cells where top - 2nd score < this threshold are labelled 'ambiguous'.

    Returns
    -------
    labels : np.ndarray[str]   shape (n_cells,)
    top_scores : np.ndarray[float]   shape (n_cells,)
    margins : np.ndarray[float]   shape (n_cells,)
    """
    import numpy as np

    # Determine n_cells from any available score array.
    n_cells = 0
    for arr in score_dict.values():
        n_cells = len(arr)
        break
    if n_cells == 0:
        return (
            np.array([], dtype=object),
            np.array([], dtype=np.float64),
            np.array([], dtype=np.float64),
        )

    # Build score matrix for milestones that have scores.
    scored_names = [m for m in milestone_names if m in score_dict]
    if not scored_names:
        labels = np.full(n_cells, "ambiguous", dtype=object)
        top_scores = np.full(n_cells, np.nan)
        margins = np.full(n_cells, np.nan)
        return labels, top_scores, margins

    score_mat = np.column_stack(
        [score_dict[m] for m in scored_names]
    )  # (n_cells, n_scored)

    top_idx = np.argmax(score_mat, axis=1)          # (n_cells,)
    top_scores = score_mat[np.arange(n_cells), top_idx]

    # Margin: top - second-best (or just top if only one milestone scored).
    if score_mat.shape[1] >= 2:
        # Mask the top column, take max of remainder.
        tmp = score_mat.copy()
        tmp[np.arange(n_cells), top_idx] = -np.inf
        second_scores = np.max(tmp, axis=1)
        margins = top_scores - second_scores
    else:
        margins = np.full(n_cells, np.inf)

    # Build label array.
    name_arr = np.array(scored_names, dtype=object)
    labels = name_arr[top_idx]

    # Apply thresholds.
    ambiguous_mask = (
        np.isnan(top_scores)
        | (top_scores < min_top_score)
        | (margins < min_score_margin)
    )
    labels[ambiguous_mask] = "ambiguous"
    top_scores = np.where(ambiguous_mask, np.nan, top_scores)
    margins = np.where(ambiguous_mask, np.nan, margins)

    return labels, top_scores, margins
