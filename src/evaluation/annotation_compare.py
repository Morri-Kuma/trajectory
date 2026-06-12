"""The six annotation comparisons (see scanvi_annotation_branch.md §5).

All functions take a query AnnData carrying two annotation columns:
    label_a_key  -> marker-silver label (query-space states)
    label_b_key  -> scANVI/surrogate label (reference-space states)
and return pandas objects / dicts that the pipeline writes to disk.
"""
from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# 1. Cell-type composition comparison
# ---------------------------------------------------------------------------
def composition_compare(adata, label_a_key, label_b_key, time_key) -> Dict:
    """Per-timepoint composition under each annotation + total-variation
    distance between the two compositions."""
    obs = adata.obs
    out_rows = []
    tv_per_time = {}
    for t, sub in obs.groupby(time_key, observed=True):
        ca = sub[label_a_key].value_counts(normalize=True)
        cb = sub[label_b_key].value_counts(normalize=True)
        states = sorted(set(ca.index) | set(cb.index))
        a = np.array([ca.get(s, 0.0) for s in states])
        b = np.array([cb.get(s, 0.0) for s in states])
        tv = 0.5 * np.abs(a - b).sum()
        tv_per_time[str(t)] = float(tv)
        for s in states:
            out_rows.append(
                {"time": str(t), "state": s,
                 "frac_marker_silver": float(ca.get(s, 0.0)),
                 "frac_scanvi": float(cb.get(s, 0.0))}
            )
    df = pd.DataFrame(out_rows)
    overall_tv = float(np.mean(list(tv_per_time.values()))) if tv_per_time else float("nan")
    return {"table": df, "tv_per_time": tv_per_time, "overall_tv": overall_tv}


# ---------------------------------------------------------------------------
# 2/3. Agreement / contingency matrix + scalar metrics
# ---------------------------------------------------------------------------
def agreement_matrix(
    adata, label_a_key, label_b_key, correspondence_map: Optional[Dict[str, str]] = None
) -> Dict:
    """Contingency table (A x B), NMI, ARI, and correspondence-aware recall.

    NMI/ARI measure statistical association and do not require matched label
    spaces. Correspondence recall is computed only where the map defines a
    reference->query correspondence.
    """
    from sklearn.metrics import (
        normalized_mutual_info_score,
        adjusted_rand_score,
    )

    a = adata.obs[label_a_key].astype(str).values
    b = adata.obs[label_b_key].astype(str).values

    contingency = pd.crosstab(pd.Series(a, name="marker_silver"),
                              pd.Series(b, name="scanvi"))
    metrics = {
        "nmi": float(normalized_mutual_info_score(a, b)),
        "ari": float(adjusted_rand_score(a, b)),
        "n_cells": int(len(a)),
    }

    if correspondence_map:
        # map scANVI (reference-space) labels into query-space via the map
        mapped_b = np.array([correspondence_map.get(x, "unmapped") for x in b])
        # recall per query state that has at least one reference correspondent
        target_states = set(correspondence_map.values())
        recalls = {}
        for s in sorted(target_states):
            mask = a == s
            n = int(mask.sum())
            if n == 0:
                continue
            recalls[s] = float(np.mean(mapped_b[mask] == s))
        metrics["correspondence_recall_per_state"] = recalls
        metrics["correspondence_recall_macro"] = (
            float(np.mean(list(recalls.values()))) if recalls else float("nan")
        )
        # overall fraction of cells whose mapped scANVI label equals marker label
        defined = mapped_b != "unmapped"
        metrics["correspondence_overall_agreement"] = (
            float(np.mean(mapped_b[defined] == a[defined])) if defined.any() else float("nan")
        )
    return {"contingency": contingency, "metrics": metrics}


# ---------------------------------------------------------------------------
# 4. Marker-gene validation
# ---------------------------------------------------------------------------
def marker_validation(adata, marker_sets: Dict[str, List[str]], label_key: str) -> pd.DataFrame:
    """For each (marker program, assigned state) compute AUROC of the marker
    score separating that state from the rest. Higher AUROC = the annotation's
    state is better aligned with the canonical program.
    """
    from sklearn.metrics import roc_auc_score
    from scipy import sparse

    X = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
    genes = {g: i for i, g in enumerate(adata.var_names)}
    labels = adata.obs[label_key].astype(str).values

    rows = []
    for program, mgenes in marker_sets.items():
        cols = [genes[g] for g in mgenes if g in genes]
        if not cols:
            continue
        score = X[:, cols].mean(axis=1)
        for state in np.unique(labels):
            y = (labels == state).astype(int)
            if y.sum() == 0 or y.sum() == len(y):
                continue
            try:
                auc = float(roc_auc_score(y, score))
            except ValueError:
                auc = float("nan")
            rows.append({"annotation": label_key, "program": program,
                         "state": state, "auroc": auc, "n_state": int(y.sum())})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 5. Pseudotime comparison
# ---------------------------------------------------------------------------
def pseudotime_compare(pt_a: np.ndarray, pt_b: np.ndarray) -> Dict:
    """Spearman + Pearson correlation between two pseudotime vectors."""
    from scipy.stats import spearmanr, pearsonr

    pt_a = np.asarray(pt_a, dtype=float)
    pt_b = np.asarray(pt_b, dtype=float)
    ok = np.isfinite(pt_a) & np.isfinite(pt_b)
    if ok.sum() < 3:
        return {"spearman": float("nan"), "pearson": float("nan"), "n": int(ok.sum())}
    sr = spearmanr(pt_a[ok], pt_b[ok]).statistic
    pr = pearsonr(pt_a[ok], pt_b[ok]).statistic
    return {"spearman": float(sr), "pearson": float(pr), "n": int(ok.sum())}


# ---------------------------------------------------------------------------
# 6. Disagreement cells (biological-interpretation effect)
# ---------------------------------------------------------------------------
def disagreement_cells(
    adata, label_a_key, label_b_key, correspondence_map: Optional[Dict[str, str]] = None
) -> pd.DataFrame:
    """Cells where the two annotations disagree (after correspondence mapping),
    with available QC columns for inspection."""
    a = adata.obs[label_a_key].astype(str).values
    b = adata.obs[label_b_key].astype(str).values
    mapped_b = (np.array([correspondence_map.get(x, "unmapped") for x in b])
                if correspondence_map else b)
    disagree = mapped_b != a
    cols = {"marker_silver": a, "scanvi": b, "scanvi_mapped": mapped_b}
    for qc in ["scanvi_confidence", "pct_counts_mt", "n_genes_by_counts", "total_counts"]:
        if qc in adata.obs:
            cols[qc] = adata.obs[qc].values
    df = pd.DataFrame(cols, index=adata.obs_names)
    return df[disagree].copy()
