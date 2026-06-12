"""
eval_representation_temporal_signal.py
======================================

Representation-space temporal & state diagnostics (Work-Plan Step 6).

These are computed on OBSERVED cells in a representation space (X_rep or a raw
scFM embedding) *before* any trajectory model is trained. They answer whether
the representation itself retains temporal order and biological-state structure.

Diagnostics
-----------
    temporal_signal_ratio (TVR)        variance explained by timepoint centroids
    adjacent_time_centroid_distance    mean distance between consecutive-time centroids
    non_adjacent_time_centroid_distance mean distance between non-adjacent centroids
    time_prediction_macro_f1           CV macro-F1 of classifying timepoint from rep
    time_prediction_r2                 CV R^2 of regressing continuous time from rep
    state_centroid_separation          mean pairwise distance between state centroids
    state_silhouette_score             silhouette of state labels in rep space
    forward_mass_fraction              transition mass moving forward in time
    backward_mass_fraction             transition mass moving backward in time
    forbidden_edge_mass                mass on edges forbidden by a lineage graph
                                       (only if a reference graph is supplied)

All functions are numpy/sklearn based (no torch) and operate on plain arrays,
so they are unit-testable on small synthetic inputs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import numpy as np


# ---------------------------------------------------------------------------
# Centroid helpers
# ---------------------------------------------------------------------------

def _group_centroids(X: np.ndarray, labels: np.ndarray):
    """Return (unique_labels, centroids) ordered by unique_labels."""
    uniq = np.unique(labels)
    cents = np.vstack([X[labels == u].mean(axis=0) for u in uniq])
    return uniq, cents


def temporal_signal_ratio(X: np.ndarray, times: np.ndarray) -> float:
    """
    Fraction of total representation variance explained by timepoint centroids
    (between-group SS / total SS), a.k.a. time-variance ratio (TVR) in [0, 1].

    Higher => more of the representation's spread is organized by time.
    """
    X = np.asarray(X, dtype=np.float64)
    times = np.asarray(times, dtype=float)
    grand = X.mean(axis=0)
    ss_total = float(np.sum((X - grand) ** 2))
    if ss_total == 0:
        return 0.0
    ss_between = 0.0
    for t in np.unique(times):
        mask = times == t
        n_t = int(mask.sum())
        cen = X[mask].mean(axis=0)
        ss_between += n_t * float(np.sum((cen - grand) ** 2))
    return float(ss_between / ss_total)


def adjacent_vs_nonadjacent_centroid_distance(X: np.ndarray, times: np.ndarray):
    """Mean centroid distance between adjacent vs non-adjacent timepoints."""
    X = np.asarray(X, dtype=np.float64)
    times = np.asarray(times, dtype=float)
    uniq = np.sort(np.unique(times))
    cents = {t: X[times == t].mean(axis=0) for t in uniq}
    adj, non_adj = [], []
    for i, ti in enumerate(uniq):
        for j, tj in enumerate(uniq):
            if j <= i:
                continue
            dist = float(np.linalg.norm(cents[ti] - cents[tj]))
            if j == i + 1:
                adj.append(dist)
            else:
                non_adj.append(dist)
    return (
        float(np.mean(adj)) if adj else None,
        float(np.mean(non_adj)) if non_adj else None,
    )


def time_prediction_scores(X: np.ndarray, times: np.ndarray, seed: int = 0):
    """
    Cross-validated predictability of time from the representation.

    Returns (macro_f1, r2). macro_f1 treats each unique timepoint as a class;
    r2 regresses the continuous time value. Uses simple, dependency-light
    estimators with k-fold CV (k <= 5, bounded by class counts).
    """
    from sklearn.linear_model import LogisticRegression, Ridge
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.metrics import f1_score, r2_score

    X = np.asarray(X, dtype=np.float64)
    times = np.asarray(times, dtype=float)
    classes, y = np.unique(times, return_inverse=True)

    macro_f1 = None
    if len(classes) >= 2:
        min_count = np.min(np.bincount(y))
        n_splits = int(max(2, min(5, min_count)))
        if min_count >= 2:
            skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
            try:
                y_pred = cross_val_predict(
                    LogisticRegression(max_iter=1000),
                    X, y, cv=skf,
                )
                macro_f1 = float(f1_score(y, y_pred, average="macro"))
            except Exception:
                macro_f1 = None

    r2 = None
    if len(times) >= 4 and len(classes) >= 2:
        n_splits = int(max(2, min(5, len(times) // 2)))
        try:
            t_pred = cross_val_predict(Ridge(alpha=1.0), X, times, cv=n_splits)
            r2 = float(r2_score(times, t_pred))
        except Exception:
            r2 = None
    return macro_f1, r2


def state_centroid_separation(X: np.ndarray, states: np.ndarray) -> Optional[float]:
    """Mean pairwise distance between state centroids (higher => more separated)."""
    X = np.asarray(X, dtype=np.float64)
    uniq, cents = _group_centroids(X, np.asarray(states))
    if len(uniq) < 2:
        return None
    dists = []
    for i in range(len(uniq)):
        for j in range(i + 1, len(uniq)):
            dists.append(float(np.linalg.norm(cents[i] - cents[j])))
    return float(np.mean(dists))


def state_silhouette(X: np.ndarray, states: np.ndarray) -> Optional[float]:
    """Silhouette score of state labels in representation space."""
    from sklearn.metrics import silhouette_score

    states = np.asarray(states)
    if len(np.unique(states)) < 2 or len(states) <= len(np.unique(states)):
        return None
    try:
        return float(silhouette_score(np.asarray(X, dtype=np.float64), states))
    except Exception:
        return None


def transition_mass_fractions(
    X: np.ndarray,
    states: np.ndarray,
    times: np.ndarray,
    *,
    state_time_order: Optional[dict] = None,
    forbidden_edges: Optional[set] = None,
    temperature: float = 1.0,
):
    """
    Estimate directionality of representation-space transitions.

    For each pair of consecutive timepoints, assign each source cell's
    *predicted next position* (here: the source cell itself, since this is a
    model-free diagnostic on observed cells, mapped to the nearest next-time
    state centroid) and accumulate soft mass over (source_state -> target_state)
    edges. Mass is split into:
      - forward  : target state's typical time > source state's typical time
      - backward : target state's typical time < source state's typical time
      - forbidden: edge not allowed by ``forbidden_edges`` (if provided)

    ``state_time_order`` maps state -> representative (e.g. median) time; if
    None it is computed from the data. Returns dict with the three fractions.

    This is a representation diagnostic, not a trained-model transition matrix;
    it tests whether nearest-centroid structure across time moves forward.
    """
    X = np.asarray(X, dtype=np.float64)
    states = np.asarray(states)
    times = np.asarray(times, dtype=float)

    if state_time_order is None:
        state_time_order = {
            s: float(np.median(times[states == s])) for s in np.unique(states)
        }

    uniq_t = np.sort(np.unique(times))
    fwd = bwd = forbidden = total = 0.0

    for ti, tj in zip(uniq_t[:-1], uniq_t[1:]):
        src_mask = times == ti
        tgt_mask = times == tj
        if not src_mask.any() or not tgt_mask.any():
            continue
        tgt_states, tgt_cents = _group_centroids(X[tgt_mask], states[tgt_mask])
        src_states = states[src_mask]
        src_X = X[src_mask]
        # soft assignment of each source cell to next-time state centroids
        d2 = _sq_dists(src_X, tgt_cents)
        w = _softmax(-d2 / max(temperature, 1e-8), axis=1)
        for r in range(src_X.shape[0]):
            s_state = src_states[r]
            s_time = state_time_order.get(s_state, float(ti))
            for c, t_state in enumerate(tgt_states):
                mass = float(w[r, c])
                total += mass
                t_time = state_time_order.get(t_state, float(tj))
                if forbidden_edges is not None and (s_state, t_state) in forbidden_edges:
                    forbidden += mass
                if t_time > s_time:
                    fwd += mass
                elif t_time < s_time:
                    bwd += mass

    if total == 0:
        return {"forward_mass_fraction": None, "backward_mass_fraction": None,
                "forbidden_edge_mass": None}
    return {
        "forward_mass_fraction": float(fwd / total),
        "backward_mass_fraction": float(bwd / total),
        "forbidden_edge_mass": (float(forbidden / total)
                                if forbidden_edges is not None else None),
    }


def _sq_dists(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    a2 = np.sum(a * a, axis=1)[:, None]
    b2 = np.sum(b * b, axis=1)[None, :]
    return np.maximum(a2 + b2 - 2.0 * a @ b.T, 0.0)


def _softmax(z: np.ndarray, axis: int) -> np.ndarray:
    z = z - np.max(z, axis=axis, keepdims=True)
    e = np.exp(z)
    return e / np.sum(e, axis=axis, keepdims=True)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def compute_all_diagnostics(
    X: np.ndarray,
    times: np.ndarray,
    states: Optional[np.ndarray] = None,
    *,
    forbidden_edges: Optional[set] = None,
    seed: int = 0,
) -> dict:
    """Compute the full temporal/state diagnostic suite for one representation."""
    adj, non_adj = adjacent_vs_nonadjacent_centroid_distance(X, times)
    macro_f1, r2 = time_prediction_scores(X, times, seed=seed)
    out = {
        "temporal_signal_ratio": temporal_signal_ratio(X, times),
        "adjacent_time_centroid_distance": adj,
        "non_adjacent_time_centroid_distance": non_adj,
        "time_prediction_macro_f1": macro_f1,
        "time_prediction_r2": r2,
        "n_cells": int(np.asarray(X).shape[0]),
        "n_timepoints": int(len(np.unique(np.asarray(times)))),
    }
    if states is not None:
        out["state_centroid_separation"] = state_centroid_separation(X, states)
        out["state_silhouette_score"] = state_silhouette(X, states)
        out.update(transition_mass_fractions(
            X, states, times, forbidden_edges=forbidden_edges
        ))
        out["n_states"] = int(len(np.unique(np.asarray(states))))
    return out


def _forbidden_edges_from_reference_graph(graph_path: str) -> set:
    """
    Build the set of FORBIDDEN (source_state, target_state) edges from a
    reference lineage graph JSON: any ordered pair of distinct states whose
    edge is NOT present (in either direction per the graph's directedness) is
    treated as forbidden. Returns a set of tuples.
    """
    with open(graph_path, encoding="utf-8") as f:
        graph = json.load(f)
    nodes = set()
    allowed = set()
    edges = graph.get("edges", graph.get("links", []))
    for e in edges:
        s = str(e.get("source") or e.get("from"))
        t = str(e.get("target") or e.get("to"))
        nodes.update([s, t])
        allowed.add((s, t))
    for n in graph.get("nodes", []):
        nodes.add(str(n.get("id", n)) if isinstance(n, dict) else str(n))
    forbidden = set()
    for s in nodes:
        for t in nodes:
            if s != t and (s, t) not in allowed:
                forbidden.add((s, t))
    return forbidden


def run_temporal_signal_evaluation(
    X: np.ndarray,
    times: np.ndarray,
    output_dir: str,
    states: Optional[np.ndarray] = None,
    *,
    reference_graph_path: Optional[str] = None,
    representation_id: Optional[str] = None,
    representation_metadata: Optional[dict] = None,
    seed: int = 0,
) -> dict:
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    forbidden = None
    if reference_graph_path and states is not None:
        forbidden = _forbidden_edges_from_reference_graph(reference_graph_path)

    diag = compute_all_diagnostics(
        X, times, states, forbidden_edges=forbidden, seed=seed
    )
    diag["representation_id"] = representation_id
    diag["metric_space"] = "representation"
    diag["representation_metadata"] = representation_metadata or {}
    diag["reference_graph_used"] = bool(forbidden is not None)

    with open(out_dir / "representation_temporal_signal.json", "w", encoding="utf-8") as f:
        json.dump(diag, f, indent=2)
    return diag


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Representation temporal/state diagnostics.")
    p.add_argument("--representation-npy", required=True,
                   help="(n_cells, d) observed representation matrix.")
    p.add_argument("--times-npy", required=True,
                   help="(n_cells,) per-cell timepoints.")
    p.add_argument("--states-npy", default=None,
                   help="(n_cells,) per-cell state labels (optional).")
    p.add_argument("--reference-graph", default=None,
                   help="Reference lineage graph JSON for forbidden_edge_mass.")
    p.add_argument("--representation-id", default=None)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--seed", type=int, default=0)
    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    X = np.load(args.representation_npy, allow_pickle=True)
    times = np.load(args.times_npy, allow_pickle=True)
    states = (np.load(args.states_npy, allow_pickle=True)
              if args.states_npy else None)
    run_temporal_signal_evaluation(
        X, times, args.output_dir, states=states,
        reference_graph_path=args.reference_graph,
        representation_id=args.representation_id,
        seed=args.seed,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
