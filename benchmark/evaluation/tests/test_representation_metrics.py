"""
Unit tests for the representation-space evaluators on small synthetic arrays.

Covers:
  - rep_ forecast metric functions (numpy backend): self-distance ~0, symmetry,
    separation increases with shift, Hausdorff not normalized by dims.
  - temporal-signal diagnostics: TVR bounds, adjacency distances, time/state
    predictability, transition mass directionality, forbidden-edge mass.
  - forbidden-edge parsing from a reference graph JSON.
"""

import json

import numpy as np
import pytest

from benchmark.evaluation.eval_representation_forecast import (
    METRIC_KEYS,
    compute_rep_metrics,
    numpy_energy_distance,
    numpy_gaussian_mmd,
    numpy_hausdorff,
    numpy_sinkhorn_wasserstein,
)
from benchmark.evaluation.eval_representation_temporal_signal import (
    _forbidden_edges_from_reference_graph,
    adjacent_vs_nonadjacent_centroid_distance,
    compute_all_diagnostics,
    state_centroid_separation,
    state_silhouette,
    temporal_signal_ratio,
    time_prediction_scores,
    transition_mass_fractions,
)


def _cloud(n, d, center, seed, scale=0.1):
    rng = np.random.default_rng(seed)
    return rng.normal(loc=center, scale=scale, size=(n, d))


# ---------------------------------------------------------------------------
# Forecast metrics
# ---------------------------------------------------------------------------

def test_gaussian_mmd_self_near_zero_and_positive_for_shift():
    X = _cloud(40, 4, 0.0, seed=1)
    Y = _cloud(40, 4, 0.0, seed=2)
    Z = _cloud(40, 4, 5.0, seed=3)
    assert abs(numpy_gaussian_mmd(X, X)) < 1e-6
    assert numpy_gaussian_mmd(X, Z) > numpy_gaussian_mmd(X, Y)


def test_energy_distance_self_near_zero_and_monotone():
    X = _cloud(40, 4, 0.0, seed=1)
    near = _cloud(40, 4, 0.3, seed=2)
    far = _cloud(40, 4, 4.0, seed=3)
    assert abs(numpy_energy_distance(X, X)) < 1e-6
    assert numpy_energy_distance(X, far) > numpy_energy_distance(X, near)


def test_hausdorff_not_normalized_by_dims():
    rng = np.random.default_rng(7)
    X3 = rng.random((10, 3))
    Y3 = rng.random((10, 3))
    # duplicate columns -> 6 dims; Hausdorff must NOT divide by dims
    X6 = np.concatenate([X3, X3], axis=1)
    Y6 = np.concatenate([Y3, Y3], axis=1)
    v3 = numpy_hausdorff(X3, Y3)
    v6 = numpy_hausdorff(X6, Y6)
    # distance grows because dims double (sqrt(2) factor), but it's not divided
    assert v6 > v3


def test_sinkhorn_wasserstein_self_near_zero_and_increasing():
    X = _cloud(30, 3, 0.0, seed=1)
    near = _cloud(30, 3, 0.5, seed=2)
    far = _cloud(30, 3, 3.0, seed=3)
    w_self = numpy_sinkhorn_wasserstein(X, X, blur=0.5)
    w_near = numpy_sinkhorn_wasserstein(X, near, blur=0.5)
    w_far = numpy_sinkhorn_wasserstein(X, far, blur=0.5)
    assert w_self < 1e-3
    assert w_far > w_near > w_self


def test_compute_rep_metrics_numpy_backend_keys():
    X = _cloud(20, 5, 0.0, seed=1)
    Y = _cloud(20, 5, 1.0, seed=2)
    m = compute_rep_metrics(X, Y, backend="numpy")
    for k in METRIC_KEYS:
        assert k in m and m[k] is not None
    assert m["backend"] == "numpy_exact"


# ---------------------------------------------------------------------------
# Temporal-signal diagnostics
# ---------------------------------------------------------------------------

def test_temporal_signal_ratio_bounds():
    # perfectly time-separated clouds with tiny within-variance -> TVR ~ 1
    X = np.vstack([_cloud(20, 2, 0.0, 1, scale=1e-3),
                   _cloud(20, 2, 10.0, 2, scale=1e-3)])
    times = np.array([0.0] * 20 + [1.0] * 20)
    tvr = temporal_signal_ratio(X, times)
    assert 0.9 < tvr <= 1.0

    # no time structure: both timepoints drawn from same distribution -> low TVR
    Xn = _cloud(40, 2, 0.0, 5, scale=1.0)
    times_n = np.array([0.0] * 20 + [1.0] * 20)
    assert temporal_signal_ratio(Xn, times_n) < 0.3


def test_adjacent_vs_nonadjacent_distance():
    # centroids at t=0,1,2 along a line; adjacent dist < non-adjacent dist
    X = np.vstack([
        _cloud(10, 1, 0.0, 1, scale=1e-3),
        _cloud(10, 1, 1.0, 2, scale=1e-3),
        _cloud(10, 1, 2.0, 3, scale=1e-3),
    ])
    times = np.array([0.0] * 10 + [1.0] * 10 + [2.0] * 10)
    adj, non_adj = adjacent_vs_nonadjacent_centroid_distance(X, times)
    assert adj is not None and non_adj is not None
    assert non_adj > adj


def test_time_prediction_scores_separable():
    X = np.vstack([_cloud(25, 3, 0.0, 1, scale=0.2),
                   _cloud(25, 3, 8.0, 2, scale=0.2)])
    times = np.array([0.0] * 25 + [1.0] * 25)
    f1, r2 = time_prediction_scores(X, times, seed=0)
    assert f1 is not None and f1 > 0.9


def test_state_separation_and_silhouette():
    X = np.vstack([_cloud(20, 2, 0.0, 1, scale=0.1),
                   _cloud(20, 2, 5.0, 2, scale=0.1)])
    states = np.array(["a"] * 20 + ["b"] * 20)
    sep = state_centroid_separation(X, states)
    sil = state_silhouette(X, states)
    assert sep > 4.0
    assert sil is not None and sil > 0.5


def test_transition_mass_forward_dominant_when_states_progress_with_time():
    # state a at t0, state b at t1, b's typical time later -> forward dominant
    X = np.vstack([_cloud(15, 2, 0.0, 1, scale=0.1),
                   _cloud(15, 2, 5.0, 2, scale=0.1)])
    states = np.array(["a"] * 15 + ["b"] * 15)
    times = np.array([0.0] * 15 + [1.0] * 15)
    res = transition_mass_fractions(X, states, times)
    assert res["forward_mass_fraction"] is not None
    assert res["forward_mass_fraction"] > res["backward_mass_fraction"]


def test_forbidden_edge_mass_uses_supplied_set():
    X = np.vstack([_cloud(15, 2, 0.0, 1, scale=0.1),
                   _cloud(15, 2, 5.0, 2, scale=0.1)])
    states = np.array(["a"] * 15 + ["b"] * 15)
    times = np.array([0.0] * 15 + [1.0] * 15)
    # forbid the a->b edge: all transition mass should be flagged forbidden
    res = transition_mass_fractions(X, states, times, forbidden_edges={("a", "b")})
    assert res["forbidden_edge_mass"] is not None
    assert res["forbidden_edge_mass"] > 0.9


def test_compute_all_diagnostics_keys():
    X = np.vstack([_cloud(20, 3, 0.0, 1, scale=0.2),
                   _cloud(20, 3, 5.0, 2, scale=0.2)])
    times = np.array([0.0] * 20 + [1.0] * 20)
    states = np.array(["a"] * 20 + ["b"] * 20)
    out = compute_all_diagnostics(X, times, states)
    for k in ("temporal_signal_ratio", "adjacent_time_centroid_distance",
              "time_prediction_macro_f1", "state_silhouette_score",
              "forward_mass_fraction", "backward_mass_fraction"):
        assert k in out


def test_forbidden_edges_from_reference_graph(tmp_path):
    graph = {
        "nodes": [{"id": "A"}, {"id": "B"}, {"id": "C"}],
        "edges": [{"source": "A", "target": "B"}, {"source": "B", "target": "C"}],
    }
    p = tmp_path / "graph.json"
    p.write_text(json.dumps(graph))
    forbidden = _forbidden_edges_from_reference_graph(str(p))
    # allowed: A->B, B->C. Everything else among {A,B,C} (excluding self) forbidden.
    assert ("A", "B") not in forbidden
    assert ("B", "C") not in forbidden
    assert ("B", "A") in forbidden
    assert ("A", "C") in forbidden
    assert ("A", "A") not in forbidden  # self-loops excluded
