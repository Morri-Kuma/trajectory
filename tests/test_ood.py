import numpy as np
from src.annotation import ood_score_and_flag

def test_far_query_flagged_ood():
    rng = np.random.default_rng(0)
    ref = rng.normal(0, 1, (200, 5))
    far = rng.normal(20, 1, (50, 5))     # clearly off the reference manifold
    res = ood_score_and_flag(ref, far, k=15, quantile=0.95)
    assert res["frac_ood"] > 0.95
    assert res["ood_score"].min() > 1.0

def test_inlier_query_not_flagged():
    rng = np.random.default_rng(1)
    ref = rng.normal(0, 1, (300, 5))
    inl = rng.normal(0, 1, (60, 5))      # same distribution as reference
    res = ood_score_and_flag(ref, inl, k=15, quantile=0.95)
    assert res["frac_ood"] < 0.20

def test_ood_gate_rejects_far_query_with_enrichment():
    from src.annotation import ood_gate
    rng = np.random.default_rng(2)
    ref = rng.normal(0, 1, (300, 5))
    far = rng.normal(20, 1, (80, 5))
    states = np.array(["a"] * 40 + ["b"] * 40)
    g = ood_gate(ref, far, k=15, quantile=0.95, state_labels=states)
    assert g["gate_decision"] == "reject"
    assert g["enrichment_vs_ref"] > 3.0
    assert set(g["ood_by_state"]) == {"a", "b"}
    assert all(0.0 <= v <= 1.0 for v in g["ood_by_state"].values())

def test_ood_gate_accepts_inlier_query():
    from src.annotation import ood_gate
    rng = np.random.default_rng(3)
    ref = rng.normal(0, 1, (300, 5))
    inl = rng.normal(0, 1, (80, 5))
    g = ood_gate(ref, inl, k=15, quantile=0.95)
    assert g["gate_decision"] == "accept"
    assert g["enrichment_vs_ref"] < 3.0
