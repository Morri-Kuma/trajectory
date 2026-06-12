import numpy as np
from benchmark.evaluation.lineage_robustness import (
    edge_confusion, pr_curve_points, bootstrap_edge_metric,
)

NODES = ["A", "B", "C", "D"]
REF = [("A", "B"), ("B", "C"), ("C", "D")]  # a 3-edge chain, like the silver graphs


def test_edge_confusion_perfect():
    c = edge_confusion(REF, REF, NODES)
    assert c["tp"] == 3 and c["fp"] == 0 and c["fn"] == 0
    assert c["n_pairs"] == 12  # 4*3 ordered off-diagonal pairs
    assert c["jaccard"] == 1.0 and c["f1"] == 1.0


def test_edge_confusion_partial():
    pred = [("A", "B"), ("A", "C")]  # 1 correct, 1 wrong, 2 missed
    c = edge_confusion(pred, REF, NODES)
    assert c["tp"] == 1 and c["fp"] == 1 and c["fn"] == 2
    assert c["tn"] == 12 - 1 - 1 - 2
    assert 0.0 < c["jaccard"] < 1.0


def test_pr_curve_and_auprc_separates():
    # candidate edges scored; true edges get high scores
    y_true = np.array([1, 1, 1, 0, 0, 0, 0, 0])
    y_score = np.array([0.9, 0.8, 0.7, 0.2, 0.1, 0.15, 0.05, 0.3])
    pr = pr_curve_points(y_true, y_score)
    assert pr["auprc"] > 0.8
    assert pr["precision"].ndim == 1 and pr["recall"].ndim == 1


def test_bootstrap_ci_brackets_point():
    rng = np.random.default_rng(0)
    y_true = np.array([1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
    y_score = y_true * 0.6 + rng.normal(0, 0.1, y_true.size)
    out = bootstrap_edge_metric(y_true, y_score, n_boot=500, seed=1)
    assert out["n_boot_valid"] > 100
    assert out["lo"] <= out["mean"] <= out["hi"]
    assert 0.0 <= out["lo"] <= 1.0 and 0.0 <= out["hi"] <= 1.0
