import numpy as np
from src.annotation import knn_surrogate_annotate

def test_separable_labels():
    rng = np.random.default_rng(0)
    ref = np.vstack([rng.normal(0, 0.2, (50, 2)), rng.normal(8, 0.2, (50, 2))])
    labels = np.array(["A"] * 50 + ["B"] * 50)
    q = np.array([[0.0, 0.0], [8.0, 8.0]])
    res = knn_surrogate_annotate(ref, labels, q, k=5, confidence_threshold=0.5)
    assert list(res["label"]) == ["A", "B"]
    assert res["confidence"].min() >= 0.5

def test_low_confidence_becomes_ood():
    ref = np.vstack([np.zeros((5, 2)), np.full((5, 2), 0.01)])
    lab = np.array(["A"] * 5 + ["B"] * 5)
    q = np.array([[0.005, 0.0]])
    res = knn_surrogate_annotate(ref, lab, q, k=10, confidence_threshold=0.9)
    assert bool(res["unknown_flag"][0]) is True
    assert res["label"][0] == "unknown_or_ood"
