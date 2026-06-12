"""KNN-on-PCA surrogate annotator.

A dependency-light stand-in for scANVI used ONLY in test_mode / when
scvi-tools is unavailable. It trains a k-NN classifier on the reference PCA
embedding + labels and predicts query labels with a confidence score equal to
the fraction of neighbours sharing the predicted label. Low-confidence cells are
flagged 'unknown_or_ood'. It is explicitly tagged so it is never mistaken for a
manuscript result.
"""
from __future__ import annotations

from typing import Dict

import numpy as np


def knn_surrogate_annotate(
    ref_pca: np.ndarray,
    ref_labels: np.ndarray,
    query_pca: np.ndarray,
    k: int = 15,
    confidence_threshold: float = 0.5,
) -> Dict[str, np.ndarray]:
    """Return dict with 'label', 'confidence', 'unknown_flag'.

    ref_pca/query_pca must live in the SAME embedding space (e.g. PCA fit on a
    shared gene space, or query projected onto the reference PCA).
    """
    from sklearn.neighbors import NearestNeighbors

    classes = np.unique(ref_labels)
    k = int(min(k, ref_pca.shape[0]))
    nn = NearestNeighbors(n_neighbors=k).fit(ref_pca)
    _, idx = nn.kneighbors(query_pca)
    neigh_labels = ref_labels[idx]  # (n_query, k)

    pred = np.empty(query_pca.shape[0], dtype=object)
    conf = np.zeros(query_pca.shape[0], dtype=float)
    for i in range(query_pca.shape[0]):
        vals, counts = np.unique(neigh_labels[i], return_counts=True)
        j = int(np.argmax(counts))
        pred[i] = vals[j]
        conf[i] = counts[j] / k

    unknown = conf < confidence_threshold
    label = pred.copy()
    label[unknown] = "unknown_or_ood"
    return {
        "label": label,
        "confidence": conf,
        "unknown_flag": unknown,
        "method": np.array(["knn_surrogate"] * query_pca.shape[0]),
    }
