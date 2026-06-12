"""Annotation branch: marker-silver (existing, read from obs) vs scANVI/scArches
reference-mapping (new). A KNN-on-PCA surrogate stands in for scANVI when
scvi-tools/torch are unavailable, so the comparison pipeline runs anywhere.
A geometric OOD score (ood.py) recalibrates the miscalibrated scANVI confidence.
"""
from .knn_surrogate import knn_surrogate_annotate
from .dispatch import annotate_query
from .ood import ood_score_and_flag, knn_reference_distance, ood_gate

__all__ = ["knn_surrogate_annotate", "annotate_query",
           "ood_score_and_flag", "knn_reference_distance", "ood_gate"]
