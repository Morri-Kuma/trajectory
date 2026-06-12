"""
benchmark.representations
=========================

Supplementary *representation-dynamics* module for the scTimeBench-aligned
benchmark. See ``benchmark/docs/scfm_representation_dynamics_plan.md`` for the
full work plan and scientific rationale.

This package is intentionally limited to four representation arms:

    rep_hvg_pca50
    rep_geneformer_cls_pca50
    rep_scgpt_cls_pca50
    rep_scfoundation_pca50

It does NOT implement GeneCompass, UCE, hybrid representations, or random
projection controls.

Reference-first guarantee
-------------------------
Every single-cell foundation-model (scFM) extractor in this package is a *thin
adapter* over the official package / checkpoint for that model. The adapters do
not re-implement tokenization, vocabulary matching, model loading, pooling, or
embedding-layer logic. When an official dependency or checkpoint is not present
locally, the adapter raises a clear, actionable error instead of substituting an
invented approximation. The authoritative references are recorded in
``benchmark.representations.references.MODEL_REFERENCES``.
"""

from .references import MODEL_REFERENCES  # noqa: F401

__all__ = ["MODEL_REFERENCES"]
