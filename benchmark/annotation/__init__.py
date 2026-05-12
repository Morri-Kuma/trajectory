"""benchmark.annotation
=======================
Milestone-based annotation module for scTimeBench framework v2.

This package provides a three-mode annotation layer that assigns milestone
labels to cells in both the reference (observed) and projected datasets.
The three label modes are:

  consensus
      state_key: consensus_milestone_label
      Primary reporting mode.  Labels are agreed upon by at least two of
      the sensitivity modes below.

  embedding_based
      state_key: milestone_embedding_label
      Sensitivity mode.  Intended for scGPT or similar foundation-model
      embeddings followed by a lightweight milestone classifier.

  classifier_based
      state_key: milestone_classifier_label
      Sensitivity mode.  Intended for label-transfer tools such as
      CellTypist, scANVI, TOSICA, SingleR, or scArches.

Script inventory
----------------
  build_marker_seed_labels        -- derive seed labels from marker genes (Step 3)
  build_embedding_milestone_model -- train scGPT embedding milestone classifier
  build_hvg_pca_embedding_milestone_model -- train HVG PCA embedding classifier
  build_classifier_milestone_model -- run / wrap external label-transfer tool
  build_consensus_milestone_labels -- merge modes into consensus_milestone_label
  annotate_projected_cells        -- apply trained models to projected outputs
  validate_milestone_markers      -- validate milestone_markers.yaml (Step 2)
  build_milestone_providers       -- build ground-truth providers (Step 4)

Config inventory
----------------
  milestone_markers.yaml          -- dataset-specific milestone marker sets and
                                     reference graph edges (Step 2)

Status: scaffold (Steps 1-4).  No heavy optional dependencies are imported here.
"""

__all__ = [
    "build_marker_seed_labels",
    "build_embedding_milestone_model",
    "build_hvg_pca_embedding_milestone_model",
    "build_classifier_milestone_model",
    "build_consensus_milestone_labels",
    "annotate_projected_cells",
    "validate_milestone_markers",
    "build_milestone_providers",
]
