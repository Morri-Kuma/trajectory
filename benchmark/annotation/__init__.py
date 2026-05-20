"""Current annotation utilities for the trajectory benchmark.

The active cell-state system is the frozen ``official_silver`` provider layer.
Observed-cell labels are built from curated milestone markers, sample/time
anchoring where available, and trajectory-aware resolution into
``final_milestone_label_coarse``.
"""

__all__ = [
    "build_marker_seed_labels",
    "build_trajectory_aware_labels",
    "build_milestone_providers",
    "validate_milestone_markers",
    "validate_milestone_configs",
]
