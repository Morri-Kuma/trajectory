"""Ground-truth provider registry for benchmark state systems."""

from .loader import GroundTruthSpec, load_ground_truth

__all__ = ["GroundTruthSpec", "load_ground_truth"]
