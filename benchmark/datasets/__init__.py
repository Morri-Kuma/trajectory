# benchmark/datasets/__init__.py
#
# Legacy package kept for existing imports. New scTimeBench-style dataset
# registries live under benchmark.shared.dataset.registry.

from benchmark.shared.dataset.registry import GSE178325Dataset, GSE230659Dataset

__all__ = ["GSE178325Dataset", "GSE230659Dataset"]
