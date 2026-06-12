"""pisde_adapter.py — adapter for PI-SDE (Jiang & Wan 2024, Bioinformatics).

PI-SDE is a physics-informed neural SDE for learning cellular dynamics from
time-series scRNA-seq. Projection-capable -> all three benchmark dimensions.
Method internals live in benchmark/methods/PISDE/run.py.
"""
from __future__ import annotations

from ._generative_adapter import GenerativeProjectionAdapter


class PISDEAdapter(GenerativeProjectionAdapter):
    _method_id = "pisde"
    _params_key = "pisde_params"
    _run_module = "benchmark.methods.PISDE.run"
