"""squiddiff_adapter.py — adapter for Squidiff (He et al. 2026, Nature Methods).

Squidiff is a conditional diffusion model predicting cellular development and
perturbation responses. Projection-capable -> all three benchmark dimensions. As it
has no native gene-expression decoder, run.py trains/simulates in PCA space and
inverse-transforms projected PCs back to expression (the MIOFlow/PRESCIENT convention
used in scTimeBench). Method internals live in benchmark/methods/Squidiff/run.py.
"""
from __future__ import annotations

from ._generative_adapter import GenerativeProjectionAdapter


class SquidiffAdapter(GenerativeProjectionAdapter):
    _method_id = "squiddiff"
    _params_key = "squiddiff_params"
    _run_module = "benchmark.methods.Squidiff.run"
