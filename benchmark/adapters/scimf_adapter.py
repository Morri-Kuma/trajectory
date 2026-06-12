"""scimf_adapter.py — adapter for scIMF (Jiang et al. 2026, PLOS Comput Biol).

scIMF is an interacting mean-field neural SDE with a Transformer encoder; in
scTimeBench it is the best overall forecaster. Projection-capable -> eligible for
Forecast Accuracy, Embedding Coherence, and Lineage Fidelity. Orchestration is the
shared scNODE-style flow; the method internals live in benchmark/methods/scIMF/run.py.
"""
from __future__ import annotations

from ._generative_adapter import GenerativeProjectionAdapter


class ScIMFAdapter(GenerativeProjectionAdapter):
    _method_id = "scimf"
    _params_key = "scimf_params"
    _run_module = "benchmark.methods.scIMF.run"
