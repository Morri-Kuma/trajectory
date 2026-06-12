"""benchmark/methods/PISDE/run.py — PI-SDE benchmark runner (contract template).

PI-SDE (Jiang & Wan 2024, Bioinformatics 40(Suppl 2):ii120-ii127; "A physics-informed
neural SDE network for learning cellular dynamics from time-series scRNA-seq data") is a
physics-informed neural SDE, projection-capable -> all three benchmark dimensions.

This is the integration CONTRACT. The adapter (benchmark/adapters/pisde_adapter.py)
calls the same five functions as benchmark/methods/scNODE/run.py — see that file as the
structural template and the scIMF run.py docstring for the full I/O contract:

    prepare_data, train_or_load, run_forecast_accuracy, run_embedding_coherence,
    run_lineage_fidelity

TO IMPLEMENT (Shirokane): vendor the PI-SDE source under
benchmark/methods/PISDE/PISDE_module/ and wire train_or_load() + the SDE forward
simulation; reuse scNODE/run.py's output/metric/STM helpers (method-agnostic once the
model can project cells to held-out time points). Config params block: pisde_params.
Note: PI-SDE was CPU-bound in scTimeBench (GPU/OOM) — budget runtime accordingly.
"""
from __future__ import annotations


class VendoringRequired(NotImplementedError):
    pass


_MSG = (
    "PI-SDE method source is not vendored yet. Place it under "
    "benchmark/methods/PISDE/PISDE_module/ and wire benchmark/methods/PISDE/run.py "
    "(use benchmark/methods/scNODE/run.py as the structural template)."
)


def prepare_data(train_adata, time_key, train_times):
    raise VendoringRequired(_MSG)


def train_or_load(train_data, train_tps, n_genes, cfg, model_cache):
    raise VendoringRequired(_MSG)


def run_forecast_accuracy(*args, **kwargs):
    raise VendoringRequired(_MSG)


def run_embedding_coherence(*args, **kwargs):
    raise VendoringRequired(_MSG)


def run_lineage_fidelity(*args, **kwargs):
    raise VendoringRequired(_MSG)
