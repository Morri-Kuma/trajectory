"""benchmark/methods/Squidiff/run.py — Squidiff benchmark runner (contract template).

Squidiff (He, Zhu, et al. 2026, Nature Methods 23(1):65-77; "Squidiff: predicting
cellular development and responses to perturbations using a diffusion model") is a
conditional diffusion model, projection-capable -> all three benchmark dimensions.

This is the integration CONTRACT. The adapter (benchmark/adapters/squiddiff_adapter.py)
calls the same five functions as benchmark/methods/scNODE/run.py (see that file as the
structural template; see scIMF run.py for the full I/O contract):

    prepare_data, train_or_load, run_forecast_accuracy, run_embedding_coherence,
    run_lineage_fidelity

TO IMPLEMENT (Shirokane): vendor the Squidiff source under
benchmark/methods/Squidiff/Squidiff_module/ and wire train_or_load() + the diffusion
sampling forward pass. IMPORTANT: Squidiff has no native gene-expression decoder — fit
a PCA on the observed training data, train/sample in PCA space, and inverse-transform
projected PCs back to expression before writing projected_expression.npy (the
MIOFlow/PRESCIENT/PI-SDE convention used throughout scTimeBench). Reuse scNODE/run.py's
output/metric/STM helpers. Config params block: squiddiff_params.
"""
from __future__ import annotations


class VendoringRequired(NotImplementedError):
    pass


_MSG = (
    "Squidiff method source is not vendored yet. Place it under "
    "benchmark/methods/Squidiff/Squidiff_module/ and wire "
    "benchmark/methods/Squidiff/run.py (use benchmark/methods/scNODE/run.py as the "
    "structural template; train/sample in PCA space, inverse-transform to expression)."
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
