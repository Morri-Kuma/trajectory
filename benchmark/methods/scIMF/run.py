"""benchmark/methods/scIMF/run.py — scIMF benchmark runner (contract template).

scIMF (Jiang, Li, et al. 2026, PLOS Comput Biol 22(1):e1013916; "Learning collective
multicellular dynamics with an interacting mean field neural SDE model") is a
Transformer + interacting-mean-field neural-SDE generative model, and scTimeBench's
best overall forecaster.

This file is the integration CONTRACT for scIMF. The thin adapter
(benchmark/adapters/scimf_adapter.py via _generative_adapter.py) calls the five
functions below; they must produce the same standardized outputs as
benchmark/methods/scNODE/run.py (use that file as the structural template):

    prepare_data(train_adata, time_key, train_times)
        -> (train_data, train_tps, train_unique_tps, ...)   # per-timepoint expression arrays
    train_or_load(train_data, train_tps, n_genes, cfg, model_cache) -> model
    run_forecast_accuracy(model, adata_full, time_key, all_unique_tps, heldout_tps,
                          n_sim_cells, output_dir, metric_sample_cells, seed)
        -> writes projected_expression.npy, forecast_metrics.json,
           per_timepoint_forecast_metrics.csv
    run_embedding_coherence(model, adata_full, time_key, all_unique_tps, heldout_tps,
                            n_sim_cells, cell_state_key, output_dir)
        -> writes embedding.npy, projected_embedding.npy, next_timepoint_embedding.npy,
           embedding_metrics.json, projected_cluster_labels.csv
    run_lineage_fidelity(model, adata_train, time_key, train_unique_tps,
                         cell_state_key, output_dir)
        -> writes state_transition_matrix.csv, lineage_graph_edges.csv

TO IMPLEMENT (Shirokane, where the source + GPU live):
  1. Vendor the scIMF source under benchmark/methods/scIMF/scIMF_module/ (from the
     authors' repository or scTimeBench's methods/ submodule, mirroring
     benchmark/methods/scNODE/scNODE_module/).
  2. Wire train_or_load() to train the scIMF SDE, and the model.simulate(...) forward
     pass used by the run_* functions; reuse scNODE/run.py's output-writing,
     forecast-metric, embedding-coherence, and STM->graph helpers verbatim (they are
     method-agnostic once the model can project cells to held-out time points).

Config params block: scimf_params (num training epochs, latent dim, n_sim_cells, seed,
use_cuda, ...). Forecast metrics are recomputed centrally by the unified evaluator, so
the runner only needs to emit projected_expression.npy in gene-expression space.
"""
from __future__ import annotations


class VendoringRequired(NotImplementedError):
    pass


_MSG = (
    "scIMF method source is not vendored yet. Place the scIMF implementation under "
    "benchmark/methods/scIMF/scIMF_module/ and wire benchmark/methods/scIMF/run.py to "
    "it (use benchmark/methods/scNODE/run.py as the structural template). See this "
    "file's module docstring for the 5-function contract."
)


def prepare_data(train_adata, time_key, train_times):  # noqa: D401
    raise VendoringRequired(_MSG)


def train_or_load(train_data, train_tps, n_genes, cfg, model_cache):
    raise VendoringRequired(_MSG)


def run_forecast_accuracy(*args, **kwargs):
    raise VendoringRequired(_MSG)


def run_embedding_coherence(*args, **kwargs):
    raise VendoringRequired(_MSG)


def run_lineage_fidelity(*args, **kwargs):
    raise VendoringRequired(_MSG)
