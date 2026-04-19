# Benchmark

scTimeBench-aligned benchmark for human chemical reprogramming trajectory analysis.

## Framework

This benchmark follows the three core evaluation dimensions defined by scTimeBench:

| Dimension | Status | Eligible methods |
|---|---|---|
| Forecast Accuracy | Inactive (current stage) | Future generative models only |
| Embedding Coherence | Inactive (current stage) | Future generative models only |
| Lineage Fidelity | **Active** | WOT, CellRank2 |

WOT and CellRank2 are evaluated on **Lineage Fidelity only**, consistent with the scTimeBench treatment of OT-based methods that cannot generate projected cells at unseen future time points. No OT-projection workaround is used.

## Directory structure

```
benchmark/
├── README.md                    ← this file
├── configs/
│   ├── benchmark_master.yaml    ← master config (dimensions, methods, outputs)
│   ├── method_capabilities.yaml ← per-method capability flags
│   ├── scenario_observed.yaml   ← scenarios A, B, C (observed time axis)
│   └── scenario_pseudotime.yaml ← scenarios D, E, F (pseudotime axis)
├── adapters/
│   ├── base_adapter.py          ← abstract base class with capability gating
│   ├── wot_adapter.py           ← WOT adapter (Lineage Fidelity only)
│   ├── cellrank2_adapter.py     ← CellRank2 adapter (Lineage Fidelity only)
│   └── future_model_adapter.py  ← template for future generative models
├── evaluation/
│   ├── eval_dispatch.py         ← capability-gated benchmark dispatcher
│   ├── eval_lineage.py          ← Lineage Fidelity evaluator (ACTIVE)
│   ├── eval_forecast.py         ← Forecast Accuracy evaluator (inactive stub)
│   └── eval_embedding.py        ← Embedding Coherence evaluator (inactive stub)
├── docs/
│   └── framework_summary.md     ← framework design notes
├── reports/                     ← aggregated benchmark reports (generated)
└── results/                     ← per-method per-scenario outputs (generated)
```

## Running the benchmark

```bash
python benchmark/evaluation/eval_dispatch.py \
    --method wot \
    --scenario A \
    --adata data/processed/adata_benchmark.h5ad \
    --output-dir benchmark/results/wot/A
```

Replace `--method wot` with `cellrank2` to run CellRank2.

## Benchmark prerequisites before full Lineage Fidelity runs

Per the framework (v2 §9.3), the following must be explicitly frozen before
Lineage Fidelity metrics can be computed:

1. A defined cell-state system (cell-state annotation for `adata_benchmark.h5ad`)
2. A frozen benchmark reference lineage graph
3. A fixed rule for aggregating method predictions to the state level

Until these are provided, `eval_lineage.py` records a `deferred` status in
`lineage_metrics.json` rather than producing invalid results.

## Adding a future generative model

1. Copy `benchmark/adapters/future_model_adapter.py` and rename it.
2. Set `supports_unseen_timepoint_projection = True`.
3. Implement `_run_forecast_accuracy_impl()`, `_run_embedding_coherence_impl()`,
   and `_run_lineage_fidelity_impl()`.
4. Add the new adapter to `METHOD_REGISTRY` in `eval_dispatch.py`.
5. Add the method to `method_capabilities.yaml`.
