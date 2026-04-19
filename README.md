# scTimeBench-Aligned Benchmark for Human Chemical iPSC Reprogramming

**Author**: Kuma — Graduate School of Frontier Sciences, The University of Tokyo
**Dataset**: GSE230659 (Liuyang et al. 2023 *Cell Stem Cell*, 75,194 cells, 15 timepoints)
**Framework**: Experimental Framework v2 (see `experimental framework v2.md`)

---

## Project overview

This repository benchmarks trajectory inference methods for human chemical iPSC reprogramming using a framework aligned with **scTimeBench**. The benchmark evaluates methods across three core dimensions: Forecast Accuracy, Embedding Coherence, and Lineage Fidelity.

The current first-stage comparison focuses on **WOT** and **CellRank2**. Because these methods do not generate projected cells at unseen future time points, they are evaluated on **Lineage Fidelity only**, consistent with the scTimeBench treatment of OT-based methods.

---

## Three core benchmark dimensions

| Dimension | Status | Eligible methods |
|---|---|---|
| Forecast Accuracy | Inactive (current stage) | Future generative models only |
| Embedding Coherence | Inactive (current stage) | Future generative models only |
| Lineage Fidelity | **Active** | WOT, CellRank2 |

No OT-projection workaround is used to force WOT or CellRank2 into unsupported dimensions.

---

## Benchmark scenarios

| ID | Time axis | Type |
|---|---|---|
| A | Observed time | Interpolation |
| B | Observed time | Extrapolation |
| C | Observed time | Interpolation + extrapolation |
| D | Pseudotime | Interpolation |
| E | Pseudotime | Extrapolation |
| F | Pseudotime | Interpolation + extrapolation |

For WOT and CellRank2, scenarios are entered for Lineage Fidelity only.

---

## Repository structure

```
trajectory/
├── experimental framework v2.md     ← source of truth for all benchmark logic
├── README.md                        ← this file
├── data/                            ← all dataset assets (do not modify)
│   ├── gse230659(human/             ← primary dataset (raw scRNA-seq, 15 samples)
│   ├── gse178325_human/             ← external validation dataset
│   ├── gse280956(human/             ← additional dataset
│   ├── GSE298212(human/             ← additional dataset
│   ├── FCR_iPSC(mouse/              ← mouse dataset
│   └── processed/                   ← processed h5ad files
└── benchmark/                       ← benchmark framework
    ├── README.md                    ← benchmark-specific documentation
    ├── configs/
    │   ├── benchmark_master.yaml    ← master config (dimensions, methods, outputs)
    │   ├── method_capabilities.yaml ← per-method capability flags
    │   ├── scenario_observed.yaml   ← scenarios A, B, C
    │   └── scenario_pseudotime.yaml ← scenarios D, E, F
    ├── adapters/
    │   ├── base_adapter.py          ← abstract base class with capability gating
    │   ├── wot_adapter.py           ← WOT (Lineage Fidelity only)
    │   ├── cellrank2_adapter.py     ← CellRank2 (Lineage Fidelity only)
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

---

## Running the benchmark

```bash
# Run WOT through the dispatcher (Lineage Fidelity only)
python benchmark/evaluation/eval_dispatch.py \
    --method wot \
    --scenario A \
    --adata data/processed/adata_benchmark.h5ad \
    --output-dir benchmark/results/wot/A

# Run CellRank2 through the dispatcher (Lineage Fidelity only)
python benchmark/evaluation/eval_dispatch.py \
    --method cellrank2 \
    --scenario A \
    --adata data/processed/adata_benchmark.h5ad \
    --output-dir benchmark/results/cellrank2/A
```

The dispatcher reads capability flags from each adapter and routes the method
only to its eligible evaluators. For WOT and CellRank2, only `eval_lineage.py`
is called. Forecast Accuracy and Embedding Coherence are skipped with an
explanatory log message.

---

## Prerequisites before Lineage Fidelity runs

Per the framework (v2 §9.3), the following must be explicitly frozen before
Lineage Fidelity metrics can be computed:

1. A defined cell-state system (cell-state annotations for `adata_benchmark.h5ad`)
2. A frozen benchmark reference lineage graph
3. A fixed rule for aggregating method predictions to the state level

Until these are provided, `eval_lineage.py` records a `deferred` status
in `lineage_metrics.json` rather than producing invalid results.

---

## Dataset summary

| Property | Value |
|---|---|
| Primary GEO accession | GSE230659 |
| Reference | Liuyang et al. 2023 *Cell Stem Cell* |
| Cells (post-QC) | 75,194 |
| Genes (HVG) | 2,000 |
| Timepoints | 15 (day 0.5 – 30) |
| Terminal state | hCiPSC (day 30) |

---

## Adding a future generative model

When a model that can generate projected cells at unseen time points is added:

1. Copy `benchmark/adapters/future_model_adapter.py` and rename it.
2. Set `supports_unseen_timepoint_projection = True`.
3. Implement `_run_forecast_accuracy_impl()`, `_run_embedding_coherence_impl()`,
   and `_run_lineage_fidelity_impl()`.
4. Add the new adapter to `METHOD_REGISTRY` in `eval_dispatch.py`.
5. Add the method to `benchmark/configs/method_capabilities.yaml`.

At that point, `eval_dispatch.py` will automatically activate Forecast Accuracy
and Embedding Coherence for the new model.

---

## Framework reference

All benchmark logic is governed by `experimental framework v2.md` in the
repository root. That file is the authoritative design document and takes
precedence over all other documentation.
