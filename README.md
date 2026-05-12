# scTimeBench-Aligned iPSC Trajectory Benchmark

**Author**: Kuma, Graduate School of Frontier Sciences, The University of Tokyo  
**Primary domain**: human chemical iPSC reprogramming  
**Framework**: Experimental Framework v2 (see `experimental framework v2.md`)

---

## Project Overview

This repository implements a benchmark for trajectory inference and forecasting methods in human chemical iPSC reprogramming. It is aligned with the evaluation philosophy of **scTimeBench**, but adapts the benchmark target, state system, reference graph, and reporting workflow to iPSC reprogramming.

The project is not a direct copy of scTimeBench. Instead, it uses scTimeBench's core evaluation dimensions as the base and adds iPSC-specific benchmark assets:

- human chemical reprogramming datasets, including GSE230659 and GSE178325-derived validation inputs;
- scGPT-derived pseudo-state annotations as the benchmark state abstraction;
- frozen state-level reference lineage graphs;
- method capability gating, so each method is evaluated only on supported tasks;
- per-method, per-scenario outputs for forecast accuracy, embedding coherence, and lineage fidelity.

---

## Current Status

The current benchmark has moved beyond the initial WOT/CellRank2 lineage-only stage.

| Component | Status |
|---|---|
| GSE230659 scGPT-v1 benchmark | Completed for formal A/B/C result sets where available |
| GSE178325_0618 HVG2000 benchmark | Completed for MIOFlow, PRESCIENT, and scNODE across scenarios A/B/C |
| Forecast Accuracy | Active for generative / projected-cell methods |
| Embedding Coherence | Active for methods with projected embeddings |
| Lineage Fidelity | Active for all methods with state-level transition output |
| WOT and CellRank2 | Evaluated on Lineage Fidelity only |

Recent summarized outputs are available under:

- `benchmark/results/summary_gse178325_0618_forecast_embedding.csv`
- `benchmark/results/summary_lineage_metrics_preferred.csv`
- `benchmark/reports/formal_benchmark_summary.csv`

---

## Evaluation Dimensions

| Dimension | Status | Eligible methods |
|---|---|---|
| Forecast Accuracy | Active | Methods that generate projected expression at held-out or future time points, such as MIOFlow, PRESCIENT, and scNODE |
| Embedding Coherence | Active | Methods that provide projected embeddings or projected cells that can be embedded |
| Lineage Fidelity | Active | Methods that provide or can be converted to state-level transition predictions |

Capability gating is intentional. WOT and CellRank2 are not forced into forecast or embedding metrics because they do not natively generate unseen-timepoint projected cells in the same sense as generative trajectory models.

---

## Benchmark Scenarios

| ID | Time axis | Type |
|---|---|---|
| A | Observed time | Interpolation / observed-time evaluation |
| B | Observed time | Extrapolation / held-out future-time evaluation |
| C | Observed time | Mixed interpolation and extrapolation |
| D | Pseudotime | Interpolation |
| E | Pseudotime | Extrapolation |
| F | Pseudotime | Mixed interpolation and extrapolation |

The currently reported iPSC formal results focus on observed-time scenarios A/B/C. Pseudotime scenarios remain part of the framework design and can be activated when the corresponding inputs and reference state system are frozen.

---

## Repository Structure

```text
trajectory/
|-- experimental framework v2.md      # benchmark design reference
|-- README.md                         # project-level overview
|-- data/                             # raw and processed dataset assets
|-- benchmark/
|   |-- README.md                     # benchmark-specific usage notes
|   |-- configs/                      # method, scenario, and run configs
|   |-- adapters/                     # method adapters and capability gates
|   |-- evaluation/                   # forecast, embedding, and lineage evaluators
|   |-- ground_truth/                 # reference state labels and lineage graphs
|   |-- inputs/                       # benchmark-ready h5ad inputs
|   |-- methods/                      # vendored or wrapped method implementations
|   |-- reports/                      # aggregated reports
|   |-- results/                      # per-method per-scenario outputs
|   `-- shared/                       # dataset and utility code
|-- scripts/                          # data preparation and reporting scripts
|-- logs/                             # run logs
`-- results/                          # auxiliary analysis outputs
```

---

## Methods

The current benchmark includes:

| Method | Forecast Accuracy | Embedding Coherence | Lineage Fidelity |
|---|---:|---:|---:|
| MIOFlow | yes | yes | yes |
| PRESCIENT | yes | yes | yes |
| scNODE | yes | yes | yes |
| WOT | no | no | yes |
| CellRank2 | no | no | yes |

Method capability flags are configured in `benchmark/configs/method_capabilities.yaml`.

---

## Running the Benchmark

Representative observed-time configs are stored in `benchmark/configs/`. For example:

```bash
python benchmark/evaluation/eval_dispatch.py \
    --config benchmark/configs/scnode_gse178325_observed_0618_hvg2000_A.yaml
```

For method-specific batch runs, use the root-level helper scripts such as:

```bash
bash run_scnode_gse178325_0618_hvg_array.sh
bash run_mioflow_gse178325_0618_hvg_array.sh
bash run_prescient_gse178325_0618_hvg_array.sh
```

The dispatcher reads method capabilities and skips unsupported metrics with explicit metadata rather than producing forced or invalid outputs.

---

## Ground Truth and Lineage Fidelity

Lineage Fidelity depends on three frozen assets:

1. a benchmark cell-state system, currently based on scGPT pseudo-state annotations;
2. a reference state-level lineage graph;
3. a fixed rule for aggregating method predictions to state-level transitions.

These assets are stored under `benchmark/ground_truth/` and referenced from the benchmark configs. Per-run lineage outputs include `state_transition_matrix.csv`, `lineage_graph_edges.csv`, `lineage_metrics.json`, and diagnostic metadata where available.

---

## Relationship to scTimeBench

This project is best understood as a domain-specific extension of scTimeBench ideas:

- scTimeBench provides the general benchmark vocabulary and metric families.
- trajectory specializes the benchmark to iPSC reprogramming.
- trajectory uses scGPT-derived pseudo-states and iPSC-specific reference graphs.
- trajectory keeps method capability gating explicit so transition-only methods and generative methods are compared on appropriate evaluation surfaces.

The goal is therefore not strict code-structure parity with scTimeBench, but a scTimeBench-aligned benchmark that answers a narrower biological question.

---

## Framework Reference

`experimental framework v2.md` is the main design document for benchmark logic. If implementation notes and README text diverge, treat the framework document and current configs/results as the more specific source of truth.
