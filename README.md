# scTimeBench-Aligned iPSC Trajectory Benchmark

**Author**: Kuma, Graduate School of Frontier Sciences, The University of Tokyo  
**Primary domain**: human iPSC reprogramming trajectories  
**Framework**: official silver trajectory benchmark

---

## Project Overview

This repository implements a benchmark for trajectory inference and forecasting
methods in human iPSC reprogramming. It is aligned with the evaluation
philosophy of scTimeBench, while adapting the datasets, state systems,
reference graphs, and reporting workflow to iPSC reprogramming.

The project is not a direct copy of scTimeBench. It uses scTimeBench's core
evaluation dimensions as the base and adds domain-specific benchmark assets:

- frozen silver-standard milestone providers for GSE178325 and GSE230659;
- an OSKM-specific GSE242424 silver provider built from author-cluster-matched
  annotations;
- method capability gating, so each method is evaluated only on supported
  tasks;
- per-method, per-scenario outputs for forecast accuracy, embedding coherence,
  and lineage fidelity;
- explicit separation of current reportable results from pilot, smoke, and
  archived/debug outputs.

## Current Status

The benchmark uses only the current official silver label providers for
reportable runs.

| Component | Status |
|---|---|
| GSE178325 marker-FM transition silver benchmark | Formal A/B/C result sets available for MIOFlow, PRESCIENT, and scNODE |
| GSE230659 marker-FM transition silver benchmark | Formal A/B/C result sets available for MIOFlow, PRESCIENT, and scNODE; WOT and CellRank2 have lineage-only scenario A checks |
| GSE242424 OSKM silver benchmark | Formal A/B/C result sets available for MIOFlow, PRESCIENT, and scNODE on the 59,187-cell author-cluster-matched subset |
| Forecast Accuracy | Active for generative/projected-cell methods |
| Embedding Coherence | Active for methods with projected embeddings or projected cells |
| Lineage Fidelity | Active for methods with state-level transition output |
Current report entry points:

- `benchmark/reports/official_silver/official_silver_model_rankings.md`
- `benchmark/reports/official_silver/official_silver_combined_method_summary_scnode_mioflow_prescient.csv`
- `benchmark/reports/official_silver/official_silver_embedding_method_summary.csv`
- `benchmark/reports/official_silver/official_silver_lineage_method_summary.csv`
- `benchmark/reports/gse242424_oskm_silver/gse242424_oskm_silver_report.md`
- `benchmark/results/result_manifest.yaml`

## Evaluation Dimensions

| Dimension | Status | Eligible methods |
|---|---|---|
| Forecast Accuracy | Active | Methods that generate projected expression at held-out or future time points, such as MIOFlow, PRESCIENT, and scNODE |
| Embedding Coherence | Active | Methods that provide projected embeddings or projected cells that can be embedded |
| Lineage Fidelity | Active | Methods that provide or can be converted to state-level transition predictions |

Capability gating is intentional. WOT and CellRank2 are not forced into
forecast or embedding metrics because they do not natively generate
unseen-timepoint projected cells in the same sense as generative trajectory
models.

## Benchmark Scenarios

| ID | Time axis | Type |
|---|---|---|
| A | Observed time | Interpolation / observed-time evaluation |
| B | Observed time | Extrapolation / held-out future-time evaluation |
| C | Observed time | Mixed interpolation and extrapolation |
| D | Pseudotime | Interpolation |
| E | Pseudotime | Extrapolation |
| F | Pseudotime | Mixed interpolation and extrapolation |

Current formal reports focus on observed-time scenarios A/B/C. Pseudotime
scenarios remain part of the framework design and can be activated when the
corresponding inputs and reference state systems are frozen.

## Repository Structure

```text
trajectory/
|-- README.md                         # project-level overview
|-- data/                             # raw and processed dataset assets
|-- benchmark/
|   |-- README.md                     # benchmark-specific usage notes
|   |-- annotation/                   # milestone/provider annotation utilities
|   |-- configs/                      # method, scenario, runtime, and run configs
|   |-- adapters/                     # method adapters and capability gates
|   |-- evaluation/                   # forecast, embedding, and lineage evaluators
|   |-- ground_truth/                 # frozen state labels and reference graphs
|   |-- inputs/                       # benchmark-ready h5ad inputs
|   |-- methods/                      # vendored or wrapped method implementations
|   |-- reports/                      # aggregated reports
|   |-- results/                      # per-method per-scenario outputs
|   `-- shared/                       # dataset and utility code
|-- scripts/                          # data preparation and reporting scripts
|-- logs/                             # run logs
`-- reference/                        # paper and external-reference notes
```

## Methods

The current benchmark includes:

| Method | Forecast Accuracy | Embedding Coherence | Lineage Fidelity |
|---|---:|---:|---:|
| MIOFlow | yes | yes | yes |
| PRESCIENT | yes | yes | yes |
| scNODE | yes | yes | yes |
| WOT | no | no | yes |
| CellRank2 | no | no | yes |

Method capability flags are configured in
`benchmark/configs/method_capabilities.yaml`.

## Running the Benchmark

Representative observed-time configs are stored in `benchmark/configs/` and
`benchmark/configs/runtime/`.

Example GSE242424 OSKM silver run:

```bash
python benchmark/methods/scNODE/run.py \
    --config benchmark/configs/scnode_gse242424_oskm_silver_A_hvg2000_formal.yaml
```

Example marker-FM transition silver runtime config:

```bash
python benchmark/methods/scNODE/run.py \
    --config benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_A_hvg2000_formal.yaml
```

For method-specific batch runs, use the root-level helper scripts such as:

```bash
bash run_marker_fm_transition_silver_primary_array.sh
bash run_gse242424_oskm_silver_formal_array.sh
```

The dispatcher and method runners read method capabilities and skip unsupported
metrics with explicit metadata rather than producing forced or invalid outputs.

## Ground Truth and Lineage Fidelity

Lineage Fidelity depends on three frozen assets:

1. a benchmark cell-state system;
2. a reference state-level lineage graph;
3. a fixed rule for aggregating method predictions to state-level transitions.

Current primary providers are registered in `benchmark/ground_truth/registry.yaml`:

- `gse178325_marker_fm_transition_silver_v1`
- `gse230659_marker_fm_transition_silver_v1`
- `gse242424_oskm_reprogramming_silver_v1`

The legacy `scgpt_v1` providers are retained for backward compatibility and
historical comparison only. New annotation systems should be added as new
versioned providers instead of changing existing frozen provider directories.

Per-run lineage outputs include `state_transition_matrix.csv`,
`lineage_graph_edges.csv`, `lineage_metrics.json`, and diagnostic metadata where
available.

## Relationship to scTimeBench

This project is best understood as a domain-specific extension of scTimeBench
ideas:

- scTimeBench provides the general benchmark vocabulary and metric families.
- trajectory specializes the benchmark to iPSC reprogramming.
- trajectory uses frozen silver-standard state systems and iPSC-specific
  reference graphs.
- trajectory keeps method capability gating explicit so transition-only methods
  and generative methods are compared on appropriate evaluation surfaces.

The goal is therefore not strict code-structure parity with scTimeBench, but a
scTimeBench-aligned benchmark that answers a narrower biological question.

## Framework Reference

`experimental framework v2.md` is the main design document for benchmark logic.
If implementation notes and README text diverge, treat the framework document
and current configs/results as the more specific source of truth.
