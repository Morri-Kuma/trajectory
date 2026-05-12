# Benchmark

This directory contains the scTimeBench-aligned benchmark framework used by the trajectory project for human chemical iPSC reprogramming.

The benchmark currently supports all three major evaluation dimensions for methods that provide the required outputs, while retaining capability gating for methods that only support a subset of tasks.

## Evaluation Dimensions

| Dimension | Status | Eligible methods |
|---|---|---|
| Forecast Accuracy | Active | Generative or projected-cell methods such as MIOFlow, PRESCIENT, and scNODE |
| Embedding Coherence | Active | Methods with projected embeddings or projected cells that can be embedded |
| Lineage Fidelity | Active | Methods with state-level transition predictions or convertible transition outputs |

WOT and CellRank2 are evaluated on Lineage Fidelity only. They are not forced into Forecast Accuracy or Embedding Coherence because they do not natively generate unseen-timepoint projected cells in the same output format as the generative methods.

## Directory Structure

```text
benchmark/
|-- README.md                       # this file
|-- adapters/                       # method wrappers and capability gates
|-- configs/                        # scenario, method, and run configs
|-- datasets/                       # dataset construction helpers
|-- docs/                           # framework notes and diagnostics
|-- evaluation/                     # metric dispatch and evaluators
|-- ground_truth/                   # frozen state labels and reference graphs
|-- inputs/                         # benchmark-ready input h5ad files
|-- methods/                        # vendored or wrapped method implementations
|-- reports/                        # aggregated summaries and figures
|-- results/                        # per-method per-scenario outputs
|-- scgpt/                          # scGPT annotation/provider utilities
`-- shared/                         # shared dataset and utility code
```

## Scenarios

| ID | Time axis | Type |
|---|---|---|
| A | Observed time | Interpolation / observed-time evaluation |
| B | Observed time | Extrapolation / held-out future-time evaluation |
| C | Observed time | Mixed interpolation and extrapolation |
| D | Pseudotime | Interpolation |
| E | Pseudotime | Extrapolation |
| F | Pseudotime | Mixed interpolation and extrapolation |

The currently reported formal iPSC runs focus on observed-time scenarios A/B/C.

## Current Result Sets

The canonical index for reportable, validation, pilot, and archived/debug result sets is:

- `results/result_manifest.yaml`

Important summary files:

- `results/summary_gse178325_0618_forecast_embedding.csv`
- `results/summary_lineage_metrics_preferred.csv`
- `reports/formal_benchmark_summary.csv`
- `reports/official/gse230659_projected_hvg_annotation_audit.md`

Per-run outputs usually include:

- `run_metadata.json`
- `forecast_metrics.json`
- `embedding_metrics.json`
- `lineage_metrics.json`
- `per_timepoint_forecast_metrics.csv`
- `per_timepoint_embedding_metrics.csv`
- `state_transition_matrix.csv`
- `lineage_graph_edges.csv`

Not every method produces every file. Missing files should be interpreted through the method capability configuration rather than treated automatically as failed runs.

## Running a Config

Representative configs are stored in `configs/`.

```bash
python benchmark/evaluation/eval_dispatch.py \
    --config benchmark/configs/scnode_gse178325_observed_0618_hvg2000_A.yaml
```

Root-level helper scripts are available for array-style runs, for example:

```bash
bash run_scnode_gse178325_0618_hvg_array.sh
bash run_mioflow_gse178325_0618_hvg_array.sh
bash run_prescient_gse178325_0618_hvg_array.sh
```

## Capability Gating

Method capabilities are declared in `configs/method_capabilities.yaml` and reflected in the adapter layer.

The dispatcher uses these flags to decide which evaluators to run:

- Forecast Accuracy requires projected expression at held-out or future time points.
- Embedding Coherence requires projected embeddings or projectable cells.
- Lineage Fidelity requires state-level transition scores, transition matrices, or outputs that can be aggregated into transitions.

Unsupported metrics are skipped explicitly instead of being approximated with method-inappropriate workarounds.

## Ground Truth Requirements

Lineage Fidelity requires:

1. a frozen cell-state system;
2. a frozen reference lineage graph;
3. a fixed aggregation rule from method outputs to state-level transitions.

Current iPSC runs use scGPT-derived pseudo-state annotations and reference graph providers under `ground_truth/`.

## Adding a Method

1. Add or vendor the method implementation under `methods/`.
2. Create or update an adapter under `adapters/`.
3. Register the method in the dispatcher.
4. Add capability flags to `configs/method_capabilities.yaml`.
5. Add one or more run configs under `configs/`.
6. Verify that unsupported metrics are skipped and supported metrics produce the expected files.
