# Benchmark

This directory contains the scTimeBench-aligned benchmark framework used by the
trajectory project for human iPSC reprogramming.

The benchmark supports all three major evaluation dimensions for methods that
provide the required outputs, while retaining capability gating for methods that
only support a subset of tasks.

## Evaluation Dimensions

| Dimension | Status | Eligible methods |
|---|---|---|
| Forecast Accuracy | Active | Generative or projected-cell methods such as MIOFlow, PRESCIENT, and scNODE |
| Embedding Coherence | Active | Methods with projected embeddings or projected cells that can be embedded |
| Lineage Fidelity | Active | Methods with state-level transition predictions or convertible transition outputs |

WOT and CellRank2 are evaluated on Lineage Fidelity only. They are not forced
into Forecast Accuracy or Embedding Coherence because they do not natively
generate unseen-timepoint projected cells in the same output format as the
generative methods.

## Directory Structure

```text
benchmark/
|-- README.md                       # this file
|-- adapters/                       # method wrappers and capability gates
|-- annotation/                     # annotation/provider-building utilities
|-- configs/                        # scenario, method, runtime, and run configs
|-- datasets/                       # legacy dataset construction helpers
|-- docs/                           # framework notes and diagnostics
|-- evaluation/                     # metric dispatch and evaluators
|-- ground_truth/                   # frozen state labels and reference graphs
|-- inputs/                         # benchmark-ready input h5ad files
|-- methods/                        # vendored or wrapped method implementations
|-- reports/                        # aggregated summaries and figures
|-- results/                        # per-method per-scenario outputs
|-- scgpt/                          # legacy scGPT annotation/provider utilities
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

The currently reported formal runs focus on observed-time scenarios A/B/C.

## Current Result Sets

The canonical index for reportable, validation, pilot, and archived/debug result
sets is:

- `results/result_manifest.yaml`

Current reportable summaries:

- `reports/official_silver/official_silver_model_rankings.md`
- `reports/gse242424_oskm_silver/gse242424_oskm_silver_report.md`
- `reports/core_summary.csv`
- `reports/embedding_summary.csv`
- `reports/lineage_summary.csv`
- `reports/official/gse230659_projected_hvg_annotation_audit.md`

Legacy summaries retained for comparison:

- `results/summary_gse178325_0618_forecast_embedding.csv`
- `results/summary_lineage_metrics_preferred.csv`
- `reports/formal_benchmark_summary.csv`

Per-run outputs usually include:

- `run_metadata.json`
- `forecast_metrics.json`
- `embedding_metrics.json`
- `lineage_metrics.json`
- `per_timepoint_forecast_metrics.csv`
- `per_timepoint_embedding_metrics.csv`
- `state_transition_matrix.csv`
- `lineage_graph_edges.csv`

Not every method produces every file. Missing files should be interpreted
through the method capability configuration rather than treated automatically as
failed runs.

## Current Ground Truth

Current reportable runs use frozen silver-standard providers under
`ground_truth/providers/`:

| Provider | Dataset | Result class |
|---|---|---|
| `gse178325_marker_fm_transition_silver_v1` | GSE178325 | `marker_fm_transition_silver_formal` |
| `gse230659_marker_fm_transition_silver_v1` | GSE230659 | `marker_fm_transition_silver_formal` |
| `gse242424_oskm_reprogramming_silver_v1` | GSE242424 | `gse242424_oskm_silver_formal` |

The GSE242424 provider is built from the 59,187-cell author-cluster-matched
subset, not the full local 156,969-cell GSE242424 input. The legacy scGPT
pseudostate providers remain available for backward compatibility and historical
reports, but they are not the current primary benchmark reference.

## Running a Config

Representative configs are stored in `configs/` and `configs/runtime/`.

```bash
python benchmark/methods/scNODE/run.py \
    --config benchmark/configs/scnode_gse242424_oskm_silver_A_hvg2000_formal.yaml
```

```bash
python benchmark/methods/scNODE/run.py \
    --config benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_A_hvg2000_formal.yaml
```

Root-level helper scripts are available for array-style runs, for example:

```bash
bash run_marker_fm_transition_silver_primary_array.sh
bash run_marker_fm_transition_silver_summarize_and_validate.sh
bash run_gse242424_oskm_silver_formal_array.sh
```

## Capability Gating

Method capabilities are declared in `configs/method_capabilities.yaml` and
reflected in the adapter layer.

The dispatcher uses these flags to decide which evaluators to run:

- Forecast Accuracy requires projected expression at held-out or future time
  points.
- Embedding Coherence requires projected embeddings or projectable cells.
- Lineage Fidelity requires state-level transition scores, transition matrices,
  or outputs that can be aggregated into transitions.

Unsupported metrics are skipped explicitly instead of being approximated with
method-inappropriate workarounds.

## Ground Truth Requirements

Lineage Fidelity requires:

1. a frozen cell-state system;
2. a frozen reference lineage graph;
3. a fixed aggregation rule from method outputs to state-level transitions.

Provider details are documented in `ground_truth/README.md` and registered in
`ground_truth/registry.yaml`.

## Adding a Method

1. Add or vendor the method implementation under `methods/`.
2. Create or update an adapter under `adapters/`.
3. Register the method in the dispatcher or method runner.
4. Add capability flags to `configs/method_capabilities.yaml`.
5. Add one or more run configs under `configs/` or `configs/runtime/`.
6. Verify that unsupported metrics are skipped and supported metrics produce the
   expected files.
