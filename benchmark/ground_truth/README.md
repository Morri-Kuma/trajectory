# Ground Truth Providers

This directory defines the benchmark's replaceable annotation and lineage
reference layer. The active registry is `registry.yaml`.

Frozen reportable providers expose the following assets:

- `state_labels.tsv`: per-cell state labels.
- `state_metadata.tsv`: per-state metadata and cell counts.
- `reference_graph.json`: state-to-state reference lineage graph.
- `reference_graph_edges.csv`: tabular copy of the reference graph edges.
- `ground_truth_metadata.json`: provenance for the provider export.
- `annotation_votes.tsv`: provider-specific label provenance where available.

Method configs should:

1. set `ground_truth.provider_id` to a registered provider;
2. provide a full `ground_truth` block with explicit paths.

Legacy `state_system` / `lineage` fallbacks have been removed from the
reportable benchmark configuration surface.

## Current Reportable Providers

The current primary benchmark providers are frozen milestone ground-truth
providers:

| Provider | Dataset | Role | Notes |
|---|---|---|---|
| `gse178325_marker_fm_transition_silver_v1` | GSE178325 | primary report | Marker-FM transition silver provider; excludes `ambiguous` and `unknown_or_ood` from official metrics. |
| `gse230659_marker_fm_transition_silver_v1` | GSE230659 | primary report | Marker-FM transition silver provider; excludes `ambiguous` and `unknown_or_ood` from official metrics. |
| `gse242424_oskm_reprogramming_ground_truth_v1` | GSE242424 | primary report | Author-cluster OSKM reprogramming ground truth provider on the 59,187-cell matched subset. |

Older pseudostate, stage-proxy, and placeholder milestone providers have been
removed from the active registry and should not be referenced by new reports or
method configs.

## Provider Policy

Do not change a frozen provider in place to match model outcomes. Add a new
versioned provider entry when the annotation policy, label universe, graph
topology, exclusion policy, or source data changes.

Lineage Fidelity depends on a fixed state system, a fixed reference graph, and a
fixed aggregation rule from method outputs to state-level transitions. Keep those
choices explicit in the provider metadata and in the run config.
