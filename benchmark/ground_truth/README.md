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

Method configs may either:

1. set `ground_truth.provider_id` to a registered provider;
2. provide a full `ground_truth` block with explicit paths; or
3. rely on legacy `state_system` / `lineage` fields.

## Current Reportable Providers

The current primary benchmark providers are frozen silver-standard milestone
providers:

| Provider | Dataset | Role | Notes |
|---|---|---|---|
| `gse178325_marker_fm_transition_silver_v1` | GSE178325 | primary report | Marker-FM transition silver provider; excludes `ambiguous` and `unknown_or_ood` from official metrics. |
| `gse230659_marker_fm_transition_silver_v1` | GSE230659 | primary report | Marker-FM transition silver provider; excludes `ambiguous` and `unknown_or_ood` from official metrics. |
| `gse242424_oskm_reprogramming_silver_v1` | GSE242424 | primary report | Author-cluster OSKM reprogramming silver provider on the 59,187-cell matched subset. |

Milestone consensus, embedding, and classifier providers for GSE178325 and
GSE230659 are retained for sensitivity analyses. The `scgpt_v1` and
`scgpt_v1_gse178325_0618` providers are deprecated legacy pseudostate
references retained for backward compatibility with older configs and reports;
they should not be used as primary benchmark ground truth.

## Provider Policy

Do not change a frozen provider in place to match model outcomes. Add a new
versioned provider entry when the annotation policy, label universe, graph
topology, exclusion policy, or source data changes.

Lineage Fidelity depends on a fixed state system, a fixed reference graph, and a
fixed aggregation rule from method outputs to state-level transitions. Keep those
choices explicit in the provider metadata and in the run config.
