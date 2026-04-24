# Ground Truth Providers

This directory defines the benchmark's replaceable annotation / lineage
reference layer.

Each provider must expose the same assets:

- `state_labels.tsv`: per-cell state labels
- `state_metadata.tsv`: per-state metadata
- `reference_graph.json`: state-to-state reference lineage graph
- `reference_graph_edges.csv`: tabular copy of the reference graph edges
- `ground_truth_metadata.json`: provenance for the provider export

The active provider registry is `registry.yaml`. Method configs may either:

1. set `ground_truth.provider_id` to a registered provider, or
2. provide a full `ground_truth` block with explicit paths, or
3. rely on legacy `state_system` / `lineage` fields.

The current provider, `scgpt_v1`, is a silver-standard working reference
derived from scGPT embeddings and consecutive-timepoint kNN transition
counting. Future providers such as CellTypist, scANVI, SingleR, or marker-based
annotation should be added as new provider entries rather than changing model
adapters.
