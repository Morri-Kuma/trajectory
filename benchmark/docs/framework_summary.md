# Framework Design Summary

scTimeBench-aligned benchmark for human chemical reprogramming to iPSCs.
Source of truth: `docs/framework/experimental_framework_v2.md`.

## Core principle

- Methods that can generate unseen future-timepoint cells → evaluate on all three dimensions.
- Methods that cannot generate unseen future-timepoint cells → evaluate on **Lineage Fidelity only**.

This rule applies directly to WOT and CellRank2, which are OT-style methods
that infer lineage structure without projecting cells to new time points.

## Three core dimensions

### 1. Forecast Accuracy
Evaluates whether a method can project cells from time `t` to unseen time `t+1`
with gene expression aligned to the observed cells at `t+1`.

**Not applicable to WOT or CellRank2.** Will be activated for future generative models.

Metrics: Wasserstein Distance, Gaussian MMD, Energy Distance MMD, Hausdorff Loss.

### 2. Embedding Coherence
Evaluates whether projected cells preserve biologically meaningful cellular structure.

**Not applicable to WOT or CellRank2.** Will be activated for future generative models.

Metrics: Adjusted Rand Index, average normalized classifier entropy.

### 3. Lineage Fidelity
Evaluates whether a method recovers a cell-state transition structure consistent
with a benchmark reference lineage.

**Applicable to WOT and CellRank2.** This is the only active dimension in the
current benchmark stage.

Metrics: AUROC, AUPRC, Jaccard Similarity, single-step recovery, multi-step recovery.
Baseline: correlation-based baseline (analogous to scTimeBench's lineage baseline).

## Six benchmark scenarios

| ID | Time axis | Type |
|---|---|---|
| A | Observed time | Interpolation |
| B | Observed time | Extrapolation |
| C | Observed time | Interpolation + extrapolation |
| D | Pseudotime | Interpolation |
| E | Pseudotime | Extrapolation |
| F | Pseudotime | Interpolation + extrapolation |

## Method capability flags

Each adapter declares `supports_unseen_timepoint_projection` and
`supports_lineage_inference`. The dispatcher (`eval_dispatch.py`) reads
these flags and routes each method only to its eligible evaluators.

**WOT and CellRank2 (current stage):**
- `supports_unseen_timepoint_projection = False`
- `supports_lineage_inference = True`
- Eligible for: Lineage Fidelity only

**Future generative models:**
- `supports_unseen_timepoint_projection = True`
- `supports_lineage_inference = True`
- Eligible for: all three dimensions

## What is intentionally not used now

The following assets are outside the active benchmark core per v2 §15:
- `celltype_labels_v1.tsv` / `pseudotime_labels_v1.tsv` / `terminal_labels_v1*.tsv`
- `reference_graph_v1.json`
- Previous T1/T2/T3/T4/T5 task structure
- Label-driven terminal identification evaluation
- OCLR-based endpoint evaluation flow

These may be reconsidered later as supplementary analyses after the core
scTimeBench-aligned Lineage Fidelity pipeline is stable and validated.

## Ranking rule

- Rank methods within each metric.
- Average metric ranks for the dimension rank within each scenario.
- For WOT and CellRank2: report Lineage Fidelity rank only.
- Do NOT fabricate a combined overall rank for dimensions a method is ineligible for.
