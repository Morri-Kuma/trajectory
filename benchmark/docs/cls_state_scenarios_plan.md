# CLS-State Scenario Design Plan

## 1. Purpose

This document proposes a supplementary benchmark module that replaces the inactive pseudotime scenarios D-F with classifier-defined ordered cell-state scenarios.

The current benchmark uses observed experimental time for scenarios A-C and keeps pseudotime scenarios D-F inactive because DPT pseudotime inferred from the full expression matrix can introduce data leakage and circularity. The proposed CLS-state design avoids treating internally inferred pseudotime as an independent time axis. Instead, it uses an independently frozen classifier-derived cell-state provider as a state-order axis.

In this document, `CLS` means **classifier-defined latent state** or **classifier-defined cell state**. It does not assume that the referenced scBERT paper uses a standard BERT `[CLS]` token. Yang et al. 2022 scBERT uses a Performer encoder, gene/expression embeddings, one-dimensional convolution, and a classifier head for cell type annotation.

## 2. Benchmark Role

The CLS-state scenarios should be treated as supplementary, not as a replacement for the primary observed-time benchmark.

| Scenario family | Role | Time/state axis | Report status |
| --- | --- | --- | --- |
| A-C | Primary benchmark | Observed experimental time, for example `abs_day` | Official primary ranking |
| D-F | Supplementary CLS-state benchmark | Frozen classifier-defined ordered states | Separate supplementary ranking |

The main question becomes:

> Given an independently frozen classifier-defined state axis, can temporal trajectory models predict held-out state distributions, preserve state structure, and recover reference state transitions?

## 3. Scenario Redefinition

The previous D-F definitions were pseudotime-based:

- D: pseudotime interpolation
- E: pseudotime extrapolation
- F: pseudotime interpolation plus extrapolation

The proposed definitions are:

- D: CLS-state interpolation
- E: CLS-state extrapolation
- F: CLS-state interpolation plus extrapolation

The state axis is discrete but ordered. Each cell receives both a state label and an order index:

```text
cls_state_label
cls_state_order_numeric
cls_state_confidence
cls_state_status
```

Example for iPSC reprogramming:

```text
hADSCs                  0
epithelial_like         1
intermediate_plastic    2
xen_like                3
hCiPS                   4
```

Example for organoid differentiation:

```text
pluripotent             0
neuroectoderm           1
NPC                     2
neuron_branch           3
glia_branch             3
mature_neuron           4
```

Branching states may share the same order index while remaining distinct `cls_state_label` values for lineage and embedding evaluation.

## 4. Core Design Rules

1. A-C remain the primary benchmark scenarios.
2. D-F become supplementary CLS-state scenarios and must be reported separately.
3. A frozen CLS-state provider must be built before model evaluation.
4. The CLS-state provider must not use embeddings, projected cells, transition matrices, or outputs from any evaluated model.
5. Held-out cells must not be used to train or tune the CLS provider in a way that leaks evaluation information.
6. Low-confidence, ambiguous, and out-of-distribution cells must be tracked and excluded from official CLS-state metrics.
7. The state order and reference graph must be frozen from literature, author annotations, or an independent reference policy. They should not be automatically inferred from evaluated model outputs.

## 5. Provider Specification

Add a new frozen provider under:

```text
benchmark/ground_truth/providers/{provider_id}/
```

Suggested provider IDs:

```text
gse230659_cls_state_scbert_lodo_v1
gse178325_cls_state_scbert_lodo_v1
gse242424_cls_state_scbert_lodo_v1
```

Recommended provider files:

```text
state_labels.tsv
state_metadata.tsv
reference_graph.json
reference_graph_edges.csv
annotation_votes.tsv
ground_truth_metadata.json
cls_provider_model_metadata.json
```

Recommended `state_labels.tsv` columns:

```text
cell_id
cls_state_label
cls_state_order_numeric
cls_state_confidence
cls_state_status
cls_state_source
```

Recommended excluded labels:

```text
ambiguous
unknown_or_ood
unknown_or_low_confidence
```

The provider should be registered in:

```text
benchmark/ground_truth/registry.yaml
```

Example registry entry:

```yaml
gse230659_cls_state_scbert_lodo_v1:
  provider_id: gse230659_cls_state_scbert_lodo_v1
  annotation_method: scbert_cls_state_provider
  status: frozen_supplementary_cls_state_provider
  state_key: cls_state_label
  time_key: cls_state_order_numeric
  label_path: benchmark/ground_truth/providers/gse230659_cls_state_scbert_lodo_v1/state_labels.tsv
  metadata_path: benchmark/ground_truth/providers/gse230659_cls_state_scbert_lodo_v1/state_metadata.tsv
  reference_graph_path: benchmark/ground_truth/providers/gse230659_cls_state_scbert_lodo_v1/reference_graph.json
  reference_edges_path: benchmark/ground_truth/providers/gse230659_cls_state_scbert_lodo_v1/reference_graph_edges.csv
  confidence_mode: all
  exclude_uncertain_states: true
  excluded_labels:
    - ambiguous
    - unknown_or_ood
    - unknown_or_low_confidence
  label_mode: cls_state_silver
  label_type: frozen_classifier_defined_state
  analysis_role: supplementary_cls_state
```

## 6. CLS Provider Construction

The recommended implementation uses a reference-guided classifier, such as scBERT, to assign state labels and confidence values.

### 6.1 Input Requirements

scBERT is designed around full-gene expression input rather than HVG-only input. Therefore, the provider-building step should prefer full-gene h5ad inputs.

Minimum input:

```text
X: full-gene expression matrix
obs[cell_id]
obs[sample_id]
obs[abs_day] or equivalent observed-time key
obs[reference_state_label] for reference/training cells
var[gene_symbol]
```

If only HVG2000 is available, the provider should be explicitly named as an HVG-limited provider, for example:

```text
gse230659_cls_state_scbert_hvg2000_sensitivity_v1
```

It should not be presented as a full scBERT-style provider.

### 6.2 Anti-Leakage Training Policies

Use one of the following policies:

| Policy | Description | Strength |
| --- | --- | --- |
| Leave-one-dataset-out | Train or fine-tune on one dataset, map another dataset | Strongest for cross-dataset validation |
| Leave-one-line/organoid-out | Train on some cell lines or organoids, map held-out line/organoid | Recommended for organoid data |
| Train-reference-only | Train on an explicitly defined reference subset, map benchmark cells | Acceptable if the reference split is fixed before evaluation |
| Same-dataset all-cell training | Train and evaluate on the same cells | Not acceptable for official CLS-state metrics |

The provider metadata must record:

```text
training_datasets
query_dataset
heldout_policy
gene_vocabulary
preprocessing
model_checkpoint
confidence_threshold
excluded_labels
date_created
code_version
```

## 7. Scenario Configuration

Create:

```text
benchmark/configs/scenario_cls.yaml
```

Suggested content:

```yaml
time_axis: cls_state_order
time_key: cls_state_order_numeric
cell_state_key: cls_state_label
provider_type: frozen_cls_state_provider
source_method: scbert_classifier
analysis_role: supplementary_cls_state

scenarios:
  D:
    id: D
    name: "CLS-state interpolation"
    description: >
      Hold out one or more internal classifier-defined ordered states.
      Train on surrounding states and evaluate predictions on the held-out
      internal states.
    dimension_eligibility:
      forecast_accuracy: true_if_method_supports_projection
      embedding_coherence: true_if_method_supports_projection
      lineage_fidelity: true_if_method_supports_lineage_inference
    scenario_params:
      train_times: [0, 1, 3, 4]
      heldout_times: [2]
      test_includes_start: true

  E:
    id: E
    name: "CLS-state extrapolation"
    description: >
      Train on early classifier-defined ordered states and evaluate on later
      unseen states.
    dimension_eligibility:
      forecast_accuracy: true_if_method_supports_projection
      embedding_coherence: true_if_method_supports_projection
      lineage_fidelity: true_if_method_supports_lineage_inference
    scenario_params:
      train_times: [0, 1, 2]
      heldout_times: [3, 4]
      test_includes_start: true

  F:
    id: F
    name: "CLS-state interpolation + extrapolation"
    description: >
      Train on a subset of classifier-defined ordered states and evaluate on
      both internal and later held-out states.
    dimension_eligibility:
      forecast_accuracy: true_if_method_supports_projection
      embedding_coherence: true_if_method_supports_projection
      lineage_fidelity: true_if_method_supports_lineage_inference
    scenario_params:
      train_times: [0, 1, 3]
      heldout_times: [2, 4]
      test_includes_start: true
```

The exact state splits should be dataset-specific and must be frozen in runtime configs.

## 8. Dataset Configuration

Add CLS-state dataset entries to:

```text
benchmark/shared/dataset/default_datasets.yaml
```

Example:

```yaml
- name: GSE230659Dataset
  tag: gse230659_cls_state_scenario_D
  id: GSE230659
  h5ad_path: benchmark/inputs/gse230659_cls_state_hvg2000/GSE230659_cls_state_HVG2000_benchmark_input.h5ad
  time_key: cls_state_order_numeric
  cell_state_key: cls_state_label
  data_preprocessing_steps:
    - name: ScenarioTimepointSplit
      time_key: cls_state_order_numeric
      train_times: [0, 1, 3, 4]
      heldout_times: [2]
```

The benchmark-ready h5ad must contain:

```text
obs[cls_state_label]
obs[cls_state_order_numeric]
obs[cls_state_confidence]
obs[cls_state_status]
```

## 9. Metric Interpretation

The existing three benchmark dimensions remain valid, but the evaluation axis changes from observed time or DPT pseudotime to frozen CLS-state order.

### 9.1 Forecast Accuracy

Forecast Accuracy compares:

```text
projected expression at held-out CLS state/order
vs
observed expression of cells assigned to the same held-out CLS state/order
```

Metrics remain:

```text
wasserstein_distance
gaussian_mmd
energy_distance_mmd
hausdorff_loss
```

This is not a comparison between predicted CLS labels and reference CLS labels. The CLS provider defines the evaluation target groups.

### 9.2 Embedding Coherence

Embedding Coherence compares:

```text
Leiden clusters from projected embedding
vs
frozen cls_state_label
```

Metrics remain:

```text
adjusted_rand_index
mean_normalized_entropy
weighted_mean_normalized_entropy
```

Important caveat: if the CLS provider is itself based on a model embedding, avoid evaluating that same embedding family as if it were independent biological truth.

### 9.3 Lineage Fidelity

Lineage Fidelity compares:

```text
method state-transition matrix
vs
frozen CLS-state reference graph
```

Metrics remain:

```text
auroc
auprc
jaccard_similarity_topk
single_step_recovery
multi_step_recovery
```

The reference graph must be frozen independently from evaluated model outputs.

## 10. Optional CLS-Specific Diagnostics

The following diagnostics may be added to the supplementary report but should not enter the first official CLS-state ranking:

```text
cls_state_distribution_distance
cls_unknown_rate
cls_low_confidence_rate
cls_confidence_mean
cls_confidence_by_state
cls_order_monotonicity
cls_state_transition_entropy
soft_cls_probability_kl
```

These diagnostics answer different questions from the three core benchmark dimensions. They should be reported as interpretability and provider-quality checks.

## 11. Minimal Viable Implementation

Recommended MVP:

1. Select one dataset, preferably GSE230659 or GSE242424.
2. Build or import a full-gene benchmark h5ad.
3. Build one frozen CLS-state provider.
4. Create `scenario_cls.yaml`.
5. Create CLS-state runtime configs for:
   - scNODE D/E/F
   - MIOFlow D/E/F
   - PRESCIENT D/E/F
   - WOT D/E/F lineage-only
   - CellRank2 D/E/F lineage-only
6. Run a small smoke test on capped cells.
7. Run full formal supplementary jobs.
8. Generate a separate report under:

```text
benchmark/reports/cls_state/
```

Recommended report files:

```text
cls_state_method_summary.md
cls_state_forecast_summary.csv
cls_state_embedding_summary.csv
cls_state_lineage_summary.csv
cls_state_provider_qc.md
```

## 12. Implementation Checklist

### Provider Layer

- [ ] Define provider ID and provider directory.
- [ ] Decide the training/query anti-leakage policy.
- [ ] Generate `state_labels.tsv`.
- [ ] Generate `state_metadata.tsv`.
- [ ] Generate `reference_graph.json`.
- [ ] Generate `reference_graph_edges.csv`.
- [ ] Generate `annotation_votes.tsv`.
- [ ] Generate `ground_truth_metadata.json`.
- [ ] Register provider in `benchmark/ground_truth/registry.yaml`.

### Input Layer

- [ ] Build full-gene h5ad where feasible.
- [ ] Add CLS-state labels and order fields to h5ad obs.
- [ ] Produce HVG2000 derivative only after full-gene provider assignment, if needed by current model runners.
- [ ] Add dataset entries to `benchmark/shared/dataset/default_datasets.yaml`.

### Scenario Layer

- [ ] Add `benchmark/configs/scenario_cls.yaml`.
- [ ] Freeze dataset-specific D/E/F train and held-out state orders.
- [ ] Ensure D-F no longer refer to DPT pseudotime in active configs.

### Method Layer

- [ ] Confirm all adapters read `time_key=cls_state_order_numeric`.
- [ ] Confirm all adapters read `cell_state_key=cls_state_label`.
- [ ] Keep capability gating unchanged:
  - scNODE, MIOFlow, PRESCIENT: Forecast, Embedding, Lineage
  - WOT, CellRank2: Lineage only

### Evaluation Layer

- [ ] Reuse existing forecast evaluator.
- [ ] Reuse existing embedding milestone evaluator.
- [ ] Reuse existing lineage evaluator.
- [ ] Add CLS-specific diagnostics only as supplementary outputs.

### Reporting Layer

- [ ] Create `benchmark/reports/cls_state/`.
- [ ] Keep CLS-state results separate from `official_silver`.
- [ ] Report provider ID, model checkpoint, confidence threshold, excluded labels, and leakage policy in every summary.

## 13. Main Risks and Mitigations

| Risk | Why it matters | Mitigation |
| --- | --- | --- |
| CLS provider trained and evaluated on the same cells | Creates circularity and overstates performance | Use leave-one-dataset, leave-one-line, or frozen reference-only policy |
| HVG2000 input weakens scBERT assumptions | scBERT was designed for full-gene input | Build full-gene provider first; label HVG-only provider as sensitivity |
| CLS states lack biological direction | D-F require interpolation/extrapolation along an order | Freeze `cls_state_order_numeric` from literature or author annotations |
| Reference graph is inferred from evaluated model outputs | Creates circular lineage evaluation | Build graph from provider/literature before benchmark runs |
| Unknown states are forced into known labels | Inflates metrics and hides OOD biology | Use confidence threshold and explicit excluded labels |
| Branching states collapse into one order | May erase lineage topology | Use shared order indices for parallel branches while preserving distinct labels |

## 14. Recommended Interpretation

CLS-state D-F should be described as:

> A supplementary reference-guided state-axis benchmark that evaluates temporal models under a frozen classifier-defined cell-state ordering.

Avoid claiming:

> The CLS-state axis is independent biological ground truth.

Recommended wording:

> We evaluate D-F against a frozen classifier-defined state provider, not against experimentally lineage-traced ground truth. The provider is used to define ordered state groups and a reference transition graph. Results are interpreted as performance against this fixed state system and are reported separately from the observed-time primary benchmark.

## 15. Summary

The CLS-state design is feasible and potentially innovative if it is implemented as a frozen supplementary provider layer. Its main value is that it replaces full-data DPT pseudotime with an externally constrained, confidence-aware, classifier-defined state system. This can reduce pseudotime leakage and provide a bridge between marker-defined silver labels, foundation-model-based annotation, and trajectory method evaluation.

The design should not replace A-C. It should extend the benchmark with a separate D-F supplementary module focused on classifier-defined ordered state dynamics.
