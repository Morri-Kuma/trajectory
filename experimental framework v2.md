# Experimental Framework v2
## scTimeBench-aligned benchmark framework for human chemical reprogramming trajectory analysis

## 1. Project overview

This repository is being developed into a benchmark-oriented single-cell trajectory analysis project for **human chemical reprogramming to iPSCs**.

The benchmark framework is now designed to follow the **same core evaluation dimensions as scTimeBench**, rather than the previously drafted five-task project-specific design. The three core benchmark dimensions are:

1. **Forecast Accuracy**
2. **Embedding Coherence**
3. **Lineage Fidelity**

The current first-stage comparison focuses on **WOT** and **CellRank2**. Because these methods do **not** directly generate unseen future-timepoint cells, they are evaluated only on **Lineage Fidelity**, consistent with the logic used in scTimeBench for OT-based methods. Forecast Accuracy and Embedding Coherence remain part of the framework, but they are activated only for future models that can generate projected cells at unseen time points.

---

## 2. Benchmark goal

Under a unified benchmark structure, compare temporal / trajectory inference methods using the **same three core dimensions defined by scTimeBench**, while adapting the benchmark content to the biological problem of **chemical reprogramming toward pluripotency**.

The benchmark therefore asks three top-level questions:

### 2.1 Forecast Accuracy
Can a method project cells from timepoint `t` to an unseen timepoint `t+1` such that the predicted gene expression is well aligned with the observed cells at `t+1`?

### 2.2 Embedding Coherence
Can a method preserve biologically meaningful cellular structure when cells are projected forward in time?

### 2.3 Lineage Fidelity
Can a method recover a cell-state transition structure that is consistent with a benchmark reference lineage?

---

## 3. Core benchmark policy

### 3.1 Same top-level dimensions as scTimeBench
The framework must always retain the following three top-level benchmark dimensions:

- Forecast Accuracy
- Embedding Coherence
- Lineage Fidelity

No additional project-specific task may be promoted to the same top level unless the benchmark is intentionally redesigned in the future.

### 3.2 Method eligibility follows scTimeBench logic
Methods are evaluated only on the dimensions they are biologically and technically able to support.

#### Rule A: Methods that can generate unseen future-timepoint cells
These methods are eligible for:
- Forecast Accuracy
- Embedding Coherence
- Lineage Fidelity

#### Rule B: Methods that cannot generate unseen future-timepoint cells
These methods are eligible for:
- Lineage Fidelity only

This is the rule that applies to **WOT** and **CellRank2** in the current benchmark stage. No workaround, projection adapter, pseudo-generated future cells, or transport-based synthetic projection layer should be introduced merely to force such methods into Forecast Accuracy or Embedding Coherence. This directly follows the scTimeBench treatment of OT-based methods, which were excluded from the first two benchmark dimensions because they cannot predict new cells at unseen time points.

### 3.3 Core benchmark and supplementary materials must be separated
At the current stage, the benchmark should focus on making the three core dimensions logically clean and code-stable.

**Updated metric policy, 2026-05-08.** The benchmark still keeps the same three scTimeBench-level metric families, but the annotation-dependent parts of the benchmark now use **marker-defined milestone labels** rather than scGPT pseudostates as the primary biological state system.

The active metric design is:

1. **Forecast Accuracy**
   - label-free / annotation-independent
   - retains the scTimeBench expression-distribution metrics without milestone labels

2. **Embedding Coherence**
   - uses the frozen `official_silver` milestone labels for ARI
   - clusters method-projected embeddings with the scTimeBench-style
     unsupervised Leiden route
   - compares those clusters against `final_milestone_label_coarse`
   - does not use projected marker labels, copied milestone labels, or
     consensus labels as clustering input

3. **Lineage Fidelity**
   - uses `official_silver` milestone transition matrices and
     marker-defined milestone reference graphs
   - aggregates model transitions by `final_milestone_label_coarse`

The active annotation-dependent state system is `official_silver`, with
`state_key: final_milestone_label_coarse`. Earlier `consensus`,
`embedding_based`, and `classifier_based` provider modes are historical design
iterations and are not part of the currently executable primary benchmark.

Items such as:
- terminal success labels
- pseudotime v1 files
- custom coarse cell-state labels
- existing reference graph drafts
- endpoint diagnostics
- early-fate diagnostics
- biology interpretation modules
- robustness / runtime reporting

are **not part of the active core benchmark flow right now**. They may be reintroduced later as supplementary analyses after the main benchmark dimensions are implemented and validated. This is a deliberate change from the earlier project-specific five-task structure.

**Deprecated note.** The paragraph above described the initial framework-rebuild stage before milestone annotation was promoted into the metrics design. It is retained for provenance only. In the updated metrics plan, marker-defined milestone labels are no longer considered an optional biology interpretation layer; they are the primary state system for Embedding Coherence and Lineage Fidelity.

---

## 4. Benchmark scenarios

The benchmark keeps the same scenario logic as a temporal benchmarking framework, but scenario activation depends on method eligibility.

### 4.1 Observed-time scenarios

#### Scenario A: Observed-time interpolation
- Time axis: observed sampling time
- Goal: project cells to held-out internal time points

#### Scenario B: Observed-time extrapolation
- Time axis: observed sampling time
- Goal: project cells to later unseen time points

#### Scenario C: Observed-time interpolation + extrapolation
- Time axis: observed sampling time
- Goal: benchmark performance under the hardest observed-time setting

### 4.2 Pseudotime scenarios（on hold）

Pseudotime scenarios are specified for framework completeness, but they are **out of scope for the current GSE230659 manuscript-stage report**. The primary dataset has no biological or technical replicates across time points, so pseudotime derived from the same time-confounded expression matrix would be difficult to interpret as an independent temporal axis. Scenarios D-F should therefore be reported as planned supplementary work unless an alternative dataset or independent pseudotime provider is introduced.

#### Scenario D: Pseudotime interpolation
- Time axis: pseudotime bins

#### Scenario E: Pseudotime extrapolation
- Time axis: pseudotime bins

#### Scenario F: Pseudotime interpolation + extrapolation
- Time axis: pseudotime bins

### 4.3 Scenario activation policy
For a given method:

- If the method supports unseen-timepoint cell generation, it may enter Scenarios A–F for all three core dimensions.
- If the method does not support unseen-timepoint cell generation, it may only enter **Lineage Fidelity** evaluation under the relevant scenario setup.

For the current benchmark stage:
- **WOT**: Lineage Fidelity only
- **CellRank2**: Lineage Fidelity only
- **scNODE, PRESCIENT, MIOFlow**: active for Scenarios A-C across Forecast Accuracy, Embedding Coherence, and Lineage Fidelity
- **Scenarios D-F**: inactive for the current manuscript-stage report

### 4.4 Primary dataset limitations

GSE230659 is a single-cell-line, single-library-per-timepoint series. Time point and library are therefore perfectly confounded: every observed time point is represented by one library from cell line 0618, with no independent donors, no technical replicates, and no statistical design that can separate temporal biology from library-specific technical variation.

This limitation affects all metric interpretation. A method that recovers the expected lineage may be learning a true reprogramming trajectory, library-specific artifacts, or both. Formal results on this dataset should therefore be described as performance against a frozen silver-standard temporal reference, not as proof that a method has recovered replicate-validated biology.

The minimum sensitivity plan for submission is:
- report leave-one-timepoint-out metric stability for Scenarios A-C where feasible
- explicitly inspect the Stage I Day 4 sample because it has the smallest post-QC cell count
- explicitly inspect the hCiPSC endpoint because it has the highest low-quality-cell rate before filtering
- repeat Lineage Fidelity with `edge_confidence_mode: high_only` as a supplementary threshold sensitivity analysis

---

## 5. Unified input specification

All methods should consume the same standardized dataset object at the benchmark level.

### Main input object
```text
adata_benchmark.h5ad
```

### Required minimum structure

#### Required `obs` fields
- `cell_id`
- `sample_id`
- `time_label`

#### Required `var` fields
- `gene_symbol`

#### Required matrix
- gene expression matrix `X`

### Recommended optional fields
- dataset-level metadata
- batch metadata
- donor metadata
- stage metadata
- preprocessing metadata

At the current stage, the framework should **not hardwire repository-specific custom label files into the core execution path**. The core benchmark should first run on minimal required inputs plus whatever method-specific metadata is strictly necessary. Existing custom label assets may be reconsidered later only after the core evaluation pipeline is stable.

---

## 6. Unified method capability declaration

To keep the codebase rigorous and prevent invalid evaluations, each method adapter must explicitly declare what it supports.

Each method must expose the following capability flags:

- `supports_unseen_timepoint_projection`
- `supports_lineage_inference`

### Interpretation

#### If `supports_unseen_timepoint_projection = True`
The method is eligible for:
- Forecast Accuracy
- Embedding Coherence
- Lineage Fidelity

#### If `supports_unseen_timepoint_projection = False`
The method is **not** eligible for:
- Forecast Accuracy
- Embedding Coherence

It may still be eligible for:
- Lineage Fidelity

### Current method settings
- `WOT`:
  - `supports_unseen_timepoint_projection = False`
  - `supports_lineage_inference = True`
- `CellRank2`:
  - `supports_unseen_timepoint_projection = False`
  - `supports_lineage_inference = True`

This explicit gating replaces the previous bad idea of forcing OT methods into unsupported benchmark dimensions. It also ensures clean evaluator dispatch and avoids logic drift later.

---

## 7. Core task 1 — Forecast Accuracy

### 7.1 Definition
Forecast Accuracy evaluates how well a method projects cells from timepoint `t` to unseen timepoint `t+1` in gene expression space. This follows scTimeBench directly.

### 7.2 Eligibility
A method is evaluated on Forecast Accuracy **only if** it directly generates projected future cells or projected future gene expression profiles.

Therefore:
- eligible: future generative / forecasting models
- not eligible: WOT, CellRank2, and similar OT-style models that do not generate unseen future-timepoint cells directly

### 7.3 Metrics
Use the same metric family as scTimeBench:
- Wasserstein Distance
- Gaussian MMD
- Energy Distance MMD
- Hausdorff Loss

### 7.3.1 Metric implementation contract
Forecast Accuracy metrics must be computed by the unified benchmark evaluator, not by individual method wrappers. Method code may emit native diagnostic metrics for debugging, but those values must not be used for official ranking unless they are recomputed by the shared evaluator.

The official implementation follows the scTimeBench `OTLossMetric` family:
- Wasserstein Distance: `geomloss.SamplesLoss("sinkhorn", p=2, blur=0.05, scaling=0.5, debias=True, backend="tensorized")`, divided by the number of genes.
- Gaussian MMD: `geomloss.SamplesLoss("gaussian", blur=1.0, debias=True, backend="tensorized")`, divided by the number of genes.
- Energy Distance MMD: `geomloss.SamplesLoss("energy", blur=1.0, debias=True, backend="tensorized")`, divided by the number of genes.
- Hausdorff Loss: bidirectional nearest-neighbor Hausdorff distance computed with `torch.cdist`, without gene-count normalization by default.

For every eligible method, the evaluator must compute metrics per evaluation time point and aggregate by the mean across evaluation time points. Dispatcher logic must not bypass the unified evaluator simply because an adapter has already written `forecast_metrics.json`; adapter-produced values should be preserved only as native diagnostics, for example `method_native_forecast_metrics.json`.

This rule prevents metric drift across models. In particular, scNODE, PRESCIENT, and MIOFlow must be compared on the same expression arrays with the same metric definitions, normalization policy, sampling policy, and aggregation policy.

### 7.4 Aggregation
Within each scenario:
1. compute the four forecast metrics
2. rank methods within each metric
3. average the metric ranks to obtain the Forecast Accuracy rank

### 7.5 Current benchmark status
This task is currently **inactive** for the WOT vs CellRank2 stage, because neither method qualifies for this benchmark dimension under the scTimeBench rule. It will be activated when additional models with true unseen-timepoint projection capability are added.

---

## 8. Core task 2 — Embedding Coherence

### 8.1 Definition
Embedding Coherence evaluates whether projected cells preserve biologically meaningful structure in embedding space after forward temporal projection. This follows scTimeBench directly.

### 8.2 Eligibility
A method is evaluated on Embedding Coherence **only if** it produces projected future cells or projected future expression profiles at unseen time points.

Therefore:
- eligible: future generative / forecasting models
- not eligible: WOT, CellRank2, and similar OT-style methods without direct unseen-timepoint cell generation

### 8.3 Metrics
Use the same two metric families as scTimeBench:
- Adjusted Rand Index (ARI)
- average normalized classifier entropy

The metric definitions remain scTimeBench-style, but the label system used by the official evaluator is updated from `scgpt_pseudostate_provisional` to milestone labels.

### 8.3.1 Active official-silver labels

Embedding Coherence is currently computed under one active primary annotation
mode:

```text
label_mode: official_silver
state_key: final_milestone_label_coarse
role: primary report
```

For projection-capable methods, the evaluator must derive clusters from the
method output itself. The active official-silver implementation builds a kNN
graph from `projected_embedding.npy`, performs Leiden clustering, and compares
those unsupervised clusters against the frozen reference labels in
`final_milestone_label_coarse`.

Projected milestone labels are not used as clustering input. Earlier
`consensus`, `embedding_based`, and `classifier_based` projected-cell
annotation modes are superseded design iterations and are not part of the
current official-silver execution path.

### 8.3.2 Official-silver diagnostics

Diagnostics may be reported alongside, but not mixed into, the official
Embedding Coherence rank:

- ambiguous / low-confidence fraction in the official reference labels
- excluded-label fraction for `ambiguous` and `unknown_or_ood`
- label distribution by time point and scenario
- Leiden cluster count and cluster-size balance from projected embeddings
- confidence distribution from the official-silver annotation workflow

These diagnostics answer whether Embedding Coherence conclusions are stable
given the frozen `official_silver` reference labels.

### 8.4 Aggregation
Within each scenario:
1. rank methods on ARI
2. rank methods on entropy
3. average the ranks to obtain the Embedding Coherence rank

Rank aggregation must be performed within each `(dataset, scenario,
result_class, ground_truth_provider, label_mode)` group. The active primary
report uses `label_mode: official_silver` and
`state_key: final_milestone_label_coarse`.

### 8.5 Current benchmark status
Embedding Coherence is active for projection-capable official-silver runs from
scNODE, PRESCIENT, and MIOFlow. It remains inactive for WOT, CellRank2, and
other methods that do not generate projected future cells.

For the current projection-capable formal report, Embedding Coherence ARI values
are low in absolute magnitude and should not be overinterpreted without a null.
Until a within-dataset random-projection or label-permutation null is available,
ARI should be framed as a relative coherence diagnostic against the frozen
`official_silver` milestone system, not as an absolute biological-validity
score.

---

## 9. Core task 3 — Lineage Fidelity

### 9.1 Definition
Lineage Fidelity evaluates whether a method can recover a cell-state lineage structure that matches a benchmark reference lineage. This is the only active core dimension for the current WOT vs CellRank2 stage. It also directly follows scTimeBench.

### 9.2 Eligibility
A method is eligible for Lineage Fidelity if it can produce:
- cell-to-cell transition probabilities
- cell-state transition probabilities
- lineage graph structure
- or equivalent forward transition information

Current eligible methods:
- WOT
- CellRank2
- scNODE
- PRESCIENT
- MIOFlow

### 9.3 Reference requirement
A formal Lineage Fidelity benchmark requires:
1. a defined cell-state system
2. a benchmark reference lineage graph
3. a method for assigning or aggregating predicted next-state transitions

At the present stage, these should be treated as **benchmark prerequisites** rather than assumed repository defaults. Existing project-side draft label files are **not yet automatically wired into the active benchmark flow**. They may be adopted later only after the core evaluator logic is stabilized and the reference definition is explicitly frozen.

### 9.3.1 Ground-truth provider strategy
The benchmark should treat the definition of "ground truth" as a **replaceable provider layer**, not as a hard-coded property of any single annotation method.

For Forecast Accuracy, the ground truth is comparatively direct:
- observed held-out cells at the target time point
- observed gene expression used as the reference distribution

For Embedding Coherence and Lineage Fidelity, the ground truth is more model-dependent because it requires:
- a cell-state annotation system
- per-cell state labels
- per-state metadata
- a reference lineage graph over those states

Therefore, each annotation / reference system should be represented as a modular **ground-truth provider**. A provider must produce the following standard assets:

```text
state_labels.tsv
state_metadata.tsv
reference_graph.json
reference_graph_edges.csv
ground_truth_metadata.json
```

Each provider must also declare:
- `provider_id`
- `annotation_method`
- `state_key`
- `label_path`
- `metadata_path`
- `reference_graph_path`
- `reference_edges_path`
- `confidence_mode`
- `exclude_uncertain_states`
- `status`
- provenance notes

### 9.3.2 Active official-silver provider strategy

The active annotation-dependent benchmark uses a frozen **official silver**
milestone provider layer. It is not the older scGPT pseudostate graph and it no
longer builds the former `consensus`, `embedding_based`, and
`classifier_based` provider triplet for the primary benchmark.

For GSE178325, the active primary milestone graph is:

```text
hADSCs -> epithelial_like -> intermediate_plastic -> xen_like -> hCiPS
```

For GSE230659, the active primary milestone graph is:

```text
hADSCs -> epithelial_like -> intermediate_plastic -> hCiPS
```

The active registered official-silver providers for these two datasets are:

```text
gse178325_marker_fm_transition_silver_v1
gse230659_marker_fm_transition_silver_v1
```

Both providers use:

```text
label_mode: official_silver
state_key: final_milestone_label_coarse
label_type: frozen_silver_standard
excluded official-metric labels: ambiguous, unknown_or_ood
```

The per-cell labels are assigned by the actually executed two-stage workflow:

1. `build_marker_seed_labels.py` performs Stage 1 annotation. It assigns
   initial high-confidence milestone labels using available sample/time
   provenance, manual or literature-supported marker rules, and marker-score
   evidence. Cells that are not confidently assigned in this stage are routed
   to `ambiguous`.
2. `build_trajectory_aware_labels.py` performs Stage 2 annotation. It revisits
   Stage 1 `ambiguous` cells using marker-score strength, score margins,
   adjacent-transition logic, and optional embedding-centroid evidence. Cells
   can be resolved to coarse milestones, adjacent transition states,
   `unknown_or_ood`, or remain `ambiguous`.
3. `build_milestone_providers.py` freezes the official provider assets using
   `final_milestone_label_coarse` as the state label. The current implementation
   only builds `official_silver` providers; it does not generate the older
   `consensus`, `embedding_based`, or `classifier_based` provider modes.

The active GSE178325 / GSE230659 silver-standard construction is:

```mermaid
flowchart TD
    A["Observed cells in GSE178325 / GSE230659"] --> B["Stage 1: build_marker_seed_labels.py"]
    B --> C["High-confidence sample/time, manual/literature marker, or marker-score milestone"]
    B --> D["Uncertain cells: ambiguous"]
    C --> E["Seed-labeled cells"]
    D --> F["Stage 2: build_trajectory_aware_labels.py"]
    F --> G["Marker scores + adjacent-transition logic"]
    F --> H["Optional embedding-centroid evidence"]
    G --> I["Resolved milestone / transition / unknown_or_ood"]
    H --> I
    F --> J["Still unresolved: ambiguous"]
    E --> K["final_milestone_label_coarse"]
    I --> K
    J --> K
    K --> L["build_milestone_providers.py"]
    L --> M["Frozen official_silver provider"]
```

GSE242424 uses a separate active provider,
`gse242424_oskm_reprogramming_silver_v1`, based on author-cluster matched OSKM
reprogramming states. It is also registered as `official_silver` with
`state_key: final_milestone_label_coarse`, but it is not produced by the
GSE178325/GSE230659 marker-FM transition workflow above.

### 9.3.3 Historical provider designs

Earlier framework drafts used `scgpt_v1`,
`consensus_milestone_label`, `milestone_embedding_label`, and
`milestone_classifier_label` as candidate silver-standard routes. Those designs
are retained only as chronological context in the experimental record. They are
not active provider-generation instructions for the current benchmark.

When the active provider is selected in runtime configs or summaries, the
benchmark should use the registry-backed provider metadata:

```yaml
ground_truth:
  label_mode: official_silver
  state_key: final_milestone_label_coarse
  exclude_uncertain_states: true
```

All summary and ranking tables must therefore carry the ground-truth provider
identity and `label_mode: official_silver`. Rankings should be computed within
the same `(dataset, scenario, result_class, ground_truth_provider, label_mode)`
group, so that scores from different datasets or provider versions are not
mixed as if they shared the same reference scale.

### 9.4 Predicted lineage construction
For each eligible method:
- infer transitions across benchmark time intervals
- aggregate these transitions to the state level
- construct a predicted state-transition matrix
- convert the matrix into a predicted lineage graph using a fixed benchmark rule

Under the active official-silver plan, this construction is performed with one
primary state key:

```text
official_silver:
  cell_state_key = final_milestone_label_coarse
```

The official state-transition matrix for manuscript-facing Lineage Fidelity is
the `official_silver` milestone transition matrix. Historical
`consensus`/`embedding_based`/`classifier_based` matrices are not part of the
current primary benchmark.

### 9.5 Metrics
Use the same metric family as scTimeBench:
- AUROC
- AUPRC
- Jaccard Similarity
- single-step lineage recovery
- multi-step lineage recovery

Metric definitions are unchanged from the scTimeBench-style lineage evaluator. The change is the biological state system: official positives now come from the marker-defined milestone graph rather than from the legacy scGPT pseudostate graph.

### 9.6 Baseline
Formal Lineage Fidelity must include a **scTimeBench-style correlation baseline**, not only a model-vs-model comparison.
The baseline should be implemented as a direct reproduction of the `Correlation` method in the scTimeBench source code (`methods/correlation/run.py` with `configs/correlation/spearman_max.yaml`), so that baseline interpretation is aligned with the reference benchmark rather than being a project-specific approximation.

Required formal implementation:

1. Treat the baseline as a method-like lineage predictor that consumes the same scenario-filtered AnnData used by the evaluated method.
2. Iterate over adjacent available training timepoint pairs `(t, t_next)` in sorted order; do not use held-out target timepoints unless they are part of that scenario's training/evaluation cell universe by design.
3. For each adjacent pair, extract cell-level expression matrices `X_t` and `X_t_next`.
4. Use **Spearman rank correlation** by ranking each cell's gene-expression vector across genes (`rankdata(..., axis=1, method="average")`) before z-scoring.
5. Compute all cell-to-cell correlations between `t` and `t_next` as `corr = (z_t @ z_t_next.T) / n_features`, replacing non-finite values with zero and clipping scores to `[-1, 1]`.
6. For each source cell at `t`, score every target state present at `t_next` by the **maximum** correlation to cells in that target state (`averaging_method: maximum`, matching `spearman_max.yaml`).
7. Assign the source cell one vote for the target state with the largest score, accumulating votes into a source-state by target-state matrix.
8. After all adjacent timepoint pairs are processed, row-normalize the vote matrix so each represented source-state row sums to 1; source states with no valid source cells should remain all-zero rows.
9. Evaluate this weighted baseline state-transition matrix against the same provider reference graph and with the same Lineage Fidelity metrics used for WOT, CellRank2, and scNODE.
10. For thresholded graph metrics, follow the scTimeBench principle of deriving a predicted graph from the weighted matrix by an automatic PR/precision+recall threshold, while retaining AUROC/AUPRC on the weighted matrix.

This replaces the current state-level Pearson mean-expression baseline for formal reporting. The older baseline may remain only as a diagnostic legacy comparison and must be labeled explicitly as `state_mean_pearson_baseline`, not as the formal scTimeBench baseline.

### 9.6.1 Baseline interpretation

The scTimeBench-style Spearman correlation baseline is not merely a sanity check in GSE230659. In Scenarios A and C it outperforms several trajectory methods on Lineage AUROC, which means the benchmark reference graph is strongly recoverable from static gene-expression similarity alone. This pattern should be reported as a primary result, not buried as a control.

The interpretation is that Lineage Fidelity on these datasets measures alignment with a coarse marker-defined milestone graph under strong time-library confounding. A trajectory model that does not beat the correlation baseline may still have useful generative or dynamic properties, but it has not demonstrated additional lineage-recovery value beyond expression-state similarity under this reference. Scenario B remains especially difficult because early-only training leaves later reprogramming milestones unobserved and the single-library-per-timepoint design makes extrapolation vulnerable to library-specific shifts.

**Deprecated note.** Earlier versions described this as alignment with a coarse scGPT-derived state graph. That interpretation applies only to legacy scGPT-v1 provider runs.

### 9.7 Aggregation
Within each scenario:
1. compute AUROC
2. compute AUPRC
3. compute Jaccard Similarity
4. compute single-step recovery
5. compute multi-step recovery
6. rank methods within each metric
7. average the metric ranks to obtain the Lineage Fidelity rank

Rank aggregation must be performed within each `(dataset, scenario,
result_class, ground_truth_provider, label_mode)` group. The active primary
Lineage Fidelity rank uses `label_mode: official_silver` and
`state_key: final_milestone_label_coarse`.

### 9.8 Current benchmark status
Lineage Fidelity is active for all lineage-capable methods. WOT and CellRank2 remain Lineage-Fidelity-only methods because they do not generate projected future cells. scNODE, PRESCIENT, and MIOFlow are lineage-capable and also eligible for Forecast Accuracy and Embedding Coherence.

---

## 10. Method-specific evaluation rules

### 10.1 WOT
WOT is treated as an OT-based lineage inference method.

It is:
- excluded from Forecast Accuracy
- excluded from Embedding Coherence
- included in Lineage Fidelity only

### 10.2 CellRank2
CellRank2 is treated the same way at the current benchmark stage.

It is:
- excluded from Forecast Accuracy
- excluded from Embedding Coherence
- included in Lineage Fidelity only

The current CellRank2 adapter uses WOT transport maps through `RealTimeKernel.from_wot()`. The result should therefore be described as CellRank2 fate/lineage post-processing on WOT-derived transport maps, not as an independent velocity-kernel CellRank2 benchmark. A native CellRank2 evaluation with an RNA-velocity or other independent kernel can be added as a separate method configuration if the required inputs are available.

### 10.3 Future generative / forecasting models
If a future model can directly generate unseen timepoint cells or projected future expression profiles, then it will be evaluated on:
- Forecast Accuracy
- Embedding Coherence
- Lineage Fidelity

### 10.4 Future OT-like methods
If a future method is similar to WOT and still cannot generate unseen future cells directly, it should again follow the same rule:
- Lineage Fidelity only

No exception layer should be created just to force benchmark symmetry.

### 10.5 Projection-method run limitations

Full formal PRESCIENT A-C runs completed as CPU formal runs because CUDA was unavailable on the HPC execution path. These results are valid outputs of the recorded protocol, but should be labeled CPU full-formal results. Poor PRESCIENT performance, especially in Scenario B, should not be overinterpreted as a definitive method failure until GPU-enabled repeat runs or seed replicates are available.

All stochastic projection-capable methods require uncertainty estimates before final submission. The minimum acceptable plan is either at least three independent seeds per method-scenario combination or a bootstrap over held-out cells for the reported metrics. Until then, small rank differences should be treated as descriptive rather than statistically resolved.

---

## 11. Ranking strategy

Use **rank aggregation**, following the scTimeBench spirit of ranking within metrics and then aggregating.

### 11.1 Within each active core dimension
- rank methods within each metric
- average ranks across the metrics belonging to that dimension

### 11.2 Scenario-level reporting
For each scenario, report:
- Forecast Accuracy rank, if applicable
- Embedding Coherence rank, if applicable
- Lineage Fidelity rank, if applicable

### 11.3 Overall reporting rule
Overall reporting must respect method eligibility.

#### For methods that support all three dimensions
Report:
- three task ranks
- combined overall core rank

#### For methods that support only Lineage Fidelity
Report:
- Lineage Fidelity rank only

Do **not** fabricate a combined overall score for methods that were not eligible for the first two dimensions. In other words, WOT and CellRank2 should not receive a fake three-task overall rank under the current benchmark stage.

This point is essential for fairness and for preventing misleading comparisons.

---

## 12. Output specification

### 12.1 Forecast Accuracy outputs
Only for eligible methods:
- `projected_expression.npy`
- `forecast_metrics.json`
- `per_timepoint_forecast_metrics.csv`

`forecast_metrics.json` and `per_timepoint_forecast_metrics.csv` are official benchmark outputs only when written by the unified evaluator described in Section 7.3.1. If a method wrapper writes its own forecast diagnostics, those files must be treated as method-native diagnostics and preserved separately, for example:
- `method_native_forecast_metrics.json`
- `method_native_per_timepoint_forecast_metrics.csv`

### 12.2 Embedding Coherence outputs
Only for eligible methods:
- `projected_embedding.npy`
- `embedding_metrics.json`
- `projected_cluster_labels.csv`

Official-silver metadata should record:
- `label_mode: official_silver`
- `state_key: final_milestone_label_coarse`
- `cluster_source: leiden_from_projected_embedding`
- excluded-label handling for `ambiguous` and `unknown_or_ood`

The older label-mode-specific filenames
`embedding_metrics_consensus.json`, `embedding_metrics_embedding_based.json`,
and `embedding_metrics_classifier_based.json` are historical names and are not
required for the active official-silver benchmark.

### 12.3 Lineage Fidelity outputs
For all lineage-eligible methods:
- `state_transition_matrix.csv`
- `lineage_graph_edges.csv`
- `lineage_metrics.json`

Official-silver metadata should record:
- `label_mode: official_silver`
- `state_key: final_milestone_label_coarse`
- provider ID and reference graph path
- excluded-label handling for `ambiguous` and `unknown_or_ood`

The older consensus/embedding/classifier suffix-specific lineage outputs are
historical names and are not required for the active official-silver benchmark.

### 12.4 Run metadata
For all methods:
- `run_metadata.json`

This file should include:
- method name
- dataset
- scenario
- ground-truth provider identity
- ground-truth provider status
- annotation label mode (`official_silver`)
- state key used for annotation / aggregation
- reference graph path and confidence mode
- capability flags
- benchmark dimensions actually executed
- runtime
- status
- notes

---

## 13. Repository structure

```text
benchmark/
├─ README.md
├─ configs/
│  ├─ benchmark_master.yaml
│  ├─ scenario_observed.yaml
│  ├─ scenario_pseudotime.yaml
│  └─ method_capabilities.yaml
├─ adapters/
│  ├─ base_adapter.py
│  ├─ wot_adapter.py
│  ├─ cellrank2_adapter.py
│  └─ future_model_adapter.py
├─ evaluation/
│  ├─ eval_forecast.py
│  ├─ eval_embedding.py
│  ├─ eval_lineage.py
│  └─ eval_dispatch.py
├─ reports/
├─ results/
└─ docs/
```

The ground-truth provider layer is now represented explicitly under:

```text
benchmark/ground_truth/
  registry.yaml
  loader.py
  providers/
    gse178325_marker_fm_transition_silver_v1/
      state_labels.tsv
      state_metadata.tsv
      annotation_votes.tsv
      reference_graph.json
      reference_graph_edges.csv
      ground_truth_metadata.json
    gse230659_marker_fm_transition_silver_v1/
    gse242424_oskm_reprogramming_silver_v1/
```

Repository-specific draft label assets remain outside the active benchmark path unless they are promoted into a registered ground-truth provider with frozen provenance and standard output files.

**Active official-silver annotation structure.** The tree above is retained as
the original framework scaffold. The current executable benchmark uses the
following annotation layer and frozen providers:

```text
benchmark/annotation/
  build_marker_seed_labels.py
  build_trajectory_aware_labels.py
  build_milestone_providers.py

benchmark/ground_truth/providers/
  gse178325_marker_fm_transition_silver_v1/
    state_labels.tsv
    state_metadata.tsv
    annotation_votes.tsv
    reference_graph.json
    reference_graph_edges.csv
    ground_truth_metadata.json
  gse230659_marker_fm_transition_silver_v1/
  gse242424_oskm_reprogramming_silver_v1/
```

New milestone label assets should first be generated under
`benchmark/annotation/` and then frozen under `benchmark/ground_truth/providers/`
with standard provider metadata. The active provider builder accepts
`label_mode: official_silver` and freezes `final_milestone_label_coarse`; it no
longer creates the older consensus/embedding/classifier provider triplet.

---

## 14. Current implementation priority

The previous execution order was:

### Step 1
Implement the capability-gated benchmark dispatcher.

### Step 2
Finalize `eval_lineage.py` as the first active core evaluator.

### Step 3
Make sure WOT and CellRank2 can both enter Lineage Fidelity cleanly.

### Step 4
Keep `eval_forecast.py` and `eval_embedding.py` implemented as framework modules, but do not force-run them for WOT / CellRank2.

### Step 5
Only after future generative models are added, activate Forecast Accuracy and Embedding Coherence.

This order keeps the logic strict and avoids pretending that unsupported methods can do tasks they were never meant to do.

**Deprecated note.** The step list above is retained for the initial framework rebuild. It is no longer the complete implementation priority because projection-capable methods have already been added and the metrics plan has shifted to milestone annotation.

The updated execution order should be:

### Step 1
Define official-silver marker-supported milestones for each dataset.

- GSE178325: `hADSCs`, `epithelial_like`, `intermediate_plastic`, `xen_like`, `hCiPS`
- GSE230659: `hADSCs`, `epithelial_like`, `intermediate_plastic`, `hCiPS`

### Step 2
Run `build_marker_seed_labels.py` to build Stage 1 labels from sample/time
provenance, manual or literature-supported marker rules, and marker-score
diagnostics. Cells that are not confidently assigned remain `ambiguous`.

### Step 3
Run `build_trajectory_aware_labels.py` to resolve Stage 1 `ambiguous` cells
with marker-score strength, adjacent-transition logic, and optional
embedding-centroid evidence.

### Step 4
Freeze the final annotation columns:

- `final_milestone_label_coarse`
- `final_milestone_label_expanded`
- `final_milestone_confidence`
- `final_milestone_source`

### Step 5
Run `build_milestone_providers.py` to freeze official-silver ground-truth
providers:

- `gse178325_marker_fm_transition_silver_v1`
- `gse230659_marker_fm_transition_silver_v1`
- `gse242424_oskm_reprogramming_silver_v1`

### Step 6
Run Embedding Coherence and Lineage Fidelity under `label_mode:
official_silver`, while keeping Forecast Accuracy unchanged and label-free.

### Step 7
Regenerate reports so that every annotation-dependent metric row carries
`provider_id`, `label_mode`, `state_key`, and `reference_graph_path`.

---

## 15. What is intentionally not used now

The following repository-specific assets are **temporarily outside the active core benchmark**:

- `celltype_labels_v1.tsv`
- `pseudotime_labels_v1.tsv`
- `reference_graph_v1.json`
- `terminal_labels_v1.tsv`
- `terminal_labels_v1_1.tsv`
- previous project-specific T1 / T2 / T3 / T4 / T5 task structure

These files may later be reviewed and reintroduced as:
- supplementary descriptions
- dataset-specific annotation assets
- sensitivity analyses
- auxiliary biological interpretation layers

But they should not be allowed to complicate the first pass of the scTimeBench-aligned core benchmark.

**Deprecated note.** This section describes the earlier framework-stabilization
phase. Under the active official-silver plan, marker-defined milestone labels
are the primary state system for annotation-dependent metrics, but the former
consensus-provider and scGPT-v1 pseudostate systems are historical designs
rather than current executable providers.

---

## 16. Summary

This framework is now aligned with scTimeBench at the level of **core evaluation logic**:

- **Forecast Accuracy**
- **Embedding Coherence**
- **Lineage Fidelity**

At the same time, it follows the same method-eligibility rule as scTimeBench:

- methods that can generate unseen future cells: evaluate all three dimensions
- methods that cannot generate unseen future cells, such as WOT-like methods: evaluate **Lineage Fidelity only**

Under this rule, **WOT** and **CellRank2** remain Lineage-Fidelity-only methods,
while **scNODE**, **PRESCIENT**, and **MIOFlow** are official
projection-capable methods eligible for all three dimensions. The current
annotation-dependent benchmark uses registry-backed `official_silver`
providers and `final_milestone_label_coarse`; earlier `scgpt_v1` runs are
historical rather than active primary evidence.

The updated primary metrics plan is:

```text
Forecast Accuracy:
  label-free scTimeBench expression-distribution metrics

Embedding Coherence:
  scTimeBench ARI / classifier entropy
  projected embeddings clustered with Leiden
  clusters compared against final_milestone_label_coarse

Lineage Fidelity:
  scTimeBench lineage metrics
  computed on milestone transition matrices against marker-defined milestone graphs
  aggregated by final_milestone_label_coarse under official_silver
```

The primary annotation-dependent state system is now
`final_milestone_label_coarse` under `label_mode: official_silver`. The older
`consensus_milestone_label`, `milestone_embedding_label`, and
`milestone_classifier_label` systems are historical design iterations and are
not active primary providers.

PRESCIENT first completed A/B/C CPU low-memory reduced-validation runs, then completed full formal A/B/C HVG2000 runs without per-timepoint training subsampling. MIOFlow full formal A/B/C HVG2000 runs are included in the primary projection-capable table.

No OT projection workaround is introduced. No repository-specific label asset is forcibly embedded into the active core pipeline. Instead, annotation-derived labels and lineage references are handled through an explicit ground-truth provider layer. The current benchmark is runnable for the formal A/B/C observed-time scenarios, with formal reporting separated from smoke tests, HPC validation runs, pilot backups, and reduced-validation runs by explicit `result_class` metadata.

The experimental record below is retained as a reproducibility log. Manuscript-facing text should use the stable framework sections above; the chronological record should be moved to a supplementary appendix or `CHANGELOG.md` before submission.

---

## 17. Experimental record

Entries before the 2026.5.11 official-silver completion record are historical.
When those entries use words such as "current" or "active", they refer to the
state of the project on that date, not to the current executable benchmark.

**2026-05-08 - Metrics plan update to consensus milestone annotation**
The metrics plan was revised while retaining the scTimeBench top-level metric families. Forecast Accuracy remains label-free and keeps the centralized scTimeBench expression-distribution metrics. Embedding Coherence and Lineage Fidelity now use marker-defined milestone labels as the primary biological state system. Two milestone annotation model families are required: `embedding_based` annotation, where scGPT or similar embeddings are followed by a milestone classifier trained from marker-defined seed cells, and `classifier_based` annotation, where CellTypist, scANVI, TOSICA, SingleR, or scArches-style label transfer predicts milestone labels. The `consensus` milestone provider is the primary reporting mode, while `embedding_based` and `classifier_based` providers are required sensitivity modes. The older scGPT-v1 pseudostate provider is deprecated as the primary reference and retained only for reproducibility and legacy sensitivity analysis.

---

**2026-04-18 - Repository rebuild, WOT pipeline, and scGPT integration**
Repo rebuilt from the earlier project-specific T1-T5 layout into a scTimeBench-aligned benchmark structure with configs, adapters, evaluators, dispatcher, reports, and results directories. The preprocessing path produced `adata_benchmark.h5ad` (75,194 x 2,000) with the required benchmark fields. WOT execution was implemented and debugged against the correct `OTModel -> compute_all_transport_maps -> from_directory -> get_coupling` API. scGPT was selected as the zero-shot embedding source for a working state system, using pretrained whole-human scGPT embeddings to define provisional pseudostates for Lineage Fidelity.

---

**2026-04-19 - scGPT-v1 state system and CellRank2 integration**
The `scgpt_v1` silver-standard provider was established with 14 pseudostates (`PS_00`-`PS_13`) and a frozen reference graph with 42 directed edges (3 high, 28 medium, 11 low). The active formal setting uses `edge_confidence_mode: medium_and_above`, giving 31 reference edges. WOT and CellRank2 configs were wired to `abs_day`, `scgpt_pseudostate_provisional`, and the scGPT-v1 reference. CellRank2 was implemented as a real adapter using WOT transport maps and `RealTimeKernel.from_wot()`, followed by sparse state-level aggregation.

---

**2026-04-20 - Scenario B baseline leak diagnosis**
A baseline leak was found and fixed: the dispatcher had passed the full AnnData object into `run_lineage_evaluation()` instead of the adapter's scenario-filtered training universe. After the fix, the correlation baseline is computed on the same cell universe used by the method. `summarize_lineage.py` was updated to rank only the five framework Lineage Fidelity metrics, keep `jaccard_similarity_topk` diagnostic-only, and collapse the correlation baseline to one row per `(scenario, result_class, ground_truth_provider)`.

---

**2026-04-21 - Scenario C split and memory-efficient runs**
Scenario C was defined as observed-time interpolation plus extrapolation: 10 training time points (0.5, 2.0, 8.0, 16.0, 16.33, 16.67, 17.0, 18.0, 22.0, 24.0), interpolation holdouts (4.0, 12.0, 20.0), and extrapolation holdouts (28.0, 30.0). Backed AnnData loading and post-filter materialization were added for large Scenario B/C jobs. Pilot C runs confirmed the split and baseline logic before full-data execution.

---

**2026-04-22 - Performance optimization and HPC validation**
WOT state aggregation was vectorized from a Python double loop to sparse `S_src @ M @ S_tgt.T`, `--skip-tmap-if-exists` was added for transport-map reuse, CellRank2 memory use was reduced by avoiding an extra full AnnData copy, and hard-coded local paths were replaced with `TRAJ_PROJECT_ROOT` / `SCGPT_REPO` overrides plus project-root discovery. Shirokane validation reproduced WOT/A, CellRank2/A, CellRank2/B, and full CellRank2/C. Full WOT/C also completed with AUROC about 0.817. Scenario B consistently remained below baseline for WOT and CellRank2, supporting the interpretation that early-only training is insufficient for full lineage recovery.

---

**2026-04-23 - scNODE adapter integration**
scNODE was added as the first projection-capable benchmark method. `scnode_adapter.py` was registered in `eval_dispatch.py`, `method_capabilities.yaml` was updated, and scNODE was marked eligible for Forecast Accuracy, Embedding Coherence, and Lineage Fidelity. A memory-feasible HVG2000 route was introduced after the naive full-gene full-cell run exceeded memory. HVG selection was made robust by replacing Scanpy's failing Seurat-style routine with a finite sparse mean/variance dispersion selector.

---

**2026-04-24 - scNODE reduced and formal benchmark runs**
scNODE Lineage Fidelity was stabilized by replacing naive soft-assignment normalization with a stable softmax, adding deterministic nearest-centroid fallback, zero-row handling for unsupported source states, and `lineage_diagnostics.json`. Reduced A/B/C validation runs completed, followed by formal A/B/C HVG2000 runs under `scgpt_v1` with `result_class: official`, `formal_benchmark: true`, `pretrain_iters=200`, `epochs=10`, `iters=100`, `batch_size=32`, `latent_dim=50`, and `n_sim_cells=2000`. The scNODE formal runs produced Forecast Accuracy, Embedding Coherence, and Lineage Fidelity outputs for all three observed-time scenarios.

---

**2026-04-25 - Formal A/B/C benchmark synchronization**
Formal WOT, CellRank2, and scNODE results were synchronized into the local result tree. CellRank2 Scenario C was promoted from pilot backup to a full-data official result. The formal lineage comparison now contains WOT, CellRank2, and scNODE for Scenarios A, B, and C under the fixed `scgpt_v1` provider. WOT and CellRank2 remain Lineage-Fidelity-only by capability, while scNODE carries all three core dimensions.

---

**2026-04-28 - Result-class cleanup and report regeneration**
All run directories with metadata were assigned explicit `result_class` values. Formal runs are marked `official`; smoke tests, Shirokane validation runs, pilot backups, and reduced-validation runs are marked as non-official classes (`smoke`, `hpc_validation`, `pilot_backup`, `reduced_validation`). `summarize_core_results.py` and `summarize_lineage.py` now recognize these classes, so `--official-only` excludes all diagnostic runs reliably. The report assets were regenerated:

- `benchmark/reports/core_summary.csv` contains 18 traceability rows across all result classes.
- `benchmark/reports/lineage_summary.csv` contains 33 rows grouped by result class.
- `benchmark/reports/lineage_summary_official.csv` contains 15 official rows only.
- `benchmark/reports/formal_benchmark_summary.csv` contains the 9 formal model-scenario rows.
- `benchmark/reports/figures/lineage_auroc_by_model_scenario.png` and `benchmark/reports/scnode_formal_metrics_table.md` were rebuilt from the cleaned formal summary.

Current formal result snapshot:

- Scenario A Lineage AUROC: WOT 0.8559, CellRank2 0.7636, scNODE 0.7611.
- Scenario B Lineage AUROC: WOT 0.5341, CellRank2 0.5631, scNODE 0.5015.
- Scenario C Lineage AUROC: WOT 0.8174, CellRank2 0.8194, scNODE 0.7161.

**2026-04-29 - PRESCIENT integration and neutral HVG benchmark input**
The shared HVG2000 input was moved out of the scNODE result tree and renamed as a model-neutral benchmark input: `benchmark/inputs/gse230659_scgpt_hvg2000/GSE230659_scGPT_annotated_HVG2000_benchmark_input.h5ad`. Its metadata was changed from `scnode_hvg_*` to `benchmark_hvg_*`, with `benchmark_input_label` set to `GSE230659 scGPT-annotated HVG2000 benchmark input`. scNODE and PRESCIENT scripts/configs now consume this neutral input instead of a scNODE-owned path.

PRESCIENT was added as a separate benchmark method with the original author repository vendored under `benchmark/methods/PRESCIENT/prescient_module`, matching the existing scNODE pattern (`benchmark/methods/scNODE/scNODE_module`). The PRESCIENT wrapper now defaults to this project-local author-code copy, while `PRESCIENT_REPO` remains available as an explicit override. The previous accidental fallback to scNODE's internal PRESCIENT baseline copy was removed from the normal search path.

PRESCIENT A/B/C CPU low-memory runs completed on Shirokane with `result_class: reduced_validation` and `formal_benchmark: false`. These runs used capped training/reference cells and are treated as integration validation rather than official formal benchmark evidence. Outputs were synchronized locally, a malformed trailing-NUL metadata file for Scenario B was cleaned, and report summaries were regenerated. PRESCIENT reduced-validation metrics now appear in `benchmark/reports/core_summary.csv`, `benchmark/reports/lineage_summary.csv`, and `benchmark/reports/prescient_reduced_validation_metrics_table.md`, while `lineage_summary_official.csv` and `formal_benchmark_summary.csv` correctly exclude them.

PRESCIENT reduced-validation snapshot:

- Scenario A: Forecast WD 3.3461, Embedding ARI 0.1194, Lineage AUROC 0.6610.
- Scenario B: Forecast WD 3.6448, Embedding ARI 0.1276, Lineage AUROC 0.5591.
- Scenario C: Forecast WD 4.4809, Embedding ARI 0.1594, Lineage AUROC 0.6490.

**2026-04-30 - Full formal PRESCIENT completion**
Full formal PRESCIENT A/B/C HVG2000 runs completed on Shirokane with `result_class: official`, `formal_benchmark: true`, and `training_protocol: full_formal_hvg2000`. The runs used no per-timepoint training subsampling and no embedding reference-cell cap. The training universes were A: 75,194 cells, B: 27,891 cells, and C: 49,820 cells. The qsub logs show `s_vmem=128G`, internal-cache cleanup before execution, `missing: []`, and successful job completion. CUDA was unavailable in these runs, so the completed full formal PRESCIENT results are CPU full formal results.

Formal summaries were regenerated after synchronization. PRESCIENT is now included in `benchmark/reports/core_summary.csv`, `benchmark/reports/lineage_summary.csv`, `benchmark/reports/lineage_summary_official.csv`, `benchmark/reports/formal_benchmark_summary.csv`, `benchmark/reports/figures/lineage_auroc_by_model_scenario.png`, and `benchmark/reports/projection_formal_metrics_table.md`.

PRESCIENT full formal snapshot:

- Scenario A: Forecast WD 3.1054, Gaussian MMD 0.1316, Embedding ARI 0.1449, Lineage AUROC 0.6766.
- Scenario B: Forecast WD 4.0936, Gaussian MMD 0.1842, Embedding ARI 0.2345, Lineage AUROC 0.5470.
- Scenario C: Forecast WD 4.4904, Gaussian MMD 0.2081, Embedding ARI 0.2881, Lineage AUROC 0.6524.

The current next scientific step is interpretation rather than pipeline construction: compare official projection-capable methods (scNODE vs PRESCIENT), keep WOT and CellRank2 lineage-only, and explain why Scenario B remains difficult across methods.

**2026.4.30 - MIOFlow integration and formal benchmark completion**
MIOFlow was added as a projection-capable benchmark method using the project-local source under `benchmark/methods/MIOFlow`, with PCA-space training/inverse-PCA expression reconstruction kept consistent with scTimeBench rather than applying nonnegative clipping.
CPU reduced A/B/C validation runs completed first, then full formal A/B/C HVG2000 runs completed on Shirokane and were synchronized locally with all required metric/result files present.
Formal MIOFlow snapshot: A WD 5.2884, ARI 0.0585, Lineage AUROC 0.7589; B WD 8.5655, ARI 0.1896, Lineage AUROC 0.4990; C WD 4.5306, ARI 0.1097, Lineage AUROC 0.7720.
Scenario B remains the main failure case, showing residual negative extrapolation drift and near-random lineage AUROC despite improvement over the reduced validation run.
The MIOFlow integration, formal lightweight `json/csv` metrics, and ignore rules for local heavy artifacts were committed in `df5f6e8`, `e873ed2`, and `ef628f1`.

**2026-05-02 - GSE178325 metric-definition drift diagnosis**
During the first GSE178325 0618-only projection benchmark, scNODE showed extremely large Forecast WD while PRESCIENT and MIOFlow showed small WD. Direct range checks showed this was not caused by scNODE expression-scale explosion. The real issue was evaluator drift: scNODE reported its author-code GeomLoss value without gene-count normalization, while PRESCIENT and MIOFlow wrappers reported project-local per-gene/scipy-style metrics. This violates the scTimeBench framework contract because models were not being compared through one evaluator.

The framework rule is now explicit: projection-capable methods output standardized predicted expression arrays, and official Forecast Accuracy metrics are recomputed centrally with the scTimeBench `OTLossMetric` definitions. Method-native metrics may be retained only as diagnostics and must not drive benchmark interpretation or ranking.

**2026.5.03 - GSE178325 scGPT-aligned benchmark completion**
GSE178325 0618-only processing was realigned with the GSE230659/scTimeBench-style route: build raw/full-gene input first, derive scGPT embeddings/pseudostates on full genes, then create the shared HVG2000 benchmark input.
Duplicate gene symbols in GSE178325 are now merged by summed raw expression before normalization/log1p so that scGPT embedding uses unique gene-symbol features.
The active provider is `scgpt_v1_gse178325_0618`, with `scgpt_pseudostate_provisional` as the state key and a frozen reference lineage graph used for Lineage Fidelity.
scNODE, PRESCIENT, and MIOFlow A/B/C runs were rerun on the aligned input; all 9 result directories completed with Forecast Accuracy, Embedding Coherence, and Lineage Fidelity metrics.
Forecast metrics are centrally recomputed with the `scTimeBench_geomloss` backend, while embedding and lineage metrics use the scGPT-derived provider rather than the earlier stage-proxy labels.
Aggregate ranking and the GSE178325 projection benchmark report were regenerated; current aggregate leaders are scNODE for Scenarios A/B and MIOFlow for Scenario C.
Local/Shirokane scripts were updated to use the aligned input/provider paths, H100 scGPT provider construction where required, and clean output handling before reruns.
scTimeBench code review found fixed seeds/random states for reproducibility, but no standard multi-random-seed repeated-run protocol; the current single-seed design is therefore framework-consistent.

---

**2026-05-03 - Unified GSE230659 forecast re-evaluation and reviewer-response updates**
All GSE230659 projection-capable formal forecast metrics were recomputed from existing `projected_expression.npy` outputs using the centralized `scTimeBench_geomloss` evaluator in `benchmark/evaluation/eval_forecast.py`. This replaced method-native wrapper metrics in `forecast_metrics.json` and preserved older wrapper outputs in `method_native_forecast_metrics.json` / `method_native_per_timepoint_forecast_metrics.csv`.

Unified Forecast WD snapshot: scNODE A/B/C 81.7124, 414.6872, 266.7125; PRESCIENT A/B/C 114.6700, 153.1575, 224.9176; MIOFlow A/B/C 163.3884, 406.9134, 224.4649. `benchmark/reports/formal_benchmark_summary.csv` and `benchmark/reports/projection_formal_metrics_table.md` were regenerated/updated to include MIOFlow and the unified forecast metrics.

Reviewer-facing limitations were added for the GSE230659 time-library confound, scGPT-v1 silver-standard circularity, low ARI interpretability without a null, the strong Spearman-correlation baseline, WOT-backed CellRank2, CPU-only PRESCIENT formal runs, inactive pseudotime Scenarios D-F, and the need for seed/bootstrap uncertainty estimates.

---

**2026.5.04 - scTimeBench-style dataset/preprocessor split refactor**
The trajectory split logic formerly duplicated across method YAML/adapters was lifted into `benchmark/shared/dataset/preprocessors/ScenarioTimepointSplit`, with scTimeBench-style `BaseDataset`, preprocessor registry, factory helpers, and GSE230659/GSE178325 dataset registries added under `benchmark/shared/dataset`.
Existing Scenario A/B/C semantics were preserved: `train_times` and `heldout_times` remain the source of truth, while adapters now share the same split helper instead of reimplementing timepoint filtering.
Shirokane validation completed for GSE230659 WOT Scenario B smoke/fullcheck, scNODE/PRESCIENT/MIOFlow Scenario B, and scNODE/PRESCIENT/MIOFlow Scenario C; all completed with expected heldout forecast timepoints and scenario-specific lineage baselines.
GSE178325 scNODE/PRESCIENT/MIOFlow A/B/C reruns also completed after the refactor; Scenario A/B/C baseline pair counts matched the intended train-time universes (14, 7, and 6 respectively).
The refactor therefore preserves capability-gated evaluation and scGPT provider lineage metrics while aligning the data split layer more closely with scTimeBench.
Minor follow-up: MIOFlow metadata should explicitly record `train_times`/`heldout_times`, although its forecast and lineage outputs already confirm the split was applied correctly.

---

**2026.5.09 - Milestone-based metrics refactor and real sensitivity provider transition**
The benchmark metrics were refactored away from the earlier scGPT-only pseudostate silver standard toward marker-defined milestone labels while preserving the three scTimeBench-level evaluation families: Forecast Accuracy, Embedding Coherence, and Lineage Fidelity. Forecast Accuracy remains annotation-independent and continues to use the centralized expression-distribution evaluator. The annotation-dependent metric families now use milestone state systems with explicit provider IDs, label modes, state keys, and reference graphs. The active milestone graphs are `epithelial_like -> intermediate_plastic -> hCiPS` for GSE230659 and `epithelial_like -> intermediate_plastic -> xen_like -> hCiPS` for GSE178325.

Full-gene milestone annotation was completed for both datasets. GSE178325 uses full-gene milestone labels directly. GSE230659 uses full-gene marker scoring first, then transfers the resulting labels back to the HVG2000 benchmark h5ad by cell identity so that model training remains HVG2000-based while annotation uses full-gene-derived biological labels. The GSE230659 full-gene-derived label distribution before the real-sensitivity transition was epithelial-like dominant with `epithelial_like`, `intermediate_plastic`, `hCiPS`, and `ambiguous` states. The corresponding milestone providers were created under `benchmark/ground_truth/providers/` and registered as `gse230659_milestone_consensus_v1`, `gse230659_milestone_embedding_v1`, and `gse230659_milestone_classifier_v1`, with GSE178325 analogues also populated.

The legacy scGPT providers were explicitly marked as deprecated legacy pseudostate references rather than primary benchmark references. Registry and metadata fields now distinguish `label_mode`, `analysis_role`, `deprecated`, and `legacy` status. Summary scripts enforce that `legacy_scgpt_pseudostate` rows cannot become primary results. The preferred primary reporting mode is `consensus`; `embedding_based` and `classifier_based` are sensitivity modes.

The GSE230659 consensus milestone formal benchmark was executed on Shirokane using `run_milestone_gse230659_primary_array.sh`. The completed formal runs were scNODE A/B/C plus WOT A and CellRank2 A. All five runs wrote lineage metrics using `provider_id=gse230659_milestone_consensus_v1`, `label_mode=consensus`, and `cell_state_key=consensus_milestone_label`. WOT and CellRank2 remain Lineage-Fidelity-only by capability because they do not produce projected expression/embedding outputs. scNODE additionally produced Forecast Accuracy outputs for A/B/C. The milestone lineage summary after this run contained five consensus primary lineage rows and one legacy scGPT comparison row. The GSE230659 consensus lineage snapshot from that run was approximately: CellRank2 A AUROC 0.0714/AUPRC 0.1024, WOT A AUROC 0.0714/AUPRC 0.1024, scNODE A AUROC 0.2857/AUPRC 0.1270, scNODE B AUROC 0.6429/AUPRC 0.2667, and scNODE C AUROC 0.3571/AUPRC 0.1500. These values became historical after the later real-sensitivity provider switch and should be rerun before final interpretation.

Projected-cell milestone annotation and Embedding Coherence were then completed for the scNODE GSE230659 formal runs using the then-current projected marker-scoring route. This paragraph is now a **historical implementation note**, superseded by the 2026.5.10 projected-cell provider policy below. `run_milestone_gse230659_embedding_coherence_array.sh` was added to annotate projected cells and run `eval_embedding_milestone.py` for scNODE A/B/C. WOT and CellRank2 tasks intentionally write skip metadata because their result directories do not contain `projected_expression.npy`, `projected_embedding.npy`, and `projected_cluster_labels.csv`. An initial bug in projected marker scoring was fixed by correctly unpacking `resolve_milestone_genes()` into `milestone_gene_indices` and marker overlap metadata. The projected-label validator was also updated so formal projected annotation is allowed; it now rejects only contradictory metadata where `smoke_test=true` and `formal_benchmark=true` simultaneously.

Embedding Coherence was revised to follow the scTimeBench entropy logic more closely. Hard per-cell milestone labels are still used for ARI and lineage transition counting, but the primary entropy metric is now a per-cell prediction/probability entropy over milestone probabilities rather than hard-label cluster entropy. In the superseded 2026.5.09 projected-marker implementation, `annotate_projected_cells.py` wrote probability columns such as `prob_epithelial_like`, `prob_intermediate_plastic`, and `prob_hCiPS` from softmax-transformed marker scores. Under the 2026.5.10 policy, formal projected-cell entropy should instead use provider probabilities such as `embedding_prob_*`, `classifier_prob_*`, or calibrated `consensus_prob_*`. `eval_embedding_milestone.py` reads milestone probability columns and reports `mean_prediction_entropy`, `weighted_prediction_entropy`, `entropy_basis=milestone_probability_vector`, and `n_probability_classes`. The old hard-label cluster entropy is retained only as a diagnostic and may be `null` with the note `entropy undefined: only one label present` when projected hard labels collapse to a single milestone.

The report layer was updated to carry the new entropy fields. `embedding_summary.csv` now includes `mean_prediction_entropy`, `weighted_prediction_entropy`, `entropy_basis`, `n_probability_classes`, `prediction_entropy_note`, and `hard_label_entropy_note`. `core_summary.csv`, `lineage_summary.csv`, `embedding_summary.csv`, and `provider_agreement_summary.csv` regenerate with provider-aware primary/sensitivity/legacy flags. A reporting hygiene issue where CellRank2 core rows had `dataset_id=unknown` was fixed by inferring dataset ID from provider ID, source h5ad, or result path when run metadata is incomplete.

The sensitivity annotation providers were then upgraded from placeholder copies to real providers for GSE230659. `build_embedding_milestone_model.py` now trains a real `X_scGPT -> logistic regression -> milestone_embedding_label` model and writes `milestone_embedding_confidence` plus `embedding_prob_*` columns. `build_classifier_milestone_model.py` now trains a CellTypist-backed milestone classifier and writes `milestone_classifier_label`, `milestone_classifier_confidence`, and `classifier_prob_*` columns, with an expression-classifier fallback available. CellTypist on Shirokane required installing the `celltypist` package into the active `traj_env`; the required dependency set for this route is `anndata`, `numpy`, `pandas`, `scipy`, `scikit-learn`, `joblib`, `scanpy`, `celltypist`, `h5py`, and `pyyaml`. The wrapper prepares an in-memory CellTypist-compatible normalized/log1p copy so the benchmark h5ad does not need to be rewritten just to satisfy CellTypist's input-scale check.

`run_gse230659_real_sensitivity_annotation.sh` was added and completed on Shirokane. It generated `benchmark/inputs/gse230659_milestone_real_sensitivity_hvg2000/GSE230659_real_sensitivity_milestone_HVG2000_benchmark_input.h5ad` with 75,194 cells and 2,000 genes. The new real-sensitivity labels are no longer marker-copy placeholders. The synchronized label distributions were: marker-defined labels, `epithelial_like` 62,989, `ambiguous` 5,259, `hCiPS` 4,696, `intermediate_plastic` 2,250; embedding-based labels, `epithelial_like` 45,517, `intermediate_plastic` 19,020, `hCiPS` 10,657; CellTypist classifier labels, `epithelial_like` 61,156, `hCiPS` 7,298, `intermediate_plastic` 6,740. The final consensus label is now a three-provider majority vote over marker, scGPT-embedding, and CellTypist labels: `epithelial_like` 60,847, `hCiPS` 6,991, `intermediate_plastic` 5,475, and `ambiguous` 1,881. The consensus metadata reported `agreement_fraction=0.7427`, `three_provider_majority=73313`, and `discordant_abstain=1881`.

GSE230659 was then formally switched to this real-sensitivity h5ad. `build_milestone_providers.py` rebuilt the three GSE230659 milestone providers from `GSE230659_real_sensitivity_milestone_HVG2000_benchmark_input.h5ad`, and all seven GSE230659 milestone configs were patched so both `dataset.h5ad_path` and `state_system.source_h5ad` point to the new real-sensitivity input. `validate_milestone_configs.py` passed for all seven GSE230659 milestone configs after the switch. GSE178325 providers and configs were intentionally left unchanged. The current active GSE230659 providers therefore now use real sensitivity labels, while previously synchronized GSE230659 formal milestone results should be treated as historical and rerun before final interpretation under the active provider set.

Current next steps: rerun the GSE230659 formal milestone benchmark under the active real-sensitivity providers, beginning with `run_milestone_gse230659_primary_array.sh`, then rerun `run_milestone_gse230659_embedding_coherence_array.sh` for scNODE A/B/C, and finally run the scNODE sensitivity array to compare consensus, embedding-based, and classifier-based providers. After GSE230659 is regenerated, the same real-sensitivity strategy can be extended to GSE178325 once an embedding source or alternative representation is available for its embedding-based provider.

---

**2026.5.10 - Projected-cell provider policy revision**
After comparing the current projected-cell annotation behavior with the scTimeBench implementation, the projected-cell provider design was narrowed to two formal automatic annotation routes: `embedding_based` and `classifier_based`. The observed/reference h5ad may still use marker-defined labels as seed labels, biological diagnostics, and consensus-support evidence. However, projected cells should no longer treat marker scoring as a formal third annotation provider.

The formal projected-cell output should therefore contain `milestone_embedding_label`, `milestone_classifier_label`, `consensus_milestone_label`, `consensus_confidence`, and provider probability columns such as `embedding_prob_*`, `classifier_prob_*`, and optionally `consensus_prob_*`. `milestone_marker_label` is no longer required for projected cells. If it appears, it must be explicitly marked diagnostic/deprecated and must not be copied into `milestone_embedding_label` or `milestone_classifier_label`.

This revision prevents a misleading sensitivity analysis in which `consensus`, `embedding_based`, and `classifier_based` are identical only because all three were copied from projected marker scores. Future formal projected-cell metrics should annotate projected embeddings with the scGPT/embedding milestone classifier, annotate projected expression with CellTypist or an equivalent classifier/reference-transfer tool, and then form projected consensus labels from those two automatic providers.

Projected-cell formal annotation was then fixed to an HVG2000-consistent route. The full-gene scGPT embedding classifier is retained only for observed-cell provenance/sensitivity and must not annotate HVG2000-only projected expression. A new `hvg2000_pca_logistic_regression` embedding provider was trained from observed HVG2000 expression and marker seed labels, saved under `benchmark/annotation_runs/gse230659_real_sensitivity/embedding_model_pca_hvg2000/`, and paired with the existing CellTypist classifier provider. `annotate_projected_cells.py` now derives `projected_hvg_embedding.npy` from `projected_expression.npy`, writes sidecar metadata, and treats `projected_embedding.npy` as method-native latent space only. scNODE GSE230659 A/B/C projected labels and Embedding Coherence were regenerated and validated with `formal_two_provider_complete=true`. Current projected audit: A provider match 0.814, ambiguous 18.6%, consensus ARI 0.033077; B provider match 0.658, ambiguous 34.2%, consensus ARI 0.000935; C provider match 0.962, ambiguous 3.8%, consensus ARI 0.003935. The low ARI values should be interpreted as projected-label coherence diagnostics, not overall model performance; Scenario C has the most stable annotation agreement, while Scenario B has the largest provider disagreement. Audit outputs are stored in `benchmark/reports/official/gse230659_projected_hvg_annotation_audit.*`.

---

**2026.5.11 - Official silver trajectory-aware milestone benchmark completion**
The milestone annotation layer was finalized as an `official_silver` state system to make the trajectory benchmark follow the scTimeBench-style `embedding -> Leiden clusters -> ARI against reference labels` logic while avoiding reliance on unavailable native cell-type annotations. The primary label source is now `final_milestone_label_coarse`, produced by a two-stage marker/trajectory-aware annotation workflow. Stage 1 assigns high-confidence marker or sample/time labels, including `hADSCs` as the starting milestone where supported by the source dataset. Stage 2 resolves remaining ambiguous cells with trajectory-aware marker/foundation-model evidence into coarse milestones, adjacent transition labels, `unknown_or_ood`, or `ambiguous`. Diagnostic milestone and `failure_or_offtrajectory` routing were removed from the active benchmark because they were not part of the agreed label-construction plan.

`milestone_markers.yaml` was updated so the active primary milestone orders are `hADSCs -> epithelial_like -> intermediate_plastic -> xen_like -> hCiPS` for GSE178325 and `hADSCs -> epithelial_like -> intermediate_plastic -> hCiPS` for GSE230659. The excluded official-metric labels are now only `ambiguous` and `unknown_or_ood`. The annotation scripts were revised accordingly: `build_marker_seed_labels.py` assigns Stage 1 labels and provenance, `build_trajectory_aware_labels.py` produces final coarse/expanded labels and confidence fields, and `build_milestone_providers.py` exports the `gse178325_marker_fm_transition_silver_v1` and `gse230659_marker_fm_transition_silver_v1` providers.

Full annotation jobs completed on Shirokane for both datasets and were synchronized locally. GSE178325 produced a trajectory-dominant official silver label distribution of approximately 67,817 `epithelial_like`, 5,706 `hADSCs`, 4,654 `intermediate_plastic`, 2,144 `hCiPS`, and 154 `xen_like` cells. GSE230659 produced approximately 67,112 `epithelial_like`, 5,370 `hCiPS`, and 2,712 `intermediate_plastic` cells; no `hADSCs` were assigned there because the available earliest time point was 0.5 rather than 0. All synchronized annotation outputs contained no `failure_or_offtrajectory`, `out_of_trajectory`, or `marker_diagnostic` labels.

The official silver formal benchmark was rerun for scNODE, MIOFlow, and
PRESCIENT across GSE178325/GSE230659 Scenarios A/B/C, plus WOT A and CellRank2
A for GSE230659 as lineage-only baselines. The original primary array launcher
used for that historical rerun is no longer part of the active root script set;
the current project tracks these outputs through the result manifest and keeps
only the active official-silver validation and embedding-coherence rerun scripts
at the repository root.

Embedding Coherence was corrected to match the intended scTimeBench-style implementation. The first official silver embedding-coherence run was invalid because `eval_embedding_milestone.py` read `projected_cluster_labels.csv`, which had already been copied from milestone labels and therefore caused a self-comparison artifact with ARI equal to 1.0. The evaluator now ignores those projected milestone labels for clustering and instead builds a kNN graph directly from `projected_embedding.npy`, performs Leiden clustering with `n_neighbors=15` and `resolution=1.0`, and compares those unsupervised clusters against `final_milestone_label_coarse`. After rerunning `run_marker_fm_transition_silver_embedding_coherence_array.sh`, all 18 projection-capable official silver embedding results used `cluster_source=leiden_from_projected_embedding:projected_embedding.npy:n_neighbors=15:resolution=1.0` with nontrivial ARI values ranging from about 0.0093 to 0.2031.

The reporting layer now treats `official_silver` as the active primary label
mode for annotation-dependent metrics. The active root validation script is
`run_marker_fm_transition_silver_validate_configs_array.sh`; the older
consensus/scGPT summary policy and the historical combined
summarize-and-validate launcher are not part of the current executable root
script set.

The current official silver reports are regenerated under
`benchmark/reports/official_silver/` with explicit official-silver filenames:
`official_silver_embedding_method_summary.csv`,
`official_silver_embedding_rankings.csv`,
`official_silver_lineage_method_summary.csv`,
`official_silver_lineage_rankings.csv`,
`official_silver_combined_method_summary_scnode_mioflow_prescient.csv`, and
`official_silver_model_rankings.md`. Under the current capability rules, all
three benchmark metric families are complete for scNODE, MIOFlow, and PRESCIENT
across both GSE178325/GSE230659 datasets and Scenarios A/B/C. WOT and CellRank2
remain complete only for their intended lineage-only role and are not expected
to produce forecast or embedding-coherence metrics.

Official silver ranking tables were generated in `benchmark/reports/official_silver/`. The combined embedding-plus-lineage ranking over the three projection-capable methods is very close: PRESCIENT ranks first with combined rank score 1.7333, MIOFlow second with 1.7389, and scNODE third with 1.7611. This difference is small and should be interpreted cautiously. By metric family, scNODE ranks best for Embedding Coherence because it has the lowest average rank after combining ARI and entropy, while MIOFlow ranks best for Lineage Fidelity. PRESCIENT has the highest mean ARI among the three projection-capable methods and the best combined rank by a narrow margin. WOT and CellRank2 are reported only in the lineage summary because they lack projected expression/embedding outputs.

At this point the official silver benchmark loop is closed: annotation providers, formal model runs, corrected Embedding Coherence, lineage metrics, forecast metrics, validation scripts, summary tables, and model rankings have all been generated and synchronized locally. The next scientific step is interpretation and manuscript/report writing rather than another pipeline rebuild, unless additional seed replicates, bootstrap uncertainty, or alternative annotation-provider sensitivity analyses are requested.

---

**2026.5.14 - GSE242424 input integration**
GSE242424 raw scRNA matrices were audited with `GSE242423_scRNA_genes.tsv`; all 9 samples matched the 36,601-gene feature table and yielded 156,969 QC-passing cells.
Reusable MTX streaming utilities and `build_gse242424_raw_full_gene_input.py` produced the full-gene h5ad, followed by `build_gse242424_hvg2000_input.py` producing the 156,969 x 2,000 HVG2000 benchmark input with 50D PCA.
GSE242424 was registered in the dataset factory with Scenario A/B/C configs; `scTimeBench_cell_type` is currently a time/stage proxy and `iPSC abs_day=16` remains provisional.
