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

### 4.2 Pseudotime scenarios

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

### 8.4 Aggregation
Within each scenario:
1. rank methods on ARI
2. rank methods on entropy
3. average the ranks to obtain the Embedding Coherence rank

### 8.5 Current benchmark status
This task is currently **inactive** for the WOT vs CellRank2 stage. It will be activated only after adding models that qualify for future-cell generation.

For the current projection-capable formal report, Embedding Coherence ARI values are low in absolute magnitude and should not be overinterpreted without a null. The required next reporting step is to add a within-dataset random-projection or label-permutation null for each scenario, then report method ARI relative to that null distribution. Until that null is available, ARI should be framed as a relative coherence diagnostic against the scGPT-v1 state system, not as an absolute biological-validity score.

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

The current provider is:

```text
provider_id: scgpt_v1
status: silver_standard_working
annotation_method: scGPT
state_key: scgpt_pseudostate_provisional
```

This provider is derived from:
1. pretrained whole-human scGPT embeddings (`X_scGPT`)
2. Leiden clustering into 14 pseudostates (`PS_00`-`PS_13`)
3. per-state biological review into confirmed / merge-review / uncertain states
4. consecutive-timepoint kNN transition counting
5. confidence-tier filtering of reference edges

The resulting scGPT-v1 reference graph has:
- 14 nodes
- 42 total directed edges
- 3 high-confidence edges
- 28 medium-confidence edges
- 11 low-confidence edges

The default active Lineage Fidelity setting uses:

```text
edge_confidence_mode: medium_and_above
```

This means the active reference contains 31 edges (high + medium), while low-confidence edges involving uncertain states are excluded from the main metric calculation.

This scGPT-v1 provider is a **silver-standard working reference**, not final biological ground truth. Its role is to make the benchmark executable, reproducible, and comparable while preserving the option to replace or compare annotation systems later.

Because the reference graph is derived from scGPT embeddings of the same GSE230659 cells used in the benchmark, the provider carries a representational circularity risk. This is acceptable only if it is reported as a silver-standard reference and not as independent biological truth. Any method that uses scGPT-like representations, directly or indirectly, must be interpreted with this caveat. Ranking sensitivity to alternative providers such as CellTypist, SingleR, scANVI, marker-rule labels, or a consensus graph is a planned validation experiment rather than an optional cosmetic extension.

Future providers may include, for example:
- CellTypist-derived labels and lineage graph
- scANVI-derived labels and lineage graph
- SingleR-derived labels and lineage graph
- marker-rule-derived labels and lineage graph
- ensemble or consensus annotation providers

When a new provider is added, the model adapters should not need to change. The benchmark should instead switch provider metadata and paths through a `ground_truth` configuration block or provider registry:

```yaml
ground_truth:
  provider_id: scgpt_v1
  state_key: scgpt_pseudostate_provisional
  confidence_mode: medium_and_above
  exclude_uncertain_states: false
```

This enables two complementary comparisons:

1. **Fixed ground truth, compare trajectory models**
   - Example: WOT vs CellRank2 vs scNODE under `scgpt_v1`

2. **Fixed model, compare ground-truth providers**
   - Example: scNODE under `scgpt_v1` vs `celltypist_v1` vs `scanvi_v1`

All summary and ranking tables must therefore carry the ground-truth provider identity. Rankings should be computed within the same `(scenario, result_class, ground_truth_provider)` group, so that scores from different annotation systems are not mixed as if they shared the same reference scale.

### 9.4 Predicted lineage construction
For each eligible method:
- infer transitions across benchmark time intervals
- aggregate these transitions to the state level
- construct a predicted state-transition matrix
- convert the matrix into a predicted lineage graph using a fixed benchmark rule

### 9.5 Metrics
Use the same metric family as scTimeBench:
- AUROC
- AUPRC
- Jaccard Similarity
- single-step lineage recovery
- multi-step lineage recovery

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

The interpretation is that Lineage Fidelity on this dataset measures alignment with a coarse scGPT-derived state graph under strong time-library confounding. A trajectory model that does not beat the correlation baseline may still have useful generative or dynamic properties, but it has not demonstrated additional lineage-recovery value beyond expression-state similarity under this reference. Scenario B remains especially difficult because early-only training leaves later reprogramming states unobserved and the single-library-per-timepoint design makes extrapolation vulnerable to library-specific shifts.

### 9.7 Aggregation
Within each scenario:
1. compute AUROC
2. compute AUPRC
3. compute Jaccard Similarity
4. compute single-step recovery
5. compute multi-step recovery
6. rank methods within each metric
7. average the metric ranks to obtain the Lineage Fidelity rank

### 9.8 Current benchmark status
This is the **only active core benchmark dimension** for the current WOT vs CellRank2 comparison.

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

### 12.3 Lineage Fidelity outputs
For all lineage-eligible methods:
- `state_transition_matrix.csv`
- `lineage_graph_edges.csv`
- `lineage_metrics.json`

### 12.4 Run metadata
For all methods:
- `run_metadata.json`

This file should include:
- method name
- dataset
- scenario
- ground-truth provider identity
- ground-truth provider status
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
    scgpt_v1/
      state_labels.tsv
      state_metadata.tsv
      reference_graph.json
      reference_graph_edges.csv
      ground_truth_metadata.json
```

Repository-specific draft label assets remain outside the active benchmark path unless they are promoted into a registered ground-truth provider with frozen provenance and standard output files.

---

## 14. Current implementation priority

The current execution order should be:

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

---

## 16. Summary

This framework is now aligned with scTimeBench at the level of **core evaluation logic**:

- **Forecast Accuracy**
- **Embedding Coherence**
- **Lineage Fidelity**

At the same time, it follows the same method-eligibility rule as scTimeBench:

- methods that can generate unseen future cells: evaluate all three dimensions
- methods that cannot generate unseen future cells, such as WOT-like methods: evaluate **Lineage Fidelity only**

Under this rule, **WOT** and **CellRank2** remain Lineage-Fidelity-only methods, while **scNODE**, **PRESCIENT**, and **MIOFlow** are official projection-capable methods evaluated across all three dimensions under the fixed `scgpt_v1` ground-truth provider. PRESCIENT first completed A/B/C CPU low-memory reduced-validation runs, then completed full formal A/B/C HVG2000 runs without per-timepoint training subsampling. MIOFlow full formal A/B/C HVG2000 runs are included in the primary projection-capable table.

No OT projection workaround is introduced. No repository-specific label asset is forcibly embedded into the active core pipeline. Instead, annotation-derived labels and lineage references are handled through an explicit ground-truth provider layer. The current benchmark is runnable for the formal A/B/C observed-time scenarios, with formal reporting separated from smoke tests, HPC validation runs, pilot backups, and reduced-validation runs by explicit `result_class` metadata.

The experimental record below is retained as a reproducibility log. Manuscript-facing text should use the stable framework sections above; the chronological record should be moved to a supplementary appendix or `CHANGELOG.md` before submission.

---

## 17. Experimental record

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
