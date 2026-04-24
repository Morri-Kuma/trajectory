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
Include a correlation-based baseline analogous to scTimeBench’s lineage baseline, so that lineage reconstruction quality is not judged only relative to other complex methods.

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

### 10.3 Future generative / forecasting models
If a future model can directly generate unseen timepoint cells or projected future expression profiles, then it will be evaluated on:
- Forecast Accuracy
- Embedding Coherence
- Lineage Fidelity

### 10.4 Future OT-like methods
If a future method is similar to WOT and still cannot generate unseen future cells directly, it should again follow the same rule:
- Lineage Fidelity only

No exception layer should be created just to force benchmark symmetry.

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

Under this rule, the current first-stage benchmark for **WOT vs CellRank2** activates only **Lineage Fidelity**. Forecast Accuracy and Embedding Coherence remain in the framework, but they are deferred until future models that actually support unseen-timepoint generation are added.

No OT projection workaround is introduced. No repository-specific label asset is forcibly embedded into the active core pipeline at this stage. Instead, annotation-derived labels and lineage references are handled through an explicit ground-truth provider layer. The first priority is to make the scTimeBench-aligned core benchmark logic clean, strict, modular, and runnable.

---

## 17. Experimental record

**2026-04-18 — Repository rebuild, WOT pipeline, and scGPT integration**
Repo rebuilt from scratch (v1→v2): legacy T1–T5 structure deleted; benchmark scaffold created (configs, adapters, evaluators, dispatcher). Preprocessing script updated (v3→v4) to produce `adata_benchmark.h5ad` (75,194 × 2,000; all required obs/var/obsm/uns fields validated). `GSE230659Dataset` pkl class, builder script, and `wot_gse230659_observed.yaml` config created. WOT `run.py` written; three bugs patched: GBK encoding on Windows, wrong `TransportMapModel.compute()` API (replaced with correct `OTModel → compute_all_transport_maps → from_directory → get_coupling` sequence), wrong coupling accessor (`.x` → `.X`). scGPT role defined as zero-shot embedding tool: embed all cells into 512-dim space → cluster → freeze pseudo-states for Lineage Fidelity. Pretrained whole-human model validated (`vocab.json`, `args.json`, `best_model.pt` 205 MB). Benchmark h5ad chosen for smoke test (92.3% vocab overlap); raw h5ad for production (84.8%, 27,267 genes, raw counts). Smoke-test script written at `benchmark/scgpt/smoke_test.py` (1,000-cell balanced subset → `embed_data` → UMAP/Leiden → save). Windows patch applied to `scgpt/tasks/cell_emb.py`: `os.sched_getaffinity` guarded with `hasattr` (`num_workers=0` on Windows). Status: WOT and scGPT smoke tests both pending local execution.

---

**2026-04-19 — scGPT v1 state system integration and CellRank2 adapter**
scGPT v1 silver-standard state system established: 14 pseudostates (PS_00–PS_13) from Leiden clustering on `X_scGPT` (512-dim); reference graph frozen with 42 edges at three confidence tiers (3 high / 28 medium / 11 low). Parallel YAML configs created for WOT and CellRank2 (scgpt_v1 paths; `abs_day` time key; `scgpt_pseudostate_provisional` state key; `edge_confidence_mode: medium_and_above`, 31 edges). `eval_lineage.py`: reference graph loader implemented with confidence-tier filtering; honest `has_predictions` status field replacing always-"completed" placeholder. `eval_dispatch.py`: `--method-config` CLI flag added; `exclude_uncertain_states` wired end-to-end; scenario config (time key, state key, cellrank2 params) injected into adapter; UTF-8 encoding fix. WOT `run.py`: `day_field` hardcoding fixed (`time_label` → configurable); WOT confirmed producing non-empty outputs with scGPT-v1 config. CellRank2 adapter fully implemented (replacing scaffold): WOT transport → `RealTimeKernel.from_wot(path=...)` → state-level aggregation via sparse S@T@S.T. Three CellRank2 API errors patched: `from_wot` takes a directory path (not object); time key must be categorical; categorical must be sorted numerically (not lexicographically). Status: CellRank2 end-to-end run pending.

---

**2026-04-20 — Scenario B baseline leak diagnosis and cleanup**
Baseline leak confirmed: `eval_dispatch.py` passed full adata (75,194 cells) to `run_lineage_evaluation` instead of the adapter's scenario-filtered view, so the Scenario B baseline incorrectly included held-out timepoints; all pilot-pass baseline CSVs were byte-identical across scenarios. Fix applied in dispatcher: `eval_adata = getattr(adapter, "adata", adata)` — baseline now computed on the same cell universe the method trained on; Scenario A unchanged (AUROC 0.6743 verified). `cellrank2_adapter.py` Step 0 now filters to `train_times` from config and reassigns `self.adata` (75,194 → 27,891 cells for Scenario B). `scripts/recompute_scenario_b_baseline.py` written and applied to both WOT-B and CR2-B (`subsample=0`); baseline CSVs are now byte-identical between methods (MD5 verified). WOT-B JSON files: NUL-byte padding stripped (Windows sync artefact); `result_class: "pilot"` and `sampling` block added to `run_metadata.json`. `summarize_lineage.py` rewritten: ranks on five framework metrics only (`jaccard_topk` diagnostic-only), one baseline per `(scenario, result_class)` with `baseline_consistent` flag, `--official-only` flag. Policy established: the correlation baseline is a property of the scenario split, not the method's subsample. Status: Scenario C ready to run.

---

**2026-04-21 — Scenario C (interpolation + extrapolation) implemented and pilot-run**
Scenario C split defined: 10 training time points (days 0.5, 2.0, 8.0, 16.0, 16.33, 16.67, 17.0, 18.0, 22.0, 24.0; 49,820 cells); three interpolation holdouts (days 4.0, 12.0, 20.0); two extrapolation holdouts (days 28.0, 30.0). YAML configs created for WOT (`wot_gse230659_observed_scgpt_v1_scenarioC.yaml`) and CellRank2 (`cellrank2_gse230659_observed_scgpt_v1_scenarioC.yaml`); A/B configs untouched. Three code fixes to support memory-efficient runs: backed-mode loading added to `eval_dispatch.py` and `run.py` (load with `backed='r'`; materialize after scenario filter); `cellrank2_adapter.py` Step 0 fixed to use `.to_memory()` on backed slice; shutil.rmtree replaced with graceful per-file fallback (mounted-fs permission issue). Both methods run as pilot (300 cells/timepoint, 3,000 cells): WOT-C AUROC 0.686 / AUPRC 0.336 / SSR 0.387; CellRank2-C AUROC 0.694 / AUPRC 0.343 / SSR 0.387; both exceed baseline (AUROC 0.675). Baseline confirmed scenario-specific (computed on 3,000 pilot training cells). `lineage_summary.csv` regenerated covering A/B (official) + C (pilot). Full-data Scenario C runs pending on local machine (`conda activate traj_env`; same YAML configs; no further code changes needed).

---

**2026-04-22 — Performance audit, code optimizations, WOT-C full run, and Shirokane HPC validation**
Full-data WOT Scenario C completed on local machine (4.5 h, AUROC 0.817 / AUPRC 0.403 / SSR 0.452), confirming both methods score ~0.818 AUROC on Scenario C with full training data. Performance audit conducted across all source files; four optimizations implemented: (1) `_aggregate_to_state_level()` in `WOT/run.py` vectorized — Python O(n²) double loop replaced with sparse S_src @ M @ S_tgt.T, eliminating 80 M+ iterations on the largest transport maps; (2) `--skip-tmap-if-exists` CLI flag added to `WOT/run.py` for tmap reuse across reruns; (3) redundant `adata.copy()` eliminated in `CellRank2Adapter._run_lineage_fidelity_impl()` — reuses existing copy for both WOT and CellRank2 stages, reducing peak RAM from 3× to 2× dataset size; (4) all hard-coded Windows paths (`C:\Users\37620\trajectory`, `C:\...\scGPT`) removed from 9 source files and replaced with `TRAJ_PROJECT_ROOT` / `SCGPT_REPO` env-var overrides plus upward-search fallback, making the codebase portable to Linux/HPC. SGE job scripts written (`run_wot_validate.sh`, `run_cellrank_validate.sh`, `run_cellrank_validate_array.sh`) and project deployed to Shirokane HPC (`/home/xzy0723/projects/trajectory`). Four Shirokane jobs completed successfully (Apr 22 21:14–21:55 JST): WOT/A reproduced exactly (AUROC 0.8559, bit-for-bit match; tmap cache hit on 2nd run reduced runtime from 676 s to 24 s); CellRank2/A reproduced (AUROC 0.7634, Δ ≤ 0.001 vs local); CellRank2/B reproduced (AUROC 0.5665); CellRank2/C run at full scale for the first time (49,820 cells, 376 s, AUROC 0.8194 vs pilot 0.694). One SGE job failed due to a config filename typo in the shell script (non-critical; corrected). Confirmed: Scenario B below baseline for both methods (WOT AUROC 0.534, CR2 0.566 vs baseline 0.678) — early-only training is insufficient for full lineage recovery, a stable scientific finding. Outstanding: WOT/B and WOT/C not yet re-run on Shirokane; `result_class` label missing from CR2/C Shirokane result; summary CSV needs regeneration to cover all 6 runs.
---

**2026-04-23 - scNODE adapter integration, HPC scripts, and full-cell HVG run**
Added the first scNODE benchmark path: created `scnode_adapter.py`, registered `scnode` in `eval_dispatch.py`, injected `scnode_params` / `scenario_params` / `dataset_id`, and updated `method_capabilities.yaml`.
Replaced the WOT scaffold entry by wiring `WOTAdapter` to the real `benchmark/methods/WOT/run.py` logic.
Wrote Shirokane qsub scripts for scNODE smoke and full runs; the smoke test on GSE230659 human data completed successfully with all required outputs.
The first naive full-cell run was killed by memory pressure, so an HVG2000 reduced-training route was designed and scripted.
The original HVG selection failed because Scanpy's Seurat-style HVG routine hit `inf` values; this was fixed by replacing it with a finite sparse mean/variance dispersion selector.
The rerun completed full-cell training plus forecast/embedding output generation on Shirokane, but lineage evaluation still produced NaNs in the state-transition matrix, so `lineage_metrics.json` is still missing.

---

**2026-04-24 - scNODE A/B/C reduced runs and ground-truth provider modularization**
Fixed scNODE Lineage Fidelity NaNs by replacing the naive soft assignment normalization with a stable softmax, adding zero-row handling for unsupported source states, and writing `lineage_diagnostics.json`.
Re-ran scNODE full-cell HVG2000 Scenario A on Shirokane using the cached trained model; `lineage_metrics.json` completed successfully and the state-transition matrix contains no non-finite values.
Added Shirokane array script support for scNODE Scenario B and Scenario C HVG2000 reduced runs; both completed with Forecast Accuracy and Lineage Fidelity outputs.
Implemented the first modular ground-truth provider layer under `benchmark/ground_truth/`, with `scgpt_v1` registered as the current silver-standard working provider.
The provider layer standardizes `state_labels.tsv`, `state_metadata.tsv`, `reference_graph.json`, `reference_graph_edges.csv`, and `ground_truth_metadata.json`, allowing future annotation systems such as CellTypist, scANVI, SingleR, marker-rule, or consensus providers to replace scGPT without changing model adapters.
`eval_dispatch.py`, `eval_lineage.py`, and `summarize_lineage.py` now carry provider metadata so results can be grouped and ranked by `(scenario, result_class, ground_truth_provider)`.
