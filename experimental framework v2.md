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

At the current stage, label directories and repository-specific draft label assets do not need to be inserted into the active benchmark path yet.

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

No OT projection workaround is introduced. No repository-specific label asset is forcibly embedded into the active core pipeline at this stage. The first priority is to make the scTimeBench-aligned core benchmark logic clean, strict, and runnable.

---

## 17. Experimental record

### Session log — 2026-04-18

**Objective:** Build the end-to-end pipeline from raw h5ad to a runnable WOT Lineage Fidelity smoke test, and eliminate all API errors blocking execution.

---

#### 17.1 Repository refactoring (v1 → v2)

The repository was rebuilt from scratch to match this framework. All legacy content was deleted: the five-task (T1–T5) project-specific structure, old evaluators, old label assets, old results and figures, reference PDFs, and HPC artifacts. The `data/` directory was preserved intact. New framework scaffold was created:

- `benchmark/configs/benchmark_master.yaml` — rewritten to three core dimensions
- `benchmark/configs/method_capabilities.yaml` — explicit capability flags for WOT and CellRank2
- `benchmark/configs/scenario_observed.yaml` and `scenario_pseudotime.yaml` — Scenarios A–F
- `benchmark/adapters/base_adapter.py` — abstract base with capability-gated dispatch
- `benchmark/adapters/wot_adapter.py`, `cellrank2_adapter.py`, `future_model_adapter.py`
- `benchmark/evaluation/eval_lineage.py` — active evaluator (AUROC, AUPRC, Jaccard, single-step, multi-step recovery)
- `benchmark/evaluation/eval_forecast.py`, `eval_embedding.py` — inactive stubs, framework-ready
- `benchmark/evaluation/eval_dispatch.py` — capability-gated dispatcher

---

#### 17.2 Preprocessing script refactor — `scripts/01_load_and_qc.py` (v3 → v4)

The preprocessing script was converted from a "npz/csv + reload snippet" workflow into a canonical single-step AnnData builder. Key changes:

- Added `_find_project_root()` with Windows hardcoded path fallback and upward-walk auto-detection
- Added `anndata` and `scipy.sparse` imports
- Output changed from npz/csv to `data/processed/adata_benchmark.h5ad` directly
- New required obs columns added during Pass 1 cell metadata accumulation: `cell_id`, `dataset_id`, `sample_id`, `stage`, `day_within_stage`, `abs_day`, `time_label`, `scTimeBench_timepoint`, `scTimeBench_cell_type`
- Cell-type proxy mapping added: StageI → `somatic`, StageII → `early_transition`, StageIII → `intermediate`, hCiPSCs → `pluripotent`
- var frame extended with `gene_symbol` and `is_hvg` fields
- New Section 9 constructs and writes `AnnData` with `X` (sparse CSR log-normalized HVG matrix), `obsm["X_pca"]`, and `uns` metadata
- New Section 12 reloads the h5ad and validates all required fields

Result confirmed: **75,194 cells × 2,000 genes**, all required obs/var/obsm/uns fields present, validation passed.

---

#### 17.3 Dataset pkl pipeline

Three new files created:

- **`benchmark/datasets/gse230659_dataset.py`** — `GSE230659Dataset` class. Pickle-safe by design: `__getstate__`/`__setstate__` store only paths, not AnnData objects. `load_data()` returns `(train_adata, test_adata)` where train = all 75,194 cells across all 15 timepoints (required by WOT for consecutive pair transport) and test = 9,142 cells at day 30 (terminal hCiPSC, scaffold split). Helper methods `time_points()` and `cell_counts_per_timepoint()` use backed read to avoid loading the full expression matrix.

- **`benchmark/datasets/build_gse230659_dataset_pkl.py`** — Builder script. Instantiates the dataset, smoke-tests `load_data()`, pickles with `HIGHEST_PROTOCOL`, and verifies round-trip. All assertions passed; output pkl is 307 bytes.

- **`benchmark/configs/wot_gse230659_observed.yaml`** — WOT run configuration for Scenario A. Specifies pkl path, h5ad fallback path, OT hyperparameters (`epsilon=0.05`, `lambda1=1.0`, `lambda2=50.0`, `local_pca=30`, `growth_iters=3`, `growth_rate_source="uniform"`), cell-state key (`scTimeBench_cell_type`), and output directory.

---

#### 17.4 WOT run script — `benchmark/methods/WOT/run.py`

Full pipeline script created: YAML config → dataset pkl → `load_data()` → WOT transport inference → state-level aggregation → eval_lineage dispatch → output files. Three bugs were encountered and patched sequentially.

**Bug 1 — YAML UTF-8/GBK encoding error (Windows)**
- Error: `UnicodeDecodeError: 'gbk' codec can't decode byte` during YAML load
- Fix: both `open(path)` calls in `_load_yaml()` changed to `open(path, "r", encoding="utf-8-sig")`

**Bug 2 — WOT API mismatch: `TransportMapModel.compute()` does not exist**
- Error: `type object 'TransportMapModel' has no attribute 'compute'`
- Investigation: inspected `wot.ot`, `wot.tmap`, and `wot.ot.OTModel` in the sandbox
- Correct WOT 1.0.x API: `wot.ot.OTModel(adata, day_field=..., growth_rate_field=..., ...)` → `compute_all_transport_maps(tmap_out=prefix, overwrite=True)` → `wot.tmap.TransportMapModel.from_directory(prefix)` → `get_coupling(t_src, t_tgt)`
- Fix applied to `_run_wot_transport()`: replaced the entire wrong `TransportMapModel.compute()` call with the correct three-step sequence

**Bug 3 — State aggregation using wrong coupling accessor and wrong attribute**
- Error (anticipated): `_aggregate_to_state_level()` was using `tmap_model.get_transport_map()` (does not exist) and extracting `.x` (wot.Dataset attribute, removed in v1)
- Fix: replaced with `get_coupling(t_src, t_tgt)` (returns AnnData), extracted `tmap.X` (AnnData expression matrix), and used `tmap.obs_names`/`tmap.var_names` for authoritative row/column cell ID ordering

---

#### 17.5 Current pipeline state as of session end

| Component | Status |
|---|---|
| `data/processed/adata_benchmark.h5ad` | ✓ Built and validated (75,194 × 2,000) |
| `benchmark/datasets/gse230659_observed.dataset.pkl` | ✓ Built and round-trip verified (307 bytes) |
| `benchmark/configs/wot_gse230659_observed.yaml` | ✓ Complete |
| `benchmark/methods/WOT/run.py` | ✓ All known API bugs patched |
| WOT transport smoke test | **Pending** — run on local machine with `traj_env` |
| Lineage metrics | Deferred — `reference_graph_path: null` in YAML; `eval_lineage.py` will record status = `"deferred"` until a frozen reference graph is provided |
| `benchmark/evaluation/eval_forecast.py` | Inactive stub |
| `benchmark/evaluation/eval_embedding.py` | Inactive stub |

**Next step:** run the smoke test locally:
```
conda activate traj_env
cd C:\Users\37620\trajectory
python benchmark\methods\WOT\run.py --config benchmark\configs\wot_gse230659_observed.yaml
```
Expected output under `benchmark/results/wot/scenario_A/`: `state_transition_matrix.csv`, `lineage_graph_edges.csv`, `run_metadata.json`. Lineage metrics will be deferred until the reference lineage graph is frozen (§9.3).

---

## 18. Experimental record — scGPT integration planning and smoke test

### Session log — 2026-04-18

**Objective:** Reason about how scGPT should be used in this project, inspect available input data, design the embedding and pseudo-state pipeline, write a first smoke-test script, and resolve the first Windows-specific runtime error encountered during execution.

---

#### 18.1 Methodological role of scGPT in this project

scGPT is being used as a **representation-learning and annotation-assistance tool**, not as a source of biological hard ground truth. The intended workflow is: embed all cells across all timepoints into a shared 512-dimensional space using the pretrained whole-human model, cluster globally in that space, define frozen pseudo-states (silver-standard labels), and use those labels as the cell-state system for downstream Lineage Fidelity evaluation. This is a supplementary annotation asset and does not modify the active benchmark path or WOT configuration.

---

#### 18.2 scGPT codebase study

The local scGPT repository at `C:\Users\37620\Documents\GitHub\scGPT` was studied in full before any implementation was proposed. Key findings:

**API:** The primary entry point for zero-shot cell embedding is `scg.tasks.embed_data(adata, model_dir, gene_col=..., batch_size=..., device=..., use_fast_transformer=...)`. This function reads all three required model files, filters genes to the vocabulary intersection, tokenizes expression values on-the-fly using 51-bin quantile discretization inside `DataCollator`, runs the 12-layer transformer encoder, extracts the `<cls>` position output as the cell embedding, L2-normalizes it, and stores the result in `adata.obsm["X_scGPT"]` (shape: n_cells × 512).

**Vocabulary:** The vocab is a JSON mapping of HGNC gene symbols to integer token IDs. Genes not present in the vocab are silently dropped before embedding. The `gene_col` parameter must point to a column (or `"index"`) in `adata.var` containing HGNC approved gene symbols — Ensembl IDs produce zero overlap.

**Model directory requirements:** Exactly three files must be present with these exact names: `vocab.json`, `args.json`, `best_model.pt`. These are hardcoded in `embed_data`.

**`args.json` keys used by the code:** `embsize`, `nheads`, `d_hid`, `nlayers`, `n_layers_cls`, `dropout`, `pad_token`, `pad_value`. The whole-human model uses `embsize=512`, `nlayers=12`, `nheads=8`.

**No fine-tuning or classification head needed:** Zero-shot clustering — `embed_data` then standard scanpy `neighbors`/`umap`/`leiden` — is the correct path for our goal.

---

#### 18.3 Pretrained model directory validation

`C:\Users\37620\trajectory\models\scgpt_whole_human` confirmed complete:

| File | Size | Content |
|---|---|---|
| `vocab.json` | 1.3 MB | 60,697 tokens: 60,694 HGNC gene symbols + `<pad>`, `<cls>`, `<eoc>` (special tokens already present) |
| `args.json` | 1.3 KB | Full model config: `embsize=512`, `nlayers=12`, `nheads=8`, `n_bins=51`, `max_seq_len=1200`, `pad_value=-2` |
| `best_model.pt` | 205 MB | PyTorch state dict for the whole-human foundation model |

This is the CellXGene May 2023 census pretrained model — trained on diverse human tissue data with `no_cls=true` (no supervised classification objective), appropriate for general-purpose embedding.

---

#### 18.4 Data input inspection

Two h5ad files were inspected. The file `human_germ.h5ad` referenced in the original plan was not found in any accessible directory and could not be assessed.

**`data/processed/adata_benchmark.h5ad`**
- Shape: 75,194 cells × 2,000 genes (HVG-reduced, log-normalized)
- obs fields (14): `cell_id`, `dataset_id`, `sample_id`, `stage`, `stage_day_label`, `day_within_stage`, `abs_day`, `time_label`, `scTimeBench_timepoint`, `scTimeBench_cell_type`, `n_genes_by_counts`, `total_counts`, `pct_counts_mt`, `pct_counts_ribo`
- var fields (10): `gene_symbol`, `gene_id`, `gene_name`, `highly_variable`, `is_hvg`, `mt`, `ribo`, `n_cells`, `mean_counts`, `var_counts`. HGNC symbols in both var.index and `var["gene_symbol"]`.
- X: float32 sparse, log1p-normalized (non-integer values, range ~0–6.4)
- Vocab overlap: **1,846 / 2,000 (92.3%)** using `gene_col="gene_symbol"`
- Role: smoke-test input (smaller, faster, all benchmark obs fields present)

**`data/processed/20260329_1343_GSE230659_raw.h5ad`**
- Shape: 75,194 cells × 27,267 genes (full post-QC transcriptome, no HVG subsetting)
- obs fields (11): `sample_id`, `stage`, `stage_day_label`, `day_within_stage`, `abs_day`, `n_genes_by_counts`, `total_counts`, `total_counts_mt`, `pct_counts_mt`, `total_counts_ribo`, `pct_counts_ribo`. Cell barcodes in obs.index encode sample/stage/day.
- var fields (9): HGNC symbols in var.index; `gene_ids` column contains Ensembl IDs. No separate named gene-symbol column.
- X: float32 sparse, **raw integer counts** (values are whole numbers; median ~5.5 per non-zero entry; median non-zero genes per cell: 4,303)
- Vocab overlap: **23,112 / 27,267 (84.8%)** using `gene_col="index"`
- Role: production embedding input (raw counts ideal for scGPT's internal 51-bin discretization; 12× more vocab-matching genes than benchmark h5ad)

**Decision:** No new `adata_scgpt_input.h5ad` needs to be built. The raw h5ad is directly usable for the production embedding run. The benchmark h5ad is used for the smoke test only.

---

#### 18.5 Smoke-test script

Script written at: `benchmark/scgpt/smoke_test.py`

Pipeline implemented:
1. Auto-detect project root (same `_find_project_root()` logic as all other scripts)
2. Pre-flight checks: verify h5ad, model files, CUDA, output directory
3. Load `adata_benchmark.h5ad`; build balanced ~1,000-cell subset: 67 cells per timepoint × 15 timepoints, fixed seed (numpy seed 42)
4. Report vocabulary overlap between `var["gene_symbol"]` and `models/scgpt_whole_human/vocab.json`
5. Run `scg.tasks.embed_data(...)` with `gene_col="gene_symbol"`, `device="cuda"`, `use_fast_transformer=False`, `batch_size=32`
6. Validate embedding: shape check, NaN/Inf check, L2 norm check (expected ≈ 1.0)
7. Run `sc.pp.neighbors(use_rep="X_scGPT")` → `sc.tl.umap()` → `sc.tl.leiden(resolution=0.5, key_added="leiden_scgpt")`
8. Save: `smoke_test_subset.h5ad`, `smoke_test_metadata.json`, `umap_by_timepoint.png`, `umap_by_leiden.png`

Output directory: `benchmark/results/scgpt/smoke_test/`

Run command:
```
conda activate scgpt_env
cd C:\Users\37620\trajectory
python benchmark\scgpt\smoke_test.py
```

---

#### 18.6 Windows compatibility patch — `scgpt/tasks/cell_emb.py`

During the first local run, the smoke test crashed inside the scGPT library at the DataLoader construction step.

**Root cause:** `os.sched_getaffinity(0)` is a POSIX/Linux system call that returns the set of CPU cores the process is pinned to. It does not exist on Windows — `os` has no attribute by that name on Windows. The call was used unconditionally to compute `num_workers` for the PyTorch DataLoader.

**Location:** `C:\Users\37620\Documents\GitHub\scGPT\scgpt\tasks\cell_emb.py`, line 112 (original).

**Original code:**
```python
num_workers=min(len(os.sched_getaffinity(0)), batch_size),
```

**Patch applied:**
```python
if hasattr(os, "sched_getaffinity"):
    _num_workers = min(len(os.sched_getaffinity(0)), batch_size)
else:
    _num_workers = 0  # Windows: run data loading in the main process
```

`num_workers=0` runs DataLoader iteration in the main process with no subprocesses. On Windows this avoids the additional multiprocessing complexity (Windows uses spawn rather than fork, requiring `if __name__ == "__main__":` guards that the scGPT library does not provide). For the smoke test workload, the performance difference is negligible.

**Other occurrences of the same pattern in the repo (not in the active smoke-test path):**
- `scgpt/trainer.py` line 553 — reached only by fine-tuning workflows, not by `embed_data`
- `tutorials/Tutorial_Annotation.ipynb` — notebook-local code, not part of the installed package
- `tutorials/build_atlas_index_faiss.py` — already has a `try/except` guard; already Windows-safe

No other Windows-specific risks were identified in the `get_batch_cell_embeddings` → `embed_data` call chain.

---

#### 18.7 Pipeline state as of session end (scGPT track)

| Component | Status |
|---|---|
| scGPT pretrained model directory | ✓ Complete and validated |
| `scgpt_env` environment | ✓ Confirmed working (torch 2.3.0, CUDA available) |
| Input data decision | ✓ Benchmark h5ad for smoke test; raw h5ad for production |
| `benchmark/scgpt/smoke_test.py` | ✓ Written |
| Windows compatibility patch (`cell_emb.py`) | ✓ Applied |
| Smoke test execution | **Pending** — Windows patch applied; rerun required to confirm embedding completes |
| Production embedding (full 75,194 cells, raw h5ad) | Not started — awaiting smoke test confirmation |
| Pseudo-state clustering and label freezing | Not started |
| `scGPT_pseudostate` wired into benchmark | **Intentionally deferred** — supplementary asset only at this stage |
