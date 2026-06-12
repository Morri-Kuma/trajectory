# scTimeBench-Aligned iPSC Trajectory Benchmark

**Author**: Kuma, Graduate School of Frontier Sciences, The University of Tokyo  
**Primary domain**: human iPSC reprogramming trajectories  
**Framework**: official silver trajectory benchmark

---

## Project Overview

This repository implements a benchmark for trajectory inference and forecasting
methods in human iPSC reprogramming. It is aligned with the evaluation
philosophy of scTimeBench, while adapting the datasets, state systems,
reference graphs, and reporting workflow to iPSC reprogramming.

The project is not a direct copy of scTimeBench. It uses scTimeBench's core
evaluation dimensions as the base and adds domain-specific benchmark assets:

- frozen silver-standard milestone providers for GSE178325 and GSE230659;
- an OSKM-specific GSE242424 ground truth provider built from author-cluster-matched
  annotations;
- method capability gating, so each method is evaluated only on supported
  tasks;
- per-method, per-scenario outputs for forecast accuracy, embedding coherence,
  and lineage fidelity;
- explicit separation of current reportable results from pilot, smoke, and
  archived/debug outputs.

## Current Status

The benchmark uses only the current registered ground-truth providers for
reportable runs.

| Component | Status |
|---|---|
| GSE178325 marker-FM transition silver benchmark | Formal A/B/C result sets available for MIOFlow, PRESCIENT, and scNODE |
| GSE230659 marker-FM transition silver benchmark | Formal A/B/C result sets available for MIOFlow, PRESCIENT, and scNODE; WOT and CellRank2 have lineage-only scenario A checks |
| GSE242424 OSKM ground truth benchmark | Formal A/B/C result sets available for MIOFlow, PRESCIENT, and scNODE on the 59,187-cell author-cluster-matched subset |
| Forecast Accuracy | Active for generative/projected-cell methods |
| Embedding Coherence | Active for methods with projected embeddings or projected cells |
| Lineage Fidelity | Active for methods with state-level transition output |
Current report entry points:

- `benchmark/reports/official_silver/official_silver_model_rankings.md`
- `benchmark/reports/official_silver/official_silver_combined_method_summary_scnode_mioflow_prescient.csv`
- `benchmark/reports/official_silver/official_silver_embedding_method_summary.csv`
- `benchmark/reports/official_silver/official_silver_lineage_method_summary.csv`
- `benchmark/reports/gse242424_oskm_ground_truth/gse242424_oskm_ground_truth_report.md`
- `benchmark/results/result_manifest.yaml`

## Evaluation Dimensions

| Dimension | Status | Eligible methods |
|---|---|---|
| Forecast Accuracy | Active | Methods that generate projected expression at held-out or future time points, such as MIOFlow, PRESCIENT, and scNODE |
| Embedding Coherence | Active | Methods that provide projected embeddings or projected cells that can be embedded |
| Lineage Fidelity | Active | Methods that provide or can be converted to state-level transition predictions |

Capability gating is intentional. WOT and CellRank2 are not forced into
forecast or embedding metrics because they do not natively generate
unseen-timepoint projected cells in the same sense as generative trajectory
models.

## Benchmark Scenarios

| ID | Time axis | Type |
|---|---|---|
| A | Observed time | Interpolation / observed-time evaluation |
| B | Observed time | Extrapolation / held-out future-time evaluation |
| C | Observed time | Mixed interpolation and extrapolation |
| D | Pseudotime | Interpolation |
| E | Pseudotime | Extrapolation |
| F | Pseudotime | Mixed interpolation and extrapolation |

Current formal reports focus on observed-time scenarios A/B/C. Pseudotime
scenarios remain part of the framework design and can be activated when the
corresponding inputs and reference state systems are frozen.

## Repository Structure

```text
trajectory/
|-- README.md                              # project-level overview
|-- benchmark/                            # benchmark package, configs, inputs, methods, reports, and results
|   |-- README.md                         # benchmark-specific usage notes
|   |-- adapters/                         # method adapters and capability gates
|   |-- annotation/                       # milestone/provider annotation utilities
|   |-- configs/                          # method, scenario, runtime, and representation configs
|   |-- evaluation/                       # forecast, embedding, lineage, and representation evaluators
|   |-- ground_truth/                     # frozen state labels and reference graphs
|   |-- inputs/                           # benchmark-ready h5ad inputs and QC summaries
|   |-- methods/                          # vendored or wrapped method implementations
|   |-- reports/                          # aggregated reports, rankings, and figures
|   |-- representations/                  # scFM representation extraction and input builders
|   |-- results/                          # per-method per-scenario outputs
|   `-- shared/                           # shared dataset registry and utility code
|-- data/                                 # raw and processed dataset assets outside benchmark-ready inputs
|-- docs/                                 # project-level documents, plans, reports, and meeting materials
|   |-- framework/                        # main benchmark design document
|   |-- meetings/                         # meeting notes and presentation materials
|   |-- method_notes/                     # method-specific notes and test plans
|   `-- reports/                          # audit and data reports
|-- jobs/                                 # cluster and batch execution entry points
|   `-- shirokane/                        # Shirokane qsub run/submit scripts, docs, and portable bundles
|-- logs/                                 # local and cluster run logs
|-- models/                               # downloaded model checkpoints and model assets
|-- reference/                            # papers and external-reference notes
`-- scripts/                              # standalone data preparation, validation, and plotting scripts
```

## Methods

The current benchmark includes:

| Method | Forecast Accuracy | Embedding Coherence | Lineage Fidelity |
|---|---:|---:|---:|
| MIOFlow | yes | yes | yes |
| PRESCIENT | yes | yes | yes |
| scNODE | yes | yes | yes |
| WOT | no | no | yes |
| CellRank2 | no | no | yes |

Method capability flags are configured in
`benchmark/configs/method_capabilities.yaml`.

## Running the Benchmark

Representative observed-time configs are stored in `benchmark/configs/` and
`benchmark/configs/runtime/`.

Example GSE242424 OSKM ground truth run:

```bash
python benchmark/methods/scNODE/run.py \
    --config benchmark/configs/scnode_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml
```

Example marker-FM transition silver runtime config:

```bash
python benchmark/methods/scNODE/run.py \
    --config benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_A_hvg2000_formal.yaml
```

For Shirokane/qsub batch runs, use the helper scripts in `jobs/shirokane/`,
for example:

```bash
bash jobs/shirokane/submit_silver_forecast_exact.sh
bash jobs/shirokane/submit_silver_lineage_graphsim.sh
```

The dispatcher and method runners read method capabilities and skip unsupported
metrics with explicit metadata rather than producing forced or invalid outputs.

## Ground Truth and Lineage Fidelity

Lineage Fidelity depends on three frozen assets:

1. a benchmark cell-state system;
2. a reference state-level lineage graph;
3. a fixed rule for aggregating method predictions to state-level transitions.

Current primary providers are registered in `benchmark/ground_truth/registry.yaml`:

- `gse178325_marker_fm_transition_silver_v1`
- `gse230659_marker_fm_transition_silver_v1`
- `gse242424_oskm_reprogramming_ground_truth_v1`

The legacy `scgpt_v1` providers are retained for backward compatibility and
historical comparison only. New annotation systems should be added as new
versioned providers instead of changing existing frozen provider directories.

Per-run lineage outputs include `state_transition_matrix.csv`,
`lineage_graph_edges.csv`, `lineage_metrics.json`, and diagnostic metadata where
available.

## Relationship to scTimeBench

This project is best understood as a domain-specific extension of scTimeBench
ideas:

- scTimeBench provides the general benchmark vocabulary and metric families.
- trajectory specializes the benchmark to iPSC reprogramming.
- trajectory uses frozen ground-truth state systems and iPSC-specific
  reference graphs.
- trajectory keeps method capability gating explicit so transition-only methods
  and generative methods are compared on appropriate evaluation surfaces.

The goal is therefore not strict code-structure parity with scTimeBench, but a
scTimeBench-aligned benchmark that answers a narrower biological question.

## Framework Reference

`docs/framework/experimental_framework_v2.md` is the main design document for benchmark logic.
If implementation notes and README text diverge, treat the framework document
and current configs/results as the more specific source of truth.

---

## Annotation-robustness MSc branch (added 2026-06)

A focused, three-month, computation-only slice sits on top of this benchmark:
**does the cell-annotation strategy change the inferred trajectory?** It compares
the existing marker-"silver" labels against **scANVI/scArches reference-mapping**
from the labelled GSE242424 dataset, across GSE178325 and GSE230659.

New top-level documents:

- `research_design_three_month_project.md` — full project design (MVP + manuscript versions).
- `scanvi_annotation_branch.md` — scANVI/scArches design + the six comparisons.
- `datasets/dataset_inventory.{md,csv}` — 3 required + 5 verified additional datasets.
- `reference/literature_review_extra.md` — existing + 8 newly verified references.
- `figure_plan.md` — five main figures + supplements.
- `docs/reports/project_audit_report.md` — critical audit of this repository.
- `deletion_proposal.md` — cleanup candidates (awaiting approval; nothing deleted).
- `TODO.md` — running task list.

New modular code (`src/`) and a config-driven entry point (`config.yaml`,
`mode: test|server`):

```text
src/{data,preprocessing,annotation,trajectory,evaluation,plotting,utils}/
```

### Run the local smoke test (CPU, no big downloads)

```bash
pip install -r requirements-annotation-branch.txt   # core deps only
python scripts/run_smoke_test.py                     # synthetic data -> full pipeline
# outputs: results/test_outputs/annotation_branch/{gse178325,gse230659}/ + summary.json
```

The smoke test generates tiny synthetic reference/query data, runs the whole
annotation→trajectory→comparison pipeline with a KNN-on-PCA surrogate for
scANVI, and asserts all tables/figures are produced. It proves the pipeline
**executes**; real results come from `mode: server` with scvi-tools + GPU on the
full `.h5ad` inputs. See `research_design_three_month_project.md` §16.

### Server run (Shirokane)

The full scANVI run is a 3-stage `qsub` chain (prep → train → map+compare):

```bash
bash jobs/shirokane/submit_scanvi_annotation_branch.sh
```

See `jobs/shirokane/docs/scanvi_annotation_branch_jobs.md` for env overrides (conda envs, GPU queue, input paths) and the raw-counts / barcode-join pre-checks.

### Tests & local reproducibility

Unit tests cover the dependency-light core (config loading, exact count
reconstruction, annotation-comparison metrics, trajectory utilities, surrogate
annotation, preprocessing):

```bash
pip install -r requirements-annotation-branch.txt pytest
pytest          # 24 tests, ~3s (scoped to tests/ via pytest.ini)
```

Run the full annotation-comparison + trajectory pipeline on the **real prepped
inputs** locally (surrogate annotation, CPU, subsampled — no scvi/GPU), once the
prep step has produced `results/annotation_branch/inputs/*.h5ad`:

```bash
python -c "import yaml;c=yaml.safe_load(open('config.yaml'));c['mode']='local_real';open('/tmp/lr.yaml','w').write(yaml.safe_dump(c))"
python -m src.pipeline --config /tmp/lr.yaml   # writes results/test_outputs/local_real_run/
```

`environment.yml` pins the full stack (core + scanpy + scvi-tools).

### First server run: finding & manuscript

The cross-protocol scANVI run completed: reference self-recovery **0.990**, but the
query mapping **collapsed** — 77–96% of chemical-reprogramming cells were confidently
assigned to the OSKM-only `hOSK` state, which scANVI confidence did **not** flag
(≤0.3% OOD). A geometric out-of-distribution score (`src/annotation/ood.py`) instead
flags **39–44%** of query cells (`scripts/demo_ood_recalibration.py`). This reframes
the project as a cautionary + diagnostic methods study.

- Read-out: `docs/scanvi_run_results_interpretation.md`
- Draft article: `manuscript/manuscript_draft.{md,docx}`
- Remaining validation (positive control, same-protocol transfer, scANVI-latent OOD via `scripts/score_ood.py`): `docs/validation_runbook.md`

Tests now: **24** (`pytest`). Repo: `LICENSE` (MIT), `CITATION.cff`, `environment.yml`.
