# Shirokane submission checklist — v2.1 extensions

Two new tracks were added in v2.1: **pseudotime Scenarios D/E/F** on the two chemical
primaries, and **two new iPSC datasets** (GSE298212 human blood reprogramming,
GSE218855 mouse MEF reprogramming) across Scenarios A/B/C. Everything except the
full-data runs is complete locally (registry, configs, builders, jobs, smoke test).
Upload the repo to Shirokane and submit in the order below.

First set the three env vars in your shell (or edit the defaults at the top of each
`.sh`): `TRAJ_PROJECT_ROOT`, `CONDA_SH`, `CONDA_ENV`.

---

## Track 1 — Pseudotime Scenarios D/E/F  (READY NOW, no new data)

Uses the existing GSE230659 / GSE178325 HVG2000 inputs and frozen providers.

```bash
# 1. Build the DPT pseudotime inputs (2 tasks: GSE230659, GSE178325)
qsub jobs/shirokane/run_pseudotime_axis_build.sh

# 2. After (1) finishes, run the D/E/F benchmark (22 tasks).
#    -hold_jid uses the job NAME set by '#$ -N pt_axis_build' so the array
#    really waits for the axis build to finish (do NOT submit them at the same
#    time without the hold — the benchmark needs the pseudotime h5ad to exist).
qsub -hold_jid pt_axis_build jobs/shirokane/run_pseudotime_benchmark_array.sh
```

Produces forecast/embedding/lineage results under
`benchmark/results/{scnode,prescient,mioflow}/gse2306**_marker_fm_silver_{D,E,F}_pseudotime_formal/`
and lineage-only results for wot/cellrank2 Scenario D.

## Track 2 — New datasets GSE298212 + GSE218855  (needs a one-time data prep)

Prerequisite on Shirokane (one remaining manual step — staging the raw data, which
lives on the cluster):

1. Place the raw data under `data/gse298212/` and `data/gse218855/`
   (accessions + provenance in `datasets/dataset_inventory.md`):
   - GSE298212: one 10x directory per sample (PBMC_EPC, S1D1, S1D3, S1D6, S1D8).
   - GSE218855: the processed `.h5ad` (e.g. `GSE218855_FCR_scRNA_raw_qc.h5ad`).

The rest is now in the repo and ready:
- Input builders exist: `scripts/build_gse298212_chemical_input.py` and
  `scripts/build_gse218855_fcr_input.py` (the preprocess job calls them). If your
  GEO download names samples/timepoints differently, adjust `SAMPLE_TO_ABS_DAY` /
  `DAY_TO_ABS_DAY` (and `--time-key` for GSE218855) at the top of each builder.
- Milestone markers exist: GSE298212 is a new section in
  `benchmark/annotation/milestone_markers.yaml`; GSE218855 (mouse) is in
  `benchmark/annotation/milestone_markers_mouse.yaml`. Both are provisional
  marker-defined milestone sets (canonical markers) — refine against author
  annotations if you obtain them.

Then:

```bash
# 3. Build HVG2000 inputs + official_silver providers (2 tasks)
qsub jobs/shirokane/run_v21_new_datasets_preprocess.sh

# 4. After (3) finishes, run A/B/C benchmark for both datasets (22 tasks)
qsub -hold_jid v21_prep jobs/shirokane/run_v21_new_datasets_benchmark_array.sh
```

Post-run integration (to get them into the official ranking tables — same two gates
the pseudotime track needed; the job for (a) is already written):

```bash
# a. official embedding metric for the projection runs (18 tasks, turnkey)
qsub -hold_jid v21_bench jobs/shirokane/run_v21_new_datasets_embedding_coherence_array.sh
# b. add gse298212_* and gse218855_* entries to benchmark/results/result_manifest.yaml
#    (mirror the existing entries; only list scenarios whose result dirs exist), then:
qsub -hold_jid v21_embed jobs/shirokane/run_silver_lineage_graphsim_summary.sh
```
Forecast metrics are already exact-mode from the dispatcher, so only the embedding
metric + the manifest edit are needed.

## Track 3 — Regenerate reports (after Track 1 and/or Track 2 finish)

```bash
qsub -hold_jid run_pseudotime_benchmark,run_v21_new_datasets_benchmark \
     jobs/shirokane/run_silver_lineage_graphsim_summary.sh        # lineage tables
qsub jobs/shirokane/run_silver_forecast_exact_summary.sh          # forecast tables
python -m benchmark.evaluation.summarize_official_silver          # combined rankings
```

---

## What each new file is

| File | Purpose |
|---|---|
| `scripts/generate_v21_benchmark_configs.py` | Emits all 44 new runtime configs (re-runnable, idempotent) |
| `scripts/build_pseudotime_axis.py` | Builds the DPT-bin pseudotime h5ad for D/E/F |
| `scripts/smoke_test_v21_extensions.py` | Local code-path smoke test (run before submitting) |
| `benchmark/configs/runtime/*_{D,E,F}_pseudotime_formal.yaml` | 22 pseudotime configs |
| `benchmark/configs/runtime/*_gse{298212,218855}_*_hvg2000_formal.yaml` | 22 new-dataset configs |
| `jobs/shirokane/run_pseudotime_axis_build.sh` | Track 1 step 1 |
| `jobs/shirokane/run_pseudotime_benchmark_array.sh` | Track 1 step 2 (22 tasks) |
| `jobs/shirokane/run_v21_new_datasets_preprocess.sh` | Track 2 step 3 |
| `jobs/shirokane/run_v21_new_datasets_benchmark_array.sh` | Track 2 step 4 (22 tasks) |

Validate the smoke test first: `PYTHONPATH=. python scripts/smoke_test_v21_extensions.py`
(expects all `[PASS]`).
