# Geneformer Representation-Dynamics Work Plan

## 1. Current Status

The goal is to reproduce the completed scGPT representation-dynamics workflow
with Geneformer:

```text
full-gene labeled source
  -> frozen Geneformer cell embedding
  -> scenario-specific train-only StandardScaler/PCA50
  -> X_rep benchmark input
  -> MIOFlow / scNODE / PRESCIENT A/B/C
  -> metrics + supplementary representation report
```

Available implementation pieces:

- `benchmark/representations/extract_scfm_embeddings.py` already supports
  `--model geneformer`.
- `benchmark/representations/scfm_extractors.py` already wraps the official
  `geneformer.TranscriptomeTokenizer` and `geneformer.EmbExtractor`.
- `benchmark/representations/build_representation_inputs.py` already supports
  `rep_geneformer_cls_pca50`.
- `benchmark/configs/representation/` already contains Geneformer configs for
  MIOFlow / scNODE / PRESCIENT across scenarios A/B/C.
- The three trajectory runners already support representation-mode input via
  `X_rep`.

Local status checks:

- The source file exists:
  `benchmark/inputs/representation/gse230659/source/GSE230659_full_gene_with_final_labels.h5ad`
  with shape `75194 x 27219`.
- This source file does not contain `var["ensembl_id"]` or `obs["n_counts"]`,
  but it does contain `var["gene_ids"]` and `obs["total_counts"]`. The current
  extractor can pass these to the official tokenizer with
  `--ensembl-col gene_ids --counts-col total_counts`.
- Local scGPT A/B/C `X_rep` inputs and scGPT trajectory results exist.
- No local `rep_geneformer_cls_pca50` input h5ad or Geneformer result directory
  was found.
- Existing scGPT result directories contain the normal runner metrics
  (`forecast_metrics.json`, `embedding_metrics.json`, `lineage_metrics.json`),
  but not the plan-aligned representation-prefixed outputs
  (`representation_forecast_metrics.json`,
  `representation_temporal_signal.json`). Geneformer should produce these, and
  scGPT should be backfilled with the same evaluator path for strict reporting.

## 2. Deliverables

The Geneformer MVP should produce:

1. Geneformer embedding h5ad:
   `benchmark/inputs/representation/gse230659/source/GSE230659_geneformer_cls_full_gene.h5ad`
2. Geneformer A/B/C X_rep inputs:
   `benchmark/inputs/representation/gse230659/rep_geneformer_cls_pca50/{A,B,C}/GSE230659_rep_geneformer_cls_pca50_{A,B,C}_X_rep.h5ad`
3. Nine trajectory result directories:
   `benchmark/results/representation_dynamics/{mioflow,scnode,prescient}/gse230659_rep_geneformer_cls_pca50_{A,B,C}`
4. Each result directory should contain at minimum:
   `run_metadata.json`, `projected_expression.npy`,
   `forecast_metrics.json`, `embedding_metrics.json`,
   `lineage_metrics.json`.
5. Each result directory should also contain representation-space evaluator
   outputs:
   `representation_forecast_metrics.json`,
   `representation_forecast_per_timepoint.csv`,
   `representation_temporal_signal.json`.
6. The supplementary report should be regenerated under:
   `benchmark/reports/representation_dynamics/`.

## 3. Work Phases

### Phase 0. Validate Geneformer Environment And Input

Purpose: confirm that the current source h5ad can be used with the official
Geneformer tokenizer without changing the Python extractor.

Checks:

- A Shirokane Geneformer environment exists, for example
  `GENEFORMER_CONDA_ENV=geneformer`.
- This import succeeds:
  `python -c "from geneformer import TranscriptomeTokenizer, EmbExtractor"`.
- A Geneformer V2 checkpoint exists, preferably:
  `/home/xzy0723/projects/trajectory/models/Geneformer-V2-104M`.
- The extraction should use:
  `--ensembl-col gene_ids`,
  `--counts-col total_counts`,
  `--model-version V2`,
  `--emb-mode cls`,
  `--emb-layer -1`.

Acceptance criteria:

- The extraction job prints the Geneformer package path, checkpoint path,
  input h5ad shape, and sample `gene_ids` / `total_counts` values before
  running inference.
- If the available checkpoint is V1, stop the CLS plan. V1 requires explicit
  `emb_mode=cell` and a renamed representation id; do not silently mix it with
  `rep_geneformer_cls_pca50`.

### Phase 1. Add Geneformer HPC Scripts

Copy the scGPT shell workflow and make minimal Geneformer-specific changes.
Avoid changing the existing Python core unless the smoke run exposes a real
package-version incompatibility.

Add:

- `jobs/shirokane/run_scfm_geneformer_extract_embedding.sh`
- `jobs/shirokane/run_scfm_geneformer_build_xrep_array.sh`
- `jobs/shirokane/run_scfm_geneformer_trajectory_array.sh`
- `jobs/shirokane/submit_scfm_geneformer_representation_pipeline.sh`

`jobs/shirokane/run_scfm_scgpt_prepare_source.sh` can be reused initially because it prepares
a full-gene source h5ad with final milestone labels; it is not actually
scGPT-specific. It can be renamed/generalized later.

Default Geneformer extraction parameters:

```bash
INPUT_H5AD=benchmark/inputs/representation/gse230659/source/GSE230659_full_gene_with_final_labels.h5ad
OUTPUT_H5AD=benchmark/inputs/representation/gse230659/source/GSE230659_geneformer_cls_full_gene.h5ad
REPRESENTATION_ID=rep_geneformer_cls_pca50
ENSEMBL_COL=gene_ids
COUNTS_COL=total_counts
MODEL_VERSION=V2
EMB_MODE=cls
EMB_LAYER=-1
BATCH_SIZE=64
NPROC=4
```

Acceptance criteria:

- The extraction script writes `adata.obsm["X_geneformer_cls"]`.
- The sidecar metadata records `model_name=Geneformer`, checkpoint,
  paper/code references, matched gene count, embedding dim, `emb_mode=cls`,
  and `model_version=V2`.

### Phase 2. Build Scenario-Specific Train-Only X_rep Inputs

Reuse `benchmark.representations.build_representation_inputs`.

Use the same train times as the scGPT build script:

- A:
  `0.5 2.0 4.0 8.0 12.0 16.0 16.33 16.67 17.0 18.0 20.0 22.0 24.0 28.0 30.0`
- B:
  `0.5 2.0 4.0 8.0 12.0 16.0`
- C:
  `0.5 2.0 8.0 16.0 16.33 16.67 17.0 18.0 22.0 24.0`

Acceptance criteria:

- All three h5ad files contain `adata.obsm["X_rep"]` with shape
  `75194 x 50`.
- `adata.uns["active_representation_id"] == "rep_geneformer_cls_pca50"`.
- `representation_inputs_metadata["rep_geneformer_cls_pca50"]["reducer_fit_scope"] == "train_only"`.
- No `all_cells_exploratory` reducer output is used in formal supplementary
  summaries.

### Phase 3. Smoke Test Then Run The 9 Trajectory Tasks

Use the existing Geneformer configs:

| Task | Method | Scenario | Config |
| --- | --- | --- | --- |
| 1 | MIOFlow | A | `mioflow_gse230659_rep_geneformer_cls_pca50_A.yaml` |
| 2 | MIOFlow | B | `mioflow_gse230659_rep_geneformer_cls_pca50_B.yaml` |
| 3 | MIOFlow | C | `mioflow_gse230659_rep_geneformer_cls_pca50_C.yaml` |
| 4 | scNODE | A | `scnode_gse230659_rep_geneformer_cls_pca50_A.yaml` |
| 5 | scNODE | B | `scnode_gse230659_rep_geneformer_cls_pca50_B.yaml` |
| 6 | scNODE | C | `scnode_gse230659_rep_geneformer_cls_pca50_C.yaml` |
| 7 | PRESCIENT | A | `prescient_gse230659_rep_geneformer_cls_pca50_A.yaml` |
| 8 | PRESCIENT | B | `prescient_gse230659_rep_geneformer_cls_pca50_B.yaml` |
| 9 | PRESCIENT | C | `prescient_gse230659_rep_geneformer_cls_pca50_C.yaml` |

Execution:

1. Smoke run first with `TRAJ_TASK_RANGE=2` (MIOFlow scenario B).
2. After smoke succeeds, run all tasks with `TRAJ_TASK_RANGE=1-9`.
3. Before each task, validate that the config's `dataset.h5ad_path` exists.

Acceptance criteria:

- All nine result directories contain `run_metadata.json`.
- `run_metadata.json` records `representation.representation_id` as
  `rep_geneformer_cls_pca50`.
- Projected output is representation-space output, not gene-expression-space
  output. The final dimension should be 50, unless a method has a documented
  internal latent-output convention.

### Phase 4. Add Representation-Space Metrics

To match `scfm_representation_dynamics_plan.md`, Geneformer should not stop at
the normal runner metrics. It should also emit `rep_*` forecast metrics and
observed-representation diagnostics.

Recommended helper:

```text
benchmark/evaluation/eval_representation_run_from_config.py
```

Inputs:

- representation config yaml
- result directory

Helper responsibilities:

1. Read `dataset.h5ad_path` from the config.
2. Extract observed representation from `adata.obsm["X_rep"]`.
3. Extract observed times from `adata.obs[time_key]`.
4. Extract states from `adata.obs["final_milestone_label_coarse"]`.
5. Read `projected_expression.npy` from the result directory; in a
   representation run, treat it as projected representation.
6. Call `eval_representation_forecast.py` and write
   `representation_forecast_metrics.json`.
7. Call `eval_representation_temporal_signal.py` and write
   `representation_temporal_signal.json`.

Notes:

- Keep the normal `forecast_metrics.json`; it is useful for continuity with the
  completed scGPT workflow.
- Use `representation_forecast_metrics.json` as the primary input to the
  supplementary representation report.
- Run the same helper on the existing scGPT result directories, so scGPT and
  Geneformer use the same reporting path.

Acceptance criteria:

- All nine Geneformer result directories contain
  `representation_forecast_metrics.json`.
- All nine scGPT result directories are backfilled with the same outputs.
- `representation_forecast_metrics.json["metric_space"] == "representation"`.
- Primary forecast keys are `rep_wasserstein_distance`,
  `rep_gaussian_mmd`, `rep_energy_distance_mmd`, and
  `rep_hausdorff_loss`.

### Phase 5. Regenerate Supplementary Reports

Run:

```bash
python -m benchmark.reports.representation_dynamics.make_representation_reports \
  --results-root benchmark/results/representation_dynamics \
  --report-dir benchmark/reports/representation_dynamics
```

Interpretation rules:

- Geneformer vs scGPT can be compared directly under representation-space
  metrics because both use PCA50 active input space with the same cells,
  splits, methods, and seeds.
- Geneformer vs HVG-PCA is a strict controlled baseline only after
  `rep_hvg_pca50` is run through the same representation-input pipeline.
- Do not merge Geneformer/scGPT representation-space forecast rankings into
  the official expression-space benchmark.

Acceptance criteria:

- `representation_forecast_rankings.csv` includes
  `rep_geneformer_cls_pca50` and `rep_scgpt_cls_pca50`.
- `representation_metadata_manifest.csv` includes Geneformer checkpoint,
  upstream release, matched gene count, and embedding dim.
- The summary is clearly marked supplementary and does not overwrite the
  official silver report.

## 4. Tests And Validation

Local code tests:

```bash
python -m pytest benchmark/representations/tests benchmark/evaluation/tests/test_representation_metrics.py
```

HPC artifact validation:

```bash
python - <<'PY'
from pathlib import Path
import anndata as ad

for scenario in "ABC":
    p = Path(f"benchmark/inputs/representation/gse230659/rep_geneformer_cls_pca50/{scenario}/GSE230659_rep_geneformer_cls_pca50_{scenario}_X_rep.h5ad")
    a = ad.read_h5ad(p, backed="r")
    print(scenario, a.shape, a.obsm["X_rep"].shape, a.uns.get("active_representation_id"))
    a.file.close()
PY
```

Result directory validation:

```bash
python - <<'PY'
from pathlib import Path

required = [
    "run_metadata.json",
    "forecast_metrics.json",
    "embedding_metrics.json",
    "lineage_metrics.json",
    "representation_forecast_metrics.json",
    "representation_temporal_signal.json",
]
root = Path("benchmark/results/representation_dynamics")
missing = []
for method in ["mioflow", "scnode", "prescient"]:
    for scenario in "ABC":
        d = root / method / f"gse230659_rep_geneformer_cls_pca50_{scenario}"
        for name in required:
            if not (d / name).exists():
                missing.append(str(d / name))
print("missing:", missing)
raise SystemExit(1 if missing else 0)
PY
```

## 5. Risks And Controls

- Geneformer input fields:
  The source h5ad lacks official fixed columns `ensembl_id` and `n_counts`.
  Use the extractor's supported mapping
  `--ensembl-col gene_ids --counts-col total_counts`; do not mutate the source
  h5ad unless a package-version issue requires it.
- Checkpoint version:
  `emb_mode=cls` requires Geneformer V2. If only V1 is available, switch to
  `emb_mode=cell` and rename the representation id before running benchmarks.
- Metric path drift:
  Existing scGPT results mainly use normal runner metrics. If Geneformer only
  repeats that, it will not fully match the representation-dynamics plan.
  Phase 4 backfills representation-prefixed metrics for both Geneformer and
  scGPT.
- Baseline interpretation:
  Old HVG2000 expression-space formal results are not a strict
  representation-input baseline. Strict comparison requires running
  `rep_hvg_pca50` through the same pipeline.

## 6. Immediate Execution Order

1. Add the four Geneformer shell scripts and submit pipeline.
2. Run Geneformer extraction smoke on Shirokane and verify
   `X_geneformer_cls` plus metadata.
3. Build Geneformer A/B/C `X_rep` inputs.
4. Run MIOFlow B smoke, then the full 9 trajectory tasks.
5. Add and run the representation evaluator backfill helper.
6. Backfill existing scGPT result directories with the same
   `representation_*` metrics.
7. Regenerate `benchmark/reports/representation_dynamics/`.
