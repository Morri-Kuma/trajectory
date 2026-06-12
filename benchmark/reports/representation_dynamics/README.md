# Representation-Dynamics Supplementary Report

This directory holds the **supplementary** representation-dynamics benchmark
outputs described in `benchmark/docs/scfm_representation_dynamics_plan.md`. These
results are reported **separately** from the official expression-space silver
benchmark and must never be mixed into the official silver report (Risk 1 / 4).

The benchmark compares exactly four representation arms (and no others):

| ID | Representation |
| --- | --- |
| `rep_hvg_pca50` | HVG2000 → PCA50 (expression-derived baseline) |
| `rep_geneformer_cls_pca50` | Geneformer CLS embedding → PCA50 |
| `rep_scgpt_cls_pca50` | scGPT CLS embedding → PCA50 |
| `rep_scfoundation_pca50` | scFoundation cell embedding → PCA50 |

## Output files (conventions)

| File | Produced by | Contents |
| --- | --- | --- |
| `representation_method_summary.md` | `make_representation_reports.py` | Human-readable summary across methods × representations × scenarios. |
| `representation_forecast_rankings.csv` | `make_representation_reports.py` | Aggregated `rep_*` forecast metrics, ranked per scenario. |
| `representation_temporal_signal.csv` | `make_representation_reports.py` | Temporal/state diagnostics per representation (observed cells). |
| `representation_lineage_metrics.csv` | `make_representation_reports.py` | forward/backward/forbidden transition mass per representation. |
| `representation_embedding_coherence.csv` | `make_representation_reports.py` | State separation / silhouette per representation. |
| `representation_metadata_manifest.csv` | `make_representation_reports.py` | One row per artifact with full scFM/representation provenance. |

## How the pieces fit together

1. **Extract frozen scFM embeddings** (Shirokane GPU host, official packages):
   `python -m benchmark.representations.extract_scfm_embeddings --model {geneformer,scgpt,scfoundation} ...`
   → writes `adata.obsm["X_geneformer_cls" | "X_scgpt_cls" | "X_scfoundation"]`
   and `adata.uns["scfm_embedding_metadata"]`.

2. **Build 50-dim representation inputs** (per scenario, **train-only** fitting):
   `python -m benchmark.representations.build_representation_inputs --representation-id <arm> --fit-scope train_only --train-times ... `
   → writes `adata.obsm["X_rep"]` + `adata.uns["representation_inputs_metadata"]`.

3. **Run trajectory models in representation space** using the configs under
   `benchmark/configs/representation/`. The runners use
   `benchmark.representations.model_input.get_model_input_matrix` so the model
   trains on `X_rep`. **Outputs are representation forecasts, not
   gene-expression forecasts.**

4. **Evaluate** with the representation evaluators (not the official
   gene-expression `eval_forecast`):
   - `benchmark/evaluation/eval_representation_forecast.py` →
     `rep_wasserstein_distance`, `rep_gaussian_mmd`, `rep_energy_distance_mmd`,
     `rep_hausdorff_loss`.
   - `benchmark/evaluation/eval_representation_temporal_signal.py` →
     `temporal_signal_ratio`, adjacent/non-adjacent centroid distances,
     `time_prediction_macro_f1`, `time_prediction_r2`,
     `state_centroid_separation`, `state_silhouette_score`,
     forward/backward/forbidden transition mass.

5. **Assemble the report** with `make_representation_reports.py`.

## Anti-leakage

For formal supplementary rankings the scaler/PCA reducers are fit on **training
cells only** (`reducer_fit_scope: train_only`). All-cell fitting is allowed only
as an explicitly-labelled `all_cells_exploratory` mode and is never used for
rankings.
