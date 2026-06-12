# Shirokane scTimeBench Embedding Coherence Run Notes

This note lists the files and commands needed to run the corrected
scTimeBench-aligned Embedding Coherence evaluation on Shirokane.

## Files To Sync

Sync these files from this workspace to the same relative paths under
`/home/xzy0723/projects/trajectory` on Shirokane:

```text
benchmark/adapters/mioflow_adapter.py
benchmark/adapters/prescient_adapter.py
benchmark/adapters/scnode_adapter.py
benchmark/evaluation/backfill_sctimebench_embeddings.py
benchmark/evaluation/eval_embedding_milestone.py
benchmark/evaluation/summarize_official_silver.py
benchmark/methods/MIOFlow/run.py
benchmark/methods/PRESCIENT/run.py
benchmark/methods/scNODE/run.py
jobs/shirokane/run_marker_fm_transition_silver_embedding_backfill_array.sh
jobs/shirokane/run_marker_fm_transition_silver_embedding_coherence_array.sh
jobs/shirokane/run_marker_fm_reference_ari_embedding_array.sh
jobs/shirokane/run_marker_fm_reference_ari_embedding_sensitivity_array.sh
jobs/shirokane/submit_marker_fm_transition_silver_embedding_coherence.sh
jobs/shirokane/submit_marker_fm_reference_ari_embedding.sh
jobs/shirokane/submit_marker_fm_reference_ari_embedding_with_sensitivity.sh
```

## Preflight

The corrected evaluator requires method outputs to contain:

```text
embedding.npy
next_timepoint_embedding.npy
```

`projected_embedding.npy` is still accepted as a fallback for projected
embeddings, but `embedding.npy` is required. If an old result directory does
not have `embedding.npy`, rerun the corresponding scNODE, MIOFlow, or
PRESCIENT method task first with the updated method runner.

For existing formal result directories, the included backfill array can create
these files from cached method artifacts. It does not retrain the methods:

```text
scNODE    trained_scnode_model.pth -> observed VAE latent embedding.npy
MIOFlow   trained_mioflow_model.pth cached PCA -> embedding.npy
PRESCIENT prescient_internal/data.pt scaler/PCA -> embedding.npy
all       projected_embedding.npy -> next_timepoint_embedding.npy
```

## Submit Main Official-Silver Embedding Coherence

```bash
cd /home/xzy0723/projects/trajectory
chmod +x jobs/shirokane/submit_marker_fm_transition_silver_embedding_coherence.sh
chmod +x jobs/shirokane/run_marker_fm_transition_silver_embedding_backfill_array.sh
chmod +x jobs/shirokane/run_marker_fm_transition_silver_embedding_coherence_array.sh
bash jobs/shirokane/submit_marker_fm_transition_silver_embedding_coherence.sh
```

By default this first submits a 18-task backfill array, then holds the
18-task Embedding Coherence evaluator array until backfill completes. Each
evaluation task validates that the output uses:

```text
embedding_metric_protocol=sctimebench
label_source=knn_transfer_from_observed_embedding
entropy_basis=classifier_probability_vector
cluster_source includes n_neighbors=15:resolution=0.5
pred_tp_avg_normalized_entropy is present
```

If you already know every result directory has the new embedding files, skip
the backfill stage with:

```bash
RUN_BACKFILL=0 bash jobs/shirokane/submit_marker_fm_transition_silver_embedding_coherence.sh
```

## Optional Reference-ARI Rerun

Run only the main resolution:

```bash
cd /home/xzy0723/projects/trajectory
chmod +x jobs/shirokane/submit_marker_fm_reference_ari_embedding.sh
chmod +x jobs/shirokane/run_marker_fm_reference_ari_embedding_array.sh
bash jobs/shirokane/submit_marker_fm_reference_ari_embedding.sh
```

Run main plus sensitivity:

```bash
cd /home/xzy0723/projects/trajectory
chmod +x jobs/shirokane/submit_marker_fm_reference_ari_embedding_with_sensitivity.sh
chmod +x jobs/shirokane/run_marker_fm_reference_ari_embedding_array.sh
chmod +x jobs/shirokane/run_marker_fm_reference_ari_embedding_sensitivity_array.sh
bash jobs/shirokane/submit_marker_fm_reference_ari_embedding_with_sensitivity.sh
```

## Regenerate Official-Silver Summary

After the 18 official-silver embedding tasks finish:

```bash
cd /home/xzy0723/projects/trajectory
source /home/xzy0723/miniconda3/etc/profile.d/conda.sh
conda activate traj_env
python -m benchmark.evaluation.summarize_official_silver --leiden-resolution 0.5
```

The generated official report should now rank Embedding Coherence using
`pred_tp_avg_normalized_entropy`, not cluster hard-label entropy.

## Expected Outputs

Each run directory should contain:

```text
benchmark/results/<method>/<run_id>/embedding_milestone_eval/embedding_metrics_official_silver.json
```

The summary step writes:

```text
benchmark/reports/official_silver/official_silver_model_rankings.md
benchmark/reports/official_silver/official_silver_embedding_rankings.csv
benchmark/reports/official_silver/official_silver_embedding_method_summary.csv
benchmark/reports/official_silver/official_silver_combined_method_summary_scnode_mioflow_prescient.csv
```
