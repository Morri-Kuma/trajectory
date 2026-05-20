# Official Benchmark Reports

This directory is the stable entry point for legacy official reports and the
GSE230659 projected-cell annotation audit. Current silver-standard aggregate
reports are stored in sibling report directories and indexed from the manifest.

The source of truth for reportable result sets is:

```text
benchmark/results/result_manifest.yaml
```

Current reportable result groups:

| Manifest group | Result set | Dataset | Role | Notes |
|---|---|---|---|---|
| `official_silver` | `gse178325_marker_fm_transition_silver_hvg2000` | GSE178325 | current primary benchmark | Frozen marker-FM transition silver provider; MIOFlow, PRESCIENT, and scNODE A/B/C |
| `official_silver` | `gse230659_marker_fm_transition_silver_hvg2000` | GSE230659 | current primary benchmark | Frozen marker-FM transition silver provider; MIOFlow, PRESCIENT, and scNODE A/B/C, plus scenario A WOT/CellRank2 lineage checks |
| `official_silver` | `gse242424_oskm_silver_hvg2000` | GSE242424 | OSKM benchmark | Author-cluster-matched 59,187-cell silver subset; MIOFlow, PRESCIENT, and scNODE A/B/C |
| `official` | `gse230659_scgpt_v1_hvg2000` | GSE230659 | legacy benchmark | Formal A/B/C scGPT-v1 pseudostate results retained for historical comparison |
| `external_validation` | `gse178325_0618_hvg2000` | GSE178325 | legacy validation | A/B/C validation runs for MIOFlow, PRESCIENT, and scNODE |

Pilot, smoke, CPU, reduced, backup, split-check, and Shirokane test outputs are
retained under their historical paths, but they should not be used in official
summaries unless explicitly selected from the manifest.

Primary report files:

```text
benchmark/reports/official_silver/official_silver_model_rankings.md
benchmark/reports/gse242424_oskm_silver/gse242424_oskm_silver_report.md
benchmark/reports/core_summary.csv
benchmark/reports/embedding_summary.csv
benchmark/reports/lineage_summary.csv
```

## GSE230659 Milestone Projected Annotation

As of 2026-05-11, the formal GSE230659 projected-cell Embedding Coherence route
uses the milestone system, not the legacy scGPT pseudostate system.

Projected cells are annotated with two automatic providers in the same HVG2000
gene universe:

```text
projected_expression.npy
-> hvg2000_pca_logistic_regression embedding provider
-> CellTypist classifier provider
-> consensus_milestone_label
```

The full-gene scGPT embedding provider remains useful for observed-cell
provenance and sensitivity analysis, but it is not a formal projected-cell
provider for HVG2000-only projected expression.

The HVG2000 PCA milestone model is a generated artifact expected at:

```text
benchmark/annotation_runs/gse230659_real_sensitivity/embedding_model_pca_hvg2000/hvg_pca_embedding_milestone_model.joblib
```

Regenerate it with `run_gse230659_hvg_pca_embedding_provider.sh` if the local
artifact is absent.

The projected annotation audit files are:

```text
benchmark/reports/official/gse230659_projected_hvg_annotation_audit.md
benchmark/reports/official/gse230659_projected_hvg_annotation_audit_summary.csv
benchmark/reports/official/gse230659_projected_hvg_annotation_label_distribution.csv
benchmark/reports/official/gse230659_projected_hvg_annotation_provider_crosstab.csv
benchmark/reports/official/gse230659_projected_hvg_annotation_embedding_metrics.csv
```

Current scNODE GSE230659 projected-cell audit summary:

| Scenario | Provider match | Ambiguous fraction | Consensus ARI | Embedding-based ARI | Classifier-based ARI |
|---|---:|---:|---:|---:|---:|
| A | 0.814 | 0.186 | 0.033077 | 0.025481 | 0.040033 |
| B | 0.658 | 0.342 | 0.000935 | 0.000015 | 0.008364 |
| C | 0.962 | 0.038 | 0.003935 | 0.002605 | 0.008095 |

Interpretation:

- The technical route is complete for scNODE A/B/C:
  `formal_two_provider_complete=True`, `gene_universe=HVG2000`, exact feature
  order matches, and projected HVG sidecar embeddings are 50-dimensional.
- B has the largest provider disagreement. The current consensus rule is
  conservative: provider disagreement becomes `ambiguous`.
- Any future confidence-based tie-break should be recorded as a versioned policy
  change or sensitivity analysis, not silently mixed with the current formal
  route.

To regenerate a lineage summary from the legacy official manifest group:

```bash
python scripts/summarize_results.py \
    --manifest benchmark/results/result_manifest.yaml \
    --groups official \
    --output-dir benchmark/reports/official
```

To include current official silver lineage sets:

```bash
python scripts/summarize_results.py \
    --manifest benchmark/results/result_manifest.yaml \
    --groups official_silver \
    --output-dir benchmark/reports/official_silver
```
