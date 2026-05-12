# Official Benchmark Reports

This directory is the stable entry point for reportable benchmark outputs.

The source of truth for which result directories are official is:

```text
benchmark/results/result_manifest.yaml
```

Current reportable result groups:

| Group | Dataset | Role | Notes |
|---|---|---|---|
| `official` | `gse230659_scgpt_v1_hvg2000` | primary benchmark | Formal A/B/C scGPT-v1 iPSC benchmark results |
| `external_validation` | `gse178325_0618_hvg2000` | external validation | A/B/C validation runs for MIOFlow, PRESCIENT, and scNODE |

Pilot, smoke, CPU, reduced, backup, split-check, and Shirokane test outputs are retained under their historical paths, but they should not be used in official summaries unless explicitly selected from the manifest.

## GSE230659 Milestone Projected Annotation

As of 2026-05-11, the formal GSE230659 projected-cell Embedding Coherence
route uses the milestone system, not the legacy scGPT pseudostate system.

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

The relevant model asset is:

```text
benchmark/annotation_runs/gse230659_real_sensitivity/embedding_model_pca_hvg2000/hvg_pca_embedding_milestone_model.joblib
```

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

- The technical route is complete for scNODE A/B/C: `formal_two_provider_complete=True`, `gene_universe=HVG2000`, exact feature-order matches, and projected HVG sidecar embeddings are 50-dimensional.
- B has the largest provider disagreement. The current consensus rule is conservative: provider disagreement becomes `ambiguous`.
- Any future confidence-based tie-break should be recorded as a versioned policy change or sensitivity analysis, not silently mixed with the current formal route.

To regenerate a lineage summary from the official manifest:

```bash
python scripts/summarize_results.py \
    --manifest benchmark/results/result_manifest.yaml \
    --groups official \
    --output-dir benchmark/reports/official
```

To include the external validation lineage set as well:

```bash
python scripts/summarize_results.py \
    --manifest benchmark/results/result_manifest.yaml \
    --groups official external_validation \
    --output-dir benchmark/reports/official
```
