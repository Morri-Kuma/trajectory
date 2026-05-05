# GSE178325 0618 Full-Gene scGPT-v1 Benchmark Report

Date: 2026-05-03

## Scope

Input: `benchmark/inputs/gse178325_scgpt_hvg2000/GSE178325_scGPT_annotated_HVG2000_benchmark_input.h5ad`
Provider: `scgpt_v1_gse178325_0618` (`scgpt_pseudostate_provisional`, 15 states, 38 medium-and-above reference edges used by evaluator).

This report replaces the earlier stage-proxy/HVG-scGPT pass. GSE178325 was processed through the aligned full-gene scGPT workflow: duplicate gene symbols are merged by raw counts before normalization/log1p, scGPT embedding is computed from the full-gene source artifact, and HVG2000 is exported from the scGPT-annotated full-gene object.

## Scenario A

### Forecast Accuracy

| Method | WD | Gaussian MMD | Energy MMD | Hausdorff |
|---|---:|---:|---:|---:|
| scNODE | 0.0202 | 6.21e-06 | 0.0005 | 17.1138 |
| PRESCIENT | 0.0263 | 5.85e-06 | 0.0011 | 17.2194 |
| MIOFlow | 0.0298 | 2.02e-05 | 0.0014 | 17.7762 |

### Embedding Coherence

| Method | ARI | Avg Normalized Entropy |
|---|---:|---:|
| scNODE | 0.1898 | 0.3382 |
| PRESCIENT | 0.1163 | 0.2857 |
| MIOFlow | 0.1091 | 0.2983 |

### Lineage Fidelity

| Method | AUROC | AUPRC | Jaccard | Top-k Jaccard | Single-step | Multi-step |
|---|---:|---:|---:|---:|---:|---:|
| MIOFlow | 0.7898 | 0.3438 | 0.1689 | 0.4074 | 0.4211 | 0.6027 |
| scNODE | 0.7391 | 0.3034 | 0.1689 | 0.2881 | 0.2895 | 0.5342 |
| PRESCIENT | 0.6728 | 0.2616 | 0.1818 | 0.2667 | 0.3684 | 0.4384 |

### Aggregate Ranking

| Rank | Method | Score | Forecast | Embedding | Lineage |
|---:|---|---:|---:|---:|---:|
| 1 | scNODE | 1.806 | 1 | 2 | 2 |
| 2 | PRESCIENT | 1.917 | 2 | 1 | 3 |
| 3 | MIOFlow | 2.222 | 3 | 3 | 1 |

## Scenario B

### Forecast Accuracy

| Method | WD | Gaussian MMD | Energy MMD | Hausdorff |
|---|---:|---:|---:|---:|
| scNODE | 0.0565 | 1.14e-05 | 0.0024 | 21.1867 |
| MIOFlow | 0.0596 | 6.93e-07 | 0.0024 | 22.1044 |
| PRESCIENT | 0.0598 | 2.72e-06 | 0.0026 | 21.8637 |

### Embedding Coherence

| Method | ARI | Avg Normalized Entropy |
|---|---:|---:|
| scNODE | 0.2104 | 0.1007 |
| PRESCIENT | 0.0573 | 0.2205 |
| MIOFlow | 0.0518 | 0.1434 |

### Lineage Fidelity

| Method | AUROC | AUPRC | Jaccard | Top-k Jaccard | Single-step | Multi-step |
|---|---:|---:|---:|---:|---:|---:|
| MIOFlow | 0.6994 | 0.3348 | 0.1642 | 0.3333 | 0.3947 | 0.5890 |
| scNODE | 0.6842 | 0.3626 | 0.1642 | 0.2459 | 0.3684 | 0.5205 |
| PRESCIENT | 0.6748 | 0.2948 | 0.1798 | 0.3103 | 0.3947 | 0.4521 |

### Aggregate Ranking

| Rank | Method | Score | Forecast | Embedding | Lineage |
|---:|---|---:|---:|---:|---:|
| 1 | scNODE | 1.639 | 1 | 1 | 2 |
| 2 | MIOFlow | 1.861 | 1 | 2 | 1 |
| 3 | PRESCIENT | 2.389 | 3 | 2 | 2 |

## Scenario C

### Forecast Accuracy

| Method | WD | Gaussian MMD | Energy MMD | Hausdorff |
|---|---:|---:|---:|---:|
| PRESCIENT | 0.0392 | 3.26e-06 | 0.0014 | 19.9903 |
| scNODE | 0.0486 | 4.40e-06 | 0.0015 | 21.2376 |
| MIOFlow | 0.0495 | 6.02e-07 | 0.0015 | 22.2244 |

### Embedding Coherence

| Method | ARI | Avg Normalized Entropy |
|---|---:|---:|
| scNODE | 0.2712 | 0.2114 |
| MIOFlow | 0.1352 | 0.1779 |
| PRESCIENT | 0.0746 | 0.2438 |

### Lineage Fidelity

| Method | AUROC | AUPRC | Jaccard | Top-k Jaccard | Single-step | Multi-step |
|---|---:|---:|---:|---:|---:|---:|
| MIOFlow | 0.7411 | 0.3987 | 0.1828 | 0.3818 | 0.4737 | 0.5753 |
| scNODE | 0.7204 | 0.4217 | 0.1828 | 0.2459 | 0.3684 | 0.5616 |
| PRESCIENT | 0.6980 | 0.2820 | 0.2000 | 0.3333 | 0.3421 | 0.4795 |

### Aggregate Ranking

| Rank | Method | Score | Forecast | Embedding | Lineage |
|---:|---|---:|---:|---:|---:|
| 1 | MIOFlow | 1.778 | 3 | 1 | 1 |
| 2 | scNODE | 1.917 | 2 | 1 | 2 |
| 3 | PRESCIENT | 2.250 | 1 | 3 | 3 |

## Interpretation

- Scenario A: scNODE ranks first overall (score 1.806).
- Scenario B: scNODE ranks first overall (score 1.639).
- Scenario C: MIOFlow ranks first overall (score 1.778).

Across the aligned full-gene scGPT run, scNODE has the best Wasserstein Distance in Scenarios A and B, while PRESCIENT has the best Wasserstein Distance in Scenario C. scNODE ranks first overall in A and B; MIOFlow ranks first overall in C by combining the best Embedding Coherence and Lineage Fidelity ranks. MIOFlow is consistently strong on Lineage Fidelity, while PRESCIENT remains competitive in Forecast Accuracy but is not top-ranked under the aggregate score in this pass.
