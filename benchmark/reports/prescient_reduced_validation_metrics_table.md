# PRESCIENT Reduced-Validation Metrics (CPU Low-Memory, HVG2000, scGPT-v1)

| Scenario | Forecast WD (lower) | Gaussian MMD (lower) | Energy MMD (lower) | Hausdorff (lower) | Embedding ARI (higher) | Entropy (lower) | Lineage AUROC (higher) | Lineage AUPRC (higher) | Baseline AUROC (higher) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 3.3461 | 0.1481 | 216.13 | 2288.56 | 0.1194 | 0.1244 | 0.6610 | 0.2220 | 0.7812 |
| B | 3.6448 | 0.1769 | 286.13 | 2980.74 | 0.1276 | 0.0910 | 0.5591 | 0.1830 | 0.5598 |
| C | 4.4809 | 0.1813 | 320.63 | 2689.28 | 0.1594 | 0.1005 | 0.6490 | 0.2201 | 0.7459 |

Note: these rows are `result_class=reduced_validation` and `formal_benchmark=false`; they are excluded from official summaries.
CPU low-memory settings use capped training/reference cells and are intended to validate PRESCIENT integration before a formal run.
