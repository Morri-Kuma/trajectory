# scNODE Formal Metrics (HVG2000, scGPT-v1)

| Scenario | Forecast WD (lower) | Gaussian MMD (lower) | Energy MMD (lower) | Hausdorff (lower) | Embedding ARI (higher) | Entropy (lower) | Lineage AUROC (higher) | Lineage AUPRC (higher) | Baseline AUROC (higher) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 156754.48 | 0.0788 | 131.92 | 2558.69 | 0.0822 | 0.1280 | 0.7611 | 0.3676 | 0.7812 |
| B | 815520.59 | 0.2391 | 579.07 | 4634.55 | 0.1162 | 0.1084 | 0.5015 | 0.2264 | 0.5598 |
| C | 526241.69 | 0.2055 | 424.89 | 4698.39 | 0.0557 | 0.1556 | 0.7161 | 0.3038 | 0.7459 |

Note: lower is better for Forecast WD, MMD, Hausdorff, and entropy; higher is better for ARI, AUROC, and AUPRC.
WOT and CellRank2 are not included in this table because they are lineage-only under the current benchmark capability rules.
