# Projection-Capable Formal Metrics (HVG2000, scGPT-v1)

| Scenario | Method | Forecast WD (lower) | Gaussian MMD (lower) | Energy MMD (lower) | Hausdorff (lower) | Embedding ARI (higher) | Entropy (lower) | Lineage AUROC (higher) | Lineage AUPRC (higher) | Baseline AUROC (higher) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A | PRESCIENT | 3.1054 | 0.1316 | 193.76 | 2849.61 | 0.1449 | 0.1534 | 0.6766 | 0.2323 | 0.7812 |
| A | scNODE | 156754.4832 | 0.0788 | 131.92 | 2558.69 | 0.0822 | 0.1280 | 0.7611 | 0.3676 | 0.7812 |
| B | PRESCIENT | 4.0936 | 0.1842 | 324.40 | 2393.48 | 0.2345 | 0.0624 | 0.5470 | 0.1900 | 0.5598 |
| B | scNODE | 815520.5920 | 0.2391 | 579.07 | 4634.55 | 0.1162 | 0.1084 | 0.5015 | 0.2264 | 0.5598 |
| C | PRESCIENT | 4.4904 | 0.2081 | 336.54 | 4526.41 | 0.2881 | 0.0876 | 0.6524 | 0.2178 | 0.7459 |
| C | scNODE | 526241.6885 | 0.2055 | 424.89 | 4698.39 | 0.0557 | 0.1556 | 0.7161 | 0.3038 | 0.7459 |

Note: lower is better for Forecast WD, MMD, Hausdorff, and entropy; higher is better for ARI, AUROC, and AUPRC.
This table includes only official projection-capable formal runs under the shared HVG2000 scGPT-v1 input.
