# Projection-Capable Formal Metrics (HVG2000, scGPT-v1)

| Scenario | Method | Forecast WD (lower) | Gaussian MMD (lower) | Energy MMD (lower) | Hausdorff (lower) | Embedding ARI (higher) | Entropy (lower) | Lineage AUROC (higher) | Lineage AUPRC (higher) | Baseline AUROC (higher) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A | MIOFlow | 163.3884 | 0.00000056 | 0.1281 | 3013.99 | 0.0585 | 0.1387 | 0.7589 | 0.3342 | 0.7812 |
| A | PRESCIENT | 114.6700 | 0.00000052 | 0.0479 | 3257.71 | 0.1449 | 0.1534 | 0.6766 | 0.2323 | 0.7812 |
| A | scNODE | 81.7124 | 0.00000052 | 0.0330 | 2558.69 | 0.0822 | 0.1280 | 0.7611 | 0.3676 | 0.7812 |
| B | MIOFlow | 406.9134 | 0.00000054 | 0.2574 | 3001.51 | 0.1896 | 0.0336 | 0.4990 | 0.1785 | 0.5598 |
| B | PRESCIENT | 153.1575 | 0.00000050 | 0.0808 | 3142.35 | 0.2345 | 0.0624 | 0.5470 | 0.1900 | 0.5598 |
| B | scNODE | 414.6872 | 0.00000050 | 0.1448 | 4634.55 | 0.1162 | 0.1084 | 0.5015 | 0.2264 | 0.5598 |
| C | MIOFlow | 224.4649 | 0.00000060 | 0.0823 | 4801.32 | 0.1097 | 0.1091 | 0.7720 | 0.3131 | 0.7459 |
| C | PRESCIENT | 224.9176 | 0.00000055 | 0.0859 | 4052.48 | 0.2881 | 0.0876 | 0.6524 | 0.2178 | 0.7459 |
| C | scNODE | 266.7125 | 0.00000055 | 0.1062 | 4698.39 | 0.0557 | 0.1556 | 0.7161 | 0.3038 | 0.7459 |

Note: lower is better for Forecast WD, MMD, Hausdorff, and entropy; higher is better for ARI, AUROC, and AUPRC.
Forecast metrics in this table were recomputed for all projection-capable methods with the unified `scTimeBench_geomloss` evaluator. Method-native forecast metrics from older wrappers are retained separately as diagnostics.
This table includes only official projection-capable formal runs under the shared HVG2000 scGPT-v1 input.
