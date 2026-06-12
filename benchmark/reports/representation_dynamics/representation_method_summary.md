# Representation-Dynamics Summary (supplementary)

Runs discovered: **9** under `benchmark\results\representation_dynamics`.

> Supplementary representation-space results. NOT the official gene-expression silver benchmark. scFM embeddings are frozen and were not used to define the official silver labels.

## Forecast rankings (lower is better)

See `representation_forecast_rankings.csv`. Primary ranking metric: `rep_wasserstein_distance` (per scenario).

| method | rep | scenario | rep_wasserstein | temporal_signal_ratio | state_silhouette |
| --- | --- | --- | --- | --- | --- |
| mioflow | rep_geneformer_cls_pca50 | A | 8.173229761597273 | 0.22088579172698217 | 0.032217279443530616 |
| mioflow | rep_geneformer_cls_pca50 | B | 24.64281586934678 | 0.2053231583964255 | 0.02924734144794944 |
| mioflow | rep_geneformer_cls_pca50 | C | 6.9069464215851895 | 0.21486283133193523 | 0.03050105622293457 |
| prescient | rep_geneformer_cls_pca50 | A | 10.918261028696355 | 0.22088579172698217 | 0.032217279443530616 |
| prescient | rep_geneformer_cls_pca50 | B | 14.606198423143004 | 0.2053231583964255 | 0.02924734144794944 |
| prescient | rep_geneformer_cls_pca50 | C | 15.492177322773497 | 0.21486283133193523 | 0.03050105622293457 |
| scnode | rep_geneformer_cls_pca50 | A | 11.154518801685414 | 0.22088579172698217 | 0.032217279443530616 |
| scnode | rep_geneformer_cls_pca50 | B | 13.789106108935167 | 0.2053231583964255 | 0.02924734144794944 |
| scnode | rep_geneformer_cls_pca50 | C | 15.280951556738746 | 0.21486283133193523 | 0.03050105622293457 |

## Interpretation guide

- scFM better forecast **and** comparable temporal signal -> useful state-aware geometry.
- scFM better state separation but worse Scenario B forecast / lower temporal signal -> identity preserved, temporal signal compressed.
- HVG-PCA best across A/B/C -> local expression variation preserves trajectory signal better than generic pretrained reps.
