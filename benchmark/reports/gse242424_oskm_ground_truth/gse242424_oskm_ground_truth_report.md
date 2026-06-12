# GSE242424 OSKM Reprogramming Ground Truth Benchmark Report

**Generated:** 2026-05-16  
**Provider:** `gse242424_oskm_reprogramming_ground_truth_v1`  
**Label mode:** `ground_truth`  
**State key:** `final_milestone_label_coarse`  
**Result class:** `gse242424_oskm_ground_truth_formal`

## Scope

Results in this report use the author-cluster-matched GSE242424 subset only: 59,187 cells with HVG2000 features. This is not the full 156,969-cell local dataset. All methods were trained and evaluated on:

`benchmark/inputs/gse242424_author_cluster_matched/GSE242424_author_cluster_matched_HVG2000_benchmark_input.h5ad`

The ground truth labels come from the frozen provider built from author ATAC-to-RNA cluster transfer annotations and the matched local scRNA-seq cells.

**Dataset:** GSE242424 OSKM fibroblast reprogramming  
**Timepoints:** `abs_day` in `{0, 2, 4, 6, 8, 10, 12, 14, 16}`  
**State system:** 10 coarse milestones, 9 directed reference edges  
**Methods evaluated:** MIOFlow, PRESCIENT, scNODE  
**Scenarios:**

- **A:** all 9 timepoints used for training/evaluation.
- **B:** train on `{0, 2, 4, 6, 8, 10}` and hold out `{12, 14, 16}` for future extrapolation.
- **C:** train on `{0, 2, 6, 10, 14}`, hold out `{4, 8, 12}` for interpolation and `{16}` for extrapolation.

## Per-Run Metrics

| method | scenario | auroc | auprc | jaccard_topk | embedding_ari | wasserstein |
| --- | --- | --- | --- | --- | --- | --- |
| mioflow | A | 0.8535 | 0.2693 | 0.3846 | 0.0660 | 0.0202 |
| mioflow | B | 0.7387 | 0.2011 | 0.2857 | 0.1569 | 0.0407 |
| mioflow | C | 0.7473 | 0.2677 | 0.3846 | 0.2447 | 0.0252 |
| prescient | A | 0.6947 | 0.1686 | 0.2000 | 0.1244 | 0.0196 |
| prescient | B | 0.7228 | 0.2087 | 0.2000 | 0.2744 | 0.0380 |
| prescient | C | 0.6984 | 0.1952 | 0.2857 | 0.1979 | 0.0254 |
| scnode | A | 0.8400 | 0.4099 | 0.3846 | 0.0344 | 0.0202 |
| scnode | B | 0.6740 | 0.1417 | 0.0588 | 0.0892 | 0.0543 |
| scnode | C | 0.7643 | 0.2181 | 0.0588 | 0.3170 | 0.0327 |

Metric notes:

- `auroc` and `auprc`: lineage-fidelity edge prediction against the frozen reference graph.
- `jaccard_topk`: top-k transition overlap with the reference graph.
- `embedding_ari`: adjusted Rand index of projected embedding clusters versus silver milestone labels.
- `wasserstein`: mean Wasserstein distance between predicted and observed cell distributions; lower is better.

## Per-Run Rankings

Ranks are computed within each scenario across the three methods. Lower rank score is better. AUROC, AUPRC, Jaccard, and ARI are ranked high-to-low; Wasserstein and entropy are ranked low-to-high.

| method | scenario | rank_auroc | rank_auprc | rank_jaccard | lineage_rank_score | rank_ari | emb_rank_score | combined_rank_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mioflow | A | 1.00 | 2.00 | 1.00 | 1.40 | 2.00 | 2.50 | 1.80 |
| mioflow | B | 1.00 | 2.00 | 1.00 | 1.40 | 2.00 | 2.00 | 1.60 |
| mioflow | C | 2.00 | 1.00 | 1.00 | 1.20 | 2.00 | 2.00 | 1.60 |
| prescient | A | 3.00 | 3.00 | 3.00 | 2.80 | 1.00 | 1.00 | 2.20 |
| prescient | B | 2.00 | 1.00 | 2.00 | 2.00 | 1.00 | 1.00 | 1.40 |
| prescient | C | 3.00 | 3.00 | 2.00 | 2.40 | 3.00 | 3.00 | 2.80 |
| scnode | A | 2.00 | 1.00 | 2.00 | 1.80 | 3.00 | 2.50 | 2.00 |
| scnode | B | 3.00 | 3.00 | 3.00 | 2.60 | 3.00 | 3.00 | 3.00 |
| scnode | C | 1.00 | 2.00 | 3.00 | 2.40 | 1.00 | 1.00 | 1.60 |

## Method Summary

| rank | method | mean_auroc | mean_auprc | mean_jaccard | mean_ari | mean_entropy | mean_wasserstein | combined_rank_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | mioflow | 0.7798 | 0.2460 | 0.3516 | 0.1559 | 0.2931 | 0.0287 | 1.6667 |
| 2 | prescient | 0.7053 | 0.1908 | 0.2286 | 0.1989 | 0.2710 | 0.0276 | 2.1333 |
| 3 | scnode | 0.7595 | 0.2566 | 0.1674 | 0.1469 | 0.2935 | 0.0357 | 2.2000 |

## Interpretation

MIOFlow ranks first overall by the combined score, driven by the highest mean AUROC and Jaccard overlap. It is strongest in Scenario A and remains competitive in B and C.

scNODE has the highest mean AUPRC and strong AUROC in Scenarios A and C, but Scenario B is much weaker, especially for Jaccard and embedding ARI. This suggests the late-timepoint extrapolation split is the most difficult setting for scNODE here.

PRESCIENT has the best mean embedding ARI and the lowest mean Wasserstein distance, but lower lineage Jaccard and AUPRC. It appears comparatively good at preserving projected state structure while recovering fewer reference graph edges.

## Caveats

1. Each method was evaluated from one formal run per scenario; run-to-run variance is not estimated here.
2. The ground truth provider is based on author cluster transfer labels, not on an independently validated lineage tree.
3. The frozen graph is a coarse milestone graph for OSKM reprogramming, including productive, partial, stalled, and off-target branches.
4. These results should be compared within this 59,187-cell matched subset, not against full-dataset runs or earlier smoke runs.

## Files

| File | Description |
| --- | --- |
| `gse242424_oskm_ground_truth_core_metrics.csv` | Per-run metric table for the 9 formal runs |
| `gse242424_oskm_ground_truth_rankings.csv` | Within-scenario ranks and composite scores |
| `gse242424_oskm_ground_truth_method_summary.csv` | Per-method means and overall rank |
| `gse242424_oskm_ground_truth_report.md` | This report |

Reference graph: `benchmark/ground_truth/providers/gse242424_oskm_reprogramming_ground_truth_v1/reference_graph.json`  
Provider registry entry: `benchmark/ground_truth/registry.yaml` -> `gse242424_oskm_reprogramming_ground_truth_v1`
