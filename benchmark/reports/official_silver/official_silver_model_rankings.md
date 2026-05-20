# Official Silver Model Rankings

Generated from `benchmark/reports/*_summary.csv` using `label_mode == official_silver`.

Ranking convention: ARI, AUROC, AUPRC, Jaccard, and recovery metrics are ranked high-to-low; entropy metrics are ranked low-to-high. Composite rank scores are averages of within-dataset/scenario ranks, so lower is better.

## Combined Summary: scNODE / MIOFlow / PRESCIENT

| overall_combined_rank | method | embedding_runs | lineage_runs | combined_rank_score | mean_ari | mean_weighted_entropy | mean_auprc | mean_auroc | mean_jaccard_topk |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | prescient | 6 | 6 | 1.7333 | 0.1150 | 0.2430 | 0.1819 | 0.4332 | 0.0667 |
| 2 | mioflow | 6 | 6 | 1.7389 | 0.0712 | 0.2928 | 0.1855 | 0.4129 | 0.1000 |
| 3 | scnode | 6 | 6 | 1.7611 | 0.0767 | 0.1740 | 0.1785 | 0.3898 | 0.0667 |

## Embedding Coherence Method Summary

| overall_embedding_rank | method | embedding_runs | mean_embedding_composite_rank_score | mean_ari | median_ari | mean_weighted_entropy |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | scnode | 6 | 1.8889 | 0.0767 | 0.0572 | 0.1740 |
| 2 | prescient | 6 | 2.0000 | 0.1150 | 0.1053 | 0.2430 |
| 3 | mioflow | 6 | 2.1111 | 0.0712 | 0.0652 | 0.2928 |

## Lineage Fidelity Method Summary

| overall_lineage_rank | method | lineage_runs | mean_lineage_composite_rank_score | mean_auprc | mean_auroc | mean_jaccard_topk | mean_single_step_recovery | mean_multi_step_recovery |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | mioflow | 6 | 1.3667 | 0.1855 | 0.4129 | 0.1000 | 0.0556 | 0.1111 |
| 2 | prescient | 6 | 1.4667 | 0.1819 | 0.4332 | 0.0667 | 0.0000 | 0.0556 |
| 3 | scnode | 6 | 1.6333 | 0.1785 | 0.3898 | 0.0667 | 0.0556 | 0.1111 |
| 4 | cellrank2 | 1 | 2.8000 | 0.1935 | 0.4359 | 0.0000 | 0.0000 | 0.0000 |
| 5 | wot | 1 | 3.2000 | 0.1842 | 0.4103 | 0.0000 | 0.0000 | 0.0000 |

## Embedding Winners By Dataset/Scenario

| dataset_id | scenario | method | adjusted_rand_index | weighted_mean_normalized_entropy | embedding_composite_rank_score |
| --- | --- | --- | --- | --- | --- |
| GSE178325 | A | mioflow | 0.1546 | 0.0855 | 1.6667 |
| GSE178325 | B | scnode | 0.0122 | 0.1224 | 1.6667 |
| GSE178325 | C | mioflow | 0.0783 | 0.1742 | 1.6667 |
| GSE178325 | C | prescient | 0.0483 | 0.1182 | 1.6667 |
| GSE230659 | A | mioflow | 0.0189 | 0.1708 | 1.6667 |
| GSE230659 | A | scnode | 0.0814 | 0.2129 | 1.6667 |
| GSE230659 | B | prescient | 0.1536 | 0.1280 | 1.6667 |
| GSE230659 | B | scnode | 0.0093 | 0.0664 | 1.6667 |
| GSE230659 | C | prescient | 0.0354 | 0.2169 | 1.6667 |
| GSE230659 | C | scnode | 0.1401 | 0.3300 | 1.6667 |

## Lineage Winners By Dataset/Scenario

| dataset_id | scenario | method | auprc | auroc | jaccard_topk | single_step_recovery | multi_step_recovery | lineage_composite_rank_score |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE178325 | A | prescient | 0.1641 | 0.4286 | 0.0000 | 0.0000 | 0.3333 | 1.0000 |
| GSE178325 | B | prescient | 0.1613 | 0.4167 | 0.0000 | 0.0000 | 0.0000 | 1.2000 |
| GSE178325 | C | prescient | 0.1510 | 0.3690 | 0.0000 | 0.0000 | 0.0000 | 1.4000 |
| GSE178325 | C | scnode | 0.1396 | 0.3095 | 0.0000 | 0.0000 | 0.3333 | 1.4000 |
| GSE230659 | A | mioflow | 0.2125 | 0.4872 | 0.2000 | 0.0000 | 0.0000 | 1.0000 |
| GSE230659 | B | scnode | 0.2688 | 0.5641 | 0.2000 | 0.3333 | 0.0000 | 1.0000 |
| GSE230659 | C | mioflow | 0.2125 | 0.4872 | 0.2000 | 0.0000 | 0.0000 | 1.0000 |

## Notes

- WOT and CellRank2 only have official_silver lineage results for GSE230659 scenario A in this run, so they are summarized in lineage tables but excluded from the combined embedding+lineage summary.
- Embedding coherence uses `leiden_from_projected_embedding:projected_embedding.npy:n_neighbors=15:resolution=1.0`.
- The complete per-run ranking tables are saved as CSV files in this same directory.
