# GSE242424 Frozen Reference Lineage Graph Notes

Date created: 2026-05-15  
Last updated: 2026-05-20  
Status: frozen official_silver provider provenance note

This note documents the frozen GSE242424 OSKM reprogramming reference graph and
its author-cluster-matched silver provider:

```text
benchmark/ground_truth/providers/gse242424_oskm_reprogramming_silver_v1/
```

The previous temporary design note has been superseded by this provider. The
graph is now registered in `benchmark/ground_truth/registry.yaml` and used by
formal configs whose `result_class` is `gse242424_oskm_silver_formal`.

## Why GSE242424 Needs Its Own Graph

GSE242424 should not reuse the chemical-reprogramming milestone graph used for
GSE230659/GSE178325. The public metadata describe this dataset as an OSKM-induced
human fibroblast reprogramming time course toward iPSC, with D0-D14 scRNA-seq
samples plus an iPSC endpoint. The biological states center on fibroblast
silencing, OSK-associated states, intermediate/pre-iPSC states, iPSC, and
off-target or stalled branches.

The earlier local GSE242424 input only had time/stage proxy labels:

```text
benchmark/inputs/gse242424_hvg2000/GSE242424_HVG2000_benchmark_input.h5ad
```

That file's `scTimeBench_cell_type` is a time proxy and is not the official
lineage state system. Formal silver runs use the author-cluster-matched subset:

```text
benchmark/inputs/gse242424_author_cluster_matched/GSE242424_author_cluster_matched_HVG2000_benchmark_input.h5ad
```

## Evidence Sources

The provider uses the authors' released analysis products as the primary source
for frozen labels and graph design:

- GEO SuperSeries GSE242424: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242424
- GEO scRNA SubSeries GSE242423: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242423
- Zenodo analysis products, DOI 10.5281/zenodo.8313962: https://zenodo.org/record/8313962
- Zenodo `clusters.tsv`: https://zenodo.org/record/8313962/files/clusters.tsv?download=1
- Author analysis code: https://github.com/kundajelab/scATAC-reprog
- Author data/browser entry point: https://kundajelab.github.io/reprogramming-browser/home.html
- Author ATAC-to-RNA cluster transfer table: https://raw.githubusercontent.com/kundajelab/scATAC-reprog/master/src/analysis/20200828_RNA_Seurat/sessions/20210725_n59378/atac.20210717_n62599.cluster.transfer.tsv
- Author cluster conversion table: https://raw.githubusercontent.com/kundajelab/scATAC-reprog/master/src/figures_factory/configs/cluster.tsv

Zenodo `clusters.tsv` provides the paper-level biological names for 15 clusters:

```text
1   Fibroblast
2   Fibroblast-like
3   Fibroblast-like
4   Fibroblast-like
5   Fibroblast-like
6   Keratinocyte-like
7   hOSK
8   xOSK
9   Intermediate
10  Partially-reprogrammed
11  Intermediate
12  Intermediate
13  Pre-iPSC
14  Pre-iPSC
15  iPSC
```

## Executed Author-Label Matching Workflow

The author-cluster-matched input was built on 2026-05-16 with:

```bash
python scripts/build_gse242424_author_cluster_matched_input.py --overwrite
```

The workflow does not download the full 14GB `scRNA.zip`. It uses the smaller
public author transfer table from GitHub, the old-cluster to paper-cluster
conversion table, and the Zenodo `clusters.tsv` cluster-name dictionary.

The local h5ad stores cells as IDs such as:

```text
GSM7763420_D2_AAACCCAAGGGAGTTC-1
```

The author's merged Seurat object uses sample-specific barcode suffixes, while
the local GEO-derived h5ad was built sample-by-sample and usually uses `-1`.
Therefore full barcode strings are not used directly. The matching rule is:

```text
barcode_core = barcode after removing the trailing -N suffix
local_join_key  = local time_label + "_" + barcode_core
author_join_key = author sample + "_" + barcode_core
```

The executed workflow then performs:

```text
author old atac_cluster
-> paper new_cluster using src/figures_factory/configs/cluster.tsv
-> author_cluster_label using Zenodo clusters.tsv
-> final_milestone_label_coarse using the frozen mapping below
```

Cells without an exact author-label match are excluded from the generated
matched subset. The original full local h5ad remains unchanged.

Outputs:

```text
benchmark/inputs/gse242424_author_cluster_matched/
  GSE242424_author_cluster_matched_HVG2000_benchmark_input.h5ad
  gse242424_author_cluster_matched_cells.tsv
  gse242424_author_cluster_match_summary.json
  gse242424_author_cluster_match_by_time.csv
  gse242424_author_cluster_counts.csv
  source_tables/
    author_atac_to_rna_cluster_transfer.tsv
    author_cluster_conversion.tsv
    author_clusters.tsv
```

Execution result:

```text
local GSE242424 h5ad cells: 156,969
author transfer table rows: 59,378
matched local cells: 59,187
excluded unmatched local cells: 97,782
transfer rows not found in local h5ad: 191
local match rate: 37.71%
```

Matched cells by time point:

```text
D0      9,759 / 16,842
D2      2,523 / 3,804
D4      5,284 / 7,358
D6      7,833 / 11,005
D8      5,497 / 6,921
D10    10,443 / 13,716
D12     7,189 / 9,994
D14     7,534 / 10,315
iPSC    3,125 / 77,014
```

## Frozen Coarse Milestone Mapping

```text
C1        -> fibroblast
C2-C5     -> fibroblast_like_stalled
C6        -> keratinocyte_like
C7        -> hOSK
C8        -> xOSK
C9        -> partial_intermediate
C10       -> partially_reprogrammed
C11-C12   -> primary_intermediate
C13-C14   -> pre_iPSC
C15       -> iPSC
```

## Frozen Reference Graph

The official graph preserves both the productive reprogramming path and the
major off-target/stalled branches:

```text
fibroblast -> fibroblast_like_stalled
fibroblast -> keratinocyte_like
fibroblast -> hOSK -> partial_intermediate -> partially_reprogrammed
fibroblast -> xOSK -> primary_intermediate -> pre_iPSC -> iPSC
```

Provider files:

```text
benchmark/ground_truth/providers/gse242424_oskm_reprogramming_silver_v1/
  state_labels.tsv
  state_metadata.tsv
  annotation_votes.tsv
  reference_graph.json
  reference_graph_edges.csv
  ground_truth_metadata.json
```

## Formal Results

Formal GSE242424 OSKM silver runs use `label_mode=official_silver`,
`state_key=final_milestone_label_coarse`, and
`provider_id=gse242424_oskm_reprogramming_silver_v1`.

Configs:

```text
benchmark/configs/mioflow_gse242424_oskm_silver_A_hvg2000_formal.yaml
benchmark/configs/mioflow_gse242424_oskm_silver_B_hvg2000_formal.yaml
benchmark/configs/mioflow_gse242424_oskm_silver_C_hvg2000_formal.yaml
benchmark/configs/prescient_gse242424_oskm_silver_A_hvg2000_formal.yaml
benchmark/configs/prescient_gse242424_oskm_silver_B_hvg2000_formal.yaml
benchmark/configs/prescient_gse242424_oskm_silver_C_hvg2000_formal.yaml
benchmark/configs/scnode_gse242424_oskm_silver_A_hvg2000_formal.yaml
benchmark/configs/scnode_gse242424_oskm_silver_B_hvg2000_formal.yaml
benchmark/configs/scnode_gse242424_oskm_silver_C_hvg2000_formal.yaml
```

Report:

```text
benchmark/reports/gse242424_oskm_silver/gse242424_oskm_silver_report.md
```

## Caveats

- The provider covers the 59,187 author-cluster-matched cells only. Do not mix
  these results with full-dataset GSE242424 runs.
- The silver provider is based on author cluster transfer labels, not on an
  independently validated lineage tree.
- The graph includes productive, partial, stalled, and off-target branches.
- `iPSC abs_day=16` is a benchmark encoding for the endpoint; it should not be
  overinterpreted as a strict continuous D16 measurement.
