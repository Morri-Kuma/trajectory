# Server dry-run notes (real data, 2k-cell subsample)

**Date:** 2026-06-09 · **Profile:** `config.yaml` → `server_dryrun` · **Result:** the
loader → preprocess → harmonize → annotate (surrogate) → 6 comparisons → figures
pipeline **ran end-to-end on the real benchmark `.h5ad` inputs** (2,000 cells per
query) in ~15 s on CPU. Outputs in `results/test_outputs/server_dryrun/{gse178325,gse230659}/`.

This was a **plumbing** test on real data (KNN surrogate, not scANVI). It is not
a biological result. Its value is that it surfaced two concrete, important issues
with the *server input choice* before any GPU time is spent.

## Observed (real, subsampled — illustrative only)
| metric | gse178325 | gse230659 |
|---|---:|---:|
| cells | 2,000 | 2,000 |
| **shared genes after ref∩query intersect** | **181** | **238** |
| pseudotime Spearman (A vs B) | 0.77 | 0.56 |
| correspondence overall agreement | 0.62 | 0.63 |
| fraction unknown/OOD | 0.47 | 0.35 |
| state-graph edge Jaccard | 0.00 | 0.00 |

## Finding 1 — the per-dataset HVG2000 inputs are gene-disjoint (blocker for scANVI)
Each benchmark HVG2000 input has its **own** top-2000 highly-variable genes, so
across datasets they barely overlap (verified):

| intersection | genes | canonical markers present |
|---|---:|---:|
| ref ∩ gse178325 (HVG2000) | 181 | few |
| ref ∩ gse230659 (HVG2000) | 238 | few |
| ref ∩ gse178325 ∩ gse230659 (HVG2000) | **79** | ~0 |
| **ref ∩ q178 ∩ q230 (FULL-GENE)** | **17,829** | **11 / 11** |

scANVI/scArches require the reference and every query on the **same** gene
space. Intersecting the HVG2000 inputs leaves ~80–240 genes and drops the
canonical markers — unusable. **The full-gene inputs share 17,829 genes and keep
all 11 markers**, so they are the correct matrix source.

## Finding 2 — labels and full-gene matrices live in different files
- **Per-cell labels** are only in the HVG2000 inputs:
  `author_cluster_label` → `benchmark/inputs/gse242424_author_cluster_matched/…HVG2000…h5ad`;
  `final_milestone_label_coarse` → `…/gse{178325,230659}_marker_fm_transition_silver_hvg2000/…h5ad`.
- **Full-gene matrices** (29,267 / 27,642 / 27,219 genes) have the genes+markers but **no labels**.
- They share cell barcodes (`obs_names`), so labels can be **joined by cell id**.

## Finding 3 — inconsistent matrix scales
`gse230659` HVG2000 `X` is raw counts (max ≈ 629) while `gse242424` and
`gse178325` HVG2000 `X` are log-normalized (max ≈ 5); full-gene inputs are
log-normalized. **scVI needs integer counts**, so the server prep must start from
the raw counts and apply one uniform normalization — not trust the heterogeneous
`X` already in these files.

## Corrected server recipe (what the real run should do)
1. Build scANVI inputs from the **full-gene** matrices, **joining** the per-cell
   labels from the HVG2000 inputs by `obs_names` (helper:
   `src/data/build_scanvi_inputs.py`).
2. **Freeze one shared gene space** on the reference (HVG selected *within* the
   ~17,829-gene intersection, or the intersection itself), and reindex queries to
   it (scArches contract; zero-fill any missing).
3. Feed **raw counts** to scVI (`layers['counts']`); apply uniform preprocessing.
4. Then run `mode: server` scANVI/scArches as designed.

The `server` profile in `config.yaml` has been annotated to point at the prep
step rather than the raw HVG2000 inputs. The `server_dryrun` profile is retained
to re-validate plumbing on real data anytime.
