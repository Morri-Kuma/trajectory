# scANVI Annotation Branch — Design & Implementation

**Date:** 2026-06-09
**Reference pipeline:** scArches scANVI surgery — https://docs.scarches.org/en/latest/scanvi_surgery_pipeline.html
**Status:** design + modular code (`src/annotation/`), runnable on small local subsets first.

---

## 1. Purpose and scope

The benchmark currently annotates the two label-free chemical-reprogramming
datasets (GSE178325, GSE230659) with an in-repo **marker "silver"** pipeline.
This branch adds an independent, learned annotation — **scANVI reference-mapping
from the labelled GSE242424 dataset via scArches** — and compares the two so we
can quantify how much trajectory results depend on annotation choice.

This document is the contract for `src/annotation/` and the comparison code in
`src/evaluation/annotation_compare.py`. It is written so the pipeline can run on
a tiny laptop subset (`test_mode`) and scale unchanged on Shirokane.

---

## 2. The honest caveat that shapes the whole branch

| | Reference (GSE242424) | Queries (GSE178325, GSE230659) |
|---|---|---|
| Reprogramming type | **OSKM (transcription-factor)** | **Chemical** |
| Starting cell | Fibroblast | hADSC (adipose-derived stromal) |
| Label source | Author ATAC→RNA transfer (per-cell) | None (in-repo marker silver) |
| Coarse label space | Fibroblast, Fibroblast-like, Keratinocyte-like, Intermediate, Partially-reprogrammed, Pre-iPSC, iPSC, hOSK, xOSK | hADSCs, epithelial_like, intermediate_plastic, xen_like, hCiPS |

The two label spaces are **related but not identical**. Start and terminal
states correspond (fibroblast↔hADSCs as the somatic origin; iPSC↔hCiPS as
pluripotent terminus); intermediate states correspond only partially. Therefore:

- scANVI predictions are a **label-correspondence hypothesis**, not a gold standard.
- Agreement is measured with **correspondence-aware** metrics (contingency, NMI, per-state recall *where a correspondence is defined*), **never** naive accuracy over mismatched label sets.
- A **curated correspondence map** (below) is fixed in advance and validated with canonical markers.

This caveat is repeated in every results table the branch produces.

### 2.1 Pre-registered correspondence map (coarse)

| Reference label (GSE242424) | Query silver label | Basis |
|---|---|---|
| Fibroblast / Fibroblast-like | hADSCs | somatic origin; COL1A2, DCN, fibroblast program |
| Keratinocyte-like / Intermediate | epithelial_like | MET / epithelial program (CDH1, EPCAM) |
| Partially-reprogrammed / Pre-iPSC | intermediate_plastic | partial pluripotency activation |
| iPSC | hCiPS | POU5F1, SOX2, NANOG, LIN28A |
| hOSK / xOSK | (no query equivalent) | OSKM-transgene-specific states; expected "unknown/OOD" in queries |
| (no reference equivalent) | xen_like | XEN program (GATA6, SOX17); chemical-specific, may map to "unknown" |

`hOSK/xOSK` are transgene-driven OSKM states with **no chemical analogue**, so
query cells should *not* map there; if they do, that is a flagged artefact. The
chemical-specific `xen_like` state has **no OSKM analogue** and is expected to
fall into low-confidence "unknown_or_ood" — itself an informative result.

---

## 3. Pipeline overview

```text
GSE242424 (labelled, OSKM)                 GSE178325 / GSE230659 (label-free, chemical)
        │                                                   │
        ▼                                                   │
[A] Build reference AnnData                                 │
    - common HVG gene space (intersect with queries)        │
    - obs['ref_label'] = author coarse labels               │
        │                                                   │
        ▼                                                   │
[B] Train scVI on reference  ──►  [C] Train scANVI (semi-supervised)
        │                                                   │
        ▼                                                   ▼
[D] Freeze reference scANVI  ──►  [E] scArches surgery: map each query
                                       (load_query_data → train query-specific weights)
                                                            │
                                                            ▼
                                       [F] Predict labels + confidence
                                            obs['scanvi_label'], obs['scanvi_confidence']
                                            low-confidence → 'unknown_or_ood'
                                                            │
                                                            ▼
                                       [G] Comparison vs marker-silver (§5)
```

Steps map to scArches' published scANVI surgery tutorial: train reference
`SCVI`→`SCANVI`, then `SCANVI.load_query_data(query, ref_model)` and
`q.train(...)`, then `q.predict()` for query labels. We add gene-space
harmonization, confidence thresholding/OOD flagging, and the comparison layer.

---

## 4. Implementation details

### 4.1 Reference preparation (`src/annotation/scanvi_reference.py`)
- Load `benchmark/inputs/gse242424_author_cluster_matched/…HVG2000…h5ad`.
- Set `ref_label = author_cluster_label` (the 9 author ATAC->RNA clusters; **chosen 2026-06-09**). Alt: `final_milestone_label_coarse` (10 states; only difference is it splits the author `Intermediate` cluster into primary/partial). Configurable via `config.yaml`.
- Restrict to a **frozen shared gene space** (scArches requires identical `var` order; reference defines it, queries reindex/zero-fill). **Verified 2026-06-09:** the per-dataset HVG2000 benchmark inputs are gene-disjoint (ref∩query ≈ 181–238; all-three = 79; canonical markers absent) and must NOT be used directly. Build inputs from the **full-gene** matrices (ref∩q178∩q230 = 17,829 genes, all 11 markers) with labels joined by barcode via `src/data/build_scanvi_inputs.py`, then freeze ~2,000 HVGs within that intersection and feed raw counts to scVI.
- Counts layer (`layers['counts']`) is the scVI input; library-size handled by scVI.
- Save `models/scanvi_reference_gse242424/` (scVI + scANVI) and the frozen gene list.

### 4.2 scVI / scANVI training (`src/annotation/scanvi_train.py`)
- `SCVI.setup_anndata(ref, layer='counts', batch_key='sample_id')`.
- scVI: `n_latent` (default 30), `n_layers` (2), trained `max_epochs` (test: 2; server: 200–400 with early stopping).
- scANVI from the scVI model with `labels_key='ref_label'`, `unlabeled_category='unknown'`; train `max_epochs` (test: 2; server: 20–100).
- Sanity gate: reference label **recovery accuracy** (scANVI re-predicting its own labels) must clear a threshold (e.g. ≥0.9 on held-out reference cells) before mapping queries.

### 4.3 Query mapping (`src/annotation/scanvi_map.py`)
- Reindex query to the frozen reference gene space; `batch_key='sample_id'`.
- `SCANVI.load_query_data(query, ref_model_path)`; `q.train(max_epochs=…, plan_kwargs={'weight_decay':0})` per the surgery recipe.
- `q.predict()` → `scanvi_label`; `q.predict(soft=True)` → posterior; `scanvi_confidence = max posterior`.
- `scanvi_unknown_flag = confidence < tau` (default `tau=0.5`, configurable) → relabel as `unknown_or_ood`.
- Export latent `q.get_latent_representation()` for UMAP (joint ref+query embedding).
- Persist query `.h5ad` with new `obs` columns + `uns['scanvi_branch']` provenance (config hash, ref id, tau, gene-space size).

### 4.4 test_mode surrogate (`src/annotation/knn_surrogate.py`)
When `scvi`/`torch` are unavailable (typical laptop), a **KNN-on-PCA** surrogate
classifier (scikit-learn) trained on the reference provides drop-in labels so the
**comparison and trajectory code can be exercised end-to-end** without a GPU. It
is clearly labelled `annotation_method='knn_surrogate'` and is **never** used for
manuscript results — only to prove the plumbing. This keeps the local smoke test
dependency-light while the real scANVI path runs on the server.

---

## 5. Comparison design (the six required comparisons)

Implemented in `src/evaluation/annotation_compare.py`; each emits a table + figure.

1. **Cell-type composition comparison.** Per-timepoint stacked composition under marker-silver vs scANVI; report total-variation distance per timepoint and overall. *Output:* `composition_compare.csv`, stacked-bar figure.

2. **UMAP / latent-space visualization.** Joint scANVI latent UMAP of reference + query, coloured by (a) dataset, (b) marker-silver label, (c) scANVI label, (d) confidence. Shows where queries land on the reference manifold. *Output:* `latent_umap_*.png`.

3. **Confusion / agreement matrix.** Contingency table marker-silver × scANVI (raw counts + row-normalized), plus **NMI** and **adjusted Rand index**, and **per-state recall against the pre-registered correspondence map** (§2.1). Because label sets differ, the contingency table is the primary object; scalar metrics are secondary. *Output:* `agreement_contingency.csv`, `agreement_metrics.json`, heatmap.

4. **Marker gene validation.** For each corresponding state, score canonical markers (pluripotency POU5F1/SOX2/NANOG/LIN28A; stromal/fibroblast COL1A2/DCN/LUM; epithelial/MET CDH1/EPCAM; XEN GATA6/SOX17) and test which annotation yields cleaner marker enrichment per assigned state (AUROC of marker score separating the state). *Output:* `marker_validation.csv`, dotplot.

5. **Pseudotime / trajectory structure comparison.** Run DPT (and PAGA) under each annotation with identical preprocessing/seed/root; compare per-cell pseudotime (Spearman), branch assignment of intermediates, and PAGA graph similarity (edge Jaccard / weighted graph distance). *Output:* `pseudotime_corr.json`, `paga_graph_compare.csv`, side-by-side trajectory plots.

6. **Effect on biological interpretation.** Identify cells where the two annotations disagree; characterise them (marker co-expression, QC, doublet score) and report the transition(s) whose biological identity is annotation-dependent. *Output:* `disagreement_cells.csv`, short interpretation note.

---

## 6. Effect on downstream benchmark (sensitivity)

Beyond annotation-vs-annotation, re-score **Lineage Fidelity** for the existing
methods (scNODE, PRESCIENT, MIOFlow, WOT, CellRank2) using each annotation as the
state system, and report **ranking deltas**. Forecast Accuracy (label-independent)
is the invariant control. *Output:* `lineage_fidelity_sensitivity.csv` with a
"rank under A / rank under B / Δ" column. (Full-data, server-side.)

---

## 7. Outputs & directory layout

```text
src/annotation/
    scanvi_reference.py     # build reference AnnData (gene-space harmonization)
    scanvi_train.py         # scVI + scANVI training (+ reference recovery gate)
    scanvi_map.py           # scArches surgery + predict + confidence/OOD
    knn_surrogate.py        # test_mode dependency-light classifier
    __init__.py
src/evaluation/
    annotation_compare.py   # the six comparisons (§5)
models/
    scanvi_reference_gse242424/    # frozen reference model + gene list (gitignored)
results/test_outputs/annotation_branch/
    composition_compare.csv, agreement_contingency.csv, agreement_metrics.json,
    marker_validation.csv, pseudotime_corr.json, *.png   # smoke-test artifacts
```

All query `.h5ad` gain: `scanvi_label`, `scanvi_confidence`, `scanvi_unknown_flag`,
and `uns['scanvi_branch']` provenance. Existing frozen providers are **not**
modified — scANVI labels are written as a *new* annotation layer/provider
(`gse178325_scanvi_ref242424_v1`, `gse230659_scanvi_ref242424_v1`), per the
project rule "add new versioned providers instead of changing frozen ones."

---

## 8. Validation strategy (does the branch itself work?)

- **Positive control (GSE175634).** Run the same reference-mapping idea where the query *does* have author labels: train a reference on part of the labelled data, map the rest, and confirm scANVI recovers known labels at the extremes. This validates the machinery independently of the OSKM↔chemical mismatch.
- **Reference recovery gate.** scANVI must re-predict GSE242424's own held-out labels at ≥0.9 before any query mapping is trusted.
- **Confidence calibration.** Report the confidence distribution and the fraction flagged `unknown_or_ood`; an all-confident or all-unknown result is a red flag to investigate, not to publish.

---

## 9. Local vs server execution

| Aspect | `test_mode` (laptop) | `server_mode` (Shirokane) |
|---|---|---|
| Cells | `max_cells_per_dataset` (e.g. 2,000) | full |
| Genes | HVG2000 | HVG2000 or full-gene |
| Annotation | KNN surrogate **or** 2-epoch CPU scANVI | full scANVI/scArches on GPU |
| Trajectory | DPT + PAGA | DPT + PAGA + Slingshot + generative methods |
| Purpose | prove the pipeline executes | produce manuscript results |

The same `src/` code runs in both; only `config.yaml` changes.

---

## 10. Limitations (stated, not hidden)

- Cross-protocol transfer (OSKM→chemical) may be imperfect; the branch reports correspondence, not accuracy.
- HVG2000 + gene-intersection reduces the feature space; full-gene runs are a server robustness check.
- scANVI confidence is model-internal, not a biological truth; OOD flags are heuristic (`tau`).
- No new wet-lab validation is possible; biological reads are hypotheses supported by markers, consistent with the project's computation-only constraint.
