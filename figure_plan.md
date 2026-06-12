# Manuscript Figure Plan

**Date:** 2026-06-09
**Project:** annotation-robustness of trajectory inference in human reprogramming
**Companion:** `research_design_three_month_project.md`, `scanvi_annotation_branch.md`

Each figure lists: the question it answers, panels, input data, method,
expected visual, the conclusion it supports, and whether it can be drafted in
**local test mode** (synthetic / small subset) or needs the **full server run**.
"Local-draft" means the *layout and code* are validated by the smoke test; the
*real numbers* always come from the server run.

---

## Figure 1 — Project design, datasets & workflow
- **Question:** What is being compared, on what data, and why?
- **Panels:**
  (a) schematic of the two annotation routes (marker-silver vs scANVI reference-mapping);
  (b) dataset table/timeline (GSE242424 reference + GSE178325/GSE230659 queries, time points);
  (c) analysis DAG (load → preprocess → annotate ×2 → trajectory ×2 → compare).
- **Input:** `datasets/dataset_inventory.*`, design docs (no cell data).
- **Method:** diagram (drawn), small metadata tables.
- **Expected visual:** one orienting schematic + dataset timeline.
- **Conclusion:** frames the central question and the OSKM↔chemical correspondence caveat.
- **Mode:** **Local-draft** (no heavy compute).

## Figure 2 — Integrated cell landscape across datasets
- **Question:** Do reference and query cells share a coherent manifold, and where do queries land?
- **Panels:**
  (a) joint scANVI-latent UMAP coloured by dataset;
  (b) by time point;
  (c) by reference label (reference cells) with query cells overlaid;
  (d) batch/QC diagnostic.
- **Input:** GSE242424 + GSE178325 + GSE230659 HVG2000 inputs.
- **Method:** scVI/scANVI latent → UMAP (server); PCA/UMAP stand-in locally (`plot_embedding`).
- **Expected visual:** integrated embedding showing start/terminal overlap and intermediate spread.
- **Conclusion:** queries are mappable onto the reference at the extremes (supports H1); intermediates are diffuse (motivates H2).
- **Mode:** **Server** for real UMAP; **Local-draft** layout via smoke test.

## Figure 3 — Annotation-strategy comparison (current vs scANVI)
- **Question:** How much, and where, do the two annotations disagree?
- **Panels:**
  (a) cell-type composition per time point, marker-silver vs scANVI (mapped);
  (b) contingency/agreement heatmap (marker-silver × scANVI) + NMI/ARI;
  (c) per-state correspondence recall bar plot;
  (d) confidence / fraction flagged `unknown_or_ood`.
- **Input:** query `.obs` with both annotations.
- **Method:** `evaluation/annotation_compare.py` (`composition_compare`, `agreement_matrix`).
- **Expected visual:** terminal/start states agree strongly; intermediates and chemical-specific `xen_like` fall to low agreement / OOD.
- **Conclusion:** quantifies annotation (dis)agreement and localises it to intermediate states (supports H1/H2).
- **Mode:** **Local-draft** (smoke test already emits all four panels); **Server** for real labels.

## Figure 4 — Trajectory & pseudotime comparison under each annotation
- **Question:** Does changing the annotation change the inferred trajectory?
- **Panels:**
  (a) pseudotime-vs-pseudotime scatter (annotation A root vs B root) with Spearman;
  (b) PAGA/state graphs side by side (A vs B) + edge-Jaccard;
  (c) pseudotime on the embedding under each annotation;
  (d) branch-assignment shift for intermediate cells.
- **Input:** preprocessed query + both annotations.
- **Method:** `trajectory/pseudotime.py` (DPT + PAGA server; geodesic/surrogate locally); `graph_similarity`.
- **Expected visual:** high global pseudotime correlation (≈0.8) with localised reordering; mostly-preserved topology with some edge differences.
- **Conclusion:** main trajectory axis is robust; intermediate ordering/branching is annotation-sensitive (supports H3).
- **Mode:** **Local-draft** (smoke test emits pseudotime scatter + graph compare); **Server** for DPT/PAGA + generative methods.

## Figure 5 — Biological interpretation of annotation-sensitive states
- **Question:** What is the biological identity of the cells where annotations disagree, and which annotation is better supported by markers?
- **Panels:**
  (a) marker dotplot/heatmap (pluripotency, stromal/fibroblast, MET, XEN) by state, both annotations;
  (b) marker-AUROC comparison (which annotation separates programs better);
  (c) characterisation of disagreement cells (marker co-expression, QC);
  (d) DE genes/TFs for the most sensitive transition (server).
- **Input:** query expression + both annotations + disagreement table.
- **Method:** `marker_validation`, `disagreement_cells`; DE via scanpy (server).
- **Expected visual:** transition cells co-express competing programs; a short list of annotation-dependent states.
- **Conclusion:** gives the disagreement a biological reading (MET / partial pluripotency) and flags which conclusions are annotation-robust.
- **Mode:** **Server** for DE; **Local-draft** for marker-AUROC and disagreement panels.

---

## Optional / supplementary figures

- **S1 — Dataset QC:** per-sample genes/UMI/mito, low-quality fractions (reuse `logs/*qc*` + `benchmark/reports/qc/`). *Local-draft.*
- **S2 — Method-benchmark sensitivity:** Lineage-Fidelity ranking of scNODE/PRESCIENT/MIOFlow/WOT/CellRank2 under annotation A vs B, with rank-delta; Forecast Accuracy as the invariant control (supports H4). *Server.*
- **S3 — Positive control (GSE175634):** the same reference-mapping where author labels exist, validating the machinery. *Server.*
- **S4 — External validation (GSE298212):** cross-protocol transfer is in-scope. The **E-MTAB-7552 organoid generalization is deferred** (decided 2026-06-09) — it needs its own organoid reference and is out of MSc core scope. *Server.*
- **S5 — Robustness/ablation:** confidence threshold `tau`, HVG count, reference label key (`author_cluster_label` vs coarse milestones), full-gene vs HVG2000. *Server.*
- **S6 — scFM-embedding annotation:** scGPT/Geneformer latent as an alternative annotation substrate (reuses existing representation inputs). *Server.*

---

## Coverage vs required figure set

| Required (brief) | This plan |
|---|---|
| Fig 1: design, datasets, workflow | Figure 1 |
| Fig 2: integrated cell landscape | Figure 2 |
| Fig 3: annotation strategies incl. scANVI | Figure 3 |
| Fig 4: trajectory & pseudotime comparison | Figure 4 |
| Fig 5: biological interpretation | Figure 5 |
| Optional (QC, benchmarking, robustness, generalization, ablation) | S1–S6 |

All five main figures map onto code that the local smoke test already exercises
for layout; only the real data values require the server run. No figure depends
on fabricated data.

---

## Revised figure set (2026-06-09, after the first scANVI server run)

The first cross-protocol run produced a confident **collapse** (see
`docs/scanvi_run_results_interpretation.md` and `manuscript/manuscript_draft.md`),
reframing the paper from "annotation changes the trajectory" to **"when is
reference-mapping annotation trustworthy for reprogramming, and how do you tell?"**
Updated main figures, grounded in the real run:

- **Fig 1 — Design & datasets:** the two annotation routes, the three datasets, the OSKM↔chemical correspondence caveat. *(in hand)*
- **Fig 2 — Reference learnable, query off-manifold:** scANVI reference self-recovery = 0.990; joint latent UMAP with chemical query cells landing off the OSKM manifold at the intermediates. *(server: recovery in hand; UMAP pending)*
- **Fig 3 — The confident collapse (core result):** (a) scANVI label distribution = the `hOSK` stripe (96.3% / 76.6%); (b) agreement heatmap vs marker-silver; (c) scANVI confidence (≤0.3% OOD) vs (d) geometric-OOD distribution (query med ~18–20 vs ref ~5). *(in hand: `results/annotation_branch/*`, `results/test_outputs/ood_recalibration/*`)*
- **Fig 4 — OOD recalibration recovers the mismatch:** geometric OOD flags 38.7% / 43.6% of query cells vs scANVI ≤0.3%; threshold/percentile sweep; per-timepoint OOD (start cells in-distribution, reprogramming intermediates OOD). *(in hand locally; server scANVI-latent version via `scripts/score_ood.py`)*
- **Fig 5 — What works (validation):** positive control (GSE175634, labels known) + same-protocol transfer (GSE230659↔GSE178325) showing sensible recovery/agreement, isolating the failure to the cross-protocol gap. *(pending one GPU run — `docs/validation_runbook.md`)*

Supplements: dataset QC; correspondence map; ablations (k, OOD quantile, HVG count, reference label key); count-reconstruction validation (corr = 1.0000).

---

## Reviewer-r3 figure-logic consolidation (proposed layout — one conclusion per figure)

Current draft carries 8 result figures + 1 supplementary (figS1). Reviewer 3 asks that
each *main* figure serve one conclusion and that detailed tables move to the Supplement.
Proposed main-text set (5 figures), mapping existing assets:

| Main fig | Conclusion | Built from current assets |
|---|---|---|
| Fig 1 — Benchmark design | datasets, tasks, scenarios, reference/OOD pipeline schematic | new schematic (Table 1 + a flow diagram) |
| Fig 2 — Reference / OOD validation | the gate flags cross-modality transfers confidence misses | current fig5 (scANVI gate) + OOD k/q sweep panel |
| Fig 3 — Primary task results | no method dominates all three tasks | current fig1 (lineage) + fig2 (forecast) + fig3 (embedding) as a 3-panel |
| Fig 4 — Cross-system instability | lineage rankings do not transfer (Kendall's W) | current fig4 (cross-system) + fig7 (rank concordance) |
| Fig 5 — Seed robustness + metric validity | single-seed leaderboards are draws; metric is discriminating | current fig8 (multi-seed CI) + figS1 (negative control) |

Supplement: pseudotime panel (current fig6), the per-method/per-dataset tables
(Tables 3–4 detail), and the full negative-control table. This is a layout/montage pass
on existing PNGs; the underlying numbers and scripts are unchanged. Figure files and their
regeneration scripts are listed in `manuscript/REPRODUCIBILITY.md`.
