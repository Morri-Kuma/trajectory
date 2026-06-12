# Three-Month Master's Research Project — Design Document

**Working title:**
**"Does how you label your cells change the trajectory you infer? An annotation-robustness study of single-cell reprogramming dynamics."**

**Author:** Kuma (Graduate School of Frontier Sciences, The University of Tokyo)
**Date:** 2026-06-09
**Status:** design proposal built on the existing `trajectory` benchmark framework

---

## 1. One-paragraph summary

Single-cell trajectory inference is only as trustworthy as the cell-state labels
it is evaluated against. Two of this project's three core datasets (GSE178325,
GSE230659) have **no author cell-type labels**, so the existing benchmark scores
trajectories against *silver-standard* labels it generated itself with a marker
pipeline. That makes the lineage results potentially circular. This project
turns that weakness into the research question: it adds a second, independent
annotation route — **scANVI/scArches reference-mapping from the labelled
GSE242424 dataset** — and asks, rigorously, **whether and where the choice of
annotation strategy changes the inferred trajectory and its biological
interpretation.** It is a self-contained, three-month, computation-only slice of
the larger benchmark, finishable on a laptop for development and the Shirokane
server for the full run, and strong enough to support a five-figure manuscript.

---

## 2. Central research question

> **In human somatic-cell reprogramming, does the cell-annotation strategy
> (dataset-specific marker "silver" labels vs. cross-dataset scANVI
> reference-mapping) materially change (a) the cell-state composition, (b) the
> inferred pseudotime/lineage topology, and (c) the biological conclusions —
> and if so, where and by how much?**

Sub-questions:
1. **Concordance.** How well do marker-silver labels and scANVI-transferred labels agree, and is disagreement random or concentrated at specific transition states?
2. **Trajectory sensitivity.** Do pseudotime orderings and state-transition graphs (DPT/PAGA/Slingshot, plus the generative methods) change when the annotation underneath them changes?
3. **Benchmark sensitivity.** Do method *rankings* on Lineage Fidelity flip depending on which annotation is treated as ground truth?
4. **Cross-protocol transfer.** Does an OSKM (GSE242424) reference annotate chemical-reprogramming (GSE178325/GSE230659) cells coherently, and what does the correspondence reveal about shared vs protocol-specific intermediate states?

---

## 3. Background and rationale

**Trajectory inference and its dependence on labels.** Methods from classical
pseudotime (DPT, Slingshot, PAGA) to generative time-series models (PRESCIENT,
MIOFlow, scNODE) and optimal-transport fate mapping (Waddington-OT, CellRank2)
reconstruct how cells move through state space over time. Benchmarks of these
methods (dynverse; scTimeBench) score them against reference cell states. When
the reference states are themselves uncertain — as in label-free reprogramming
data — the benchmark inherits that uncertainty.

**The specific gap in this project.** GSE178325 (Guan 2022) and GSE230659
(Liuyang 2023) provide raw chemical-reprogramming time courses but **no
published per-cell labels** (confirmed from the companion GitHub repos, whose
terminal cluster-ID assignments are left as empty placeholders). The project
therefore built silver labels via a marker-seed → milestone → trajectory-aware
pipeline. These are reasonable but *one particular choice*. No independent,
learned annotation has ever been compared against them. Meanwhile GSE242424
provides the missing ingredient: a reprogramming time course **with** per-cell
author labels (ATAC→RNA transfer), usable as a reference.

**Why this is a real contribution, not a chore.** "Annotation choice changes
downstream results" is suspected but rarely quantified end-to-end on the same
data with a frozen benchmark. Doing it on reprogramming — where intermediate
states are biologically contested — is novel and useful, and it directly
hardens the existing benchmark by reporting which conclusions are
annotation-invariant.

---

## 4. Hypothesis / main objective

**Primary objective.** Quantify the effect of annotation strategy on trajectory
inference and on benchmark conclusions, across the three required datasets.

**Working hypotheses.**
- **H1 (terminal robustness):** Terminal/extreme states (start fibroblast/hADSC, terminal iPSC/hCiPS) are annotation-invariant; both strategies agree strongly there.
- **H2 (intermediate fragility):** Disagreement concentrates in *intermediate/transition* states, where marker scoring and reference-mapping diverge most.
- **H3 (topology stability, ordering shift):** The coarse lineage *topology* (a mostly linear start→intermediate→terminal progression with limited branching) is preserved across annotations, but pseudotime *ordering* and branch assignment of intermediate cells shift measurably.
- **H4 (ranking stability):** Generative-method *relative rankings* on forecast accuracy are largely annotation-independent (forecast does not use the labels), whereas *lineage-fidelity* rankings are more sensitive.

These are falsifiable and each maps to a figure (§13/figure_plan.md).

---

## 5. Dataset strategy

**Core (required, retained):** GSE178325 and GSE230659 as **queries**; GSE242424
as the **labelled reference**. See `datasets/dataset_inventory.md`.

**Honest mismatch handling.** The reference is OSKM, the queries are chemical.
Their label spaces correspond only partially (start/terminal yes, intermediates
partial). The project treats scANVI output as a **label-correspondence
hypothesis**, analysed with contingency/agreement tables and marker validation,
**not** as a gold accuracy target. This limitation is stated wherever results
are reported.

**Supporting datasets (graded use):**
- **GSE175634 (human iPSC→CM, author-labelled):** *positive control* — the one case where ground truth is known, so we can validate that scANVI reference-mapping and the agreement metrics behave correctly when a real answer exists.
- **GSE298212 (human blood chemical reprog):** *external validation* of cross-protocol transfer.
- **E-MTAB-7552 (human organoid, labelled):** *generalization* beyond reprogramming (organoid fallback tier).
- **GSE218855 (mouse FCR, ships `.h5ad`):** *smoke-test input* and cross-species robustness control only.

**Fallback logic (explicit, per brief).** Human iPSC *reprogramming* series with
per-cell labels are scarce, so additional datasets descend the priority ladder:
labelled human iPSC *differentiation* (GSE175634) → other human chemical
reprogramming (GSE298212) → human organoid (E-MTAB-7552) → cross-species control
(GSE218855). Each step is justified by scarcity at the tier above.

---

## 6. Methodological workflow

```text
            ┌─────────────────────────────────────────────────────────┐
            │ 0. Config-driven entry (config.yaml): test vs server mode │
            └─────────────────────────────────────────────────────────┘
                                      │
  ┌───────────────────────┬──────────┴───────────┬────────────────────────┐
  ▼                       ▼                      ▼                         ▼
1. Data load        2. Preprocess         3a. Annotation A:        3b. Annotation B:
(.h5ad via shared   (QC, normalize,       marker-silver            scANVI reference-map
 registry)           HVG, PCA, neighbours)(existing in-repo)       (GSE242424 → queries)
                                                │                         │
                                                └──────────┬──────────────┘
                                                           ▼
                                        4. Annotation comparison
                                        (composition, latent UMAP,
                                         contingency/agreement, markers)
                                                           │
                                                           ▼
                                        5. Trajectory inference ×2 annotations
                                        (DPT/PAGA/Slingshot + generative methods)
                                                           │
                                                           ▼
                                        6. Comparison & benchmark sensitivity
                                        (pseudotime corr, graph similarity,
                                         lineage-fidelity ranking deltas)
                                                           │
                                                           ▼
                                        7. Biological interpretation
                                        (markers, transition states, regulators)
                                                           │
                                                           ▼
                                        8. Figures + reports (manuscript)
```

Steps 1–2 and the generative methods (5) reuse the existing `benchmark/`
package; steps 3b, 4, 6 are the new contribution implemented under `src/`.

---

## 7. Cell annotation strategy

**Annotation A — marker "silver" (existing).** The in-repo pipeline:
marker-seed scoring → milestone-marker labels → Stage-1/Stage-2 trajectory-aware
resolution → `final_milestone_label_coarse`. Dataset-specific, interpretable,
already frozen as official providers. Used as-is.

**Annotation B — scANVI/scArches reference-mapping (new).**
1. Build a reference from **GSE242424** with author labels (`author_cluster_label` / coarse milestones).
2. Train **scVI** on the reference, then **scANVI** (semi-supervised) to make the latent space label-aware.
3. Map each query (GSE178325, GSE230659) onto the frozen reference via **scArches** surgery; predict per-cell labels **with uncertainty**.
4. Keep `scanvi_label`, `scanvi_confidence`, and `scanvi_unknown_flag` (low-confidence → "unknown_or_ood") in `obs`.

**Comparison layer.** Because label spaces differ, agreement is reported as a
**contingency matrix + normalized mutual information + per-state precision/recall
where states correspond**, plus a curated correspondence map (e.g.
fibroblast↔hADSCs, iPSC↔hCiPS) validated by canonical markers (POU5F1/SOX2/NANOG
for pluripotency; COL1A2/DCN for stromal/fibroblast). Full design in
`scanvi_annotation_branch.md`.

---

## 8. Trajectory inference strategy

Run trajectories **twice per dataset** — once under each annotation — so that any
difference is attributable to annotation alone (everything else held fixed by
config + seed).

- **Classical, annotation-sensitive:** DPT (already computed in-repo), PAGA (state-graph topology), Slingshot (lineages/pseudotime). These consume the labels directly and are the most sensitive probes of annotation effect.
- **Generative / OT (existing benchmark methods):** scNODE, PRESCIENT, MIOFlow (forecast + embedding + lineage), WOT and CellRank2 (lineage only). These are largely label-independent for *training* but label-dependent for *lineage-fidelity scoring*, which is exactly the sensitivity we want to measure.

Pseudotime is anchored with a fixed root (start population) per dataset to make
runs comparable.

---

## 9. Benchmarking / comparison strategy

Two comparison axes, both reusing the scTimeBench-aligned metrics already in
`benchmark/evaluation/`:

1. **Annotation-vs-annotation (new).** On identical inputs: composition shift, label agreement (contingency/NMI), pseudotime correlation (Spearman of per-cell DPT order under A vs B), and PAGA/lineage graph similarity under A vs B.
2. **Method benchmark sensitivity (extends existing).** Re-score Lineage Fidelity for each method under annotation A and under annotation B; report ranking deltas (does any method's rank flip?). Forecast Accuracy serves as an annotation-invariant control.

A result is reported as **robust** if it holds under both annotations, and
**annotation-dependent** otherwise — this labelling is the manuscript's main
takeaway and a concrete service to the benchmark.

---

## 10. Biological interpretation strategy

- **Marker validation** of corresponding states (pluripotency, stromal/fibroblast, epithelial/MET, XEN-like) to judge which annotation places transition cells more defensibly.
- **Transition-state focus:** characterise the cells where A and B disagree — are they genuine intermediates (co-expressing programs), low-quality, or doublets?
- **Regulators/pathways (full version):** for the most annotation-sensitive transition (e.g. intermediate_plastic ↔ Partially-reprogrammed/Pre-iPSC), report differentially expressed TFs/pathways to give the disagreement a biological reading (MET, pluripotency activation).
- **Interpretation guardrail:** biological claims are made **only** for analyses actually run on full data (server); laptop runs prove the pipeline, not the biology.

---

## 11. Expected results

- A concordance map showing **high terminal agreement, lower intermediate agreement** (supports H1/H2).
- Latent-space UMAPs where query cells land coherently on the reference at the extremes and diffusely at intermediates.
- Pseudotime that is **highly correlated** (likely Spearman ≳ 0.8) between annotations along the main axis, with localised reordering at branch points (supports H3).
- A method-ranking table that is **stable for forecast** and **partly sensitive for lineage fidelity** (supports H4).
- A short list of transition states whose biological identity is annotation-dependent — the most novel, discussion-worthy result.

(These are *hypothesised* outcomes to be confirmed by the full run; none is a claimed finding yet.)

---

## 12. Risk management

| Risk | Likelihood | Mitigation |
|---|---|---|
| scANVI cross-protocol transfer is poor (OSKM→chemical) | Medium | Treat as a *finding*, not failure; validate with GSE175634 positive control where ground truth exists; report correspondence honestly |
| scvi-tools/GPU not available locally | High (laptop) | `test_mode` runs a fast KNN-on-PCA surrogate classifier + tiny CPU scANVI smoke; full scANVI runs on Shirokane |
| Compute/storage blow-up | Medium | HVG2000 inputs, `max_cells_per_dataset`, subsampling, reuse existing `.h5ad`; never re-download raw |
| Label-space mismatch confuses metrics | Medium | Use correspondence-aware metrics (contingency/NMI), not naive accuracy |
| Scope creep back to full benchmark | Medium | MVP gate (§14); generative-method re-runs are "full version" only |
| Batch effects masquerade as annotation effects | Medium | Fix preprocessing + seeds across A/B; integrate on reference; report batch diagnostics |

---

## 13. Three-month timeline

Assumes ~1 FTE master's student; laptop for development, Shirokane for full runs.

| Weeks | Phase | Deliverable |
|---|---|---|
| 1–2 | Setup & reproduce | `config.yaml`, `src/` scaffold, small-subset run of existing pipeline; reproduce one marker-silver provider |
| 3–4 | scANVI reference | Trained scVI/scANVI reference on GSE242424; reference UMAP + label recovery sanity check (Fig 3a) |
| 5–6 | Query mapping | scArches mapping of GSE178325/GSE230659; predicted labels + uncertainty; **annotation comparison** (Figs 2, 3) |
| 7–8 | Trajectory ×2 | DPT/PAGA/Slingshot under both annotations; pseudotime + graph comparison (Fig 4) |
| 9 | Benchmark sensitivity | Re-score lineage fidelity under A/B for existing methods; ranking deltas |
| 10 | Biology | Marker/transition-state interpretation; positive control (GSE175634); external validation (GSE298212) (Fig 5) |
| 11 | Full server run | Scale to full cells on Shirokane; finalize figures |
| 12 | Write-up | Manuscript draft, methods, reproducibility check |

Buffer is built into weeks 9–10; if scANVI transfer is weak, weeks 10–11 pivot
to the GSE175634 positive-control story without changing the thesis.

---

## 14. Minimum viable version (MVP)

The thesis stands on the **two required query datasets + the reference**, three
classical trajectory methods, and the annotation comparison:

- GSE242424 scANVI reference + scArches mapping of GSE178325 & GSE230659.
- Annotation A vs B comparison (composition, latent UMAP, contingency/NMI, markers).
- DPT + PAGA under both annotations; pseudotime correlation + graph similarity.
- Figures 1–4 + a reduced Figure 5 (markers only).
- One full-data server run for the two queries.

This MVP is achievable even if the generative-method re-scoring slips, and is
already a coherent, defensible MSc thesis.

---

## 15. Full manuscript-level version

Adds, as time allows:
- Slingshot + the generative/OT methods (scNODE, PRESCIENT, MIOFlow, WOT, CellRank2) re-scored under both annotations → method-ranking sensitivity (optional Fig 6).
- GSE175634 positive control and GSE298212 external validation as full panels.
- Organoid generalization (E-MTAB-7552) and/or scFM-embedding (scGPT/Geneformer) annotation as a robustness extension.
- Regulator/pathway analysis of annotation-sensitive transitions.

---

## 16. Reproducibility plan

- **Single entry point:** `config.yaml` with `mode: test|server`, dataset/output paths, seeds, model + trajectory params (see file).
- **Modular `src/`:** `data/`, `preprocessing/`, `annotation/`, `trajectory/`, `evaluation/`, `plotting/`, `utils/` — each step a function callable from a notebook or CLI, no hidden state.
- **Deterministic:** fixed `random_seed` threaded through numpy/torch/scanpy; versions pinned in an environment file.
- **Provenance:** every output `.h5ad`/table records the config hash, annotation strategy, and dataset id in `uns`/metadata.
- **Notebooks** (`notebooks/01…05`) are thin drivers over `src/` so an examiner can re-run the MVP on a laptop subset, then the same code scales on the server by switching `mode`.
- **No fabricated artifacts:** results files are produced only by code in the repo; documents separate *confirmed* (run) from *expected* (hypothesised) results.

---

## 17. Relationship to the existing repository

Nothing is thrown away. The existing `benchmark/` becomes the **"full server
version"** engine; this project defines the **focused, finishable MSc slice**
and adds the missing **annotation branch** plus a **reproducible `src/`/config
front-end**. The redesign reuses ~80% of existing assets (inputs, providers,
evaluators, methods) and fixes the project's single biggest scientific
vulnerability — circular ground truth — by measuring, rather than assuming,
annotation robustness.
