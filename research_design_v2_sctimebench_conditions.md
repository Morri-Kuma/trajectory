# Experimental Design (v2.1): A pluripotency-domain benchmark of temporal trajectory-inference methods, built on scTimeBench

**Status:** refinement grounded in `docs/framework/experimental_framework_v2.md`
(v2.1 scope update) and Osakwe, Huang & Li, *scTimeBench: A streamlined benchmarking
platform for single-cell time-series analysis* (bioRxiv 2026,
doi:10.64898/2026.03.16.712069). This document replaces the earlier conditions draft.
All result numbers are read from files in `benchmark/reports/`, `benchmark/results/`,
and `results/annotation_branch/`; nothing is fabricated.

**v2.1 changes (2026-06-10).** Three axes were **removed** from the design — input
representation (HVG vs scGPT vs Geneformer), annotation-provider-as-condition, and
time–library confounding. The scANVI/scArches work is **retained and integrated** as
the benchmark's reference-construction-and-validation engine. The dataset suite is
expanded to **six time-labeled iPSC scRNA-seq datasets** (plus an organoid extension).

**Status update (2026-06-11).** Pseudotime Scenarios **D/E/F are complete** for both
chemical primaries (GSE178325, GSE230659): run on Shirokane, integrated into the
official A–F rankings, and written into the manuscript. Headline result — pseudotime
denoises lineage recovery conditionally, most strikingly for PRESCIENT on GSE230659
(mean single-step lineage AUROC $0.64\to0.97$), reproducing scTimeBench's effect on
reprogramming data. The two integration datasets (GSE298212, GSE218855) are now fully
prepared locally — input builders, milestone-marker schemas (incl. a mouse set), configs,
and Shirokane jobs all exist; only raw-data staging and the actual runs remain (Shirokane).

---

## 1. Objective and relationship to scTimeBench

scTimeBench evaluates time-aware trajectory-inference methods on three core tasks —
Forecast Accuracy, Embedding Coherence, and Lineage Fidelity — across eight
developmental datasets spanning four species, and concludes that (i) strong
forecasters need not preserve biological signal, (ii) lineage fidelity is the central
weakness, with methods often no better than a Spearman-correlation baseline, and
(iii) pseudotime can denoise trajectories.

scTimeBench answers *"which method is best on average across developmental systems?"*.
We reuse its task definitions, metric families, scenario structure, and
method-eligibility rule **unchanged**, and apply them to a **single coherent
biological domain — induced pluripotency** — to answer a question scTimeBench does
not: **does a method's performance generalize across the different routes that
converge on the pluripotent state?** Chemical reprogramming, OSKM
transcription-factor reprogramming, iPSC-directed differentiation, mouse MEF→iPSC
reprogramming, and (as an extension) iPSC-derived organoid development all reach or
depart from pluripotency by mechanistically different paths. If method rankings are
stable across these routes, a practitioner can pick a method once; if they invert,
method choice must be conditioned on the system.

## 2. Two innovations distinguishing this benchmark from scTimeBench

**Innovation I — a pluripotency-domain cross-modality generalization benchmark.**
Rather than a single averaged leaderboard over heterogeneous developmental systems,
we hold the biology to one domain (iPSC) and make the **reprogramming/differentiation
modality** the primary axis, asking whether method performance *transfers* across
modalities that share the pluripotent endpoint. The only two remaining experimental
axes are the **biological system/modality** (six datasets) and the **forecasting
difficulty** (Scenarios A/B/C).

**Innovation II — a reference-construction-and-validation pipeline (the scANVI
integration).** scTimeBench sidesteps the annotation problem by curating
author-labeled datasets. But most iPSC time courses (e.g. both chemical-reprogramming
primaries) ship **no per-cell author labels**, and Embedding Coherence and Lineage
Fidelity both require a cell-state reference. We therefore build silver references by
**scANVI/scArches label transfer** from an annotated dataset and **gate every transfer
with a geometric out-of-distribution (OOD) score**: a transfer is accepted only if the
query cells lie inside the reference's latent support, and rejected otherwise. This
makes the benchmark extensible to unlabeled iPSC data while quantifying when a
transferred reference can be trusted — a methodological contribution scTimeBench lacks.

## 3. Core benchmark structure (inherited from scTimeBench, unchanged)

| Task | Metrics | Annotation-dependent? |
|---|---|---|
| Forecast Accuracy | Wasserstein, Gaussian MMD, Energy-MMD, Hausdorff (geomloss, recomputed centrally `metric_backend=scTimeBench_exact`) | No |
| Embedding Coherence | Leiden-cluster ARI vs reference labels; average normalized classifier entropy | Yes |
| Lineage Fidelity | single/multi-step (Floyd–Warshall) AUROC, AUPRC, Jaccard on a PR-thresholded state graph | Yes |

**Methods and capability gating** (scTimeBench rule: score only supported dimensions):
WOT and CellRank2 are optimal-transport, lineage-only; scNODE, PRESCIENT, and MIOFlow
are projection-capable and scored on all three tasks. CellRank2 runs as fate/lineage
post-processing on the WOT `RealTimeKernel` and is reported as such.

**Scenarios:** A observed-time interpolation, B extrapolation, C joint; metrics on
held-out time points. Pseudotime Scenarios D–F are specified but deferred (Section 8).

**Providers:** a single frozen `official_silver` provider per dataset, keyed on
`final_milestone_label_coarse`, with `ambiguous`/`unknown_or_ood` excluded; rankings
are computed within a `(dataset, scenario, result_class, provider, label_mode)` group
and a contract checker refuses to mix incompatible settings. *(The multi-provider
"annotation as a condition" comparison has been removed — Section 1 of the framework
v2.1 note.)*

## 4. The two experimental axes

**Axis 1 — biological system / modality (six datasets).** Chemical reprogramming
(GSE178325, GSE230659), OSKM reprogramming (GSE242424), iPSC-derived differentiation
(GSE175634), human blood-cell chemical reprogramming (GSE298212), and mouse MEF→iPSC
reprogramming (GSE218855), with a human cerebral-organoid time course (E-MTAB-7552) as
the stated extension. Details and run status: `datasets/dataset_inventory.{md,csv}`.

**Axis 2 — forecasting difficulty (Scenarios A/B/C).** Interpolation vs extrapolation
vs joint, the within-dataset stress test.

## 5. Reference construction and validation (the scANVI pipeline)

For datasets without author labels, silver references are produced and validated by a
three-gate procedure, using GSE242424 (author-labeled OSKM) and GSE175634 (author
`type` labels) as candidate reference atlases:

1. **Self-recovery gate.** scANVI must re-predict the reference's own author labels at
   ≥0.90 held-out accuracy before it may annotate any query. On GSE242424 this gate is
   met (recovery **0.990**).
2. **Geometric OOD gate.** Each query cell's mean distance to its 15 nearest reference
   cells in the shared scANVI latent is compared to the reference's own kNN-distance
   distribution; a transfer is rejected if a large fraction of query cells fall beyond
   the reference's 95th percentile.
3. **Within-modality rule.** Transfers are accepted only within a reprogramming
   modality.

This procedure is justified by a real, completed diagnostic. Mapping the **chemical**
queries onto the **OSKM** reference (a cross-modality transfer) **collapses**: 76–96%
of query cells are confidently assigned to a single OSKM-transgene state (`hOSK`,
<1% of the reference), agreement with the marker-silver labels is 0.01–0.04, and
scANVI's posterior confidence flags ≤0.3% of cells — yet the geometric OOD score flags
**75.9% / 73.0%** of query cells as out-of-distribution. Two controls confirm the gate
is calibrated: scANVI self-transfer on a labeled dataset recovers held-out labels at
**0.944**, and **within-modality** (chemical→chemical) transfer agrees with the
existing annotation at **0.76–0.86**. The benchmark therefore uses scANVI references
only for within-modality, OOD-passing transfers; cross-modality references remain
marker-defined. *(Files: `results/annotation_branch/`,
`results/positive_control_gse175634/`, `results/same_protocol_*/`,
`results/figures/ood_confidence_vs_geometric.png`.)*

This is the logical home of the scANVI work: not a separate paper, but the benchmark's
annotation engine plus a quantitative trust gate that lets Axis 1 scale to unlabeled
iPSC datasets.

## 6. Current real-result snapshot (honest)

**Primary chemical datasets (GSE178325 + GSE230659), projection-capable methods.**
Combined embedding+lineage rank is narrow: MIOFlow 1st (1.79), PRESCIENT 2nd (1.85),
scNODE 3rd (2.18). By task, PRESCIENT leads Embedding Coherence (mean ARI 0.135 vs
0.072, 0.054); MIOFlow leads Lineage Fidelity (mean single-step AUROC 0.662, ahead of
PRESCIENT 0.630 and scNODE 0.597, with WOT 0.560 and CellRank2 0.574 last); Forecast
has no single winner (scNODE wins Scenario A, PRESCIENT the harder B/C on GSE230659).

**Axis 2 (difficulty).** Scenario B degrades every method: forecast Wasserstein on
GSE230659 rises from ≈81 (scNODE A) to 148–404 (B), and the lowest lineage AUROC in
the primary benchmark is scNODE GSE230659 B = 0.436 (below chance). Absolute lineage
AUROC stays in a modest 0.44–0.77 band — scTimeBench's "lineage is hardest" result,
reproduced on reprogramming data.

**Axis 1 (system).** The same methods score markedly higher on OSKM (GSE242424:
single-step AUROC up to 0.85; method means 0.71–0.78) than on chemical GSE230659
(means 0.59–0.66): lineage performance is as much a property of the system and its
reference graph as of the method, and method ranking does **not** cleanly transfer
across modalities.

**Baseline.** Per scTimeBench and the project record, a Spearman-correlation baseline
is competitive; the current graph-sim rerun did not recompute it (no AnnData passed),
so a per-scenario graph-sim baseline is a required pre-submission item (Section 8).

## 7. Threats to validity and limitations

Silver references are marker-defined (or scANVI-transferred), not author-validated, so
lineage fidelity measures alignment with a coarse reference graph rather than
ground-truth biology; four of six datasets carry completed runs while GSE298212/
GSE218855 are integration-ready and the organoid extension needs its own reference;
PRESCIENT formal runs were CPU-only; and no multi-seed/bootstrap uncertainty is
attached yet, so small rank gaps (1.79 vs 1.85) are descriptive, not statistically
resolved.

## 8. Three-month plan

Month 1 — run GSE298212 and GSE218855 through the full A/B/C pipeline (build their
`official_silver` providers via the scANVI/OOD gate where labels are absent); recompute
the Spearman graph-sim baseline for every dataset/scenario; finish GSE175634 lineage.

Month 2 — add ≥3 seeds (or held-out-cell bootstrap) per projection method/scenario for
uncertainty; run the cross-modality generalization analysis (does ranking transfer?)
as the headline result; integrate the organoid extension (E-MTAB-7552) with its own
reference.

Month 3 — write-up (`manuscript/benchmark_manuscript.tex`) and a reviewer-facing
sensitivity appendix; package the cross-modality matrix and the scANVI/OOD
reference-validation figure as the headline set.

## 9. Reproducibility

One framework spec (`docs/framework/experimental_framework_v2.md` v2.1),
capability-gated dispatcher and per-task evaluators under `benchmark/evaluation/`,
frozen providers under `benchmark/ground_truth/providers/`, the scANVI/OOD utilities
under `src/annotation/` (unit-tested), explicit `result_class` separation of formal vs
validation runs, and Shirokane `qsub` launchers under `jobs/shirokane/`. Reports are
regenerated by `python -m benchmark.evaluation.summarize_official_silver`.
