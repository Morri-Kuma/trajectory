# Project Audit Report — `trajectory`

**Date:** 2026-06-09
**Auditor:** research-engineering pass (read-only inspection; no files deleted)
**Scope:** full inspection of `C:\Users\37620\trajectory`, with emphasis on scripts, notebooks, documentation, data-processing logic, results, the `reference/` folder, and all assets tied to GSE178325, GSE230659, GSE242424.

> This audit is deliberately critical. The goal is not to praise the existing
> work but to identify what is solid, what is fragile, and what should change to
> turn the repository into a coherent, feasible, manuscript-ready three-month
> master's project. Confirmed facts are separated from assumptions throughout.

---

## 0. Executive summary

The repository is **far more mature than a typical "trajectory inference" starter project**. It is a working, scTimeBench-aligned **benchmarking framework** for human chemical-reprogramming trajectories, with five inference methods, two single-cell foundation-model (scFM) representation pipelines, frozen ground-truth/silver-standard providers, cluster job scripts, and aggregated result reports. The engineering quality of `benchmark/` is high.

The central tension this audit surfaces: **the existing project is scoped like a multi-year benchmarking paper, not a three-month MSc project.** It spans 36 GB of data and 89 GB of benchmark artifacts, six evaluation scenarios, five methods, and two scFMs. A single master's student cannot defensibly *complete and write up* all of that in three months. The requested **scANVI annotation-robustness branch is therefore not a side quest — it is the opportunity to carve a clean, self-contained, publishable slice** out of the larger machine: *does the choice of cell-annotation strategy change the trajectory you infer?*

One biological nuance must be stated up front because it shapes the whole annotation branch: **the reference dataset (GSE242424) is OSKM transcription-factor reprogramming of fibroblasts, whereas the two query datasets (GSE178325, GSE230659) are chemical reprogramming of hADSCs.** Their label spaces are related but not identical. This makes the scANVI reference-mapping comparison scientifically interesting *and* imposes a real limitation that must be handled with a label-correspondence analysis rather than a naive accuracy score.

---

## 1. Current project structure

Confirmed top-level layout (git-tracked unless noted):

```text
trajectory/
├── README.md                    # scTimeBench-aligned benchmark overview (good, current)
├── .gitignore                   # excludes data/, models/, logs/, *.h5ad, heavy intermediates
├── benchmark/                   # the real engine (89 GB on disk; code is light, artifacts heavy)
│   ├── adapters/                # method wrappers: wot, cellrank2, prescient, mioflow, scnode, future_model
│   ├── annotation/              # marker-seed + milestone + trajectory-aware silver labelling
│   ├── configs/                 # method × scenario × dataset YAMLs (+ runtime/, representation/)
│   ├── evaluation/              # forecast / embedding / lineage / representation evaluators
│   ├── ground_truth/            # frozen providers + registry.yaml  (+ _expanded_official_backup/)
│   ├── inputs/                  # benchmark-ready .h5ad (HVG2000 + full-gene) + representation X_rep
│   ├── methods/                 # vendored MIOFlow / PRESCIENT / WOT / scNODE
│   ├── representations/         # scGPT / Geneformer embedding extraction
│   ├── reports/                 # official_silver/, gse242424_oskm_ground_truth/, figures/, qc/
│   ├── results/                 # per-method per-scenario outputs (+ smoke/, _invalid_lineage_label_space/)
│   └── shared/                  # dataset registry + constants
├── data/                        # NOT git-tracked; 36 GB raw + processed (see §6)
│   ├── gse175634/               # iPSC→cardiomyocyte (cardiac), partial local metadata
│   ├── gse178325_human/         # Guan 2022 chemical reprogramming (rna_seq = older nested dup; rna_seq_10x canonical)
│   ├── gse230659(human/         # Liuyang 2023 chemical reprogramming  ← literal "(" in folder name
│   ├── gse242424/               # OSKM reprogramming (GSE242423 gene file)
│   └── processed/               # 3 full-gene benchmark .h5ad (2.4–2.6 GB each)
├── docs/                        # framework_v2 design doc, meetings, method_notes, reports
├── jobs/shirokane/              # ~45 qsub run/submit scripts for the Shirokane cluster
├── logs/                        # NOT git-tracked; hundreds of run logs + a nested .git repo
├── models/                      # NOT git-tracked; Geneformer-V2-104M + scGPT_whole_human (607 MB)
├── reference/                   # papers (PDFs) + paper_references.md (18 refs)
└── scripts/                     # standalone dataset builders, validators, plotting
```

The design source of truth is `docs/framework/experimental_framework_v2.md` (1,223 lines): three evaluation dimensions (Forecast Accuracy, Embedding Coherence, Lineage Fidelity), six scenarios (A–C observed-time, D–F pseudotime, the latter on hold), method capability gating, and frozen ground-truth provider strategy.

**Structural assessment.** `benchmark/` is well-modularised and the separation of adapters / evaluation / ground_truth / reports is genuinely good. The weaknesses are at the *edges*: raw `data/` is messy and inconsistently named, `logs/` is enormous and contains a stray nested git repo, `benchmark/archive/` alone is 32 GB, and there is **no top-level reproducible entry point, no `config.yaml`, no `TODO.md`, no `src/` for the lighter analysis code** — the analysis logic lives inside the benchmark package and in one-off `scripts/`. For a master's project that must be re-runnable by an examiner, that is a gap.

---

## 2. Main research question implied by the existing files

The files imply a **methods-benchmarking** question, not a single biological one:

> *Among modern single-cell time-series trajectory/forecasting methods (WOT, CellRank2, PRESCIENT, MIOFlow, scNODE) — and modern scFM embeddings (scGPT, Geneformer) used as their input representation — which most faithfully reconstruct the dynamics of human somatic-cell reprogramming, judged by forecast accuracy, embedding coherence, and lineage fidelity against frozen silver-standard / ground-truth state graphs?**

This is coherent and aligned with scTimeBench. But it is **broad**, and its weakest link is the *ground truth*: GSE178325 and GSE230659 have **no author per-cell labels** (confirmed in the cell-annotation READMEs and the raw-data audit), so the project had to **manufacture** silver-standard milestone labels via a marker-seed → milestone → trajectory-aware pipeline. The whole lineage-fidelity result therefore rests on the quality and arbitrariness of that home-grown annotation. **That dependency is exactly what the scANVI branch should stress-test**, which is why the requested branch is the most valuable next step rather than adding a sixth method.

---

## 3. Current datasets used

| Accession | System | Role in repo | Per-cell labels? | Benchmark input present |
|---|---|---|---|---|
| **GSE178325** (Guan 2022, *Nature*) | Human chemical reprogramming, hADSC/HEF/ASF → hCiPS | Primary query / silver-standard | **No author labels** → silver milestones built in-repo (`final_milestone_label_coarse`, 5 states + ambiguous) | `…/gse178325_marker_fm_transition_silver_hvg2000/…h5ad` (80,475 cells, HVG2000) + full-gene |
| **GSE230659** (Liuyang 2023, *Cell Stem Cell*) | Optimized human chemical reprogramming, hADSC → hCiPS | Primary query / silver-standard | **No author labels** → silver milestones (4 states + ambiguous) | `…/gse230659_marker_fm_transition_silver_hvg2000/…h5ad` (75,194 cells) + full-gene + scGPT/Geneformer X_rep |
| **GSE242424** (OSKM reprogramming) | Human fibroblast → iPSC, **TF (OSKM)** reprogramming | **Reference / ground truth** | **Yes** — author ATAC-cluster labels transferred to RNA (`author_cluster_label` 9 states; `final_milestone_label_coarse` 10 states) | `…/gse242424_author_cluster_matched/…h5ad` (59,187 cells) |
| GSE175634 (cardiac) | Human iPSC → cardiomyocyte differentiation | Secondary benchmark dataset | Author `type` labels used directly (6 states) | `…/gse175634_cardiac_author_hvg2000/…h5ad` (292 MB) |

Additional datasets are physically present under `data/` but **not** wired into the formal benchmark (from the raw-data audit v2): GSE298212 (human blood-cell chemical reprogramming), GSE218855 / FCR (mouse fast chemical reprogramming, `.h5ad` available), GSE247600 (iPSC→dopaminergic neuron, very few cells), GSE280956 (ChIP-seq, not expression), GSE136314 (empty locally).

**Key confirmed fact for the annotation branch:** GSE242424's reference labels come from the kundajelab/scATAC-reprog ATAC→RNA transfer table and describe an **OSKM** system (Fibroblast, Keratinocyte-like, Intermediate, Partially-reprogrammed, Pre-iPSC, iPSC, hOSK, xOSK, Fibroblast-like). The query datasets' silver labels describe a **chemical** system (hADSCs, epithelial_like, intermediate_plastic, xen_like, hCiPS). The terminal (iPSC ≈ hCiPS) and start (Fibroblast ≈ hADSCs) states correspond; the intermediates only partially correspond. This must be modelled as **label correspondence**, not identity.

---

## 4. Current computational methods used

**Trajectory / forecasting methods (5), capability-gated** (`benchmark/configs/method_capabilities.yaml`):

| Method | Forecast | Embedding | Lineage | Notes (confirmed from configs/adapters) |
|---|:--:|:--:|:--:|---|
| scNODE | ✅ | ✅ | ✅ | Generative neural-ODE; runs in HVG/PCA space |
| PRESCIENT | ✅ | ✅ | ✅ | Generative; trained in PCA space, inverse-transformed |
| MIOFlow | ✅ | ✅ | ✅ | Neural-ODE / OT flow in PCA space |
| WOT | ❌ | ❌ | ✅ | Optimal transport; lineage-only by design |
| CellRank2 | ❌ | ❌ | ✅ | Fate mapping; lineage-only at current stage |

**Representations:** scGPT (`scgpt_whole_human`) and Geneformer-V2-104M zero-shot embeddings (CLS / PCA50) built as `X_rep` inputs.
**Annotation:** in-repo marker-seed → milestone-marker → Stage-1/Stage-2 trajectory-aware silver labelling (no learned reference-mapping method is currently used — this is the gap the branch fills).
**Evaluation:** forecast (exact/representation), embedding coherence (incl. official-silver label diagnostics), lineage graph-similarity, plus QC and visual-summary scripts.
**Infrastructure:** Scanpy/AnnData-style `.h5ad`, YAML-driven configs, Shirokane qsub array jobs.

**Assessment.** The method roster is appropriate and current (all five are real, peer-reviewed/2021–2024). The conspicuous **methodological gap is annotation**: there is no learned, transferable cell-type classifier anywhere in the pipeline — exactly what scANVI/scArches provides. Adding it is well-justified and not redundant.

---

## 5. Current weaknesses and gaps

1. **Ground truth is self-referential (highest-priority scientific risk).** Two of three primary datasets have no author labels; lineage fidelity is scored against silver labels the project invented. Without testing annotation robustness, a reviewer can dismiss the lineage results as circular. → *scANVI branch directly addresses this.*
2. **Scope is too large for three months.** Five methods × six scenarios × multiple datasets × two scFMs is a PhD-scale benchmark. The MSc deliverable needs a defensible *slice* with a single clear question.
3. **No reproducible top-level entry point.** No `config.yaml`, no `src/`, no `notebooks/`, no `TODO.md`. Re-running the project requires tribal knowledge of `benchmark/` internals and cluster scripts. An examiner cannot reproduce it locally.
4. **Local-vs-server split is implicit.** Everything assumes Shirokane. There is no documented `test_mode` / `max_cells` / `use_small_subset` path for a Windows laptop, which is the stated execution constraint.
5. **Raw `data/` hygiene.** Inconsistent naming (`gse230659(human` with a literal parenthesis; `gse242424` holding a `GSE242423` gene file), a duplicated `gse178325_human/rna_seq` copy (the flattened `rna_seq_10x` is canonical and used by the builder), and mixed conventions. The parenthesis in particular breaks naive shell/glob code.
6. **Heavy, partly-orphaned artifacts.** `benchmark/archive/` = 32 GB; `logs/` holds hundreds of `.log` files **and a stray nested `.git`**; `benchmark/results/{smoke,_invalid_lineage_label_space}` and `ground_truth/_expanded_official_backup` are non-reportable. These bloat the working copy and confuse "what is current."
7. **Pseudotime scenarios (D–F) are designed but inactive.** Documented as "on hold." Fine, but the README/framework still advertise six scenarios, which overstates what is actually delivered.
8. **Reference list is solid but not yet trajectory-broad.** `reference/paper_references.md` has 18 well-chosen refs (datasets + the 5 methods + scTimeBench + Scanpy/AnnData), but lacks general trajectory-inference benchmarking context (e.g. dynverse), scANVI/scArches, and organoid time-series references that the redesign needs.
9. **No annotation-uncertainty propagation.** "ambiguous"/"unknown_or_ood" cells are excluded from metrics, but the downstream sensitivity of trajectories to that exclusion is untested.

---

## 6. Files that appear redundant, broken, obsolete, or unnecessary

**No file is deleted in this audit.** These are *candidates* to be itemised in `deletion_proposal.md` for explicit approval. Sizes confirmed via `du`/`ls`.

| Path | Category | Observation | Suggested action |
|---|---|---|---|
| `benchmark/archive/` (**32 GB**) | Heavy / obsolete | Largest single consumer of disk; `rerun_backups/` already gitignored | Archive off-repo (external drive / cold storage), then remove locally |
| `data/gse178325_human/rna_seq/` | Duplicate | Older nested layout; the flattened `rna_seq_10x/` is canonical (build script reads it, holds QC) and `rna_seq/` is unreferenced | Remove `rna_seq/`, keep `rna_seq_10x/` (verified 2026-06-09) |
| `logs/` (hundreds of `.log`) + `logs/.git/` | Transient / broken | Gitignored, but a **nested git repository** sitting inside `logs/` is almost certainly accidental | Remove `logs/.git`; rotate/compress old logs |
| `benchmark/ground_truth/_expanded_official_backup/` (38 MB) | Backup | Explicit backup of superseded providers | Keep one cold copy, drop from working tree |
| `benchmark/results/_invalid_lineage_label_space/` (372 KB) | Invalid | Name marks it as invalid output | Remove after confirming nothing references it |
| `benchmark/results/smoke/` (0 B) | Empty | Empty smoke dir | Remove |
| `data/trajectory_raw_data_audit_report_v2.md` vs `docs/reports/trajectory_raw_data_audit_report_v2.md` | Possible duplicate | The two files **differ** (not byte-identical) — one is likely stale | Diff, keep canonical copy under `docs/reports/`, archive the other |
| `scripts/__pycache__/`, `.pytest_cache/`, scattered `__pycache__/` | Transient | Build caches | Safe to remove (regenerated) |
| `data/gse230659(human/` (literal `(`) | Fragile naming | Not deletion — a *rename* candidate to `gse230659_human/` | Rename + update any path references |

**Net:** roughly **32+ GB** of local disk is recoverable with zero scientific loss, the bulk from `benchmark/archive/` and the redundant `rna_seq` copy. None of it is git-tracked, so cleanup is low-risk — but per project rules, nothing is removed without your approval.

---

## 7. What I recommend changing (preview of the redesign)

Spelled out fully in `research_design_three_month_project.md`, but in one line: **keep the benchmark engine as the "full server version," and define the MSc deliverable as a focused annotation-robustness study** — *current marker-silver annotation vs scANVI reference-mapping, and its effect on inferred trajectories* — across the three required datasets, with the other public datasets used for external validation. This reuses ~80% of what already exists, fixes the project's biggest scientific vulnerability (circular ground truth), and is finishable in three months on a laptop for development plus the Shirokane server for the full run.

**Nothing in this audit has modified or deleted any file.**
