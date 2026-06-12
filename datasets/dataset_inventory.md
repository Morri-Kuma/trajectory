# Dataset Inventory — iPSC trajectory benchmark (v2.1)

**Date:** 2026-06-10 (reframed to the scTimeBench-aligned conditions benchmark)
**Companion:** `datasets/dataset_inventory.csv` (machine-readable)

**Inclusion principle (per project brief).** The benchmark uses **at least six
time-labeled iPSC scRNA-seq datasets**. Priority order: (1) time-labeled human iPSC
reprogramming; (2) human iPSC-derived differentiation; (3) organoid time-series (the
allowed extension); (4) cross-species reprogramming control. The six-plus suite below
spans the **pluripotency domain across reprogramming modalities**, which is the axis
the benchmark differentiates methods over. Every accession is public and verified
against its source paper and/or locally inspected files; nothing is fabricated. Cell
counts marked as benchmark inputs are computed from the local HVG2000 `.h5ad`.

**Suite at a glance (≥6 satisfied).** Four datasets carry completed formal results
(GSE178325, GSE230659, GSE242424 fully; GSE175634 forecast+embedding); two more are
integration-ready with data already acquired (GSE298212 human blood reprogramming,
GSE218855 mouse reprogramming); the organoid extension (E-MTAB-7552) and a small
neural set (GSE247600) round out the modality coverage.

---

## Core iPSC suite

### GSE178325 — human chemical reprogramming (Guan 2022, *Nature*) · **primary**
Human hADSC/HEF/ASF → hCiPS, chemical. ~15 stage-day points (D0…StageIV D10, hCiPSC).
**80,475** cells (HVG2000 benchmark input). Labels: in-repo silver milestones (5 coarse
states). 10x MTX, integer raw counts. **Run status: formal results — forecast +
embedding + lineage (A/B/C).** Role: primary chemical-reprogramming time course.
URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178325

### GSE230659 — optimized human chemical reprogramming (Liuyang 2023, *Cell Stem Cell*) · **primary**
Human hADSC → hCiPS, optimized chemical. 15 stage-day points. **75,194** cells (HVG2000).
Labels: in-repo silver milestones (4 coarse states). 10x MTX, raw counts. **Run status:
formal results — forecast + embedding + lineage (A/B/C).** Role: cleanest
chemical-reprogramming course; near-linear trajectory.
URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230659

### GSE242424 — human OSKM (TF) reprogramming · **system condition + scANVI reference**
Human fibroblast → iPSC, **OSKM transcription-factor** reprogramming. D0…D14, iPSC (9).
**59,187** author-cluster-matched cells (156,969 full-gene). Labels: ✅ author ATAC→RNA
cluster labels (9 states) — the **only** dataset with per-cell author labels, so it is
the **scANVI/scArches reference** for label transfer *and* an OSKM system condition.
10x MTX, raw counts. **Run status: formal results — forecast + embedding + lineage
(A/B/C).** URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242424

### GSE175634 — human iPSC → cardiomyocyte differentiation (Elorbany 2022, *PLoS Genet*) · **system condition**
Human iPSC → cardiomyocyte, multi-line. 7 timepoints over ~16 days; author states
IPSC→MES→CMES→PROG→CM with a PROG→CF branch. **224,202** cells (HVG2000). Labels: ✅
author `type` labels (6 states); frozen provider `gse175634_cardiac_silver_v1`. **Run
status: partial — forecast + embedding done; graph-sim lineage pending.** Role:
labelled iPSC-differentiation system (a known branching reference) for cross-modality
generalization. URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175634

### GSE298212 — human blood-cell chemical reprogramming (Peng 2025, *Cell Stem Cell*) · **system condition**
Human PBMC/EPC → hCiPS, chemical. PBMC_EPC, S1D1/3/6/8 (5). **40,794** cells / 5 samples
(hg38, CellRanger 5.0.1). Labels: sample-level only; silver milestones to build. **Run
status: integration-ready — data acquired; configs to build; no runs yet.** Role:
second human chemical-reprogramming system (different starting cell type) for
cross-system generalization. URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298212

### GSE218855 — mouse fast chemical reprogramming (Chen 2023, *Nat Cell Biol*) · **cross-species condition**
Mouse MEF → iPSC, fast chemical reprogramming. D0/4/8/12, iPSC (5). **46,361** cells
(processed `.h5ad` present). Labels: sample/time/condition only; silver milestones to
build. **Run status: integration-ready — `.h5ad` present; no runs yet.** Role:
cross-species robustness (does method ranking survive a species change?).
URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218855

---

## Extensions (modality coverage, user-allowed)

### E-MTAB-7552 — human cerebral organoid development (Kanton 2019, *Nature*) · **organoid extension**
Human ESC (H9) + **iPSC (409b2)** → cerebral organoid; 0d/4d/10d/15d/1mo/2mo/4mo (7).
**~43,498** cells. Labels: ✅ author cell-type + developmental-stage. 10x (ArrayExpress).
**Run status: integration-ready — needs its own organoid reference; no runs yet.** Role:
the stated organoid extension and the longest developmental trajectory; tests whether
the framework generalizes past reprogramming.
URL: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-7552

### GSE247600 — human iPSC → midbrain dopaminergic neuron differentiation · **candidate (reserve)**
Human iPSC-derived DA neuron differentiation. D21/30/45/65 (4). **~1,022** cells (small).
Labels: none per-cell. Reserve only (low cell count, Ensembl-only IDs).
URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247600

---

## Additional verified candidate (not yet acquired)

- **GSE115943 — mouse MEF→iPSC, Waddington-OT (Schiebinger 2019, *Cell*).** The
  canonical dense reprogramming time course (~315k cells, sampled every 12 h over 18
  days) and the dataset WOT was designed on; a natural cross-species reprogramming
  reference that ties directly to the WOT method in our roster. Real and public
  (verified); acquisition + config is a future-work item.
  URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115943

---

## Summary table

| Accession | Species | System | Modality | Time pts | Cells | Labels | Role | Run status |
|---|---|---|---|---|---|---|---|---|
| GSE178325 | Human | hADSC→hCiPS | chemical reprog | ~15 | 80,475 | silver | primary | results (all 3) |
| GSE230659 | Human | hADSC→hCiPS | chemical reprog | 15 | 75,194 | silver | primary | results (all 3) |
| GSE242424 | Human | fib→iPSC | OSKM reprog | 9 | 59,187 | **author** | system + scANVI ref | results (all 3) |
| GSE175634 | Human | iPSC→CM | differentiation | 7 | 224,202 | author | system | forecast+embed |
| GSE298212 | Human | blood→hCiPS | chemical reprog | 5 | 40,794 | sample | system | integration-ready |
| GSE218855 | Mouse | MEF→iPSC | chemical reprog | 5 | 46,361 | sample | cross-species | integration-ready |
| E-MTAB-7552 | Human | iPSC/ESC organoid | organoid | 7 | ~43,498 | author | extension | integration-ready |
| GSE247600 | Human | iPSC→DA neuron | differentiation | 4 | ~1,022 | none | reserve | candidate |

Six iPSC datasets (rows 1–6) form the benchmark suite; E-MTAB-7552 adds the organoid
extension. All accessions are public and verified.
