# Literature Review (Extended) — trajectory project

**Date:** 2026-06-09
**Purpose:** consolidate the project's existing references with newly added,
independently verified references that support the redesigned three-month
project (trajectory inference + the scANVI annotation-robustness branch).

**Verification policy.** Every newly added reference below was checked via at
least one reliable public source (publisher page, PubMed/PMC, or DOI resolver)
during this pass. DOIs were confirmed to resolve to the stated title, authors,
journal, and year. Items I could not fully verify are explicitly marked
**UNCERTAIN** and are *not* relied on as core references. No reference here is
fabricated; where only the paper (not a data accession) was verified, that is
stated.

---

## Part A — Existing project references (carried over)

These come from `reference/paper_references.md` and remain valid. Summarised
here so this file is self-contained; see that file for the annotated original.

**Datasets / biology**
- [A1] Guan J, Wang G, Wang J, *et al.* Chemical reprogramming of human somatic cells to pluripotent stem cells. *Nature* 605, 325–331 (2022). DOI: 10.1038/s41586-022-04593-5. — **GSE178325** source.
- [A2] Liuyang S, Wang G, Wang Y, *et al.* Highly efficient and rapid generation of human pluripotent stem cells by chemical reprogramming. *Cell Stem Cell* 30, 450–459.e9 (2023). DOI: 10.1016/j.stem.2023.02.008. — **GSE230659** source.
- [A11] Chen X, *et al.* A fast chemical reprogramming system promotes cell identity transition through a diapause-like state. *Nat Cell Biol* 25, 1146–1156 (2023). DOI: 10.1038/s41556-023-01193-x. — FCR mouse background.
- [A12] Peng F, *et al.* Chemical reprogramming of human blood cells to pluripotent stem cells. *Cell Stem Cell* 32, 1192–1199.e11 (2025). DOI: 10.1016/j.stem.2025.07.003. — GSE298212 background.
- [A13] Wang Y, *et al.* A rapid chemical reprogramming system to generate human pluripotent stem cells. *Nat Chem Biol* 21, 1030–1038 (2025). DOI: 10.1038/s41589-024-01799-8.

**Methods (trajectory / forecasting)**
- [A5] Zhang J, *et al.* scNODE: generative model for temporal single-cell transcriptomic data prediction. *Bioinformatics* 40(Suppl_2), ii146–ii154 (2024). DOI: 10.1093/bioinformatics/btae393.
- [A6] Huguet G, *et al.* Manifold interpolating optimal-transport flows for trajectory inference (MIOFlow). *NeurIPS 35*, 29705–29718 (2022). arXiv: 10.48550/arXiv.2206.14928.
- [A7] Yeo GHT, Saksena SD, Gifford DK. Generative modeling of single-cell time series with PRESCIENT. *Nat Commun* 12, 3222 (2021). DOI: 10.1038/s41467-021-23518-w.
- [A8] Schiebinger G, *et al.* Optimal-transport analysis of single-cell gene expression identifies developmental trajectories in reprogramming (Waddington-OT). *Cell* 176, 928–943.e22 (2019). DOI: 10.1016/j.cell.2019.01.006.
- [A9] Weiler P, *et al.* CellRank 2: unified fate mapping in multiview single-cell data. *Nat Methods* 21, 1196–1205 (2024). DOI: 10.1038/s41592-024-02303-9.
- [A10] Lange M, *et al.* CellRank for directed single-cell fate mapping. *Nat Methods* 19, 159–170 (2022). DOI: 10.1038/s41592-021-01346-6.

**Foundation models / benchmark / tooling**
- [A3] Cui H, *et al.* scGPT: toward building a foundation model for single-cell multi-omics. *Nat Methods* 21, 1470–1480 (2024). DOI: 10.1038/s41592-024-02201-0.
- [A4] Osakwe A, Huang EH, Li Y. scTimeBench: a streamlined benchmarking platform for single-cell time-series analysis. *bioRxiv* (2026). DOI: 10.64898/2026.03.16.712069. — **preprint**; cite only when stating scTimeBench-alignment.
- [A14] Saelens W, Cannoodt R, Todorov H, Saeys Y. A comparison of single-cell trajectory inference methods (dynverse). *Nat Biotechnol* 37, 547–554 (2019). DOI: 10.1038/s41587-019-0071-9.
- [A15] Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol* 19, 15 (2018). DOI: 10.1186/s13059-017-1382-0.
- [A16] Virshup I, *et al.* anndata: access and store annotated data matrices. *JOSS* 9(101), 4371 (2024). DOI: 10.21105/joss.04371.
- [A17] Domínguez Conde C, *et al.* Cross-tissue immune cell analysis reveals tissue-specific features in humans (CellTypist). *Science* 376, eabl5197 (2022). DOI: 10.1126/science.abl5197.
- [A18] Moon KR, *et al.* Visualizing structure and transitions in high-dimensional biological data (PHATE). *Nat Biotechnol* 37, 1482–1492 (2019). DOI: 10.1038/s41587-019-0336-3.

---

## Part B — Newly added & verified references

The redesigned project adds a **reference-mapping annotation branch** (scANVI/
scArches) and a **classical-trajectory comparison layer** (DPT/PAGA/Slingshot)
on top of the existing generative-method benchmark. The references below fill
those gaps. Each entry lists: title, authors, year, journal/server, DOI/URL,
main method/focus, why relevant, and what part of the project it supports.

### B1 — scVI (foundation for the annotation branch)
- **Title:** Deep generative modeling for single-cell transcriptomics
- **Authors:** Lopez R, Regier J, Cole MB, Jordan MI, Yosef N
- **Year / venue:** 2018, *Nature Methods* 15(12):1053–1058
- **DOI / URL:** 10.1038/s41592-018-0229-2 · https://www.nature.com/articles/s41592-018-0229-2
- **Main method / focus:** Variational-autoencoder generative model of scRNA-seq counts; batch-corrected latent space, the basis for scANVI/scArches.
- **Why relevant:** The annotation branch trains a latent model on the GSE242424 reference; scVI is the unsupervised backbone whose latent space is later made label-aware by scANVI.
- **Supports:** annotation branch, batch integration, latent-space (UMAP) figure.

### B2 — scANVI (the reference-annotation method)
- **Title:** Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models
- **Authors:** Xu C, Lopez R, Mehlman E, Regier J, Jordan MI, Yosef N
- **Year / venue:** 2021, *Molecular Systems Biology* 17(1):e9620
- **DOI / URL:** 10.15252/msb.20209620 · https://www.embopress.org/doi/full/10.15252/msb.20209620
- **Main method / focus:** Semi-supervised extension of scVI (scANVI) that propagates a partially observed cell-label set through the latent space to annotate unlabelled cells, with calibrated uncertainty.
- **Why relevant:** This is the exact method the requested branch uses to transfer GSE242424 author labels onto GSE178325/GSE230659. Currently missing from the project's references.
- **Supports:** annotation branch (core), benchmarking of annotation strategies, cell-annotation comparison figure.

### B3 — scArches (query-to-reference surgery)
- **Title:** Mapping single-cell data to reference atlases by transfer learning
- **Authors:** Lotfollahi M, Naghipourfar M, Luecken MD, *et al.* (Theis FJ, senior)
- **Year / venue:** 2022, *Nature Biotechnology* 40(1):121–130
- **DOI / URL:** 10.1038/s41587-021-01001-7 · https://www.nature.com/articles/s41587-021-01001-7
- **Main method / focus:** "Architectural surgery" — freeze a trained reference (e.g. scANVI) and fine-tune lightweight query-specific weights to map new data without retraining or sharing raw reference data.
- **Why relevant:** This is the *operational pattern* of the requested branch (the scArches scANVI surgery pipeline named in the task). It defines how GSE178325/GSE230659 are mapped onto the frozen GSE242424 reference.
- **Supports:** annotation branch (implementation pattern), reproducibility/portability of the reference model to the server.

### B4 — DPT (pseudotime the benchmark already uses)
- **Title:** Diffusion pseudotime robustly reconstructs lineage branching
- **Authors:** Haghverdi L, Büttner M, Wolf FA, Buettner F, Theis FJ
- **Year / venue:** 2016, *Nature Methods* 13:845–848
- **DOI / URL:** 10.1038/nmeth.3971 · https://www.nature.com/articles/nmeth.3971
- **Main method / focus:** Diffusion-map-based pseudotime ordering with branching detection.
- **Why relevant:** The repo's logs show DPT pseudotime is already computed (e.g. GSE230659); DPT is the natural pseudotime axis on which to test whether annotation strategy changes trajectory structure (scenarios D–F).
- **Supports:** trajectory inference strategy, pseudotime comparison figure.

### B5 — PAGA (topology-level trajectory comparison)
- **Title:** PAGA: graph abstraction reconciles clustering with trajectory inference through a topology-preserving map of single cells
- **Authors:** Wolf FA, Hamey FK, Plass M, *et al.* (Theis FJ, senior)
- **Year / venue:** 2019, *Genome Biology* 20:59
- **DOI / URL:** 10.1186/s13059-019-1663-x · https://link.springer.com/article/10.1186/s13059-019-1663-x
- **Main method / focus:** Partition-based graph abstraction; a coarse-grained connectivity graph over cell groups that preserves global topology.
- **Why relevant:** PAGA gives a state-level transition graph directly comparable to the project's reference lineage graphs, and lets us quantify how the *inferred topology* shifts between marker-silver and scANVI annotations.
- **Supports:** trajectory/lineage comparison figure, lineage-fidelity interpretation.

### B6 — Slingshot (independent lineage/pseudotime baseline)
- **Title:** Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics
- **Authors:** Street K, Risso D, Fletcher RB, Das D, Ngai J, Yosef N, Purdom E, Dudoit S
- **Year / venue:** 2018, *BMC Genomics* 19:477
- **DOI / URL:** 10.1186/s12864-018-4772-0 · https://pmc.ncbi.nlm.nih.gov/articles/PMC6007078/
- **Main method / focus:** Cluster-based minimum-spanning-tree lineages with principal-curve pseudotime; robust branching trajectory inference.
- **Why relevant:** Provides a second, annotation-dependent classical trajectory method so the robustness claim does not hinge on DPT alone.
- **Supports:** trajectory inference strategy, robustness analysis (optional figure).

### B7 — Trajectory-inference review (framing / benchmarking context)
- **Title:** Recent advances in trajectory inference from single-cell omics data
- **Authors:** Deconinck L, Cannoodt R, Saelens W, Deplancke B, Saeys Y
- **Year / venue:** 2021, *Current Opinion in Systems Biology* 27:100344
- **DOI / URL:** 10.1016/j.coisb.2021.05.005 · https://www.sciencedirect.com/science/article/pii/S2452310021000299
- **Main method / focus:** Review of modern TI methods, uncertainty reporting, and use of complementary information (time, unspliced RNA).
- **Why relevant:** Up-to-date framing for the introduction and for justifying method selection and uncertainty-aware evaluation.
- **Supports:** background/rationale, benchmarking strategy.

### B8 — Organoid time-course dataset+method (fallback-scope support)
- **Title:** Organoid single-cell genomic atlas uncovers human-specific features of brain development
- **Authors:** Kanton S, Boyle MJ, He Z, *et al.* (Treutlein B / Camp JG, senior)
- **Year / venue:** 2019, *Nature* 574:418–422
- **DOI / URL:** 10.1038/s41586-019-1654-9 · https://pubmed.ncbi.nlm.nih.gov/31619793/
- **Main method / focus:** Time-resolved scRNA-seq atlas of cerebral organoid development with pseudotemporal trajectories.
- **Why relevant:** Concrete, real human-organoid time-series option for the dataset-expansion fallback (priority 3 in the redesign) if additional iPSC-reprogramming series prove insufficient. **Data accession (E-MTAB-7552 / associated repositories) to be confirmed in `datasets/dataset_inventory.md` before any use.**
- **Supports:** dataset expansion (fallback tier), external-validation generalization.

---

## Part C — Coverage check vs. redesign needs

| Project need | Covered by |
|---|---|
| Generative time-series forecasting methods | A5 (scNODE), A6 (MIOFlow), A7 (PRESCIENT) |
| OT / fate-mapping lineage methods | A8 (WOT), A9/A10 (CellRank) |
| Classical pseudotime / topology methods | **B4 (DPT), B5 (PAGA), B6 (Slingshot)** |
| Reference-mapping annotation (the branch) | **B1 (scVI), B2 (scANVI), B3 (scArches)** |
| Marker/classifier annotation baselines | A17 (CellTypist) + in-repo marker-silver pipeline |
| Benchmarking philosophy / TI context | A4 (scTimeBench), A14 (dynverse), **B7 (review)** |
| Datasets (required) | A1 (GSE178325), A2 (GSE230659) + GSE242424 (see inventory) |
| Datasets (expansion / fallback) | A11–A13 (reprogramming), **B8 (organoid)**, + inventory |
| Tooling | A15 (Scanpy), A16 (AnnData), A3 (scGPT) |

**Newly added & verified references: 8 (B1–B8) — exceeds the required 5.**
No reference in Parts A or B is fabricated. The only item with an unverified
*data accession* (not the paper) is B8 (Kanton 2019); the publication itself is
verified, and its accession is flagged for confirmation in the dataset step.
