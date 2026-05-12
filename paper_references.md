# References for the trajectory manuscript

This list is organized for the current `trajectory` benchmark manuscript. The
first section contains the references that should be cited if the manuscript
describes the active benchmark design, datasets, annotation sources, and
methods. The second section contains optional references for datasets or
background material that are present in the workspace but are not part of the
primary formal benchmark.

## Core references

[1] Guan J, Wang G, Wang J, et al. Chemical reprogramming of human somatic cells to pluripotent stem cells[J]. Nature, 2022, 605: 325-331. DOI: [10.1038/s41586-022-04593-5](https://doi.org/10.1038/s41586-022-04593-5).

[2] Liuyang S, Wang G, Wang Y, et al. Highly efficient and rapid generation of human pluripotent stem cells by chemical reprogramming[J]. Cell Stem Cell, 2023, 30: 450-459.e9. DOI: [10.1016/j.stem.2023.02.008](https://doi.org/10.1016/j.stem.2023.02.008).

[3] Cui H, Wang C, Maan H, et al. scGPT: toward building a foundation model for single-cell multi-omics using generative AI[J]. Nature Methods, 2024, 21: 1470-1480. DOI: [10.1038/s41592-024-02201-0](https://doi.org/10.1038/s41592-024-02201-0).

[4] Osakwe A, Huang E H, Li Y. scTimeBench: a streamlined benchmarking platform for single-cell time-series analysis[EB/OL]. bioRxiv, 2026. DOI: [10.64898/2026.03.16.712069](https://doi.org/10.64898/2026.03.16.712069).  
Note: this is a preprint as of 2026-05-12; cite it only if the manuscript explicitly states that the benchmark is scTimeBench-aligned.

[5] Zhang J, Larschan E, Bigness J, Singh R. scNODE: generative model for temporal single cell transcriptomic data prediction[J]. Bioinformatics, 2024, 40(Supplement_2): ii146-ii154. DOI: [10.1093/bioinformatics/btae393](https://doi.org/10.1093/bioinformatics/btae393).

[6] Huguet G, Magruder D S, Tong A, et al. Manifold interpolating optimal-transport flows for trajectory inference[C]//Advances in Neural Information Processing Systems 35. 2022: 29705-29718. URL: [NeurIPS paper](https://proceedings.neurips.cc/paper_files/paper/2022/file/bfc03f077688d8885c0a9389d77616d0-Paper-Conference.pdf). arXiv DOI: [10.48550/arXiv.2206.14928](https://doi.org/10.48550/arXiv.2206.14928).

[7] Yeo G H T, Saksena S D, Gifford D K. Generative modeling of single-cell time series with PRESCIENT enables prediction of cell trajectories with interventions[J]. Nature Communications, 2021, 12: 3222. DOI: [10.1038/s41467-021-23518-w](https://doi.org/10.1038/s41467-021-23518-w).

[8] Schiebinger G, Shu J, Tabaka M, et al. Optimal-transport analysis of single-cell gene expression identifies developmental trajectories in reprogramming[J]. Cell, 2019, 176(4): 928-943.e22. DOI: [10.1016/j.cell.2019.01.006](https://doi.org/10.1016/j.cell.2019.01.006).

[9] Weiler P, Lange M, Klein M, Pe'er D, Theis F. CellRank 2: unified fate mapping in multiview single-cell data[J]. Nature Methods, 2024, 21: 1196-1205. DOI: [10.1038/s41592-024-02303-9](https://doi.org/10.1038/s41592-024-02303-9).

[10] Lange M, Bergen V, Klein M, et al. CellRank for directed single-cell fate mapping[J]. Nature Methods, 2022, 19: 159-170. DOI: [10.1038/s41592-021-01346-6](https://doi.org/10.1038/s41592-021-01346-6).  
Note: include this if the text discusses CellRank beyond the CellRank 2 implementation.

## Optional references for background or non-primary datasets

[11] Chen X, Lu Y, Wang L, et al. A fast chemical reprogramming system promotes cell identity transition through a diapause-like state[J]. Nature Cell Biology, 2023, 25(8): 1146-1156. DOI: [10.1038/s41556-023-01193-x](https://doi.org/10.1038/s41556-023-01193-x).  
Use if the manuscript mentions the local FCR mouse dataset or fast chemical reprogramming background.

[12] Peng F, Wang Y, Cheng L, et al. Chemical reprogramming of human blood cells to pluripotent stem cells[J]. Cell Stem Cell, 2025, 32(8): 1192-1199.e11. DOI: [10.1016/j.stem.2025.07.003](https://doi.org/10.1016/j.stem.2025.07.003).  
Use if the manuscript mentions GSE298212 / blood-cell hCiPS chemical reprogramming.

[13] Wang Y, Peng F, Yang Z, et al. A rapid chemical reprogramming system to generate human pluripotent stem cells[J]. Nature Chemical Biology, 2025, 21: 1030-1038. DOI: [10.1038/s41589-024-01799-8](https://doi.org/10.1038/s41589-024-01799-8).  
Use if the manuscript discusses newer rapid human hCiPS chemical reprogramming beyond GSE230659/GSE178325.

[14] Saelens W, Cannoodt R, Todorov H, Saeys Y. A comparison of single-cell trajectory inference methods[J]. Nature Biotechnology, 2019, 37: 547-554. DOI: [10.1038/s41587-019-0071-9](https://doi.org/10.1038/s41587-019-0071-9).  
Use if the introduction includes a general trajectory-inference benchmark context.

[15] Wolf F A, Angerer P, Theis F J. SCANPY: large-scale single-cell gene expression data analysis[J]. Genome Biology, 2018, 19: 15. DOI: [10.1186/s13059-017-1382-0](https://doi.org/10.1186/s13059-017-1382-0).  
Use if the methods section states that preprocessing, neighborhood graphs, Leiden clustering, UMAP, or DPT were run through Scanpy.

[16] Virshup I, Rybakov S, Theis F J, Angerer P, Wolf F A. anndata: access and store annotated data matrices[J]. Journal of Open Source Software, 2024, 9(101): 4371. DOI: [10.21105/joss.04371](https://doi.org/10.21105/joss.04371).  
Use if the manuscript explicitly mentions AnnData/h5ad as the benchmark data container.

[17] Dominguez Conde C, Xu C, Jarvis L B, et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans[J]. Science, 2022, 376(6594): eabl5197. DOI: [10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197).  
Use if the manuscript explicitly mentions CellTypist as a classifier-based annotation route.

[18] Moon K R, van Dijk D, Wang Z, et al. Visualizing structure and transitions in high-dimensional biological data[J]. Nature Biotechnology, 2019, 37: 1482-1492. DOI: [10.1038/s41587-019-0336-3](https://doi.org/10.1038/s41587-019-0336-3).  
Use if the manuscript discusses PHATE or MIOFlow's PHATE-distance / manifold-distance component in detail.
