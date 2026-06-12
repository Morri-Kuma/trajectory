# Silver Standard Feasibility Rationale for the Trajectory Benchmark

说明：项目路径中使用了 `sliver_standard` 目录名，但本文统一采用领域中更常见的 `silver standard` / “银标准”表述。本文面向审稿风险控制，论证当前 `official_silver` / `marker_fm_transition_silver` 设计为何可行、如何避免过度声明，以及论文中应该如何表述。

## 结论摘要

当前项目中的 silver standard 不应被描述为“独立验证的真实谱系”或“gold-standard biology”。它更适合被定义为：一个由化学重编程文献、人工整理 marker、样本时间信息和轨迹约束共同冻结的、可复现的 benchmark reference state system。其用途是让不同 trajectory 方法在同一套状态标签和 reference graph 下进行公平比较，而不是证明某个模型发现了未经实验验证的真实谱系。

这个设计是可辩护的，原因有五点。第一，time-resolved scRNA-seq 是 destructive snapshot，单细胞跨时间的真实亲缘关系通常缺失，因此真实数据 benchmark 本来就常依赖 proxy reference 或 silver standard [R1, R2]。第二，trajectory inference benchmark 文献已经明确区分 real-data gold/silver standard，并用 reference trajectory/graph 对方法进行多指标比较 [R1]。第三，当前项目的 milestone order 不是任意聚类结果，而是来自人类化学重编程主文献与优化协议：GSE178325 使用含 XEN-like 阶段的四阶段路线，GSE230659 使用更直接的优化路线 [R11, R12]。第四，细胞类型/状态注释领域普遍承认不存在绝对 gold standard，因此需要不确定性、共识、拒判、classifier/reference-transfer 和 OOD 处理；本项目将 `ambiguous` 和 `unknown_or_ood` 排除在 official metrics 外，符合这一原则 [R3-R8, R10]。第五，项目已把早期 scGPT-only pseudostate provider 标为 legacy，当前 official silver 不再把单一 foundation-model embedding 当作主 reference，从而降低 circularity 风险；这也符合近期关于 scFM embedding 在 dynamics reconstruction 中表现不稳定的 benchmark 证据 [R9]。

## 当前设计定义

当前 active silver standard 是两个 frozen provider：

| Dataset | Provider | State key | Primary milestones | Reference graph |
| --- | --- | --- | --- | --- |
| GSE178325 | `gse178325_marker_fm_transition_silver_v1` | `final_milestone_label_coarse` | `hADSCs`, `epithelial_like`, `intermediate_plastic`, `xen_like`, `hCiPS` | `hADSCs -> epithelial_like -> intermediate_plastic -> xen_like -> hCiPS` |
| GSE230659 | `gse230659_marker_fm_transition_silver_v1` | `final_milestone_label_coarse` | `hADSCs`, `epithelial_like`, `intermediate_plastic`, `hCiPS` | `hADSCs -> epithelial_like -> intermediate_plastic -> hCiPS` |

项目证据：marker 和 graph 定义见 `benchmark/annotation/milestone_markers.yaml`；provider 注册和冻结状态见 `benchmark/ground_truth/registry.yaml`；每个 provider 的 `ground_truth_metadata.json` 记录了 `label_mode=official_silver`、`label_type=frozen_silver_standard`、`status=frozen_silver_standard_milestone_provider`、label counts、排除标签和 provenance。

GSE178325 provider 覆盖 80,475 个细胞，coarse labels 为 `epithelial_like` 67,817、`hADSCs` 5,706、`intermediate_plastic` 4,654、`hCiPS` 2,144、`xen_like` 154。GSE230659 provider 覆盖 75,194 个细胞，coarse labels 为 `epithelial_like` 67,112、`hCiPS` 5,370、`intermediate_plastic` 2,712；该数据集中没有分配 `hADSCs`，因为可用 earliest time point 是 0.5 而不是 day 0。项目证据来自两个 provider 的 `ground_truth_metadata.json`。

## 逐条论证

### 1. 在 time-series scRNA-seq 中，silver standard 是必要的，而不是偷换概念

Time-series scRNA-seq 通常在每个时间点破坏性测量不同细胞，不能追踪同一个细胞从早期到晚期的真实轨迹；因此跨时间 lineage relationship 在数据层面缺失 [R2]。Sha et al. 在 TIGON 工作中直接指出，time-series scRNA-seq 只提供 unpaired snapshots，细胞间跨时间 trajectory 缺失，需要通过 computational trajectory inference 或 optimal transport 一类方法重建 [R2]。Zhou et al. 对 scFM embedding 做 dynamics benchmark 时也将该问题表述为从 time-resolved single-cell transcriptomics 重建细胞轨迹，核心困难包括 embedding 和跨时间 transport map 估计 [R9]。

因此，审稿人若要求“每个细胞的真实未来命运”作为 benchmark 标准，是在要求当前实验设计通常不提供的 gold standard。可行且规范的做法是：把评估对象限定为“是否恢复与冻结 reference milestone graph 一致的状态转移结构”，并明确该 reference 是 silver standard，而非 lineage tracing gold standard [R1, R2]。

### 2. 真实数据 benchmark 使用 silver standard 有先例

Saelens et al. 对 45 个 trajectory inference 方法、110 个真实数据集和 229 个 synthetic 数据集进行比较，并明确把 real datasets 中无法独立于表达矩阵得到 reference trajectory 的情况归为 silver standard；他们仍然使用这些 reference 对 topology、branch assignment、cell position 和 feature dynamics 等方面进行评价 [R1]。这说明，真实 scRNA-seq trajectory benchmark 中使用 silver standard 并不是不合格设计，而是领域已有实践。

本项目的做法与该原则一致：不声称 `final_milestone_label_coarse` 是绝对 biological truth，而是把它作为统一 benchmark provider，所有方法在同一 provider、同一 reference graph、同一 label exclusion policy 下评估。项目证据：`benchmark/ground_truth/registry.yaml` 要求 provider 暴露 `state_labels.tsv`、`state_metadata.tsv`、`reference_graph.json` 和 `reference_graph_edges.csv`；runtime configs 中显式记录 `provider_id`、`state_key` 和 `reference_graph_path`。

### 3. Milestone graph 有明确生物学来源，不是 post hoc 任意拼接

GSE178325 的 `hADSCs -> epithelial_like -> intermediate_plastic -> xen_like -> hCiPS` 设计来自 Guan et al. 对人类体细胞化学重编程的轨迹分析：该研究通过构建 intermediate plastic state 实现人类 chemically induced pluripotent stem cells，并将早期 chemical-induced dedifferentiation 和 intermediate plastic state 作为关键过程 [R11]。Guan et al. 还把 XEN-like state 与化学重编程通向 pluripotency 的桥接过程联系起来 [R11]。

GSE230659 的 `hADSCs -> epithelial_like -> intermediate_plastic -> hCiPS` 设计来自 Liuyang et al. 的优化协议：该研究建立更高效、更快速的化学重编程系统，并报告优化协议使过程更直接，强调早期 proliferation 和 oxidative phosphorylation 活动 [R12]。因此，GSE230659 不强制加入 `xen_like` 作为 primary graph node 是可解释的：该数据集对应 optimized / shortcut-like protocol，而不是原始四阶段路线。

项目证据：`benchmark/annotation/milestone_markers.yaml` 将两篇主文献和 companion code 作为 marker provenance，并分别定义 GSE178325 和 GSE230659 的 primary milestones、marker sets、reference graph edges 和 excluded labels。

### 4. 使用 coarse milestone 而非细粒度 cell type，可降低 annotation 噪声和审稿风险

细胞注释文献普遍强调：cell type / cell state 的边界并非总有 gold standard，连续状态、技术噪声、参考集粒度和专家知识都会造成不确定性 [R6]。PopV 因此提出 ensemble consensus 与 uncertainty score，用于标出难注释细胞群，而不是假设单一自动注释器永远正确 [R6]。scANVI 也把 annotation 和 harmonization 置于 probabilistic framework 中，强调 biological/measurement noise 和已有 label uncertainty [R4]。

本项目选择 coarse milestones，是对这种不确定性的保守处理。`epithelial_like`、`intermediate_plastic`、`xen_like`、`hCiPS` 等节点是面向 reprogramming trajectory 的粗粒度状态，而不是试图区分所有亚型。coarse graph 更适合 benchmark lineage fidelity：它检验方法是否恢复大方向上的 milestone transition，而不是把低置信、连续过渡细胞强行离散成精细 cell types [R1, R6]。

项目证据：Stage 2 生成 `final_milestone_label_coarse` 作为 official metric key，同时保留 `final_milestone_label_expanded` 作为 supplementary；`build_trajectory_aware_labels.py` 明确说明 coarse label drives official scTimeBench-style metrics，expanded label 只作补充。

### 5. 不确定和 OOD 细胞没有被强行纳入 official metrics

自动注释领域的一个核心风险是把所有细胞强制分配到参考标签。scDetect 引入 rejection / unknown 选项，避免把低置信细胞错误归入已知类型 [R5]；TOSICA 也支持 unknown cell type discovery / unknown labeling，并在 reference-query 不完全重合时识别 query 中未覆盖的状态 [R7]。PopV 进一步强调不确定性和低共识区域应进入 manual scrutiny，而不是被当作高置信标签 [R6]。

本项目与这些原则一致：Stage 1 将低置信 marker cells 路由到 `ambiguous`；Stage 2 允许 `unknown_or_ood` 和 `ambiguous` 作为 valid final labels；official provider 明确把 `ambiguous` 和 `unknown_or_ood` 排除在 official metrics 外。这样做的意义是，benchmark 只在 reference graph 定义清楚的细胞状态上比较方法，避免用最不可靠的细胞惩罚或奖励模型 [R5-R7]。

项目证据：`build_marker_seed_labels.py`、`build_trajectory_aware_labels.py` 和 `build_milestone_providers.py` 均记录了 `ambiguous` / `unknown_or_ood` policy；两个 official provider 的 metadata 都含有 `excluded_from_official_metrics: ["ambiguous", "unknown_or_ood"]`。

### 6. 多证据/共识原则有统计学和单细胞注释文献基础

Dawid and Skene 的经典 latent-class/EM 框架表明，在真实标签不可直接观测时，可以利用多个 observer 的观测来估计 error rates 和 latent true response，并可形成 majority 或 weighted consensus [R3]。这为“没有 gold truth 时使用多来源证据构造可审计标签”提供了统计学基础。

单细胞领域也采用类似思想。PopV 用多个 annotation algorithms 的 ontology-aware voting 形成 consensus label 和 consensus score [R6]；scDetect 使用 rank-based ensemble learning 与 majority-vote probability prediction [R5]；scANVI 将已有 annotations 纳入半监督 probabilistic label transfer [R4]；CellTypist 在跨组织免疫图谱中以 curated reference 和机器学习分类器进行系统性注释，支持 classifier/reference-transfer 作为自动注释证据的一类可行路线 [R10]。因此，项目采用 marker、sample time、trajectory constraints，并保留 source/confidence/provenance 字段，是与 annotation literature 一致的工程化实现 [R3-R6, R10]。

项目证据：official provider 的 `annotation_votes.tsv` 保留 Stage 1、Stage 2、final label、confidence、source、membership 等列；`ground_truth_metadata.json` 记录 `confidence_key`、`source_key`、`membership_key`、`generated_by` 和 provenance。

### 7. 当前设计降低了 scGPT / foundation-model circularity

早期项目使用 scGPT-derived pseudostate provider 作为工作银标准，但该路线有 representational circularity 风险：如果 reference graph 由同一类 embedding 产生，那么使用相近 representation 的方法可能被不公平奖励。项目已经在 `registry.yaml` 中将 `scgpt_v1` 和 GSE178325 scGPT providers 标记为 `deprecated_legacy_pseudostate_provider`、`legacy_reference_only`，并把 primary result 切换到 marker-fm transition silver provider。

这个修订有文献支持。Zhou et al. 系统 benchmark 了 zero-shot scFM embeddings 在 cellular dynamics reconstruction 中的表现，发现 scFM embeddings 不一定优于 HVG-PCA baseline，并可能压缩微弱 temporal signal [R9]。这说明，把 scGPT/scFM embedding 单独作为 biological ground truth 会过度乐观；更稳妥的做法是让 foundation-model evidence 只作为辅助或 sensitivity evidence，而不是主 reference。

项目证据：`benchmark/reports/validate_legacy_scgpt_policy.py` 要求 scGPT provider 不得成为 primary；`benchmark/ground_truth/registry.yaml` 中 official silver providers 的 `annotation_method` 是 `marker_fm_transition_silver_provider`，而非 scGPT-only pseudostate provider。

### 8. Evaluation 与 label construction 被解耦，避免 self-comparison artifact

Embedding coherence 若直接拿 projected milestone labels 当 cluster labels，会产生自我比较假象。项目已经修正该点：official silver embedding coherence 从 `projected_embedding.npy` 重新构建 kNN graph 并运行 Leiden clustering，然后将 unsupervised clusters 与 `final_milestone_label_coarse` 比较。这样，ARI 衡量的是 projected embedding 是否自然分出与 silver milestones 一致的结构，而不是标签文件与自己比较。

这一点符合 trajectory benchmark 的基本逻辑：Saelens et al. 用多个维度评估 predicted trajectory 与 reference 的一致性，而不是让方法直接复用 reference label [R1]；scTimeBench 也将 forecast accuracy、embedding coherence 和 lineage fidelity 分成不同任务 [R13]。

项目证据：`docs/framework/experimental_framework_v2.md` 记录了 2026-05-11 的修正；`benchmark/reports/official_silver/official_silver_model_rankings.md` 说明 embedding coherence 使用 `leiden_from_projected_embedding:projected_embedding.npy:n_neighbors=15:resolution=1.0`。

### 9. Silver standard 适合用于公平比较，但必须伴随限制说明

Hippen et al. 显示，实验因素会影响 computational deconvolution / reference-based inference 的表现，包括 dissociation、protocol 和参考 profile 构建方式 [R8]。这类结果提醒我们：任何基于 single-cell reference 的 benchmark 都可能携带实验设计、样本组成和 protocol 的偏差。

因此，本项目应该在 manuscript 中明确限制：official silver 衡量的是方法对 frozen milestone reference graph 的恢复能力，不是独立实验验证的真实重编程谱系；GSE230659 存在 time/library confounding 风险；GSE178325 的 `xen_like` label count 很小，应避免对该边单独过度解释；不同 provider 或 bootstrap sensitivity 应作为补充分析，而不是装饰性附录 [R1, R6, R8]。

## 审稿人可能 challenge 与建议回应

| Challenge | 建议回应 | 依据 |
| --- | --- | --- |
| “这不是 gold standard，为什么能评估方法？” | 同意它不是 gold standard；本文只声称它是 frozen silver-standard reference。真实 time-series scRNA-seq 缺少逐细胞跨时间 lineage，领域 benchmark 已经使用 real-data silver standard。 | [R1, R2] |
| “标签是否只是 marker heuristic？” | 标签不是无约束 heuristic，而是文献定义 milestone order + curated marker sets + sample/time root + trajectory-aware resolution，并冻结为 provider。coarse label 用于 official metrics，expanded label 仅补充。 | [R11, R12], 项目 provider metadata |
| “为什么 GSE230659 没有 XEN-like？” | GSE230659 对应 optimized direct/shortcut-like chemical reprogramming route，文献报告该协议更直接；因此不强制 XEN-like 进入 primary graph。 | [R12], `milestone_markers.yaml` |
| “自动注释没有 ground truth，会传播错误。” | 正因如此，设计采用 coarse milestones、confidence/source tracking、`ambiguous` 和 `unknown_or_ood` 排除策略；这符合自动注释领域关于 uncertainty 和 unknown/rejection 的建议。 | [R4-R7] |
| “使用 foundation model 会不会循环论证？” | 早期 scGPT-only provider 已降级为 legacy；official silver 不以 scGPT pseudostates 作为主 reference。近期 benchmark 也显示 scFM embeddings 对 dynamics 并不稳定，不能直接当 truth。 | [R9], `registry.yaml` |
| “模型是否在和自己生成的 label 比较？” | Embedding coherence 已改为从 projected embedding 重新 Leiden clustering，再与 frozen milestone labels 比较；lineage fidelity 使用统一 provider graph。 | [R1, R13], official_silver report |

## 论文建议表述

建议使用：

> We therefore evaluate annotation-dependent metrics against a frozen silver-standard milestone reference, not against independently lineage-traced biological ground truth. The reference is constructed from literature-supported chemical reprogramming milestones, curated marker programs, sample-time anchoring for the starting state where available, and trajectory-aware resolution of ambiguous cells. Ambiguous and out-of-trajectory cells are excluded from official metrics. All method rankings are interpreted as performance against this frozen reference graph.

建议避免：

> Our method recovers the true biological lineage.

> The silver labels are ground truth cell types.

> scGPT-defined pseudostates provide the primary biological reference.

## 最小补充验证建议

1. 报告 provider identity：所有 summary rows 必须带 `provider_id`、`label_mode`、`state_key`、`reference_graph_path`。依据：benchmark provider layer 与 scTimeBench-style task separation [R13]。
2. 报告 sensitivity：至少保留 legacy / alternative provider 结果为 supplementary，不与 official silver primary rows 混排。依据：trajectory benchmark 中方法表现依赖 dataset/topology，没有 one-size-fits-all [R1]。
3. 报告 uncertainty：对 `ambiguous`、`unknown_or_ood`、低 confidence 或低 consensus 区域进行计数，并明确 official metrics 的 exclusion policy。依据：popV、scANVI、TOSICA、scDetect 关于 uncertainty / unknown 的设计 [R4-R7]。
4. 报告 limitations：特别说明 single-library-per-timepoint、GSE178325 `xen_like` 细胞数较少、GSE230659 day 0 缺失导致 `hADSCs` 不进入 observed labels。依据：reference-based computational analysis 对 experimental factors 敏感 [R8]。

## 文献与项目证据

### Local reference PDFs

[R1] Saelens W, Cannoodt R, Todorov H, Saeys Y. A comparison of single-cell trajectory inference methods. Nature Biotechnology. 2019;37:547-554. DOI: 10.1038/s41587-019-0071-9. Local file: `reference/sliver_standard/Saelens 等 - 2019 - A comparison of single-cell trajectory inference methods.pdf`.

[R2] Sha Y, Qiu Y, Zhou P, Nie Q. Reconstructing growth and dynamic trajectories from single-cell transcriptomics data. Nature Machine Intelligence. 2024;6:25-39. DOI: 10.1038/s42256-023-00763-w. Local file: `reference/sliver_standard/Sha 等 - 2023 - Reconstructing growth and dynamic trajectories from single-cell transcriptomics data.pdf`.

[R3] Dawid A P, Skene A M. Maximum likelihood estimation of observer error-rates using the EM algorithm. Journal of the Royal Statistical Society: Series C. 1979;28:20-28. Local file: `reference/sliver_standard/EM.pdf`.

[R4] Xu C, Lopez R, Mehlman E, Regier J, Jordan M I, Yosef N. Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models. Molecular Systems Biology. 2021;17:e9620. DOI: 10.15252/msb.20209620. Local file: `reference/sliver_standard/msb.20209620.pdf`.

[R5] Shen Y, Chu Q, Timko M P, Fan L. scDetect: a rank-based ensemble learning algorithm for cell type identification of single-cell RNA sequencing in cancer. Bioinformatics. 2021;37:4115-4122. Local file: `reference/sliver_standard/btab410.pdf`.

[R6] Ergen C, Xing G, Xu C, et al. Consensus prediction of cell type labels in single-cell data with popV. Nature Genetics. 2024;56:2731-2738. DOI: 10.1038/s41588-024-01993-3. Local file: `reference/sliver_standard/s41588-024-01993-3.pdf`.

[R7] Chen J, et al. Transformer for one stop interpretable cell type annotation. Nature Communications. 2023;14:223. DOI: 10.1038/s41467-023-35923-4. Local file: `reference/sliver_standard/s41467-023-35923-4.pdf`.

[R8] Hippen A A, Omran D K, Weber L M, et al. Performance of computational algorithms to deconvolve heterogeneous bulk ovarian tumor tissue depends on experimental factors. Genome Biology. 2023;24:239. DOI: 10.1186/s13059-023-03077-7. Local file: `reference/sliver_standard/13059_2023_Article_3077.pdf`.

[R9] Zhou X, Wang Z, Ling Y, et al. Benchmarking zero-shot single-cell foundation model embeddings for cellular dynamics reconstruction. bioRxiv. 2026. DOI: 10.64898/2026.03.10.710748. Local file: `reference/metrics/Zhou 等 - Benchmarking zero-shot single-cell foundation model embeddings for cellular dynamics reconstruction.pdf`.

[R10] Dominguez Conde C, Xu C, Jarvis L B, et al. Cross-tissue immune cell analysis reveals tissue-specific features in humans. Science. 2022;376:eabl5197. DOI: 10.1126/science.abl5197. Local file: `reference/sliver_standard/science.abl5197.pdf`. This supports the feasibility of CellTypist-style classifier annotation as a sensitivity or auxiliary route.

### Additional / project-curated references

[R11] Guan J, Wang G, Wang J, et al. Chemical reprogramming of human somatic cells to pluripotent stem cells. Nature. 2022;605:325-331. DOI: 10.1038/s41586-022-04593-5.

[R12] Liuyang S, Wang G, Wang Y, et al. Highly efficient and rapid generation of human pluripotent stem cells by chemical reprogramming. Cell Stem Cell. 2023;30:450-459.e9. DOI: 10.1016/j.stem.2023.02.008.

[R13] Osakwe A, Huang E H, Li Y. scTimeBench: a streamlined benchmarking platform for single-cell time-series analysis. bioRxiv. 2026. DOI: 10.64898/2026.03.16.712069.

### Project evidence

- `benchmark/annotation/milestone_markers.yaml`: curated milestone order, marker sets, graph edges, excluded labels, and literature provenance.
- `benchmark/annotation/build_marker_seed_labels.py`: Stage 1 marker/sample-time conservative routing with `ambiguous`.
- `benchmark/annotation/build_trajectory_aware_labels.py`: Stage 2 trajectory-aware resolution, `unknown_or_ood`, final coarse/expanded labels, and rule prohibiting evaluated-model embeddings.
- `benchmark/annotation/build_milestone_providers.py`: provider export, `official_silver` mode, state labels, annotation votes, reference graphs, metadata.
- `benchmark/ground_truth/registry.yaml`: provider registry, legacy scGPT deprecation, official silver provider status.
- `benchmark/ground_truth/providers/gse178325_marker_fm_transition_silver_v1/ground_truth_metadata.json`: GSE178325 frozen provider counts and provenance.
- `benchmark/ground_truth/providers/gse230659_marker_fm_transition_silver_v1/ground_truth_metadata.json`: GSE230659 frozen provider counts and provenance.
- `benchmark/reports/official_silver/official_silver_model_rankings.md`: official silver reporting and corrected embedding coherence cluster source.
