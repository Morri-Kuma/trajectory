# manuscript_outline

## Title

### 暂定中文题目

面向人类 iPSC 重编程轨迹的能力感知单细胞时间序列方法评测框架

### 暂定英文题目

A capability-aware benchmark for single-cell trajectory forecasting and lineage recovery in human iPSC reprogramming

### 备选英文题目

1. A silver-standard benchmark of trajectory forecasting methods in human iPSC reprogramming
2. Capability-aware evaluation of single-cell trajectory methods across human iPSC reprogramming datasets
3. Benchmarking trajectory forecasting, embedding coherence, and lineage recovery in human reprogramming time series

### 推荐文章类型

这篇文章最适合写成 benchmark / resource / methods application 类型，而不是新算法论文。核心贡献不是提出一个全新的轨迹推断模型，而是构建一个面向人类 iPSC 重编程过程的、可复现的、按方法能力分层的评测框架，并系统比较已有方法的适用范围。

## Three-sentence story

### 三句话主线

人类 iPSC 重编程是一个动态、异质且可能存在分支命运选择的细胞状态转换过程，单细胞时间序列数据为研究这一过程提供了重要窗口，但不同轨迹方法的输出能力并不相同。  

本项目构建了一个 scTimeBench-aligned 的 capability-aware benchmark，将 MIOFlow、PRESCIENT、scNODE、WOT 和 CellRank2 按其是否能生成 unseen-timepoint projected cells 分层，并分别评估 forecast accuracy、embedding coherence 和 lineage fidelity。  

在 GSE178325、GSE230659 和 GSE242424 的 official silver 结果中，不同方法在不同数据集、不同 scenario、不同指标上的优势并不一致，说明轨迹方法不存在单一赢家，合理的评测必须同时考虑方法能力、数据集背景和任务定义。

### 一句话摘要

我们建立了一个面向人类 iPSC 重编程单细胞时间序列的 official-silver benchmark，用统一的 scenario、state provider、reference graph 和 capability gating 评估现有轨迹方法在预测未来细胞状态、保持投影结构和恢复谱系关系方面的表现。

### 文章要回答的核心问题

1. 现有 trajectory forecasting / lineage inference 方法在 human iPSC reprogramming 时间序列中表现如何？
2. 生成式方法和 lineage-only 方法能否放在同一个评测框架下公平比较？
3. 哪些方法更擅长 projected-cell embedding coherence，哪些方法更擅长 state-level lineage recovery？
4. 同一方法在 GSE178325、GSE230659 和 GSE242424 上是否稳定？
5. extrapolation scenario 是否比 observed-time 或 mixed scenario 更难？

### 文章的主要贡献

1. 整理并冻结了面向 iPSC 重编程轨迹评测的 official silver provider，包括 `gse178325_marker_fm_transition_silver_v1`、`gse230659_marker_fm_transition_silver_v1` 和 `gse242424_oskm_reprogramming_silver_v1`。
2. 建立 observed-time scenarios A/B/C：A 为 observed-time / interpolation-like evaluation，B 为 future-time extrapolation，C 为 interpolation + extrapolation 混合评测。
3. 建立 capability-aware evaluation：MIOFlow、PRESCIENT、scNODE 可评估 forecast accuracy、embedding coherence、lineage fidelity；WOT 和 CellRank2 只评估 lineage fidelity。
4. 统一输出 per-method, per-scenario 的 projected expression、projected embedding、state transition matrix、lineage graph 和 metric summaries。
5. 发现 benchmark 结果具有明显 task dependence：PRESCIENT 在 GSE178325/GSE230659 combined official silver ranking 中 overall combined rank 第一，scNODE 在 embedding coherence composite 中第一，MIOFlow 在 lineage fidelity composite 中第一；GSE242424 OSKM silver 中 MIOFlow overall 第一。

## Main figures

### Figure 1. Benchmark framework overview

**目标：** 让读者一眼看懂整个 benchmark 是怎么工作的。

**建议内容：**

- 输入数据：GSE178325、GSE230659、GSE242424。
- 预处理：raw count matrix、QC、HVG2000、benchmark-ready h5ad。
- official silver provider：marker-FM transition silver provider 和 OSKM author-cluster-matched silver provider。
- 三个 observed-time scenario：A、B、C。
- 方法：MIOFlow、PRESCIENT、scNODE、WOT、CellRank2。
- capability gating：projection-capable methods 与 lineage-only methods 分开评估。
- 输出指标：forecast accuracy、embedding coherence、lineage fidelity。

**图注要点：** 强调这是一个 capability-aware benchmark，不强行让 WOT/CellRank2 做 projected-cell forecast。

### Figure 2. Dataset curation and frozen silver state systems

**目标：** 展示数据集来源、细胞数量、时间点和 reference state system。

**建议内容：**

- GSE178325：human reprogramming validation / benchmark dataset，使用 marker-FM transition silver provider。
- GSE230659：human iPSC reprogramming benchmark dataset，使用 marker-FM transition silver provider。
- GSE242424：OSKM fibroblast reprogramming，使用 author-cluster-matched subset，59,187 cells，HVG2000 features，9 个 timepoints `{0, 2, 4, 6, 8, 10, 12, 14, 16}`，10 个 coarse milestones，9 条 directed reference edges。
- 每个 dataset 对应的 provider、state key、reference graph。

**可以放的 panel：**

- Panel A：数据集表格。
- Panel B：每个数据集的 timepoint cell count barplot。
- Panel C：reference graph schematic。
- Panel D：benchmark scenario split 示意图。

### Figure 3. Official silver benchmark across GSE178325 and GSE230659

**目标：** 展示当前 official silver 主要结果。

**核心数字：**

- Combined summary 中，PRESCIENT overall combined rank 第一，combined rank score 1.7333。
- MIOFlow 第二，combined rank score 1.7389。
- scNODE 第三，combined rank score 1.7611。
- scNODE embedding coherence rank 第一，mean embedding composite rank score 1.8889。
- MIOFlow lineage fidelity rank 第一，mean lineage composite rank score 1.3667。

**建议图形：**

- 方法 × 指标 heatmap。
- 每个方法的 combined rank score barplot。
- embedding 和 lineage 分开画两个 panel，避免把不同任务混成一个单一结论。

**主结论：** 在 GSE178325/GSE230659 official silver benchmark 中，PRESCIENT 的综合排名略优，但三个 projection-capable 方法差距很小；不同维度的 winner 不同，说明不能只用单一指标选择 trajectory 方法。

### Figure 4. Dataset- and scenario-specific winners

**目标：** 展示不同 dataset/scenario 下方法表现不稳定，这是 benchmark 的重要发现。

**可写结果：**

- GSE178325 lineage fidelity 中，PRESCIENT 在 A/B/C 三个 scenario 中整体最强或并列最强。
- GSE230659 lineage fidelity 中，MIOFlow 在 A 和 C 最强，scNODE 在 B 最强。
- Embedding coherence winners 在不同 dataset/scenario 间切换：GSE178325 A 由 MIOFlow 胜出，GSE178325 B 由 scNODE 胜出，GSE230659 B/C 中 PRESCIENT 或 scNODE 更突出。

**建议图形：**

- dataset × scenario × method 的 rank tile plot。
- 每个 scenario 的 winner 标记。
- A/B/C scenario 的难度比较。

**主结论：** 方法表现依赖 dataset 和 scenario；future extrapolation 不应被 observed-time interpolation 结果替代。

### Figure 5. GSE242424 OSKM silver benchmark as an independent reprogramming setting

**目标：** 用 OSKM reprogramming 数据说明 benchmark 可以迁移到另一类 reprogramming 系统。

**核心数字：**

- 使用 author-cluster-matched GSE242424 subset：59,187 cells，HVG2000 features。
- MIOFlow overall rank 第一，combined rank score 1.6667。
- PRESCIENT 第二，combined rank score 2.1333。
- scNODE 第三，combined rank score 2.2000。
- MIOFlow mean AUROC 0.7798，mean Jaccard 0.3516。
- scNODE mean AUPRC 0.2566 最高。
- PRESCIENT mean ARI 0.1989 和 mean Wasserstein 0.0276 最好。

**建议图形：**

- GSE242424 method summary barplot。
- A/B/C scenario rank plot。
- lineage AUROC/AUPRC/Jaccard 与 embedding ARI/Wasserstein 分开显示。

**主结论：** OSKM silver benchmark 中 MIOFlow lineage recovery 最稳定，PRESCIENT 更擅长 preserving projected state structure，scNODE 在部分 scenario 和 AUPRC 上有优势，但 Scenario B extrapolation 明显更难。

### Supplementary figures

1. QC summary：每个数据集的 cells、detected genes、UMI counts、mitochondrial percentage。
2. Provider construction workflow：marker-FM transition silver provider 的构建流程。
3. Method capability table：每个方法支持哪些 evaluation dimensions。
4. Per-scenario forecast metrics：Wasserstein、Gaussian MMD、Energy MMD、Hausdorff。
5. Per-method lineage graph examples：predicted lineage graph vs reference graph。
6. Sensitivity analyses：如果后续补多 seed，可以放 run-to-run variance。

## Results outline

### Result 1. A capability-aware benchmark for human iPSC reprogramming trajectories

**本节要回答的问题：** 为什么需要一个新的 benchmark，而不是直接套用已有 trajectory benchmark？

**写作要点：**

人类 iPSC 重编程过程具有明显的时间性和状态转换特征，但公开单细胞数据通常缺少 per-cell fate tracing ground truth。现有方法的输出能力也不一致：有些方法能生成未来时间点的 projected cells，有些方法主要产生 transition / fate 信息。如果直接把所有方法放在同一个指标下排名，会对 lineage-only 方法不公平，也会对 generative 方法的 forecast 能力评估不足。

本项目因此采用 capability-aware evaluation。Projection-capable 方法，包括 MIOFlow、PRESCIENT 和 scNODE，被评估 forecast accuracy、embedding coherence 和 lineage fidelity。Lineage-only 方法，包括 WOT 和 CellRank2，只在 lineage fidelity 上评估。这个设计避免了人为构造不适合某些方法的输出，也使评测结论更贴近方法真实用途。

**本节可引用的项目文件：**

- `README.md`
- `benchmark/configs/method_capabilities.yaml`
- `benchmark/results/result_manifest.yaml`

**本节小结句：**

This design turns method heterogeneity from a nuisance into an explicit benchmark variable.

### Result 2. Frozen official silver providers define reproducible state systems

**本节要回答的问题：** 没有真实 lineage tracing ground truth 时，如何定义可复现的 benchmark target？

**写作要点：**

项目使用 frozen official silver providers 定义 benchmark state systems 和 reference graphs。当前 reportable providers 包括：

- `gse178325_marker_fm_transition_silver_v1`
- `gse230659_marker_fm_transition_silver_v1`
- `gse242424_oskm_reprogramming_silver_v1`

对于 GSE178325 和 GSE230659，文章应说明它们使用 marker-FM transition silver provider，并以 `final_milestone_label_coarse` 等 state key 作为 coarse milestone labels。对于 GSE242424，文章应说明 silver labels 来自 author ATAC-to-RNA cluster transfer annotations 和 matched local scRNA-seq cells，构成 OSKM reprogramming-specific provider。

这一节要避免把 silver provider 写成真正 biological ground truth。更准确的说法是：frozen silver-standard state system、reference state graph、benchmark reference graph。

**本节小结句：**

By freezing state providers and reference graphs, the benchmark makes lineage-fidelity evaluation reproducible while keeping the biological uncertainty explicit.

### Result 3. Official silver rankings reveal task-dependent performance across GSE178325 and GSE230659

**本节要回答的问题：** 在主要 official silver benchmark 中，哪个方法表现最好？

**核心结果：**

在 GSE178325 和 GSE230659 的 six projection-capable official silver runs 中，PRESCIENT 的 combined rank score 为 1.7333，整体排名第一；MIOFlow 为 1.7389，排名第二；scNODE 为 1.7611，排名第三。三者差距非常小，因此不应写成 PRESCIENT 绝对优越，而应写成 projection-capable 方法在综合排名上表现接近。

分指标看，scNODE 在 embedding coherence summary 中排名第一，mean embedding composite rank score 为 1.8889，mean weighted entropy 为 0.1740。MIOFlow 在 lineage fidelity summary 中排名第一，mean lineage composite rank score 为 1.3667，mean Jaccard top-k 为 0.1000。PRESCIENT 的 mean ARI 为 0.1150，在三者中最高，同时 overall combined rank 第一。

**建议写法：**

这些结果支持一个重要结论：trajectory benchmark 不应该只追求一个全局冠军。不同方法在 projected embedding、state transition recovery 和综合 rank 上各有优势，实际使用时应根据任务目标选择方法。

**本节小结句：**

The official silver benchmark shows no universal winner: PRESCIENT is strongest by combined rank, scNODE by embedding coherence, and MIOFlow by lineage fidelity.

### Result 4. Dataset and scenario splits expose method-specific robustness

**本节要回答的问题：** 方法表现是否稳定？A/B/C scenario 是否改变排名？

**核心结果：**

在 GSE178325 lineage fidelity 中，PRESCIENT 在 scenario A 和 B 排名第一，在 C 中也与 scNODE 接近或并列强势。GSE230659 中则不同：MIOFlow 在 scenario A 和 C 的 lineage composite rank score 为 1.0000，表现最好；scNODE 在 scenario B 中 AUROC 0.5641、AUPRC 0.2688、single-step recovery 0.3333，排名第一。

这说明同一方法在不同 dataset 和 scenario 下的表现会变化。Scenario B 作为 held-out future extrapolation 更接近实际预测任务，不能用 A 的 observed-time 表现替代。Scenario C 同时含 interpolation 和 extrapolation，更适合观察方法对混合任务的适应能力。

**建议写法：**

不要把 scenario A/B/C 只写成技术拆分。它们代表不同科学问题：

- A：方法能否在 observed temporal manifold 上保持一致？
- B：方法能否预测未来未见时间点？
- C：方法能否同时处理插值和外推？

**本节小结句：**

Scenario-specific rankings reveal that extrapolation and mixed-time prediction are not interchangeable with observed-time evaluation.

### Result 5. GSE242424 OSKM silver benchmark validates the framework in an independent reprogramming system

**本节要回答的问题：** 这个 benchmark 是否只适用于 GSE178325/GSE230659？

**核心结果：**

GSE242424 OSKM silver benchmark 使用 author-cluster-matched subset，共 59,187 cells，HVG2000 features，9 个 timepoints，10 个 coarse milestones 和 9 条 directed reference edges。A/B/C 三个 scenario 均完成 MIOFlow、PRESCIENT 和 scNODE formal runs。

在 GSE242424 中，MIOFlow overall rank 第一，combined rank score 1.6667；PRESCIENT 第二，combined rank score 2.1333；scNODE 第三，combined rank score 2.2000。MIOFlow 的 mean AUROC 为 0.7798，mean Jaccard 为 0.3516，说明它在 lineage recovery 上较稳定。PRESCIENT 的 mean ARI 为 0.1989，mean Wasserstein 为 0.0276，说明它在 projected state structure 和 distribution matching 上更有优势。scNODE 的 mean AUPRC 为 0.2566，是三者中最高，但在 Scenario B extrapolation 中表现较弱。

**本节小结句：**

The OSKM silver benchmark suggests that the same evaluation framework can be transferred to an independent reprogramming context, while preserving dataset-specific method rankings.

### Result 6. Capability gating prevents misleading comparisons between projection-capable and lineage-only methods

**本节要回答的问题：** WOT 和 CellRank2 应该如何纳入文章？

**写作要点：**

WOT 和 CellRank2 不应放入 forecast accuracy 或 embedding coherence 排名，因为它们不直接生成 unseen-timepoint projected cells。当前 official silver result 中，WOT 和 CellRank2 只有 GSE230659 scenario A lineage-only checks，因此可以作为 lineage-fidelity comparison / sanity check，而不是综合 ranking 的组成部分。

在 GSE230659 scenario A lineage-only 结果中，CellRank2 AUROC 0.4359、AUPRC 0.1935；WOT AUROC 0.4103、AUPRC 0.1842。它们低于 projection-capable 方法在同一 scenario 的 lineage composite rank，但这些结果更适合作为 lineage-only 方法在该 silver graph 上的参考，而不是用来判定其整体方法优劣。

**本节小结句：**

Capability gating makes the benchmark fairer by reporting each method only on tasks it is designed to support.

## Methods outline

### 1. Data sources and benchmark datasets

需要写清楚三个 reportable dataset groups：

- GSE178325 marker-FM transition silver HVG2000 benchmark。
- GSE230659 marker-FM transition silver HVG2000 benchmark。
- GSE242424 OSKM silver HVG2000 benchmark。

每个数据集需要交代：

- 生物学背景。
- 原始数据类型。
- 时间点或 stage 信息。
- 是否有 per-cell ground truth。
- 最终进入 benchmark 的细胞子集。
- 使用的 provider 和 reference graph。

### 2. Preprocessing and benchmark input construction

建议写：

- 从 raw count matrix 构建 AnnData / h5ad。
- 基础 QC：过滤低质量细胞、检查 detected genes、UMI counts、mitochondrial percentage。
- 标准化和 log transformation。
- HVG2000 feature selection。
- 保留 timepoint / scenario / dataset metadata。
- 生成 benchmark-ready input files。

如果不同数据集的处理流程不同，应分别列出，避免写得过于抽象。

### 3. Silver-standard state providers

需要写清楚：

- 为什么使用 silver-standard 而不是 true ground truth。
- provider 是 frozen/versioned 的。
- 每个 provider 包含 state labels、state metadata、reference graph edges 和 provider metadata。
- 旧版 scGPT-v1 provider 仅作为历史比较，不作为当前主结果。

推荐术语：

- frozen silver provider
- silver-standard milestone labels
- reference state graph
- benchmark reference graph

避免术语：

- true lineage ground truth
- experimentally verified fate map
- real lineage tree

### 4. Benchmark scenarios

当前正式结果聚焦 observed-time A/B/C：

- Scenario A：observed-time / interpolation-like setting。
- Scenario B：future-time extrapolation，hold out later timepoints。
- Scenario C：mixed interpolation and extrapolation。

GSE242424 需要具体写：

- A：all 9 timepoints used for training/evaluation。
- B：train on `{0, 2, 4, 6, 8, 10}` and hold out `{12, 14, 16}`。
- C：train on `{0, 2, 6, 10, 14}` and hold out `{4, 8, 12}` plus `{16}`。

### 5. Methods and adapters

Projection-capable methods：

- MIOFlow：评估 forecast accuracy、embedding coherence、lineage fidelity。
- PRESCIENT：评估 forecast accuracy、embedding coherence、lineage fidelity。
- scNODE：评估 forecast accuracy、embedding coherence、lineage fidelity。

Lineage-only methods：

- WOT：只评估 lineage fidelity。
- CellRank2：只评估 lineage fidelity。

需要说明每个方法通过 adapter 统一输出：

- projected expression 或 projected embedding。
- state transition matrix。
- lineage graph edges。
- run metadata。

### 6. Forecast accuracy metrics

用于 projection-capable methods。

建议指标：

- Wasserstein distance。
- Gaussian MMD。
- Energy distance MMD。
- Hausdorff loss。

解释方向：

- lower is better。
- 衡量 projected cells 与 held-out / observed cells 的表达分布距离。
- 不适用于 WOT 和 CellRank2。

### 7. Embedding coherence metrics

用于有 projected embedding 或可从 projected cells 构建 embedding 的方法。

当前 official silver ranking 使用：

- adjusted Rand index。
- mean normalized entropy。
- weighted mean normalized entropy。
- cluster source: `leiden_from_projected_embedding:projected_embedding.npy:n_neighbors=15:resolution=1.0`。

解释方向：

- ARI 越高越好。
- entropy 越低越好。
- 评估 projected embedding 是否保留 silver milestone label structure。

### 8. Lineage fidelity metrics

用于所有能产生 state-level transition output 的方法。

指标：

- AUROC。
- AUPRC。
- Jaccard top-k。
- single-step recovery。
- multi-step recovery。

解释方向：

- 衡量 predicted transition graph 与 frozen reference graph 的一致性。
- 对 WOT 和 CellRank2，这是当前唯一正式适用的指标族。

### 9. Ranking and summary rules

需要写：

- 在每个 dataset/scenario 内对方法排名。
- ARI、AUROC、AUPRC、Jaccard、recovery metrics high-to-low。
- entropy 和 distance metrics low-to-high。
- composite rank score 是 within-dataset/scenario ranks 的平均值。
- 对不支持某类输出的方法，不计算该类指标，也不强行放入 overall projection-capable summary。

### 10. Software and reproducibility

需要写：

- 代码组织在 `benchmark/`。
- reportable result sets 由 `benchmark/results/result_manifest.yaml` 定义。
- official silver summaries 在 `benchmark/reports/official_silver/` 和 `benchmark/reports/gse242424_oskm_silver/`。
- configs 位于 `benchmark/configs/` 和 `benchmark/configs/runtime/`。
- 每个 run 输出 `run_metadata.json`、metrics JSON/CSV、transition matrix 和 lineage graph files。

## Introduction outline

### Paragraph 1. Biological motivation

人类体细胞重编程为 iPSC 是理解细胞命运转换、可塑性和再生医学的重要模型。化学重编程和 OSKM reprogramming 均表现出复杂的时间依赖状态变化，单细胞转录组数据可以解析群体平均实验无法观察到的中间状态和分支结构。

### Paragraph 2. Computational problem

已有许多单细胞 trajectory inference 和 forecasting 方法，但它们的目标和输出形式不同。有些方法重点推断 fate transition 或 transport coupling，有些方法可以生成未来时间点的细胞状态。缺少统一且公平的评测框架时，很难判断方法在真实 reprogramming 时间序列中的适用范围。

### Paragraph 3. Benchmark gap

通用 trajectory benchmark 提供了重要参考，但人类 iPSC 重编程具有特定的数据结构、时间点设计和 ground-truth 缺失问题。因此需要一个 domain-specific benchmark，将 state provider、scenario split、method capability 和 evaluation metrics 一起冻结。

### Paragraph 4. Our work

本研究构建了一个 scTimeBench-aligned official silver benchmark，系统评估 MIOFlow、PRESCIENT、scNODE、WOT 和 CellRank2 在 human iPSC reprogramming datasets 上的表现。我们显示，不同方法在 forecast、embedding coherence 和 lineage fidelity 上表现不同，且 ranking 随 dataset 和 scenario 改变。

## Discussion outline

### Main message 1. No universal best method

综合结果说明没有单一 trajectory 方法在所有数据集、scenario 和指标上都最好。PRESCIENT 在 GSE178325/GSE230659 combined rank 中略优，scNODE 在 embedding coherence 中更强，MIOFlow 在 lineage fidelity 和 GSE242424 OSKM silver benchmark 中更稳定。

### Main message 2. Capability-aware evaluation is necessary

WOT 和 CellRank2 与 MIOFlow/PRESCIENT/scNODE 的能力不同。把所有方法强行放进 forecast ranking 会产生误导。capability gating 是本文重要方法学贡献之一。

### Main message 3. Silver-standard evaluation is useful but limited

由于缺少真实 lineage tracing，本文使用 frozen silver providers 和 reference graphs。这样可以支持可复现 benchmark，但不能把结果解释为对真实细胞命运树的完全验证。

### Main message 4. Future work

未来可以补充：

- 多 seed repeated runs，估计 run-to-run variance。
- 更多 human reprogramming datasets。
- 更严格的 external biological validation。
- pseudotime scenarios D/E/F。
- perturbation-aware 或 lineage tracing datasets。

## Key limitations

1. **Silver standard is not true lineage tracing.** 当前 state labels 和 reference graphs 是 frozen silver providers，不是真实谱系追踪结果。
2. **Run-to-run variance is not yet fully estimated.** 当前多个 formal reports 是每个 method/scenario 一个 formal run；投稿前最好补 3 个以上 seeds。
3. **Timepoint and batch may be confounded.** 尤其公开 reprogramming 数据中，timepoint、sample、library 和 condition 可能不完全正交。
4. **WOT and CellRank2 have limited current official silver coverage.** 当前 official silver 中它们主要有 GSE230659 scenario A lineage-only checks，不能与 projection-capable methods 做完整 combined comparison。
5. **Forecast metrics depend on representation and feature space.** HVG2000、PCA/inverse transform、projected embedding construction 都会影响 forecast 和 embedding metrics。
6. **GSE242424 uses a matched subset.** OSKM silver benchmark 使用 59,187-cell author-cluster-matched subset，不代表本地完整 156,969-cell dataset。
7. **Pseudotime scenarios are not current main results.** D/E/F 仍是 framework design，不应写成已完成主结果。

## Target journal type

### 比较合适的文章定位

1. **Benchmark/resource paper**
   最符合当前项目状态。重点是数据整理、评测框架、方法比较和可复现结果。

2. **Methods application paper**
   如果强调 capability-aware evaluation 和 official silver provider design，也可以投偏计算方法应用的期刊。

3. **Short report / application note**
   如果短期内没有时间补多 seed 和更多验证，可以先压缩成短文或 application note。

### 潜在目标期刊方向

- Bioinformatics
- GigaScience
- BMC Bioinformatics
- Briefings in Bioinformatics，要求会更高
- Genome Biology / Nature Communications，除非补充更多生物学发现和验证，否则难度较大
- single-cell / stem-cell domain journal 的 methods/resource section

### 当前投稿前最重要的补强

1. 冻结 manuscript 使用的 official result set，只引用 `official_silver` manifest 中的 reportable runs。
2. 补充 MIOFlow、PRESCIENT、scNODE 的 multi-seed stability analysis。
3. 统一所有文档里关于 legacy scGPT-v1 与 current official silver provider 的表述。
4. 补齐 Figure 1-5，尤其是 framework overview、dataset/provider summary 和 main ranking heatmap。
5. 明确代码和数据可复现路径，避免结果文件分散造成审稿人困惑。

## Candidate abstract draft

Single-cell time-series data provide a powerful view of human induced pluripotent stem cell reprogramming, but benchmarking trajectory methods in this setting remains challenging because methods differ in whether they infer transitions, generate future cells, or preserve projected cellular structure. Here, we present a capability-aware official-silver benchmark for human iPSC reprogramming trajectories. The benchmark freezes dataset-specific milestone providers, reference state graphs, observed-time scenario splits, and metric families, and evaluates each method only on tasks supported by its output type. Across GSE178325 and GSE230659 marker-FM transition silver benchmarks, PRESCIENT achieved the best combined rank among projection-capable methods, while scNODE ranked best for embedding coherence and MIOFlow ranked best for lineage fidelity. In an independent GSE242424 OSKM silver benchmark, MIOFlow achieved the best combined score, whereas PRESCIENT better preserved projected state structure and scNODE showed strong AUPRC in selected scenarios. These results show that trajectory method performance is task-, scenario-, and dataset-dependent, supporting capability-aware benchmarking as a practical strategy for evaluating single-cell reprogramming trajectories without overstating silver-standard reference graphs as true lineage ground truth.

## Writing checklist

- [ ] 确定主标题和文章类型。
- [ ] 确定 Figure 1-5 的最终内容。
- [ ] 只使用 `benchmark/results/result_manifest.yaml` 中的 current official_silver reportable runs。
- [ ] 把 legacy scGPT-v1 结果放入 supplement 或 historical comparison，不放主结果。
- [ ] 给每个 dataset/provider 写一句准确描述。
- [ ] 给每个 metric family 写清楚 eligible methods。
- [ ] 补 multi-seed stability 或在 limitations 中明确说明尚未完成。
- [ ] 写第一版 Results，不急着写 Introduction。
- [ ] 最后写 Abstract 和 Discussion。
