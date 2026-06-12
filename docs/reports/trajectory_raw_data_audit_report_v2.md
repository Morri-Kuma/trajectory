# Trajectory data audit report v2

范围：`C:\Users\37620\trajectory\data`，明确排除 `data\processed`。

这版按数据集逐一回答同一套 10 个问题；没有把 `gse230659(human\rna_seq` 单独作为主分析对象。所有“无法确认”的项目都表示本地原始文件里缺少相应信息，而不是没有生物学意义。

## 全局结论

| 数据集 | 是否是单细胞表达矩阵 | 本地可直接用于轨迹输入 | 主要原因 |
|---|---|---|---|
| `gse230659(human` | 是，10x scRNA | 是，适合 | human chemical reprogramming to hCiPSC，时间/阶段清楚，raw counts 完整 |
| `gse178325_human` | 是，10x scRNA | 有条件适合 | human reprogramming 数据丰富，但样本来源/日期/阶段复杂，需要先整理 batch/donor/condition |
| `GSE298212(human` | 是，10x scRNA | 有条件适合 | human blood-cell chemical reprogramming，时间短，只有 5 个样本 |
| `FCR_iPSC(mouse` | 是，scRNA + 多组学 | 不适合作为本项目主 input | 小鼠 FCR 数据，物种和本项目 human benchmark 不一致 |
| `gse242424` | 不完整 scRNA | 暂不适合 | 缺 features/genes 文件，barcode 数异常大，无法确认基因名/QC |
| `gse247600` | 是，CSV single-cell/nucleus expression | 不适合作为本项目主 input | 神经分化/PD 模型，不是 hCiPSC reprogramming；样本细胞数很少 |
| `gse175634` | 不完整 MTX | 暂不适合 | 只有 count matrix，没有本地 cell/gene metadata |
| `gse280956(human` | 否 | 不适合 | ChIP-seq bigWig，不是表达矩阵 |
| `gse136314` | 无本地文件 | 不能评估 | 空目录 |

## 1. `gse230659(human`

| 问题 | 回答 |
|---|---|
| 1. 生物学背景 | human chemical reprogramming，目标是从人类体细胞诱导 hCiPSC。GEO series title 为 highly efficient and rapid generation of human pluripotent stem cells by chemical reprogramming。 |
| 2. 时间信息 | Stage I: D0.5/D2/D4/D8；Stage II: D4/D8；Stage III: D0.33/D0.67/D1/D2/D4/D6/D8/D12；终点为 hCiPSC。时间不是单一线性 day，而是 stage + day 的组合。 |
| 3. 原始格式和矩阵类型 | 15 个 scRNA 10x 文件夹，每个样本有 `matrix.mtx.gz`、`barcodes.tsv.gz`、`features.tsv.gz`。MatrixMarket header 为 integer，逐项检查为 raw count matrix。没有 `.h5ad/.h5/loom/Seurat object` 原始输入。features 第 1 列是 Ensembl ID，第 2 列是 gene symbol。样本内 barcode 唯一；跨样本 raw barcode 有 787 个重复，所以合并前必须加 sample 前缀。 |
| 4. 批次信息 | 每个样本就是一个 library/batch；本地 metadata 显示 scRNA 样本 source/cell line 都是 `0618`，platform GPL24676，10x Genomics，Cell Ranger v3.1.0，hg19。不同时间点与 batch/library 完全重合；没有独立 donor/replicate 字段。 |
| 5. 细胞数量和时间分布 | 本地 raw barcode 总数 79,907。最少样本是 StageI D4，仅 835 cells；StageIII D0.67 为 3,045；StageIII D6 为 3,623。 |
| 6. QC 指标 | 已全量逐细胞计算：每细胞 detected genes、UMI/counts、mt%、ribosomal%。低质量比例使用 `n_genes<200 OR counts<500 OR mt>20%`。doublet 信息本地没有，需要额外跑 Scrublet/DoubletFinder 等。 |
| 7. 原始细胞类型注释 | 没有 per-cell cell type annotation。只有 sample/stage/time 级别信息。 |
| 8. Ground truth / 外部参考 | 有实验定义 stage 和 hCiPSC endpoint；没有 per-cell fate label。marker genes、lineage model、生物学假设需从论文/外部知识整理，原始文件中未提供冻结版 lineage graph。 |
| 9. RNA velocity | 没有 `.loom`，没有 spliced/unspliced matrix。10x 3' scRNA 理论上可从 BAM/velocyto 重新生成 velocity，但本地文件不够。当前不适合直接 scVelo/CellRank velocity。 |
| 10. 是否适合作为本项目 input | 适合，是 human hCiPSC chemical reprogramming 时间序列 raw count。主要风险是 timepoint 与 batch 完全重合、缺少 replicate、缺少 per-cell ground truth。 |

细胞/QC 摘要：

| sample | cells | median genes | median counts | mt% median | ribo% median | low-quality% |
|---|---:|---:|---:|---:|---:|---:|
| StageI_Day0.5 | 6,331 | 5,451 | 27,275 | 3.36 | 20.62 | 0.57 |
| StageI_Day2 | 4,337 | 5,676 | 31,231 | 4.22 | 20.16 | 0.83 |
| StageI_Day4 | 835 | 6,080 | 47,313 | 1.09 | 27.08 | 7.43 |
| StageI_Day8 | 5,919 | 3,602 | 19,061 | 7.30 | 38.75 | 8.94 |
| StageII_Day4 | 6,059 | 4,222 | 17,086 | 5.44 | 27.06 | 1.01 |
| StageII_Day8 | 5,783 | 3,754 | 16,479 | 4.54 | 32.40 | 2.89 |
| StageIII_Day0.33 | 6,281 | 3,990 | 13,273 | 4.60 | 27.18 | 1.42 |
| StageIII_Day0.67 | 3,045 | 5,064 | 19,965 | 3.42 | 26.80 | 1.94 |
| StageIII_Day1 | 4,950 | 4,295 | 16,888 | 3.84 | 29.00 | 1.03 |
| StageIII_Day2 | 7,036 | 3,017 | 9,222 | 4.77 | 28.18 | 1.18 |
| StageIII_Day4 | 4,638 | 4,334 | 14,102 | 8.37 | 17.85 | 2.54 |
| StageIII_Day6 | 3,623 | 4,578 | 18,220 | 8.09 | 22.44 | 2.26 |
| StageIII_Day8 | 4,812 | 4,880 | 22,143 | 7.54 | 22.26 | 4.90 |
| StageIII_Day12 | 5,684 | 4,938 | 25,717 | 6.70 | 18.13 | 6.18 |
| hCiPSC | 10,574 | 3,362 | 10,610 | 6.52 | 22.40 | 13.39 |

## 2. `gse178325_human`

| 问题 | 回答 |
|---|---|
| 1. 生物学背景 | human chemical reprogramming / hCiPSC 相关数据，包含 HEF、hADSC、hASF、stage I-IV、hCiPSC、H1 等样本，适合作外部验证或扩展数据。 |
| 2. 时间信息 | 文件名/样本名编码多套时间与阶段：Stage I/II/III/IV，S1D0-D2、S1D16、S2D24、S2D4/D8/D12、S3D4/D8、S4D1/D2/D4/D10，以及 hCiPSC/H1/起始细胞样本。还有日期/批次后缀 0330、0605、0618、0809、1117、1230。 |
| 3. 原始格式和矩阵类型 | 本地有 `.tar.gz` 原始补充包，也有已解压的 10x `matrix.mtx.gz`、`barcodes.tsv.gz`、`features.tsv.gz`。矩阵为 integer raw counts。features 第 1 列 Ensembl ID，第 2 列 gene symbol。没有 `.h5ad/.h5/loom/Seurat object`。 |
| 4. 批次信息 | 每个 GSM/library 是 batch。时间点、stage、日期和起始细胞来源不完全正交，存在明显 batch/condition 混合；部分时间点只来自一个日期/来源。donor/replicate 需要从 GEO/论文进一步整理，本地文件名只能确认日期/来源样式。 |
| 5. 细胞数量和时间分布 | 去除 `rna_seq_10x` 重复副本后，唯一 10x 样本 32 个，总 221,916 cells。细胞数最低样本是 0330-stageII-JNKIN8，2,919 cells；另有 0618_S2D8 为 3,005 cells、0618_S3D8 为 3,605 cells。 |
| 6. QC 指标 | 已对 32 个唯一样本全量计算 detected genes、UMI/counts、mt%、ribosomal%、低质量比例。部分 hCiPSC/起始细胞样本低质量比例很高，例如 ADSC0618 59.46%、HEF1117 64.78%、ASF1230 46.10%，需要特别过滤/核查。doublet 本地无字段。 |
| 7. 原始细胞类型注释 | 没有 per-cell cell type annotation。样本名有细胞来源和阶段标签。 |
| 8. Ground truth / 外部参考 | 有实验 stage/time、起始细胞和 hCiPSC/H1 对照；没有 per-cell fate label 或本地 lineage graph。marker genes/lineage hypothesis 需外部整理。 |
| 9. RNA velocity | 无 `.loom`，无 spliced/unspliced matrix。只有 gene count MTX，不适合直接 RNA velocity。 |
| 10. 是否适合作为本项目 input | 有条件适合。它很有价值，但 batch/condition/time 结构复杂，不应未经整理直接和主数据混合。更适合作 external validation 或独立 benchmark。 |

细胞/QC 摘要：

| sample | cells | median genes | median counts | mt% median | ribo% median | low-quality% |
|---|---:|---:|---:|---:|---:|---:|
| HEFs-0330 | 4,335 | 6,018 | 41,546 | 4.48 | 18.06 | 1.78 |
| hADSCs-0618 | 4,507 | 5,781 | 38,420 | 5.90 | 19.17 | 4.37 |
| 0618-stageI | 4,931 | 4,099 | 22,901 | 4.41 | 34.65 | 10.65 |
| 0618-stageII | 5,981 | 3,254 | 12,831 | 4.97 | 35.44 | 5.23 |
| 0618-stageIII | 6,930 | 3,309 | 10,404 | 10.20 | 22.24 | 7.27 |
| 0330-stageII | 4,501 | 4,306 | 18,291 | 5.73 | 28.62 | 2.84 |
| 0330-stageII-5azaC | 3,740 | 4,421 | 18,864 | 4.77 | 28.83 | 2.22 |
| 0330-stageII-JNKIN8 | 2,919 | 5,214 | 28,587 | 6.99 | 25.16 | 9.18 |
| 0330-stageII-Tranyl | 3,911 | 4,360 | 18,798 | 5.96 | 28.30 | 2.89 |
| ASF | 5,475 | 5,470 | 34,060 | 5.25 | 20.42 | 2.10 |
| 0809_P3 | 6,129 | 6,367 | 44,101 | 5.10 | 18.77 | 1.75 |
| S1D16 | 4,653 | 4,261 | 24,289 | 9.08 | 27.84 | 29.23 |
| S2D24 | 6,706 | 4,028 | 17,187 | 7.57 | 28.01 | 6.73 |
| 0618_S1D4 | 4,415 | 5,474 | 37,442 | 8.27 | 23.21 | 1.29 |
| 0618_S2D4 | 5,605 | 4,142 | 20,949 | 5.80 | 36.28 | 5.14 |
| 0618_S2D8 | 3,005 | 4,895 | 24,387 | 4.74 | 31.15 | 4.16 |
| 0618_S2D12 | 4,841 | 3,632 | 16,299 | 4.11 | 37.31 | 3.08 |
| 0618_S3D4 | 4,520 | 4,094 | 15,390 | 8.67 | 23.05 | 5.11 |
| 0618_S3D8 | 3,605 | 4,166 | 15,938 | 10.10 | 20.74 | 8.40 |
| 0618_S4D1 | 7,865 | 4,305 | 13,682 | 8.76 | 18.01 | 3.13 |
| 0618_S4D2 | 10,649 | 3,164 | 8,562 | 10.06 | 17.88 | 7.41 |
| 0618_S4D4 | 9,587 | 1,151 | 2,655 | 6.80 | 25.04 | 2.61 |
| 0618_S4D10 | 3,868 | 5,343 | 29,424 | 6.29 | 23.62 | 13.73 |
| ADSC0618 | 18,814 | 472 | 1,050 | 30.47 | 19.09 | 59.46 |
| ADSC0809 | 5,524 | 4,054 | 14,659 | 2.94 | 25.76 | 2.75 |
| HEF1117 | 18,823 | 538 | 1,204 | 31.29 | 16.97 | 64.78 |
| ASF1230 | 20,414 | 604 | 1,266 | 17.49 | 24.12 | 46.10 |
| H1 | 16,530 | 666 | 1,297 | 8.18 | 30.15 | 3.73 |
| S1D0 | 5,797 | 5,357 | 30,142 | 4.04 | 22.48 | 1.12 |
| S1D0.5 | 4,130 | 6,371 | 42,587 | 4.94 | 21.37 | 1.40 |
| S1D1 | 4,346 | 5,872 | 32,771 | 7.22 | 19.65 | 0.76 |
| S1D2 | 4,860 | 5,177 | 30,282 | 5.77 | 21.15 | 2.35 |

## 3. `GSE298212(human`

| 问题 | 回答 |
|---|---|
| 1. 生物学背景 | human blood-cell chemical reprogramming to hCiPS cells。GEO design 写明 5 个 scRNA samples: PBMC EPC, S1D1, S1D3, S1D6, S1D8。 |
| 2. 时间信息 | PBMC_EPC 起始/参考样本，加 S1D1、S1D3、S1D6、S1D8。时间轴较短，覆盖早期 stage 1。 |
| 3. 原始格式和矩阵类型 | 5 个 10x `matrix.mtx.gz` + `barcodes.tsv.gz` + `features.tsv.gz`。矩阵是 integer raw counts。features 第 1 列 Ensembl ID，第 2 列 gene symbol。没有 `.h5ad/.h5/loom/Seurat object`。 |
| 4. 批次信息 | 每个样本就是一个 library/batch；platform GPL24676，10x Genomics，Cell Ranger v5.0.1，hg38。time 与 batch 完全重合。series matrix 有 library name，但无独立 donor/replicate 字段。 |
| 5. 细胞数量和时间分布 | 总 40,794 cells。D1 只有 3,669 cells，D6 4,877 cells；D3 最高 17,460 cells。 |
| 6. QC 指标 | 已全量计算 detected genes、UMI/counts、mt%、ribosomal%、低质量比例。doublet 无本地字段。 |
| 7. 原始细胞类型注释 | 有 sample 级 cell type，例如 EPC/self-derived S1Dx；无 per-cell cell type annotation。 |
| 8. Ground truth / 外部参考 | 有实验时间点和 PBMC/EPC 起点；没有 per-cell fate label、lineage graph、终点判定标准文件。marker genes 需外部整理。 |
| 9. RNA velocity | 无 `.loom` 和 spliced/unspliced matrix；本地不支持直接 velocity。 |
| 10. 是否适合作为本项目 input | 有条件适合。它是 human chemical reprogramming，但时间点少、缺终末 hCiPSC 长时间轴，更适合作补充/外部测试。 |

| sample | cells | median genes | median counts | mt% median | ribo% median | low-quality% |
|---|---:|---:|---:|---:|---:|---:|
| PBMC_EPC | 9,007 | 1,531 | 2,852 | 2.97 | 18.93 | 1.42 |
| S1D1 | 3,669 | 4,400 | 18,304 | 5.05 | 23.28 | 4.33 |
| S1D3 | 17,460 | 1,076 | 1,773 | 6.78 | 8.53 | 4.30 |
| S1D6 | 4,877 | 7,760 | 43,061 | 6.53 | 14.22 | 7.96 |
| S1D8 | 5,781 | 7,884 | 38,813 | 4.23 | 13.43 | 8.80 |

## 4. `FCR_iPSC(mouse`

| 问题 | 回答 |
|---|---|
| 1. 生物学背景 | mouse fast chemical reprogramming from MEF to iPSC，包含 scRNA、ATAC、CUT&Tag、bulk RNA 多组学。 |
| 2. 时间信息 | scRNA: D0/D4/D8/D12/iPSC；ATAC/CUT&Tag/bulk RNA 包含 D0/D2/D4/D6/D8/D10/D12/iPSC，bulk RNA 另有 ESC。 |
| 3. 原始格式和矩阵类型 | scRNA 同时有 10x MTX raw counts 和 `.h5ad`。`GSE218855_FCR_scRNA_raw_qc.h5ad` 中 `X` 和 `layers['counts']` 都是整数 count-like values。`.h5ad` var names 为 mouse gene symbols。bulk RNA 是 FPKM，不是 raw count。ATAC/CUT&Tag 是 bigWig/peaks，不是表达矩阵。 |
| 4. 批次信息 | `.h5ad` 有 `sample/time/line/condition`；scRNA sample 与 time 重合。ATAC/CUT&Tag/bulk RNA 有 replicate suffix，例如 `_1/_2`。无 human donor 信息。 |
| 5. 细胞数量和时间分布 | 10x raw MTX 合计 54,941 barcodes；`.h5ad` QC 后为 46,361 cells。D8/D12/iPSC 相对较少但仍有 7,535-8,790 cells。 |
| 6. QC 指标 | `.h5ad` 已有 `n_genes_by_counts`、`total_counts`、`pct_counts_mt`。本地没有 doublet 字段；ribosomal% 不在 obs 中，但可从 gene symbols 重新计算。 |
| 7. 原始细胞类型注释 | `.h5ad` 没有 per-cell cell type，只提供 sample/time/line/condition。bulk series 有 sample 级 cell type，例如 MEF/reprogramming intermediate/iPSC/mESC。 |
| 8. Ground truth / 外部参考 | 有实验时间点和 iPSC/ESC endpoint；没有 per-cell fate label。FCR 论文提供 marker/假设，但本地 raw 文件没有冻结 lineage model。 |
| 9. RNA velocity | 无 `.loom`，无 spliced/unspliced matrix。10x data 理论上可从 BAM 生成 velocity，本地文件不够。 |
| 10. 是否适合作为本项目 input | 不适合作为 human trajectory 项目主 input，因为是 mouse。可用于跨物种/方法 sanity check。 |

`.h5ad` QC 摘要：

| sample | cells | median genes | median counts | mt% median | low-quality% |
|---|---:|---:|---:|---:|---:|
| FCR_D0 | 9,918 | 4,457 | 17,091 | 3.67 | 0.00 |
| FCR_D4 | 12,459 | 3,247 | 9,539 | 7.04 | 0.00 |
| FCR_D8 | 7,659 | 3,474 | 10,889 | 6.18 | 0.00 |
| FCR_D12 | 7,535 | 4,203 | 14,131 | 4.53 | 0.00 |
| FCR_iPSC | 8,790 | 4,178 | 19,126 | 4.91 | 0.00 |

## 5. `gse242424`

| 问题 | 回答 |
|---|---|
| 1. 生物学背景 | human fibroblast reprogramming / pluripotency-related dataset；从文件名看是 D0-D14/iPSC time course。 |
| 2. 时间信息 | D0, D2, D4, D6, D8, D10, D12, D14, iPSC。 |
| 3. 原始格式和矩阵类型 | 本地只有 `*.matrix.mtx.gz` 和 `*.barcodes.tsv.gz`，没有对应 `features.tsv/genes.tsv`。MTX header 为 integer raw count matrix，但基因名/ID 无法本地确认。没有 `.h5ad/.h5/Seurat/CSV/loom`。 |
| 4. 批次信息 | 每个 sample/timepoint 是一个 batch；本地无 donor/replicate/library/platform 元数据。time 与 batch 完全重合。 |
| 5. 细胞数量和时间分布 | 本地 barcode 数异常巨大：总 22,236,487 barcode entries，单样本约 2.17-2.87M。这不像过滤后的 cells，更像 unfiltered matrix/droplet barcode 列数；不能直接当细胞数。 |
| 6. QC 指标 | 由于本地缺 features/gene annotation，不能计算 mt%、ribo%；barcode 数过大且疑似 unfiltered，不建议直接计算 per-cell QC。doublet 无字段。 |
| 7. 原始细胞类型注释 | 本地没有 per-cell annotation。 |
| 8. Ground truth / 外部参考 | 有实验 timepoint/iPSC endpoint 线索；本地无 marker、fate label、lineage graph。 |
| 9. RNA velocity | 无 `.loom`，无 spliced/unspliced matrix。 |
| 10. 是否适合作为本项目 input | 暂不适合。需要补齐 features/genes 文件和过滤后的 cell barcode/metadata 后才能评估。 |

本地 MTX barcode 分布：

| time | barcode entries |
|---|---:|
| D0 | 2,517,452 |
| D2 | 2,468,344 |
| D4 | 2,413,356 |
| D6 | 2,173,461 |
| D8 | 2,384,334 |
| D10 | 2,323,127 |
| D12 | 2,869,519 |
| D14 | 2,364,973 |
| iPSC | 2,721,921 |

## 6. `gse247600`

| 问题 | 回答 |
|---|---|
| 1. 生物学背景 | iPSC-derived neuronal / midbrain dopaminergic neuron differentiation context，文件名包含 HET clone 和 D21/D30/D45/D65；不是 hCiPSC reprogramming。 |
| 2. 时间信息 | D21, D30, D45, D65。HET1/HET2 clone 多时间点；D45 有 nuclei 和 whole-cell protocol difference。 |
| 3. 原始格式和矩阵类型 | `.csv.gz` gene-by-cell matrices；每个样本有普通 `genematrix` 和 `incl_introns_genematrix`。值为 integer counts。行名是 Ensembl IDs；没有 gene symbol annotation。没有 `.h5ad/.h5/matrix.mtx/Seurat/loom`。barcode 在每个 CSV 内唯一。 |
| 4. 批次信息 | batch/library 来自 sample 文件名，例如 HET1D21_S2。clone HET1/HET2、time、protocol partially confounded：D45 只有 HET1 且分 nuclei/whole。无 donor/replicate metadata 文件。 |
| 5. 细胞数量和时间分布 | 普通 genematrix 合计 1,022 cells；若把 incl_introns 也算入会重复成 2,044 columns。D45 whole 只有 80 cells，是明显低细胞数样本。 |
| 6. QC 指标 | 已计算 detected genes 和 total counts。因为只有 Ensembl ID、无 gene symbol/annotation，本地不能可靠计算 mt%/ribo%。doublet 无字段。 |
| 7. 原始细胞类型注释 | 无 per-cell cell type annotation。 |
| 8. Ground truth / 外部参考 | 有实验 day/clone/protocol；无 fate label、lineage graph、终点标准文件。 |
| 9. RNA velocity | 有 intron-inclusive matrix，但没有 matched spliced/unspliced layers；没有 `.loom`。不能直接用于 scVelo。 |
| 10. 是否适合作为本项目 input | 不适合作为本 hCiPSC trajectory 项目 input，生物学体系不同且细胞数很少。 |

普通 genematrix QC：

| sample | cells | median genes | median counts | low-quality% |
|---|---:|---:|---:|---:|
| HET1D21_S2 | 118 | 6,200 | 173,558 | 1.69 |
| HET2D21_S1 | 115 | 6,022 | 113,740 | 1.74 |
| HET1D30_S5 | 147 | 6,966 | 455,851 | 0.00 |
| HET2D30_S6 | 150 | 6,839 | 516,430 | 0.00 |
| HET1D45NUCLEI_S8 | 143 | 10,919 | 1,128,753 | 0.00 |
| HET1D45WHOLE_S7 | 80 | 8,247 | 1,599,288 | 0.00 |
| HET1D65_S9 | 137 | 7,120 | 642,250 | 0.73 |
| HET2D65_S10 | 132 | 6,788 | 462,487 | 0.00 |

## 7. `gse175634`

| 问题 | 回答 |
|---|---|
| 1. 生物学背景 | human iPSC-to-cardiomyocyte differentiation dataset，GEO 描述为多 cell line、多 time point 的心肌分化。 |
| 2. 时间信息 | GEO 报道 7 个 timepoints / 16-day differentiation；但本地只有一个 count matrix，无法把 cell 映射到 timepoint。 |
| 3. 原始格式和矩阵类型 | 本地只有 `GSE175634_cell_counts.mtx`。MatrixMarket header: integer raw count matrix，38,943 genes x 230,786 cells，317,416,972 nonzero entries。没有本地 barcode/gene metadata，没有 `.h5ad/.h5/Seurat/CSV/loom`。 |
| 4. 批次信息 | 本地无法确认每个 cell/sample 来自哪个 batch、donor、replicate、library 或 platform；metadata 文件缺失。 |
| 5. 细胞数量和时间分布 | 矩阵列数 230,786，可视为 cell/barcode 数；但无法按时间点或样本分布。 |
| 6. QC 指标 | 可从 mtx 计算 detected genes/counts，但由于缺 gene annotation，无法本地计算 mt%/ribo%；也无法按 sample/time 分层。doublet 无字段。 |
| 7. 原始细胞类型注释 | 本地没有。 |
| 8. Ground truth / 外部参考 | 可能有实验 timepoint/line metadata，但本地缺失；无 fate label/lineage graph。 |
| 9. RNA velocity | 无 `.loom`，无 spliced/unspliced matrix。 |
| 10. 是否适合作为本项目 input | 暂不适合。即使补齐 metadata，它也是 cardiomyocyte differentiation，不是 hCiPSC reprogramming。 |

## 8. `gse280956(human`

| 问题 | 回答 |
|---|---|
| 1. 生物学背景 | human mesenchymal stem cell hydrogel / H3K18la ChIP-seq，cartilage regeneration context。 |
| 2. 时间信息 | 文件名有 D1/D2，但 series matrix 显示是 hydrogel dynamic condition/Input/IP，不是 scRNA temporal trajectory。 |
| 3. 原始格式和矩阵类型 | 本地是 `.bigwig` ChIP-seq signal 和 series matrix；没有表达矩阵、没有 single-cell barcode。 |
| 4. 批次信息 | 有 Input/IP 和 hydrogel condition；不是 scRNA batch。 |
| 5. 细胞数量和时间分布 | 不适用，没有细胞矩阵。 |
| 6. QC 指标 | 不适用，不能计算 scRNA per-cell genes/UMI/mt/ribo/doublet。 |
| 7. 原始细胞类型注释 | sample-level source 是 hMSC；无 per-cell annotation。 |
| 8. Ground truth / 外部参考 | 有 ChIP-seq biological condition；无 trajectory ground truth。 |
| 9. RNA velocity | 不适用。 |
| 10. 是否适合作为本项目 input | 不适合。它不是 RNA expression trajectory input。 |

## 9. `gse136314`

| 问题 | 回答 |
|---|---|
| 1. 生物学背景 | 本地空目录，无法从本地文件确认。 |
| 2. 时间信息 | 无本地文件。 |
| 3. 原始格式和矩阵类型 | 无本地文件。 |
| 4. 批次信息 | 无本地文件。 |
| 5. 细胞数量和时间分布 | 无本地文件。 |
| 6. QC 指标 | 无本地文件。 |
| 7. 原始细胞类型注释 | 无本地文件。 |
| 8. Ground truth / 外部参考 | 无本地文件。 |
| 9. RNA velocity | 无本地文件。 |
| 10. 是否适合作为本项目 input | 不能评估。 |

## RNA velocity 总结

| 数据集 | loom | spliced/unspliced | 是否直接适合 velocity |
|---|---|---|---|
| `gse230659(human` | 无 | 无 | 否 |
| `gse178325_human` | 无 | 无 | 否 |
| `GSE298212(human` | 无 | 无 | 否 |
| `FCR_iPSC(mouse` | 无 | 无 | 否 |
| `gse242424` | 无 | 无 | 否 |
| `gse247600` | 无 | 无；只有 incl_introns 总矩阵 | 否 |
| `gse175634` | 无 | 无 | 否 |
| `gse280956(human` | 不适用 | 不适用 | 否 |

## 本轮生成的审计文件

| 文件 | 内容 |
|---|---|
| `logs\raw_data_samples_summary.json` | 所有本地 MTX/CSV 样本格式、维度、barcode、gene ID 类型概览 |
| `logs\gse230659_primary_qc.json` | `gse230659` 15 个 scRNA 样本逐细胞 QC 汇总 |
| `logs\gse178325_qc_summary.json` | `gse178325_human` 32 个唯一 scRNA 样本逐细胞 QC 汇总 |
| `logs\other_dataset_qc_summary.json` | `GSE298212` 和 FCR `.h5ad` QC 汇总 |
| `logs\gse247600_qc_summary.json` | `gse247600` CSV 普通 genematrix QC 汇总 |

## 外部来源

- NCBI GEO GSE230659: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE230659
- NCBI GEO GSE178325: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178325
- NCBI GEO GSE298212: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE298212
- NCBI GEO GSE218855: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218855
- NCBI GEO GSE242424: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242424
- NCBI GEO GSE247600: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247600
- NCBI GEO GSE175634: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE175634
- NCBI GEO GSE280956: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280956
