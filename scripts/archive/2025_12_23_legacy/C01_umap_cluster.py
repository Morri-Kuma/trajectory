import scanpy as sc
import numpy as np
from pathlib import Path
from datetime import datetime

# =========================
# Config
# =========================
DATASET_TAG = "gse178325"

# ✅ 输入：你这次 QC 的输出文件
IN_H5AD = r"E:\scgpt\data\processed_gse178325\merged\GSE178325_merged_qc.h5ad"

# ✅ 输出：按数据集分开存，避免混到旧的 C_traj
BASE_RESULTS = Path(r"E:\scgpt\results\C_traj") / DATASET_TAG


def get_run_dir_always_new():
    """
    ✅ 每次都新建一个 run_dir（绝不复用 latest），避免覆盖历史结果。
    同时更新 _latest_run.txt 方便你快速定位最近一次输出。
    """
    BASE_RESULTS.mkdir(parents=True, exist_ok=True)
    latest = BASE_RESULTS / "_latest_run.txt"

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = BASE_RESULTS / run_id
    run_dir.mkdir(parents=True, exist_ok=False)

    latest.write_text(str(run_dir), encoding="utf-8")
    return run_dir


RUN = get_run_dir_always_new()
FIGDIR = RUN / "figures"
ADIR = RUN / "adata"
FIGDIR.mkdir(parents=True, exist_ok=True)
ADIR.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = str(FIGDIR)
sc.settings.autoshow = False

print("[INFO] reading:", IN_H5AD)
adata = sc.read_h5ad(IN_H5AD)

# --- 必要字段检查/准备 ---
for k in ["day", "stage"]:
    if k not in adata.obs:
        raise ValueError(f"adata.obs 里必须有 '{k}' 字段")

adata.obs["day_numeric"] = adata.obs["day"].astype(float)

# --- 确保 counts layer 存在（后面 HVG 建议用它）---
if "counts" not in adata.layers:
    # 如果没有 counts layer，就把当前 X 当作 counts（你这份一般是有的）
    adata.layers["counts"] = adata.X.copy()

# --- 规范化 & log1p（只在没做过时做）---
# 注意：normalize/log1p 会把 X 变成 float，建议用 float32 降内存
if "log1p" not in adata.uns:
    # ✅ 用 counts layer 生成一个 float32 的工作矩阵作为 X（不污染 counts）
    adata.X = adata.layers["counts"].astype(np.float32)

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

# --- HVG：优先 donor batch-aware；并且用 raw counts 来估计 HVG（layer='counts'）---
hvg_kwargs = dict(n_top_genes=3000, flavor="seurat_v3", layer="counts")
if "donor" in adata.obs and adata.obs["donor"].nunique() > 1:
    hvg_kwargs["batch_key"] = "donor"

# seurat_v3 偶尔在内存/版本上会出问题，给一个稳健兜底
try:
    sc.pp.highly_variable_genes(adata, **hvg_kwargs)
except Exception as e:
    print("[WARN] seurat_v3 HVG failed, fallback to cell_ranger. Error:", repr(e))
    hvg_kwargs = dict(n_top_genes=3000, flavor="cell_ranger")
    if "donor" in adata.obs and adata.obs["donor"].nunique() > 1:
        hvg_kwargs["batch_key"] = "donor"
    sc.pp.highly_variable_genes(adata, **hvg_kwargs)

adata = adata[:, adata.var["highly_variable"]].copy()

# --- PCA / neighbors / UMAP / Leiden ---
# 大数据避免 scale+zero_center=True（会吃内存）
sc.tl.pca(adata, n_comps=50, svd_solver="randomized", zero_center=False)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50)
sc.tl.umap(adata, min_dist=0.3)
sc.tl.leiden(adata, resolution=0.8, key_added="leiden")

# --- 画图：按 stage/day/cluster ---
sc.pl.umap(
    adata,
    color=["stage", "day_numeric"],
    ncols=2,
    wspace=0.4,
    save=f"_{DATASET_TAG}_stage_day.png"
)
sc.pl.umap(
    adata,
    color=["leiden"],
    save=f"_{DATASET_TAG}_cluster.png"
)

# --- 保存 ---
out = ADIR / f"C01_umap_cluster_{DATASET_TAG}.h5ad"
adata.write_h5ad(out)

# 记录一下参数，方便复现实验
(RUN / "run_config.txt").write_text(
    "\n".join([
        f"IN_H5AD={IN_H5AD}",
        "normalize_total=1e4",
        "log1p=True",
        "HVG=n_top_genes=3000",
        "PCA=n_comps=50 randomized zero_center=False",
        "neighbors=n_neighbors=20 n_pcs=50",
        "umap=min_dist=0.3",
        "leiden=resolution=0.8",
    ]),
    encoding="utf-8"
)

print("[OK] wrote:", out)
print("[OK] figures at:", FIGDIR)
print("[OK] run dir:", RUN)
