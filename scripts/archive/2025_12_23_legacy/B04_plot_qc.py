import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime

ROOT = Path(r"E:\scgpt")
DATASET_TAG = "gse178325"

# ✅ 新数据的输入目录
INDIR = ROOT / rf"data\processed_{DATASET_TAG}\merged"
IN_MERGED = INDIR / "GSE178325_merged_rawcounts.h5ad"
IN_QC = INDIR / "GSE178325_merged_qc.h5ad"

# ✅ 用到秒，避免同一天多次运行撞目录
ts = datetime.now().strftime("%Y%m%d_%H%M%S")
OUTDIR = INDIR / f"qc_plots_{DATASET_TAG}_{ts}"
OUTDIR.mkdir(parents=True, exist_ok=False)

# 画图用的 QC 指标（你 B02 已经算好了）
METRICS = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]

# ✅ stage 顺序：补齐 Somatic / StageIV / H1
STAGE_ORDER = ["Somatic", "StageI", "StageII", "StageIII", "StageIV", "hCiPSC", "H1"]

SCATTER_MAX_POINTS = 50000
RANDOM_SEED = 0


def ensure_qc_metrics_exist(adata, tag=""):
    """B04 不再现算 QC（避免 backed 模式触发大矩阵计算）；只检查是否存在。"""
    need = set(METRICS)
    miss = sorted(list(need - set(adata.obs.columns)))
    if miss:
        raise ValueError(
            f"[ERROR] missing QC columns {miss} in {tag}. "
            f"Please make sure you ran B02 (calculate_qc_metrics) and saved the h5ad."
        )
    return adata


def make_ordered_categories(adata):
    """给 stage/day 做排序用的 category。"""
    if "stage" in adata.obs.columns:
        adata.obs["stage"] = adata.obs["stage"].astype(str)
        stages = [s for s in STAGE_ORDER if s in set(adata.obs["stage"])]
        others = sorted(set(adata.obs["stage"]) - set(stages))
        stage_cat = stages + others
        adata.obs["stage_cat"] = pd.Categorical(
            adata.obs["stage"], categories=stage_cat, ordered=True
        )

    if "day" in adata.obs.columns:
        adata.obs["day"] = pd.to_numeric(adata.obs["day"], errors="coerce")
        days = sorted([d for d in adata.obs["day"].dropna().unique()])
        adata.obs["day_cat"] = pd.Categorical(adata.obs["day"], categories=days, ordered=True)

    return adata


def _boxplot_by_group(df, group_col, metric, title, out_png):
    groups = df[group_col].dropna().unique().tolist()
    if isinstance(df[group_col].dtype, pd.CategoricalDtype) and df[group_col].dtype.ordered:
        groups = [g for g in df[group_col].dtype.categories if g in set(groups)]
    else:
        groups = sorted(groups, key=lambda x: str(x))

    data, labels = [], []
    for g in groups:
        vals = df.loc[df[group_col] == g, metric].astype(float).values
        if len(vals) == 0:
            continue
        data.append(vals)
        labels.append(str(g))

    plt.figure(figsize=(max(10, 0.45 * len(labels)), 5))
    plt.boxplot(data, showfliers=False)
    plt.xticks(range(1, len(labels) + 1), labels, rotation=45, ha="right")
    plt.ylabel(metric)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def _scatter_by_stage(df, x, y, title, out_png, max_points=SCATTER_MAX_POINTS):
    rng = np.random.default_rng(RANDOM_SEED)

    if len(df) > max_points:
        idx = rng.choice(len(df), size=max_points, replace=False)
        d = df.iloc[idx].copy()
    else:
        d = df.copy()

    stage_col = "stage"
    stages = d[stage_col].dropna().unique().tolist()
    stages_sorted = [s for s in STAGE_ORDER if s in set(stages)] + sorted(set(stages) - set(STAGE_ORDER))

    plt.figure(figsize=(7, 5))
    for s in stages_sorted:
        sub = d[d[stage_col] == s]
        if len(sub) == 0:
            continue
        plt.scatter(sub[x], sub[y], s=3, alpha=0.35, label=str(s))

    plt.xlabel(x)
    plt.ylabel(y)
    plt.title(title)
    plt.legend(markerscale=3, fontsize=9, frameon=False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()


def plot_suite(adata, tag):
    adata = ensure_qc_metrics_exist(adata, tag=tag)
    adata = make_ordered_categories(adata)
    df = adata.obs.copy()

    # 1) boxplot by stage
    for m in METRICS:
        _boxplot_by_group(
            df=df,
            group_col="stage_cat",
            metric=m,
            title=f"{tag}: {m} by stage",
            out_png=OUTDIR / f"{tag}_box_by_stage_{m}.png",
        )

    # 2) boxplot by day
    for m in METRICS:
        _boxplot_by_group(
            df=df,
            group_col="day_cat",
            metric=m,
            title=f"{tag}: {m} by day",
            out_png=OUTDIR / f"{tag}_box_by_day_{m}.png",
        )

    # 3) scatter plots
    _scatter_by_stage(
        df=df,
        x="total_counts",
        y="n_genes_by_counts",
        title=f"{tag}: total_counts vs n_genes_by_counts (colored by stage)",
        out_png=OUTDIR / f"{tag}_scatter_counts_vs_genes_by_stage.png",
    )
    _scatter_by_stage(
        df=df,
        x="total_counts",
        y="pct_counts_mt",
        title=f"{tag}: total_counts vs pct_counts_mt (colored by stage)",
        out_png=OUTDIR / f"{tag}_scatter_counts_vs_mt_by_stage.png",
    )

    # 4) pct_counts_mt histogram
    plt.figure(figsize=(7, 4))
    vals = df["pct_counts_mt"].astype(float).values
    vals = vals[np.isfinite(vals)]
    plt.hist(vals, bins=80)
    plt.xlabel("pct_counts_mt")
    plt.ylabel("cells")
    plt.title(f"{tag}: pct_counts_mt histogram")
    plt.tight_layout()
    plt.savefig(OUTDIR / f"{tag}_hist_pct_counts_mt.png", dpi=200)
    plt.close()

    print(f"[OK] QC plots written to: {OUTDIR} ({tag})")


def main():
    # ✅ 用 backed 读，避免把超大 X 全载入内存（B04 只用 obs）
    print("[INFO] reading (backed):", IN_MERGED)
    ad_pre = sc.read_h5ad(IN_MERGED, backed="r")

    print("[INFO] reading (backed):", IN_QC)
    ad_post = sc.read_h5ad(IN_QC, backed="r")

    plot_suite(ad_pre, tag="pre")
    plot_suite(ad_post, tag="post")

    # 关闭文件句柄（有时 Windows 下不关会影响后续操作）
    try:
        ad_pre.file.close()
        ad_post.file.close()
    except Exception:
        pass

    print("\n[DONE] All QC plots generated.")
    print(f"Output folder: {OUTDIR}")


if __name__ == "__main__":
    main()
