# B03_qc_filter.py
import pandas as pd
import scanpy as sc
from pathlib import Path

ROOT = Path(r"E:\scgpt")
DATASET_TAG = "gse178325"

INDIR = ROOT / rf"data\processed_{DATASET_TAG}\merged"

IN_MERGED = INDIR / "GSE178325_merged_rawcounts.h5ad"
OUT_QC = INDIR / "GSE178325_merged_qc.h5ad"
OUT_STATS = INDIR / "qc_filter_stats_gse178325.tsv"

# 阈值
MIN_UMI = 1000
MIN_GENES = 500
MAX_MT = 20.0

# 内存紧张就用 True（你现在就是这个情况）
USE_BACKED = True


def first_nonnull(s):
    s2 = s.dropna()
    return s2.iloc[0] if len(s2) else None


def main():
    adata = sc.read_h5ad(IN_MERGED, backed="r" if USE_BACKED else None)

    # --- 过滤前统计（按 GSM）---
    before = (
        adata.obs.groupby("gsm", observed=True)
        .size()
        .rename("cells_before")
        .to_frame()
    )

    # --- 检查 QC 列 ---
    required = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
    miss = [c for c in required if c not in adata.obs.columns]
    if miss:
        raise ValueError(f"[ERROR] Missing QC columns in merged.h5ad: {miss} (run B02 first)")

    # --- 过滤 mask（只用 obs，不会加载大矩阵）---
    mask = (
        (adata.obs["total_counts"] >= MIN_UMI) &
        (adata.obs["n_genes_by_counts"] >= MIN_GENES) &
        (adata.obs["pct_counts_mt"] <= MAX_MT)
    )

    # ✅ backed 模式下：统计用 obs 子集；写文件用 copy(filename=...)
    if USE_BACKED:
        obs_after = adata.obs.loc[mask].copy()

        # 直接把过滤后的 AnnData 写到新文件（不会走 .copy() 到内存）
        _ = adata[mask].copy(filename=str(OUT_QC))

        # 关闭原文件句柄
        try:
            adata.file.close()
        except Exception:
            pass
    else:
        adata_qc = adata[mask].copy()
        obs_after = adata_qc.obs
        adata_qc.write_h5ad(OUT_QC)

    # --- 过滤后统计（按 GSM）---
    after = (
        obs_after.groupby("gsm", observed=True)
        .size()
        .rename("cells_after")
        .to_frame()
    )

    stats = before.join(after, how="left").fillna(0)
    stats["cells_after"] = stats["cells_after"].astype(int)
    stats["filtered_out"] = stats["cells_before"] - stats["cells_after"]
    stats["keep_rate"] = (stats["cells_after"] / stats["cells_before"]).round(4)

    # 加上 meta：stage/day/donor
    meta_cols = ["stage", "day", "donor"]
    meta_missing = [c for c in meta_cols if c not in obs_after.columns]
    if meta_missing:
        raise ValueError(f"[ERROR] Missing meta columns in obs: {meta_missing}")

    meta = (
        obs_after.groupby("gsm", observed=True)[meta_cols]
        .agg(first_nonnull)
    )

    stats = meta.join(stats, how="right")

    # 写 stats 文件（带阈值注释）
    with open(OUT_STATS, "w", encoding="utf-8") as w:
        w.write(f"# MIN_UMI={MIN_UMI}\tMIN_GENES={MIN_GENES}\tMAX_MT={MAX_MT}\n")
    stats.to_csv(OUT_STATS, sep="\t", mode="a")

    print(f"[INFO] before cells: {int(stats['cells_before'].sum())}")
    print(f"[INFO] after  cells: {int(stats['cells_after'].sum())}")
    print(f"[OK] wrote qc: {OUT_QC}")
    print(f"[OK] wrote stats: {OUT_STATS}\n")

    print(stats.sort_values("keep_rate").head(10))


if __name__ == "__main__":
    main()
