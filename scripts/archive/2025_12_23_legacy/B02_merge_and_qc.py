# B02_merge_and_qc.py
import gc
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
from scipy import sparse

ROOT = Path(r"E:\scgpt")
DATASET_TAG = "gse178325"

INDIR = ROOT / rf"data\processed_{DATASET_TAG}\per_gsm_h5ad"
OUTDIR = ROOT / rf"data\processed_{DATASET_TAG}\merged"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_MERGED = OUTDIR / "GSE178325_merged_rawcounts.h5ad"
TMPDIR = OUTDIR / "_tmp_chunks"
TMPDIR.mkdir(parents=True, exist_ok=True)

CHUNK_SIZE = 8  # 8~12 之间一般比较稳


def downcast_sparse_counts(a):
    """尽量把 counts 压到 uint16（不够就 uint32），并保证 CSR"""
    X = a.layers.get("counts", a.X)
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)
    else:
        X = X.tocsr()

    # 只看非零 data 的最大值
    maxv = int(X.data.max()) if X.data.size else 0
    if maxv <= np.iinfo(np.uint16).max:
        X = X.astype(np.uint16)
    else:
        X = X.astype(np.uint32)

    a.X = X
    a.layers["counts"] = X
    return a


def concat_list(adatas, keys):
    """concat 一个列表，尽量省内存"""
    # 注意：这里用 join="inner" 避免 outer 对齐带来的额外内存
    adata = ad.concat(
        adatas,
        join="inner",
        label="gsm_batch",
        keys=keys,
    )
    # 保证 X 是 raw counts
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X
    adata.X = adata.layers["counts"]
    return adata


def main():
    files = sorted(INDIR.glob("GSM*.h5ad"))
    if not files:
        raise FileNotFoundError(f"No GSM*.h5ad found in {INDIR}")

    print(f"[INFO] found {len(files)} GSM files")

    tmp_files = []

    # ---------- chunk merge ----------
    for i in range(0, len(files), CHUNK_SIZE):
        chunk = files[i:i+CHUNK_SIZE]
        adatas = []
        keys = []

        print(f"\n[INFO] processing chunk {i//CHUNK_SIZE} ({i}..{i+len(chunk)-1})")

        for f in chunk:
            a = sc.read_h5ad(f)

            # 基础检查
            if "stage" not in a.obs.columns or "day" not in a.obs.columns:
                raise ValueError(f"[ERROR] Missing stage/day in {f}")

            # 确保 counts layer 存在
            if "counts" not in a.layers:
                a.layers["counts"] = a.X

            # downcast + csr
            a = downcast_sparse_counts(a)

            gsm = str(a.obs["gsm"].iloc[0])
            keys.append(gsm)
            adatas.append(a)

        merged_chunk = concat_list(adatas, keys)

        # 可选：把一些列设成 category（省内存）
        for c in ["gsm", "donor", "stage", "gsm_batch", "dataset"]:
            if c in merged_chunk.obs.columns:
                merged_chunk.obs[c] = merged_chunk.obs[c].astype("category")

        tmp = TMPDIR / f"chunk_{i//CHUNK_SIZE:02d}.h5ad"
        merged_chunk.write_h5ad(tmp)
        tmp_files.append(tmp)
        print(f"[OK] wrote chunk: {tmp}  shape={merged_chunk.shape}")

        # 释放内存
        del merged_chunk
        for a in adatas:
            del a
        del adatas
        gc.collect()

    # ---------- final merge (merge chunk files) ----------
    print(f"\n[INFO] merging {len(tmp_files)} chunks into final ...")
    adatas = []
    keys = []
    for tf in tmp_files:
        a = sc.read_h5ad(tf)
        if "counts" not in a.layers:
            a.layers["counts"] = a.X
        a = downcast_sparse_counts(a)

        # 用 chunk 文件名当 key 也行，这里不重要
        keys.append(tf.stem)
        adatas.append(a)

    adata = concat_list(adatas, keys)

    # QC metrics（仍然不做过滤）
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    gsm_counts = adata.obs.groupby("gsm", observed=True).size().sort_values(ascending=False)
    print("[INFO] per-gsm cell counts:")
    print(gsm_counts)

    cols = ["gsm", "stage", "day", "donor", "total_counts", "n_genes_by_counts", "pct_counts_mt"]
    print("\n[INFO] merged shape:", adata.shape)
    print(adata.obs[cols].head())

    adata.write_h5ad(OUT_MERGED)
    print(f"\n[OK] wrote merged: {OUT_MERGED}")

    # 如果你想删除 tmp chunk 文件，解除注释：
    # for tf in tmp_files:
    #     tf.unlink(missing_ok=True)
    # TMPDIR.rmdir()

if __name__ == "__main__":
    main()
