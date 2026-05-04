#!/usr/bin/env python3
"""
Build the GSE178325 0618 post-QC raw/full-gene h5ad used by the scGPT path.

This script intentionally mirrors the GSE230659 scGPT source artifact:
post-QC cells, filtered full-gene sparse log-normalized X, GSE230659-style
obs/var columns, and no benchmark-only framework columns. The downstream
scGPT step then produces adata_scgpt_full.h5ad and adata_scgpt_annotated.h5ad.
"""

from __future__ import annotations

import gc
import gzip
import os
import time as _time
from datetime import datetime
from pathlib import Path
from typing import Optional

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import psutil
except ImportError:  # pragma: no cover - optional runtime telemetry only
    psutil = None


TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")
T0 = _time.time()


def _find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
        raise FileNotFoundError(f"TRAJ_PROJECT_ROOT={env!r} does not exist.")
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "data").exists():
            return candidate
    return here.parent


PROJECT_ROOT = _find_project_root()
DATA_DIR = PROJECT_ROOT / "data" / "gse178325_human" / "rna_seq_10x"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed" / "gse178325_human"
OUTPUT_H5AD = PROCESSED_DIR / "GSE178325_0618_raw_full_gene_benchmark_input.h5ad"
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

QC_PARAMS = {
    "min_genes_per_cell": 200,
    "max_genes_per_cell": 8000,
    "min_cells_per_gene": 3,
    "max_pct_mito": 20.0,
    "min_counts_per_cell": 500,
}
NORM_TARGET_SUM = 1e4
DATASET_ID = "GSE178325"
BATCH_ID = "0618"
STAGE_OFFSETS = {"StageI": 0.0, "StageII": 4.0, "StageIII": 16.0, "StageIV": 24.0}
HCIPSC_ABSOLUTE_DAY = 36.0

SAMPLES = [
    ("GSM5683317", "StageI", 0.0, "GSM5683317_StageI_Day0-0618"),
    ("GSM5683318", "StageI", 0.5, "GSM5683318_StageI_Day0.5-0618"),
    ("GSM5683319", "StageI", 1.0, "GSM5683319_StageI_Day1-0618"),
    ("GSM5683320", "StageI", 2.0, "GSM5683320_StageI_Day2-0618"),
    ("GSM5534144", "StageI", 4.0, "GSM5534144_StageI_Day4-0618"),
    ("GSM5534145", "StageII", 4.0, "GSM5534145_StageII_Day4-0618"),
    ("GSM5534146", "StageII", 8.0, "GSM5534146_StageII_Day8-0618"),
    ("GSM5534147", "StageII", 12.0, "GSM5534147_StageII_Day12-0618"),
    ("GSM5534148", "StageIII", 4.0, "GSM5534148_StageIII_Day4-0618"),
    ("GSM5534149", "StageIII", 8.0, "GSM5534149_StageIII_Day8-0618"),
    ("GSM5534150", "StageIV", 1.0, "GSM5534150_StageIV_Day1-0618"),
    ("GSM5534151", "StageIV", 2.0, "GSM5534151_StageIV_Day2-0618"),
    ("GSM5534152", "StageIV", 4.0, "GSM5534152_StageIV_Day4-0618"),
    ("GSM5534153", "StageIV", 10.0, "GSM5534153_StageIV_Day10-0618"),
    ("GSM5534154", "hCiPSCs", None, "GSM5534154_hCiPSCs-0618"),
]


def rss_gb() -> float:
    if psutil is None:
        return float("nan")
    return psutil.Process().memory_info().rss / 1e9


def sample_path(gsm: str) -> Path:
    return DATA_DIR / gsm


def abs_day_for(stage: str, day: Optional[float]) -> float:
    if stage == "hCiPSCs":
        return HCIPSC_ABSOLUTE_DAY
    if day is None:
        raise ValueError(f"Missing day for non-terminal stage {stage}.")
    return STAGE_OFFSETS[stage] + float(day)


def stage_day_label(stage: str, day: Optional[float]) -> str:
    if stage == "hCiPSCs":
        return "hCiPSCs"
    if float(day).is_integer():
        return f"{stage}_Day{int(day)}"
    return f"{stage}_Day{day:g}"


def read_10x_mtx(sample_dir: Path):
    bc_path = sample_dir / "barcodes.tsv.gz"
    ft_path = sample_dir / "features.tsv.gz"
    mtx_path = sample_dir / "matrix.mtx.gz"
    if not bc_path.exists() or not ft_path.exists() or not mtx_path.exists():
        raise FileNotFoundError(f"Missing 10x files under {sample_dir}")

    barcodes = pd.read_csv(bc_path, sep="\t", header=None, names=["barcode"])
    features = pd.read_csv(
        ft_path,
        sep="\t",
        header=None,
        names=["gene_id", "gene_name", "feature_type"],
    )

    rows_list, cols_list, data_list = [], [], []
    with gzip.open(str(mtx_path), "rt") as f:
        line = f.readline()
        while line.startswith("%"):
            line = f.readline()
        n_genes, n_cells, _ = map(int, line.strip().split())

        block_rows, block_cols, block_vals = [], [], []
        for line in f:
            r, c, v = line.split()
            block_rows.append(int(r) - 1)
            block_cols.append(int(c) - 1)
            block_vals.append(int(v))
            if len(block_rows) >= 5_000_000:
                rows_list.append(np.asarray(block_rows, dtype=np.int32))
                cols_list.append(np.asarray(block_cols, dtype=np.int32))
                data_list.append(np.asarray(block_vals, dtype=np.float32))
                block_rows, block_cols, block_vals = [], [], []
        if block_rows:
            rows_list.append(np.asarray(block_rows, dtype=np.int32))
            cols_list.append(np.asarray(block_cols, dtype=np.int32))
            data_list.append(np.asarray(block_vals, dtype=np.float32))

    row = np.concatenate(rows_list) if len(rows_list) > 1 else rows_list[0]
    col = np.concatenate(cols_list) if len(cols_list) > 1 else cols_list[0]
    data = np.concatenate(data_list) if len(data_list) > 1 else data_list[0]
    return row, col, data, n_genes, n_cells, barcodes, features


def compute_qc(row, col, data, n_cells, mt_mask, ribo_mask) -> pd.DataFrame:
    total_counts = np.bincount(col, weights=data.astype(np.float64), minlength=n_cells)
    n_genes_detected = np.bincount(col, minlength=n_cells)
    mt_entry = mt_mask[row]
    ribo_entry = ribo_mask[row]
    total_counts_mt = np.bincount(
        col[mt_entry], weights=data[mt_entry].astype(np.float64), minlength=n_cells
    )
    total_counts_ribo = np.bincount(
        col[ribo_entry], weights=data[ribo_entry].astype(np.float64), minlength=n_cells
    )
    pct_counts_mt = np.where(total_counts > 0, 100.0 * total_counts_mt / total_counts, 0.0)
    pct_counts_ribo = np.where(
        total_counts > 0, 100.0 * total_counts_ribo / total_counts, 0.0
    )
    return pd.DataFrame(
        {
            "n_genes_by_counts": n_genes_detected.astype(np.int32),
            "total_counts": total_counts.astype(np.float32),
            "total_counts_mt": total_counts_mt.astype(np.float32),
            "pct_counts_mt": pct_counts_mt.astype(np.float32),
            "total_counts_ribo": total_counts_ribo.astype(np.float32),
            "pct_counts_ribo": pct_counts_ribo.astype(np.float32),
        }
    )


print(f"[{TIMESTAMP}] GSE178325 raw/full-gene builder")
print(f"  Project root : {PROJECT_ROOT}")
print(f"  Data dir     : {DATA_DIR}")
print(f"  Output h5ad  : {OUTPUT_H5AD}")

ref_features = pd.read_csv(
    sample_path(SAMPLES[0][0]) / "features.tsv.gz",
    sep="\t",
    header=None,
    names=["gene_id", "gene_name", "feature_type"],
)
GENE_IDS = ref_features["gene_id"].astype(str).to_numpy()
GENE_NAMES = ref_features["gene_name"].astype(str).to_numpy()
FEATURE_TYPES = ref_features["feature_type"].astype(str).to_numpy()
N_GENES_REF = len(GENE_NAMES)
MT_MASK = np.array([g.startswith("MT-") for g in GENE_NAMES], dtype=bool)
RIBO_MASK = np.array([g.startswith("RPS") or g.startswith("RPL") for g in GENE_NAMES])

GENE_SYMBOL_CODES, UNIQUE_GENE_SYMBOLS = pd.factorize(GENE_NAMES, sort=False)
N_SYMBOLS = len(UNIQUE_GENE_SYMBOLS)
symbol_gene_ids = []
symbol_feature_types = []
symbol_mt = np.zeros(N_SYMBOLS, dtype=bool)
symbol_ribo = np.zeros(N_SYMBOLS, dtype=bool)
for symbol_idx, symbol in enumerate(UNIQUE_GENE_SYMBOLS):
    mask = GENE_SYMBOL_CODES == symbol_idx
    symbol_gene_ids.append(";".join(pd.unique(GENE_IDS[mask]).astype(str)))
    symbol_feature_types.append(";".join(pd.unique(FEATURE_TYPES[mask]).astype(str)))
    symbol_mt[symbol_idx] = bool(MT_MASK[mask].any())
    symbol_ribo[symbol_idx] = bool(RIBO_MASK[mask].any())

symbol_sum = np.zeros(N_SYMBOLS, dtype=np.float64)
symbol_n_cells = np.zeros(N_SYMBOLS, dtype=np.int64)
total_cells_passed = 0
obs_parts = []

duplicate_symbol_entries = int(N_GENES_REF - N_SYMBOLS)
print(f"  Reference features      : {N_GENES_REF:,}")
print(f"  Unique gene symbols     : {N_SYMBOLS:,}")
print(f"  Duplicate symbol entries: {duplicate_symbol_entries:,}")

print("\nPASS 1: QC and gene statistics")
for gsm, stage, day, sample_name in SAMPLES:
    t_start = _time.time()
    row, col, data, n_genes, n_cells, barcodes, features = read_10x_mtx(sample_path(gsm))
    if n_genes != N_GENES_REF:
        raise ValueError(f"Gene count mismatch in {gsm}: {n_genes} vs {N_GENES_REF}")
    qc = compute_qc(row, col, data, n_cells, MT_MASK, RIBO_MASK)
    keep = (
        (qc["n_genes_by_counts"] >= QC_PARAMS["min_genes_per_cell"])
        & (qc["n_genes_by_counts"] < QC_PARAMS["max_genes_per_cell"])
        & (qc["total_counts"] >= QC_PARAMS["min_counts_per_cell"])
        & (qc["pct_counts_mt"] < QC_PARAMS["max_pct_mito"])
    ).to_numpy()
    kept = np.where(keep)[0]
    abs_day = abs_day_for(stage, day)
    obs = pd.DataFrame(
        {
            "sample_id": sample_name,
            "stage": stage,
            "stage_day_label": stage_day_label(stage, day),
            "day_within_stage": np.nan if day is None else float(day),
            "abs_day": abs_day,
            "n_genes_by_counts": qc.loc[keep, "n_genes_by_counts"].to_numpy(),
            "total_counts": qc.loc[keep, "total_counts"].to_numpy(),
            "total_counts_mt": qc.loc[keep, "total_counts_mt"].to_numpy(),
            "pct_counts_mt": qc.loc[keep, "pct_counts_mt"].to_numpy(),
            "total_counts_ribo": qc.loc[keep, "total_counts_ribo"].to_numpy(),
            "pct_counts_ribo": qc.loc[keep, "pct_counts_ribo"].to_numpy(),
        },
        index=[f"{sample_name}_{b}" for b in barcodes.loc[keep, "barcode"].astype(str)],
    )
    obs_parts.append(obs)

    kept_rank = np.full(n_cells, -1, dtype=np.int32)
    kept_rank[kept] = np.arange(len(kept), dtype=np.int32)
    entry_mask = kept_rank[col] >= 0
    merged_counts = sp.coo_matrix(
        (
            data[entry_mask].astype(np.float32),
            (kept_rank[col[entry_mask]], GENE_SYMBOL_CODES[row[entry_mask]]),
        ),
        shape=(len(kept), N_SYMBOLS),
        dtype=np.float32,
    ).tocsr()
    symbol_sum += np.asarray(merged_counts.sum(axis=0)).ravel()
    symbol_n_cells += merged_counts.getnnz(axis=0)
    total_cells_passed += len(kept)
    print(
        f"  {gsm}: {n_cells:>6,} -> {len(kept):>6,} cells "
        f"abs_day={abs_day:<5.1f} [{_time.time() - t_start:.1f}s, {rss_gb():.1f}GB]"
    )
    del row, col, data, barcodes, features, qc, keep, kept, kept_rank, entry_mask
    del merged_counts
    gc.collect()

obs_adata = pd.concat(obs_parts)
obs_adata.index.name = "obs_names"
del obs_parts
gc.collect()

symbol_passes = symbol_n_cells >= QC_PARAMS["min_cells_per_gene"]
symbol_indices = np.where(symbol_passes)[0]
symbol_rank = np.full(N_SYMBOLS, -1, dtype=np.int32)
symbol_rank[symbol_indices] = np.arange(len(symbol_indices), dtype=np.int32)
var = pd.DataFrame(
    {
        "gene_ids": np.asarray(symbol_gene_ids, dtype=object)[symbol_indices],
        "feature_types": np.asarray(symbol_feature_types, dtype=object)[symbol_indices],
        "mt": symbol_mt[symbol_indices],
        "ribo": symbol_ribo[symbol_indices],
        "n_cells": symbol_n_cells[symbol_indices].astype(np.int32),
        "n_cells_by_counts": symbol_n_cells[symbol_indices].astype(np.int32),
        "mean_counts": (symbol_sum[symbol_indices] / total_cells_passed).astype(np.float32),
        "pct_dropout_by_counts": (
            100.0 * (1.0 - symbol_n_cells[symbol_indices] / total_cells_passed)
        ).astype(np.float32),
        "total_counts": symbol_sum[symbol_indices].astype(np.float32),
        "index": np.asarray(UNIQUE_GENE_SYMBOLS, dtype=object)[symbol_indices],
    },
    index=np.asarray(UNIQUE_GENE_SYMBOLS, dtype=object)[symbol_indices],
)
var.index.name = "var_names"
print(
    f"\nGene symbols after merged min-cell filter: "
    f"{len(symbol_indices):,}/{N_SYMBOLS:,}"
)
print(f"  Duplicate gene-symbol entries merged: {duplicate_symbol_entries:,}")

print("\nPASS 2: Build sparse full-gene log-normalized matrix")
rows_all, cols_all, data_all = [], [], []
cell_offset = 0
cell_totals = obs_adata["total_counts"].to_numpy(dtype=np.float64)
for gsm, stage, day, sample_name in SAMPLES:
    t_start = _time.time()
    row, col, data, n_genes, n_cells, barcodes, features = read_10x_mtx(sample_path(gsm))
    qc = compute_qc(row, col, data, n_cells, MT_MASK, RIBO_MASK)
    keep = (
        (qc["n_genes_by_counts"] >= QC_PARAMS["min_genes_per_cell"])
        & (qc["n_genes_by_counts"] < QC_PARAMS["max_genes_per_cell"])
        & (qc["total_counts"] >= QC_PARAMS["min_counts_per_cell"])
        & (qc["pct_counts_mt"] < QC_PARAMS["max_pct_mito"])
    ).to_numpy()
    kept = np.where(keep)[0]
    kept_rank = np.full(n_cells, -1, dtype=np.int32)
    kept_rank[kept] = np.arange(len(kept), dtype=np.int32)
    row_symbol = GENE_SYMBOL_CODES[row]
    entry_mask = (kept_rank[col] >= 0) & (symbol_rank[row_symbol] >= 0)
    target_rows = cell_offset + kept_rank[col[entry_mask]]
    target_cols = symbol_rank[row_symbol[entry_mask]]
    scale = NORM_TARGET_SUM / cell_totals[target_rows]
    values = (data[entry_mask].astype(np.float32) * scale).astype(np.float32)
    rows_all.append(target_rows.astype(np.int32))
    cols_all.append(target_cols.astype(np.int32))
    data_all.append(values)
    cell_offset += len(kept)
    print(f"  {gsm}: nnz={len(values):,} [{_time.time() - t_start:.1f}s, {rss_gb():.1f}GB]")
    del row, col, data, barcodes, features, qc, keep, kept, kept_rank, row_symbol
    del entry_mask, target_rows, target_cols, scale, values
    gc.collect()

if cell_offset != total_cells_passed:
    raise AssertionError(f"Cell count mismatch: {cell_offset} vs {total_cells_passed}")

X = sp.coo_matrix(
    (
        np.concatenate(data_all),
        (np.concatenate(rows_all), np.concatenate(cols_all)),
    ),
    shape=(total_cells_passed, len(symbol_indices)),
    dtype=np.float32,
).tocsr()
X.data = np.log1p(X.data).astype(np.float32, copy=False)
del rows_all, cols_all, data_all
gc.collect()

adata = ad.AnnData(X=X, obs=obs_adata, var=var)
adata.uns["preprocessing_summary"] = {
    "script": "scripts/build_gse178325_raw_full_gene_input.py",
    "reference_dataset_flow": "GSE230659 scGPT full-gene path",
    "dataset": DATASET_ID,
    "batch_id": BATCH_ID,
    "timestamp": TIMESTAMP,
    "qc_params": QC_PARAMS,
    "norm_target_sum": NORM_TARGET_SUM,
    "gene_filter": f"merged gene-symbol n_cells >= {QC_PARAMS['min_cells_per_gene']}",
    "duplicate_gene_symbol_policy": (
        "sum raw counts for features with the same gene symbol before "
        "total-count normalization and log1p"
    ),
    "n_reference_features": int(N_GENES_REF),
    "n_unique_gene_symbols": int(N_SYMBOLS),
    "n_duplicate_symbol_entries_merged": int(duplicate_symbol_entries),
    "x_transform": "log1p(total-count normalized sparse counts after gene-symbol merging)",
}
adata.write_h5ad(OUTPUT_H5AD)

elapsed = _time.time() - T0
print("\nDONE")
print(f"  Output      : {OUTPUT_H5AD}")
print(f"  Shape       : {adata.n_obs:,} cells x {adata.n_vars:,} genes")
print(f"  X sparse    : {sp.issparse(adata.X)}")
print(f"  Runtime     : {elapsed / 60:.1f} min")
print(f"  RSS         : {rss_gb():.2f} GB")
