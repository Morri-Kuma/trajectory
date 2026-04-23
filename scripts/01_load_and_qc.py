#!/usr/bin/env python3
"""
01_load_and_qc.py — Load GSE230659 scRNA-seq data, perform QC, and write
                     the canonical scTimeBench benchmark input object.

Project : scTimeBench-aligned benchmark for human chemical iPSC reprogramming
Dataset : GSE230659 (Liuyang et al. 2023 Cell Stem Cell)
Step    : 1 — Data loading, QC filtering, normalization, HVG selection,
              PCA, AnnData construction, h5ad output
Author  : Kuma (Graduate School of Frontier Sciences, University of Tokyo)
Env     : conda activate traj_env  (Python 3.10+)
          Requires: numpy, pandas, matplotlib, seaborn, psutil,
                    anndata, scipy

Revision history:
    v1  2026-03-29  Initial version (requires scanpy/anndata)
    v2  2026-03-29  Memory optimization — inner join, per-sample QC
    v3  2026-03-29  Pure numpy/pandas, no anndata/scipy at runtime.
                    Outputs npz/csv + reload snippet.
    v4  2026-04-17  Refactored to scTimeBench v2 framework.
                    Directly constructs and writes AnnData h5ad.
                    No longer produces reload snippet.
                    Required obs fields added for benchmark compatibility.
                    Figures and QC summary redirected to benchmark/reports/.

Usage:
    conda activate traj_env
    cd C:\\Users\\37620\\trajectory
    python scripts/01_load_and_qc.py

Primary output (canonical benchmark input):
    data/processed/adata_benchmark.h5ad

Secondary outputs (timestamped archives):
    data/processed/{TIMESTAMP}_GSE230659_qc_obs.csv
    data/processed/{TIMESTAMP}_GSE230659_qc_var.csv
    data/processed/{TIMESTAMP}_GSE230659_qc_pca.npz
    benchmark/reports/qc/{TIMESTAMP}_qc_violin_post.png
    benchmark/reports/qc/{TIMESTAMP}_qc_scatter.png
    benchmark/reports/qc/{TIMESTAMP}_cell_counts_per_time.png
    benchmark/reports/qc/{TIMESTAMP}_pca_time.png
    benchmark/reports/qc/{TIMESTAMP}_qc_summary.csv
"""

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
import os, gc, gzip, time as _time
from datetime import datetime
from pathlib import Path
from typing import Optional, Tuple, Dict, List

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import psutil
import anndata as ad
import scipy.sparse as sp
import warnings
warnings.filterwarnings("ignore")

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")
T0 = _time.time()

# ---- Project-root detection ----
# Strategy: try the known absolute path first, then search upward from this
# script's location until we find a directory that contains a 'data/' subdirectory.

def _find_project_root() -> Path:
    # 1. Explicit override: TRAJ_PROJECT_ROOT env var (Shirokane HPC / CI).
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
        raise FileNotFoundError(
            f"TRAJ_PROJECT_ROOT={env!r} does not exist.  "
            "Correct the environment variable and retry."
        )
    # 2. Walk upward from this script until a directory containing 'data/' is found.
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "data").exists():
            return candidate
    # 3. Last-resort fallback.
    return here.parent

PROJECT_ROOT = _find_project_root()

DATA_DIR      = PROJECT_ROOT / "data" / "gse230659(human" / "rna_seq"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
FIGURES_DIR   = PROJECT_ROOT / "benchmark" / "reports" / "qc"
METRICS_DIR   = PROJECT_ROOT / "benchmark" / "reports" / "qc"

# Canonical benchmark output (stable, non-timestamped)
ADATA_BENCHMARK_PATH = PROCESSED_DIR / "adata_benchmark.h5ad"

for d in [PROCESSED_DIR, FIGURES_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ---- QC Parameters ----
QC_PARAMS = {
    "min_genes_per_cell": 200,
    "max_genes_per_cell": 8000,
    "min_cells_per_gene": 3,
    "max_pct_mito": 20.0,
    "min_counts_per_cell": 500,
}

# ---- Feature selection & dim reduction ----
NORM_TARGET_SUM = 1e4
N_HVG = 2000
N_PCS = 50

# ---- Dataset ID ----
DATASET_ID = "GSE230659"

print(f"[{TIMESTAMP}] GSE230659 QC Pipeline v4 (scTimeBench-aligned)")
print(f"  Project root : {PROJECT_ROOT}")
print(f"  Data dir     : {DATA_DIR}")
print(f"  Output h5ad  : {ADATA_BENCHMARK_PATH}")
print(f"  RAM available: {psutil.virtual_memory().available / 1e9:.1f} GB")

# =============================================================================
# 1. SAMPLE MANIFEST & TIME MAPPING
# =============================================================================
STAGE_OFFSETS = {"StageI": 0.0, "StageII": 8.0, "StageIII": 16.0}
HCIPSC_ABSOLUTE_DAY = 30.0

SAMPLES = [
    ("GSM7229998_StageI_Day0.5-0618",    "StageI",   0.5),
    ("GSM7229999_StageI_Day2-0618",       "StageI",   2.0),
    ("GSM7230000_StageI_Day4-0618",       "StageI",   4.0),
    ("GSM7230001_StageI_Day8-0618",       "StageI",   8.0),
    ("GSM7230002_StageII_Day4-0618",      "StageII",  4.0),
    ("GSM7230003_StageII_Day8-0618",      "StageII",  8.0),
    ("GSM7230004_StageIII_Day0.33-0618",  "StageIII", 0.33),
    ("GSM7230005_StageIII_Day0.67-0618",  "StageIII", 0.67),
    ("GSM7230006_StageIII_Day1-0618",     "StageIII", 1.0),
    ("GSM7230007_StageIII_Day2-0618",     "StageIII", 2.0),
    ("GSM7230008_StageIII_Day4-0618",     "StageIII", 4.0),
    ("GSM7230009_StageIII_Day6-0618",     "StageIII", 6.0),
    ("GSM7230010_StageIII_Day8-0618",     "StageIII", 8.0),
    ("GSM7230011_StageIII_Day12-0618",    "StageIII", 12.0),
    ("GSM7230012_hCiPSCs-0618",           "hCiPSCs",  None),
]

# Provisional stage → cell-type proxy mapping.
# Used for scTimeBench_cell_type when no external annotation file is available.
# This is a coarse stage-level proxy; it should be replaced with a frozen
# per-cell annotation when one is available.
STAGE_TO_CELLTYPE_PROXY: Dict[str, str] = {
    "StageI":   "somatic",
    "StageII":  "early_transition",
    "StageIII": "intermediate",
    "hCiPSCs":  "pluripotent",
}
CELLTYPE_ANNOTATION_IS_PROXY = True   # flip to False when real annotation is wired in

# =============================================================================
# 2. PURE-PYTHON 10X MTX READER
# =============================================================================
def read_10x_mtx(
    sample_dir: Path,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, int, pd.DataFrame, pd.DataFrame]:
    """
    Read 10x Genomics MTX format (matrix.mtx.gz + barcodes.tsv.gz + features.tsv.gz).
    Returns COO arrays (row, col, data) plus shape and metadata DataFrames.
    MTX is 1-indexed; we convert to 0-indexed.
    Returns gene-major COO: row=gene_idx, col=cell_idx, data=count.
    """
    bc_path = sample_dir / "barcodes.tsv.gz"
    barcodes = pd.read_csv(bc_path, sep="\t", header=None, names=["barcode"])

    ft_path = sample_dir / "features.tsv.gz"
    features = pd.read_csv(ft_path, sep="\t", header=None,
                           names=["gene_id", "gene_name", "feature_type"])

    mtx_path = sample_dir / "matrix.mtx.gz"
    rows_list, cols_list, data_list = [], [], []

    with gzip.open(str(mtx_path), "rt") as f:
        line = f.readline()
        while line.startswith("%"):
            line = f.readline()
        parts = line.strip().split()
        n_genes, n_cells, n_entries = int(parts[0]), int(parts[1]), int(parts[2])

        block_size = 5_000_000
        block_rows, block_cols, block_vals = [], [], []

        for line in f:
            parts = line.split()
            block_rows.append(int(parts[0]) - 1)
            block_cols.append(int(parts[1]) - 1)
            block_vals.append(int(parts[2]))

            if len(block_rows) >= block_size:
                rows_list.append(np.array(block_rows, dtype=np.int32))
                cols_list.append(np.array(block_cols, dtype=np.int32))
                data_list.append(np.array(block_vals, dtype=np.int32))
                block_rows, block_cols, block_vals = [], [], []

        if block_rows:
            rows_list.append(np.array(block_rows, dtype=np.int32))
            cols_list.append(np.array(block_cols, dtype=np.int32))
            data_list.append(np.array(block_vals, dtype=np.int32))

    row = np.concatenate(rows_list) if len(rows_list) > 1 else rows_list[0]
    col = np.concatenate(cols_list) if len(cols_list) > 1 else cols_list[0]
    data = np.concatenate(data_list) if len(data_list) > 1 else data_list[0]

    del rows_list, cols_list, data_list, block_rows, block_cols, block_vals
    gc.collect()

    return row, col, data, n_genes, n_cells, barcodes, features


def compute_cell_qc(
    row: np.ndarray, col: np.ndarray, data: np.ndarray,
    n_genes: int, n_cells: int, mt_mask: np.ndarray,
) -> pd.DataFrame:
    """Compute per-cell QC metrics from COO arrays using numpy bincount."""
    total_counts = np.bincount(col, weights=data.astype(np.float64), minlength=n_cells)
    n_genes_detected = np.bincount(col, minlength=n_cells)
    is_mito = mt_mask[row]
    mito_counts = np.bincount(col[is_mito], weights=data[is_mito].astype(np.float64),
                              minlength=n_cells)
    pct_mito = np.where(total_counts > 0, 100.0 * mito_counts / total_counts, 0.0)

    return pd.DataFrame({
        "n_genes_by_counts": n_genes_detected.astype(np.int32),
        "total_counts": total_counts.astype(np.float32),
        "pct_counts_mt": pct_mito.astype(np.float32),
    })


# =============================================================================
# 3. PASS 1 — Per-sample QC + gene statistics accumulation
# =============================================================================
print("\n" + "=" * 70)
print("PASS 1: Per-sample QC filtering + gene statistics")
print("=" * 70)

ref_features = pd.read_csv(DATA_DIR / SAMPLES[0][0] / "features.tsv.gz",
                           sep="\t", header=None,
                           names=["gene_id", "gene_name", "feature_type"])
GENE_NAMES = ref_features["gene_name"].values
GENE_IDS   = ref_features["gene_id"].values
N_GENES_REF = len(GENE_NAMES)

MT_MASK   = np.array([g.startswith("MT-") for g in GENE_NAMES], dtype=bool)
RIBO_MASK = np.array([g.startswith("RPS") or g.startswith("RPL") for g in GENE_NAMES],
                     dtype=bool)

print(f"  Reference genes: {N_GENES_REF:,}")
print(f"  MT genes: {MT_MASK.sum()}")
print(f"  Ribosomal genes: {RIBO_MASK.sum()}")

gene_sum     = np.zeros(N_GENES_REF, dtype=np.float64)
gene_sq_sum  = np.zeros(N_GENES_REF, dtype=np.float64)
gene_n_cells = np.zeros(N_GENES_REF, dtype=np.int64)
total_cells_passed = 0

sample_meta_list: List[pd.DataFrame] = []
sample_cell_counts: Dict[str, Dict] = {}

for folder, stage, day in SAMPLES:
    t_start = _time.time()
    sample_path = DATA_DIR / folder

    row, col, data, n_genes, n_cells, barcodes, features = read_10x_mtx(sample_path)
    assert n_genes == N_GENES_REF, f"Gene count mismatch in {folder}: {n_genes} vs {N_GENES_REF}"

    qc = compute_cell_qc(row, col, data, n_genes, n_cells, MT_MASK)

    keep = (
        (qc["n_genes_by_counts"] >= QC_PARAMS["min_genes_per_cell"])
        & (qc["n_genes_by_counts"] < QC_PARAMS["max_genes_per_cell"])
        & (qc["total_counts"] >= QC_PARAMS["min_counts_per_cell"])
        & (qc["pct_counts_mt"] < QC_PARAMS["max_pct_mito"])
    ).values

    n_kept    = keep.sum()
    n_removed = n_cells - n_kept

    is_ribo = RIBO_MASK[row]
    ribo_counts_all   = np.bincount(col[is_ribo], weights=data[is_ribo].astype(np.float64),
                                    minlength=n_cells)
    total_counts_all  = np.bincount(col, weights=data.astype(np.float64), minlength=n_cells)
    pct_ribo = np.where(total_counts_all > 0, 100.0 * ribo_counts_all / total_counts_all, 0.0)

    abs_day         = STAGE_OFFSETS.get(stage, 0.0) + day if day is not None else HCIPSC_ABSOLUTE_DAY
    stage_day_label = f"{stage}_Day{day}" if day is not None else "hCiPSCs"

    # scTimeBench_cell_type: stage-based proxy (no external annotation file available yet).
    # Replace the STAGE_TO_CELLTYPE_PROXY lookup with a per-cell annotation merge
    # once a frozen annotation file has been defined.
    cell_type_proxy = STAGE_TO_CELLTYPE_PROXY.get(stage, "unknown")

    cell_meta = pd.DataFrame({
        # --- Original columns (kept from v3) ---
        "barcode":             [f"{folder}_{b}" for b in barcodes["barcode"].values[keep]],
        "sample_id":           folder,
        "stage":               stage,
        "stage_day_label":     stage_day_label,
        "day_within_stage":    day if day is not None else np.nan,
        "abs_day":             abs_day,
        "n_genes_by_counts":   qc["n_genes_by_counts"].values[keep],
        "total_counts":        qc["total_counts"].values[keep],
        "pct_counts_mt":       qc["pct_counts_mt"].values[keep],
        "pct_counts_ribo":     pct_ribo[keep].astype(np.float32),
        # --- New columns required by scTimeBench framework ---
        "cell_id":                f"{folder}_placeholder",   # overwritten below per-cell
        "dataset_id":             DATASET_ID,
        "time_label":             abs_day,
        "scTimeBench_timepoint":  abs_day,
        "scTimeBench_cell_type":  cell_type_proxy,
    })

    # Set cell_id = barcode (per-cell, same as barcode column)
    cell_meta["cell_id"] = cell_meta["barcode"]

    sample_meta_list.append(cell_meta)

    # Accumulate gene stats from kept cells only
    kept_cell_indices = set(np.where(keep)[0])
    entry_mask = np.array([c in kept_cell_indices for c in col], dtype=bool)
    kept_row   = row[entry_mask]
    kept_data  = data[entry_mask].astype(np.float64)

    gene_sum     += np.bincount(kept_row, weights=kept_data, minlength=N_GENES_REF)
    gene_sq_sum  += np.bincount(kept_row, weights=kept_data ** 2, minlength=N_GENES_REF)
    gene_n_cells += np.bincount(kept_row, minlength=N_GENES_REF)

    total_cells_passed += n_kept

    sample_cell_counts[folder] = {
        "stage": stage, "abs_day": abs_day,
        "cells_before": n_cells, "cells_after": n_kept,
    }

    elapsed = _time.time() - t_start
    mem_gb  = psutil.Process().memory_info().rss / 1e9
    print(f"  {folder}: {n_cells:>6,} → {n_kept:>6,} cells "
          f"(dropped {n_removed:>4})  abs_day={abs_day:<6.2f}  "
          f"[{elapsed:.1f}s, {mem_gb:.1f}GB RSS]")

    del row, col, data, barcodes, features, qc, kept_row, kept_data, entry_mask
    del ribo_counts_all, total_counts_all, pct_ribo, keep
    gc.collect()

# Merge all cell metadata
obs = pd.concat(sample_meta_list, ignore_index=True)
obs.index = obs["barcode"]
del sample_meta_list
gc.collect()

print(f"\nPass 1 complete: {total_cells_passed:,} cells passed QC across {len(SAMPLES)} samples")
print(f"  RAM used: {psutil.Process().memory_info().rss / 1e9:.2f} GB")

# =============================================================================
# 4. GENE FILTERING + HVG SELECTION
# =============================================================================
print("\n" + "=" * 70)
print("GENE FILTERING + HVG SELECTION")
print("=" * 70)

gene_passes   = gene_n_cells >= QC_PARAMS["min_cells_per_gene"]
n_genes_before = N_GENES_REF
n_genes_after  = gene_passes.sum()
print(f"  Gene filter (min {QC_PARAMS['min_cells_per_gene']} cells): "
      f"{n_genes_before:,} → {n_genes_after:,}")

gene_mean = gene_sum / total_cells_passed
gene_var  = (gene_sq_sum / total_cells_passed) - gene_mean ** 2
gene_var  = np.maximum(gene_var, 0)

valid    = gene_passes & (gene_mean > 0)
log_mean = np.full(N_GENES_REF, np.nan)
log_var  = np.full(N_GENES_REF, np.nan)
log_mean[valid] = np.log10(gene_mean[valid] + 1e-10)
log_var[valid]  = np.log10(gene_var[valid]  + 1e-10)

n_bins        = 20
valid_idx     = np.where(valid)[0]
lm            = log_mean[valid]
lv            = log_var[valid]
bin_edges     = np.percentile(lm, np.linspace(0, 100, n_bins + 1))
bin_edges[0]  -= 1e-10
bin_edges[-1] += 1e-10
bin_assignment = np.clip(np.digitize(lm, bin_edges) - 1, 0, n_bins - 1)

expected_var = np.zeros(len(lm))
for b in range(n_bins):
    mask = bin_assignment == b
    if mask.sum() > 0:
        expected_var[mask] = np.median(lv[mask])

residual_var = np.full(N_GENES_REF, -np.inf)
residual_var[valid_idx] = lv - expected_var
residual_var[~gene_passes] = -np.inf

hvg_indices = np.argsort(residual_var)[::-1][:N_HVG]
hvg_mask    = np.zeros(N_GENES_REF, dtype=bool)
hvg_mask[hvg_indices] = True

print(f"  HVGs selected: {hvg_mask.sum()}")
print(f"  HVG mean expression range: [{gene_mean[hvg_mask].min():.4f}, {gene_mean[hvg_mask].max():.2f}]")

# Build var DataFrame (all genes)
var = pd.DataFrame({
    # --- Required benchmark fields ---
    "gene_symbol":      GENE_NAMES,          # required by framework v2
    "gene_id":          GENE_IDS,            # required by framework v2
    "highly_variable":  hvg_mask,            # required by framework v2
    "is_hvg":           hvg_mask,            # alias, required by framework v2
    # --- Original fields (kept) ---
    "gene_name":        GENE_NAMES,
    "mt":               MT_MASK,
    "ribo":             RIBO_MASK,
    "n_cells":          gene_n_cells.astype(np.int32),
    "mean_counts":      gene_mean.astype(np.float32),
    "var_counts":       gene_var.astype(np.float32),
    "passes_filter":    gene_passes,
})
var.index = GENE_NAMES

# Save full gene metadata (archive)
var_path = PROCESSED_DIR / f"{TIMESTAMP}_GSE230659_qc_var.csv"
var.to_csv(var_path)
print(f"  Saved gene metadata → {var_path}")

# =============================================================================
# 5. PASS 2 — Build HVG expression matrix (dense)
# =============================================================================
print("\n" + "=" * 70)
print("PASS 2: Build HVG expression matrix")
print("=" * 70)

hvg_gene_indices = np.where(hvg_mask)[0]
hvg_idx_set      = set(hvg_gene_indices.tolist())
hvg_col_map      = {g: i for i, g in enumerate(hvg_gene_indices)}

X_hvg = np.zeros((total_cells_passed, N_HVG), dtype=np.float32)
print(f"  Allocated HVG matrix: {total_cells_passed:,} x {N_HVG} = "
      f"{X_hvg.nbytes / 1e9:.2f} GB")

cell_offset = 0
for folder, stage, day in SAMPLES:
    t_start     = _time.time()
    sample_path = DATA_DIR / folder

    row, col, data, n_genes, n_cells, barcodes, features = read_10x_mtx(sample_path)

    qc   = compute_cell_qc(row, col, data, n_genes, n_cells, MT_MASK)
    keep = (
        (qc["n_genes_by_counts"] >= QC_PARAMS["min_genes_per_cell"])
        & (qc["n_genes_by_counts"] < QC_PARAMS["max_genes_per_cell"])
        & (qc["total_counts"] >= QC_PARAMS["min_counts_per_cell"])
        & (qc["pct_counts_mt"] < QC_PARAMS["max_pct_mito"])
    ).values

    kept_cell_set = set(np.where(keep)[0])
    n_kept        = len(kept_cell_set)

    old_to_new = {c_old: idx for idx, c_old in enumerate(sorted(kept_cell_set))}

    for i in range(len(row)):
        g = row[i]
        c = col[i]
        if g in hvg_idx_set and c in kept_cell_set:
            X_hvg[cell_offset + old_to_new[c], hvg_col_map[g]] += data[i]

    cell_offset += n_kept

    elapsed = _time.time() - t_start
    mem_gb  = psutil.Process().memory_info().rss / 1e9
    print(f"  {folder}: {n_kept:>6,} cells loaded  [{elapsed:.1f}s, {mem_gb:.1f}GB RSS]")

    del row, col, data, barcodes, features, qc, keep, kept_cell_set, old_to_new
    gc.collect()

assert cell_offset == total_cells_passed, \
    f"Cell count mismatch: {cell_offset} vs {total_cells_passed}"

print(f"\nPass 2 complete: HVG matrix shape = {X_hvg.shape}")

# =============================================================================
# 6. NORMALIZATION + LOG TRANSFORM
# =============================================================================
print("\n--- Normalization ---")

cell_totals  = obs["total_counts"].values.astype(np.float64)
scale_factors = NORM_TARGET_SUM / cell_totals
X_norm        = X_hvg * scale_factors[:, np.newaxis]
X_log         = np.log1p(X_norm).astype(np.float32)

del X_hvg, X_norm
gc.collect()

print(f"  Normalized to {NORM_TARGET_SUM:.0f} counts/cell + log1p")
print(f"  X_log shape: {X_log.shape}, dtype: {X_log.dtype}")

# =============================================================================
# 7. PCA
# =============================================================================
print("\n--- PCA ---")

gene_means  = X_log.mean(axis=0)
X_centered  = X_log - gene_means

print(f"  Computing covariance matrix ({N_HVG} x {N_HVG})...")
cov = (X_centered.T @ X_centered) / (X_centered.shape[0] - 1)

print(f"  Eigendecomposition...")
eigenvalues, eigenvectors = np.linalg.eigh(cov)
eigenvalues  = eigenvalues[::-1]
eigenvectors = eigenvectors[:, ::-1]
eigenvalues  = eigenvalues[:N_PCS]
eigenvectors = eigenvectors[:, :N_PCS]

var_ratio = eigenvalues / cov.trace()
X_pca     = X_centered @ eigenvectors   # (n_cells x N_PCS)

print(f"  PCA: {N_PCS} components computed")
print(f"  Variance explained (top 10): {var_ratio[:10].round(4)}")
print(f"  Cumulative (top 30): {var_ratio[:30].sum():.3f}")

del X_centered, cov
gc.collect()

# =============================================================================
# 8. SAVE ARCHIVE OUTPUTS (obs CSV, var CSV, PCA npz)
# =============================================================================
print("\n--- Saving archive outputs ---")

obs_path = PROCESSED_DIR / f"{TIMESTAMP}_GSE230659_qc_obs.csv"
obs.to_csv(obs_path)
print(f"  Cell metadata (archive) → {obs_path}")

pca_path = PROCESSED_DIR / f"{TIMESTAMP}_GSE230659_qc_pca.npz"
np.savez_compressed(
    pca_path,
    X_pca=X_pca.astype(np.float32),
    components=eigenvectors.T.astype(np.float32),
    eigenvalues=eigenvalues.astype(np.float32),
    variance_ratio=var_ratio.astype(np.float32),
    gene_means=gene_means.astype(np.float32),
)
print(f"  PCA archive → {pca_path}")

# =============================================================================
# 9. CONSTRUCT AND WRITE ANNDATA (canonical benchmark input)
# =============================================================================
print("\n" + "=" * 70)
print("CONSTRUCTING AnnData: adata_benchmark.h5ad")
print("=" * 70)

# --- Build obs for AnnData ---
# Required columns per framework v2 §5 and this script's requirements:
#   cell_id, dataset_id, sample_id, stage, day_within_stage, abs_day,
#   time_label, scTimeBench_timepoint, scTimeBench_cell_type
# All of these were added in Pass 1.

obs_adata = obs[[
    "cell_id",
    "dataset_id",
    "sample_id",
    "stage",
    "stage_day_label",
    "day_within_stage",
    "abs_day",
    "time_label",
    "scTimeBench_timepoint",
    "scTimeBench_cell_type",
    # QC metrics retained for downstream reference
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
    "pct_counts_ribo",
]].copy()
obs_adata.index = obs["barcode"]     # obs_names = unique barcode strings
obs_adata.index.name = "obs_names"

# --- Build var for AnnData (HVG subset only) ---
# Required columns per framework v2: gene_symbol, gene_id, highly_variable, is_hvg
var_hvg = var.loc[var["highly_variable"], [
    "gene_symbol",
    "gene_id",
    "gene_name",
    "highly_variable",
    "is_hvg",
    "mt",
    "ribo",
    "n_cells",
    "mean_counts",
    "var_counts",
]].copy()
var_hvg.index = GENE_NAMES[hvg_gene_indices]
var_hvg.index.name = "var_names"

# --- Construct AnnData ---
# X is the log-normalized HVG matrix, stored as a sparse CSR matrix.
adata = ad.AnnData(
    X=sp.csr_matrix(X_log),
    obs=obs_adata,
    var=var_hvg,
)

# PCA embedding
adata.obsm["X_pca"] = X_pca.astype(np.float32)

# PCA metadata
adata.uns["pca"] = {
    "variance_ratio": var_ratio.astype(np.float32),
    "variance":       eigenvalues.astype(np.float32),
}

# --- Required uns fields (framework v2 §5, §12.4) ---
adata.uns["dataset_id"]  = DATASET_ID
adata.uns["time_axis"]   = "observed_time"

adata.uns["preprocessing_summary"] = {
    "script":              "scripts/01_load_and_qc.py",
    "version":             "v4",
    "framework":           "scTimeBench v2",
    "dataset":             DATASET_ID,
    "timestamp":           TIMESTAMP,
    "n_samples":           len(SAMPLES),
    "qc_params":           QC_PARAMS,
    "norm_target_sum":     NORM_TARGET_SUM,
    "n_hvg":               N_HVG,
    "n_pcs":               N_PCS,
    "cell_type_annotation": (
        "provisional_stage_proxy"
        if CELLTYPE_ANNOTATION_IS_PROXY
        else "external_annotation"
    ),
    "cell_type_proxy_note": (
        "scTimeBench_cell_type is a coarse stage-level proxy mapping: "
        "StageI→somatic, StageII→early_transition, StageIII→intermediate, "
        "hCiPSCs→pluripotent. "
        "Replace with a frozen per-cell annotation when one is available."
        if CELLTYPE_ANNOTATION_IS_PROXY
        else ""
    ),
}

print(f"  AnnData shape: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"  obs columns  : {list(adata.obs.columns)}")
print(f"  var columns  : {list(adata.var.columns)}")
print(f"  obsm keys    : {list(adata.obsm.keys())}")
print(f"  uns keys     : {list(adata.uns.keys())}")

# --- Write canonical h5ad ---
print(f"\n  Writing canonical benchmark input → {ADATA_BENCHMARK_PATH}")
adata.write_h5ad(ADATA_BENCHMARK_PATH)
print(f"  Done. File size: {ADATA_BENCHMARK_PATH.stat().st_size / 1e6:.1f} MB")

# =============================================================================
# 10. QC VISUALIZATIONS (secondary outputs)
# =============================================================================
print("\n--- Generating QC plots ---")

STAGE_COLORS = {
    "StageI": "#2196F3", "StageII": "#FF9800",
    "StageIII": "#4CAF50", "hCiPSCs": "#E91E63",
}
stage_order = ["StageI", "StageII", "StageIII", "hCiPSCs"]

# Post-QC violins
fig, axes = plt.subplots(1, 4, figsize=(20, 5))
for ax, metric, label in zip(
    axes,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
    ["Genes per cell", "UMI counts", "Mito %", "Ribo %"],
):
    sns.violinplot(
        data=obs, x="stage", y=metric, order=stage_order,
        palette=STAGE_COLORS, ax=ax, inner="box", cut=0, linewidth=0.5,
    )
    ax.set_title(label, fontsize=12)
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=30)

fig.suptitle("Post-QC Distributions by Stage", fontsize=14, y=1.02)
fig.tight_layout()
violin_path = FIGURES_DIR / f"{TIMESTAMP}_qc_violin_post.png"
fig.savefig(violin_path, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {violin_path}")

# Scatter: genes vs counts, mito%
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
sc1 = axes[0].scatter(
    obs["total_counts"], obs["n_genes_by_counts"],
    c=obs["pct_counts_mt"], cmap="RdYlBu_r", s=1, alpha=0.3, rasterized=True,
)
axes[0].set_xlabel("Total counts")
axes[0].set_ylabel("Genes detected")
axes[0].set_title("Genes vs Counts (color=mito%)")
plt.colorbar(sc1, ax=axes[0], label="Mito %")

for stg in stage_order:
    mask = obs["stage"] == stg
    axes[1].scatter(
        obs.loc[mask, "total_counts"], obs.loc[mask, "pct_counts_mt"],
        c=STAGE_COLORS[stg], s=1, alpha=0.3, label=stg, rasterized=True,
    )
axes[1].set_xlabel("Total counts")
axes[1].set_ylabel("Mito %")
axes[1].set_title("Mito% vs Counts (color=stage)")
axes[1].legend(markerscale=5, fontsize=9)
fig.tight_layout()
scatter_path = FIGURES_DIR / f"{TIMESTAMP}_qc_scatter.png"
fig.savefig(scatter_path, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {scatter_path}")

# Cell counts per timepoint
time_counts = (
    obs.groupby(["abs_day", "stage_day_label", "stage"], observed=True)
    .size()
    .reset_index(name="n_cells")
    .sort_values("abs_day")
)

fig, ax = plt.subplots(figsize=(14, 5))
bar_colors = [STAGE_COLORS.get(s, "#999") for s in time_counts["stage"]]
ax.bar(range(len(time_counts)), time_counts["n_cells"], color=bar_colors)
ax.set_xticks(range(len(time_counts)))
ax.set_xticklabels(time_counts["stage_day_label"], rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Number of cells (post-QC)")
ax.set_title("Cell counts per timepoint after QC")
for i, (_, r) in enumerate(time_counts.iterrows()):
    ax.text(i, r["n_cells"] + 50, f"{r['n_cells']:,}", ha="center", va="bottom", fontsize=8)
fig.tight_layout()
cellcount_path = FIGURES_DIR / f"{TIMESTAMP}_cell_counts_per_time.png"
fig.savefig(cellcount_path, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {cellcount_path}")

# PCA plot
fig, axes = plt.subplots(1, 3, figsize=(21, 6))
sc1 = axes[0].scatter(
    X_pca[:, 0], X_pca[:, 1], c=obs["abs_day"].values,
    cmap="viridis", s=1, alpha=0.3, rasterized=True,
)
axes[0].set_xlabel(f"PC1 ({var_ratio[0]:.1%})")
axes[0].set_ylabel(f"PC2 ({var_ratio[1]:.1%})")
axes[0].set_title("PCA — Absolute Day")
plt.colorbar(sc1, ax=axes[0], label="Absolute Day")

for stg in stage_order:
    mask = obs["stage"].values == stg
    axes[1].scatter(X_pca[mask, 0], X_pca[mask, 1],
                    c=STAGE_COLORS[stg], s=1, alpha=0.3, label=stg, rasterized=True)
axes[1].set_xlabel(f"PC1 ({var_ratio[0]:.1%})")
axes[1].set_ylabel(f"PC2 ({var_ratio[1]:.1%})")
axes[1].set_title("PCA — Stage")
axes[1].legend(markerscale=5, fontsize=9)

unique_days = sorted(obs["abs_day"].unique())
day_colors  = plt.cm.turbo(np.linspace(0, 1, len(unique_days)))
day_cmap    = dict(zip(unique_days, day_colors))
for d in unique_days:
    mask = obs["abs_day"].values == d
    lbl  = obs.loc[mask, "stage_day_label"].iloc[0] if mask.sum() > 0 else ""
    axes[2].scatter(X_pca[mask, 0], X_pca[mask, 1],
                    c=[day_cmap[d]], s=1, alpha=0.3, label=lbl, rasterized=True)
axes[2].set_xlabel(f"PC1 ({var_ratio[0]:.1%})")
axes[2].set_ylabel(f"PC2 ({var_ratio[1]:.1%})")
axes[2].set_title("PCA — Stage+Day Label")
axes[2].legend(markerscale=5, fontsize=7, loc="center left", bbox_to_anchor=(1.0, 0.5))
fig.tight_layout()
pca_fig_path = FIGURES_DIR / f"{TIMESTAMP}_pca_time.png"
fig.savefig(pca_fig_path, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {pca_fig_path}")

# =============================================================================
# 11. QC SUMMARY TABLE
# =============================================================================
qc_summary = (
    obs.groupby(["stage", "stage_day_label", "abs_day"], observed=True)
    .agg(
        n_cells=("n_genes_by_counts", "size"),
        median_genes=("n_genes_by_counts", "median"),
        median_counts=("total_counts", "median"),
        median_mito_pct=("pct_counts_mt", "median"),
        median_ribo_pct=("pct_counts_ribo", "median"),
    )
    .reset_index()
    .sort_values("abs_day")
)

csv_path = FIGURES_DIR / f"{TIMESTAMP}_qc_summary.csv"
qc_summary.to_csv(csv_path, index=False)
print(f"\n  QC summary → {csv_path}")
print("\n--- QC Summary ---")
print(qc_summary.to_string(index=False))

# =============================================================================
# 12. VALIDATION BLOCK
# =============================================================================
print("\n" + "=" * 70)
print("VALIDATION: Reloading adata_benchmark.h5ad")
print("=" * 70)

adata_check = ad.read_h5ad(ADATA_BENCHMARK_PATH)

REQUIRED_OBS = [
    "cell_id", "dataset_id", "sample_id", "stage",
    "day_within_stage", "abs_day", "time_label",
    "scTimeBench_timepoint", "scTimeBench_cell_type",
]

print(f"  n_obs (cells)  : {adata_check.n_obs:,}")
print(f"  n_vars (genes) : {adata_check.n_vars:,}")
print(f"  X dtype        : {adata_check.X.dtype}")
print(f"  X is sparse    : {sp.issparse(adata_check.X)}")

print(f"\n  Required obs columns:")
all_obs_ok = True
for col in REQUIRED_OBS:
    present = col in adata_check.obs.columns
    status  = "✓" if present else "✗ MISSING"
    print(f"    {status}  {col}")
    if not present:
        all_obs_ok = False

print(f"\n  X_pca in obsm  : {'✓' if 'X_pca' in adata_check.obsm else '✗ MISSING'}")
print(f"  X_pca shape    : {adata_check.obsm['X_pca'].shape if 'X_pca' in adata_check.obsm else 'n/a'}")

print(f"\n  uns keys       : {list(adata_check.uns.keys())}")
print(f"  dataset_id     : {adata_check.uns.get('dataset_id', 'n/a')}")
print(f"  time_axis      : {adata_check.uns.get('time_axis', 'n/a')}")

cell_type_note = adata_check.uns.get("preprocessing_summary", {}).get("cell_type_annotation", "n/a")
print(f"  cell_type_note : {cell_type_note}")

validation_ok = all_obs_ok and "X_pca" in adata_check.obsm
print(f"\n  Validation {'PASSED ✓' if validation_ok else 'FAILED ✗ — check missing fields above'}")

# =============================================================================
# 13. FINAL REPORT
# =============================================================================
elapsed_total = _time.time() - T0
mem_gb        = psutil.Process().memory_info().rss / 1e9

print("\n" + "=" * 70)
print(f"[{TIMESTAMP}] Step 01 COMPLETE")
print(f"  Primary output      : {ADATA_BENCHMARK_PATH}")
print(f"  Total cells post-QC : {total_cells_passed:,}")
print(f"  HVGs selected       : {hvg_mask.sum()}")
print(f"  PCA components      : {N_PCS}")
print(f"  Total runtime       : {elapsed_total:.0f}s ({elapsed_total/60:.1f}min)")
print(f"  Peak RSS memory     : {mem_gb:.2f} GB")
print("=" * 70)
