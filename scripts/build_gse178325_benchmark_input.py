#!/usr/bin/env python3
"""
Build the GSE178325 human benchmark input with the same processing strategy as
scripts/01_load_and_qc.py for GSE230659.

The first benchmark-ready input uses the clean 0618-batch time-resolved trajectory:
StageI day 0/0.5/1/2/4, StageII day 4/8/12, StageIII day 4/8,
StageIV day 1/2/4/10, and hCiPSCs-0618 as the terminal sample.

Outputs
-------
Primary benchmark input:
    benchmark/inputs/gse178325_human_hvg2000/
        GSE178325_human_HVG2000_benchmark_input.h5ad

Archive outputs:
    data/processed/gse178325_human/{timestamp}_GSE178325_qc_obs.csv
    data/processed/gse178325_human/{timestamp}_GSE178325_qc_var.csv
    data/processed/gse178325_human/{timestamp}_GSE178325_qc_pca.npz
    benchmark/reports/qc/gse178325_human/{timestamp}_*.png/csv
"""

from __future__ import annotations

import gc
import gzip
import os
import time as _time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
import matplotlib.pyplot as plt  # noqa: E402

try:
    import psutil
except ImportError:
    psutil = None


TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")
T0 = _time.time()


def _find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
        raise FileNotFoundError(
            f"TRAJ_PROJECT_ROOT={env!r} does not exist. Correct it and retry."
        )

    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "data").exists():
            return candidate
    return here.parent


PROJECT_ROOT = _find_project_root()
DATA_DIR = PROJECT_ROOT / "data" / "gse178325_human" / "rna_seq_10x"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed" / "gse178325_human"
INPUT_DIR = PROJECT_ROOT / "benchmark" / "inputs" / "gse178325_human_hvg2000"
FIGURES_DIR = PROJECT_ROOT / "benchmark" / "reports" / "qc" / "gse178325_human"
ADATA_BENCHMARK_PATH = INPUT_DIR / "GSE178325_human_HVG2000_benchmark_input.h5ad"

for d in [PROCESSED_DIR, INPUT_DIR, FIGURES_DIR]:
    d.mkdir(parents=True, exist_ok=True)


QC_PARAMS = {
    "min_genes_per_cell": 200,
    "max_genes_per_cell": 8000,
    "min_cells_per_gene": 3,
    "max_pct_mito": 20.0,
    "min_counts_per_cell": 500,
}

NORM_TARGET_SUM = 1e4
N_HVG = 2000
N_PCS = 50
DATASET_ID = "GSE178325"
BATCH_ID = "0618"

# GSE178325 has staged day labels where each stage's day count restarts from
# that stage's beginning. Convert them to a cumulative observed-time axis by
# adding the elapsed duration of previous stages:
#   StageI ends at day 4
#   StageII starts after StageI day 4 and spans to StageII day 12
#   StageIII starts after StageI day 4 + StageII day 12
#   StageIV starts after StageI day 4 + StageII day 12 + StageIII day 8
STAGE_OFFSETS = {
    "StageI": 0.0,
    "StageII": 4.0,
    "StageIII": 16.0,
    "StageIV": 24.0,
}
HCIPSC_ABSOLUTE_DAY = 36.0

# Main 0618-batch time-resolved trajectory. Other samples in the GEO series
# include untreated starting populations, perturbations, ESC, and samples from
# other date-coded batches such as 0330, 0605, 0809, 1117, and 1230; they are
# intentionally excluded from this first benchmark input.
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

for _gsm, _stage, _day, _sample_name in SAMPLES:
    if BATCH_ID not in _sample_name:
        raise ValueError(
            f"Non-{BATCH_ID} sample found in GSE178325 manifest: {_sample_name}"
        )


def available_ram_gb() -> float:
    if psutil is None:
        return float("nan")
    return psutil.virtual_memory().available / 1e9


def rss_gb() -> float:
    if psutil is None:
        return float("nan")
    return psutil.Process().memory_info().rss / 1e9

STAGE_TO_CELLTYPE_PROXY: Dict[str, str] = {
    "StageI": "somatic_or_initiating",
    "StageII": "early_transition",
    "StageIII": "intermediate",
    "StageIV": "late_transition",
    "hCiPSCs": "pluripotent",
}
CELLTYPE_ANNOTATION_IS_PROXY = True

STAGE_COLORS = {
    "StageI": "#2C7FB8",
    "StageII": "#F28E2B",
    "StageIII": "#59A14F",
    "StageIV": "#B07AA1",
    "hCiPSCs": "#E15759",
}
STAGE_ORDER = ["StageI", "StageII", "StageIII", "StageIV", "hCiPSCs"]


def _sample_path(gsm: str) -> Path:
    return DATA_DIR / gsm


def read_10x_mtx(
    sample_dir: Path,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int, int, pd.DataFrame, pd.DataFrame]:
    """Read 10x matrix.mtx.gz, barcodes.tsv.gz, and features.tsv.gz."""
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
        parts = line.strip().split()
        n_genes, n_cells, _n_entries = int(parts[0]), int(parts[1]), int(parts[2])

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
    return row, col, data, n_genes, n_cells, barcodes, features


def compute_cell_qc(
    row: np.ndarray,
    col: np.ndarray,
    data: np.ndarray,
    n_cells: int,
    mt_mask: np.ndarray,
) -> pd.DataFrame:
    total_counts = np.bincount(col, weights=data.astype(np.float64), minlength=n_cells)
    n_genes_detected = np.bincount(col, minlength=n_cells)
    is_mito = mt_mask[row]
    mito_counts = np.bincount(
        col[is_mito],
        weights=data[is_mito].astype(np.float64),
        minlength=n_cells,
    )
    pct_mito = np.where(total_counts > 0, 100.0 * mito_counts / total_counts, 0.0)
    return pd.DataFrame(
        {
            "n_genes_by_counts": n_genes_detected.astype(np.int32),
            "total_counts": total_counts.astype(np.float32),
            "pct_counts_mt": pct_mito.astype(np.float32),
        }
    )


def abs_day_for(stage: str, day: Optional[float]) -> float:
    if stage == "hCiPSCs":
        return HCIPSC_ABSOLUTE_DAY
    if day is None:
        raise ValueError(f"Non-terminal sample for {stage} is missing day.")
    return STAGE_OFFSETS[stage] + float(day)


def stage_day_label(stage: str, day: Optional[float]) -> str:
    if stage == "hCiPSCs":
        return "hCiPSCs"
    return f"{stage}_Day{day:g}"


print(f"[{TIMESTAMP}] GSE178325 benchmark input builder")
print(f"  Project root : {PROJECT_ROOT}")
print(f"  Data dir     : {DATA_DIR}")
print(f"  Output h5ad  : {ADATA_BENCHMARK_PATH}")
print(f"  RAM available: {available_ram_gb():.1f} GB")


print("\n" + "=" * 70)
print("PASS 1: Per-sample QC filtering + gene statistics")
print("=" * 70)

ref_features = pd.read_csv(
    _sample_path(SAMPLES[0][0]) / "features.tsv.gz",
    sep="\t",
    header=None,
    names=["gene_id", "gene_name", "feature_type"],
)
GENE_NAMES = ref_features["gene_name"].values
GENE_IDS = ref_features["gene_id"].values
N_GENES_REF = len(GENE_NAMES)
MT_MASK = np.array([str(g).startswith("MT-") for g in GENE_NAMES], dtype=bool)
RIBO_MASK = np.array(
    [str(g).startswith("RPS") or str(g).startswith("RPL") for g in GENE_NAMES],
    dtype=bool,
)

print(f"  Reference genes: {N_GENES_REF:,}")
print(f"  MT genes        : {MT_MASK.sum()}")
print(f"  Ribosomal genes : {RIBO_MASK.sum()}")

gene_sum = np.zeros(N_GENES_REF, dtype=np.float64)
gene_sq_sum = np.zeros(N_GENES_REF, dtype=np.float64)
gene_n_cells = np.zeros(N_GENES_REF, dtype=np.int64)
total_cells_passed = 0
sample_meta_list: List[pd.DataFrame] = []
sample_cell_counts: Dict[str, Dict[str, object]] = {}

for gsm, stage, day, sample_name in SAMPLES:
    t_start = _time.time()
    sample_path = _sample_path(gsm)
    row, col, data, n_genes, n_cells, barcodes, features = read_10x_mtx(sample_path)
    if n_genes != N_GENES_REF:
        raise AssertionError(f"Gene count mismatch in {gsm}: {n_genes} vs {N_GENES_REF}")

    qc = compute_cell_qc(row, col, data, n_cells, MT_MASK)
    keep = (
        (qc["n_genes_by_counts"] >= QC_PARAMS["min_genes_per_cell"])
        & (qc["n_genes_by_counts"] < QC_PARAMS["max_genes_per_cell"])
        & (qc["total_counts"] >= QC_PARAMS["min_counts_per_cell"])
        & (qc["pct_counts_mt"] < QC_PARAMS["max_pct_mito"])
    ).values
    n_kept = int(keep.sum())
    n_removed = int(n_cells - n_kept)

    is_ribo = RIBO_MASK[row]
    ribo_counts_all = np.bincount(
        col[is_ribo], weights=data[is_ribo].astype(np.float64), minlength=n_cells
    )
    total_counts_all = np.bincount(
        col, weights=data.astype(np.float64), minlength=n_cells
    )
    pct_ribo = np.where(
        total_counts_all > 0, 100.0 * ribo_counts_all / total_counts_all, 0.0
    )

    abs_day = abs_day_for(stage, day)
    label = stage_day_label(stage, day)
    barcodes_kept = barcodes["barcode"].values[keep]
    cell_ids = [f"{gsm}_{b}" for b in barcodes_kept]

    cell_meta = pd.DataFrame(
        {
            "barcode": cell_ids,
            "cell_id": cell_ids,
            "dataset_id": DATASET_ID,
            "batch_id": BATCH_ID,
            "sample_id": gsm,
            "sample_name": sample_name,
            "stage": stage,
            "stage_day_label": label,
            "day_within_stage": day if day is not None else np.nan,
            "abs_day": abs_day,
            "time_label": abs_day,
            "scTimeBench_timepoint": abs_day,
            "scTimeBench_cell_type": STAGE_TO_CELLTYPE_PROXY.get(stage, "unknown"),
            "n_genes_by_counts": qc["n_genes_by_counts"].values[keep],
            "total_counts": qc["total_counts"].values[keep],
            "pct_counts_mt": qc["pct_counts_mt"].values[keep],
            "pct_counts_ribo": pct_ribo[keep].astype(np.float32),
        }
    )
    sample_meta_list.append(cell_meta)

    entry_mask = keep[col]
    kept_row = row[entry_mask]
    kept_data = data[entry_mask].astype(np.float64)
    gene_sum += np.bincount(kept_row, weights=kept_data, minlength=N_GENES_REF)
    gene_sq_sum += np.bincount(kept_row, weights=kept_data**2, minlength=N_GENES_REF)
    gene_n_cells += np.bincount(kept_row, minlength=N_GENES_REF)

    total_cells_passed += n_kept
    sample_cell_counts[gsm] = {
        "sample_name": sample_name,
        "stage": stage,
        "abs_day": abs_day,
        "cells_before": n_cells,
        "cells_after": n_kept,
    }

    elapsed = _time.time() - t_start
    mem_gb = rss_gb()
    print(
        f"  {gsm}: {n_cells:>6,} -> {n_kept:>6,} cells "
        f"(dropped {n_removed:>4})  abs_day={abs_day:<5.1f} "
        f"[{elapsed:.1f}s, {mem_gb:.1f}GB RSS]"
    )

    del row, col, data, barcodes, features, qc, keep, kept_row, kept_data
    del entry_mask, ribo_counts_all, total_counts_all, pct_ribo
    gc.collect()

obs = pd.concat(sample_meta_list, ignore_index=True)
obs.index = obs["barcode"]
del sample_meta_list
gc.collect()

print(
    f"\nPass 1 complete: {total_cells_passed:,} cells passed QC "
    f"across {len(SAMPLES)} samples"
)


print("\n" + "=" * 70)
print("GENE FILTERING + HVG SELECTION")
print("=" * 70)

gene_passes = gene_n_cells >= QC_PARAMS["min_cells_per_gene"]
gene_mean = gene_sum / total_cells_passed
gene_var = (gene_sq_sum / total_cells_passed) - gene_mean**2
gene_var = np.maximum(gene_var, 0)

valid = gene_passes & (gene_mean > 0)
log_mean = np.full(N_GENES_REF, np.nan)
log_var = np.full(N_GENES_REF, np.nan)
log_mean[valid] = np.log10(gene_mean[valid] + 1e-10)
log_var[valid] = np.log10(gene_var[valid] + 1e-10)

n_bins = 20
valid_idx = np.where(valid)[0]
lm = log_mean[valid]
lv = log_var[valid]
bin_edges = np.percentile(lm, np.linspace(0, 100, n_bins + 1))
bin_edges[0] -= 1e-10
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
hvg_mask = np.zeros(N_GENES_REF, dtype=bool)
hvg_mask[hvg_indices] = True

print(
    f"  Gene filter (min {QC_PARAMS['min_cells_per_gene']} cells): "
    f"{N_GENES_REF:,} -> {int(gene_passes.sum()):,}"
)
print(f"  HVGs selected: {int(hvg_mask.sum())}")

var = pd.DataFrame(
    {
        "gene_symbol": GENE_NAMES,
        "gene_id": GENE_IDS,
        "highly_variable": hvg_mask,
        "is_hvg": hvg_mask,
        "gene_name": GENE_NAMES,
        "mt": MT_MASK,
        "ribo": RIBO_MASK,
        "n_cells": gene_n_cells.astype(np.int32),
        "mean_counts": gene_mean.astype(np.float32),
        "var_counts": gene_var.astype(np.float32),
        "passes_filter": gene_passes,
    },
    index=GENE_NAMES,
)
var.index.name = "var_names"

var_path = PROCESSED_DIR / f"{TIMESTAMP}_GSE178325_qc_var.csv"
var.to_csv(var_path)
print(f"  Saved gene metadata -> {var_path}")


print("\n" + "=" * 70)
print("PASS 2: Build HVG expression matrix")
print("=" * 70)

hvg_gene_indices = np.where(hvg_mask)[0]
hvg_rank = np.full(N_GENES_REF, -1, dtype=np.int32)
hvg_rank[hvg_gene_indices] = np.arange(len(hvg_gene_indices), dtype=np.int32)

X_hvg = np.zeros((total_cells_passed, N_HVG), dtype=np.float32)
print(
    f"  Allocated HVG matrix: {total_cells_passed:,} x {N_HVG} = "
    f"{X_hvg.nbytes / 1e9:.2f} GB"
)

cell_offset = 0
for gsm, stage, day, sample_name in SAMPLES:
    t_start = _time.time()
    row, col, data, n_genes, n_cells, barcodes, features = read_10x_mtx(_sample_path(gsm))
    qc = compute_cell_qc(row, col, data, n_cells, MT_MASK)
    keep = (
        (qc["n_genes_by_counts"] >= QC_PARAMS["min_genes_per_cell"])
        & (qc["n_genes_by_counts"] < QC_PARAMS["max_genes_per_cell"])
        & (qc["total_counts"] >= QC_PARAMS["min_counts_per_cell"])
        & (qc["pct_counts_mt"] < QC_PARAMS["max_pct_mito"])
    ).values
    kept_rank = np.full(n_cells, -1, dtype=np.int32)
    kept_cells = np.where(keep)[0]
    kept_rank[kept_cells] = np.arange(len(kept_cells), dtype=np.int32)

    entry_mask = (hvg_rank[row] >= 0) & (kept_rank[col] >= 0)
    target_rows = cell_offset + kept_rank[col[entry_mask]]
    target_cols = hvg_rank[row[entry_mask]]
    np.add.at(X_hvg, (target_rows, target_cols), data[entry_mask].astype(np.float32))

    cell_offset += len(kept_cells)
    elapsed = _time.time() - t_start
    mem_gb = rss_gb()
    print(
        f"  {gsm}: {len(kept_cells):>6,} cells loaded "
        f"[{elapsed:.1f}s, {mem_gb:.1f}GB RSS]"
    )

    del row, col, data, barcodes, features, qc, keep, kept_rank, kept_cells
    del entry_mask, target_rows, target_cols
    gc.collect()

if cell_offset != total_cells_passed:
    raise AssertionError(f"Cell count mismatch: {cell_offset} vs {total_cells_passed}")


print("\n--- Normalization ---")
cell_totals = obs["total_counts"].values.astype(np.float64)
scale_factors = NORM_TARGET_SUM / cell_totals
X_norm = X_hvg * scale_factors[:, np.newaxis]
X_log = np.log1p(X_norm).astype(np.float32)
del X_hvg, X_norm
gc.collect()
print(f"  Normalized to {NORM_TARGET_SUM:.0f} counts/cell + log1p")


print("\n--- PCA ---")
gene_means = X_log.mean(axis=0)
X_centered = X_log - gene_means
cov = (X_centered.T @ X_centered) / (X_centered.shape[0] - 1)
eigenvalues, eigenvectors = np.linalg.eigh(cov)
eigenvalues = eigenvalues[::-1][:N_PCS]
eigenvectors = eigenvectors[:, ::-1][:, :N_PCS]
var_ratio = eigenvalues / cov.trace()
X_pca = X_centered @ eigenvectors
del X_centered, cov
gc.collect()
print(f"  PCA: {N_PCS} components computed")
print(f"  Variance explained (top 10): {var_ratio[:10].round(4)}")


print("\n--- Saving archive outputs ---")
obs_path = PROCESSED_DIR / f"{TIMESTAMP}_GSE178325_qc_obs.csv"
obs.to_csv(obs_path)
print(f"  Cell metadata archive -> {obs_path}")

pca_path = PROCESSED_DIR / f"{TIMESTAMP}_GSE178325_qc_pca.npz"
np.savez_compressed(
    pca_path,
    X_pca=X_pca.astype(np.float32),
    components=eigenvectors.T.astype(np.float32),
    eigenvalues=eigenvalues.astype(np.float32),
    variance_ratio=var_ratio.astype(np.float32),
    gene_means=gene_means.astype(np.float32),
)
print(f"  PCA archive -> {pca_path}")


print("\n" + "=" * 70)
print("CONSTRUCTING AnnData: GSE178325 benchmark input")
print("=" * 70)

obs_adata = obs[
    [
        "cell_id",
        "dataset_id",
        "batch_id",
        "sample_id",
        "sample_name",
        "stage",
        "stage_day_label",
        "day_within_stage",
        "abs_day",
        "time_label",
        "scTimeBench_timepoint",
        "scTimeBench_cell_type",
        "n_genes_by_counts",
        "total_counts",
        "pct_counts_mt",
        "pct_counts_ribo",
    ]
].copy()
obs_adata.index = obs["barcode"]
obs_adata.index.name = "obs_names"

var_hvg = var.loc[
    var["highly_variable"],
    [
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
    ],
].copy()
var_hvg.index = GENE_NAMES[hvg_gene_indices]
var_hvg.index.name = "var_names"

adata = ad.AnnData(X=sp.csr_matrix(X_log), obs=obs_adata, var=var_hvg)
adata.obsm["X_pca"] = X_pca.astype(np.float32)
adata.uns["pca"] = {
    "variance_ratio": var_ratio.astype(np.float32),
    "variance": eigenvalues.astype(np.float32),
}
adata.uns["dataset_id"] = DATASET_ID
adata.uns["time_axis"] = "observed_time"
adata.uns["source_geo_accession"] = "GSE178325"
adata.uns["preprocessing_summary"] = {
    "script": "scripts/build_gse178325_benchmark_input.py",
    "reference_script": "scripts/01_load_and_qc.py",
    "version": "v1",
    "framework": "scTimeBench v2",
    "dataset": DATASET_ID,
    "timestamp": TIMESTAMP,
    "n_samples": len(SAMPLES),
    "batch_id": BATCH_ID,
    "sample_policy": "0618_batch_time_resolved_main_trajectory_only",
    "excluded_sample_note": (
        "Samples from date-coded batches other than 0618 are excluded. Untimed "
        "stage aggregates, perturbation samples, ESC, starting population "
        "controls, and extra hCiPSC batches are not included in this first "
        "benchmark input."
    ),
    "stage_offsets": STAGE_OFFSETS,
    "hcipsc_absolute_day": HCIPSC_ABSOLUTE_DAY,
    "qc_params": QC_PARAMS,
    "norm_target_sum": NORM_TARGET_SUM,
    "n_hvg": N_HVG,
    "n_pcs": N_PCS,
    "cell_type_annotation": (
        "provisional_stage_proxy"
        if CELLTYPE_ANNOTATION_IS_PROXY
        else "external_annotation"
    ),
}

print(f"  AnnData shape: {adata.n_obs} cells x {adata.n_vars} genes")
print(f"  obs columns  : {list(adata.obs.columns)}")
print(f"  var columns  : {list(adata.var.columns)}")
print(f"  obsm keys    : {list(adata.obsm.keys())}")

print(f"\n  Writing benchmark input -> {ADATA_BENCHMARK_PATH}")
adata.write_h5ad(ADATA_BENCHMARK_PATH)
print(f"  Done. File size: {ADATA_BENCHMARK_PATH.stat().st_size / 1e6:.1f} MB")


print("\n--- Generating QC plots ---")

fig, axes = plt.subplots(1, 4, figsize=(20, 5))
for ax, metric, label in zip(
    axes,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
    ["Genes per cell", "UMI counts", "Mito %", "Ribo %"],
):
    sns.violinplot(
        data=obs,
        x="stage",
        y=metric,
        order=STAGE_ORDER,
        palette=STAGE_COLORS,
        ax=ax,
        inner="box",
        cut=0,
        linewidth=0.5,
    )
    ax.set_title(label, fontsize=12)
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=30)
fig.suptitle("GSE178325 Post-QC Distributions by Stage", fontsize=14, y=1.02)
fig.tight_layout()
violin_path = FIGURES_DIR / f"{TIMESTAMP}_qc_violin_post.png"
fig.savefig(violin_path, dpi=150, bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(14, 5))
sc1 = axes[0].scatter(
    obs["total_counts"],
    obs["n_genes_by_counts"],
    c=obs["pct_counts_mt"],
    cmap="RdYlBu_r",
    s=1,
    alpha=0.3,
    rasterized=True,
)
axes[0].set_xlabel("Total counts")
axes[0].set_ylabel("Genes detected")
axes[0].set_title("Genes vs Counts (color=mito%)")
plt.colorbar(sc1, ax=axes[0], label="Mito %")
for stg in STAGE_ORDER:
    mask = obs["stage"] == stg
    axes[1].scatter(
        obs.loc[mask, "total_counts"],
        obs.loc[mask, "pct_counts_mt"],
        c=STAGE_COLORS[stg],
        s=1,
        alpha=0.3,
        label=stg,
        rasterized=True,
    )
axes[1].set_xlabel("Total counts")
axes[1].set_ylabel("Mito %")
axes[1].set_title("Mito% vs Counts (color=stage)")
axes[1].legend(markerscale=5, fontsize=9)
fig.tight_layout()
scatter_path = FIGURES_DIR / f"{TIMESTAMP}_qc_scatter.png"
fig.savefig(scatter_path, dpi=150, bbox_inches="tight")
plt.close(fig)

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
ax.set_title("GSE178325 cell counts per timepoint after QC")
for i, (_, r) in enumerate(time_counts.iterrows()):
    ax.text(i, r["n_cells"] + 50, f"{r['n_cells']:,}", ha="center", va="bottom", fontsize=8)
fig.tight_layout()
cellcount_path = FIGURES_DIR / f"{TIMESTAMP}_cell_counts_per_time.png"
fig.savefig(cellcount_path, dpi=150, bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(21, 6))
sc1 = axes[0].scatter(
    X_pca[:, 0],
    X_pca[:, 1],
    c=obs["abs_day"].values,
    cmap="viridis",
    s=1,
    alpha=0.3,
    rasterized=True,
)
axes[0].set_xlabel(f"PC1 ({var_ratio[0]:.1%})")
axes[0].set_ylabel(f"PC2 ({var_ratio[1]:.1%})")
axes[0].set_title("PCA - Absolute Day")
plt.colorbar(sc1, ax=axes[0], label="Absolute Day")
for stg in STAGE_ORDER:
    mask = obs["stage"].values == stg
    axes[1].scatter(
        X_pca[mask, 0],
        X_pca[mask, 1],
        c=STAGE_COLORS[stg],
        s=1,
        alpha=0.3,
        label=stg,
        rasterized=True,
    )
axes[1].set_xlabel(f"PC1 ({var_ratio[0]:.1%})")
axes[1].set_ylabel(f"PC2 ({var_ratio[1]:.1%})")
axes[1].set_title("PCA - Stage")
axes[1].legend(markerscale=5, fontsize=9)

unique_days = sorted(obs["abs_day"].unique())
day_colors = plt.cm.turbo(np.linspace(0, 1, len(unique_days)))
for d, color in zip(unique_days, day_colors):
    mask = obs["abs_day"].values == d
    lbl = obs.loc[mask, "stage_day_label"].iloc[0] if mask.sum() > 0 else ""
    axes[2].scatter(
        X_pca[mask, 0],
        X_pca[mask, 1],
        c=[color],
        s=1,
        alpha=0.3,
        label=lbl,
        rasterized=True,
    )
axes[2].set_xlabel(f"PC1 ({var_ratio[0]:.1%})")
axes[2].set_ylabel(f"PC2 ({var_ratio[1]:.1%})")
axes[2].set_title("PCA - Stage+Day Label")
axes[2].legend(markerscale=5, fontsize=7, loc="center left", bbox_to_anchor=(1.0, 0.5))
fig.tight_layout()
pca_fig_path = FIGURES_DIR / f"{TIMESTAMP}_pca_time.png"
fig.savefig(pca_fig_path, dpi=150, bbox_inches="tight")
plt.close(fig)

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
print(f"  Saved QC plots and summary under {FIGURES_DIR}")


print("\n" + "=" * 70)
print("VALIDATION: Reloading GSE178325 benchmark input")
print("=" * 70)
adata_check = ad.read_h5ad(ADATA_BENCHMARK_PATH)
REQUIRED_OBS = [
    "cell_id",
    "dataset_id",
    "sample_id",
    "stage",
    "day_within_stage",
    "abs_day",
    "time_label",
    "scTimeBench_timepoint",
    "scTimeBench_cell_type",
]
all_obs_ok = all(col in adata_check.obs.columns for col in REQUIRED_OBS)
print(f"  n_obs (cells)  : {adata_check.n_obs:,}")
print(f"  n_vars (genes) : {adata_check.n_vars:,}")
print(f"  X dtype        : {adata_check.X.dtype}")
print(f"  X is sparse    : {sp.issparse(adata_check.X)}")
print(f"  X_pca shape    : {adata_check.obsm['X_pca'].shape}")
print(f"  Required obs   : {'OK' if all_obs_ok else 'MISSING'}")
print(f"  dataset_id     : {adata_check.uns.get('dataset_id', 'n/a')}")
validation_ok = all_obs_ok and "X_pca" in adata_check.obsm
print(f"  Validation     : {'PASSED' if validation_ok else 'FAILED'}")

elapsed_total = _time.time() - T0
mem_gb = rss_gb()
print("\n" + "=" * 70)
print(f"[{TIMESTAMP}] GSE178325 builder COMPLETE")
print(f"  Primary output      : {ADATA_BENCHMARK_PATH}")
print(f"  Total cells post-QC : {total_cells_passed:,}")
print(f"  HVGs selected       : {int(hvg_mask.sum())}")
print(f"  PCA components      : {N_PCS}")
print(f"  Total runtime       : {elapsed_total:.0f}s ({elapsed_total / 60:.1f}min)")
print(f"  Peak RSS memory     : {mem_gb:.2f} GB")
print("=" * 70)
