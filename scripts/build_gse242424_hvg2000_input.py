#!/usr/bin/env python3
"""
Build the GSE242424 HVG2000 benchmark input from the validated full-gene h5ad.

Reads  : data/processed/gse242424_human/GSE242424_raw_full_gene_benchmark_input.h5ad
Writes : benchmark/inputs/gse242424_hvg2000/GSE242424_HVG2000_benchmark_input.h5ad
         benchmark/inputs/gse242424_hvg2000/GSE242424_HVG2000_summary.json

HVG method : Seurat-style binned residual variance in log10-mean / log10-var
             space, computed directly on the log-normalized sparse CSR matrix.
             This mirrors the approach in scripts/build_gse178325_benchmark_input.py.

PCA        : 50 components via numpy eigh on the sparse covariance matrix
             (no dense materialisation of the full cells x HVGs array).
             Projection is done in 20k-cell batches.

Usage
-----
    conda activate traj_env
    cd C:\\Users\\37620\\trajectory

    # Full run (requires full-gene h5ad, ~156,969 cells):
    python scripts/build_gse242424_hvg2000_input.py

    # Smoke test with explicit path override:
    python scripts/build_gse242424_hvg2000_input.py --input path/to/full_gene.h5ad

Caveats
-------
* iPSC abs_day=16 is PROVISIONAL.  Revise before using abs_day as a
  formal trajectory time axis.
* HVG scoring is on log1p-normalised expression (not raw counts), which
  is equivalent to the scanpy seurat flavour and consistent with the
  GSE178325 builder convention.

Revision history
    v1  2026-05-14  Initial implementation.
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import sys
import time as _time
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import psutil
except ImportError:
    psutil = None


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

TIMESTAMP  = datetime.now().strftime("%Y%m%d_%H%M")
T0         = _time.time()
DATASET_ID = "GSE242424"
N_HVG      = 2000
N_PCS      = 50
N_HVG_BINS = 20          # bins for residual-variance HVG selection


# ---------------------------------------------------------------------------
# Project root
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

INPUT_H5AD  = (
    PROJECT_ROOT / "data" / "processed" / "gse242424_human"
    / "GSE242424_raw_full_gene_benchmark_input.h5ad"
)
OUTPUT_DIR  = PROJECT_ROOT / "benchmark" / "inputs" / "gse242424_hvg2000"
OUTPUT_H5AD = OUTPUT_DIR / "GSE242424_HVG2000_benchmark_input.h5ad"
OUTPUT_JSON = OUTPUT_DIR / "GSE242424_HVG2000_summary.json"


# ---------------------------------------------------------------------------
# Required obs columns (must be present in the full-gene h5ad)
# ---------------------------------------------------------------------------

REQUIRED_OBS = [
    "cell_id",
    "dataset_id",
    "sample_id",
    "time_label",
    "abs_day",
    "stage",
    "stage_day_label",
    "scTimeBench_timepoint",
    "scTimeBench_cell_type",
    "n_genes_by_counts",
    "total_counts",
    "pct_counts_mt",
]

REQUIRED_VAR = [
    "gene_id",
    "gene_symbol",
    "feature_type",
    "n_cells_by_counts",
    "highly_variable",
    "is_hvg",
]


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------

def rss_gb() -> float:
    if psutil is None:
        return float("nan")
    return psutil.Process().memory_info().rss / 1e9


def _check_obs(adata: ad.AnnData) -> list[str]:
    return [c for c in REQUIRED_OBS if c not in adata.obs.columns]


def _check_var(adata: ad.AnnData) -> list[str]:
    return [c for c in REQUIRED_VAR if c not in adata.var.columns]


# ---------------------------------------------------------------------------
# HVG selection: binned residual variance (Seurat-style)
# Works directly on log-normalized sparse CSR (no raw counts needed).
# ---------------------------------------------------------------------------

def select_hvg(
    X: sp.csr_matrix,
    n_top: int = N_HVG,
    n_bins: int = N_HVG_BINS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Return (hvg_idx, residual_var, log_mean, log_var) for all genes.

    hvg_idx : int array of length n_top, sorted ascending (stable column order).
    Genes with zero mean are assigned residual_var = -inf and never selected.
    """
    X = X.tocsr().astype(np.float64)
    n_cells, n_genes = X.shape

    # Sparse mean and variance
    mean    = np.asarray(X.mean(axis=0)).ravel()               # (n_genes,)
    mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()   # (n_genes,)
    var     = np.maximum(mean_sq - mean ** 2, 0.0)

    # Log10-space for binning (avoids scale issues across magnitude ranges)
    valid    = mean > 0
    log_mean = np.full(n_genes, np.nan, dtype=np.float64)
    log_var  = np.full(n_genes, np.nan, dtype=np.float64)
    log_mean[valid] = np.log10(mean[valid])
    log_var[valid]  = np.log10(var[valid] + 1e-10)

    valid_idx    = np.where(valid)[0]
    lm           = log_mean[valid]
    lv           = log_var[valid]
    bin_edges    = np.percentile(lm, np.linspace(0, 100, n_bins + 1))
    bin_edges[0]  -= 1e-10
    bin_edges[-1] += 1e-10
    bin_assign   = np.clip(np.digitize(lm, bin_edges) - 1, 0, n_bins - 1)

    expected_var = np.zeros(len(lm), dtype=np.float64)
    for b in range(n_bins):
        mask = bin_assign == b
        if mask.sum() > 0:
            expected_var[mask] = np.median(lv[mask])

    residual_var              = np.full(n_genes, -np.inf, dtype=np.float64)
    residual_var[valid_idx]   = lv - expected_var

    hvg_idx = np.sort(np.argsort(residual_var)[::-1][:n_top])
    return hvg_idx, residual_var, log_mean, log_var


# ---------------------------------------------------------------------------
# PCA: sparse covariance + numpy eigh + batched projection
# Peak memory: O(N_HVG^2) for covariance + O(batch x N_HVG) per batch.
# The full (n_cells x N_HVG) dense matrix is never materialised.
# ---------------------------------------------------------------------------

def compute_pca(
    X_hvg: sp.csr_matrix,
    n_components: int = N_PCS,
    batch_size: int = 20_000,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Return (X_pca, components, eigenvalues, variance_ratio, gene_means).

    X_pca            : (n_cells, n_components) float32
    components       : (n_components, N_HVG) float32  -- eigenvectors transposed
    eigenvalues      : (n_components,) float32
    variance_ratio   : (n_components,) float32
    gene_means       : (N_HVG,) float32 -- used to centre before projection
    """
    X_hvg    = X_hvg.tocsr().astype(np.float64)
    n_cells, n_hvg = X_hvg.shape
    n_components   = min(n_components, n_hvg - 1)

    # Sparse mean
    mean_hvg = np.asarray(X_hvg.mean(axis=0)).ravel()  # (n_hvg,)

    # Covariance without materialising dense X:
    # Cov = (X - mean).T @ (X - mean) / (n-1)
    #     = (X.T @ X  -  n * outer(mean, mean)) / (n - 1)
    XtX      = np.asarray((X_hvg.T @ X_hvg).todense())  # (n_hvg, n_hvg)
    cov      = (XtX - n_cells * np.outer(mean_hvg, mean_hvg)) / (n_cells - 1)
    del XtX
    gc.collect()

    print(f"    Covariance ({n_hvg}x{n_hvg}) computed  RSS={rss_gb():.2f}GB")

    eigenvalues, eigenvectors = np.linalg.eigh(cov)          # ascending order
    eigenvalues  = eigenvalues[::-1][:n_components]           # descending
    eigenvectors = eigenvectors[:, ::-1][:, :n_components]    # (n_hvg, n_comp)
    eigenvalues  = np.maximum(eigenvalues, 0.0)
    trace        = max(float(np.trace(cov)), 1e-12)
    var_ratio    = eigenvalues / trace
    del cov
    gc.collect()

    # Batched projection: (X[batch] - mean) @ eigenvectors
    X_pca = np.empty((n_cells, n_components), dtype=np.float32)
    for i in range(0, n_cells, batch_size):
        j    = min(i + batch_size, n_cells)
        Xb   = np.asarray(X_hvg[i:j].todense()).astype(np.float64)
        Xb  -= mean_hvg
        X_pca[i:j] = (Xb @ eigenvectors).astype(np.float32)
    del Xb
    gc.collect()

    print(f"    PCA projection done  RSS={rss_gb():.2f}GB")

    return (
        X_pca,
        eigenvectors.T.astype(np.float32),   # (n_comp, n_hvg)
        eigenvalues.astype(np.float32),
        var_ratio.astype(np.float32),
        mean_hvg.astype(np.float32),
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(input_h5ad: Path = INPUT_H5AD) -> None:
    print(f"[{TIMESTAMP}] GSE242424 HVG2000 benchmark input builder  v1")
    print(f"  Input h5ad   : {input_h5ad}")
    print(f"  Output h5ad  : {OUTPUT_H5AD}")
    print(f"  N_HVG={N_HVG}  N_PCS={N_PCS}  bins={N_HVG_BINS}")
    print()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Load full-gene h5ad
    # ------------------------------------------------------------------
    if not input_h5ad.exists():
        sys.exit(
            f"ERROR: full-gene h5ad not found: {input_h5ad}\n"
            "  Run scripts/build_gse242424_raw_full_gene_input.py first."
        )

    print(f"Loading full-gene h5ad ...")
    adata_full = ad.read_h5ad(input_h5ad)
    print(f"  Shape : {adata_full.shape}")
    print(f"  X sparse : {sp.issparse(adata_full.X)}  dtype : {adata_full.X.dtype}")
    print(f"  RSS : {rss_gb():.2f} GB")

    # ------------------------------------------------------------------
    # 2. Validate required columns
    # ------------------------------------------------------------------
    missing_obs = _check_obs(adata_full)
    missing_var = _check_var(adata_full)
    if missing_obs or missing_var:
        msg = []
        if missing_obs:
            msg.append(f"obs missing: {missing_obs}")
        if missing_var:
            msg.append(f"var missing: {missing_var}")
        sys.exit("ERROR: full-gene h5ad is missing required columns:\n  " + "\n  ".join(msg))

    if adata_full.n_vars < N_HVG:
        sys.exit(
            f"ERROR: full-gene h5ad has only {adata_full.n_vars} genes; "
            f"cannot select {N_HVG} HVGs."
        )

    # ------------------------------------------------------------------
    # 3. HVG selection
    # ------------------------------------------------------------------
    print(f"\nSelecting top {N_HVG} HVGs (binned residual variance) ...")
    t_hvg = _time.time()
    hvg_idx, residual_var, log_mean, log_var = select_hvg(
        adata_full.X, n_top=N_HVG, n_bins=N_HVG_BINS,
    )
    print(
        f"  HVG selection done in {_time.time() - t_hvg:.1f}s  "
        f"RSS={rss_gb():.2f}GB"
    )
    assert len(hvg_idx) == N_HVG, f"Expected {N_HVG} HVGs, got {len(hvg_idx)}"

    hvg_names = adata_full.var.index[hvg_idx].tolist()
    print(f"  Top 5 HVGs: {hvg_names[:5]}")

    # ------------------------------------------------------------------
    # 4. PCA on the HVG subset
    # ------------------------------------------------------------------
    print(f"\nComputing PCA ({N_PCS} components) ...")
    t_pca = _time.time()
    X_hvg_sp = adata_full.X[:, hvg_idx]
    if not sp.issparse(X_hvg_sp):
        X_hvg_sp = sp.csr_matrix(X_hvg_sp)

    X_pca, components, eigenvalues, var_ratio, gene_means = compute_pca(
        X_hvg_sp, n_components=N_PCS,
    )
    print(
        f"  PCA done in {_time.time() - t_pca:.1f}s  "
        f"X_pca shape={X_pca.shape}"
    )
    print(
        f"  Var explained top-10 PCs: "
        f"{[round(float(v), 4) for v in var_ratio[:10]]}"
    )

    # ------------------------------------------------------------------
    # 5. Build obs / var for output
    # ------------------------------------------------------------------
    obs_out = adata_full.obs[REQUIRED_OBS].copy()

    # Build var for the HVG subset
    var_full   = adata_full.var
    var_hvg    = var_full.iloc[hvg_idx][
        ["gene_id", "gene_symbol", "feature_type", "n_cells_by_counts"]
    ].copy()
    var_hvg["highly_variable"]   = True
    var_hvg["is_hvg"]            = True
    var_hvg["hvg_residual_var"]  = residual_var[hvg_idx].astype(np.float32)
    var_hvg["log_mean"]          = log_mean[hvg_idx].astype(np.float32)
    var_hvg["log_var"]           = log_var[hvg_idx].astype(np.float32)
    # Carry over mt/ribo flags if present
    for flag in ("mt", "ribo", "passes_filter"):
        if flag in var_full.columns:
            var_hvg[flag] = var_full.iloc[hvg_idx][flag].values
    var_hvg.index.name = "var_names"

    # ------------------------------------------------------------------
    # 6. Assemble output AnnData
    # ------------------------------------------------------------------
    X_out = X_hvg_sp.tocsr().astype(np.float32)

    adata_out = ad.AnnData(X=X_out, obs=obs_out, var=var_hvg)
    adata_out.obsm["X_pca"] = X_pca

    adata_out.uns["pca"] = {
        "variance_ratio": var_ratio,
        "variance":       eigenvalues,
        "components":     components,   # (n_pcs, n_hvg)
        "gene_means":     gene_means,
    }
    adata_out.uns["dataset_id"]      = DATASET_ID
    adata_out.uns["time_axis"]       = "observed_time"
    adata_out.uns["source_h5ad"]     = str(input_h5ad)
    adata_out.uns["ipsc_abs_day_note"] = (
        "iPSC abs_day=16 is PROVISIONAL. "
        "Review before any formal benchmark run using abs_day for trajectory models."
    )
    adata_out.uns["preprocessing_summary"] = {
        "script":                "scripts/build_gse242424_hvg2000_input.py",
        "source_script":         "scripts/build_gse242424_raw_full_gene_input.py",
        "dataset":               DATASET_ID,
        "timestamp":             TIMESTAMP,
        "source_h5ad":           str(input_h5ad),
        "source_shape":          list(adata_full.shape),
        "hvg_selection_method":  "binned_residual_variance_seurat_style",
        "hvg_selection_details": (
            f"Top {N_HVG} genes by log10-var / expected-log10-var residual "
            f"in {N_HVG_BINS} log10-mean bins. "
            "Computed on log1p-normalised sparse CSR (no raw counts re-read)."
        ),
        "n_hvg":                 N_HVG,
        "n_pcs":                 N_PCS,
        "x_transform":           "log1p(total-count normalised to 1e4) -- inherited from source",
        "intended_use":          "model training / benchmark trajectory inference",
    }

    # ------------------------------------------------------------------
    # 7. Write h5ad
    # ------------------------------------------------------------------
    print(f"\nWriting output h5ad ...")
    adata_out.write_h5ad(OUTPUT_H5AD)
    size_mb = OUTPUT_H5AD.stat().st_size / 1e6
    print(f"  {OUTPUT_H5AD}")
    print(f"  File size  : {size_mb:.1f} MB")
    print(f"  Shape      : {adata_out.n_obs:,} x {adata_out.n_vars:,}")

    # ------------------------------------------------------------------
    # 8. Write summary JSON
    # ------------------------------------------------------------------
    summary = {
        "dataset_id":              DATASET_ID,
        "timestamp":               TIMESTAMP,
        "source_h5ad":             str(input_h5ad),
        "source_shape":            list(adata_full.shape),
        "output_h5ad":             str(OUTPUT_H5AD),
        "output_shape":            [int(adata_out.n_obs), int(adata_out.n_vars)],
        "n_hvg":                   N_HVG,
        "n_pcs":                   N_PCS,
        "hvg_selection_method":    "binned_residual_variance_seurat_style",
        "x_is_sparse":             sp.issparse(adata_out.X),
        "x_pca_shape":             list(adata_out.obsm["X_pca"].shape),
        "obs_columns":             list(adata_out.obs.columns),
        "var_columns":             list(adata_out.var.columns),
        "abs_day_values":          sorted(adata_out.obs["abs_day"].unique().tolist()),
        "ipsc_abs_day_note":       "iPSC abs_day=16 is PROVISIONAL",
        "top5_hvg_genes":          hvg_names[:5],
        "pca_var_explained_top10": [round(float(v), 6) for v in var_ratio[:10]],
        "runtime_s":               round(_time.time() - T0, 1),
    }
    with open(OUTPUT_JSON, "w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"  Summary JSON -> {OUTPUT_JSON}")

    # ------------------------------------------------------------------
    # 9. Inline validation
    # ------------------------------------------------------------------
    print("\nValidating ...")
    adata_check = ad.read_h5ad(OUTPUT_H5AD)
    checks = {
        "n_vars == 2000":       adata_check.n_vars == N_HVG,
        "X_pca present":        "X_pca" in adata_check.obsm,
        "X_pca n_components":   adata_check.obsm["X_pca"].shape[1] == N_PCS,
        "X is sparse":          sp.issparse(adata_check.X),
        "dataset_id correct":   adata_check.uns.get("dataset_id") == DATASET_ID,
        "all required obs OK":  all(c in adata_check.obs.columns for c in REQUIRED_OBS),
        "all required var OK":  all(c in adata_check.var.columns for c in REQUIRED_VAR),
        "no NaN in X.data":     np.isfinite(adata_check.X.data).all(),
    }
    all_ok = True
    for name, result in checks.items():
        status = "PASS" if result else "FAIL"
        print(f"  [{status}] {name}")
        if not result:
            all_ok = False

    abs_days = sorted(adata_check.obs["abs_day"].unique().tolist())
    print(f"  abs_day values : {abs_days}")

    elapsed = _time.time() - T0
    print()
    print("=" * 70)
    print(f"[{TIMESTAMP}] DONE -- {'ALL CHECKS PASSED' if all_ok else 'SOME CHECKS FAILED'}")
    print(f"  Output  : {OUTPUT_H5AD}")
    print(f"  Shape   : {adata_check.n_obs:,} x {adata_check.n_vars:,}")
    print(f"  X_pca   : {adata_check.obsm['X_pca'].shape}")
    print(f"  Runtime : {elapsed:.0f}s ({elapsed/60:.1f} min)")
    print(f"  RSS     : {rss_gb():.2f} GB")
    print("=" * 70)

    if not all_ok:
        sys.exit(1)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build GSE242424 HVG2000 benchmark input from full-gene h5ad."
    )
    parser.add_argument(
        "--input",
        metavar="PATH",
        default=None,
        help="Override path to the full-gene h5ad (default: project-root auto-detect).",
    )
    args = parser.parse_args()
    input_path = Path(args.input) if args.input else INPUT_H5AD
    main(input_h5ad=input_path)
