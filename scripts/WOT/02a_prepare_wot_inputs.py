#!/usr/bin/env python3
"""
02a_prepare_wot_inputs.py — Prepare all inputs required by WaddingtonOT.

Project : Comparative Study of Trajectory Inference Models for Chemical iPSC Reprogramming
Dataset : GSE230659 (human, Liuyang et al. 2023 Cell Stem Cell)
Step    : 2a — Reconstruct AnnData from Step 01 outputs, write WOT-format input files
Author  : Kuma
Env     : conda activate traj_env

Prerequisites:
    Step 01 outputs in data/processed/:
        YYYYMMDD_HHMM_GSE230659_qc_obs.csv
        YYYYMMDD_HHMM_GSE230659_qc_var.csv
        YYYYMMDD_HHMM_GSE230659_qc_hvg_matrix.npz
        YYYYMMDD_HHMM_GSE230659_qc_pca.npz

    pip install wot anndata scipy  (all available in traj_env)

Usage:
    conda activate traj_env
    cd C:\\Users\\37620\\trajectory
    python scripts/WOT/02a_prepare_wot_inputs.py

Outputs (all timestamped):
    data/processed/  YYYYMMDD_HHMM_GSE230659_wot_ready.h5ad  — AnnData for WOT
    data/processed/  YYYYMMDD_HHMM_cell_days.txt              — WOT cell days file
    results/figures/ YYYYMMDD_HHMM_timepoint_cells.png        — Cells/timepoint summary
"""

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
import os, gc, glob
from datetime import datetime
from pathlib import Path

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")

_candidates = [
    Path(r"C:\Users\37620\trajectory"),
    Path("/sessions/cool-admiring-hawking/mnt/trajectory"),
    Path(__file__).resolve().parents[2],
]
PROJECT_ROOT = next((p for p in _candidates if p.exists()), _candidates[-1])
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
FIGURES_DIR   = PROJECT_ROOT / "results" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ---- Locate Step 01 outputs (most recent timestamp) ----
def find_latest(pattern: str) -> Path:
    matches = sorted(PROCESSED_DIR.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No files matching '{pattern}' in {PROCESSED_DIR}.\n"
            f"Run scripts/C_traj/01_load_and_qc.py first."
        )
    return matches[-1]  # most recent by filename timestamp

OBS_PATH    = find_latest("*_GSE230659_qc_obs.csv")
VAR_PATH    = find_latest("*_GSE230659_qc_var.csv")
MATRIX_PATH = find_latest("*_GSE230659_qc_hvg_matrix.npz")
PCA_PATH    = find_latest("*_GSE230659_qc_pca.npz")

# ---- WOT hyperparameters (baseline run) ----
# Using PCA space (30 dims) as recommended for WOT on scRNA-seq
N_PCS_WOT = 30   # PCA dims to pass to WOT
EPSILON   = 0.05  # entropic regularization
LAMBDA1   = 1.0   # unbalanced OT source marginal penalty
LAMBDA2   = 50.0  # unbalanced OT target marginal penalty

# =============================================================================
# 1. IMPORTS
# =============================================================================
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

print(f"[{TIMESTAMP}] WOT Input Preparation")
print(f"  Project root : {PROJECT_ROOT}")
print(f"  obs          : {OBS_PATH.name}")
print(f"  var          : {VAR_PATH.name}")
print(f"  matrix       : {MATRIX_PATH.name}")
print(f"  pca          : {PCA_PATH.name}")

# =============================================================================
# 2. LOAD STEP 01 OUTPUTS
# =============================================================================
print("\n--- Loading Step 01 outputs ---")

obs = pd.read_csv(OBS_PATH, index_col=0)
var = pd.read_csv(VAR_PATH, index_col=0)

hvg_data  = np.load(MATRIX_PATH, allow_pickle=True)
X_log     = hvg_data["X_log"]           # (n_cells, 2000) float32
hvg_names = hvg_data["hvg_gene_names"]  # (2000,) str

pca_data  = np.load(PCA_PATH, allow_pickle=True)
X_pca     = pca_data["X_pca"]           # (n_cells, 50) float32
var_ratio = pca_data["variance_ratio"]  # (50,)

# Subset var to HVGs only (in same order as matrix columns)
var_hvg = var[var["highly_variable"]].copy()
var_hvg = var_hvg.reindex(hvg_names)

print(f"  Cells: {obs.shape[0]:,}")
print(f"  HVGs:  {len(hvg_names):,}")
print(f"  PCs:   {X_pca.shape[1]}")
print(f"  Time points: {sorted(obs['abs_day'].unique())}")

# =============================================================================
# 3. BUILD ANNDATA
# =============================================================================
print("\n--- Building AnnData ---")

adata = ad.AnnData(
    X   = sp.csr_matrix(X_log),
    obs = obs.copy(),
    var = var_hvg.copy(),
)

# Store PCA in obsm
adata.obsm["X_pca"] = X_pca.astype(np.float32)

# Store PCA metadata in uns
adata.uns["pca"] = {
    "variance_ratio": var_ratio,
    "variance": pca_data["eigenvalues"],
    "params": {"n_comps": len(var_ratio), "n_hvg": len(hvg_names)},
}

# ── F1: add 'day' column as float alias of 'abs_day' ─────────────────────────
# This ensures downstream steps (02b, 02c) can find a 'day' column without
# having to rediscover 'abs_day' and rename it at runtime.
adata.obs["day"] = adata.obs["abs_day"].astype(float)
print(f"  Added obs['day'] = obs['abs_day'].astype(float)")

print(f"  AnnData: {adata.n_obs:,} x {adata.n_vars:,}")
print(f"  obs columns: {list(adata.obs.columns)}")
print(f"  obsm keys:   {list(adata.obsm.keys())}")

# ── F2: consistency assertions ─────────────────────────────────────────────────
print("\n--- Consistency checks ---")

assert X_log.shape[0] == obs.shape[0], (
    f"Row count mismatch: X_log.shape[0]={X_log.shape[0]} "
    f"vs obs.shape[0]={obs.shape[0]}"
)
assert X_pca.shape[0] == obs.shape[0], (
    f"Row count mismatch: X_pca.shape[0]={X_pca.shape[0]} "
    f"vs obs.shape[0]={obs.shape[0]}"
)
assert X_log.shape[1] == len(hvg_names), (
    f"Column count mismatch: X_log.shape[1]={X_log.shape[1]} "
    f"vs len(hvg_names)={len(hvg_names)}"
)

# var_hvg index must exactly match hvg_names
_var_idx = var_hvg.index.tolist()
assert _var_idx == list(hvg_names), (
    f"var_hvg.index does not match hvg_names.\n"
    f"First mismatch at position "
    f"{next(i for i,(a,b) in enumerate(zip(_var_idx, hvg_names)) if a!=b)}"
)

# var_hvg must not have unexpected null rows after reindexing
_null_rows = var_hvg.index[var_hvg.isnull().all(axis=1)]
if len(_null_rows) > 0:
    raise ValueError(
        f"var_hvg has {len(_null_rows)} all-null rows after reindex. "
        f"First: {_null_rows[0]}. Check that hvg_names are present in var."
    )

print(f"  X_log shape    : {X_log.shape}  ✓")
print(f"  X_pca shape    : {X_pca.shape}  ✓")
print(f"  var_hvg shape  : {var_hvg.shape}  ✓")
print(f"  var_hvg index  : matches hvg_names  ✓")
print(f"  var_hvg nulls  : none  ✓")

# Verify the critical 'abs_day' column exists and time range
assert "abs_day" in adata.obs.columns, "CRITICAL: abs_day column missing from obs"
print(f"  Time range     : [{adata.obs['abs_day'].min()}, {adata.obs['abs_day'].max()}]  ✓")

# =============================================================================
# 4. WRITE WOT-FORMAT CELL DAYS FILE
# =============================================================================
# WOT requires a two-column TSV: id (barcode) + day (numeric)
print("\n--- Writing cell_days.txt ---")

cell_days = adata.obs[["abs_day"]].copy()
cell_days.index.name = "id"
cell_days.columns = ["day"]

cell_days_path = PROCESSED_DIR / f"{TIMESTAMP}_cell_days.txt"
cell_days.to_csv(cell_days_path, sep="\t")
print(f"  Saved → {cell_days_path}")
print(f"  Sample:\n{cell_days.head(3).to_string()}")
print(f"\n  Cells per timepoint:")
for day, n in cell_days["day"].value_counts().sort_index().items():
    print(f"    Day {day:5.2f}: {n:>6,} cells")

# =============================================================================
# 5. FLAG LOW-COVERAGE TIMEPOINTS FOR WOT
# =============================================================================
print("\n--- Timepoint coverage analysis ---")
time_summary = (
    adata.obs.groupby(["abs_day", "stage_day_label"], observed=True)
    .size()
    .reset_index(name="n_cells")
    .sort_values("abs_day")
)

# WOT quality flag: flag pairs where min/max < 0.20
time_summary["n_cells_next"] = time_summary["n_cells"].shift(-1)
time_summary["ratio"] = (
    time_summary[["n_cells", "n_cells_next"]].min(axis=1) /
    time_summary[["n_cells", "n_cells_next"]].max(axis=1)
)

print("\n  Transport pair quality:")
for _, row in time_summary.dropna(subset=["n_cells_next"]).iterrows():
    flag = " ⚠  LOW" if row["ratio"] < 0.20 else ""
    print(f"    {row['stage_day_label']:<28} → next: "
          f"{row['n_cells']:>5,} ↔ {row['n_cells_next']:>5.0f}  "
          f"ratio={row['ratio']:.2f}{flag}")

# =============================================================================
# 6. SAVE ANNDATA FOR WOT
# =============================================================================
print("\n--- Saving WOT-ready AnnData ---")

h5ad_path = PROCESSED_DIR / f"{TIMESTAMP}_GSE230659_wot_ready.h5ad"
adata.write_h5ad(h5ad_path, compression="gzip")
print(f"  Saved → {h5ad_path}")

# =============================================================================
# 7. DIAGNOSTIC PLOT — Cells per timepoint with coverage flags
# =============================================================================
print("\n--- Generating diagnostic plot ---")

STAGE_COLORS = {
    "StageI":   "#2196F3",
    "StageII":  "#FF9800",
    "StageIII": "#4CAF50",
    "hCiPSCs":  "#E91E63",
}

fig, ax = plt.subplots(figsize=(15, 5))
colors = [STAGE_COLORS.get(s, "#999")
          for s in adata.obs.groupby("abs_day", observed=True)["stage"].first()
          .sort_index().values]
x = range(len(time_summary))
bars = ax.bar(x, time_summary["n_cells"], color=colors, edgecolor="white", linewidth=0.5)

# Annotate bars
for i, (_, row) in enumerate(time_summary.iterrows()):
    ax.text(i, row["n_cells"] + 80, f"{row['n_cells']:,}",
            ha="center", va="bottom", fontsize=8)
    if pd.notna(row.get("ratio")) and row.get("ratio", 1.0) < 0.20:
        ax.text(i, row["n_cells"] / 2, "⚠", ha="center", va="center",
                fontsize=16, color="red")

ax.set_xticks(x)
ax.set_xticklabels(time_summary["stage_day_label"], rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Number of cells (post-QC)")
ax.set_title("Cell counts per timepoint — WOT readiness check\n"
             "(⚠ = transport pair ratio < 0.20, may need attention)")
ax.axhline(1000, color="red", linestyle="--", alpha=0.4, label="1,000 cell threshold")

for stage, color in STAGE_COLORS.items():
    ax.bar(0, 0, color=color, label=stage)
ax.legend(loc="upper right", fontsize=9)

fig.tight_layout()
fig_path = FIGURES_DIR / f"{TIMESTAMP}_timepoint_cells_wot.png"
fig.savefig(fig_path, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {fig_path}")

# =============================================================================
# 8. SUMMARY
# =============================================================================
print(f"\n{'='*60}")
print(f"[{TIMESTAMP}] Step 02a COMPLETE")
print(f"  AnnData       : {h5ad_path.name}")
print(f"  cell_days.txt : {cell_days_path.name}")
print(f"  AnnData shape : {adata.n_obs:,} cells x {adata.n_vars:,} HVGs")
print(f"  obs['day']    : added  (float alias of abs_day)")
print(f"  PCA dims used : {N_PCS_WOT} (of {X_pca.shape[1]} computed)")
print(f"\n  WOT baseline parameters:")
print(f"    epsilon = {EPSILON}  (entropic regularization)")
print(f"    lambda1 = {LAMBDA1}  (source marginal penalty)")
print(f"    lambda2 = {LAMBDA2}  (target marginal penalty)")
print(f"\n  Next: run 02b_estimate_growth_rates.py")
print(f"{'='*60}")
