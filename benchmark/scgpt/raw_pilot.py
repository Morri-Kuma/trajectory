"""
raw_pilot.py
scGPT pilot embedding — GSE230659 raw full-gene h5ad (medium-scale)

Purpose
-------
A controlled pilot embedding run on the full-expression raw-count file
(20260329_1343_GSE230659_raw.h5ad) at ~7,500-cell scale. This is the
intermediate step between the 1,005-cell smoke test and the full 75,194-cell
production run.

This stage answers:
  - Can the raw full-gene h5ad embed cleanly with scGPT on this machine?
  - Are embeddings stable and non-degenerate at realistic scale?
  - What GPU load / runtime is expected before we attempt the full dataset?
  - Do UMAP / Leiden structures look usable for the later pseudo-state step?

What this stage does NOT do:
  - Does not define or freeze pseudo-states
  - Does not build a reference graph
  - Does not touch the active WOT benchmark configuration

Pipeline
--------
  1. Load 20260329_1343_GSE230659_raw.h5ad
  2. Build balanced pilot subset: up to 500 cells per timepoint (~7,500 total)
  3. Report vocabulary overlap (raw h5ad has ~27K genes vs vocab of ~60K)
  4. Run scGPT embedding on GPU (whole-human pretrained model)
  5. Run sc.pp.neighbors / sc.tl.umap / sc.tl.leiden
  6. Save h5ad, metadata JSON, UMAP figures, cluster-by-timepoint table

Environment
-----------
  conda activate scgpt_env
  cd C:\\Users\\37620\\trajectory
  python benchmark\\scgpt\\raw_pilot.py

Outputs (under benchmark/results/scgpt/raw_pilot/)
---------------------------------------------------
  raw_pilot_subset.h5ad              — AnnData with X_scGPT, UMAP, Leiden
  raw_pilot_metadata.json            — full run provenance and QC numbers
  umap_by_timepoint.png              — UMAP coloured by abs_day
  umap_by_leiden.png                 — UMAP coloured by Leiden cluster
  cluster_timepoint_summary.csv      — cells per (Leiden cluster × timepoint)
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch


# ---------------------------------------------------------------------------
# 0. Project-root detection
# ---------------------------------------------------------------------------

def _find_project_root() -> Path:
    import os
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
print(f"[raw_pilot] Project root: {PROJECT_ROOT}")

# Make the local scGPT repo importable.
# Set SCGPT_REPO env var to override the default path (required on Shirokane / HPC).
import os as _os
SCGPT_REPO = Path(_os.environ.get("SCGPT_REPO", r"C:\Users\37620\Documents\GitHub\scGPT"))
del _os
if SCGPT_REPO.exists() and str(SCGPT_REPO) not in sys.path:
    sys.path.insert(0, str(SCGPT_REPO))

import scgpt as scg  # noqa: E402


# ---------------------------------------------------------------------------
# 1. Paths and parameters
# ---------------------------------------------------------------------------

RAW_H5AD_PATH = (
    PROJECT_ROOT / "data" / "processed" / "20260329_1343_GSE230659_raw.h5ad"
)
MODEL_DIR  = PROJECT_ROOT / "models" / "scgpt_whole_human"
OUTPUT_DIR = PROJECT_ROOT / "benchmark" / "results" / "scgpt" / "raw_pilot"

# Embedding parameters
GENE_COL             = "index"   # HGNC symbols live in raw h5ad var.index
CELLS_PER_TIMEPOINT  = 500       # 500 × 15 timepoints = 7,500 cells
TIME_KEY             = "abs_day" # raw h5ad uses abs_day, not time_label
BATCH_SIZE           = 64        # reduce to 32 if CUDA OOM
USE_FAST_TRANSFORMER = False      # flash_attn absent in scgpt_env
LEIDEN_RESOLUTION    = 0.5

# ---------------------------------------------------------------------------
# 2. Pre-flight checks
# ---------------------------------------------------------------------------

print("\n[raw_pilot] === Pre-flight checks ===")

if not RAW_H5AD_PATH.exists():
    raise FileNotFoundError(
        f"Raw h5ad not found:\n  {RAW_H5AD_PATH}\n"
        "Confirm the filename matches the file on disk."
    )
print(f"  Input h5ad : {RAW_H5AD_PATH}  ✓")

for fname in ["vocab.json", "args.json", "best_model.pt"]:
    if not (MODEL_DIR / fname).exists():
        raise FileNotFoundError(f"Missing model file: {MODEL_DIR / fname}")
print(f"  Model dir  : {MODEL_DIR}  ✓")

# GPU / CUDA status
if not torch.cuda.is_available():
    print("  WARNING: CUDA not available — falling back to CPU. This will be slow.")
    _device = "cpu"
    gpu_name     = "CPU"
    cuda_version = "N/A"
else:
    _device      = "cuda"
    gpu_name     = torch.cuda.get_device_name(0)
    cuda_version = torch.version.cuda
    gpu_mem_total_gb = torch.cuda.get_device_properties(0).total_memory / 1e9
    print(f"  GPU        : {gpu_name}  ✓")
    print(f"  CUDA       : {cuda_version}")
    print(f"  VRAM total : {gpu_mem_total_gb:.1f} GB")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"  Output dir : {OUTPUT_DIR}  ✓")

# ---------------------------------------------------------------------------
# 3. Load raw h5ad and build balanced pilot subset
# ---------------------------------------------------------------------------

print(f"\n[raw_pilot] === Loading {RAW_H5AD_PATH.name} ===")
adata_full = sc.read_h5ad(RAW_H5AD_PATH)
print(f"  Full dataset: {adata_full.n_obs} cells × {adata_full.n_vars} genes")

if TIME_KEY not in adata_full.obs.columns:
    raise KeyError(
        f"Time column {TIME_KEY!r} not found in obs.\n"
        f"Available: {list(adata_full.obs.columns)}"
    )

timepoints = sorted(adata_full.obs[TIME_KEY].unique())
print(f"  Timepoints  : {len(timepoints)} ({timepoints[0]} → {timepoints[-1]})")

# Balanced sampling: up to CELLS_PER_TIMEPOINT from each timepoint
print(f"\n[raw_pilot] === Building pilot subset ({CELLS_PER_TIMEPOINT} cells/timepoint) ===")

np.random.seed(42)
subset_idx = []

for tp in timepoints:
    tp_mask    = adata_full.obs[TIME_KEY] == tp
    tp_indices = np.where(tp_mask)[0]
    n_sample   = min(CELLS_PER_TIMEPOINT, len(tp_indices))
    chosen     = np.random.choice(tp_indices, size=n_sample, replace=False)
    subset_idx.extend(chosen.tolist())
    print(f"  abs_day {tp:6.2f} : sampled {n_sample:4d} / {len(tp_indices):5d} cells")

subset_idx = sorted(subset_idx)
adata = adata_full[subset_idx].copy()
del adata_full  # free memory before loading the 205 MB model

print(f"\n  Pilot subset: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"  Timepoints covered: {adata.obs[TIME_KEY].nunique()} / {len(timepoints)}")

# ---------------------------------------------------------------------------
# 4. Vocabulary overlap check
# ---------------------------------------------------------------------------

print("\n[raw_pilot] === Vocabulary overlap check ===")

vocab_path = MODEL_DIR / "vocab.json"
with open(vocab_path, "r", encoding="utf-8") as f:
    vocab = json.load(f)

special_tokens = {"<pad>", "<cls>", "<eoc>"}
vocab_genes    = set(vocab.keys()) - special_tokens

# gene_col="index" means gene symbols are in the var.index
gene_symbols = set(adata.var.index)
overlap      = gene_symbols & vocab_genes
n_overlap    = len(overlap)
n_total      = len(gene_symbols)
overlap_frac = n_overlap / n_total if n_total > 0 else 0.0

print(f"  Vocab gene count (excl. special tokens): {len(vocab_genes):,}")
print(f"  Genes in raw h5ad var.index            : {n_total:,}")
print(f"  Genes matching vocab                   : {n_overlap:,} / {n_total:,}  ({100*overlap_frac:.1f}%)")

# Also report how many vocab-matching genes are above minimal expression
# (n_cells column records how many cells express each gene)
if "n_cells" in adata.var.columns:
    expressed_mask   = adata.var.index.isin(vocab_genes)
    n_expressed_vocab = expressed_mask.sum()
    print(f"  Vocab-matching genes (any expression in pilot subset): {n_expressed_vocab:,}")

if overlap_frac < 0.3:
    print("  WARNING: <30% overlap — check that var.index contains HGNC symbols.")
else:
    print("  Overlap is strong.  ✓")

# ---------------------------------------------------------------------------
# 5. Capture GPU memory baseline
# ---------------------------------------------------------------------------

if _device == "cuda":
    torch.cuda.reset_peak_memory_stats()
    gpu_mem_before_mb = torch.cuda.memory_allocated() / 1e6
    print(f"\n[raw_pilot] GPU memory before embedding: {gpu_mem_before_mb:.1f} MB")

# ---------------------------------------------------------------------------
# 6. Run scGPT embedding
# ---------------------------------------------------------------------------

print(f"\n[raw_pilot] === Running scGPT embedding ===")
print(f"  gene_col={GENE_COL!r}  batch_size={BATCH_SIZE}  device={_device}")
print("  use_fast_transformer=False  (flash_attn absent)")

t_embed_start = time.time()

adata = scg.tasks.embed_data(
    adata,
    MODEL_DIR,
    gene_col             = GENE_COL,
    max_length           = 1200,
    batch_size           = BATCH_SIZE,
    device               = _device,
    use_fast_transformer = USE_FAST_TRANSFORMER,
    return_new_adata     = False,   # add X_scGPT in-place
)

embed_runtime_sec = time.time() - t_embed_start
print(f"  Embedding complete in {embed_runtime_sec:.1f}s  "
      f"({embed_runtime_sec/60:.1f} min)")

# Capture peak GPU memory
if _device == "cuda":
    gpu_mem_after_mb = torch.cuda.memory_allocated() / 1e6
    gpu_mem_peak_mb  = torch.cuda.max_memory_allocated() / 1e6
    print(f"  GPU memory after embedding : {gpu_mem_after_mb:.1f} MB")
    print(f"  GPU peak memory            : {gpu_mem_peak_mb:.1f} MB")
else:
    gpu_mem_before_mb = gpu_mem_after_mb = gpu_mem_peak_mb = -1.0

# Validate embedding
if "X_scGPT" not in adata.obsm:
    raise RuntimeError("X_scGPT not in adata.obsm — embed_data failed silently.")

emb      = adata.obsm["X_scGPT"]
emb_shape = emb.shape
norms     = np.linalg.norm(emb, axis=1)
norm_mean = float(norms.mean())
norm_std  = float(norms.std())
has_nan   = bool(np.isnan(emb).any())
has_inf   = bool(np.isinf(emb).any())

print(f"\n  Embedding shape  : {emb_shape}")
print(f"  L2 norm (mean±std): {norm_mean:.4f} ± {norm_std:.4f}  (expected ≈ 1.0000 ± 0.0000)")
print(f"  NaN: {has_nan}  |  Inf: {has_inf}")

if has_nan or has_inf:
    print("  WARNING: degenerate values in embedding — inspect before proceeding.")
if abs(norm_mean - 1.0) > 0.01 or norm_std > 0.01:
    print("  WARNING: L2 norms deviate from 1.0 — embeddings may not be L2-normalized.")

# ---------------------------------------------------------------------------
# 7. Downstream Scanpy: neighbors → UMAP → Leiden
# ---------------------------------------------------------------------------

print("\n[raw_pilot] === Neighbors / UMAP / Leiden ===")
sc.settings.verbosity = 1

sc.pp.neighbors(adata, use_rep="X_scGPT", n_neighbors=30, random_state=42)
print("  neighbors (n_neighbors=30): done")

sc.tl.umap(adata, random_state=42)
print("  UMAP: done")

sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, random_state=42,
             key_added="leiden_scgpt")
n_leiden = adata.obs["leiden_scgpt"].nunique()
print(f"  Leiden (res={LEIDEN_RESOLUTION}): {n_leiden} clusters")

# ---------------------------------------------------------------------------
# 8. Cluster × timepoint summary table
# ---------------------------------------------------------------------------

summary = (
    adata.obs
    .groupby(["leiden_scgpt", TIME_KEY], observed=True)
    .size()
    .unstack(fill_value=0)
)
summary.index.name   = "leiden_cluster"
summary.columns.name = "abs_day"

# Add totals
summary["total_cells"] = summary.sum(axis=1)
summary.loc["total"]   = summary.sum(axis=0)

# ---------------------------------------------------------------------------
# 9. Save outputs
# ---------------------------------------------------------------------------

print("\n[raw_pilot] === Saving outputs ===")

# 9a. Pilot AnnData
h5ad_out = OUTPUT_DIR / "raw_pilot_subset.h5ad"
adata.write_h5ad(h5ad_out)
print(f"  Saved h5ad               : {h5ad_out}")

# 9b. Cluster × timepoint summary CSV
csv_out = OUTPUT_DIR / "cluster_timepoint_summary.csv"
summary.to_csv(csv_out)
print(f"  Saved cluster summary    : {csv_out}")

# 9c. Metadata JSON
metadata = {
    "run_timestamp"             : time.strftime("%Y-%m-%dT%H:%M:%S"),
    "script"                    : "benchmark/scgpt/raw_pilot.py",
    # --- inputs ---
    "input_h5ad"                : str(RAW_H5AD_PATH),
    "model_dir"                 : str(MODEL_DIR),
    "gene_col"                  : GENE_COL,
    # --- hardware ---
    "device"                    : _device,
    "gpu_name"                  : gpu_name,
    "cuda_version"              : cuda_version,
    "gpu_mem_total_gb"          : round(gpu_mem_total_gb, 2) if _device == "cuda" else -1,
    "gpu_mem_before_embed_mb"   : round(gpu_mem_before_mb, 1),
    "gpu_mem_after_embed_mb"    : round(gpu_mem_after_mb, 1),
    "gpu_mem_peak_embed_mb"     : round(gpu_mem_peak_mb, 1),
    # --- embedding config ---
    "use_fast_transformer"      : USE_FAST_TRANSFORMER,
    "batch_size"                : BATCH_SIZE,
    "max_length"                : 1200,
    # --- subset ---
    "cells_per_timepoint_target": CELLS_PER_TIMEPOINT,
    "n_cells_pilot"             : int(adata.n_obs),
    "n_genes_raw"               : int(n_total),
    "n_timepoints"              : int(adata.obs[TIME_KEY].nunique()),
    # --- vocab ---
    "n_genes_in_vocab"          : int(n_overlap),
    "vocab_overlap_frac"        : round(overlap_frac, 4),
    # --- embedding quality ---
    "embedding_shape"           : list(emb_shape),
    "embedding_norm_mean"       : round(norm_mean, 6),
    "embedding_norm_std"        : round(norm_std, 6),
    "embedding_has_nan"         : has_nan,
    "embedding_has_inf"         : has_inf,
    # --- clustering ---
    "leiden_resolution"         : LEIDEN_RESOLUTION,
    "n_leiden_clusters"         : int(n_leiden),
    # --- runtime ---
    "embed_runtime_sec"         : round(embed_runtime_sec, 1),
    "embed_runtime_min"         : round(embed_runtime_sec / 60, 2),
}
meta_out = OUTPUT_DIR / "raw_pilot_metadata.json"
with open(meta_out, "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)
print(f"  Saved metadata JSON      : {meta_out}")

# 9d. UMAP figures
sc.settings.figdir = str(OUTPUT_DIR)

# Convert abs_day to string for categorical colouring
adata.obs["abs_day_str"] = adata.obs[TIME_KEY].astype(str)

sc.pl.umap(
    adata,
    color  = "abs_day_str",
    title  = f"scGPT pilot — by timepoint  (n={adata.n_obs})",
    palette= "tab20",
    show   = False,
    save   = "_by_timepoint.png",
)
print(f"  Saved UMAP (timepoint)   : {OUTPUT_DIR}/umap_by_timepoint.png")

sc.pl.umap(
    adata,
    color  = "leiden_scgpt",
    title  = f"scGPT pilot — Leiden res={LEIDEN_RESOLUTION}  ({n_leiden} clusters)",
    show   = False,
    save   = "_by_leiden.png",
)
print(f"  Saved UMAP (Leiden)      : {OUTPUT_DIR}/umap_by_leiden.png")

# ---------------------------------------------------------------------------
# 10. Final summary
# ---------------------------------------------------------------------------

print("\n[raw_pilot] === Pilot summary ===")
print(f"  Input file      : {RAW_H5AD_PATH.name}")
print(f"  Cells embedded  : {adata.n_obs}")
print(f"  Genes (raw)     : {n_total:,}")
print(f"  Vocab overlap   : {n_overlap:,} / {n_total:,}  ({100*overlap_frac:.1f}%)")
print(f"  Embedding shape : {emb_shape}")
print(f"  Norm mean±std   : {norm_mean:.4f} ± {norm_std:.4f}")
print(f"  Leiden clusters : {n_leiden}")
print(f"  Embed runtime   : {embed_runtime_sec:.1f}s ({embed_runtime_sec/60:.1f} min)")
if _device == "cuda":
    print(f"  GPU peak VRAM   : {gpu_mem_peak_mb:.0f} MB")

print(f"\n  Outputs → {OUTPUT_DIR}")
for f in sorted(OUTPUT_DIR.iterdir()):
    size_kb = f.stat().st_size / 1024
    print(f"    {f.name:<45} {size_kb:>8.1f} KB")

status = "PASSED ✓" if not (has_nan or has_inf) else "COMPLETED WITH WARNINGS ⚠"
print(f"\n[raw_pilot] {status}")
