"""
full_embed.py
scGPT full embedding — GSE230659 (all 75,194 cells, raw full-gene h5ad)

Purpose
-------
Compute the complete scGPT embedding for all cells in the GSE230659 dataset
using the raw full-expression h5ad as input. This script produces the stable
embedding-layer artifact that all later annotation and pseudo-state scripts
will read from without re-running the GPU step.

What this script does
---------------------
  1.  Pre-flight checks (paths, model files, CUDA)
  2.  Load 20260329_1343_GSE230659_raw.h5ad  (75,194 cells × 27,267 genes)
  3.  Vocabulary overlap report
  4.  Full scGPT embedding on GPU  →  obsm["X_scGPT"]  (75,194 × 512)
  5.  Checkpoint save of the raw embedding array  (safety net)
  6.  Embedding validation (shape, L2 norms, NaN/Inf)
  7.  sc.pp.neighbors  →  sc.tl.umap  →  sc.tl.leiden
      Leiden result stored in  obs["leiden_scgpt_res0.5"]
  8.  Save stable artifact:  adata_scgpt_full.h5ad
  9.  Save cluster × timepoint matrix CSV
  10. Save full_embed_metadata.json
  11. Save UMAP figures

What this script does NOT do
-----------------------------
  - Does not assign provisional pseudo-state names
  - Does not read or write cluster_map_provisional.json
  - Does not merge or split clusters
  - Does not touch the active WOT / CellRank2 benchmark configuration

Environment
-----------
  conda activate scgpt_env
  cd C:\\Users\\37620\\trajectory
  python benchmark\\scgpt\\full_embed.py

Expected runtime
----------------
  Embedding  : ~8–12 min  (extrapolated from 58 s / 7,500-cell pilot)
  Neighbors  : ~1–2 min
  UMAP       : ~3–5 min
  Total      : ~15–20 min

Outputs  (benchmark/results/scgpt/full/)
-----------------------------------------
  adata_scgpt_full.h5ad            stable embedding artifact (main output)
  full_embed_metadata.json         run provenance and QC numbers
  cluster_timepoint_matrix.csv     Leiden cluster × abs_day cell counts
  umap_by_timepoint.png            UMAP coloured by abs_day
  umap_by_leiden.png               UMAP coloured by leiden_scgpt_res0.5
  [checkpoint] X_scGPT_checkpoint.npy   raw embedding array (deleted on
                                         success; kept if downstream fails)
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
    known = Path(r"C:\Users\37620\trajectory")
    if known.exists() and (known / "data").exists():
        return known
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "data").exists():
            return candidate
    return here.parent


PROJECT_ROOT = _find_project_root()
print(f"[full_embed] Project root: {PROJECT_ROOT}")

# Make the local scGPT repo importable (same pattern as smoke_test / raw_pilot)
SCGPT_REPO = Path(r"C:\Users\37620\Documents\GitHub\scGPT")
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
OUTPUT_DIR = PROJECT_ROOT / "benchmark" / "results" / "scgpt" / "full"

# Embedding / clustering parameters
GENE_COL             = "index"   # HGNC symbols are in raw h5ad var.index
TIME_KEY             = "abs_day" # time column in raw h5ad obs
LEIDEN_RESOLUTION    = 0.5
LEIDEN_KEY           = f"leiden_scgpt_res{LEIDEN_RESOLUTION}"  # resolution-explicit column
N_NEIGHBORS          = 30        # larger than default; appropriate for 75K cells
BATCH_SIZE           = 64        # confirmed safe from pilot (peak VRAM ~1,360 MB)
USE_FAST_TRANSFORMER = False      # flash_attn absent in scgpt_env
RANDOM_SEED          = 42

# Checkpoint file: raw embedding array saved immediately after GPU step.
# Deleted on successful completion; kept if downstream steps fail so the
# ~10-minute GPU run does not need to be repeated.
CHECKPOINT_NPY = OUTPUT_DIR / "X_scGPT_checkpoint.npy"

# ---------------------------------------------------------------------------
# 2. Pre-flight checks
# ---------------------------------------------------------------------------

print("\n[full_embed] === Pre-flight checks ===")

if not RAW_H5AD_PATH.exists():
    raise FileNotFoundError(
        f"Raw h5ad not found:\n  {RAW_H5AD_PATH}\n"
        "Confirm the filename matches the file on disk."
    )
print(f"  Input h5ad  : {RAW_H5AD_PATH}  ✓")

for fname in ["vocab.json", "args.json", "best_model.pt"]:
    fp = MODEL_DIR / fname
    if not fp.exists():
        raise FileNotFoundError(
            f"Required model file missing: {fp}\n"
            "Confirm that models/scgpt_whole_human is complete."
        )
print(f"  Model dir   : {MODEL_DIR}  ✓")

if not torch.cuda.is_available():
    # Warn loudly — running 75 K cells on CPU would take many hours
    print(
        "\n  WARNING: CUDA not available.\n"
        "  Embedding 75,194 cells on CPU may take several hours.\n"
        "  Proceeding anyway — press Ctrl-C to abort.\n"
    )
    _device      = "cpu"
    gpu_name     = "CPU"
    cuda_version = "N/A"
    gpu_mem_total_gb = -1.0
else:
    _device          = "cuda"
    gpu_name         = torch.cuda.get_device_name(0)
    cuda_version     = torch.version.cuda
    gpu_mem_total_gb = torch.cuda.get_device_properties(0).total_memory / 1e9
    print(f"  GPU         : {gpu_name}  ✓")
    print(f"  CUDA        : {cuda_version}")
    print(f"  VRAM total  : {gpu_mem_total_gb:.1f} GB")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"  Output dir  : {OUTPUT_DIR}  ✓")

# ---------------------------------------------------------------------------
# 3. Load full raw h5ad
# ---------------------------------------------------------------------------

print(f"\n[full_embed] === Loading {RAW_H5AD_PATH.name} ===")
t_load_start = time.time()
adata = sc.read_h5ad(RAW_H5AD_PATH)
load_time = time.time() - t_load_start

print(f"  Loaded in {load_time:.1f}s")
print(f"  Shape      : {adata.n_obs} cells × {adata.n_vars} genes")
print(f"  obs columns: {list(adata.obs.columns)}")

if TIME_KEY not in adata.obs.columns:
    raise KeyError(
        f"Time column {TIME_KEY!r} not found in obs.\n"
        f"Available columns: {list(adata.obs.columns)}"
    )

timepoints = sorted(adata.obs[TIME_KEY].unique())
print(f"  Timepoints : {len(timepoints)}  "
      f"(abs_day {timepoints[0]} → {timepoints[-1]})")

tp_counts = adata.obs[TIME_KEY].value_counts().sort_index()
for tp, n in tp_counts.items():
    print(f"    abs_day {tp:6.2f} : {n:6,} cells")

# ---------------------------------------------------------------------------
# 4. Vocabulary overlap check
# ---------------------------------------------------------------------------

print("\n[full_embed] === Vocabulary overlap check ===")

with open(MODEL_DIR / "vocab.json", "r", encoding="utf-8") as f:
    vocab = json.load(f)

special_tokens = {"<pad>", "<cls>", "<eoc>"}
vocab_genes    = set(vocab.keys()) - special_tokens

# gene_col="index": HGNC symbols live in var.index for the raw h5ad
gene_symbols = set(adata.var.index)
overlap      = gene_symbols & vocab_genes
n_overlap    = len(overlap)
n_total      = len(gene_symbols)
overlap_frac = n_overlap / n_total if n_total > 0 else 0.0

print(f"  Vocab gene count (excl. special tokens) : {len(vocab_genes):,}")
print(f"  Genes in raw h5ad var.index             : {n_total:,}")
print(f"  Genes matching vocab                    : "
      f"{n_overlap:,} / {n_total:,}  ({100*overlap_frac:.1f}%)")

if overlap_frac < 0.3:
    raise RuntimeError(
        f"Vocabulary overlap is only {100*overlap_frac:.1f}%. "
        "Check that gene_col='index' points to HGNC gene symbols, not Ensembl IDs."
    )
print("  Overlap is strong  ✓")

# ---------------------------------------------------------------------------
# 5. Capture GPU memory baseline and run embedding
# ---------------------------------------------------------------------------

if _device == "cuda":
    torch.cuda.reset_peak_memory_stats()
    gpu_mem_before_mb = torch.cuda.memory_allocated() / 1e6
    print(f"\n[full_embed] GPU memory before embedding : {gpu_mem_before_mb:.1f} MB")

print(f"\n[full_embed] === Running scGPT embedding ===")
print(f"  Cells      : {adata.n_obs:,}")
print(f"  gene_col   : {GENE_COL!r}")
print(f"  batch_size : {BATCH_SIZE}")
print(f"  device     : {_device}")
print(f"  use_fast_transformer : {USE_FAST_TRANSFORMER}  (flash_attn absent)")

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
print(f"\n  Embedding complete in {embed_runtime_sec:.1f}s "
      f"({embed_runtime_sec / 60:.1f} min)")

if _device == "cuda":
    gpu_mem_after_mb = torch.cuda.memory_allocated() / 1e6
    gpu_mem_peak_mb  = torch.cuda.max_memory_allocated() / 1e6
    print(f"  GPU memory after    : {gpu_mem_after_mb:.1f} MB")
    print(f"  GPU peak VRAM       : {gpu_mem_peak_mb:.1f} MB")
else:
    gpu_mem_before_mb = gpu_mem_after_mb = gpu_mem_peak_mb = -1.0

# ---------------------------------------------------------------------------
# 6. Checkpoint — save raw embedding array immediately after GPU step
#    so that neighbors / UMAP / Leiden failures do not require re-embedding
# ---------------------------------------------------------------------------

if "X_scGPT" not in adata.obsm:
    raise RuntimeError(
        "X_scGPT not found in adata.obsm after embed_data.\n"
        "The embedding call returned without populating obsm."
    )

print(f"\n[full_embed] Writing embedding checkpoint → {CHECKPOINT_NPY.name}")
np.save(CHECKPOINT_NPY, adata.obsm["X_scGPT"])
print(f"  Checkpoint written  ({CHECKPOINT_NPY.stat().st_size / 1e6:.0f} MB)")

# ---------------------------------------------------------------------------
# 7. Validate the embedding
# ---------------------------------------------------------------------------

emb       = adata.obsm["X_scGPT"]
emb_shape = emb.shape
norms     = np.linalg.norm(emb, axis=1)
norm_mean = float(norms.mean())
norm_std  = float(norms.std())
has_nan   = bool(np.isnan(emb).any())
has_inf   = bool(np.isinf(emb).any())

print(f"\n[full_embed] Embedding validation:")
print(f"  Shape           : {emb_shape}")
print(f"  L2 norm mean±std: {norm_mean:.6f} ± {norm_std:.6f}  (expected 1.0 ± 0.0)")
print(f"  NaN: {has_nan}  |  Inf: {has_inf}")

if has_nan or has_inf:
    print("  WARNING: degenerate values detected — inspect before proceeding.")
if abs(norm_mean - 1.0) > 0.01 or norm_std > 0.01:
    print("  WARNING: L2 norms deviate from expected 1.0.")

# ---------------------------------------------------------------------------
# 8. Downstream Scanpy: neighbors → UMAP → Leiden
# ---------------------------------------------------------------------------

print(f"\n[full_embed] === Neighbors / UMAP / Leiden ===")
sc.settings.verbosity = 1

t_downstream_start = time.time()

sc.pp.neighbors(
    adata,
    use_rep     = "X_scGPT",
    n_neighbors = N_NEIGHBORS,
    random_state= RANDOM_SEED,
)
print(f"  neighbors (n_neighbors={N_NEIGHBORS}): done")

sc.tl.umap(adata, random_state=RANDOM_SEED)
print("  UMAP: done")

# Resolution-explicit key: preserves ability to add other resolutions later
# without overwriting this run's result
sc.tl.leiden(
    adata,
    resolution  = LEIDEN_RESOLUTION,
    random_state= RANDOM_SEED,
    key_added   = LEIDEN_KEY,
)
n_leiden = adata.obs[LEIDEN_KEY].nunique()
print(f"  Leiden (res={LEIDEN_RESOLUTION}, key='{LEIDEN_KEY}'): {n_leiden} clusters")

downstream_runtime_sec = time.time() - t_downstream_start
print(f"  Downstream steps: {downstream_runtime_sec:.1f}s")

# ---------------------------------------------------------------------------
# 9. Save stable full embedding artifact
# ---------------------------------------------------------------------------

print(f"\n[full_embed] === Saving outputs ===")

h5ad_out = OUTPUT_DIR / "adata_scgpt_full.h5ad"
print(f"  Saving adata_scgpt_full.h5ad  (this may take ~1–2 min for a large file) ...")
t_save_start = time.time()
adata.write_h5ad(h5ad_out)
save_time = time.time() - t_save_start
file_size_mb = h5ad_out.stat().st_size / 1e6
print(f"  Saved in {save_time:.1f}s  →  {h5ad_out.name}  ({file_size_mb:.0f} MB)")

# Checkpoint no longer needed — embedding is safely inside the h5ad
if CHECKPOINT_NPY.exists():
    CHECKPOINT_NPY.unlink()
    print("  Embedding checkpoint removed (data is in h5ad)")

# ---------------------------------------------------------------------------
# 10. Cluster × timepoint matrix
# ---------------------------------------------------------------------------

ct_matrix = (
    adata.obs
    .groupby([LEIDEN_KEY, TIME_KEY], observed=True)
    .size()
    .unstack(fill_value=0)
)
ct_matrix.index.name   = "leiden_cluster"
ct_matrix.columns.name = "abs_day"
ct_matrix["total_cells"] = ct_matrix.sum(axis=1)
ct_matrix.loc["total"]   = ct_matrix.sum(axis=0)

csv_out = OUTPUT_DIR / "cluster_timepoint_matrix.csv"
ct_matrix.to_csv(csv_out)
print(f"  Saved cluster_timepoint_matrix.csv")

# ---------------------------------------------------------------------------
# 11. Metadata JSON
# ---------------------------------------------------------------------------

total_runtime_sec = embed_runtime_sec + downstream_runtime_sec + save_time

metadata = {
    "run_timestamp"             : time.strftime("%Y-%m-%dT%H:%M:%S"),
    "script"                    : "benchmark/scgpt/full_embed.py",
    # --- inputs ---
    "input_h5ad"                : str(RAW_H5AD_PATH),
    "model_dir"                 : str(MODEL_DIR),
    "gene_col"                  : GENE_COL,
    # --- hardware ---
    "device"                    : _device,
    "gpu_name"                  : gpu_name,
    "cuda_version"              : cuda_version,
    "gpu_mem_total_gb"          : round(gpu_mem_total_gb, 2),
    "gpu_mem_before_embed_mb"   : round(gpu_mem_before_mb, 1),
    "gpu_mem_after_embed_mb"    : round(gpu_mem_after_mb, 1),
    "gpu_mem_peak_embed_mb"     : round(gpu_mem_peak_mb, 1),
    # --- embedding config ---
    "use_fast_transformer"      : USE_FAST_TRANSFORMER,
    "batch_size"                : BATCH_SIZE,
    "max_length"                : 1200,
    "random_seed"               : RANDOM_SEED,
    # --- data ---
    "n_cells"                   : int(adata.n_obs),
    "n_genes_raw"               : int(n_total),
    "n_timepoints"              : int(len(timepoints)),
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
    "leiden_key"                : LEIDEN_KEY,
    "n_leiden_clusters"         : int(n_leiden),
    "n_neighbors"               : N_NEIGHBORS,
    # --- runtime ---
    "embed_runtime_sec"         : round(embed_runtime_sec, 1),
    "embed_runtime_min"         : round(embed_runtime_sec / 60, 2),
    "downstream_runtime_sec"    : round(downstream_runtime_sec, 1),
    "save_runtime_sec"          : round(save_time, 1),
    "total_runtime_sec"         : round(total_runtime_sec, 1),
    "total_runtime_min"         : round(total_runtime_sec / 60, 2),
    # --- output ---
    "output_h5ad"               : str(h5ad_out),
    "output_h5ad_size_mb"       : round(file_size_mb, 1),
}

meta_out = OUTPUT_DIR / "full_embed_metadata.json"
with open(meta_out, "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)
print(f"  Saved full_embed_metadata.json")

# ---------------------------------------------------------------------------
# 12. UMAP figures
# ---------------------------------------------------------------------------

sc.settings.figdir = str(OUTPUT_DIR)

# Timepoint figure — convert abs_day to string for categorical palette
adata.obs["abs_day_str"] = adata.obs[TIME_KEY].astype(str)

sc.pl.umap(
    adata,
    color   = "abs_day_str",
    title   = f"scGPT full embedding — by timepoint  (n={adata.n_obs:,})",
    palette = "tab20",
    show    = False,
    save    = "_by_timepoint.png",
)
print("  Saved umap_by_timepoint.png")

# Leiden figure
sc.pl.umap(
    adata,
    color   = LEIDEN_KEY,
    title   = f"scGPT full embedding — Leiden res={LEIDEN_RESOLUTION}  "
               f"({n_leiden} clusters)",
    show    = False,
    save    = "_by_leiden.png",
)
print("  Saved umap_by_leiden.png")

# ---------------------------------------------------------------------------
# 13. Final summary
# ---------------------------------------------------------------------------

print(f"\n[full_embed] === Run complete ===")
print(f"  Input file      : {RAW_H5AD_PATH.name}")
print(f"  Cells embedded  : {adata.n_obs:,}")
print(f"  Genes (raw)     : {n_total:,}")
print(f"  Vocab overlap   : {n_overlap:,} / {n_total:,}  ({100*overlap_frac:.1f}%)")
print(f"  Embedding shape : {emb_shape}")
print(f"  Norm mean±std   : {norm_mean:.6f} ± {norm_std:.6f}")
print(f"  Leiden clusters : {n_leiden}  (key: '{LEIDEN_KEY}')")
print(f"  GPU peak VRAM   : {gpu_mem_peak_mb:.0f} MB"
      if _device == "cuda" else "  Device: CPU")
print(f"  Embed runtime   : {embed_runtime_sec:.1f}s  ({embed_runtime_sec/60:.1f} min)")
print(f"  Total runtime   : {total_runtime_sec:.1f}s  ({total_runtime_sec/60:.1f} min)")
print(f"\n  Outputs → {OUTPUT_DIR}")
for f in sorted(OUTPUT_DIR.iterdir()):
    if f.is_file():
        print(f"    {f.name:<45}  {f.stat().st_size / 1e6:>7.1f} MB")

status = "PASSED ✓" if not (has_nan or has_inf) else "COMPLETED WITH WARNINGS ⚠"
print(f"\n[full_embed] {status}")
print(
    "\n  Next step:\n"
    f"    1. Review cluster_timepoint_matrix.csv and UMAP figures\n"
    f"    2. Update benchmark/scgpt/cluster_map_provisional.json if cluster\n"
    f"       count at full scale ({n_leiden}) differs from the pilot (12)\n"
    f"    3. Run benchmark/scgpt/annotate_provisional.py to add pseudo-state\n"
    f"       labels without re-running this embedding"
)
