"""
smoke_test.py
scGPT smoke test — GSE230659 (adata_benchmark.h5ad subset)

Purpose
-------
End-to-end infrastructure validation for the scGPT embedding pipeline.
This script does NOT produce pseudo-states or modify any benchmark configuration.
It is strictly a pipeline sanity check.

Pipeline
--------
  1. Load adata_benchmark.h5ad
  2. Build a balanced ~1000-cell subset spanning all 15 timepoints
  3. Report vocabulary overlap between var["gene_symbol"] and scGPT vocab
  4. Run scGPT zero-shot embedding (GPU, use_fast_transformer=False)
  5. Run sc.pp.neighbors / sc.tl.umap / sc.tl.leiden on the embedding
  6. Save subset h5ad, metadata JSON, and UMAP figures

Environment
-----------
  conda activate scgpt_env
  cd C:\\Users\\37620\\trajectory
  python benchmark\\scgpt\\smoke_test.py

Outputs (under benchmark/results/scgpt/smoke_test/)
----------------------------------------------------
  smoke_test_subset.h5ad       — AnnData with X_scGPT, UMAP, Leiden labels
  smoke_test_metadata.json     — run metadata and QC numbers
  umap_by_timepoint.png        — UMAP coloured by time_label
  umap_by_leiden.png           — UMAP coloured by Leiden cluster
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import scanpy as sc
import torch


# ---------------------------------------------------------------------------
# 0. Project-root detection (same logic as other scripts in this repo)
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
print(f"[smoke_test] Project root: {PROJECT_ROOT}")

# Make the local scGPT repo importable if the package is not installed.
# Set SCGPT_REPO env var to override the default path (required on Shirokane / HPC).
import os as _os
SCGPT_REPO = Path(_os.environ.get("SCGPT_REPO", r"C:\Users\37620\Documents\GitHub\scGPT"))
del _os
if SCGPT_REPO.exists() and str(SCGPT_REPO) not in sys.path:
    sys.path.insert(0, str(SCGPT_REPO))

import scgpt as scg  # noqa: E402  (import after sys.path adjustment)


# ---------------------------------------------------------------------------
# 1. Paths — edit these only if your layout differs
# ---------------------------------------------------------------------------

H5AD_PATH  = PROJECT_ROOT / "data" / "processed" / "adata_benchmark.h5ad"
MODEL_DIR  = PROJECT_ROOT / "models" / "scgpt_whole_human"
OUTPUT_DIR = PROJECT_ROOT / "benchmark" / "results" / "scgpt" / "smoke_test"

# scGPT embedding parameters
GENE_COL            = "gene_symbol"   # var column with HGNC gene symbols
CELLS_PER_TIMEPOINT = 67              # ~67 × 15 timepoints ≈ 1,000 cells total
BATCH_SIZE          = 32              # conservative; reduce to 8 if OOM
DEVICE              = "cuda"
USE_FAST_TRANSFORMER = False          # flash_attn absent in this environment
LEIDEN_RESOLUTION   = 0.5

# ---------------------------------------------------------------------------
# 2. Pre-flight checks
# ---------------------------------------------------------------------------

print("\n[smoke_test] === Pre-flight checks ===")

if not H5AD_PATH.exists():
    raise FileNotFoundError(f"adata_benchmark.h5ad not found at:\n  {H5AD_PATH}")
print(f"  Input h5ad : {H5AD_PATH}  ✓")

for required_file in ["vocab.json", "args.json", "best_model.pt"]:
    fp = MODEL_DIR / required_file
    if not fp.exists():
        raise FileNotFoundError(
            f"Required model file missing: {fp}\n"
            "Confirm that the scgpt_whole_human directory is complete."
        )
print(f"  Model dir  : {MODEL_DIR}  ✓ (vocab.json, args.json, best_model.pt present)")

if DEVICE == "cuda":
    if not torch.cuda.is_available():
        print("  WARNING: CUDA requested but torch.cuda.is_available() is False.")
        print("           Falling back to CPU. This will be slow.")
        _device = "cpu"
    else:
        _device = "cuda"
        print(f"  CUDA       : available  ✓  ({torch.cuda.get_device_name(0)})")
else:
    _device = DEVICE

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"  Output dir : {OUTPUT_DIR}  ✓")

# ---------------------------------------------------------------------------
# 3. Load full AnnData and build balanced subset
# ---------------------------------------------------------------------------

print("\n[smoke_test] === Loading adata_benchmark.h5ad ===")
adata_full = sc.read_h5ad(H5AD_PATH)
print(f"  Full dataset: {adata_full.n_obs} cells × {adata_full.n_vars} genes")

# Confirm gene_col exists
if GENE_COL not in adata_full.var.columns:
    raise KeyError(
        f"gene_col={GENE_COL!r} not found in adata.var.columns.\n"
        f"Available columns: {list(adata_full.var.columns)}"
    )

# Build a balanced subset: CELLS_PER_TIMEPOINT cells sampled from each timepoint.
# Sampling is done with a fixed random seed for reproducibility.
print(f"\n[smoke_test] === Building balanced subset (~{CELLS_PER_TIMEPOINT} cells/timepoint) ===")

np.random.seed(42)
time_key   = "time_label"
timepoints = sorted(adata_full.obs[time_key].unique())
subset_idx = []

for tp in timepoints:
    tp_mask = adata_full.obs[time_key] == tp
    tp_indices = np.where(tp_mask)[0]
    n_sample   = min(CELLS_PER_TIMEPOINT, len(tp_indices))
    chosen     = np.random.choice(tp_indices, size=n_sample, replace=False)
    subset_idx.extend(chosen.tolist())
    print(f"  timepoint {tp:6.2f} : sampled {n_sample:3d} / {len(tp_indices):5d} cells")

subset_idx = sorted(subset_idx)
adata = adata_full[subset_idx].copy()
print(f"\n  Subset total: {adata.n_obs} cells × {adata.n_vars} genes")
print(f"  Timepoints covered: {len(adata.obs[time_key].unique())} / {len(timepoints)}")

# ---------------------------------------------------------------------------
# 4. Vocabulary overlap check
# ---------------------------------------------------------------------------

print("\n[smoke_test] === Vocabulary overlap check ===")

vocab_path = MODEL_DIR / "vocab.json"
with open(vocab_path, "r", encoding="utf-8") as f:
    vocab = json.load(f)

# Exclude special tokens from the gene count
special_tokens = {"<pad>", "<cls>", "<eoc>"}
vocab_genes    = set(vocab.keys()) - special_tokens

gene_symbols   = set(adata.var[GENE_COL])
overlap        = gene_symbols & vocab_genes
n_overlap      = len(overlap)
n_total        = len(gene_symbols)
overlap_frac   = n_overlap / n_total if n_total > 0 else 0.0

print(f"  Vocab gene count (excl. special tokens): {len(vocab_genes):,}")
print(f"  Genes in subset var[{GENE_COL!r}]:        {n_total}")
print(f"  Genes matching vocab:                   {n_overlap} / {n_total}  ({100*overlap_frac:.1f}%)")

if overlap_frac < 0.3:
    print(
        "  WARNING: Less than 30% of genes match the scGPT vocabulary.\n"
        "  Check that gene_col contains HGNC gene symbols (not Ensembl IDs)."
    )
elif overlap_frac < 0.6:
    print("  NOTE: Moderate overlap. Embeddings will be computed but token sequences are short.")
else:
    print("  Overlap is good. Embeddings should be informative.  ✓")

# ---------------------------------------------------------------------------
# 5. Run scGPT embedding
# ---------------------------------------------------------------------------

print(f"\n[smoke_test] === Running scGPT embedding (device={_device}, batch_size={BATCH_SIZE}) ===")
print("  Note: use_fast_transformer=False (flash_attn not installed in scgpt_env)")

t0 = time.time()

adata = scg.tasks.embed_data(
    adata,
    MODEL_DIR,
    gene_col         = GENE_COL,
    max_length       = 1200,
    batch_size       = BATCH_SIZE,
    device           = _device,
    use_fast_transformer = USE_FAST_TRANSFORMER,
    return_new_adata = False,   # add X_scGPT in-place on adata.obsm
)

embed_runtime = time.time() - t0
print(f"  Embedding done in {embed_runtime:.1f}s")

# Validate the embedding
if "X_scGPT" not in adata.obsm:
    raise RuntimeError(
        "X_scGPT not found in adata.obsm after embed_data — embedding failed."
    )

emb = adata.obsm["X_scGPT"]
emb_shape   = emb.shape
norms       = np.linalg.norm(emb, axis=1)
norm_mean   = float(norms.mean())
norm_std    = float(norms.std())
has_nan     = bool(np.isnan(emb).any())
has_inf     = bool(np.isinf(emb).any())

print(f"  Embedding shape : {emb_shape}")
print(f"  L2 norm (mean ± std): {norm_mean:.4f} ± {norm_std:.4f}  (should be ~1.0)")
print(f"  NaN present: {has_nan}  |  Inf present: {has_inf}")

if has_nan or has_inf:
    print("  WARNING: NaN or Inf in embeddings. Check model loading and input data.")

# ---------------------------------------------------------------------------
# 6. Downstream Scanpy steps: neighbors → UMAP → Leiden
# ---------------------------------------------------------------------------

print("\n[smoke_test] === Running neighbors / UMAP / Leiden ===")

sc.settings.verbosity = 1

sc.pp.neighbors(adata, use_rep="X_scGPT", n_neighbors=15, random_state=42)
print("  neighbors: done")

sc.tl.umap(adata, random_state=42)
print("  UMAP: done")

sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, random_state=42, key_added="leiden_scgpt")
n_leiden_clusters = adata.obs["leiden_scgpt"].nunique()
print(f"  Leiden (res={LEIDEN_RESOLUTION}): {n_leiden_clusters} clusters")

# ---------------------------------------------------------------------------
# 7. Save outputs
# ---------------------------------------------------------------------------

print("\n[smoke_test] === Saving outputs ===")

# 7a. Subset AnnData
h5ad_out = OUTPUT_DIR / "smoke_test_subset.h5ad"
adata.write_h5ad(h5ad_out)
print(f"  Saved subset h5ad   : {h5ad_out}")

# 7b. Metadata JSON
metadata = {
    "run_timestamp"        : time.strftime("%Y-%m-%dT%H:%M:%S"),
    "input_h5ad"           : str(H5AD_PATH),
    "model_dir"            : str(MODEL_DIR),
    "gene_col"             : GENE_COL,
    "device"               : _device,
    "use_fast_transformer" : USE_FAST_TRANSFORMER,
    "batch_size"           : BATCH_SIZE,
    "cells_per_timepoint"  : CELLS_PER_TIMEPOINT,
    "n_cells_subset"       : int(adata.n_obs),
    "n_genes_input"        : int(n_total),
    "n_genes_in_vocab"     : int(n_overlap),
    "vocab_overlap_frac"   : round(overlap_frac, 4),
    "embedding_shape"      : list(emb_shape),
    "embedding_norm_mean"  : round(norm_mean, 4),
    "embedding_norm_std"   : round(norm_std, 4),
    "embedding_has_nan"    : has_nan,
    "embedding_has_inf"    : has_inf,
    "leiden_resolution"    : LEIDEN_RESOLUTION,
    "n_leiden_clusters"    : int(n_leiden_clusters),
    "embed_runtime_sec"    : round(embed_runtime, 1),
}
meta_out = OUTPUT_DIR / "smoke_test_metadata.json"
with open(meta_out, "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)
print(f"  Saved metadata JSON : {meta_out}")

# 7c. UMAP figures
sc.settings.figdir = str(OUTPUT_DIR)

# Figure 1: coloured by timepoint
sc.pl.umap(
    adata,
    color    = time_key,
    title    = "scGPT smoke test — by timepoint (adata_benchmark subset)",
    palette  = "tab20",
    show     = False,
    save     = "_by_timepoint.png",
)
# scanpy appends the save string to "umap"; rename to something clearer
src = OUTPUT_DIR / "umap_by_timepoint.png"
if not src.exists():
    # scanpy saves as umap + save string — handle both naming conventions
    candidate = OUTPUT_DIR / "umap_by_timepoint.png"
    alt       = OUTPUT_DIR / f"umap{time_key}_by_timepoint.png"
    if not candidate.exists() and alt.exists():
        alt.rename(candidate)
print(f"  Saved UMAP (timepoint): {OUTPUT_DIR / 'umap_by_timepoint.png'}")

# Figure 2: coloured by Leiden cluster
sc.pl.umap(
    adata,
    color    = "leiden_scgpt",
    title    = f"scGPT smoke test — Leiden (res={LEIDEN_RESOLUTION})",
    show     = False,
    save     = "_by_leiden.png",
)
print(f"  Saved UMAP (Leiden)  : {OUTPUT_DIR / 'umap_by_leiden.png'}")

# ---------------------------------------------------------------------------
# 8. Final summary
# ---------------------------------------------------------------------------

print("\n[smoke_test] === Smoke test complete ===")
print(f"  Cells embedded  : {adata.n_obs}")
print(f"  Vocab overlap   : {n_overlap}/{n_total} ({100*overlap_frac:.1f}%)")
print(f"  Embedding shape : {emb_shape}")
print(f"  Norm mean       : {norm_mean:.4f}  (expected ≈ 1.0)")
print(f"  Leiden clusters : {n_leiden_clusters}")
print(f"  Embed runtime   : {embed_runtime:.1f}s")
print(f"\n  Outputs written to: {OUTPUT_DIR}")
print("  Files:")
for f in sorted(OUTPUT_DIR.iterdir()):
    size_kb = f.stat().st_size / 1024
    print(f"    {f.name:<40} {size_kb:>8.1f} KB")

print("\n[smoke_test] PASSED ✓" if not (has_nan or has_inf) else "\n[smoke_test] COMPLETED WITH WARNINGS ⚠")
