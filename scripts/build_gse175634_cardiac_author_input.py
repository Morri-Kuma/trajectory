#!/usr/bin/env python3
"""
Build the GSE175634 cardiac-differentiation author-annotation HVG2000 benchmark input.

# ─────────────────────────────────────────────────────────────────────────────
# DATASET ID NOTE
# ─────────────────────────────────────────────────────────────────────────────
# The target GEO accession is GSE175634.
# The user originally requested "GSE174534" — that appears to be a typo
# (digits transposed: 174534 vs 175634). All source files in data/gse175634/
# are named GSE175634_*, confirming GSE175634 is the correct ID.
# This script and all downstream configs use GSE175634 as the canonical ID.
# ─────────────────────────────────────────────────────────────────────────────

Pipeline:
  1. Load raw-count MTX + gene/cell index TSVs from data/gse175634/.
  2. Load per-cell metadata (diffday, type, individual, etc.).
  3. QC filter cells (min/max genes, max mito pct, min counts).
  4. Normalise (library-size 1e4, log1p) and select top-2000 HVGs.
  5. Store raw counts in adata.layers["counts"].
  6. Compute abs_day from diffday string (day0→0, day1→1, …).
  7. Copy author cell-type label (type) to final_milestone_label_coarse,
     treating UNK cells as "uncertain" (retained but flagged for exclusion
     from official metrics).
  8. Add standard obs columns: cell_id, sample_id, time_label.
  9. Write HVG2000 h5ad to
       benchmark/inputs/gse175634_cardiac_author_hvg2000/
         GSE175634_cardiac_author_HVG2000_benchmark_input.h5ad

Silver-standard annotation step: SKIPPED.
  GSE175634 already has complete per-cell state labels from the authors
  (the `type` column). No marker-seed or trajectory-aware annotation pass
  is required.

Run from the repository root:
    python scripts/build_gse175634_cardiac_author_input.py [--overwrite]

Shirokane (runs ~10-20 min with 32 G RAM; the MTX files are ~4 GB each):
    python scripts/build_gse175634_cardiac_author_input.py
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import time as _time
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

# ─── Constants ────────────────────────────────────────────────────────────────

DATASET_ID = "GSE175634"
# NOTE: folder is gse175634 (lowercase); GEO accession is GSE175634.
DATA_FOLDER = "gse175634"
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")
T0 = _time.time()

N_TOP_GENES = 2000
NORM_TARGET_SUM = 1e4

QC_PARAMS = {
    "min_genes_per_cell": 200,
    "max_genes_per_cell": 10_000,
    "min_counts_per_cell": 500,
    "max_pct_mito": 25.0,
    "min_cells_per_gene": 3,
}

# Author cell-state type labels and their mapping to final_milestone_label_coarse
AUTHOR_TYPE_TO_MILESTONE: dict[str, str] = {
    "IPSC":  "IPSC",
    "MES":   "MES",
    "CMES":  "CMES",
    "PROG":  "PROG",
    "CM":    "CM",
    "CF":    "CF",
    "UNK":   "UNK",   # retained as uncertain; excluded from official metrics
}

# diffday string → numeric abs_day
DIFFDAY_TO_ABS_DAY: dict[str, float] = {
    "day0":  0.0,
    "day1":  1.0,
    "day3":  3.0,
    "day5":  5.0,
    "day7":  7.0,
    "day11": 11.0,
    "day15": 15.0,
}

# ─── Project-root detection ───────────────────────────────────────────────────

def find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
        raise FileNotFoundError(f"TRAJ_PROJECT_ROOT={env!r} does not exist.")
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "benchmark").is_dir() and (candidate / "scripts").is_dir():
            return candidate
    return here.parent


# ─── Sanity checks ────────────────────────────────────────────────────────────

def validate_inputs(data_dir: Path) -> None:
    """Raise FileNotFoundError if any required raw file is missing."""
    required = [
        "GSE175634_cell_counts.mtx",
        "gene_indices_counts.tsv",
        "GSE175634_cell_indices.tsv",
        "GSE175634_cell_metadata.tsv",
    ]
    missing = [f for f in required if not (data_dir / f).exists()]
    if missing:
        raise FileNotFoundError(
            f"Missing required files in {data_dir}:\n  " + "\n  ".join(missing)
        )
    print(f"[validate_inputs] all required raw files present in {data_dir}")


def validate_alignment(cell_index: pd.DataFrame, metadata: pd.DataFrame) -> None:
    """Check that cell indices and metadata have the same cells in the same order."""
    if len(cell_index) != len(metadata):
        raise ValueError(
            f"Cell index ({len(cell_index)}) and metadata ({len(metadata)}) "
            "have different row counts."
        )
    # cell names should match
    ci_names = cell_index["cell_name"].values
    md_names = metadata["cell"].values
    mismatches = (ci_names != md_names).sum()
    if mismatches > 0:
        raise ValueError(
            f"Cell name mismatch between index and metadata: {mismatches} / {len(ci_names)} differ."
        )
    print(f"[validate_alignment] cell index and metadata aligned ({len(cell_index)} cells).")


def validate_output(adata: ad.AnnData) -> None:
    """Lightweight QC checks on the built AnnData."""
    required_obs = ["cell_id", "sample_id", "time_label", "abs_day",
                    "final_milestone_label_coarse", "diffday", "individual", "type"]
    missing_obs = [c for c in required_obs if c not in adata.obs.columns]
    if missing_obs:
        raise ValueError(f"Output adata missing obs columns: {missing_obs}")

    if "counts" not in adata.layers:
        raise ValueError("Output adata missing layers['counts']")

    n_unk = (adata.obs["final_milestone_label_coarse"] == "UNK").sum()
    n_total = adata.n_obs
    label_counts = adata.obs["final_milestone_label_coarse"].value_counts().to_dict()

    # abs_day should cover all expected values
    expected_days = set(DIFFDAY_TO_ABS_DAY.values())
    observed_days = set(adata.obs["abs_day"].unique().tolist())
    missing_days = expected_days - observed_days
    if missing_days:
        print(f"  [WARN] Some expected abs_day values not present after QC: {missing_days}")

    print(f"  [OK] obs columns: {required_obs}")
    print(f"  [OK] layers['counts'] present")
    print(f"  [OK] n_obs={n_total}, n_vars={adata.n_vars}")
    print(f"  [OK] label counts: {label_counts}")
    print(f"  [OK] UNK cells: {n_unk} / {n_total} "
          f"({100 * n_unk / n_total:.1f}%) — flagged for exclusion from official metrics")
    print(f"  [OK] abs_day range: {adata.obs['abs_day'].min()} – {adata.obs['abs_day'].max()}")


# ─── Build helpers ────────────────────────────────────────────────────────────

def load_sparse_matrix(mtx_path: Path) -> sp.csr_matrix:
    print(f"  loading MTX: {mtx_path}  ({mtx_path.stat().st_size / 1e9:.2f} GB)")
    mat = sio.mmread(str(mtx_path))
    return sp.csr_matrix(mat.T)   # transpose: cells × genes


def qc_filter(adata: ad.AnnData, params: dict) -> ad.AnnData:
    """Apply QC filters and return filtered AnnData."""
    n_before = adata.n_obs
    print(f"  QC: {n_before} cells before filtering")

    # Mito genes (human ENSEMBL / gene-symbol prefix MT-)
    mito_mask = adata.var_names.str.startswith("MT-")
    adata.obs["pct_counts_mt"] = (
        np.array(adata[:, mito_mask].X.sum(axis=1)).flatten()
        / np.array(adata.X.sum(axis=1)).flatten().clip(1)
        * 100
    )
    adata.obs["n_genes"] = np.array((adata.X > 0).sum(axis=1)).flatten()
    adata.obs["n_counts"] = np.array(adata.X.sum(axis=1)).flatten()

    keep = (
        (adata.obs["n_genes"] >= params["min_genes_per_cell"])
        & (adata.obs["n_genes"] <= params["max_genes_per_cell"])
        & (adata.obs["n_counts"] >= params["min_counts_per_cell"])
        & (adata.obs["pct_counts_mt"] <= params["max_pct_mito"])
    )
    adata = adata[keep].copy()

    # Gene filter: min cells
    gene_cells = np.array((adata.X > 0).sum(axis=0)).flatten()
    adata = adata[:, gene_cells >= params["min_cells_per_gene"]].copy()

    print(f"  QC: {adata.n_obs} cells remaining after filtering "
          f"(removed {n_before - adata.n_obs}); {adata.n_vars} genes")
    return adata


def normalize_and_hvg(adata: ad.AnnData, n_top_genes: int) -> ad.AnnData:
    """Normalise counts, log1p, and select HVGs."""
    try:
        import scanpy as sc
    except ImportError:
        raise ImportError(
            "scanpy is required for normalisation and HVG selection.\n"
            "Install with: pip install scanpy"
        )

    # Stash raw counts before normalisation
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=NORM_TARGET_SUM)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"  HVG selection: {adata.n_vars} genes selected (n_top_genes={n_top_genes})")
    return adata


# ─── Main ─────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build GSE175634 cardiac author HVG2000 benchmark input."
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Allow overwriting existing output h5ad."
    )
    parser.add_argument(
        "--project-root", type=Path, default=None,
        help="Override project root directory (default: auto-detected)."
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    project_root = args.project_root or find_project_root()
    data_dir = project_root / "data" / DATA_FOLDER
    output_dir = (
        project_root / "benchmark" / "inputs" / "gse175634_cardiac_author_hvg2000"
    )
    output_h5ad = output_dir / "GSE175634_cardiac_author_HVG2000_benchmark_input.h5ad"
    qc_dir = output_dir / "build_qc"

    print("=" * 64)
    print(f"Build: GSE175634 cardiac author HVG2000 benchmark input")
    print(f"  Dataset ID  : {DATASET_ID}")
    print(f"  Project root: {project_root}")
    print(f"  Data dir    : {data_dir}")
    print(f"  Output h5ad : {output_h5ad}")
    print("=" * 64)

    # ── guard ──────────────────────────────────────────────────────────────────
    if output_h5ad.exists() and not args.overwrite:
        raise FileExistsError(
            f"{output_h5ad} already exists. Pass --overwrite to rebuild."
        )

    validate_inputs(data_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    # ── load gene index ────────────────────────────────────────────────────────
    print("[step 1] Loading gene index …")
    gene_idx = pd.read_csv(
        data_dir / "gene_indices_counts.tsv", sep="\t",
        dtype={"gene_index": int, "gene_name": str},
    )
    gene_names = gene_idx["gene_name"].values
    print(f"  {len(gene_names)} genes")

    # ── load cell index ────────────────────────────────────────────────────────
    print("[step 2] Loading cell index …")
    cell_idx = pd.read_csv(
        data_dir / "GSE175634_cell_indices.tsv", sep="\t",
        dtype={"cell_index": int, "cell_name": str},
    )
    print(f"  {len(cell_idx)} cells")

    # ── load metadata ──────────────────────────────────────────────────────────
    print("[step 3] Loading cell metadata …")
    metadata = pd.read_csv(
        data_dir / "GSE175634_cell_metadata.tsv", sep="\t",
        dtype=str,
    )
    print(f"  {len(metadata)} rows, columns: {list(metadata.columns)}")

    validate_alignment(cell_idx, metadata)

    # ── build AnnData ──────────────────────────────────────────────────────────
    print("[step 4] Loading expression matrix …")
    X_csr = load_sparse_matrix(data_dir / "GSE175634_cell_counts.mtx")
    print(f"  matrix shape: {X_csr.shape}  (cells × genes)")

    if X_csr.shape != (len(cell_idx), len(gene_names)):
        raise ValueError(
            f"Matrix shape {X_csr.shape} does not match "
            f"({len(cell_idx)} cells, {len(gene_names)} genes)"
        )

    obs = metadata.copy().reset_index(drop=True)
    obs.index = cell_idx["cell_name"].values

    # Parse abs_day from diffday
    obs["diffday_lower"] = obs["diffday"].str.lower().str.strip()
    obs["abs_day"] = obs["diffday_lower"].map(DIFFDAY_TO_ABS_DAY)
    unknown_days = obs.loc[obs["abs_day"].isna(), "diffday"].unique().tolist()
    if unknown_days:
        print(f"  [WARN] unmapped diffday values → abs_day NaN: {unknown_days}")

    # Standard annotation columns
    obs["cell_id"] = obs["cell"].astype(str)
    obs["sample_id"] = obs["sample"].astype(str)
    obs["time_label"] = obs["diffday"].astype(str)

    # Map author type → final_milestone_label_coarse
    obs["final_milestone_label_coarse"] = (
        obs["type"].map(AUTHOR_TYPE_TO_MILESTONE).fillna("UNK")
    )
    unmapped = obs.loc[
        ~obs["type"].isin(AUTHOR_TYPE_TO_MILESTONE), "type"
    ].unique().tolist()
    if unmapped:
        print(f"  [WARN] author type values not in mapping (treated as UNK): {unmapped}")

    var = pd.DataFrame({"gene_name": gene_names}, index=gene_names)

    adata = ad.AnnData(X=X_csr, obs=obs, var=var)
    adata.uns["dataset_id"] = DATASET_ID
    adata.uns["data_folder"] = DATA_FOLDER
    adata.uns["build_timestamp"] = TIMESTAMP
    adata.uns["silver_standard_skipped"] = (
        "GSE175634 uses author-provided cell-state annotations (type column). "
        "The silver-standard annotation generation step is not required."
    )
    print(f"  Built AnnData: {adata.shape}")
    gc.collect()

    # ── QC filter ──────────────────────────────────────────────────────────────
    print("[step 5] QC filtering …")
    adata = qc_filter(adata, QC_PARAMS)

    # ── normalise + HVG ───────────────────────────────────────────────────────
    print("[step 6] Normalising and selecting HVGs …")
    adata = normalize_and_hvg(adata, N_TOP_GENES)

    # ── QC reports ────────────────────────────────────────────────────────────
    print("[step 7] Writing build QC reports …")
    label_counts = adata.obs["final_milestone_label_coarse"].value_counts()
    label_counts.to_csv(qc_dir / "label_counts.csv", header=["n_cells"])

    day_counts = adata.obs["diffday"].value_counts().sort_index()
    day_counts.to_csv(qc_dir / "timepoint_counts.csv", header=["n_cells"])

    summary = {
        "dataset_id": DATASET_ID,
        "build_timestamp": TIMESTAMP,
        "n_cells_raw": int(len(cell_idx)),
        "n_cells_qc": int(adata.n_obs),
        "n_cells_removed": int(len(cell_idx) - adata.n_obs),
        "n_genes_raw": int(len(gene_names)),
        "n_genes_hvg": int(adata.n_vars),
        "qc_params": QC_PARAMS,
        "norm_target_sum": NORM_TARGET_SUM,
        "n_hvg": N_TOP_GENES,
        "label_counts": label_counts.to_dict(),
        "timepoint_counts": day_counts.to_dict(),
        "obs_columns": list(adata.obs.columns),
        "output_h5ad": str(output_h5ad),
        "silver_standard_skipped": True,
        "silver_standard_note": (
            "Author-provided type labels used directly as final_milestone_label_coarse. "
            "UNK cells are retained in the h5ad but will be excluded from official metrics."
        ),
        "runtime_s": round(_time.time() - T0, 2),
    }
    with (qc_dir / "build_summary.json").open("w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=True)

    # ── validate + write ──────────────────────────────────────────────────────
    print("[step 8] Validating output …")
    validate_output(adata)

    print(f"[step 9] Writing h5ad → {output_h5ad}")
    adata.write_h5ad(output_h5ad)

    elapsed = round(_time.time() - T0, 1)
    print()
    print(f"Done in {elapsed} s.")
    print(f"Output: {output_h5ad}")
    print()
    print("Next step: run scripts/build_gse175634_cardiac_author_provider.py")


if __name__ == "__main__":
    main()
