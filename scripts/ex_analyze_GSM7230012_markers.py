#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
ex_analyze_GSM7230012_markers.py

Purpose
-------
Analyze single-cell POU5F1(OCT4), SOX2, and NANOG expression
within one hCiPSC sample: GSM7230012_hCiPSCs-0618.

What this script does
---------------------
1. Load the latest *_GSE230659_wot_fates.h5ad
2. Subset to one sample_id
3. Extract marker expression:
      - POU5F1 (fallback alias: OCT4)
      - SOX2
      - NANOG
4. Compute within-sample marker thresholds
5. Assign operational cell states:
      - tri_very_high
      - tri_high
      - bi_high
      - single_high
      - all_low
6. Build a local neighborhood graph + UMAP + Leiden clusters
7. Summarize marker patterns at:
      - cell level
      - cluster level
8. Save figures, tables, and an annotated h5ad

Important
---------
This is an exploratory state-stratification script INSIDE the hCiPSC sample.
It does NOT replace the paper's formal colony-level OCT4+ definition.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# =============================================================================
# CONFIG
# =============================================================================
PROJECT_ROOT = Path(r"C:\Users\37620\trajectory")
METRICS_DIR = PROJECT_ROOT / "results" / "metrics"
OUT_BASE = PROJECT_ROOT / "results" / "marker_states"

INPUT_H5AD = None
INPUT_PATTERN = "*_GSE230659_wot_fates.h5ad"

TARGET_SAMPLE = "GSM7230012_hCiPSCs-0618"
SAMPLE_ID_COL = "sample_id"

# Marker aliases: first existing one will be used
MARKER_ALIASES = {
    "POU5F1": ["POU5F1", "OCT4"],
    "SOX2": ["SOX2"],
    "NANOG": ["NANOG"],
}

# Threshold strategy:
# use non-zero expression distribution INSIDE THIS SAMPLE
HIGH_QUANTILE = 0.50       # operational "high"
VERY_HIGH_QUANTILE = 0.75  # operational "very high"

# Neighborhood / embedding
N_NEIGHBORS = 20
LEIDEN_RESOLUTION = 0.40
USE_EXISTING_PCA_IF_AVAILABLE = True
PCA_COMPONENTS = 30

# Figure settings
FIG_DPI = 180
POINT_SIZE = 12
ALPHA = 0.85

# If your h5ad uses .raw or a specific layer, change here
USE_RAW = False
EXPR_LAYER = None  # e.g. "log1p"

# =============================================================================
# HELPERS
# =============================================================================
def find_latest(directory: Path, pattern: str) -> Path:
    files = sorted(directory.glob(pattern), key=lambda p: p.stat().st_mtime)
    if not files:
        raise FileNotFoundError(f"No file matching '{pattern}' found in {directory}")
    return files[-1]


def choose_gene_name(adata: sc.AnnData, aliases: List[str]) -> str:
    for g in aliases:
        if g in adata.var_names:
            return g
    raise KeyError(f"None of these aliases were found in adata.var_names: {aliases}")


def extract_gene_vector(
    adata: sc.AnnData,
    gene_name: str,
    use_raw: bool = False,
    layer: str | None = None,
) -> np.ndarray:
    if use_raw:
        if adata.raw is None:
            raise ValueError("USE_RAW=True, but adata.raw is None.")
        if gene_name not in adata.raw.var_names:
            raise KeyError(f"{gene_name} not found in adata.raw.var_names.")
        x = adata.raw[:, gene_name].X
    elif layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers.")
        idx = adata.var_names.get_loc(gene_name)
        x = adata.layers[layer][:, idx]
    else:
        idx = adata.var_names.get_loc(gene_name)
        x = adata.X[:, idx]

    if sp.issparse(x):
        x = x.toarray().ravel()
    else:
        x = np.asarray(x).ravel()

    return x.astype(float)


def safe_zscore(x: np.ndarray) -> np.ndarray:
    mu = np.nanmean(x)
    sd = np.nanstd(x)
    if sd < 1e-12:
        return np.zeros_like(x, dtype=float)
    return (x - mu) / sd


def nonzero_quantile_threshold(x: np.ndarray, q: float) -> float:
    nz = x[x > 0]
    if nz.size == 0:
        return 0.0
    return float(np.quantile(nz, q))


def assign_state(n_high: int, n_very_high: int) -> str:
    if n_very_high == 3:
        return "tri_very_high"
    if n_high == 3:
        return "tri_high"
    if n_high == 2:
        return "bi_high"
    if n_high == 1:
        return "single_high"
    return "all_low"


def summarize_marker(vec: np.ndarray) -> Dict[str, float]:
    return {
        "n_cells": int(vec.size),
        "n_nonzero": int((vec > 0).sum()),
        "pct_nonzero": float((vec > 0).mean() * 100.0),
        "min": float(np.min(vec)),
        "q25": float(np.quantile(vec, 0.25)),
        "median": float(np.quantile(vec, 0.50)),
        "q75": float(np.quantile(vec, 0.75)),
        "max": float(np.max(vec)),
        "mean": float(np.mean(vec)),
    }


def ensure_embedding_and_clusters(adata: sc.AnnData) -> None:
    # Use existing PCA if available
    if USE_EXISTING_PCA_IF_AVAILABLE and "X_pca" in adata.obsm:
        pass
    else:
        sc.pp.pca(adata, n_comps=min(PCA_COMPONENTS, adata.n_vars))

    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X_pca")
    sc.tl.umap(adata)

    try:
        sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, key_added="leiden_markers")
    except Exception:
        # fallback: one cluster only
        adata.obs["leiden_markers"] = pd.Categorical(["0"] * adata.n_obs)


def make_output_dirs(base: Path, target_sample: str) -> Tuple[Path, Path]:
    safe_name = target_sample.replace("/", "_").replace("\\", "_")
    out_dir = base / safe_name
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)
    return out_dir, fig_dir


def save_histogram(x: np.ndarray, title: str, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.hist(x, bins=40)
    ax.set_title(title)
    ax.set_xlabel("expression")
    ax.set_ylabel("cell count")
    fig.tight_layout()
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def save_scatter(df: pd.DataFrame, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(6, 5))
    sca = ax.scatter(
        df["expr_POU5F1"].values,
        df["expr_SOX2"].values,
        c=df["expr_NANOG"].values,
        s=10,
        alpha=0.8,
    )
    ax.set_xlabel("POU5F1/OCT4")
    ax.set_ylabel("SOX2")
    ax.set_title("POU5F1 vs SOX2, colored by NANOG")
    cbar = fig.colorbar(sca, ax=ax)
    cbar.set_label("NANOG")
    fig.tight_layout()
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


def save_umap_color(adata: sc.AnnData, color: str, title: str, out_png: Path) -> None:
    fig = sc.pl.umap(
        adata,
        color=color,
        return_fig=True,
        show=False,
        title=title,
        size=POINT_SIZE,
        alpha=ALPHA,
    )
    fig.savefig(out_png, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================
def main() -> None:
    out_dir, fig_dir = make_output_dirs(OUT_BASE, TARGET_SAMPLE)

    if INPUT_H5AD is None:
        h5ad_path = find_latest(METRICS_DIR, INPUT_PATTERN)
    else:
        h5ad_path = Path(INPUT_H5AD)

    print("=" * 80)
    print("LOAD H5AD")
    print("=" * 80)
    print(f"Input h5ad: {h5ad_path}")

    adata = sc.read_h5ad(str(h5ad_path))
    print(f"adata shape: {adata.shape}")

    if SAMPLE_ID_COL not in adata.obs.columns:
        raise KeyError(
            f"'{SAMPLE_ID_COL}' not in adata.obs.columns. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    unique_samples = adata.obs[SAMPLE_ID_COL].astype(str).unique().tolist()
    if TARGET_SAMPLE not in unique_samples:
        raise ValueError(
            f"TARGET_SAMPLE='{TARGET_SAMPLE}' not found in adata.obs['{SAMPLE_ID_COL}'].\n"
            f"First 20 sample IDs:\n{unique_samples[:20]}"
        )

    ad = adata[adata.obs[SAMPLE_ID_COL].astype(str) == TARGET_SAMPLE].copy()

    print("=" * 80)
    print("SUBSET SAMPLE")
    print("=" * 80)
    print(f"Target sample : {TARGET_SAMPLE}")
    print(f"Subset shape  : {ad.shape}")

    if ad.n_obs == 0:
        raise ValueError("No cells left after subsetting.")

    # choose actual gene names
    gene_map = {
        canonical: choose_gene_name(ad, aliases)
        for canonical, aliases in MARKER_ALIASES.items()
    }
    print("Marker mapping:")
    for k, v in gene_map.items():
        print(f"  {k} -> {v}")

    # expression
    expr_pou5f1 = extract_gene_vector(ad, gene_map["POU5F1"], USE_RAW, EXPR_LAYER)
    expr_sox2 = extract_gene_vector(ad, gene_map["SOX2"], USE_RAW, EXPR_LAYER)
    expr_nanog = extract_gene_vector(ad, gene_map["NANOG"], USE_RAW, EXPR_LAYER)

    # thresholds inside this sample
    thr = {
        "POU5F1_high": nonzero_quantile_threshold(expr_pou5f1, HIGH_QUANTILE),
        "SOX2_high": nonzero_quantile_threshold(expr_sox2, HIGH_QUANTILE),
        "NANOG_high": nonzero_quantile_threshold(expr_nanog, HIGH_QUANTILE),
        "POU5F1_very_high": nonzero_quantile_threshold(expr_pou5f1, VERY_HIGH_QUANTILE),
        "SOX2_very_high": nonzero_quantile_threshold(expr_sox2, VERY_HIGH_QUANTILE),
        "NANOG_very_high": nonzero_quantile_threshold(expr_nanog, VERY_HIGH_QUANTILE),
    }

    # cell-level states
    pou5f1_high = expr_pou5f1 >= thr["POU5F1_high"]
    sox2_high = expr_sox2 >= thr["SOX2_high"]
    nanog_high = expr_nanog >= thr["NANOG_high"]

    pou5f1_vhigh = expr_pou5f1 >= thr["POU5F1_very_high"]
    sox2_vhigh = expr_sox2 >= thr["SOX2_very_high"]
    nanog_vhigh = expr_nanog >= thr["NANOG_very_high"]

    n_high = pou5f1_high.astype(int) + sox2_high.astype(int) + nanog_high.astype(int)
    n_vhigh = pou5f1_vhigh.astype(int) + sox2_vhigh.astype(int) + nanog_vhigh.astype(int)

    z_pou5f1 = safe_zscore(expr_pou5f1)
    z_sox2 = safe_zscore(expr_sox2)
    z_nanog = safe_zscore(expr_nanog)
    pluri_3gene_score = np.mean(np.vstack([z_pou5f1, z_sox2, z_nanog]), axis=0)

    cell_state = [assign_state(int(a), int(b)) for a, b in zip(n_high, n_vhigh)]

    # save to obs
    ad.obs["expr_POU5F1"] = expr_pou5f1
    ad.obs["expr_SOX2"] = expr_sox2
    ad.obs["expr_NANOG"] = expr_nanog
    ad.obs["z_POU5F1"] = z_pou5f1
    ad.obs["z_SOX2"] = z_sox2
    ad.obs["z_NANOG"] = z_nanog
    ad.obs["pluri_3gene_score"] = pluri_3gene_score
    ad.obs["n_markers_high"] = n_high
    ad.obs["n_markers_very_high"] = n_vhigh
    ad.obs["marker_state"] = pd.Categorical(cell_state)

    # local embedding + clusters
    print("=" * 80)
    print("EMBEDDING / CLUSTERING")
    print("=" * 80)
    ensure_embedding_and_clusters(ad)

    # cell-level table
    df_cells = ad.obs.copy()
    df_cells.insert(0, "cell_id", ad.obs_names.astype(str))
    cell_out = out_dir / "cell_level_marker_expression.tsv"
    df_cells.to_csv(cell_out, sep="\t", index=False)

    # threshold summary
    threshold_rows = []
    for k, v in thr.items():
        threshold_rows.append({"parameter": k, "value": float(v)})
    threshold_df = pd.DataFrame(threshold_rows)
    threshold_out = out_dir / "marker_thresholds.tsv"
    threshold_df.to_csv(threshold_out, sep="\t", index=False)

    # marker summary
    marker_summary = pd.DataFrame(
        [
            {"marker": "POU5F1", **summarize_marker(expr_pou5f1)},
            {"marker": "SOX2", **summarize_marker(expr_sox2)},
            {"marker": "NANOG", **summarize_marker(expr_nanog)},
        ]
    )
    marker_summary_out = out_dir / "marker_summary.tsv"
    marker_summary.to_csv(marker_summary_out, sep="\t", index=False)

    # cluster summary
    cluster_key = "leiden_markers"
    cluster_rows = []
    for cl in ad.obs[cluster_key].astype(str).unique():
        idx = ad.obs[cluster_key].astype(str) == cl
        sub = ad.obs.loc[idx]
        state_counts = sub["marker_state"].value_counts().to_dict()

        cluster_rows.append(
            {
                "cluster": cl,
                "n_cells": int(idx.sum()),
                "mean_POU5F1": float(sub["expr_POU5F1"].mean()),
                "mean_SOX2": float(sub["expr_SOX2"].mean()),
                "mean_NANOG": float(sub["expr_NANOG"].mean()),
                "mean_pluri_3gene_score": float(sub["pluri_3gene_score"].mean()),
                "pct_tri_very_high": float((sub["marker_state"] == "tri_very_high").mean() * 100),
                "pct_tri_high": float((sub["marker_state"] == "tri_high").mean() * 100),
                "pct_bi_high": float((sub["marker_state"] == "bi_high").mean() * 100),
                "pct_single_high": float((sub["marker_state"] == "single_high").mean() * 100),
                "pct_all_low": float((sub["marker_state"] == "all_low").mean() * 100),
                "dominant_state": max(state_counts, key=state_counts.get),
            }
        )

    cluster_df = pd.DataFrame(cluster_rows).sort_values(
        ["mean_pluri_3gene_score", "n_cells"], ascending=[False, False]
    )
    cluster_out = out_dir / "cluster_level_marker_summary.tsv"
    cluster_df.to_csv(cluster_out, sep="\t", index=False)

    # save annotated h5ad
    h5ad_out = out_dir / "GSM7230012_hCiPSCs-0618_marker_annotated.h5ad"
    ad.write_h5ad(h5ad_out)

    # figures
    print("=" * 80)
    print("SAVE FIGURES")
    print("=" * 80)

    save_histogram(expr_pou5f1, "POU5F1/OCT4 expression", fig_dir / "hist_POU5F1.png")
    save_histogram(expr_sox2, "SOX2 expression", fig_dir / "hist_SOX2.png")
    save_histogram(expr_nanog, "NANOG expression", fig_dir / "hist_NANOG.png")

    save_scatter(df_cells, fig_dir / "scatter_POU5F1_vs_SOX2_color_NANOG.png")

    save_umap_color(ad, "expr_POU5F1", "UMAP: POU5F1/OCT4", fig_dir / "umap_POU5F1.png")
    save_umap_color(ad, "expr_SOX2", "UMAP: SOX2", fig_dir / "umap_SOX2.png")
    save_umap_color(ad, "expr_NANOG", "UMAP: NANOG", fig_dir / "umap_NANOG.png")
    save_umap_color(ad, "pluri_3gene_score", "UMAP: 3-gene pluripotency score", fig_dir / "umap_pluri_3gene_score.png")
    save_umap_color(ad, "marker_state", "UMAP: marker_state", fig_dir / "umap_marker_state.png")
    save_umap_color(ad, "leiden_markers", "UMAP: leiden_markers", fig_dir / "umap_leiden_markers.png")

    # console summary
    print("=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print("Thresholds:")
    print(threshold_df.to_string(index=False))

    print("\nMarker-state counts:")
    print(ad.obs["marker_state"].value_counts().to_string())

    print("\nTop clusters by mean_pluri_3gene_score:")
    print(cluster_df.head(10).to_string(index=False))

    print("\nSaved files:")
    print(f"  {cell_out}")
    print(f"  {threshold_out}")
    print(f"  {marker_summary_out}")
    print(f"  {cluster_out}")
    print(f"  {h5ad_out}")
    print(f"  figures -> {fig_dir}")


if __name__ == "__main__":
    main()