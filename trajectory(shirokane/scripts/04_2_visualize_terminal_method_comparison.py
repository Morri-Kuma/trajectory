#!/usr/bin/env python3
"""
05_visualize_terminal_only_methods_v2.py
=======================================
Visualize ONLY the terminal-cell universe and compare how the dynamic method
vs the biology-based TF method define iPSC cells.

Designed for outputs from 04_cellrank2_v12_gpt.py.

Figures
-------
1. terminal_only_overlap_umap.png
   Terminal cells only, colored as:
   - dynamic_only
   - bio_only
   - overlap_dyn_and_bio
   - terminal_neither

2. terminal_only_two_method_panels.png
   Left: dynamic method positive/negative in terminal cells
   Right: biological method categories in terminal cells

3. terminal_only_summary.png
   Count summary + overlap metrics
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


_CANDIDATES = [
    Path("/home/xzy0723/projects/trajectory"),
    Path(r"C:\Users\37620\trajectory"),
    Path("/sessions/cool-admiring-hawking/mnt/trajectory"),
    Path(__file__).resolve().parents[1],
]
PROJECT_ROOT = next((p for p in _CANDIDATES if p.exists()), Path("/home/xzy0723/projects/trajectory"))
FIGURES_DIR = PROJECT_ROOT / "results" / "figures"
METRICS_DIR = PROJECT_ROOT / "results" / "metrics"
CR2_DIR = PROJECT_ROOT / "results" / "cellrank2"


def find_latest_h5ad() -> Path:
    candidates = []
    for base in (METRICS_DIR, CR2_DIR):
        if base.exists():
            candidates.extend(base.glob("*_GSE230659_cellrank2.h5ad"))
    if not candidates:
        raise FileNotFoundError("Could not find '*_GSE230659_cellrank2.h5ad'. Pass --h5ad explicitly.")
    return sorted(candidates, key=lambda p: p.stat().st_mtime)[-1]


def require_columns(obs: pd.DataFrame, cols: Iterable[str]) -> None:
    missing = [c for c in cols if c not in obs.columns]
    if missing:
        raise ValueError(f"Missing required columns in adata.obs: {missing}")


def ensure_umap(adata: sc.AnnData) -> np.ndarray:
    if "X_umap" in adata.obsm:
        return np.asarray(adata.obsm["X_umap"])
    if "X_pca" in adata.obsm and adata.obsm["X_pca"].shape[1] >= 2:
        print("[WARN] X_umap not found. Using first 2 PCA dims.")
        return np.asarray(adata.obsm["X_pca"][:, :2])
    raise ValueError("Neither X_umap nor X_pca available.")


def scatter_by_category(ax, coords, labels, order, title, size=28.0, alpha=0.9):
    labels = pd.Series(labels).astype(str)
    plotted = 0
    for cat in order:
        mask = labels.values == cat
        if not np.any(mask):
            continue
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            s=size, alpha=alpha, label=f"{cat} (n={int(mask.sum())})",
            rasterized=True,
        )
        plotted += int(mask.sum())
    ax.set_title(title, fontsize=11)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(frameon=False, fontsize=8, loc="best")
    return plotted


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=str, default=None)
    parser.add_argument("--outdir", type=str, default=str(FIGURES_DIR))
    parser.add_argument("--prefix", type=str, default=None)
    args = parser.parse_args()

    h5ad_path = Path(args.h5ad) if args.h5ad else find_latest_h5ad()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Reading {h5ad_path}")
    adata = sc.read_h5ad(str(h5ad_path))
    require_columns(
        adata.obs,
        [
            "terminal_cell",
            "dyn_ips_terminal_cell",
            "bio_ips_terminal_cell",
            "bio_ips_category",
            "ips_final_mask",
        ],
    )

    coords = ensure_umap(adata)
    term = adata.obs["terminal_cell"].astype(bool).values
    dyn = adata.obs["dyn_ips_terminal_cell"].astype(bool).values
    bio = adata.obs["bio_ips_terminal_cell"].astype(bool).values

    prefix = args.prefix or h5ad_path.stem.replace("_GSE230659_cellrank2", "")

    # Restrict EVERYTHING to terminal cells only
    term_coords = coords[term]
    dyn_t = dyn[term]
    bio_t = bio[term]
    bio_cat = adata.obs.loc[term, "bio_ips_category"].astype(str)

    # Proper terminal-only overlap label
    term_overlap = np.where(
        dyn_t & bio_t, "overlap_dyn_and_bio",
        np.where(dyn_t & ~bio_t, "dynamic_only",
                 np.where(~dyn_t & bio_t, "bio_only", "terminal_neither"))
    )

    n_term = int(term.sum())
    n_dyn = int(dyn_t.sum())
    n_bio = int(bio_t.sum())
    n_overlap = int((dyn_t & bio_t).sum())
    n_dyn_only = int((dyn_t & ~bio_t).sum())
    n_bio_only = int((~dyn_t & bio_t).sum())
    n_neither = int((~dyn_t & ~bio_t).sum())
    union = n_dyn + n_bio - n_overlap
    jaccard = n_overlap / union if union > 0 else np.nan

    # Figure 1: terminal cells only overlap map
    fig1, ax1 = plt.subplots(figsize=(7.2, 6.0))
    scatter_by_category(
        ax1,
        term_coords,
        term_overlap,
        ["overlap_dyn_and_bio", "dynamic_only", "bio_only", "terminal_neither"],
        "Terminal cells only: dynamic vs biological definitions",
        size=42.0,
        alpha=0.95,
    )
    fig1.tight_layout()
    p1 = outdir / f"{prefix}_terminal_only_overlap_umap.png"
    fig1.savefig(p1, dpi=180, bbox_inches="tight")
    plt.close(fig1)
    print(f"[FIG] {p1}")

    # Figure 2: two method panels
    fig2, axes = plt.subplots(1, 2, figsize=(12.5, 5.4))
    scatter_by_category(
        axes[0],
        term_coords,
        np.where(dyn_t, "dynamic_iPSC", "dynamic_non_iPSC"),
        ["dynamic_iPSC", "dynamic_non_iPSC"],
        "Dynamic method within terminal cells",
        size=42.0,
        alpha=0.95,
    )
    scatter_by_category(
        axes[1],
        term_coords,
        bio_cat,
        ["POU5F1-only", "single_high", "bi_high", "tri_high", "tri_very_high", "terminal_non_ips_bio"],
        "Biological TF method within terminal cells",
        size=42.0,
        alpha=0.95,
    )
    fig2.tight_layout()
    p2 = outdir / f"{prefix}_terminal_only_two_method_panels.png"
    fig2.savefig(p2, dpi=180, bbox_inches="tight")
    plt.close(fig2)
    print(f"[FIG] {p2}")

    # Figure 3: summary plot
    fig3, axes3 = plt.subplots(1, 2, figsize=(12.5, 4.8))
    labels_left = ["dynamic", "biology", "overlap", "dyn_only", "bio_only", "neither"]
    vals_left = [n_dyn, n_bio, n_overlap, n_dyn_only, n_bio_only, n_neither]
    axes3[0].bar(labels_left, vals_left)
    axes3[0].set_title(f"Terminal-cell counts (n={n_term})")
    axes3[0].set_ylabel("Cell count")
    for i, v in enumerate(vals_left):
        axes3[0].text(i, v + 0.3, str(v), ha="center", va="bottom", fontsize=9)

    bio_counts = bio_cat.value_counts().reindex(
        ["POU5F1-only", "single_high", "bi_high", "tri_high", "tri_very_high", "terminal_non_ips_bio"],
        fill_value=0,
    )
    axes3[1].bar(bio_counts.index.tolist(), bio_counts.values.tolist())
    axes3[1].set_title("Biological method categories")
    axes3[1].set_ylabel("Cell count")
    axes3[1].tick_params(axis="x", rotation=35)
    for i, v in enumerate(bio_counts.values.tolist()):
        axes3[1].text(i, v + 0.3, str(int(v)), ha="center", va="bottom", fontsize=9)

    fig3.text(
        0.5, -0.03,
        f"Jaccard={jaccard:.4f} | overlap/dynamic={n_overlap / n_dyn if n_dyn else np.nan:.4f} | overlap/biology={n_overlap / n_bio if n_bio else np.nan:.4f}",
        ha="center", va="top", fontsize=10,
    )
    fig3.tight_layout()
    p3 = outdir / f"{prefix}_terminal_only_summary.png"
    fig3.savefig(p3, dpi=180, bbox_inches="tight")
    plt.close(fig3)
    print(f"[FIG] {p3}")

    # Summary table
    summary = pd.DataFrame(
        [
            ("n_terminal", n_term),
            ("n_dynamic", n_dyn),
            ("n_biological", n_bio),
            ("n_overlap", n_overlap),
            ("n_dynamic_only", n_dyn_only),
            ("n_bio_only", n_bio_only),
            ("n_neither", n_neither),
            ("jaccard", jaccard),
        ],
        columns=["metric", "value"],
    )
    p4 = outdir / f"{prefix}_terminal_only_summary.tsv"
    summary.to_csv(p4, sep="\t", index=False)
    print(f"[TSV] {p4}")

    print("\nDone.")


if __name__ == "__main__":
    main()
