"""Plotting helpers (matplotlib, headless Agg backend). Each returns the saved
path. UMAP is used when available; otherwise the first two PCs stand in.
"""
from __future__ import annotations

from typing import Optional

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from ..utils.io import ensure_dir  # noqa: E402


def _save(fig, path):
    ensure_dir(path.rsplit("/", 1)[0] if "/" in path else ".")
    fig.savefig(path, dpi=120, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_composition(df, path, title="Composition: marker-silver vs scANVI"):
    """df columns: time, state, frac_marker_silver, frac_scanvi."""
    times = sorted(df["time"].unique())
    states = sorted(df["state"].unique())
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    for ax, col, name in zip(axes, ["frac_marker_silver", "frac_scanvi"],
                             ["marker-silver", "scANVI"]):
        bottom = np.zeros(len(times))
        for s in states:
            vals = [float(df[(df.time == t) & (df.state == s)][col].sum()) for t in times]
            ax.bar(range(len(times)), vals, bottom=bottom, label=s)
            bottom += np.array(vals)
        ax.set_xticks(range(len(times)))
        ax.set_xticklabels(times, rotation=45, ha="right", fontsize=7)
        ax.set_title(name)
    axes[0].set_ylabel("fraction")
    axes[1].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
    fig.suptitle(title)
    return _save(fig, path)


def plot_agreement_heatmap(contingency, path, title="Annotation agreement"):
    fig, ax = plt.subplots(figsize=(max(5, contingency.shape[1]),
                                    max(4, contingency.shape[0] * 0.6)))
    mat = contingency.values.astype(float)
    mat_norm = mat / np.maximum(mat.sum(axis=1, keepdims=True), 1)
    im = ax.imshow(mat_norm, aspect="auto", cmap="viridis")
    ax.set_xticks(range(contingency.shape[1]))
    ax.set_xticklabels(contingency.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(contingency.shape[0]))
    ax.set_yticklabels(contingency.index, fontsize=7)
    ax.set_xlabel("scANVI label")
    ax.set_ylabel("marker-silver label")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="row-normalized")
    return _save(fig, path)


def plot_embedding(adata, color_key, path, title=None, rep: Optional[str] = None):
    """2D scatter coloured by an obs key. Uses obsm['X_umap'] if present, else
    first two PCs."""
    if rep and rep in adata.obsm:
        xy = adata.obsm[rep][:, :2]
    elif "X_umap" in adata.obsm:
        xy = adata.obsm["X_umap"][:, :2]
    else:
        xy = adata.obsm["X_pca"][:, :2]
    labels = adata.obs[color_key].astype(str).values
    fig, ax = plt.subplots(figsize=(6, 5))
    for s in sorted(np.unique(labels)):
        m = labels == s
        ax.scatter(xy[m, 0], xy[m, 1], s=6, alpha=0.7, label=s)
    ax.set_title(title or color_key)
    ax.set_xlabel("dim 1")
    ax.set_ylabel("dim 2")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=6, markerscale=2)
    return _save(fig, path)


def plot_pseudotime_scatter(pt_a, pt_b, path, title="Pseudotime: A vs B"):
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(pt_a, pt_b, s=6, alpha=0.5)
    ax.plot([0, 1], [0, 1], "r--", lw=1)
    ax.set_xlabel("pseudotime (marker-silver annotation)")
    ax.set_ylabel("pseudotime (scANVI annotation)")
    ax.set_title(title)
    return _save(fig, path)
