#!/usr/bin/env python3
"""
Plot original-vs-scNODE PCA scatter plots for every dataset with a scNODE
Scenario-A result.

For each dataset, PCA is fit on a combined sample of:
  - original benchmark-input HVG expression
  - scNODE projected_expression

The original and scNODE samples are then shown as two separate PC1/PC2 plots
using the same PCA basis and color map. This is intended as the primary visual
comparison for expression-space distribution differences.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from plot_scnode_vs_original_umaps import (
    color_map,
    discover_specs,
    find_project_root,
    load_original_sample,
    load_scnode_sample,
)


def fit_pca_projection(
    X_orig: np.ndarray,
    X_scn: np.ndarray,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    X = np.vstack([X_orig, X_scn]).astype(np.float32, copy=False)
    pca = PCA(n_components=2, svd_solver="randomized", random_state=seed)
    coords = pca.fit_transform(X)
    return coords[: len(X_orig)], coords[len(X_orig):], pca.explained_variance_ratio_


def plot_pca_panel(
    coords: np.ndarray,
    labels: np.ndarray,
    title: str,
    out_path: Path,
    label_to_color: dict[str, tuple],
    explained_variance: np.ndarray,
    point_size: float,
) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 7.2), dpi=180)
    labels = labels.astype(str)

    for label in sorted(label_to_color):
        mask = labels == label
        if not np.any(mask):
            continue
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=point_size,
            c=[label_to_color[label]],
            label=label,
            linewidths=0,
            alpha=0.72,
        )

    for label in sorted(label_to_color):
        mask = labels == label
        if not np.any(mask):
            continue
        center = np.median(coords[mask], axis=0)
        ax.text(
            center[0],
            center[1],
            label,
            fontsize=8,
            ha="center",
            va="center",
            color="black",
            bbox={
                "boxstyle": "round,pad=0.18",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.72,
            },
        )

    ax.set_title(title, fontsize=12)
    ax.set_xlabel(f"PC1 ({explained_variance[0] * 100:.1f}% var.)")
    ax.set_ylabel(f"PC2 ({explained_variance[1] * 100:.1f}% var.)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        markerscale=3.0,
        fontsize=7,
        title="cell label",
        title_fontsize=8,
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-original-cells", type=int, default=20000)
    parser.add_argument("--max-scnode-cells", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("benchmark/reports/visualizations/scnode_vs_original_pca"),
    )
    args = parser.parse_args()

    root = find_project_root()
    out_dir = args.output_dir if args.output_dir.is_absolute() else root / args.output_dir
    specs = discover_specs(root)
    if not specs:
        raise RuntimeError("No scNODE Scenario-A dataset specs found.")

    manifest = []
    for i, spec in enumerate(specs):
        print("=" * 80)
        print(f"[{i + 1}/{len(specs)}] {spec.dataset_id} ({spec.run_id})")
        print(f"  h5ad      : {spec.h5ad_path}")
        print(f"  scNODE dir: {spec.scnode_dir}")
        print(f"  label key : {spec.label_key}")

        X_orig, y_orig = load_original_sample(
            spec,
            max_cells=args.max_original_cells,
            seed=args.seed + i * 100,
        )
        X_scn, y_scn = load_scnode_sample(
            spec,
            max_cells=args.max_scnode_cells,
            seed=args.seed + i * 100,
        )
        print(f"  original sample: {X_orig.shape}")
        print(f"  scNODE sample  : {X_scn.shape}")

        pc_orig, pc_scn, explained = fit_pca_projection(
            X_orig,
            X_scn,
            seed=args.seed + i,
        )
        colors = color_map(np.concatenate([y_orig.astype(str), y_scn.astype(str)]))

        prefix = out_dir / spec.run_id
        original_png = prefix.with_name(prefix.name + "_original_cell_labels_pca.png")
        scnode_png = prefix.with_name(prefix.name + "_scnode_cell_labels_pca.png")

        plot_pca_panel(
            pc_orig,
            y_orig,
            f"{spec.dataset_id} original benchmark input",
            original_png,
            colors,
            explained,
            point_size=3.5,
        )
        plot_pca_panel(
            pc_scn,
            y_scn,
            f"{spec.dataset_id} scNODE inferred cells",
            scnode_png,
            colors,
            explained,
            point_size=4.0,
        )

        manifest.append(
            {
                "dataset_id": spec.dataset_id,
                "run_id": spec.run_id,
                "config_path": str(spec.config_path.relative_to(root)),
                "h5ad_path": str(spec.h5ad_path.relative_to(root)),
                "scnode_dir": str(spec.scnode_dir.relative_to(root)),
                "label_key": spec.label_key,
                "n_original_plotted": int(len(y_orig)),
                "n_scnode_plotted": int(len(y_scn)),
                "explained_variance_ratio_pc1": float(explained[0]),
                "explained_variance_ratio_pc2": float(explained[1]),
                "original_png": str(original_png.relative_to(root)),
                "scnode_png": str(scnode_png.relative_to(root)),
                "labels": sorted(pd.Series(np.concatenate([y_orig, y_scn])).unique().tolist()),
            }
        )
        print(f"  wrote: {original_png}")
        print(f"  wrote: {scnode_png}")

    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("=" * 80)
    print(f"Wrote manifest: {manifest_path}")


if __name__ == "__main__":
    main()
