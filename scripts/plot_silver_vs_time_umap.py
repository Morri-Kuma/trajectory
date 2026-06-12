#!/usr/bin/env python3
"""Plot time and silver-standard annotations in expression-derived UMAP space.

The saved X_umap coordinates in GSE178325 and GSE230659 were constructed from
X_scGPT, a representation also available to the silver annotation workflow.
This script intentionally computes a new UMAP from HVG2000 expression PCA so
the visualization does not reuse that stored representation.
"""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from umap import UMAP

from plot_silver_vs_time_pca import (
    DATASETS,
    STATE_COLORS,
    DatasetSpec,
    display_state_label,
    normalize_log1p,
    observed_state_order,
    padded_limits,
    value_counts,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT / "benchmark" / "reports" / "official_silver" / "umap_annotation_comparison"
)
TIME_KEY = "abs_day"
N_PCS = 50
N_NEIGHBORS = 30
MIN_DIST = 0.3
RANDOM_STATE = 0
UMAP_INIT = "random"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Destination directory for PNG, PDF, and JSON summary outputs.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=240,
        help="PNG rendering resolution (default: 240).",
    )
    parser.add_argument(
        "--datasets",
        nargs="*",
        default=[spec.dataset_id for spec in DATASETS],
        choices=[spec.dataset_id for spec in DATASETS],
        help="Dataset identifiers to plot (default: all three).",
    )
    return parser.parse_args()


def compute_expression_umap(
    adata: ad.AnnData, spec: DatasetSpec
) -> tuple[np.ndarray, np.ndarray]:
    matrix = adata.X
    if spec.normalize_log1p_for_pca:
        matrix = normalize_log1p(matrix)

    pca = PCA(n_components=N_PCS, svd_solver="arpack", random_state=RANDOM_STATE)
    pcs = pca.fit_transform(matrix).astype(np.float32, copy=False)
    reducer = UMAP(
        n_components=2,
        n_neighbors=N_NEIGHBORS,
        min_dist=MIN_DIST,
        metric="euclidean",
        init=UMAP_INIT,
        random_state=RANDOM_STATE,
        n_jobs=1,
        low_memory=True,
    )
    coordinates = reducer.fit_transform(pcs).astype(np.float32, copy=False)
    return coordinates, pca.explained_variance_ratio_


def format_day(day: float) -> str:
    if float(day).is_integer():
        return str(int(day))
    return f"{day:g}"


def plot_paired_umap(
    spec: DatasetSpec,
    coordinates: np.ndarray,
    time_values: np.ndarray,
    silver_values: np.ndarray,
    output_dir: Path,
    dpi: int,
) -> tuple[Path, Path]:
    n_cells = coordinates.shape[0]
    timepoints = np.sort(np.unique(time_values.astype(float)))
    states = observed_state_order(silver_values, spec.state_order)
    state_index = {state: index for index, state in enumerate(states)}
    silver_codes = np.asarray([state_index[str(value)] for value in silver_values])

    x_limits = padded_limits(coordinates[:, 0])
    y_limits = padded_limits(coordinates[:, 1])
    point_size = 1.8 if n_cells > 50000 else 2.4
    draw_order = np.random.default_rng(RANDOM_STATE).permutation(n_cells)

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(15.8, 6.7),
        sharex=True,
        sharey=True,
        constrained_layout=False,
    )
    fig.suptitle(
        f"{spec.dataset_id}: experimental time versus silver-standard annotation",
        fontsize=15,
        fontweight="bold",
        y=0.98,
    )
    fig.text(
        0.5,
        0.937,
        (
            f"Expression-derived UMAP; n = {n_cells:,} cells | "
            f"50-PC input, n_neighbors = {N_NEIGHBORS}, min_dist = {MIN_DIST}"
        ),
        ha="center",
        va="center",
        fontsize=10,
        color="#444444",
    )

    boundaries = np.r_[
        timepoints[0] - 0.5 * (timepoints[1] - timepoints[0])
        if len(timepoints) > 1
        else timepoints[0] - 0.5,
        (timepoints[:-1] + timepoints[1:]) / 2,
        timepoints[-1] + 0.5 * (timepoints[-1] - timepoints[-2])
        if len(timepoints) > 1
        else timepoints[0] + 0.5,
    ]
    time_cmap = mpl.colormaps["viridis"].resampled(len(timepoints))
    time_norm = mpl.colors.BoundaryNorm(boundaries, time_cmap.N)
    scatter_time = axes[0].scatter(
        coordinates[draw_order, 0],
        coordinates[draw_order, 1],
        c=time_values[draw_order].astype(float),
        s=point_size,
        cmap=time_cmap,
        norm=time_norm,
        alpha=0.62,
        linewidths=0,
        rasterized=True,
    )
    axes[0].set_title("Observed experimental time", fontsize=12, pad=12)
    colorbar = fig.colorbar(scatter_time, ax=axes[0], fraction=0.044, pad=0.025)
    colorbar.set_label("Day (abs_day)", fontsize=10)
    colorbar.set_ticks(timepoints)
    colorbar.set_ticklabels([format_day(day) for day in timepoints])
    colorbar.ax.tick_params(labelsize=8)

    state_colors = [STATE_COLORS.get(state, "#333333") for state in states]
    state_cmap = mpl.colors.ListedColormap(state_colors)
    state_norm = mpl.colors.BoundaryNorm(np.arange(len(states) + 1) - 0.5, len(states))
    axes[1].scatter(
        coordinates[draw_order, 0],
        coordinates[draw_order, 1],
        c=silver_codes[draw_order],
        s=point_size,
        cmap=state_cmap,
        norm=state_norm,
        alpha=0.62,
        linewidths=0,
        rasterized=True,
    )
    axes[1].set_title(spec.silver_panel_title, fontsize=12, pad=12)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="none",
            markerfacecolor=STATE_COLORS.get(state, "#333333"),
            markeredgecolor="none",
            markersize=7,
            label=display_state_label(state),
        )
        for state in states
    ]
    axes[1].legend(
        handles=handles,
        title="Silver label",
        loc="upper left",
        bbox_to_anchor=(1.01, 1.0),
        frameon=False,
        borderaxespad=0,
        fontsize=9,
        title_fontsize=10,
        handletextpad=0.4,
        labelspacing=0.62,
    )

    for axis in axes:
        axis.set_xlim(x_limits)
        axis.set_ylim(y_limits)
        axis.set_xlabel("UMAP1")
        axis.grid(False)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.set_aspect("equal", adjustable="box")
    axes[0].set_ylabel("UMAP2")
    fig.text(
        0.01,
        0.012,
        "All cells shown; ambiguous labels are retained where present. Stored X_scGPT/X_umap is not used.",
        fontsize=8.5,
        color="#555555",
    )
    fig.subplots_adjust(left=0.065, right=0.82, bottom=0.1, top=0.88, wspace=0.18)

    stem = f"{spec.dataset_id.lower()}_experimental_time_vs_silver_umap"
    png_path = output_dir / f"{stem}.png"
    pdf_path = output_dir / f"{stem}.pdf"
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png_path, pdf_path


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    selected = set(args.datasets)
    summaries: list[dict[str, object]] = []

    for spec in DATASETS:
        if spec.dataset_id not in selected:
            continue
        print(f"Reading {spec.dataset_id}: {spec.input_path}")
        adata = ad.read_h5ad(spec.input_path)
        missing_columns = {TIME_KEY, spec.silver_key} - set(adata.obs.columns)
        if missing_columns:
            raise KeyError(
                f"{spec.dataset_id} missing obs columns: {sorted(missing_columns)}"
            )

        time_values = adata.obs[TIME_KEY].to_numpy(dtype=float)
        silver_values = adata.obs[spec.silver_key].astype(str).to_numpy()
        coordinates, explained_variance = compute_expression_umap(adata, spec)
        png_path, pdf_path = plot_paired_umap(
            spec, coordinates, time_values, silver_values, output_dir, args.dpi
        )
        summaries.append(
            {
                "dataset_id": spec.dataset_id,
                "input_h5ad": str(spec.input_path),
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "time_key": TIME_KEY,
                "silver_key": spec.silver_key,
                "timepoint_counts": value_counts(time_values),
                "silver_label_counts": value_counts(silver_values),
                "expression_preprocessing": spec.pca_preprocessing,
                "embedding_method": (
                    "UMAP computed from 50-component PCA of the HVG2000 expression "
                    "matrix; stored X_scGPT and X_umap are not used."
                ),
                "pca_50_explained_variance_sum": float(explained_variance.sum()),
                "umap_parameters": {
                    "n_pcs": N_PCS,
                    "n_neighbors": N_NEIGHBORS,
                    "min_dist": MIN_DIST,
                    "metric": "euclidean",
                    "init": UMAP_INIT,
                    "random_state": RANDOM_STATE,
                },
                "png_output": str(png_path),
                "pdf_output": str(pdf_path),
            }
        )
        print(f"Wrote {png_path}")
        print(f"Wrote {pdf_path}")
        del coordinates, adata
        gc.collect()

    summary_path = output_dir / "umap_annotation_comparison_summary.json"
    summary_path.write_text(
        json.dumps({"datasets": summaries}, indent=2, ensure_ascii=True),
        encoding="utf-8",
    )
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
