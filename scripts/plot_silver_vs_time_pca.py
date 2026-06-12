#!/usr/bin/env python3
"""Plot experimental-time and silver-standard annotations in shared PCA space.

Each dataset is shown as a paired panel: cells colored by observed experimental
time on the left and by its available transition-aware silver-standard label on
the right. GSE178325 and GSE230659 use the expanded supplementary annotation;
GSE242424 has one milestone label field and uses that field directly.
The PCA coordinates are computed once per dataset from its HVG2000 expression
matrix, so the two panels differ only in their annotation colors.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from sklearn.utils.sparsefuncs import inplace_row_scale


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT / "benchmark" / "reports" / "official_silver" / "pca_annotation_comparison"
)
TIME_KEY = "abs_day"


@dataclass(frozen=True)
class DatasetSpec:
    dataset_id: str
    input_path: Path
    pca_preprocessing: str
    normalize_log1p_for_pca: bool
    silver_key: str
    silver_panel_title: str
    state_order: tuple[str, ...]


DATASETS = (
    DatasetSpec(
        dataset_id="GSE178325",
        input_path=PROJECT_ROOT
        / "benchmark"
        / "inputs"
        / "gse178325_marker_fm_transition_silver_hvg2000"
        / "GSE178325_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad",
        pca_preprocessing=(
            "PCA computed from stored HVG2000 matrix; input metadata records "
            "log1p(total-count normalized expression)."
        ),
        normalize_log1p_for_pca=False,
        silver_key="final_milestone_label_expanded",
        silver_panel_title="Silver standard (expanded milestones)",
        state_order=(
            "stage_00_hADSCs",
            "stage_01_hADSCs_to_epithelial_like",
            "stage_02_epithelial_like",
            "stage_03_epithelial_like_to_intermediate_plastic",
            "stage_04_intermediate_plastic",
            "stage_05_intermediate_plastic_to_xen_like",
            "stage_06_xen_like",
            "stage_07_xen_like_to_hCiPS",
            "stage_08_hCiPS",
            "ambiguous",
        ),
    ),
    DatasetSpec(
        dataset_id="GSE230659",
        input_path=PROJECT_ROOT
        / "benchmark"
        / "inputs"
        / "gse230659_marker_fm_transition_silver_hvg2000"
        / "GSE230659_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad",
        pca_preprocessing=(
            "PCA computed after total-count normalization to 1e4 and log1p; "
            "the stored HVG2000 matrix is on count scale."
        ),
        normalize_log1p_for_pca=True,
        silver_key="final_milestone_label_expanded",
        silver_panel_title="Silver standard (expanded milestones)",
        state_order=(
            "stage_00_hADSCs",
            "stage_01_hADSCs_to_epithelial_like",
            "stage_02_epithelial_like",
            "stage_03_epithelial_like_to_intermediate_plastic",
            "stage_04_intermediate_plastic",
            "stage_05_intermediate_plastic_to_hCiPS",
            "stage_06_hCiPS",
            "ambiguous",
        ),
    ),
    DatasetSpec(
        dataset_id="GSE242424",
        input_path=PROJECT_ROOT
        / "benchmark"
        / "inputs"
        / "gse242424_author_cluster_matched"
        / "GSE242424_author_cluster_matched_HVG2000_benchmark_input.h5ad",
        pca_preprocessing=(
            "PCA computed from stored HVG2000 matrix; input metadata records "
            "log1p(total-count normalized expression)."
        ),
        normalize_log1p_for_pca=False,
        silver_key="final_milestone_label_coarse",
        silver_panel_title="Silver standard (milestones)",
        state_order=(
            "fibroblast",
            "fibroblast_like_stalled",
            "keratinocyte_like",
            "hOSK",
            "partial_intermediate",
            "partially_reprogrammed",
            "xOSK",
            "primary_intermediate",
            "pre_iPSC",
            "iPSC",
        ),
    ),
)


STATE_COLORS = {
    "hADSCs": "#4e79a7",
    "epithelial_like": "#f28e2b",
    "intermediate_plastic": "#59a14f",
    "xen_like": "#b07aa1",
    "hCiPS": "#e15759",
    "ambiguous": "#9d9d9d",
    "stage_00_hADSCs": "#4e79a7",
    "stage_01_hADSCs_to_epithelial_like": "#76b7b2",
    "stage_02_epithelial_like": "#f28e2b",
    "stage_03_epithelial_like_to_intermediate_plastic": "#edc948",
    "stage_04_intermediate_plastic": "#59a14f",
    "stage_05_intermediate_plastic_to_xen_like": "#d4a6c8",
    "stage_06_xen_like": "#b07aa1",
    "stage_07_xen_like_to_hCiPS": "#ff9d9a",
    "stage_08_hCiPS": "#e15759",
    "stage_05_intermediate_plastic_to_hCiPS": "#ff9d9a",
    "stage_06_hCiPS": "#e15759",
    "fibroblast": "#4e79a7",
    "fibroblast_like_stalled": "#bab0ac",
    "keratinocyte_like": "#edc948",
    "hOSK": "#76b7b2",
    "partial_intermediate": "#59a14f",
    "partially_reprogrammed": "#8cd17d",
    "xOSK": "#f28e2b",
    "primary_intermediate": "#ff9d9a",
    "pre_iPSC": "#af7aa1",
    "iPSC": "#e15759",
}


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


def normalize_log1p(matrix: sp.spmatrix | np.ndarray) -> sp.spmatrix | np.ndarray:
    """Return a total-count normalized, log1p transformed matrix."""
    if sp.issparse(matrix):
        transformed = matrix.tocsr(copy=True).astype(np.float32)
        totals = np.asarray(transformed.sum(axis=1)).ravel()
        scaling = np.divide(
            1e4,
            totals,
            out=np.zeros_like(totals, dtype=np.float32),
            where=totals > 0,
        )
        inplace_row_scale(transformed, scaling)
        transformed.data = np.log1p(transformed.data).astype(np.float32, copy=False)
        return transformed

    transformed = np.asarray(matrix, dtype=np.float32).copy()
    totals = transformed.sum(axis=1)
    scaling = np.divide(
        1e4,
        totals,
        out=np.zeros_like(totals, dtype=np.float32),
        where=totals > 0,
    )
    transformed *= scaling[:, None]
    np.log1p(transformed, out=transformed)
    return transformed


def compute_pca(adata: ad.AnnData, spec: DatasetSpec) -> tuple[np.ndarray, np.ndarray]:
    matrix = adata.X
    if spec.normalize_log1p_for_pca:
        matrix = normalize_log1p(matrix)
    model = PCA(n_components=2, svd_solver="arpack", random_state=0)
    coordinates = model.fit_transform(matrix).astype(np.float32, copy=False)
    return coordinates, model.explained_variance_ratio_


def format_day(day: float) -> str:
    if float(day).is_integer():
        return str(int(day))
    return f"{day:g}"


def observed_state_order(values: np.ndarray, preferred: tuple[str, ...]) -> list[str]:
    seen = set(values.astype(str))
    ordered = [state for state in preferred if state in seen]
    ordered.extend(sorted(seen - set(ordered)))
    return ordered


def display_state_label(state: str) -> str:
    if not state.startswith("stage_"):
        return state
    parts = state.split("_", 2)
    stage_number = parts[1]
    semantic_label = parts[2].replace("_to_", " -> ")
    return f"S{stage_number} {semantic_label}"


def padded_limits(values: np.ndarray) -> tuple[float, float]:
    low, high = np.quantile(values, [0.002, 0.998])
    span = high - low
    pad = max(span * 0.06, 0.25)
    return float(low - pad), float(high + pad)


def plot_paired_pca(
    spec: DatasetSpec,
    coordinates: np.ndarray,
    time_values: np.ndarray,
    silver_values: np.ndarray,
    explained_variance: np.ndarray,
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
    rng = np.random.default_rng(0)
    draw_order = rng.permutation(n_cells)

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
        f"HVG2000 expression PCA; n = {n_cells:,} cells | identical coordinates in both panels",
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
        axis.set_xlabel(f"PC1 ({explained_variance[0] * 100:.1f}% variance)")
        axis.grid(False)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.set_aspect("equal", adjustable="box")
    axes[0].set_ylabel(f"PC2 ({explained_variance[1] * 100:.1f}% variance)")
    fig.text(
        0.01,
        0.012,
        "All cells shown; ambiguous labels are retained where present.",
        fontsize=8.5,
        color="#555555",
    )
    fig.subplots_adjust(left=0.065, right=0.82, bottom=0.1, top=0.88, wspace=0.18)

    stem = f"{spec.dataset_id.lower()}_experimental_time_vs_silver_pca"
    png_path = output_dir / f"{stem}.png"
    pdf_path = output_dir / f"{stem}.pdf"
    fig.savefig(png_path, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png_path, pdf_path


def value_counts(values: np.ndarray) -> dict[str, int]:
    labels, counts = np.unique(values.astype(str), return_counts=True)
    return {str(label): int(count) for label, count in zip(labels, counts)}


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    selected = set(args.datasets)
    summaries: list[dict[str, object]] = []

    for spec in DATASETS:
        if spec.dataset_id not in selected:
            continue
        if not spec.input_path.exists():
            raise FileNotFoundError(f"Missing registered input: {spec.input_path}")

        print(f"Reading {spec.dataset_id}: {spec.input_path}")
        adata = ad.read_h5ad(spec.input_path)
        missing_columns = {TIME_KEY, spec.silver_key} - set(adata.obs.columns)
        if missing_columns:
            raise KeyError(
                f"{spec.dataset_id} missing obs columns: {sorted(missing_columns)}"
            )

        time_values = adata.obs[TIME_KEY].to_numpy(dtype=float)
        silver_values = adata.obs[spec.silver_key].astype(str).to_numpy()
        coordinates, explained_variance = compute_pca(adata, spec)
        png_path, pdf_path = plot_paired_pca(
            spec,
            coordinates,
            time_values,
            silver_values,
            explained_variance,
            output_dir,
            args.dpi,
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
                "pca_preprocessing": spec.pca_preprocessing,
                "pca_explained_variance_ratio": [
                    float(explained_variance[0]),
                    float(explained_variance[1]),
                ],
                "png_output": str(png_path),
                "pdf_output": str(pdf_path),
            }
        )
        print(f"Wrote {png_path}")
        print(f"Wrote {pdf_path}")

    summary_path = output_dir / "pca_annotation_comparison_summary.json"
    summary_path.write_text(
        json.dumps({"datasets": summaries}, indent=2, ensure_ascii=True),
        encoding="utf-8",
    )
    print(f"Wrote {summary_path}")


if __name__ == "__main__":
    main()
