#!/usr/bin/env python
"""
Create a presentation-ready overview figure for GSE230659 scGPT annotation.

The figure explains two ideas:
1. How cell annotators such as scGPT process time-course single-cell data.
2. How annotation after pooled representation learning differs from simply
   grouping cells by experimental time point.

Inputs:
  benchmark/results/scgpt/full/adata_scgpt_annotated.h5ad

Outputs:
  benchmark/reports/figures/gse230659_scgpt_annotation_workflow.png
  benchmark/reports/figures/gse230659_scgpt_annotation_workflow.pdf
  benchmark/reports/figures/gse230659_timepoint_state_family_fraction.csv
"""

from __future__ import annotations

from pathlib import Path
import textwrap

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
INPUT_H5AD = ROOT / "benchmark" / "results" / "scgpt" / "full" / "adata_scgpt_annotated.h5ad"
FIG_DIR = ROOT / "benchmark" / "reports" / "figures"
OUT_PNG = FIG_DIR / "gse230659_scgpt_annotation_workflow.png"
OUT_PDF = FIG_DIR / "gse230659_scgpt_annotation_workflow.pdf"
OUT_FRAC = FIG_DIR / "gse230659_timepoint_state_family_fraction.csv"

STATE_ORDER = [
    "somatic",
    "stageII",
    "mid_transition",
    "stageIII_transition",
    "stageIII_peak",
    "late_transition",
    "late_branch",
    "pre_pluripotent",
    "pluripotent",
    "uncertain",
]

STATE_LABELS = {
    "somatic": "Somatic",
    "stageII": "StageII",
    "mid_transition": "Mid\ntransition",
    "stageIII_transition": "StageIII\ntransition",
    "stageIII_peak": "StageIII\npeak",
    "late_transition": "Late\ntransition",
    "late_branch": "Late\nbranch",
    "pre_pluripotent": "Pre-\npluripotent",
    "pluripotent": "Pluripotent",
    "uncertain": "Uncertain",
}

STATE_COLORS = {
    "somatic": "#2f6db5",
    "stageII": "#7fc97f",
    "mid_transition": "#2f8f4e",
    "stageIII_transition": "#f46d43",
    "stageIII_peak": "#d73027",
    "late_transition": "#984ea3",
    "late_branch": "#c51b7d",
    "pre_pluripotent": "#8c6bb1",
    "pluripotent": "#253494",
    "uncertain": "#bdbdbd",
}

TIME_CMAP = mpl.colormaps["turbo"]


def read_plot_table() -> pd.DataFrame:
    if not INPUT_H5AD.exists():
        raise FileNotFoundError(f"Missing input h5ad: {INPUT_H5AD}")

    adata = ad.read_h5ad(INPUT_H5AD, backed="r")
    if "X_umap" not in adata.obsm:
        adata.file.close()
        raise KeyError("Input h5ad does not contain obsm['X_umap'].")

    obs_cols = [
        "sample_id",
        "stage",
        "stage_day_label",
        "day_within_stage",
        "abs_day",
        "leiden_scgpt_res0.5",
        "scgpt_pseudostate_provisional",
        "scgpt_state_family",
        "scgpt_state_status",
    ]
    available = [c for c in obs_cols if c in adata.obs.columns]
    obs = adata.obs[available].copy()
    umap = np.asarray(adata.obsm["X_umap"])
    adata.file.close()

    obs["UMAP_1"] = umap[:, 0]
    obs["UMAP_2"] = umap[:, 1]
    obs["abs_day"] = pd.to_numeric(obs["abs_day"], errors="coerce")
    obs["time_label"] = obs.apply(format_time_label, axis=1)
    state_family = obs["scgpt_state_family"].astype(str)
    obs["state_family_label"] = state_family.map(STATE_LABELS).fillna(state_family)
    return obs


def format_time_label(row: pd.Series) -> str:
    stage = str(row.get("stage", ""))
    day = row.get("day_within_stage")
    if stage == "hCiPSCs":
        return "hCiPSCs"
    if pd.isna(day):
        return stage
    day_text = f"{float(day):g}"
    return f"{stage} D{day_text}"


def ordered_time_labels(df: pd.DataFrame) -> list[str]:
    return (
        df[["abs_day", "time_label"]]
        .drop_duplicates()
        .sort_values("abs_day")["time_label"]
        .tolist()
    )


def downsample_for_plot(df: pd.DataFrame, max_per_time: int = 2500) -> pd.DataFrame:
    parts = []
    for _, sub in df.groupby("time_label", sort=False):
        if len(sub) > max_per_time:
            parts.append(sub.sample(max_per_time, random_state=42))
        else:
            parts.append(sub)
    return pd.concat(parts, axis=0).sample(frac=1.0, random_state=42)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=16,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def style_umap(ax: plt.Axes) -> None:
    ax.set_xlabel("UMAP_1", fontsize=9)
    ax.set_ylabel("UMAP_2", fontsize=9)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color("#9a9a9a")


def plot_umap_by_time(ax: plt.Axes, df_plot: pd.DataFrame, time_labels: list[str]) -> None:
    time_to_idx = {label: i for i, label in enumerate(time_labels)}
    denom = max(len(time_labels) - 1, 1)
    colors = {label: TIME_CMAP(i / denom) for label, i in time_to_idx.items()}
    for label in time_labels:
        sub = df_plot[df_plot["time_label"] == label]
        ax.scatter(
            sub["UMAP_1"],
            sub["UMAP_2"],
            s=1.6,
            c=[colors[label]],
            linewidths=0,
            alpha=0.75,
            label=label,
            rasterized=True,
        )
    ax.set_title("Pooled GSE230659 human cells by time point", fontsize=11)
    style_umap(ax)
    leg = ax.legend(
        bbox_to_anchor=(0.5, -0.12),
        loc="upper center",
        frameon=False,
        markerscale=5,
        fontsize=5.8,
        handletextpad=0.2,
        borderaxespad=0,
        ncol=3,
    )
    for handle in leg.legend_handles:
        handle.set_alpha(1)


def plot_umap_by_family(ax: plt.Axes, df_plot: pd.DataFrame) -> None:
    for family in STATE_ORDER:
        sub = df_plot[df_plot["scgpt_state_family"] == family]
        if len(sub) == 0:
            continue
        ax.scatter(
            sub["UMAP_1"],
            sub["UMAP_2"],
            s=1.7,
            c=STATE_COLORS[family],
            linewidths=0,
            alpha=0.78,
            label=STATE_LABELS[family].replace("\n", " "),
            rasterized=True,
        )
    ax.set_title("scGPT-derived cell-state families", fontsize=11)
    style_umap(ax)
    leg = ax.legend(
        bbox_to_anchor=(0.5, -0.12),
        loc="upper center",
        frameon=False,
        markerscale=5,
        fontsize=5.8,
        handletextpad=0.25,
        borderaxespad=0,
        ncol=2,
    )
    for handle in leg.legend_handles:
        handle.set_alpha(1)


def state_fraction_table(df: pd.DataFrame) -> pd.DataFrame:
    counts = pd.crosstab(df["time_label"], df["scgpt_state_family"])
    counts = counts.reindex(index=ordered_time_labels(df), columns=STATE_ORDER, fill_value=0)
    frac = counts.div(counts.sum(axis=1), axis=0)
    frac.index.name = "timepoint"
    return frac


def plot_heatmap(ax: plt.Axes, frac: pd.DataFrame) -> None:
    sns.heatmap(
        frac,
        ax=ax,
        cmap="RdBu_r",
        vmin=0,
        vmax=max(0.65, float(frac.max().max())),
        linewidths=0.35,
        linecolor="white",
        cbar_kws={"label": "Fraction of cells", "shrink": 0.8},
    )
    ax.set_title("State composition after restoring labels to time points", fontsize=11)
    ax.set_xlabel("scGPT state family", fontsize=9)
    ax.set_ylabel("Observed time point", fontsize=9)
    ax.set_xticklabels([STATE_LABELS.get(x.get_text(), x.get_text()) for x in ax.get_xticklabels()], rotation=45, ha="right", fontsize=7)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=7)


def workflow_box(ax: plt.Axes, xy: tuple[float, float], w: float, h: float, title: str, body: str, color: str) -> None:
    box = FancyBboxPatch(
        xy,
        w,
        h,
        boxstyle="round,pad=0.018,rounding_size=0.02",
        linewidth=1.0,
        edgecolor=color,
        facecolor=mpl.colors.to_rgba(color, 0.08),
    )
    ax.add_patch(box)
    ax.text(xy[0] + w / 2, xy[1] + h - 0.08, title, ha="center", va="top", fontsize=10, fontweight="bold", color=color)
    ax.text(
        xy[0] + w / 2,
        xy[1] + h - 0.18,
        "\n".join(textwrap.wrap(body, width=24)),
        ha="center",
        va="top",
        fontsize=8,
        color="#333333",
        linespacing=1.25,
    )


def workflow_arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=14,
            linewidth=1.3,
            color="#555555",
        )
    )


def plot_workflow(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("scGPT annotation workflow", fontsize=10.5, loc="left", pad=8)

    steps = [
        ("Input", "15 observed GSE230659 human time points", "#2f6db5"),
        ("Pool all cells", "Cells from every time point are combined", "#2f8f4e"),
        ("scGPT annotation", "Embedding, neighbor graph, Leiden clustering, state labels", "#984ea3"),
        ("Restore labels", "Annotated cells are mapped back to their original time point", "#d73027"),
        ("Further analysis", "Use restored labels to compare state composition across time", "#444444"),
    ]
    y_positions = [0.78, 0.62, 0.46, 0.30, 0.12]
    box_x, box_w, box_h = 0.12, 0.76, 0.125
    for (title, body, color), y in zip(steps, y_positions):
        box = FancyBboxPatch(
            (box_x, y),
            box_w,
            box_h,
            boxstyle="round,pad=0.014,rounding_size=0.018",
            linewidth=1.0,
            edgecolor=color,
            facecolor=mpl.colors.to_rgba(color, 0.08),
        )
        ax.add_patch(box)
        ax.text(
            box_x + 0.04,
            y + box_h * 0.72,
            title,
            ha="left",
            va="center",
            fontsize=8.8,
            fontweight="bold",
            color=color,
        )
        ax.text(
            box_x + 0.04,
            y + box_h * 0.36,
            "\n".join(textwrap.wrap(body, width=42)),
            ha="left",
            va="center",
            fontsize=6.8,
            color="#333333",
            linespacing=1.15,
        )

    for y0, y1 in zip(y_positions[:-1], y_positions[1:]):
        workflow_arrow(ax, (0.50, y0 - 0.006), (0.50, y1 + box_h + 0.006))


def plot_time_vs_annotation_concept(ax: plt.Axes) -> None:
    ax.axis("off")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title("Automatic annotation vs time-point grouping", fontsize=10.5, loc="left", pad=8)

    left_x, right_x = 0.06, 0.56
    col_w = 0.36
    header_y = 0.83
    body_y = 0.22
    body_h = 0.48

    def column(x: float, title: str, body: list[str], color: str) -> None:
        ax.text(x, header_y, title, fontsize=9.5, fontweight="bold", color=color, ha="left", va="center")
        box = FancyBboxPatch(
            (x, body_y),
            col_w,
            body_h,
            boxstyle="round,pad=0.016,rounding_size=0.02",
            linewidth=1.1,
            edgecolor=color,
            facecolor=mpl.colors.to_rgba(color, 0.08),
        )
        ax.add_patch(box)
        y = body_y + body_h - 0.085
        for line in body:
            wrapped = "\n".join(textwrap.wrap(line, width=24))
            ax.text(x + 0.04, y, wrapped, fontsize=7.0, color="#333333", ha="left", va="top", linespacing=1.12)
            y -= 0.115 + 0.04 * wrapped.count("\n")

    column(
        left_x,
        "Time-point grouping",
        [
            "Groups cells by sample time",
            "Uses experimental labels",
            "May mix multiple states",
        ],
        "#2f6db5",
    )
    column(
        right_x,
        "Automatic annotation",
        [
            "Pools all cells first",
            "Learns expression states",
            "Clusters by similarity",
        ],
        "#984ea3",
    )


def make_figure() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.linewidth": 0.8,
        }
    )

    df = read_plot_table()
    time_labels = ordered_time_labels(df)
    df_plot = downsample_for_plot(df)
    frac = state_fraction_table(df)
    frac.to_csv(OUT_FRAC)

    fig = plt.figure(figsize=(17.0, 6.7), dpi=220)
    gs = fig.add_gridspec(
        nrows=1,
        ncols=4,
        width_ratios=[0.95, 1.05, 1.05, 1.34],
        wspace=0.62,
    )

    ax_d = fig.add_subplot(gs[0, 0])
    ax_a1 = fig.add_subplot(gs[0, 1])
    ax_a2 = fig.add_subplot(gs[0, 2])
    ax_c = fig.add_subplot(gs[0, 3])

    plot_workflow(ax_d)
    plot_umap_by_time(ax_a1, df_plot, time_labels)
    plot_umap_by_family(ax_a2, df_plot)
    plot_heatmap(ax_c, frac)

    add_panel_label(ax_d, "A")
    add_panel_label(ax_a1, "B")
    add_panel_label(ax_c, "C")

    fig.suptitle(
        "An example of annotated by realtime points vs. auto-annotation tool",
        fontsize=14,
        y=0.985,
    )
    fig.savefig(OUT_PNG, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT_PNG}")
    print(f"Wrote {OUT_PDF}")
    print(f"Wrote {OUT_FRAC}")


if __name__ == "__main__":
    make_figure()
