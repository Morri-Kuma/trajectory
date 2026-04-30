#!/usr/bin/env python
"""
Build presentation-ready stage report assets from benchmark summaries.

Inputs:
  benchmark/reports/core_summary.csv

Outputs:
  benchmark/reports/formal_benchmark_summary.csv
  benchmark/reports/figures/lineage_auroc_by_model_scenario.png
  benchmark/reports/scnode_formal_metrics_table.md
  benchmark/reports/projection_formal_metrics_table.md
  benchmark/reports/prescient_reduced_validation_metrics_table.md
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
REPORTS = ROOT / "benchmark" / "reports"
FIG_DIR = REPORTS / "figures"


def normalize_scenario(value: object) -> str | None:
    match = re.match(r"([ABC])_", str(value))
    return match.group(1) if match else None


def load_presentation_rows() -> pd.DataFrame:
    core = pd.read_csv(REPORTS / "core_summary.csv")
    core["scenario_normalized"] = core["scenario"].map(normalize_scenario)

    rows = core[
        (core["result_class"] == "official")
        & (core["scenario_normalized"].isin(["A", "B", "C"]))
        & (
            core["method"].isin(["wot", "cellrank2"])
            | (
                core["method"].isin(["scnode", "prescient"])
                & core["scenario"].str.contains("hvg2000_formal", na=False)
            )
        )
    ].copy()

    method_label = {
        "wot": "WOT",
        "cellrank2": "CellRank2",
        "scnode": "scNODE",
        "prescient": "PRESCIENT",
    }
    rows["method_label"] = rows["method"].map(method_label).fillna(rows["method"])
    rows["evaluation_scope"] = np.where(
        rows["method"].isin(["wot", "cellrank2"]),
        "Lineage Fidelity only (method capability)",
        "Forecast Accuracy + Embedding Coherence + Lineage Fidelity",
    )
    return rows.sort_values(["scenario_normalized", "method_label"])


def write_formal_summary(rows: pd.DataFrame) -> Path:
    columns = [
        "method_label",
        "method",
        "scenario_normalized",
        "scenario",
        "result_class",
        "ground_truth_provider",
        "forecast_status",
        "forecast_wasserstein",
        "forecast_gaussian_mmd",
        "forecast_energy_distance_mmd",
        "forecast_hausdorff_loss",
        "embedding_status",
        "embedding_ari",
        "embedding_entropy",
        "lineage_status",
        "lineage_auroc",
        "lineage_auprc",
        "lineage_jaccard",
        "lineage_single_step_recovery",
        "lineage_multi_step_recovery",
        "baseline_auroc",
        "baseline_auprc",
        "evaluation_scope",
        "source",
    ]
    out = REPORTS / "formal_benchmark_summary.csv"
    rows[columns].to_csv(out, index=False)
    return out


def write_scnode_table(rows: pd.DataFrame) -> Path:
    scnode = rows[rows["method"] == "scnode"].copy()
    columns = [
        "scenario_normalized",
        "forecast_wasserstein",
        "forecast_gaussian_mmd",
        "forecast_energy_distance_mmd",
        "forecast_hausdorff_loss",
        "embedding_ari",
        "embedding_entropy",
        "lineage_auroc",
        "lineage_auprc",
        "baseline_auroc",
    ]
    scnode = scnode[columns].sort_values("scenario_normalized")
    for column in columns[1:]:
        scnode[column] = pd.to_numeric(scnode[column], errors="coerce")

    out = REPORTS / "scnode_formal_metrics_table.md"
    lines = [
        "# scNODE Formal Metrics (HVG2000, scGPT-v1)",
        "",
        "| Scenario | Forecast WD (lower) | Gaussian MMD (lower) | Energy MMD (lower) | Hausdorff (lower) | Embedding ARI (higher) | Entropy (lower) | Lineage AUROC (higher) | Lineage AUPRC (higher) | Baseline AUROC (higher) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in scnode.iterrows():
        lines.append(
            f"| {row['scenario_normalized']} | "
            f"{row['forecast_wasserstein']:.2f} | "
            f"{row['forecast_gaussian_mmd']:.4f} | "
            f"{row['forecast_energy_distance_mmd']:.2f} | "
            f"{row['forecast_hausdorff_loss']:.2f} | "
            f"{row['embedding_ari']:.4f} | "
            f"{row['embedding_entropy']:.4f} | "
            f"{row['lineage_auroc']:.4f} | "
            f"{row['lineage_auprc']:.4f} | "
            f"{row['baseline_auroc']:.4f} |"
        )
    lines.extend(
        [
            "",
            "Note: lower is better for Forecast WD, MMD, Hausdorff, and entropy; higher is better for ARI, AUROC, and AUPRC.",
            "WOT and CellRank2 are not included in this table because they are lineage-only under the current benchmark capability rules.",
        ]
    )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def write_projection_formal_table(rows: pd.DataFrame) -> Path:
    projection = rows[rows["method"].isin(["scnode", "prescient"])].copy()
    columns = [
        "method_label",
        "scenario_normalized",
        "forecast_wasserstein",
        "forecast_gaussian_mmd",
        "forecast_energy_distance_mmd",
        "forecast_hausdorff_loss",
        "embedding_ari",
        "embedding_entropy",
        "lineage_auroc",
        "lineage_auprc",
        "baseline_auroc",
    ]
    projection = projection[columns].sort_values(["scenario_normalized", "method_label"])
    for column in columns[2:]:
        projection[column] = pd.to_numeric(projection[column], errors="coerce")

    out = REPORTS / "projection_formal_metrics_table.md"
    lines = [
        "# Projection-Capable Formal Metrics (HVG2000, scGPT-v1)",
        "",
        "| Scenario | Method | Forecast WD (lower) | Gaussian MMD (lower) | Energy MMD (lower) | Hausdorff (lower) | Embedding ARI (higher) | Entropy (lower) | Lineage AUROC (higher) | Lineage AUPRC (higher) | Baseline AUROC (higher) |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in projection.iterrows():
        lines.append(
            f"| {row['scenario_normalized']} | "
            f"{row['method_label']} | "
            f"{row['forecast_wasserstein']:.4f} | "
            f"{row['forecast_gaussian_mmd']:.4f} | "
            f"{row['forecast_energy_distance_mmd']:.2f} | "
            f"{row['forecast_hausdorff_loss']:.2f} | "
            f"{row['embedding_ari']:.4f} | "
            f"{row['embedding_entropy']:.4f} | "
            f"{row['lineage_auroc']:.4f} | "
            f"{row['lineage_auprc']:.4f} | "
            f"{row['baseline_auroc']:.4f} |"
        )
    lines.extend(
        [
            "",
            "Note: lower is better for Forecast WD, MMD, Hausdorff, and entropy; higher is better for ARI, AUROC, and AUPRC.",
            "This table includes only official projection-capable formal runs under the shared HVG2000 scGPT-v1 input.",
        ]
    )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def write_prescient_reduced_table() -> Path:
    core = pd.read_csv(REPORTS / "core_summary.csv")
    rows = core[
        (core["method"] == "prescient")
        & (core["result_class"] == "reduced_validation")
        & (core["scenario"].str.contains("hvg2000_cpu", na=False))
    ].copy()
    rows["scenario_normalized"] = rows["scenario"].map(normalize_scenario)
    rows = rows.sort_values("scenario_normalized")

    numeric_columns = [
        "forecast_wasserstein",
        "forecast_gaussian_mmd",
        "forecast_energy_distance_mmd",
        "forecast_hausdorff_loss",
        "embedding_ari",
        "embedding_entropy",
        "lineage_auroc",
        "lineage_auprc",
        "baseline_auroc",
    ]
    for column in numeric_columns:
        rows[column] = pd.to_numeric(rows[column], errors="coerce")

    out = REPORTS / "prescient_reduced_validation_metrics_table.md"
    lines = [
        "# PRESCIENT Reduced-Validation Metrics (CPU Low-Memory, HVG2000, scGPT-v1)",
        "",
        "| Scenario | Forecast WD (lower) | Gaussian MMD (lower) | Energy MMD (lower) | Hausdorff (lower) | Embedding ARI (higher) | Entropy (lower) | Lineage AUROC (higher) | Lineage AUPRC (higher) | Baseline AUROC (higher) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, row in rows.iterrows():
        lines.append(
            f"| {row['scenario_normalized']} | "
            f"{row['forecast_wasserstein']:.4f} | "
            f"{row['forecast_gaussian_mmd']:.4f} | "
            f"{row['forecast_energy_distance_mmd']:.2f} | "
            f"{row['forecast_hausdorff_loss']:.2f} | "
            f"{row['embedding_ari']:.4f} | "
            f"{row['embedding_entropy']:.4f} | "
            f"{row['lineage_auroc']:.4f} | "
            f"{row['lineage_auprc']:.4f} | "
            f"{row['baseline_auroc']:.4f} |"
        )
    lines.extend(
        [
            "",
            "Note: these rows are `result_class=reduced_validation` and `formal_benchmark=false`; they are excluded from official summaries.",
            "CPU low-memory settings use capped training/reference cells and are intended to validate PRESCIENT integration before a formal run.",
        ]
    )
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def write_lineage_figure(rows: pd.DataFrame) -> Path:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    plot_df = rows[
        ["scenario_normalized", "method_label", "lineage_auroc", "baseline_auroc"]
    ].copy()
    plot_df["lineage_auroc"] = pd.to_numeric(plot_df["lineage_auroc"], errors="coerce")
    plot_df["baseline_auroc"] = pd.to_numeric(plot_df["baseline_auroc"], errors="coerce")

    scenarios = ["A", "B", "C"]
    methods = ["WOT", "CellRank2", "scNODE", "PRESCIENT"]
    colors = {
        "WOT": "#377eb8",
        "CellRank2": "#4daf4a",
        "scNODE": "#984ea3",
        "PRESCIENT": "#e41a1c",
    }
    x = np.arange(len(scenarios))
    width = 0.18

    fig, ax = plt.subplots(figsize=(9.2, 5.2), dpi=180)
    for i, method in enumerate(methods):
        values = []
        for scenario in scenarios:
            sub = plot_df[
                (plot_df["scenario_normalized"] == scenario)
                & (plot_df["method_label"] == method)
            ]
            values.append(float(sub["lineage_auroc"].iloc[0]) if len(sub) else np.nan)
        bars = ax.bar(
            x + (i - 1.5) * width,
            values,
            width,
            label=method,
            color=colors[method],
            edgecolor="white",
            linewidth=0.8,
        )
        for bar, value in zip(bars, values):
            if np.isfinite(value):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    value + 0.015,
                    f"{value:.3f}",
                    ha="center",
                    va="bottom",
                    fontsize=8,
                )

    baseline_values = []
    for scenario in scenarios:
        sub = plot_df[plot_df["scenario_normalized"] == scenario]["baseline_auroc"].dropna()
        baseline_values.append(float(sub.iloc[0]) if len(sub) else np.nan)
    ax.plot(
        x,
        baseline_values,
        color="#555555",
        marker="o",
        linewidth=1.8,
        label="Correlation baseline",
    )
    for xi, value in zip(x, baseline_values):
        if np.isfinite(value):
            ax.text(
                xi,
                value - 0.045,
                f"baseline {value:.3f}",
                ha="center",
                va="top",
                fontsize=8,
                color="#444444",
            )

    ax.set_title("Lineage Fidelity AUROC by Model and Scenario", fontsize=13, pad=12)
    ax.set_ylabel("AUROC (higher is better)")
    ax.set_xlabel("Scenario")
    ax.set_xticks(x)
    ax.set_xticklabels(
        ["A\nObserved-time", "B\nExtrapolation", "C\nInterpolation + extrapolation"]
    )
    ax.set_ylim(0, 1.02)
    ax.grid(axis="y", linestyle="--", alpha=0.25)
    ax.legend(frameon=False, ncol=2, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    out = FIG_DIR / "lineage_auroc_by_model_scenario.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    rows = load_presentation_rows()
    outputs = [
        write_formal_summary(rows),
        write_lineage_figure(rows),
        write_scnode_table(rows),
        write_projection_formal_table(rows),
        write_prescient_reduced_table(),
    ]
    for output in outputs:
        print(f"wrote {output}")


if __name__ == "__main__":
    main()
