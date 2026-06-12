#!/usr/bin/env python
"""fig8: multi-seed robustness of lineage fidelity on the primary datasets.

Reads benchmark/reports/official_silver/multiseed_ci.csv (5 seeds) and overlays:
  - bars = multi-seed mean single-step lineage AUROC, whiskers = 95% CI,
  - open markers = the single-seed (seed 42) value the manuscript reports,
  - dashed line = the scTimeBench Spearman correlation baseline for that dataset.
Shows that the single-seed leader is not robust and that the baseline is the floor.
"""
from __future__ import annotations

import csv
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.environ.get("TRAJ_PROJECT_ROOT", os.getcwd())
REP = os.path.join(ROOT, "benchmark", "reports", "official_silver")
OUT = os.path.join(ROOT, "manuscript", "figures_v2")
GEN = ["scnode", "prescient", "mioflow"]
PRETTY = {"scnode": "scNODE", "prescient": "PRESCIENT", "mioflow": "MIOFlow"}
DATASETS = ["GSE178325", "GSE230659"]
SCEN = ["A", "B", "C"]
SC_COLOR = {"A": "#2c7fb8", "B": "#d95f0e", "C": "#31a354"}


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    # multi-seed mean + CI
    mean = defaultdict(dict)
    ci = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(REP, "multiseed_ci.csv"))):
        if r["metric"] != "lineage_single_step_auroc":
            continue
        mean[(r["dataset"], r["scenario"])][r["method"]] = num(r["mean"])
        ci[(r["dataset"], r["scenario"])][r["method"]] = num(r["ci95_half_width"])
    # single-seed canonical (seed 42)
    canon = defaultdict(dict)
    for r in csv.DictReader(open(os.path.join(REP, "official_silver_lineage_rankings.csv"))):
        if r["dataset_id"] in DATASETS and r["scenario"] in SCEN and r["method"] in GEN:
            canon[(r["dataset_id"], r["scenario"])][r["method"]] = num(r["single_step_auc_roc"])
    # Spearman baseline
    base = {}
    for r in csv.DictReader(open(os.path.join(REP, "spearman_baseline.csv"))):
        base[r["dataset"]] = num(r["single_step_auroc"])

    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.0), sharey=True)
    for ax, ds in zip(axes, DATASETS):
        x = np.arange(len(GEN))
        w = 0.26
        for i, sc in enumerate(SCEN):
            m = [mean[(ds, sc)].get(g, np.nan) for g in GEN]
            e = [ci[(ds, sc)].get(g, 0.0) for g in GEN]
            xs = x + (i - 1) * w
            ax.bar(xs, m, w, yerr=e, capsize=3, color=SC_COLOR[sc],
                   edgecolor="black", linewidth=0.4, label=f"Scenario {sc}",
                   error_kw=dict(lw=1.0))
            cs = [canon[(ds, sc)].get(g, np.nan) for g in GEN]
            ax.scatter(xs, cs, marker="o", facecolors="none", edgecolors="black",
                       s=26, zorder=5, linewidths=1.0)
        if ds in base:
            ax.axhline(base[ds], ls="--", lw=1.0, color="grey")
            ax.text(len(GEN) - 1.4, base[ds] + 0.01, "Spearman baseline",
                    fontsize=7, color="grey")
        ax.set_xticks(x)
        ax.set_xticklabels([PRETTY[g] for g in GEN], fontsize=8, rotation=12)
        ax.set_title(ds, fontsize=10)
        ax.grid(axis="y", ls=":", lw=0.5, alpha=0.6)
        ax.set_ylim(0, 1.0)
    axes[0].set_ylabel("Single-step lineage AUROC", fontsize=9)
    axes[0].legend(fontsize=7, frameon=False, loc="upper left")
    fig.suptitle("Multi-seed robustness (bars = mean ± 95% CI over 5 seeds; "
                 "open circles = single-seed value)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    figpath = os.path.join(OUT, "fig8_multiseed_robustness.png")
    fig.savefig(figpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote", figpath)


if __name__ == "__main__":
    main()
