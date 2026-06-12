#!/usr/bin/env python
"""Regenerate manuscript figures fig1/fig2/fig3 from the official_silver ranking CSVs.

These three panels show the two primary chemical-reprogramming datasets
(GSE178325, GSE230659) by method and observed-time scenario (A/B/C):
  fig1_lineage_auroc.png   single-step graph-sim lineage AUROC (incl. WOT/CellRank2 @A)
  fig2_forecast_wd.png     exact OT-loss Wasserstein (lower better; per-panel y-scale)
  fig3_embedding_ari.png   embedding-coherence Leiden ARI vs reference milestones
Reads benchmark/reports/official_silver/*_rankings.csv; writes manuscript/figures_v2/.
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
os.makedirs(OUT, exist_ok=True)

PRIMARIES = ["GSE178325", "GSE230659"]
SCEN = ["A", "B", "C"]
SC_COLOR = {"A": "#2c7fb8", "B": "#d95f0e", "C": "#31a354"}
GEN = ["scnode", "prescient", "mioflow"]
LIN_METHODS = ["scnode", "prescient", "mioflow", "wot", "cellrank2"]
PRETTY = {"scnode": "scNODE", "prescient": "PRESCIENT", "mioflow": "MIOFlow",
          "wot": "WOT", "cellrank2": "CellRank2"}


def load(name):
    with open(os.path.join(REP, name)) as f:
        return list(csv.DictReader(f))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def lookup(rows, valcol):
    d = defaultdict(dict)  # d[(dataset, method)][scenario] = value
    for r in rows:
        d[(r["dataset_id"], r["method"])][r["scenario"]] = num(r.get(valcol))
    return d


def grouped_panel(ax, d, methods, title, ylabel, chance=None, ymax=None):
    x = np.arange(len(methods))
    w = 0.26
    for i, sc in enumerate(SCEN):
        vals = [d.get((title_ds(title), m), {}).get(sc, np.nan) for m in methods]
        ax.bar(x + (i - 1) * w, vals, w, label=f"Scenario {sc}",
               color=SC_COLOR[sc], edgecolor="black", linewidth=0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([PRETTY[m] for m in methods], rotation=20, ha="right", fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=9)
    if chance is not None:
        ax.axhline(chance, ls="--", lw=0.9, color="grey")
    if ymax is not None:
        ax.set_ylim(0, ymax)
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.6)


def title_ds(title):
    return title.split()[0]


def make(valcol, csvname, methods, ylabel, outfile, chance=None, share_y=True):
    d = lookup(load(csvname), valcol)
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))
    ymax = None
    if share_y:
        allv = [v for (ds, m), sv in d.items() if ds in PRIMARIES
                for v in sv.values() if not np.isnan(v)]
        ymax = max(allv) * 1.18 if allv else None
    for ax, ds in zip(axes, PRIMARIES):
        panel_ymax = ymax
        if not share_y:
            vv = [v for (dd, m), sv in d.items() if dd == ds
                  for v in sv.values() if not np.isnan(v)]
            panel_ymax = max(vv) * 1.18 if vv else None
        grouped_panel(ax, d, methods, f"{ds}", ylabel, chance=chance, ymax=panel_ymax)
    axes[0].legend(fontsize=8, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, outfile), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {outfile}")


def make_cross_system():
    """fig4: lineage single-step AUROC (mean over A/B/C) across all marker-silver
    systems plus the OSKM ground-truth system, by method. Shows the ranking flip."""
    import json
    import glob
    lin = load("official_silver_lineage_rankings.csv")
    sil = defaultdict(lambda: defaultdict(list))  # ds -> method -> [auroc]
    for r in lin:
        if r["method"] in GEN and r["scenario"] in SCEN:
            v = num(r.get("single_step_auc_roc"))
            if v is not None:
                sil[r["dataset_id"]][r["method"]].append(v)
    # OSKM (GSE242424) uses the ground-truth protocol column 'auroc'
    osk = defaultdict(list)
    for p in glob.glob(os.path.join(ROOT, "benchmark/results/*/gse242424_oskm_ground_truth_*_hvg2000_formal/lineage_metrics.json")):
        meth = p.split(os.sep)[-2].split(os.sep)[0]
        meth = os.path.normpath(p).split(os.sep)[-3]
        try:
            d = json.load(open(p))
        except Exception:
            continue
        v = d.get("auroc")
        if meth in GEN and v is not None:
            osk[meth].append(v)
    order = ["GSE178325", "GSE230659", "GSE218855", "GSE298212"]
    labels = {"GSE178325": "hADSC-1", "GSE230659": "hADSC-2",
              "GSE218855": "mouse", "GSE298212": "blood", "GSE242424": "OSKM*"}
    systems = [d for d in order if d in sil] + ["GSE242424"]
    means = {}
    for d in order:
        if d in sil:
            means[d] = {m: float(np.mean(sil[d][m])) for m in GEN if sil[d].get(m)}
    means["GSE242424"] = {m: float(np.mean(osk[m])) for m in GEN if osk.get(m)}
    fig, ax = plt.subplots(figsize=(8.4, 3.8))
    x = np.arange(len(systems))
    w = 0.26
    mcol = {"scnode": "#1b9e77", "prescient": "#7570b3", "mioflow": "#d95f02"}
    for i, m in enumerate(GEN):
        vals = [means[d].get(m, np.nan) for d in systems]
        ax.bar(x + (i - 1) * w, vals, w, label=PRETTY[m], color=mcol[m],
               edgecolor="black", linewidth=0.4)
    ax.axhline(0.5, ls="--", lw=0.9, color="grey")
    ax.set_xticks(x)
    ax.set_xticklabels([labels[d] for d in systems], fontsize=9)
    ax.set_ylabel("Single-step lineage AUROC (mean A/B/C)", fontsize=9)
    ax.set_ylim(0, 1.0)
    ax.legend(fontsize=8, frameon=False, ncol=3, loc="upper center")
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.6)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig4_cross_system.png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("wrote fig4_cross_system.png (OSKM* = ground-truth protocol)")


if __name__ == "__main__":
    make("single_step_auc_roc", "official_silver_lineage_rankings.csv", LIN_METHODS,
         "Single-step AUROC", "fig1_lineage_auroc.png", chance=0.5, share_y=True)
    make("wasserstein_distance", "official_silver_forecast_rankings.csv", GEN,
         "Wasserstein (OT loss)", "fig2_forecast_wd.png", share_y=False)
    make("adjusted_rand_index", "official_silver_embedding_rankings.csv", GEN,
         "Leiden ARI", "fig3_embedding_ari.png", chance=0.0, share_y=True)
    make_cross_system()
