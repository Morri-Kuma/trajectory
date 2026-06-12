#!/usr/bin/env python
"""Rank-concordance analysis: quantify how (un)stable the method leaderboard is
across biological systems --- the statistical core of the 'rankings do not
transfer' claim.

For each scTimeBench task (lineage single-step AUROC, forecast Wasserstein,
embedding ARI) we rank the three projection-capable methods (scNODE, PRESCIENT,
MIOFlow) within each official-silver dataset (mean over observed-time scenarios
A/B/C), then report Kendall's W (coefficient of concordance, 0=no agreement,
1=identical rankings) and the full matrix of pairwise Spearman rank correlations
between datasets. Writes a CSV and a bump-chart figure. Pure local; reads the
regenerated official_silver ranking CSVs.
"""
from __future__ import annotations

import csv
import os
import statistics as st
from collections import defaultdict
from itertools import combinations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.environ.get("TRAJ_PROJECT_ROOT", os.getcwd())
REP = os.path.join(ROOT, "benchmark", "reports", "official_silver")
OUT_FIG = os.path.join(ROOT, "manuscript", "figures_v2")
os.makedirs(OUT_FIG, exist_ok=True)

GEN = ["scnode", "prescient", "mioflow"]
PRETTY = {"scnode": "scNODE", "prescient": "PRESCIENT", "mioflow": "MIOFlow"}
DATASETS = ["GSE178325", "GSE230659", "GSE218855", "GSE298212"]
DS_SHORT = {"GSE178325": "hADSC-1", "GSE230659": "hADSC-2",
            "GSE218855": "mouse", "GSE298212": "blood"}
ABC = ["A", "B", "C"]
COLOR = {"scnode": "#1b9e77", "prescient": "#7570b3", "mioflow": "#d95f02"}

# task: (csv, value column, higher_is_better)
TASKS = [
    ("Lineage AUROC", "official_silver_lineage_rankings.csv", "single_step_auc_roc", True),
    ("Forecast WD", "official_silver_forecast_rankings.csv", "wasserstein_distance", False),
    ("Embedding ARI", "official_silver_embedding_rankings.csv", "adjusted_rand_index", True),
]


def load(name):
    with open(os.path.join(REP, name)) as f:
        return list(csv.DictReader(f))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def mean_table(rows, valcol):
    agg = defaultdict(lambda: defaultdict(list))
    for r in rows:
        if r["method"] in GEN and r["dataset_id"] in DATASETS and r["scenario"] in ABC:
            v = num(r.get(valcol))
            if v is not None:
                agg[r["dataset_id"]][r["method"]].append(v)
    return {d: {m: st.mean(agg[d][m]) for m in GEN if agg[d].get(m)} for d in DATASETS if d in agg}


def ranks_from(table, higher_better):
    out = {}
    for d, mv in table.items():
        order = sorted(mv, key=lambda m: (-mv[m] if higher_better else mv[m]))
        out[d] = {m: order.index(m) + 1 for m in order}
    return out


def kendalls_w(ranks, methods, datasets):
    k, n = len(methods), len(datasets)
    Rj = {m: sum(ranks[d][m] for d in datasets) for m in methods}
    Rbar = sum(Rj.values()) / k
    S = sum((Rj[m] - Rbar) ** 2 for m in methods)
    return 12 * S / (n ** 2 * (k ** 3 - k))


def spearman(ra, rb, methods):
    a = np.array([ra[m] for m in methods], float)
    b = np.array([rb[m] for m in methods], float)
    return float(np.corrcoef(a, b)[0, 1])


def main():
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8))
    summary_rows = []
    for ax, (tname, csvname, valcol, hib) in zip(axes, TASKS):
        table = mean_table(load(csvname), valcol)
        dss = [d for d in DATASETS if d in table]
        rk = ranks_from(table, hib)
        W = kendalls_w(rk, GEN, dss)
        # bump chart: x = datasets, y = rank (1 top)
        x = np.arange(len(dss))
        for m in GEN:
            ys = [rk[d][m] for d in dss]
            ax.plot(x, ys, "-o", color=COLOR[m], lw=2, ms=7, label=PRETTY[m])
        ax.set_xticks(x)
        ax.set_xticklabels([DS_SHORT[d] for d in dss], fontsize=8, rotation=15)
        ax.set_yticks([1, 2, 3])
        ax.set_yticklabels(["1st", "2nd", "3rd"], fontsize=8)
        ax.invert_yaxis()
        ax.set_title(f"{tname}\nKendall's $W$ = {W:.2f}", fontsize=9)
        ax.grid(axis="y", ls=":", lw=0.5, alpha=0.6)
        summary_rows.append({"task": tname, "kendalls_w": round(W, 3),
                             "datasets": ";".join(dss)})
        for d1, d2 in combinations(dss, 2):
            summary_rows.append({"task": tname,
                                 "pair": f"{DS_SHORT[d1]}~{DS_SHORT[d2]}",
                                 "spearman_rho": round(spearman(rk[d1], rk[d2], GEN), 3)})
    axes[0].legend(fontsize=7, frameon=False, loc="center left")
    fig.suptitle("Method-ranking (in)stability across biological systems", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    figpath = os.path.join(OUT_FIG, "fig7_rank_concordance.png")
    fig.savefig(figpath, dpi=200, bbox_inches="tight")
    plt.close(fig)

    outcsv = os.path.join(REP, "rank_concordance.csv")
    keys = ["task", "kendalls_w", "datasets", "pair", "spearman_rho"]
    with open(outcsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        for r in summary_rows:
            w.writerow({k: r.get(k, "") for k in keys})
    print("wrote", figpath)
    print("wrote", outcsv)
    for r in summary_rows:
        print("  ", r)


if __name__ == "__main__":
    main()
