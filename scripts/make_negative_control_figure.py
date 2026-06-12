#!/usr/bin/env python
"""Supplementary figure for the lineage negative control (all 72 runs).
Two panels (observed A/B/C, pseudotime D/E/F): real single-step AUROC per
method x system x scenario vs the label-permutation null. One conclusion: real
scores clear the empirical null (~0.62) only where genuine lineage signal exists.
Output: manuscript/figures_v2/figS1_negative_control.png
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

df = pd.read_csv("results/test_outputs/negative_control/lineage_negative_control.csv")
OUT = "manuscript/figures_v2/figS1_negative_control.png"
datasets = sorted(df["dataset"].unique())
dpos = {d: i for i, d in enumerate(datasets)}
mcol = {"scnode": "#1b9e77", "prescient": "#d95f02", "mioflow": "#7570b3"}
moff = {"scnode": -0.22, "prescient": 0.0, "mioflow": 0.22}

fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), sharey=True)
for ax, axis in zip(axes, ["observed", "pseudotime"]):
    sub = df[df["axis"] == axis]
    nullmean = sub["perm_null_mean"].mean()
    ax.axhspan(sub["perm_null_mean"].min(), sub["perm_null_p95"].max(), color="0.90", zorder=0)
    ax.axhline(nullmean, color="0.45", ls="--", lw=1, zorder=1,
               label=f"perm-null mean ≈ {nullmean:.2f}")
    ax.axhline(0.5, color="k", ls=":", lw=1, zorder=1, label="nominal chance 0.5")
    for _, r in sub.iterrows():
        x = dpos[r["dataset"]] + moff.get(r["method"], 0.0)
        face = mcol[r["method"]] if r["real_above_null_p95"] else "none"
        ax.scatter(x, r["real_auroc"], s=46, facecolors=face, edgecolors=mcol[r["method"]],
                   linewidths=1.6, zorder=3)
    ax.set_xticks(range(len(datasets)))
    ax.set_xticklabels([d.replace("GSE", "") for d in datasets], fontsize=8)
    ax.set_title(f"{axis} time (n={len(sub)}; filled = real > null p95)", fontsize=10)
    ax.set_xlabel("system (GSE)"); ax.set_ylim(0.35, 1.04)
    ax.legend(fontsize=7, loc="lower right")
axes[0].set_ylabel("single-step lineage AUROC")
# method legend
from matplotlib.lines import Line2D
handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=c, markeredgecolor=c,
                  markersize=8, label=m) for m, c in mcol.items()]
axes[0].legend(handles=handles, fontsize=7, loc="upper left", title="method", title_fontsize=7)
fig.suptitle("Lineage negative control: real prediction vs label-permutation null "
             "(2000 perms; grey band = null mean→p95)", fontsize=11)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT, dpi=200, bbox_inches="tight"); print("wrote", OUT)
