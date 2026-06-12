#!/usr/bin/env python
"""fig4: lineage single-step AUROC (mean A/B/C) across all marker-silver systems
plus the OSKM ground-truth system, by method. Shows the cross-system ranking flip."""
import csv, json, glob, os
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.environ.get("TRAJ_PROJECT_ROOT", os.getcwd())
REP = os.path.join(ROOT, "benchmark", "reports", "official_silver")
OUT = os.path.join(ROOT, "manuscript", "figures_v2")
GEN = ["scnode", "prescient", "mioflow"]
PRETTY = {"scnode": "scNODE", "prescient": "PRESCIENT", "mioflow": "MIOFlow"}
SCEN = ["A", "B", "C"]
def num(x):
    try: return float(x)
    except: return None
lin = list(csv.DictReader(open(os.path.join(REP, "official_silver_lineage_rankings.csv"))))
sil = defaultdict(lambda: defaultdict(list))
for r in lin:
    if r["method"] in GEN and r["scenario"] in SCEN:
        v = num(r.get("single_step_auc_roc"))
        if v is not None: sil[r["dataset_id"]][r["method"]].append(v)
osk = defaultdict(list)
for p in glob.glob(os.path.join(ROOT, "benchmark/results/*/gse242424_oskm_ground_truth_*_hvg2000_formal/lineage_metrics.json")):
    meth = os.path.normpath(p).split(os.sep)[-3]
    try: d = json.load(open(p))
    except: continue
    v = d.get("auroc")
    if meth in GEN and v is not None: osk[meth].append(v)
order = ["GSE178325", "GSE230659", "GSE218855", "GSE298212"]
labels = {"GSE178325":"hADSC-1","GSE230659":"hADSC-2","GSE218855":"mouse","GSE298212":"blood","GSE242424":"OSKM*"}
systems = [d for d in order if d in sil] + ["GSE242424"]
means = {d: {m: float(np.mean(sil[d][m])) for m in GEN if sil[d].get(m)} for d in order if d in sil}
means["GSE242424"] = {m: float(np.mean(osk[m])) for m in GEN if osk.get(m)}
fig, ax = plt.subplots(figsize=(8.4, 3.8))
x = np.arange(len(systems)); w = 0.26
mcol = {"scnode":"#1b9e77","prescient":"#7570b3","mioflow":"#d95f02"}
for i, m in enumerate(GEN):
    vals = [means[d].get(m, np.nan) for d in systems]
    ax.bar(x + (i-1)*w, vals, w, label=PRETTY[m], color=mcol[m], edgecolor="black", linewidth=0.4)
ax.axhline(0.5, ls="--", lw=0.9, color="grey")
ax.set_xticks(x); ax.set_xticklabels([labels[d] for d in systems], fontsize=9)
ax.set_ylabel("Single-step lineage AUROC (mean A/B/C)", fontsize=9)
ax.set_ylim(0, 1.0)
ax.legend(fontsize=8, frameon=False, ncol=3, loc="upper center")
ax.grid(axis="y", ls=":", lw=0.5, alpha=0.6)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig4_cross_system.png"), dpi=200, bbox_inches="tight")
print("wrote fig4_cross_system.png ; means:")
for d in systems: print(" ", labels[d], {m: round(means[d].get(m, float('nan')),3) for m in GEN})
