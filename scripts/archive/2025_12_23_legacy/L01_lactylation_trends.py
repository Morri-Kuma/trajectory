import os
import re
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------
# 0) locate run + inputs
# -----------------------
BASE = r"E:\scgpt\results\c_traj"
LATEST = os.path.join(BASE, "_latest_run.txt")
with open(LATEST, "r", encoding="utf-8") as f:
    RUNPATH = f.read().strip()

C4 = os.path.join(RUNPATH, "adata", "C4_wot_labels_scoreTarget.h5ad")
if not os.path.exists(C4):
    C4 = os.path.join(RUNPATH, "adata", "C4_wot_labels.h5ad")
if not os.path.exists(C4):
    C4 = os.path.join(RUNPATH, "adata", "C2_pseudotime.h5ad")

FULL = r"E:\scgpt\data\processed\GSE230617_qc.h5ad"
if not os.path.exists(FULL):
    raise FileNotFoundError(f"FULL expression h5ad not found: {FULL}")

OUTDIR = os.path.join(RUNPATH, "lactylation_v2")
FIGDIR = os.path.join(OUTDIR, "figures")
TABDIR = os.path.join(OUTDIR, "tables")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(TABDIR, exist_ok=True)

print("[INFO] C4/labels =", C4)
print("[INFO] FULL expr  =", FULL)
print("[INFO] OUTDIR     =", OUTDIR)

# -----------------------
# 1) load
# -----------------------
ad_lab = sc.read_h5ad(C4)
ad_full = sc.read_h5ad(FULL)

ad_lab.obs_names = ad_lab.obs_names.astype(str)
ad_full.obs_names = ad_full.obs_names.astype(str)

print("[INFO] ad_lab  n_cells =", ad_lab.n_obs, "n_genes =", ad_lab.n_vars)
print("[INFO] ad_full n_cells =", ad_full.n_obs, "n_genes =", ad_full.n_vars)
print("[INFO] ad_lab  obs_names head:", ad_lab.obs_names[:5].tolist())
print("[INFO] ad_full obs_names head:", ad_full.obs_names[:5].tolist())

# -----------------------
# 2) alignment
#   priority:
#   (A) exact obs_names overlap
#   (B) join by gsm+barcode etc (if exists)
#   (C) FALLBACK: positional alignment if lab obs_names are 0..n-1 and n_obs equal
# -----------------------
# (A) direct overlap
inter = pd.Index(ad_full.obs_names).intersection(pd.Index(ad_lab.obs_names))
if len(inter) > 1000:
    print(f"[INFO] direct obs_names overlap = {len(inter)}; aligning by obs_names.")
    ad_full = ad_full[inter].copy()
    ad_lab = ad_lab[inter].copy()
else:
    # (C) positional fallback
    can_positional = False
    if ad_full.n_obs == ad_lab.n_obs:
        try:
            idx = pd.to_numeric(pd.Index(ad_lab.obs_names), errors="raise").astype(int).to_numpy()
            # check it's a permutation of 0..n-1 (or at least within range and unique)
            if idx.min() == 0 and idx.max() == ad_full.n_obs - 1 and len(np.unique(idx)) == ad_full.n_obs:
                can_positional = True
        except Exception:
            can_positional = False

    if can_positional:
        print("[WARN] No reliable cell-id key found; using POSITIONAL alignment (lab obs_names are 0..n-1).")
        # reorder FULL to lab order by integer indices
        ad_full = ad_full[idx, :].copy()
        # keep lab as-is
        # add a stable key for future use
        ad_full.obs["cell_key_full"] = ad_full.obs_names.astype(str)
        ad_lab.obs["cell_key_full"] = ad_full.obs_names.astype(str)
    else:
        raise RuntimeError(
            "No overlapping cell IDs and positional alignment not possible.\n"
            "This means C4 lost barcode/cell_id and is not 0..n-1 in order.\n"
            "Fix: preserve barcode (or a gsm+barcode cell_key) in C1/C2/C4 outputs."
        )

print("[INFO] aligned ad_full n_cells =", ad_full.n_obs)
print("[INFO] aligned ad_lab  n_cells =", ad_lab.n_obs)

# -----------------------
# 3) copy embeddings + labels (LAB -> FULL)
# -----------------------
if "X_umap" not in ad_lab.obsm:
    raise RuntimeError("X_umap missing in label h5ad. Use C1/C2/C4 outputs that contain UMAP.")
ad_full.obsm["X_umap"] = ad_lab.obsm["X_umap"].copy()

# copy useful obs columns
for col in [
    "stage", "day", "day_numeric", "gsm", "sample_name",
    "pseudotime", "p_iPSC", "right_traj",
    "pluri_score", "failed", "is_failed", "fail", "leiden"
]:
    if col in ad_lab.obs.columns:
        ad_full.obs[col] = ad_lab.obs[col].copy()

if "day_numeric" not in ad_full.obs.columns and "day" in ad_full.obs.columns:
    ad_full.obs["day_numeric"] = pd.to_numeric(ad_full.obs["day"], errors="coerce")

# -----------------------
# 4) make gsm_short (robust ordering)
# -----------------------
if "gsm" in ad_full.obs.columns:
    s = ad_full.obs["gsm"].astype(str)
    gsm_short = s.str.extract(r"(GSM\d+)", expand=False)
    # if most extracted, use it; else try from full obs_names prefix
    if gsm_short.notna().mean() < 0.8:
        gsm_short = pd.Series(ad_full.obs_names).str.extract(r"(GSM\d+)", expand=False).values
    ad_full.obs["gsm_short"] = gsm_short.astype(str)
else:
    ad_full.obs["gsm_short"] = pd.Series(ad_full.obs_names).str.extract(r"(GSM\d+)", expand=False).astype(str).values

GSM_ORDER_SHORT = [
    "GSM7229998","GSM7229999","GSM7230000","GSM7230001","GSM7230002",
    "GSM7230003","GSM7230004","GSM7230005","GSM7230006","GSM7230007",
    "GSM7230008","GSM7230009","GSM7230010","GSM7230011","GSM7230012"
]
if ad_full.obs["gsm_short"].isin(GSM_ORDER_SHORT).mean() >= 0.9:
    gsm_order = GSM_ORDER_SHORT
else:
    gsm_order = sorted(ad_full.obs["gsm_short"].unique().tolist())
    print("[WARN] gsm_short not matching expected list; using sorted unique gsm_short.")

ad_full.obs["gsm_ordered"] = pd.Categorical(ad_full.obs["gsm_short"], categories=gsm_order, ordered=True)

# nice label for x-axis
if "stage" in ad_full.obs.columns and "day_numeric" in ad_full.obs.columns:
    ad_full.obs["timepoint_label"] = (
        ad_full.obs["gsm_short"].astype(str) + "_" +
        ad_full.obs["stage"].astype(str) + "_Day" +
        ad_full.obs["day_numeric"].astype(str)
    )
else:
    ad_full.obs["timepoint_label"] = ad_full.obs["gsm_short"].astype(str)

# -----------------------
# 5) gene sets
# -----------------------
GENESETS = {
    "lactate_metabolism": [
        "LDHA","LDHB","SLC16A1","SLC16A3","SLC16A7","SLC5A12","SLC5A8",
        "SLC2A1","HK2","PFKP","ALDOA","GAPDH","ENO1","PKM","PDK1","PDK3","HIF1A"
    ],
    "lactylation_machinery": [
        "EP300","CREBBP","HDAC1","HDAC2","HDAC3","SIRT2"
    ],
}

var_u = set(pd.Index(ad_full.var_names.astype(str)).str.upper())

def keep_genes(gs):
    keep, miss = [], []
    for g in gs:
        if g.upper() in var_u:
            if g in ad_full.var_names:
                keep.append(g)
            else:
                idx = np.where(np.char.upper(ad_full.var_names.astype(str)) == g.upper())[0]
                keep.append(ad_full.var_names[idx[0]])
        else:
            miss.append(g)
    return keep, miss

for k, gs in list(GENESETS.items()):
    keep, miss = keep_genes(gs)
    GENESETS[k] = keep
    print(f"[INFO] {k}: keep={len(keep)} miss={len(miss)}")
    if miss:
        print("  missing:", ",".join(miss))

# -----------------------
# 6) score
# -----------------------
for k, genes in GENESETS.items():
    score_name = f"{k}_score"
    if len(genes) >= 2:
        sc.tl.score_genes(ad_full, gene_list=genes, score_name=score_name, use_raw=False)
    elif len(genes) == 1:
        g = genes[0]
        X = ad_full[:, [g]].X
        ad_full.obs[score_name] = X.toarray().ravel() if hasattr(X, "toarray") else np.asarray(X).ravel()
    else:
        ad_full.obs[score_name] = np.nan

# -----------------------
# 7) tables
# -----------------------
def summarize_by(col):
    df = ad_full.obs[["gsm_ordered", "timepoint_label", col]].copy().dropna()
    g = (df.groupby(["gsm_ordered", "timepoint_label"], observed=True)[col]
           .agg(["count","median","mean","std"])
           .reset_index()
           .sort_values(["gsm_ordered"]))
    return g

for col in ["lactate_metabolism_score", "lactylation_machinery_score"]:
    g = summarize_by(col)
    g.to_csv(os.path.join(TABDIR, f"{col}_by_timepoint.tsv"), sep="\t", index=False)
    print("[OK] wrote", f"{col}_by_timepoint.tsv")

if "right_traj" in ad_full.obs.columns:
    df = ad_full.obs[["gsm_ordered","timepoint_label","right_traj",
                      "lactate_metabolism_score","lactylation_machinery_score"]].copy()
    df = df.dropna(subset=["right_traj"])
    df["right_traj"] = df["right_traj"].astype(int)
    out = (
        df.groupby(["gsm_ordered","timepoint_label","right_traj"], observed=True)
          .agg(n=("right_traj","size"),
               lactate_med=("lactate_metabolism_score","median"),
               mach_med=("lactylation_machinery_score","median"))
          .reset_index()
          .sort_values(["gsm_ordered","right_traj"])
    )
    out.to_csv(os.path.join(TABDIR, "scores_by_timepoint_right_traj.tsv"), sep="\t", index=False)
    print("[OK] wrote scores_by_timepoint_right_traj.tsv")

# -----------------------
# 8) plots
# -----------------------
sc.set_figure_params(dpi=150)

for col in ["lactate_metabolism_score", "lactylation_machinery_score"]:
    sc.pl.umap(ad_full, color=col, show=False, size=5)
    outp = os.path.join(FIGDIR, f"umap_{col}.png")
    plt.tight_layout()
    plt.savefig(outp, dpi=200)
    plt.close()
    print("[OK] wrote", outp)

if "p_iPSC" in ad_full.obs.columns:
    sc.pl.umap(ad_full, color="p_iPSC", show=False, size=5)
    outp = os.path.join(FIGDIR, "umap_p_iPSC.png")
    plt.tight_layout()
    plt.savefig(outp, dpi=200)
    plt.close()
    print("[OK] wrote", outp)

def box_by_time(col):
    df = ad_full.obs[["gsm_ordered", "timepoint_label", col]].dropna().copy()
    order = list(df["gsm_ordered"].cat.categories)
    labels, data = [], []
    for o in order:
        sub = df[df["gsm_ordered"] == o]
        if sub.shape[0] == 0:
            continue
        labels.append(sub["timepoint_label"].iloc[0])
        data.append(sub[col].astype(float).values)
    plt.figure(figsize=(12,4))
    plt.boxplot(data, tick_labels=labels, showfliers=False)
    plt.xticks(rotation=45, ha="right")
    plt.ylabel(col)
    plt.title(f"{col} by timepoint")
    plt.tight_layout()
    outp = os.path.join(FIGDIR, f"box_{col}_by_timepoint.png")
    plt.savefig(outp, dpi=200)
    plt.close()
    print("[OK] wrote", outp)

for col in ["lactate_metabolism_score","lactylation_machinery_score"]:
    box_by_time(col)

if "p_iPSC" in ad_full.obs.columns:
    x = pd.to_numeric(ad_full.obs["p_iPSC"], errors="coerce").to_numpy()
    for col in ["lactate_metabolism_score","lactylation_machinery_score"]:
        y = pd.to_numeric(ad_full.obs[col], errors="coerce").to_numpy()
        m = np.isfinite(x) & np.isfinite(y)
        plt.figure(figsize=(5,4))
        plt.scatter(x[m], y[m], s=2, alpha=0.3)
        plt.xlabel("p_iPSC")
        plt.ylabel(col)
        plt.title(f"{col} vs p_iPSC")
        plt.tight_layout()
        outp = os.path.join(FIGDIR, f"scatter_{col}_vs_p_iPSC.png")
        plt.savefig(outp, dpi=200)
        plt.close()
        print("[OK] wrote", outp)

# -----------------------
# 9) save
# -----------------------
OUT_H5AD = os.path.join(OUTDIR, "C4_plus_lactylation_scores_fullgenes.h5ad")
ad_full.write_h5ad(OUT_H5AD)
print("[OK] wrote", OUT_H5AD)
