# -*- coding: utf-8 -*-
"""
Fig2A-like plot (Yan 2021 style):
Left : dynamic genes heatmap (bins are forced by stage order -> within-stage pseudotime)
Middle: cluster trend curves
Right : GO(BP) enrichment per cluster (Enrichr via gseapy)

Default input:
  E:\scgpt\data\gse230659(human\rna_seq

Default output:
  E:\scgpt\results\fig2a_like_gse230659

Run (CMD):
  C:/Users/37620/anaconda3/envs/scgpt/python.exe e:/scgpt/scripts/2025_12_23/make_fig2a_stage_binned.py
"""

import os
import re
import glob
import sys
import warnings

import numpy as np
import pandas as pd

import scanpy as sc
import anndata as ad

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans


# ----------------------------
# Self-check: python + gseapy
# ----------------------------
print("[INFO] sys.executable =", sys.executable)

HAS_GSEAPY = False
gp = None
try:
    import gseapy as gp  # type: ignore
    HAS_GSEAPY = True
except Exception as e:
    HAS_GSEAPY = False
    gp = None
    print("[WARN] gseapy import failed:", repr(e))

print("[INFO] HAS_GSEAPY     =", HAS_GSEAPY)


STAGE_ORDER = ["StageI", "StageII", "StageIII", "hCiPSCs"]


def parse_meta(folder_name: str):
    """
    Parse sample folder names like:
      GSM7229998_StageI_Day0.5-0618
      GSM7230012_hCiPSCs-0618
    """
    m = re.match(r"^(GSM\d+)_(StageI|StageII|StageIII)_Day([0-9.]+)-", folder_name)
    if m:
        gsm, stage, day = m.group(1), m.group(2), float(m.group(3))
        return gsm, stage, day

    m = re.match(r"^(GSM\d+)_hCiPSCs-", folder_name)
    if m:
        gsm = m.group(1)
        return gsm, "hCiPSCs", np.nan

    raise ValueError(f"Unrecognized folder name: {folder_name}")


def find_10x_dir(sample_dir: str):
    """
    Try common 10x output layouts under each sample folder.
    """
    candidates = [
        os.path.join(sample_dir, "filtered_feature_bc_matrix"),
        os.path.join(sample_dir, "outs", "filtered_feature_bc_matrix"),
        os.path.join(sample_dir, "filtered_gene_bc_matrices", "GRCh38"),
        os.path.join(sample_dir, "outs", "filtered_gene_bc_matrices", "GRCh38"),
        sample_dir,
    ]
    for c in candidates:
        if not os.path.isdir(c):
            continue
        mtx = glob.glob(os.path.join(c, "matrix.mtx*"))
        bar = glob.glob(os.path.join(c, "barcodes.tsv*"))
        feat = glob.glob(os.path.join(c, "features.tsv*")) + glob.glob(os.path.join(c, "genes.tsv*"))
        if mtx and bar and feat:
            return c
    raise FileNotFoundError(f"Cannot find 10x matrix under: {sample_dir}")


def load_all_samples(root_dir: str) -> ad.AnnData:
    sample_dirs = [d for d in glob.glob(os.path.join(root_dir, "GSM*")) if os.path.isdir(d)]
    sample_dirs = sorted(sample_dirs)
    if len(sample_dirs) == 0:
        raise FileNotFoundError(f"No GSM* folders found under: {root_dir}")

    adatas = []
    keys = []
    for sd in sample_dirs:
        name = os.path.basename(sd)
        gsm, stage, day = parse_meta(name)
        tenx_dir = find_10x_dir(sd)

        a = sc.read_10x_mtx(tenx_dir, var_names="gene_symbols", cache=False)
        a.var_names_make_unique()

        a.obs["gsm"] = gsm
        a.obs["sample_name"] = name
        a.obs["stage"] = stage
        a.obs["day"] = day

        adatas.append(a)
        keys.append(name)
        print(f"[OK] loaded {name}: {a.n_obs} cells, {a.n_vars} genes from {tenx_dir}")

    adata = ad.concat(adatas, label="sample_name", keys=keys, join="outer", fill_value=0)
    return adata


def preprocess_and_pseudotime(adata: ad.AnnData) -> ad.AnnData:
    """
    Basic preprocessing + DPT pseudotime.
    Root is earliest StageI day.
    """
    # Filter
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=50)

    # Normalize + log
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # keep log-normalized for later bin means

    # Stage ordered categorical
    adata.obs["stage"] = pd.Categorical(adata.obs["stage"], categories=STAGE_ORDER, ordered=True)

    # HVG pipeline for pseudotime (DPT on diffmap)
    sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3")
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50)
    sc.pp.neighbors(adata_hvg, n_neighbors=30, n_pcs=30)
    sc.tl.diffmap(adata_hvg)

    # Root selection: earliest StageI day
    stageI = adata_hvg.obs["stage"] == "StageI"
    if stageI.sum() == 0:
        raise ValueError("No StageI cells found; cannot set DPT root.")

    days = adata_hvg.obs.loc[stageI, "day"].astype(float).values
    min_day = np.nanmin(days)
    root_mask = stageI & (adata_hvg.obs["day"].astype(float) == min_day)
    root_candidates = np.where(root_mask.values)[0]
    if len(root_candidates) == 0:
        root_candidates = np.where(stageI.values)[0]

    adata_hvg.uns["iroot"] = int(root_candidates[0])

    sc.tl.dpt(adata_hvg, n_dcs=10)
    adata.obs["pseudotime"] = adata_hvg.obs["dpt_pseudotime"].values

    return adata


def allocate_bins_by_stage(adata: ad.AnnData, n_bins_total: int = 10):
    """
    Make bins in strict order:
      StageI -> StageII -> StageIII -> hCiPSCs

    Each stage gets at least 1 bin.
    Remaining bins distributed across StageI/II/III proportional to cell count.
    Within each stage, cells are sorted by pseudotime then split into bins.
    """
    if n_bins_total < len(STAGE_ORDER):
        raise ValueError("n_bins_total must be >= number of stages (4).")

    bins = {s: 1 for s in STAGE_ORDER}
    remaining = n_bins_total - len(STAGE_ORDER)

    counts = {s: int((adata.obs["stage"] == s).sum()) for s in STAGE_ORDER}
    weight_stages = [s for s in STAGE_ORDER if s != "hCiPSCs"]  # keep hCiPSCs fixed 1 bin
    w = np.array([max(counts[s], 1) for s in weight_stages], dtype=float)
    w = w / w.sum()

    add = np.floor(w * remaining).astype(int)
    for s, a in zip(weight_stages, add):
        bins[s] += int(a)

    # Fix rounding
    while sum(bins.values()) < n_bins_total:
        s = max(weight_stages, key=lambda x: counts[x])
        bins[s] += 1
    while sum(bins.values()) > n_bins_total:
        s = max(weight_stages, key=lambda x: bins[x])
        if bins[s] > 1:
            bins[s] -= 1
        else:
            break

    # Assign bin labels in strict stage order
    bin_labels = []
    assignments = []
    for s in STAGE_ORDER:
        n = bins[s]
        idx = np.where((adata.obs["stage"] == s).values)[0]
        idx = idx[np.argsort(adata.obs["pseudotime"].values[idx])]
        chunks = np.array_split(idx, n) if len(idx) > 0 else []
        for j, chunk in enumerate(chunks, start=1):
            label = f"{s}_bin{j:02d}"
            bin_labels.append(label)
            for k in chunk:
                assignments.append((k, label))

    bin_series = pd.Series(index=adata.obs_names, dtype="object")
    for k, label in assignments:
        bin_series.iloc[k] = label
    bin_series = bin_series.fillna("hCiPSCs_bin01")

    adata.obs["bin"] = pd.Categorical(bin_series, categories=bin_labels, ordered=True)
    return adata, bin_labels


def pick_dynamic_genes_from_bins(adata: ad.AnnData, top_n: int = 2000):
    """
    Quick proxy "dynamic genes":
      - compute mean expression per bin (on adata.raw, log-normalized)
      - score each gene by variance across bins
      - take top_n genes

    This is for visualization. (If you want tradeSeq-like stats, we can swap later.)
    """
    X = adata.raw.X
    genes = adata.raw.var_names.to_list()
    bin_cats = list(adata.obs["bin"].cat.categories)

    means = []
    for b in bin_cats:
        m = (adata.obs["bin"] == b).values
        if m.sum() == 0:
            means.append(np.zeros((len(genes),), dtype=float))
        else:
            mu = X[m].mean(axis=0)
            means.append(np.asarray(mu).ravel())

    M = np.vstack(means).T  # genes x bins
    v = M.var(axis=1)
    top_idx = np.argsort(v)[::-1][:top_n]

    dyn_genes = [genes[i] for i in top_idx]
    M_top = M[top_idx, :]
    return dyn_genes, M_top, bin_cats


def cluster_gene_patterns(M_z: np.ndarray, k: int = 3):
    """
    Cluster genes by their binned z-score patterns using KMeans,
    then reorder clusters by peak position (early -> late).
    """
    km = KMeans(n_clusters=k, random_state=0, n_init=10)
    lab = km.fit_predict(M_z)

    peak = np.array([M_z[lab == c].mean(axis=0).argmax() if np.any(lab == c) else 0 for c in range(k)])
    order = np.argsort(peak)

    new_lab = np.zeros_like(lab)
    for new_c, old_c in enumerate(order):
        new_lab[lab == old_c] = new_c
    return new_lab


def go_enrich_enrichr(gene_list, organism="Human", db="GO_Biological_Process_2023", top=6) -> pd.DataFrame:
    """
    GO enrichment via Enrichr (gseapy.enrichr). Needs internet.
    Returns: Term, Adjusted P-value, Combined Score (if available)
    """
    if not HAS_GSEAPY:
        return pd.DataFrame(columns=["Term", "Adjusted P-value", "Combined Score"])

    gene_list = [g for g in gene_list if isinstance(g, str) and g.strip() != ""]
    if len(gene_list) < 20:
        return pd.DataFrame(columns=["Term", "Adjusted P-value", "Combined Score"])

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=db,
            organism=organism,
            outdir=None,
            no_plot=True
        )
        res = enr.results
        if res is None or res.shape[0] == 0:
            return pd.DataFrame(columns=["Term", "Adjusted P-value", "Combined Score"])

        # Standard columns from gseapy/enrichr
        need_cols = []
        if "Term" in res.columns:
            need_cols.append("Term")
        else:
            return pd.DataFrame(columns=["Term", "Adjusted P-value", "Combined Score"])

        if "Adjusted P-value" in res.columns:
            padj_col = "Adjusted P-value"
        elif "P-value" in res.columns:
            padj_col = "P-value"
        else:
            return pd.DataFrame(columns=["Term", "Adjusted P-value", "Combined Score"])

        score_col = "Combined Score" if "Combined Score" in res.columns else None

        out = res.copy()
        out = out.sort_values(padj_col, ascending=True).head(top)

        out2 = pd.DataFrame({
            "Term": out["Term"].astype(str).values,
            "Adjusted P-value": out[padj_col].astype(float).values,
            "Combined Score": out[score_col].astype(float).values if score_col else np.nan
        })
        return out2

    except Exception as e:
        warnings.warn(f"[WARN] GO enrichment failed (network/Enrichr issue): {e}")
        return pd.DataFrame(columns=["Term", "Adjusted P-value", "Combined Score"])


def main(
    root_dir: str,
    out_dir: str,
    n_bins: int = 10,
    top_dyn: int = 2000,
    k: int = 3,
    go_db: str = "GO_Biological_Process_2023",
    go_top: int = 6
):
    os.makedirs(out_dir, exist_ok=True)

    # 1) Load
    adata = load_all_samples(root_dir)

    # 2) Pseudotime
    adata = preprocess_and_pseudotime(adata)

    # 3) Stage-binned bins
    adata, _ = allocate_bins_by_stage(adata, n_bins_total=n_bins)

    # 4) Dynamic genes (quick proxy)
    dyn_genes, M, bin_cats = pick_dynamic_genes_from_bins(adata, top_n=top_dyn)

    # Z-score per gene across bins
    M_z = (M - M.mean(axis=1, keepdims=True)) / (M.std(axis=1, keepdims=True) + 1e-8)

    # 5) Cluster patterns
    lab = cluster_gene_patterns(M_z, k=k)

    # sort genes by cluster then by peak bin
    peak = M_z.argmax(axis=1)
    order = np.lexsort((peak, lab))
    dyn_genes = [dyn_genes[i] for i in order]
    M_z = M_z[order, :]
    lab = lab[order]

    # Save gene clusters
    df_genes = pd.DataFrame({"gene": dyn_genes, "cluster": (lab + 1)})
    df_genes.to_csv(os.path.join(out_dir, "dynamic_genes_clusters.tsv"), sep="\t", index=False)

    # Save adata
    adata.write(os.path.join(out_dir, "adata_with_pseudotime_and_bins.h5ad"))

    # 6) GO enrichment per cluster
    go_tables = []
    go_by_cluster = {}
    for c in range(k):
        genes_c = [dyn_genes[i] for i in range(len(dyn_genes)) if lab[i] == c]
        res = go_enrich_enrichr(
            genes_c, organism="Human",
            db=go_db, top=go_top
        )
        res = res.copy()
        res["cluster"] = c + 1
        go_tables.append(res)
        go_by_cluster[c] = res

        print(f"[INFO] cluster {c+1}: genes={len(genes_c)}, go_terms={0 if res is None else res.shape[0]}")

    go_all = pd.concat(go_tables, ignore_index=True) if len(go_tables) > 0 else pd.DataFrame()
    go_all.to_csv(os.path.join(out_dir, "go_enrichment_per_cluster.tsv"), sep="\t", index=False)

    # 7) Plot: heatmap + trends + GO
    fig = plt.figure(figsize=(19, 9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.0, 1.15], wspace=0.35)

    # --- Heatmap (left) ---
    ax0 = fig.add_subplot(gs[0, 0])
    im = ax0.imshow(M_z, aspect="auto", interpolation="nearest")
    ax0.set_title("Dynamic genes (binned by stage → pseudotime)")
    ax0.set_xticks(range(len(bin_cats)))
    ax0.set_xticklabels(bin_cats, rotation=90, fontsize=8)
    ax0.set_yticks([])

    cbar = fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.02)
    cbar.set_label("Z-score (per gene)")

    # cluster separator lines
    for c in range(k - 1):
        idx = np.where(lab == c)[0]
        if len(idx) == 0:
            continue
        ax0.axhline(idx.max() + 0.5, linewidth=1)

    # --- Trends (middle) ---
    ax1 = fig.add_subplot(gs[0, 1])
    x = np.arange(len(bin_cats))

    for c in range(k):
        idx = np.where(lab == c)[0]
        if len(idx) == 0:
            continue

        step = max(1, len(idx) // 200)  # subsample thin lines
        for i in idx[::step]:
            ax1.plot(x, M_z[i], linewidth=0.5, alpha=0.12)

        ax1.plot(x, M_z[idx].mean(axis=0), linewidth=2.5,
                 label=f"Cluster {c+1} (n={len(idx)})")

    ax1.set_title("Cluster trends")
    ax1.set_xticks(x)
    ax1.set_xticklabels(bin_cats, rotation=90, fontsize=8)
    ax1.set_ylabel("Z-score")
    ax1.legend(frameon=False)

    # --- GO panel (right) ---
    sub = gs[0, 2].subgridspec(k, 1, hspace=0.55)

    def _shorten(s: str, maxlen: int = 55) -> str:
        s = str(s)
        return s if len(s) <= maxlen else (s[:maxlen - 3] + "...")

    for c in range(k):
        ax = fig.add_subplot(sub[c, 0])
        res = go_by_cluster.get(c, pd.DataFrame())

        ax.set_title(f"Cluster {c+1} GO(BP) top terms", fontsize=10)

        if res is None or res.shape[0] == 0:
            note = "No enriched terms"
            if not HAS_GSEAPY:
                note += "\n(gseapy not available)"
            ax.text(0.5, 0.5, note, ha="center", va="center")
            ax.set_axis_off()
            continue

        terms = [_shorten(t) for t in res["Term"].astype(str).values[::-1]]
        p = res["Adjusted P-value"].astype(float).values[::-1]
        score = -np.log10(np.clip(p, 1e-300, 1.0))

        ax.barh(np.arange(len(terms)), score)
        ax.set_yticks(np.arange(len(terms)))
        ax.set_yticklabels(terms, fontsize=8)
        ax.set_xlabel("-log10(adj P)")
        ax.grid(axis="x", alpha=0.2)

    fig.tight_layout()
    out_png = os.path.join(out_dir, "fig2a_like_stage_binned_with_GO.png")
    fig.savefig(out_png, dpi=300)
    plt.close(fig)

    print("[DONE]")
    print(" - figure:", out_png)
    print(" - genes :", os.path.join(out_dir, "dynamic_genes_clusters.tsv"))
    print(" - GO    :", os.path.join(out_dir, "go_enrichment_per_cluster.tsv"))
    print(" - h5ad  :", os.path.join(out_dir, "adata_with_pseudotime_and_bins.h5ad"))


if __name__ == "__main__":
    ROOT = r"E:\scgpt\data\gse230659(human\rna_seq"
    OUT = r"E:\scgpt\results\fig2a_like_gse230659"
    main(
        root_dir=ROOT,
        out_dir=OUT,
        n_bins=10,
        top_dyn=2000,
        k=3,
        go_db="GO_Biological_Process_2023",
        go_top=6
    )
