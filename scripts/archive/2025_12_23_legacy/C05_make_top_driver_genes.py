# -*- coding: utf-8 -*-
"""
Make top driver genes table by DE (high vs low) using a continuous label in adata.obs.
Default: high=top 10% label, low=bottom 10%, method=wilcoxon.

Example (PowerShell):
  python scripts\2025_12_23\C05_make_top_driver_genes.py `
    --h5ad "E:\scgpt\results\C_traj\20251224_153430\adata\C4_wot_labels_scoreTarget.h5ad" `
    --label p_iPSC `
    --out "E:\scgpt\results\tables\top_driver_genes.tsv" `
    --top-n 300
"""

import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="Input .h5ad path")
    ap.add_argument("--label", required=True, help="Column name in adata.obs, e.g. p_iPSC / pluri_score / pseudotime")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--top-n", type=int, default=300, help="Top N genes to save (default: 300)")
    ap.add_argument("--q-low", type=float, default=0.10, help="Low quantile (default: 0.10)")
    ap.add_argument("--q-high", type=float, default=0.90, help="High quantile (default: 0.90)")

    # Optional hard thresholds (override quantiles if provided)
    ap.add_argument("--low-threshold", type=float, default=None,
                    help="Hard threshold for low group (<=). Overrides q-low.")
    ap.add_argument("--high-threshold", type=float, default=None,
                    help="Hard threshold for high group (>=). Overrides q-high.")

    ap.add_argument("--method", default="wilcoxon", choices=["wilcoxon", "t-test", "logreg"],
                    help="scanpy rank_genes_groups method (default: wilcoxon)")
    return ap.parse_args()


def safe_read_logfc(res, group="high"):
    """
    Safely read log fold-change from scanpy rank_genes_groups result.
    In scanpy 1.11.x, res fields are often numpy structured arrays / recarray
    and should be accessed by indexing (no .get()).
    """
    # default: all NaN
    names = res["names"][group]
    logfc = np.full(len(names), np.nan, dtype=float)

    # Most common: res["logfoldchanges"][group]
    try:
        lf = res["logfoldchanges"][group]
        # ensure array-like and correct length
        if len(lf) == len(logfc):
            return np.asarray(lf, dtype=float)
    except Exception:
        pass

    # Sometimes direct 1D array (rare)
    try:
        lf = res["logfoldchanges"]
        if hasattr(lf, "__len__") and len(lf) == len(logfc):
            return np.asarray(lf, dtype=float)
    except Exception:
        pass

    return logfc


def main():
    args = parse_args()

    if not os.path.exists(args.h5ad):
        raise FileNotFoundError(f"H5AD not found: {args.h5ad}")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    print(f"[INFO] Reading: {args.h5ad}")
    adata = sc.read_h5ad(args.h5ad)

    if args.label not in adata.obs.columns:
        raise ValueError(
            f"Label '{args.label}' not in adata.obs. Available (first 50): {list(adata.obs.columns)[:50]}"
        )

    # make sure label is numeric
    y = pd.to_numeric(adata.obs[args.label], errors="coerce").values
    finite = np.isfinite(y)
    if finite.mean() < 0.8:
        print(f"[WARN] Only {finite.mean():.2%} of label values are finite. Check your label column.")

    # compute thresholds
    if args.low_threshold is not None:
        lo = float(args.low_threshold)
    else:
        lo = float(np.nanquantile(y, args.q_low))

    if args.high_threshold is not None:
        hi = float(args.high_threshold)
    else:
        hi = float(np.nanquantile(y, args.q_high))

    print(f"[INFO] LABEL={args.label} low<= {lo:.9g} ; high>= {hi:.9g}")

    grp = np.full(adata.n_obs, "mid", dtype=object)
    grp[y <= lo] = "low"
    grp[y >= hi] = "high"
    adata.obs["driver_group"] = grp

    n_low = int((grp == "low").sum())
    n_high = int((grp == "high").sum())
    n_mid = int((grp == "mid").sum())
    print(f"[INFO] group sizes: low={n_low}, high={n_high}, mid={n_mid}")

    if n_low < 50 or n_high < 50:
        print("[WARN] One group has <50 cells. Consider relaxing thresholds (e.g., q-low=0.2, q-high=0.8).")

    sub = adata[adata.obs["driver_group"].isin(["high", "low"])].copy()

    print(f"[INFO] Running rank_genes_groups (method={args.method}) ...")
    sc.tl.rank_genes_groups(
        sub,
        groupby="driver_group",
        groups=["high"],
        reference="low",
        method=args.method,
    )

    res = sub.uns["rank_genes_groups"]

    # IMPORTANT: res is a dict-like mapping to numpy structured arrays; no .get()
    names = res["names"]["high"]
    scores = res["scores"]["high"]
    pvals_adj = res["pvals_adj"]["high"]
    logfc = safe_read_logfc(res, group="high")

    df = pd.DataFrame({
        "gene": names,
        "score": scores,
        "pvals_adj": pvals_adj,
        "logFC": logfc,
    })

    df = df.dropna(subset=["gene"]).drop_duplicates("gene").head(args.top_n)
    df.to_csv(args.out, sep="\t", index=False)

    print(f"[OK] wrote: {args.out} (n={df.shape[0]})")
    print(df.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
