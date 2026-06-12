#!/usr/bin/env python3
"""Build the GSE218855 HVG2000 benchmark input (mouse MEF fast chemical reprogramming).

GSE218855 (Chen et al. 2023, Nat Cell Biol, doi:10.1038/s41556-023-01193-x): fast
chemical reprogramming of mouse embryonic fibroblasts (MEF) to iPSCs through a
diapause-like state. GEO ships this as per-sample 10x MatrixMarket triplets
(GSM*_<sample>-{barcodes,features,matrix}), samples D0/D4/D8/D12/iPSC — NOT a single
.h5ad. This builder reads those flat triplets directly (from data/gse218855/ and its
_raw/ subdir), so no manual reorganization is needed after scripts/download_new_datasets.sh.

Cross-species robustness dataset; UNLABELED -> milestone labels are added afterwards by
the marker-FM annotation step using milestone_markers_mouse.yaml (mouse markers; MEF
start at abs_day==0). Marker matching is case-insensitive, so uppercase markers match
mouse Title-case var_names.

    python scripts/build_gse218855_fcr_input.py [--overwrite]
"""
from __future__ import annotations

import argparse
import glob
import gzip
import json
import os
import re
import time as _time
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp
import yaml

DATASET_ID = "GSE218855"
DATA_FOLDER = "gse218855"
N_TOP_GENES = 2000
NORM_TARGET_SUM = 1e4
T0 = _time.time()

DAY_TO_ABS_DAY = {"d0": 0.0, "d4": 4.0, "d8": 8.0, "d12": 12.0, "ipsc": 16.0, "ips": 16.0}
QC_PARAMS = {"min_genes_per_cell": 200, "max_genes_per_cell": 10_000,
             "min_counts_per_cell": 500, "max_pct_mito": 25.0, "min_cells_per_gene": 3}


def find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env and Path(env).exists():
        return Path(env)
    here = Path(__file__).resolve().parent
    for c in [here, *here.parents]:
        if (c / "benchmark").is_dir() and (c / "scripts").is_dir():
            return c
    return here.parent


def _force_include_marker_genes(adata, root, markers_filename, dataset_id):
    """Return actual var_names that must be force-kept so the marker-FM milestone
    annotation always finds >=2 markers per milestone. Canonical pluripotency /
    lineage TFs are low-variance and get dropped by top-HVG selection even when
    expressed; with min_cells_per_gene=3 they survive QC, so we union them back in.
    Matching is case-insensitive (mirrors marker_utils.build_gene_index_map)."""
    path = Path(root) / "benchmark" / "annotation" / markers_filename
    cfg = (yaml.safe_load(path.read_text(encoding="utf-8-sig")) or {}).get(dataset_id, {})
    want = set()
    for ms in (cfg.get("marker_sets") or {}).values():
        for g in (ms.get("genes") or []):
            want.add(str(g).upper())
    upper2name = {}
    for name in map(str, adata.var_names):
        upper2name.setdefault(name.upper(), name)
    return [upper2name[u] for u in want if u in upper2name]


def _abs_day(token: str) -> float:
    s = re.sub(r"[^a-z0-9]", "", token.lower())
    for k, d in DAY_TO_ABS_DAY.items():
        if k == s or k in s:
            return d
    m = re.search(r"d(\d+)", s)
    return float(m.group(1)) if m else float("nan")


def _read_lines(path: str):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as f:
        return [ln.rstrip("\n").split("\t") for ln in f]


def _load_triplet(mtx, bc, feat, sample):
    X = sio.mmread(mtx).astype(np.float32)          # genes x cells
    X = sp.csr_matrix(X).T                           # cells x genes
    barcodes = [r[0] for r in _read_lines(bc)]
    feats = _read_lines(feat)
    symbols = [(r[1] if len(r) >= 2 else r[0]) for r in feats]   # symbol col if present
    a = ad.AnnData(X=X)
    if X.shape[0] != len(barcodes) or X.shape[1] != len(symbols):
        raise ValueError(f"{sample}: shape {X.shape} != ({len(barcodes)},{len(symbols)})")
    a.obs_names = [f"{sample}_{b}" for b in barcodes]
    a.var_names = pd.Index(symbols)
    a.var_names_make_unique()
    a.obs["sample_id"] = sample
    a.obs["abs_day"] = _abs_day(sample)
    return a


def _discover_samples(data_dir: Path):
    """Return {sample_token: (matrix, barcodes, features)} from flat GEO 10x files."""
    search = [str(data_dir), str(data_dir / "_raw")]
    mtx_files = []
    for d in search:
        mtx_files += glob.glob(os.path.join(d, "*matrix.mtx*"))
    out = {}
    for mtx in sorted(set(mtx_files)):
        base = os.path.basename(mtx)
        stem = re.sub(r"[-_.]?matrix\.mtx(\.gz)?$", "", base)        # strip role+ext
        token = re.sub(r"^gsm\d+[-_]?", "", stem, flags=re.I) or stem  # drop GSM id
        d = os.path.dirname(mtx)
        def sib(role_alts):
            for r in role_alts:
                for ext in (".tsv.gz", ".tsv"):
                    p = os.path.join(d, f"{stem}-{r}{ext}")
                    if os.path.exists(p):
                        return p
                    p = os.path.join(d, f"{stem}_{r}{ext}")
                    if os.path.exists(p):
                        return p
            return None
        bc = sib(["barcodes"]); feat = sib(["features", "genes"])
        if bc and feat:
            out[token] = (mtx, bc, feat)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--project-root", type=Path, default=None)
    args = ap.parse_args()
    import scanpy as sc

    root = args.project_root or find_project_root()
    data_dir = root / "data" / DATA_FOLDER
    out_dir = root / "benchmark" / "inputs" / "gse218855_marker_fm_transition_silver_hvg2000"
    out_h5ad = out_dir / "GSE218855_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"
    print("=" * 64)
    print(f"Build {DATASET_ID} HVG2000 benchmark input (mouse, per-sample 10x)")
    print(f"  data: {data_dir}\n  out : {out_h5ad}")
    print("=" * 64)
    if out_h5ad.exists() and not args.overwrite:
        raise FileExistsError(f"{out_h5ad} exists; pass --overwrite.")
    samples = _discover_samples(data_dir)
    if not samples:
        raise FileNotFoundError(
            f"No 10x triplets found under {data_dir} or {data_dir}/_raw. Run "
            f"scripts/download_new_datasets.sh first.")
    print(f"  discovered samples: {sorted(samples)}")
    out_dir.mkdir(parents=True, exist_ok=True)

    parts = []
    for token, (mtx, bc, feat) in sorted(samples.items()):
        print(f"  loading {token}: {os.path.basename(mtx)}")
        parts.append(_load_triplet(mtx, bc, feat, token))
    adata = ad.concat(parts, join="outer", index_unique=None)
    adata.obs_names_make_unique()
    if adata.obs["abs_day"].isna().any():
        bad = sorted(adata.obs.loc[adata.obs["abs_day"].isna(), "sample_id"].unique())
        raise ValueError(f"samples without abs_day: {bad}; edit DAY_TO_ABS_DAY.")
    print(f"  concatenated -> {adata.shape}")

    # QC + normalize + HVG (raw counts -> normalize+log)
    adata.var_names_make_unique()
    n0 = adata.n_obs
    mito = adata.var_names.str.upper().str.startswith("MT-")
    if mito.sum():
        adata.obs["pct_counts_mt"] = (np.asarray(adata[:, mito].X.sum(1)).ravel()
                                      / np.asarray(adata.X.sum(1)).ravel().clip(1) * 100)
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(1)).ravel()
    adata.obs["n_counts"] = np.asarray(adata.X.sum(1)).ravel()
    keep = ((adata.obs["n_genes"] >= QC_PARAMS["min_genes_per_cell"])
            & (adata.obs["n_genes"] <= QC_PARAMS["max_genes_per_cell"])
            & (adata.obs["n_counts"] >= QC_PARAMS["min_counts_per_cell"]))
    if "pct_counts_mt" in adata.obs:
        keep &= adata.obs["pct_counts_mt"] <= QC_PARAMS["max_pct_mito"]
    adata = adata[keep].copy()
    gci = np.asarray((adata.X > 0).sum(0)).ravel()
    adata = adata[:, gci >= QC_PARAMS["min_cells_per_gene"]].copy()
    print(f"  QC: {n0} -> {adata.n_obs} cells, {adata.n_vars} genes")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=NORM_TARGET_SUM)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, flavor="seurat")
    forced = _force_include_marker_genes(adata, root, "milestone_markers_mouse.yaml", DATASET_ID)
    hv = adata.var["highly_variable"].to_numpy()
    keep_g = hv | adata.var_names.isin(forced)
    n_extra = int(keep_g.sum() - hv.sum())
    adata = adata[:, keep_g].copy()
    print(f"  HVG: {int(hv.sum())} HVG + {n_extra} forced markers -> {adata.n_vars} genes")
    print(f"  forced-in markers ({len(forced)} present): {sorted(forced)}")

    adata.obs["cell_id"] = adata.obs_names.astype(str)
    adata.obs["time_label"] = adata.obs["sample_id"].astype(str)
    adata.uns["dataset_id"] = DATASET_ID
    adata.uns["species"] = "mouse"
    adata.uns["build_timestamp"] = datetime.now().strftime("%Y%m%d_%H%M")
    adata.uns["annotation_note"] = (
        "Unlabeled mouse dataset (per-sample 10x). final_milestone_label_coarse added by "
        "marker-FM annotation using milestone_markers_mouse.yaml; MEF start at abs_day==0.")
    out_dir.joinpath("build_qc").mkdir(exist_ok=True)
    json.dump({"dataset_id": DATASET_ID, "species": "mouse", "n_cells": int(adata.n_obs),
               "n_vars": int(adata.n_vars), "samples": sorted(samples),
               "abs_day_counts": adata.obs["abs_day"].value_counts().sort_index().to_dict(),
               "runtime_s": round(_time.time() - T0, 1)},
              open(out_dir / "build_qc" / "build_summary.json", "w"), indent=2)
    print(f"[write] {out_h5ad}")
    adata.write_h5ad(out_h5ad)
    print(f"Done. cells={adata.n_obs} genes={adata.n_vars}")
    print("Next: marker-FM annotation with milestone_markers_mouse.yaml.")


if __name__ == "__main__":
    main()
