#!/usr/bin/env python3
"""Build the GSE298212 HVG2000 benchmark input (human blood chemical reprogramming).

GSE298212 (Peng et al. 2025, Cell Stem Cell, doi:10.1016/j.stem.2025.07.003) is a
fully-chemical reprogramming of human peripheral/cord blood cells toward hCiPS.
The deposited time course covers the Stage-1 reprogramming samples
(PBMC_EPC, S1D1, S1D3, S1D6, S1D8).

This builder mirrors scripts/build_gse175634_cardiac_author_input.py but for an
UNLABELED dataset: it produces the HVG2000 h5ad with raw counts + an abs_day axis;
the per-cell milestone labels (final_milestone_label_coarse) are added afterwards by
the marker-FM annotation step (benchmark/annotation/build_marker_seed_labels.py ->
build_trajectory_aware_labels.py), which reads the GSE298212 entry in
benchmark/annotation/milestone_markers.yaml. The reprogramming start (blood) is
assigned by source sample / abs_day==0, following the framework's hADSCs convention.

Expected raw layout under data/gse298212/ (standard 10x; adjust SAMPLES / loader if
your GEO download differs):
    data/gse298212/<sample>/{barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz}
  or a flat set of GSE298212_<sample>_{barcodes,features,matrix}.* files.

Run on Shirokane (the raw data is not staged locally):
    python scripts/build_gse298212_chemical_input.py [--overwrite]
"""
from __future__ import annotations

import argparse
import gc
import json
import os
import time as _time
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

DATASET_ID = "GSE298212"
DATA_FOLDER = "gse298212"
N_TOP_GENES = 2000
NORM_TARGET_SUM = 1e4
T0 = _time.time()

# Sample -> abs_day (Stage-1 reprogramming time course). Adjust to match the
# actual sample names in your GEO download.
SAMPLE_TO_ABS_DAY: dict[str, float] = {
    "PBMC_EPC": 0.0,   # blood starting population (trajectory root, time==0)
    "EPC":      0.0,
    "S1D1":     1.0,
    "S1D3":     3.0,
    "S1D6":     6.0,
    "S1D8":     8.0,
}

QC_PARAMS = {
    "min_genes_per_cell": 200,
    "max_genes_per_cell": 10_000,
    "min_counts_per_cell": 500,
    "max_pct_mito": 25.0,
    "min_cells_per_gene": 3,
}


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


def _abs_day_for(sample: str) -> float:
    s = str(sample)
    if s in SAMPLE_TO_ABS_DAY:
        return SAMPLE_TO_ABS_DAY[s]
    for key, day in SAMPLE_TO_ABS_DAY.items():
        if key.lower() in s.lower():
            return day
    return float("nan")


def _manual_10x(sd: Path) -> ad.AnnData:
    """Read a 10x triplet from a directory, robust to features-file format."""
    import gzip
    import scipy.io as sio
    import scipy.sparse as sp

    def pick(patterns):
        for pat in patterns:
            hits = list(sd.glob(pat))
            if hits:
                return str(hits[0])
        return None

    mtx = pick(["*matrix.mtx.gz", "*matrix.mtx"])
    bc = pick(["*barcodes.tsv.gz", "*barcodes.tsv"])
    feat = pick(["*features.tsv.gz", "*genes.tsv.gz", "*features.tsv", "*genes.tsv"])
    if not (mtx and bc and feat):
        raise FileNotFoundError(f"{sd}: missing one of matrix/barcodes/features")

    def lines(p):
        op = gzip.open if p.endswith(".gz") else open
        with op(p, "rt") as f:
            return [ln.rstrip("\n").split("\t") for ln in f]

    X = sp.csr_matrix(sio.mmread(mtx).astype(np.float32)).T   # cells x genes
    barcodes = [r[0] for r in lines(bc)]
    feats = lines(feat)
    symbols = [(r[1] if len(r) >= 2 else r[0]) for r in feats]
    a = ad.AnnData(X=X)
    a.obs_names = barcodes
    a.var_names = pd.Index(symbols)
    a.var_names_make_unique()
    return a


def load_samples(data_dir: Path) -> ad.AnnData:
    """Load each 10x sample, tag sample_id + abs_day, and concatenate."""
    import scanpy as sc

    sample_dirs = sorted([p for p in data_dir.iterdir()
                          if p.is_dir() and p.name != "_raw"
                          and list(p.glob("*matrix.mtx*"))])
    if not sample_dirs:
        raise FileNotFoundError(
            f"No per-sample 10x directories (with a matrix.mtx) found under {data_dir}. "
            f"Run scripts/download_new_datasets.sh first."
        )
    parts = []
    for sd in sample_dirs:
        sample = sd.name
        print(f"  loading sample {sample} from {sd}")
        try:
            a = sc.read_10x_mtx(sd, var_names="gene_symbols", make_unique=True)
        except Exception as e:
            print(f"    read_10x_mtx failed ({type(e).__name__}); manual triplet read")
            a = _manual_10x(sd)
        a.obs["sample_id"] = sample
        a.obs["abs_day"] = _abs_day_for(sample)
        a.obs_names = [f"{sample}_{bc}" for bc in a.obs_names]
        parts.append(a)
    adata = ad.concat(parts, join="outer", index_unique=None)
    adata.obs_names_make_unique()
    print(f"  concatenated {len(parts)} samples -> {adata.shape}")
    return adata


def qc_normalize_hvg(adata: ad.AnnData) -> ad.AnnData:
    import scanpy as sc
    n0 = adata.n_obs
    adata.var_names_make_unique()
    mito = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
    adata.obs["pct_counts_mt"] = (
        np.asarray(adata[:, mito].X.sum(1)).ravel()
        / np.asarray(adata.X.sum(1)).ravel().clip(1) * 100
    )
    adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(1)).ravel()
    adata.obs["n_counts"] = np.asarray(adata.X.sum(1)).ravel()
    keep = ((adata.obs["n_genes"] >= QC_PARAMS["min_genes_per_cell"])
            & (adata.obs["n_genes"] <= QC_PARAMS["max_genes_per_cell"])
            & (adata.obs["n_counts"] >= QC_PARAMS["min_counts_per_cell"])
            & (adata.obs["pct_counts_mt"] <= QC_PARAMS["max_pct_mito"]))
    adata = adata[keep].copy()
    gc = np.asarray((adata.X > 0).sum(0)).ravel()
    adata = adata[:, gc >= QC_PARAMS["min_cells_per_gene"]].copy()
    print(f"  QC: {n0} -> {adata.n_obs} cells, {adata.n_vars} genes")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=NORM_TARGET_SUM)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_TOP_GENES, flavor="seurat")
    forced = _force_include_marker_genes(adata, find_project_root(), "milestone_markers.yaml", DATASET_ID)
    hv = adata.var["highly_variable"].to_numpy()
    keep_g = hv | adata.var_names.isin(forced)
    n_extra = int(keep_g.sum() - hv.sum())
    adata = adata[:, keep_g].copy()
    print(f"  HVG: {int(hv.sum())} HVG + {n_extra} forced markers -> {adata.n_vars} genes")
    print(f"  forced-in markers ({len(forced)} present): {sorted(forced)}")
    return adata


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--project-root", type=Path, default=None)
    args = ap.parse_args()

    root = args.project_root or find_project_root()
    data_dir = root / "data" / DATA_FOLDER
    out_dir = root / "benchmark" / "inputs" / "gse298212_marker_fm_transition_silver_hvg2000"
    out_h5ad = out_dir / "GSE298212_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"
    print("=" * 64)
    print(f"Build {DATASET_ID} HVG2000 benchmark input")
    print(f"  data dir : {data_dir}")
    print(f"  output   : {out_h5ad}")
    print("=" * 64)
    if out_h5ad.exists() and not args.overwrite:
        raise FileExistsError(f"{out_h5ad} exists; pass --overwrite.")
    if not data_dir.is_dir():
        raise FileNotFoundError(f"raw data dir not found: {data_dir}")
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = load_samples(data_dir)
    if adata.obs["abs_day"].isna().any():
        bad = sorted(adata.obs.loc[adata.obs["abs_day"].isna(), "sample_id"].unique())
        raise ValueError(f"samples without an abs_day mapping: {bad}; edit SAMPLE_TO_ABS_DAY.")
    adata = qc_normalize_hvg(adata)

    adata.obs["cell_id"] = adata.obs_names.astype(str)
    adata.obs["time_label"] = adata.obs["sample_id"].astype(str)
    adata.uns["dataset_id"] = DATASET_ID
    adata.uns["build_timestamp"] = datetime.now().strftime("%Y%m%d_%H%M")
    adata.uns["annotation_note"] = (
        "Unlabeled dataset. final_milestone_label_coarse is added by the marker-FM "
        "annotation step using the GSE298212 entry in milestone_markers.yaml; the blood "
        "start is assigned by abs_day==0 (source sample), per the framework hADSCs convention."
    )
    out_dir.joinpath("build_qc").mkdir(exist_ok=True)
    summary = {"dataset_id": DATASET_ID, "n_cells": int(adata.n_obs), "n_vars": int(adata.n_vars),
               "abs_day_counts": adata.obs["abs_day"].value_counts().sort_index().to_dict(),
               "sample_counts": adata.obs["sample_id"].value_counts().to_dict(),
               "qc_params": QC_PARAMS, "runtime_s": round(_time.time() - T0, 1)}
    json.dump(summary, open(out_dir / "build_qc" / "build_summary.json", "w"), indent=2)
    gc.collect()
    print(f"[write] {out_h5ad}")
    adata.write_h5ad(out_h5ad)
    print(f"Done ({summary['runtime_s']}s). cells={adata.n_obs} genes={adata.n_vars}")
    print("Next: marker-FM annotation (run_v21_new_datasets_preprocess.sh handles this).")


if __name__ == "__main__":
    main()
