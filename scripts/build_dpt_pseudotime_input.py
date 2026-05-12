#!/usr/bin/env python3
"""Build a scTimeBench-style DPT pseudotime benchmark input.

The pseudotime computation follows scTimeBench's default Pseudotime
preprocessor: Scanpy neighbors -> diffusion map -> earliest-time root -> DPT.
The continuous DPT values are then discretized into ordered equal-cell-count
bins so they can be used as the benchmark time axis for D-F scenarios.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc


def _find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
        raise FileNotFoundError(f"TRAJ_PROJECT_ROOT={env!r} does not exist.")
    here = Path(__file__).resolve()
    for candidate in [here.parent, *here.parents]:
        if (candidate / "benchmark").exists() and (candidate / "data").exists():
            return candidate
    return here.parents[1]


def _parse_args() -> argparse.Namespace:
    root = _find_project_root()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-h5ad",
        default=str(
            root
            / "benchmark/inputs/gse230659_scgpt_hvg2000/"
            / "GSE230659_scGPT_annotated_HVG2000_benchmark_input.h5ad"
        ),
    )
    parser.add_argument(
        "--output-h5ad",
        default=str(
            root
            / "benchmark/inputs/gse230659_scgpt_hvg2000/"
            / "GSE230659_scGPT_annotated_HVG2000_DPT_pseudotime_bins.h5ad"
        ),
    )
    parser.add_argument("--observed-time-key", default="abs_day")
    parser.add_argument("--pseudotime-key", default="dpt_pseudotime")
    parser.add_argument("--bin-key", default="dpt_pseudotime_bin")
    parser.add_argument("--bin-numeric-key", default="dpt_pseudotime_bin_numeric")
    parser.add_argument("--num-bins", type=int, default=None)
    parser.add_argument(
        "--representation",
        choices=["hvg", "zheng_hvg", "pca", "obsm"],
        default="obsm",
        help="Representation used before neighbors/diffmap. 'obsm' uses --obsm-key.",
    )
    parser.add_argument("--obsm-key", default="X_scGPT")
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--n-top-genes", type=int, default=1000)
    parser.add_argument("--pca-components", type=int, default=50)
    parser.add_argument("--metadata-json", default=None)
    return parser.parse_args()


def _prepare_for_dpt(adata, args):
    if args.representation == "obsm":
        if args.obsm_key not in adata.obsm:
            raise RuntimeError(
                f"Requested obsm[{args.obsm_key!r}], available keys: {list(adata.obsm.keys())}"
            )
        X = np.asarray(adata.obsm[args.obsm_key], dtype=np.float32)
        return ad.AnnData(X=X, obs=adata.obs.copy())

    if args.representation == "hvg":
        work = adata.copy()
        sc.pp.highly_variable_genes(work, n_top_genes=args.n_top_genes, inplace=True)
        return work[:, work.var.highly_variable].copy()

    if args.representation == "zheng_hvg":
        return sc.pp.recipe_zheng17(adata, n_top_genes=args.n_top_genes, copy=True)

    work = adata.copy()
    sc.tl.pca(work, n_comps=args.pca_components)
    return ad.AnnData(X=work.obsm["X_pca"], obs=adata.obs.copy())


def _assign_equal_cell_bins(values: np.ndarray, num_bins: int) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    bins = np.empty(values.shape[0], dtype=np.int32)
    for bin_id, idx in enumerate(np.array_split(order, num_bins)):
        bins[idx] = bin_id
    return bins


def _default_splits(num_bins: int) -> dict:
    bins = list(range(num_bins))
    d_holdout = [b for b in (num_bins // 4, num_bins // 2, (3 * num_bins) // 4)]
    d_holdout = sorted(set(b for b in d_holdout if 0 < b < num_bins - 1))
    e_start = max(2, int(np.floor(num_bins * 2 / 3)))
    e_holdout = bins[e_start:]
    f_internal = d_holdout[:2]
    f_holdout = sorted(set(f_internal + e_holdout))
    return {
        "D": {
            "name": "Pseudotime interpolation",
            "train_times": [float(b) for b in bins if b not in d_holdout],
            "heldout_times": [float(b) for b in d_holdout],
        },
        "E": {
            "name": "Pseudotime extrapolation",
            "train_times": [float(b) for b in bins if b not in e_holdout],
            "heldout_times": [float(b) for b in e_holdout],
        },
        "F": {
            "name": "Pseudotime interpolation + extrapolation",
            "train_times": [float(b) for b in bins if b not in f_holdout],
            "heldout_times": [float(b) for b in f_holdout],
        },
    }


def main() -> None:
    args = _parse_args()
    input_h5ad = Path(args.input_h5ad)
    output_h5ad = Path(args.output_h5ad)
    metadata_json = (
        Path(args.metadata_json)
        if args.metadata_json
        else output_h5ad.with_suffix(".metadata.json")
    )
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    metadata_json.parent.mkdir(parents=True, exist_ok=True)

    print(f"[dpt] Loading {input_h5ad}")
    adata = ad.read_h5ad(input_h5ad)
    if args.observed_time_key not in adata.obs.columns:
        raise RuntimeError(f"obs[{args.observed_time_key!r}] is missing.")

    num_bins = args.num_bins
    if num_bins is None:
        num_bins = int(adata.obs[args.observed_time_key].astype(float).nunique())
    if num_bins < 3:
        raise ValueError("--num-bins must be at least 3.")

    print(f"[dpt] Preparing representation: {args.representation}")
    work = _prepare_for_dpt(adata, args)
    print("[dpt] Computing neighbors, diffusion map, and DPT")
    use_rep = "X" if args.representation in ("obsm", "pca") else None
    sc.pp.neighbors(work, n_neighbors=args.n_neighbors, use_rep=use_rep)
    sc.tl.diffmap(work)

    observed = adata.obs[args.observed_time_key].astype(float)
    earliest_tp = float(observed.min())
    root_idx = observed[observed == earliest_tp].index[0]
    work.uns["iroot"] = work.obs_names.get_loc(root_idx)
    sc.tl.dpt(work)

    pseudotime = np.asarray(work.obs["dpt_pseudotime"], dtype=np.float32)
    if not np.all(np.isfinite(pseudotime)):
        raise RuntimeError("DPT produced non-finite pseudotime values.")
    bins = _assign_equal_cell_bins(pseudotime, num_bins=num_bins)

    adata.obs[args.pseudotime_key] = pseudotime
    adata.obs[args.bin_numeric_key] = bins.astype(float)
    adata.obs[args.bin_key] = np.array([f"PT_{b:02d}" for b in bins], dtype=object)
    adata.obs["scTimeBench_timepoint"] = adata.obs[args.bin_numeric_key]
    adata.uns["dpt_pseudotime"] = {
        "source": "trajectory/scripts/build_dpt_pseudotime_input.py",
        "sctimebench_alignment": (
            "Scanpy neighbors -> diffmap -> earliest observed time root -> dpt; "
            "continuous DPT discretized into equal-cell-count bins."
        ),
        "observed_time_key": args.observed_time_key,
        "pseudotime_key": args.pseudotime_key,
        "bin_key": args.bin_key,
        "bin_numeric_key": args.bin_numeric_key,
        "num_bins": int(num_bins),
        "representation": args.representation,
        "obsm_key": args.obsm_key if args.representation == "obsm" else None,
        "n_neighbors": int(args.n_neighbors),
        "root_observed_time": earliest_tp,
        "root_obs_name": str(root_idx),
    }

    print(f"[dpt] Writing {output_h5ad}")
    adata.write_h5ad(output_h5ad)

    counts = (
        adata.obs.groupby(args.bin_numeric_key, observed=False)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    metadata = {
        "input_h5ad": str(input_h5ad),
        "output_h5ad": str(output_h5ad),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "observed_time_key": args.observed_time_key,
        "pseudotime_key": args.pseudotime_key,
        "bin_key": args.bin_key,
        "bin_numeric_key": args.bin_numeric_key,
        "num_bins": int(num_bins),
        "representation": args.representation,
        "obsm_key": args.obsm_key if args.representation == "obsm" else None,
        "n_neighbors": int(args.n_neighbors),
        "root_observed_time": earliest_tp,
        "root_obs_name": str(root_idx),
        "bin_counts": [
            {
                "bin": float(row[args.bin_numeric_key]),
                "n_cells": int(row["n_cells"]),
            }
            for _, row in counts.iterrows()
        ],
        "scenario_splits": _default_splits(int(num_bins)),
    }
    metadata_json.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"[dpt] Wrote metadata {metadata_json}")


if __name__ == "__main__":
    main()
