#!/usr/bin/env python3
"""
Export the GSE178325 scGPT-annotated HVG2000 benchmark input.

This mirrors the existing GSE230659 artifact:
benchmark/inputs/gse230659_scgpt_hvg2000/
  GSE230659_scGPT_annotated_HVG2000_benchmark_input.h5ad

The output keeps the same obs/var/obsm preservation policy: subset genes only,
keep scGPT annotation columns, keep X_scGPT/X_umap, and add the same benchmark
HVG metadata keys.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp


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
    return here.parents[2]


def _parse_args() -> argparse.Namespace:
    root = _find_project_root()
    parser = argparse.ArgumentParser(
        description="Build GSE178325 scGPT-annotated HVG2000 benchmark input"
    )
    parser.add_argument(
        "--input-h5ad",
        default=str(
            root
            / "benchmark/results/scgpt/gse178325_0618/full/adata_scgpt_annotated.h5ad"
        ),
    )
    parser.add_argument(
        "--output-h5ad",
        default=str(
            root
            / "benchmark/inputs/gse178325_scgpt_hvg2000/"
            / "GSE178325_scGPT_annotated_HVG2000_benchmark_input.h5ad"
        ),
    )
    parser.add_argument("--n-top-genes", type=int, default=2000)
    return parser.parse_args()


def _sparse_mean_var(X) -> tuple[np.ndarray, np.ndarray]:
    if sp.issparse(X):
        X = X.tocsr()
        mean = np.asarray(X.mean(axis=0)).ravel()
        mean_sq = np.asarray(X.multiply(X).mean(axis=0)).ravel()
    else:
        X = np.asarray(X)
        mean = X.mean(axis=0)
        mean_sq = np.square(X).mean(axis=0)
    var = np.maximum(mean_sq - np.square(mean), 0.0)
    return mean.astype(np.float64), var.astype(np.float64)


def _dispersion_scores(X) -> np.ndarray:
    mean, var = _sparse_mean_var(X)
    with np.errstate(divide="ignore", invalid="ignore"):
        score = var / mean
    score[~np.isfinite(score)] = -np.inf
    score[mean <= 0] = -np.inf
    return score


def main() -> None:
    args = _parse_args()
    input_h5ad = Path(args.input_h5ad)
    output_h5ad = Path(args.output_h5ad)
    output_h5ad.parent.mkdir(parents=True, exist_ok=True)

    print(f"[hvg2000] Loading {input_h5ad}")
    adata = ad.read_h5ad(input_h5ad)
    if not adata.var_names.is_unique:
        raise ValueError(
            "Input var_names are not unique. Rebuild the scGPT full artifact "
            "from a deduplicated raw/full-gene h5ad before exporting HVG2000."
        )
    if adata.n_vars < args.n_top_genes:
        raise ValueError(f"Input has only {adata.n_vars} genes; cannot select {args.n_top_genes}.")

    scores = _dispersion_scores(adata.X)
    hvg_idx = np.argsort(scores)[::-1][: args.n_top_genes]
    hvg_idx = np.sort(hvg_idx)

    out = adata[:, hvg_idx].copy()
    out.var["highly_variable"] = True
    out.var["benchmark_hvg_score"] = scores[hvg_idx].astype(np.float32)
    out.uns["benchmark_hvg_n_top_genes"] = int(args.n_top_genes)
    out.uns["benchmark_hvg_selection_method"] = (
        "Top genes by finite sparse mean/variance dispersion from "
        "scGPT-annotated full-gene AnnData."
    )
    out.uns["benchmark_hvg_source"] = str(input_h5ad)
    out.uns["benchmark_input_id"] = "GSE178325_scGPT_annotated_HVG2000"
    out.uns["benchmark_input_label"] = "GSE178325 scGPT-annotated HVG2000 benchmark input"

    out.write_h5ad(output_h5ad)
    print("[hvg2000] Wrote")
    print(f"  {output_h5ad}")
    print(f"  shape: {out.n_obs:,} cells x {out.n_vars:,} genes")
    print(f"  obs columns: {list(out.obs.columns)}")
    print(f"  var columns: {list(out.var.columns)}")
    print(f"  obsm keys: {list(out.obsm.keys())}")


if __name__ == "__main__":
    main()
