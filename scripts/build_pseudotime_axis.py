#!/usr/bin/env python
"""Build the pseudotime benchmark input for Scenarios D/E/F.

Loads an existing observed-time HVG2000 benchmark h5ad, computes a Scanpy DPT
pseudotime axis (rooted at the earliest observed time point) using the project's
existing scTimeBench-style preprocessors, discretizes it into equal-cell-count
bins, and writes a new ``*_pseudotime_benchmark_input.h5ad`` that carries:

    obs['dpt_pseudotime']             continuous DPT
    obs['dpt_pseudotime_bin']         string bins  PT_00..PT_NN
    obs['dpt_pseudotime_bin_numeric'] float bin index 0..N-1   <- scenario time axis

The D/E/F runtime configs set ``time_key: dpt_pseudotime_bin_numeric`` and the
bin-based ``scenario_params`` from benchmark/configs/scenario_pseudotime.yaml, so
the unchanged dispatcher filters train/heldout on these bins.

Run on Shirokane (full data); locally only on a small subsample for smoke testing.

    python scripts/build_pseudotime_axis.py \
        --input  benchmark/inputs/.../GSE230659_..._HVG2000_benchmark_input.h5ad \
        --output benchmark/inputs/.../GSE230659_..._HVG2000_pseudotime_benchmark_input.h5ad \
        --time-key abs_day --num-bins 15
"""
from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad

from benchmark.shared.dataset.preprocessors.pseudotime import (
    DPTPseudotime,
    RoundPseudotimeToBins,
)


def build(input_path, output_path, time_key="abs_day", num_bins=15,
          n_top_genes=1000, n_neighbors=15):
    import numpy as np  # noqa: F401  (kept for finite checks / future use)
    import scanpy as sc

    adata = ad.read_h5ad(input_path)
    if time_key not in adata.obs.columns:
        raise SystemExit(f"obs[{time_key!r}] required for DPT root selection; "
                         f"have {list(adata.obs.columns)[:20]}")

    dataset_dict = {"time_key": time_key}

    # Build a scale-robust PCA representation for DPT if one is not present.
    # Benchmark inputs are inconsistent in X scale: some store count/library-scale
    # X (which makes scanpy's seurat-flavor HVG expm1 OVERFLOW to infinity, e.g.
    # GSE230659), others store log1p-scale X (e.g. GSE178325). We therefore avoid
    # the seurat-HVG path entirely: normalize+log only when the data looks
    # count-scale, then PCA, so the DPT manifold is consistent across datasets.
    if "X_pca" not in adata.obsm:
        work = adata.copy()
        xmax = float(work.X.max())
        if xmax > 30.0:  # count- or library-scale, not log1p
            sc.pp.normalize_total(work, target_sum=1e4)
            sc.pp.log1p(work)
            print(f"[info] X looked count-scale (max={xmax:.1f}); applied "
                  f"normalize_total(1e4)+log1p before PCA.")
        n_comps = int(min(50, work.n_vars - 1, work.n_obs - 1))
        sc.pp.pca(work, n_comps=n_comps)
        adata.obsm["X_pca"] = work.obsm["X_pca"]
        del work

    dpt = DPTPseudotime(
        dataset_dict,
        preprocess_type="obsm",
        observed_time_key=time_key,
        n_neighbors=n_neighbors,
        obsm_key="X_pca",
        copy_to_timepoint=False,
    )
    adata = dpt.preprocess(adata)

    binner = RoundPseudotimeToBins(
        dataset_dict, num_bins=num_bins, copy_to_timepoint=False,
    )
    adata = binner.preprocess(adata)

    adata.uns["pseudotime_axis"] = {
        "source_method": "scanpy_dpt",
        "preprocess_type": "obsm_pca",
        "num_bins": int(num_bins),
        "time_key_observed": time_key,
        "bin_numeric_key": "dpt_pseudotime_bin_numeric",
        "built_by": "scripts/build_pseudotime_axis.py",
    }
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)

    counts = adata.obs["dpt_pseudotime_bin_numeric"].value_counts().sort_index()
    print(f"[OK] wrote {output_path}")
    print(f"     cells={adata.n_obs} bins={num_bins} "
          f"dpt range=[{adata.obs['dpt_pseudotime'].min():.3f}, "
          f"{adata.obs['dpt_pseudotime'].max():.3f}]")
    print(f"     per-bin cell counts:\n{counts.to_string()}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--time-key", default="abs_day")
    ap.add_argument("--num-bins", type=int, default=15)
    ap.add_argument("--n-top-genes", type=int, default=1000)
    ap.add_argument("--n-neighbors", type=int, default=15)
    args = ap.parse_args()
    build(args.input, args.output, args.time_key, args.num_bins,
          args.n_top_genes, args.n_neighbors)


if __name__ == "__main__":
    main()
