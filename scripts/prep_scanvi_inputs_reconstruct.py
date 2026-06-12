#!/usr/bin/env python
"""Build scVI-ready scANVI inputs from the log-normalized full-gene .h5ad files.

The `*_raw_full_gene_benchmark_input.h5ad` files store a log-normalized X
(normalize_total(target_sum=1e4) -> log1p) and drop the raw counts, but keep the
per-cell library size in obs['total_counts']. That normalization is exactly
invertible element-wise (verified on real data: reconstructed library sizes
correlate 1.0000 with obs['total_counts']):

    counts_ij = round( expm1(X_ij) / 1e4 * total_counts_i )

This script: (1) intersects the gene space across the 3 full-gene files,
(2) selects 2000 HVGs on the reference within that intersection, (3) reconstructs
integer counts on those genes for each dataset, (4) joins the per-cell label from
the matching HVG2000 file by barcode, (5) writes scVI-ready inputs (X = counts,
layers['counts'] = counts) where the `server` profile in config.yaml expects them.

Run as a CPU job. Depends only on scanpy/anndata/scipy (the traj_env).
"""
from __future__ import annotations

import os
import numpy as np
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

TARGET = 1e4
OUT = os.environ.get("PREP_OUT_DIR", "results/annotation_branch/inputs")
N_HVG = int(os.environ.get("PREP_N_HVG", "2000"))

# (output_name, full_gene_h5ad, label_h5ad, label_key)
SPECS = [
    ("gse242424_ref_fullgene_labelled.h5ad",
     "data/processed/gse242424_human/GSE242424_raw_full_gene_benchmark_input.h5ad",
     "benchmark/inputs/gse242424_author_cluster_matched/GSE242424_author_cluster_matched_HVG2000_benchmark_input.h5ad",
     "author_cluster_label"),
    ("gse178325_query_fullgene_labelled.h5ad",
     "data/processed/gse178325_human/GSE178325_0618_raw_full_gene_benchmark_input.h5ad",
     "benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad",
     "final_milestone_label_coarse"),
    ("gse230659_query_fullgene_labelled.h5ad",
     "data/processed/gse230659_human/GSE230659_0618_raw_full_gene_benchmark_input.h5ad",
     "benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad",
     "final_milestone_label_coarse"),
]


def reconstruct_counts(A):
    """Invert normalize_total(1e4)->log1p back to integer counts (sparse-safe)."""
    X = A.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    X = X.astype(np.float64).copy()
    X.data = np.expm1(X.data)
    tc = A.obs["total_counts"].values.astype(np.float64)
    X = sp.diags(tc / TARGET) @ X
    X.data = np.rint(X.data)
    X.eliminate_zeros()
    return X.astype(np.float32).tocsr()


def vnames(fp):
    a = ad.read_h5ad(fp, backed="r")
    v = list(a.var_names)
    a.file.close()
    return v


def main():
    os.makedirs(OUT, exist_ok=True)
    ref_vars = vnames(SPECS[0][1])
    inter = set(ref_vars)
    for spec in SPECS[1:]:
        inter &= set(vnames(spec[1]))
    inter = [g for g in ref_vars if g in inter]
    print("intersection genes:", len(inter), flush=True)

    ref = ad.read_h5ad(SPECS[0][1])[:, inter].copy()
    sc.pp.highly_variable_genes(ref, n_top_genes=N_HVG)
    hvg = ref.var_names[ref.var["highly_variable"].values].tolist()
    print("HVG:", len(hvg), flush=True)
    del ref

    for out_name, fg, lab, key in SPECS:
        A = ad.read_h5ad(fg)[:, hvg].copy()
        counts = reconstruct_counts(A)
        A.X = counts.copy()
        A.layers["counts"] = counts
        L = ad.read_h5ad(lab, backed="r")
        common = A.obs_names.intersection(L.obs.index)
        A = A[common].copy()
        A.obs[key] = L.obs.loc[A.obs_names, key].astype(str).values
        L.file.close()
        p = os.path.join(OUT, out_name)
        A.write_h5ad(p)
        int_ok = bool(np.allclose(A.layers["counts"].data, np.round(A.layers["counts"].data)))
        print("wrote", p, "cells", A.n_obs, "genes", A.n_vars, key,
              A.obs[key].nunique(), "int", int_ok, flush=True)
        del A
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
