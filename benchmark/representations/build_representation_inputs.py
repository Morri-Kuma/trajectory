"""
build_representation_inputs.py
==============================

Work-Plan Step 2: build the common 50-dimensional representation inputs that the
trajectory models consume via ``adata.obsm["X_rep"]``.

Four arms (and ONLY these):

    rep_hvg_pca50              raw counts -> libsize norm -> log1p -> top 2000 HVG -> PCA50
    rep_geneformer_cls_pca50   X_geneformer_cls -> StandardScaler -> PCA50
    rep_scgpt_cls_pca50        X_scgpt_cls      -> StandardScaler -> PCA50
    rep_scfoundation_pca50     X_scfoundation   -> StandardScaler -> PCA50

Anti-leakage rule (Work-Plan Step 2 / Risk 2)
---------------------------------------------
For *formal* runs the scaler and PCA must be fit on TRAINING cells only for each
scenario, then applied to all cells. The fit scope is recorded in
``reducer_fit_scope`` (``train_only`` or ``all_cells_exploratory``). All-cell
fitting is allowed only as an explicitly-labelled exploratory mode and is never
used silently for formal configs.

The pure functions here (``build_hvg_pca``, ``build_scfm_pca``) are import-safe
without anndata/scanpy so they can be unit-tested directly on numpy arrays.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

from .metadata import CODE_VERSION, utc_now_iso
from .references import REP_TO_RAW_OBSM, REPRESENTATION_ARMS

DEFAULT_N_COMPONENTS = 50
DEFAULT_N_HVG = 2000

TRAIN_ONLY = "train_only"
ALL_CELLS_EXPLORATORY = "all_cells_exploratory"


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class RepresentationResult:
    X_rep: np.ndarray
    representation_id: str
    reducer_fit_scope: str
    n_components: int
    metadata: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------
# Core numpy/sklearn building blocks (test-friendly)
# ---------------------------------------------------------------------------

def _library_size_normalize(counts: np.ndarray, target_sum: float = 1e4) -> np.ndarray:
    """Library-size normalize raw counts to ``target_sum`` per cell, then log1p.

    Matches the reference-paper HVG preprocessing:
        raw counts -> normalize_total(target_sum) -> log1p
    """
    counts = np.asarray(counts, dtype=np.float64)
    totals = counts.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.where(totals == 0, 0.0, target_sum / totals)
    return np.log1p(counts * scale)


def _seurat_dispersion_hvg(lognorm: np.ndarray, n_top: int) -> np.ndarray:
    """
    Select ``n_top`` highly variable genes by normalized dispersion on
    log-normalized data (Seurat/scanpy ``flavor='seurat'`` style, binned by
    mean expression). Returns the selected column indices (sorted ascending).

    Used only as a fallback when scanpy is unavailable; the CLI prefers
    ``scanpy.pp.highly_variable_genes`` for exact parity with the reference.
    """
    mean = lognorm.mean(axis=0)
    var = lognorm.var(axis=0)
    n_genes = lognorm.shape[1]
    n_top = min(n_top, n_genes)
    with np.errstate(divide="ignore", invalid="ignore"):
        dispersion = np.where(mean > 0, var / mean, 0.0)
    # Bin genes by mean into 20 bins; z-score dispersion within each bin.
    n_bins = 20
    df_mean = mean.copy()
    # rank-based binning to avoid empty bins
    order = np.argsort(df_mean)
    bins = np.zeros(n_genes, dtype=int)
    bins[order] = (np.arange(n_genes) * n_bins // max(n_genes, 1)).clip(0, n_bins - 1)
    norm_disp = np.zeros(n_genes, dtype=float)
    for b in range(n_bins):
        mask = bins == b
        if not np.any(mask):
            continue
        d = dispersion[mask]
        mu, sd = d.mean(), d.std()
        sd = sd if sd > 0 else 1.0
        norm_disp[mask] = (d - mu) / sd
    top_idx = np.argsort(norm_disp)[::-1][:n_top]
    return np.sort(top_idx)


def build_hvg_pca(
    counts_all: np.ndarray,
    train_mask: Optional[np.ndarray],
    *,
    n_components: int = DEFAULT_N_COMPONENTS,
    n_hvg: int = DEFAULT_N_HVG,
    seed: int = 0,
    hvg_indices: Optional[np.ndarray] = None,
) -> RepresentationResult:
    """
    Build rep_hvg_pca50 from raw counts.

    Pipeline: library-size normalize -> log1p -> select ``n_hvg`` HVGs -> PCA.

    Anti-leakage: HVG selection and PCA are fit on ``train_mask`` rows only when
    ``train_mask`` is provided (``reducer_fit_scope='train_only'``); otherwise on
    all cells (``all_cells_exploratory``). The fitted transform is applied to all
    cells.

    Parameters
    ----------
    counts_all : (n_cells, n_genes) raw counts for ALL cells.
    train_mask : boolean array selecting training cells, or None for exploratory.
    hvg_indices : optionally supply pre-computed HVG column indices (e.g. from
        scanpy) to use instead of the internal dispersion fallback.
    """
    from sklearn.decomposition import PCA

    counts_all = np.asarray(counts_all, dtype=np.float64)
    n_cells, n_genes = counts_all.shape
    if train_mask is None:
        fit_scope = ALL_CELLS_EXPLORATORY
        fit_rows = np.ones(n_cells, dtype=bool)
    else:
        fit_scope = TRAIN_ONLY
        fit_rows = np.asarray(train_mask, dtype=bool)
        if fit_rows.sum() == 0:
            raise ValueError("train_mask selects zero cells.")

    lognorm_all = _library_size_normalize(counts_all)

    if hvg_indices is None:
        hvg_indices = _seurat_dispersion_hvg(lognorm_all[fit_rows], n_hvg)
        hvg_source = "internal_dispersion_fallback"
    else:
        hvg_indices = np.sort(np.asarray(hvg_indices, dtype=int))
        hvg_source = "provided"

    lognorm_hvg = lognorm_all[:, hvg_indices]

    n_comp = int(min(n_components, fit_rows.sum() - 1, lognorm_hvg.shape[1]))
    if n_comp < 2:
        raise ValueError(
            f"Not enough cells/genes for PCA: n_components={n_comp}."
        )
    pca = PCA(n_components=n_comp, svd_solver="randomized", random_state=seed)
    pca.fit(lognorm_hvg[fit_rows])
    X_rep = pca.transform(lognorm_hvg).astype(np.float32)

    meta = {
        "representation_id": "rep_hvg_pca50",
        "input_space": "expression_counts",
        "output_space": "representation",
        "preprocessing": "libsize_norm_1e4 -> log1p -> top_HVG -> PCA",
        "n_hvg": int(len(hvg_indices)),
        "hvg_source": hvg_source,
        "n_components": n_comp,
        "reducer_fit_scope": fit_scope,
        "n_train_cells_fit": int(fit_rows.sum()),
        "n_total_cells": int(n_cells),
        "code_version": CODE_VERSION,
        "date_created": utc_now_iso(),
    }
    return RepresentationResult(X_rep, "rep_hvg_pca50", fit_scope, n_comp, meta)


def build_scfm_pca(
    emb_all: np.ndarray,
    train_mask: Optional[np.ndarray],
    representation_id: str,
    *,
    n_components: int = DEFAULT_N_COMPONENTS,
    seed: int = 0,
) -> RepresentationResult:
    """
    Build a scFM representation: raw embedding -> StandardScaler -> PCA.

    Anti-leakage: scaler + PCA fit on ``train_mask`` rows only when provided
    (``train_only``); else all cells (``all_cells_exploratory``). Transform is
    applied to all cells.
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    emb_all = np.asarray(emb_all, dtype=np.float64)
    n_cells, emb_dim = emb_all.shape
    if train_mask is None:
        fit_scope = ALL_CELLS_EXPLORATORY
        fit_rows = np.ones(n_cells, dtype=bool)
    else:
        fit_scope = TRAIN_ONLY
        fit_rows = np.asarray(train_mask, dtype=bool)
        if fit_rows.sum() == 0:
            raise ValueError("train_mask selects zero cells.")

    scaler = StandardScaler()
    scaler.fit(emb_all[fit_rows])
    scaled_all = scaler.transform(emb_all)

    n_comp = int(min(n_components, fit_rows.sum() - 1, emb_dim))
    if n_comp < 2:
        raise ValueError(f"Not enough cells/dims for PCA: n_components={n_comp}.")
    pca = PCA(n_components=n_comp, svd_solver="randomized", random_state=seed)
    pca.fit(scaled_all[fit_rows])
    X_rep = pca.transform(scaled_all).astype(np.float32)

    meta = {
        "representation_id": representation_id,
        "input_space": "scfm_embedding",
        "output_space": "representation",
        "preprocessing": "StandardScaler -> PCA",
        "raw_embedding_dim": int(emb_dim),
        "n_components": n_comp,
        "reducer_fit_scope": fit_scope,
        "n_train_cells_fit": int(fit_rows.sum()),
        "n_total_cells": int(n_cells),
        "code_version": CODE_VERSION,
        "date_created": utc_now_iso(),
    }
    return RepresentationResult(X_rep, representation_id, fit_scope, n_comp, meta)


# ---------------------------------------------------------------------------
# AnnData-level orchestration (CLI)
# ---------------------------------------------------------------------------

def _train_mask_from_times(adata, time_key: str, train_times) -> np.ndarray:
    if time_key not in adata.obs.columns:
        raise ValueError(
            f"time_key '{time_key}' not in adata.obs (columns: "
            f"{list(adata.obs.columns)})."
        )
    cell_times = adata.obs[time_key].to_numpy().astype(float)
    train_set = {float(t) for t in train_times}
    mask = np.array([t in train_set for t in cell_times], dtype=bool)
    if mask.sum() == 0:
        raise ValueError(
            f"No cells matched train_times {sorted(train_set)} on '{time_key}'."
        )
    return mask


def build_representation(
    adata,
    representation_id: str,
    *,
    train_mask: Optional[np.ndarray] = None,
    n_components: int = DEFAULT_N_COMPONENTS,
    n_hvg: int = DEFAULT_N_HVG,
    seed: int = 0,
    use_scanpy_hvg: bool = True,
) -> RepresentationResult:
    """Build one representation arm from an AnnData object."""
    if representation_id not in REPRESENTATION_ARMS:
        raise ValueError(
            f"representation_id must be one of {REPRESENTATION_ARMS}, "
            f"got {representation_id!r}."
        )

    if representation_id == "rep_hvg_pca50":
        import scipy.sparse as sp

        counts = adata.X.toarray() if sp.issparse(adata.X) else np.asarray(adata.X)
        hvg_indices = None
        if use_scanpy_hvg:
            hvg_indices = _scanpy_hvg_indices(adata, n_hvg, train_mask)
        return build_hvg_pca(
            counts, train_mask, n_components=n_components, n_hvg=n_hvg,
            seed=seed, hvg_indices=hvg_indices,
        )

    raw_obsm = REP_TO_RAW_OBSM[representation_id]
    if raw_obsm not in adata.obsm:
        raise ValueError(
            f"adata.obsm['{raw_obsm}'] required for {representation_id} not found. "
            f"Run extract_scfm_embeddings.py first. "
            f"Available obsm keys: {list(adata.obsm.keys())}."
        )
    emb = np.asarray(adata.obsm[raw_obsm])
    return build_scfm_pca(
        emb, train_mask, representation_id, n_components=n_components, seed=seed,
    )


def _scanpy_hvg_indices(adata, n_hvg, train_mask):
    """Use scanpy HVG (seurat flavor) on libsize-norm/log1p data if available.

    Fit on training cells only when train_mask provided (anti-leakage). Returns
    None if scanpy is unavailable so the caller falls back to the numpy method.
    """
    try:
        import scanpy as sc
    except Exception:
        return None
    import scipy.sparse as sp

    sub = adata[train_mask].copy() if train_mask is not None else adata.copy()
    # Operate on a counts copy to avoid mutating caller state.
    sub = sub.copy()
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    sc.pp.highly_variable_genes(sub, n_top_genes=min(n_hvg, sub.n_vars),
                                flavor="seurat")
    hv = np.where(sub.var["highly_variable"].to_numpy())[0]
    return np.sort(hv)


def run(args) -> Path:
    try:
        import anndata
    except Exception as exc:  # pragma: no cover
        raise RuntimeError(f"anndata is required: {exc}") from exc

    adata = anndata.read_h5ad(args.input_h5ad)

    train_mask = None
    if args.fit_scope == TRAIN_ONLY:
        if not args.train_times:
            raise ValueError(
                "--fit-scope train_only requires --train-times (the scenario "
                "training timepoints). For exploratory all-cell fitting pass "
                "--fit-scope all_cells_exploratory explicitly."
            )
        train_mask = _train_mask_from_times(
            adata, args.time_key, [float(t) for t in args.train_times]
        )

    result = build_representation(
        adata,
        args.representation_id,
        train_mask=train_mask,
        n_components=args.n_components,
        n_hvg=args.n_hvg,
        seed=args.seed,
    )

    adata.obsm["X_rep"] = result.X_rep
    reps_meta = dict(adata.uns.get("representation_inputs_metadata", {}))
    reps_meta[args.representation_id] = result.metadata
    adata.uns["representation_inputs_metadata"] = reps_meta
    adata.uns["active_representation_id"] = args.representation_id

    out_path = Path(args.output_h5ad)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)

    sidecar = out_path.with_name(out_path.stem + ".representation_inputs.json")
    with open(sidecar, "w", encoding="utf-8") as f:
        json.dump(result.metadata, f, indent=2)

    print(
        f"[build_representation_inputs] {args.representation_id}: "
        f"X_rep {result.X_rep.shape} fit_scope={result.reducer_fit_scope} "
        f"-> {out_path}"
    )
    return out_path


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Build common 50-dim representation inputs (X_rep).",
    )
    p.add_argument("--input-h5ad", required=True)
    p.add_argument("--output-h5ad", required=True)
    p.add_argument("--representation-id", required=True,
                   choices=list(REPRESENTATION_ARMS))
    p.add_argument("--n-components", type=int, default=DEFAULT_N_COMPONENTS)
    p.add_argument("--n-hvg", type=int, default=DEFAULT_N_HVG)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--time-key", default="abs_day",
                   help="obs column with the timepoint for train/heldout split.")
    p.add_argument("--train-times", nargs="*", default=None,
                   help="Scenario training timepoints (required for train_only).")
    p.add_argument(
        "--fit-scope",
        default=TRAIN_ONLY,
        choices=[TRAIN_ONLY, ALL_CELLS_EXPLORATORY],
        help=(
            "train_only (formal, anti-leakage) fits scaler/PCA on training "
            "cells only; all_cells_exploratory fits on all cells and must NOT "
            "be used for formal rankings."
        ),
    )
    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        run(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
