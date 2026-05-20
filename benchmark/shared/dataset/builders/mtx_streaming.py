#!/usr/bin/env python3
"""
mtx_streaming.py -- Reusable streaming reader for 10x-style gzipped MTX files.

Public API
----------
stream_mtx_and_filter(mtx_gz, barcodes_gz, n_genes, mt_mask, ...)
    -> StreamResult

Design
------
Two streaming passes over the gzip-compressed MatrixMarket file, each driven
by the pandas C tokeniser (pd.read_csv with chunksize).  Peak RAM is
O(n_barcodes + chunk_size + filtered_nnz), NOT O(raw nnz).

  Pass 0  -- accumulate per-barcode QC metrics into O(n_barcodes) arrays
  Pass 1  -- emit COO triples for QC-passing barcodes only

The returned X_csr has shape (n_pass, n_genes) with raw integer counts.
Callers are responsible for normalization and log-transform.

Typical usage
-------------
    from benchmark.shared.dataset.builders.mtx_streaming import (
        StreamResult, stream_mtx_and_filter,
    )

    result = stream_mtx_and_filter(
        mtx_gz=Path("sample.matrix.mtx.gz"),
        barcodes_gz=Path("sample.barcodes.tsv.gz"),
        n_genes=36601,
        mt_mask=mt_mask_array,
        ribo_mask=ribo_mask_array,
        min_counts=500, min_genes=200, max_genes=8000, max_pct_mt=20.0,
    )
    # result.X_csr  -- scipy.sparse.csr_matrix (n_pass, n_genes)
    # result.barcodes -- np.ndarray of barcode strings (n_pass,)
    # result.qc_df    -- pd.DataFrame with QC columns for passing cells
"""

from __future__ import annotations

import gzip
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp


# Coordinate entries per pandas read_csv chunk.  2 M rows ~ 48 MB for int32
# triples, comfortably within typical per-sample RAM budgets.
_DEFAULT_CHUNK = 2_000_000


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class StreamResult:
    """
    All outputs from a single stream_mtx_and_filter() call.

    Attributes
    ----------
    X_csr          : scipy CSR matrix, shape (n_pass, n_genes), dtype float32
                     Contains RAW integer counts (NOT normalized).
    barcodes       : str array, shape (n_pass,) -- QC-passing barcode strings
    qc_df          : DataFrame, index=barcodes, columns:
                       n_genes_by_counts, total_counts,
                       total_counts_mt, pct_counts_mt,
                       total_counts_ribo, pct_counts_ribo
    n_genes        : validated gene dimension from the MTX header
    n_barcodes_raw : raw barcode count (MTX column dimension)
    nnz_raw        : actual non-zero entries streamed from the raw MTX
    n_pass         : number of QC-passing barcodes
    pass_fraction  : n_pass / n_barcodes_raw
    """
    X_csr: sp.csr_matrix
    barcodes: np.ndarray
    qc_df: pd.DataFrame
    n_genes: int
    n_barcodes_raw: int
    nnz_raw: int
    n_pass: int
    pass_fraction: float


# ---------------------------------------------------------------------------
# MTX header utilities
# ---------------------------------------------------------------------------

def _read_mtx_header(fh) -> tuple[int, int, int]:
    """
    Consume comment/metadata lines and the dimension line from an open MTX
    text file handle.  Returns (n_genes, n_barcodes, nnz_header).
    """
    line = fh.readline()
    while line.startswith("%"):
        line = fh.readline()
    parts = line.split()
    return int(parts[0]), int(parts[1]), int(parts[2])


def read_mtx_dims(mtx_gz: Path) -> tuple[int, int, int]:
    """
    Open a gzipped MTX, parse the header, and return (n_genes, n_barcodes, nnz).
    File is closed immediately after; no data lines are read.
    """
    with gzip.open(mtx_gz, "rt") as fh:
        return _read_mtx_header(fh)


# ---------------------------------------------------------------------------
# Barcode file utilities
# ---------------------------------------------------------------------------

def count_barcodes_gz(barcodes_gz: Path) -> int:
    """
    Return line count of a gzipped barcodes TSV.
    Fast path: zcat | wc -l (C implementation).
    Fallback: pure-Python gzip loop.
    """
    try:
        result = subprocess.run(
            ["sh", "-c", f"zcat '{barcodes_gz}' | wc -l"],
            capture_output=True, text=True, check=True,
        )
        return int(result.stdout.strip())
    except Exception:
        n = 0
        with gzip.open(barcodes_gz, "rt") as fh:
            for _ in fh:
                n += 1
        return n


def read_barcodes_gz(barcodes_gz: Path) -> np.ndarray:
    """
    Read all barcodes from a gzipped single-column TSV into a string array.
    Raises FileNotFoundError if the file is absent.
    """
    if not barcodes_gz.exists():
        raise FileNotFoundError(f"Barcodes file not found: {barcodes_gz}")
    df = pd.read_csv(barcodes_gz, sep="\t", header=None, names=["barcode"], dtype=str)
    return df["barcode"].to_numpy()


# ---------------------------------------------------------------------------
# Pass 0 -- per-barcode QC accumulation
# ---------------------------------------------------------------------------

def _pass0_qc(
    mtx_gz: Path,
    n_barcodes: int,
    mt_mask: np.ndarray,
    ribo_mask: Optional[np.ndarray],
    chunk_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Stream MTX and accumulate per-barcode QC arrays without storing COO data.

    Returns
    -------
    total_counts, n_genes_by_counts, total_counts_mt, total_counts_ribo,
    nnz_raw
    All count arrays have shape (n_barcodes,).
    """
    n = n_barcodes
    total_counts      = np.zeros(n, dtype=np.float64)
    n_genes_by_counts = np.zeros(n, dtype=np.int32)
    total_counts_mt   = np.zeros(n, dtype=np.float64)
    total_counts_ribo = np.zeros(n, dtype=np.float64)
    nnz_raw = 0

    with gzip.open(mtx_gz, "rt") as fh:
        _read_mtx_header(fh)   # consume header; validation done by caller
        for chunk in pd.read_csv(
            fh, sep=r"\s+", header=None,
            names=["r", "c", "v"],
            dtype={"r": np.int32, "c": np.int32, "v": np.int32},
            engine="c", chunksize=chunk_size,
        ):
            gene_idx = chunk["r"].values - 1   # 0-based gene index
            bc_idx   = chunk["c"].values - 1   # 0-based barcode index
            vals     = chunk["v"].values.astype(np.float64)

            total_counts      += np.bincount(bc_idx, weights=vals, minlength=n)
            n_genes_by_counts += np.bincount(bc_idx, minlength=n).astype(np.int32)

            mt_sel = mt_mask[gene_idx]
            if mt_sel.any():
                total_counts_mt += np.bincount(
                    bc_idx[mt_sel], weights=vals[mt_sel], minlength=n
                )

            if ribo_mask is not None:
                ribo_sel = ribo_mask[gene_idx]
                if ribo_sel.any():
                    total_counts_ribo += np.bincount(
                        bc_idx[ribo_sel], weights=vals[ribo_sel], minlength=n
                    )

            nnz_raw += len(chunk)

    return total_counts, n_genes_by_counts, total_counts_mt, total_counts_ribo, nnz_raw


# ---------------------------------------------------------------------------
# QC filter application
# ---------------------------------------------------------------------------

def _apply_qc_filter(
    total_counts: np.ndarray,
    n_genes_by_counts: np.ndarray,
    total_counts_mt: np.ndarray,
    total_counts_ribo: np.ndarray,
    barcodes_all: np.ndarray,
    min_counts: float,
    min_genes: int,
    max_genes: int,
    max_pct_mt: float,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """
    Apply QC thresholds and build the passing-cell QC DataFrame.

    Returns
    -------
    kept_idx  : int array of raw-barcode indices that pass QC
    kept_rank : int32 array, length n_barcodes_raw, maps raw idx -> kept rank
                (-1 for dropped barcodes)
    qc_df     : DataFrame for passing cells only, indexed by barcode string
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        pct_counts_mt = np.where(
            total_counts > 0,
            100.0 * total_counts_mt / np.maximum(total_counts, 1.0),
            0.0,
        )
        pct_counts_ribo = np.where(
            total_counts > 0,
            100.0 * total_counts_ribo / np.maximum(total_counts, 1.0),
            0.0,
        )

    kept_mask = (
        (total_counts      >= min_counts)
        & (n_genes_by_counts >= min_genes)
        & (n_genes_by_counts <= max_genes)
        & (pct_counts_mt    <= max_pct_mt)
    )
    kept_idx  = np.where(kept_mask)[0]
    kept_rank = np.full(len(total_counts), -1, dtype=np.int32)
    kept_rank[kept_idx] = np.arange(len(kept_idx), dtype=np.int32)

    qc_df = pd.DataFrame(
        {
            "n_genes_by_counts": n_genes_by_counts[kept_idx].astype(np.int32),
            "total_counts":      total_counts[kept_idx].astype(np.float32),
            "total_counts_mt":   total_counts_mt[kept_idx].astype(np.float32),
            "pct_counts_mt":     pct_counts_mt[kept_idx].astype(np.float32),
            "total_counts_ribo": total_counts_ribo[kept_idx].astype(np.float32),
            "pct_counts_ribo":   pct_counts_ribo[kept_idx].astype(np.float32),
        },
        index=barcodes_all[kept_idx],
    )
    qc_df.index.name = "barcode"

    return kept_idx, kept_rank, qc_df


# ---------------------------------------------------------------------------
# Pass 1 -- filtered COO collection
# ---------------------------------------------------------------------------

def _pass1_filtered_coo(
    mtx_gz: Path,
    kept_rank: np.ndarray,
    chunk_size: int,
) -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
    """
    Second streaming pass: collect COO triples for QC-passing barcodes only.

    Barcode indices are remapped to 0..n_pass-1 using kept_rank.
    Gene indices are kept as-is (0-based).

    Returns
    -------
    gene_idx_list  : list of int32 arrays -- 0-based gene indices
    cell_rank_list : list of int32 arrays -- remapped 0-based cell ranks
    data_list      : list of float32 arrays -- raw counts
    """
    gene_idx_list:  list[np.ndarray] = []
    cell_rank_list: list[np.ndarray] = []
    data_list:      list[np.ndarray] = []

    with gzip.open(mtx_gz, "rt") as fh:
        _read_mtx_header(fh)
        for chunk in pd.read_csv(
            fh, sep=r"\s+", header=None,
            names=["r", "c", "v"],
            dtype={"r": np.int32, "c": np.int32, "v": np.int32},
            engine="c", chunksize=chunk_size,
        ):
            gene_idx = chunk["r"].values - 1   # 0-based gene index
            bc_idx   = chunk["c"].values - 1   # 0-based barcode index
            vals     = chunk["v"].values

            cell_rank = kept_rank[bc_idx]      # -1 for dropped barcodes
            valid     = cell_rank >= 0

            if valid.any():
                gene_idx_list.append(gene_idx[valid].astype(np.int32))
                cell_rank_list.append(cell_rank[valid].astype(np.int32))
                data_list.append(vals[valid].astype(np.float32))

    return gene_idx_list, cell_rank_list, data_list


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def stream_mtx_and_filter(
    mtx_gz: Path,
    barcodes_gz: Path,
    n_genes: int,
    mt_mask: np.ndarray,
    ribo_mask: Optional[np.ndarray] = None,
    min_counts: float = 500.0,
    min_genes: int = 200,
    max_genes: int = 8000,
    max_pct_mt: float = 20.0,
    chunk_size: int = _DEFAULT_CHUNK,
) -> StreamResult:
    """
    Stream a raw barcode-level MTX file and return a QC-filtered sparse matrix.

    Performs two passes over the gzipped MTX file:
      Pass 0  per-barcode QC accumulation (O(n_barcodes) RAM)
      Pass 1  filtered COO collection     (O(filtered_nnz) RAM)

    The returned X_csr has shape (n_pass, n_genes) with raw integer counts
    (dtype float32, ready for total-count normalization and log1p by caller).

    Parameters
    ----------
    mtx_gz      : path to *.matrix.mtx.gz
    barcodes_gz : path to *.barcodes.tsv.gz
    n_genes     : expected gene count; validated against MTX header
    mt_mask     : bool array length n_genes, True for MT- genes
    ribo_mask   : bool array length n_genes, True for RPS*/RPL* genes (optional)
    min_counts  : minimum total UMI count per cell
    min_genes   : minimum number of genes detected per cell
    max_genes   : maximum number of genes detected per cell
    max_pct_mt  : maximum mitochondrial percentage per cell
    chunk_size  : coordinate entries per pandas read_csv chunk

    Returns
    -------
    StreamResult -- see dataclass docstring for field descriptions.

    Raises
    ------
    FileNotFoundError : if either input file is absent
    ValueError        : if MTX gene dimension != n_genes, or barcode file line
                        count != MTX column dimension
    """
    if not mtx_gz.exists():
        raise FileNotFoundError(f"MTX file not found: {mtx_gz}")
    if not barcodes_gz.exists():
        raise FileNotFoundError(f"Barcodes file not found: {barcodes_gz}")

    # --- Validate MTX header ---
    n_genes_mtx, n_barcodes_mtx, _nnz_hdr = read_mtx_dims(mtx_gz)
    if n_genes_mtx != n_genes:
        raise ValueError(
            f"{mtx_gz.name}: gene dimension {n_genes_mtx} != expected {n_genes}"
        )

    # --- Load and validate barcodes ---
    barcodes_all = read_barcodes_gz(barcodes_gz)
    if len(barcodes_all) != n_barcodes_mtx:
        raise ValueError(
            f"{barcodes_gz.name}: {len(barcodes_all)} lines "
            f"!= MTX column dimension {n_barcodes_mtx}"
        )

    # --- Pass 0: per-barcode QC ---
    tc, ngc, tc_mt, tc_ribo, nnz_raw = _pass0_qc(
        mtx_gz, n_barcodes_mtx, mt_mask, ribo_mask, chunk_size
    )

    # --- Apply QC thresholds ---
    kept_idx, kept_rank, qc_df = _apply_qc_filter(
        tc, ngc, tc_mt, tc_ribo, barcodes_all,
        min_counts, min_genes, max_genes, max_pct_mt,
    )
    n_pass = len(kept_idx)

    if n_pass == 0:
        return StreamResult(
            X_csr=sp.csr_matrix((0, n_genes), dtype=np.float32),
            barcodes=barcodes_all[kept_idx],
            qc_df=qc_df,
            n_genes=n_genes_mtx,
            n_barcodes_raw=n_barcodes_mtx,
            nnz_raw=nnz_raw,
            n_pass=0,
            pass_fraction=0.0,
        )

    # --- Pass 1: collect filtered COO ---
    gene_idx_list, cell_rank_list, data_list = _pass1_filtered_coo(
        mtx_gz, kept_rank, chunk_size
    )

    # --- Build CSR (rows=cells, cols=genes) ---
    if gene_idx_list:
        coo_row = np.concatenate(cell_rank_list)  # cell axis (0..n_pass-1)
        coo_col = np.concatenate(gene_idx_list)   # gene axis (0..n_genes-1)
        coo_dat = np.concatenate(data_list)
    else:
        coo_row = np.empty(0, dtype=np.int32)
        coo_col = np.empty(0, dtype=np.int32)
        coo_dat = np.empty(0, dtype=np.float32)

    X_csr = sp.coo_matrix(
        (coo_dat, (coo_row, coo_col)),
        shape=(n_pass, n_genes),
        dtype=np.float32,
    ).tocsr()

    return StreamResult(
        X_csr=X_csr,
        barcodes=barcodes_all[kept_idx],
        qc_df=qc_df,
        n_genes=n_genes_mtx,
        n_barcodes_raw=n_barcodes_mtx,
        nnz_raw=nnz_raw,
        n_pass=n_pass,
        pass_fraction=n_pass / n_barcodes_mtx,
    )
