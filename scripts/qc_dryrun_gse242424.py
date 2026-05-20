#!/usr/bin/env python3
"""
qc_dryrun_gse242424.py - Lightweight streaming QC dry-run for GSE242424

Streams each raw barcode-level MTX file and computes per-barcode QC metrics
without building an AnnData object or loading the full sparse matrix into
memory.  Per-barcode accumulators are O(n_barcodes) numpy arrays; coordinate
entries are read via pandas C-tokeniser in chunks so peak RAM stays at
O(n_barcodes + chunk_size), far below O(nnz).

Checks performed
  1. Gene dimension of each MTX == rows in GSE242423_scRNA_genes.tsv
  2. Barcode dimension of each MTX == lines in matching barcodes.tsv.gz

Outputs
  benchmark/reports/qc/gse242424_qc_dryrun_summary.csv
  Stdout: formatted summary table

Usage
  cd C:\\Users\\37620\\trajectory
  python scripts/qc_dryrun_gse242424.py

Dataset : GSE242424 / GSE242423 (iPSC reprogramming time-course)
"""

from __future__ import annotations

import gzip
import os
import subprocess
import sys
from pathlib import Path
import time as _time

import numpy as np
import pandas as pd

# pandas read_csv chunksize: 2M rows ~ 48 MB parsed int32 data
_PD_CHUNK = 2_000_000

# ---------------------------------------------------------------------------
# Project-root detection
# ---------------------------------------------------------------------------

def _find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
    # Walk upward from this script until we find a directory containing data/
    candidate = Path(__file__).resolve().parent
    for _ in range(6):
        candidate = candidate.parent
        if (candidate / "data").is_dir():
            return candidate
    raise FileNotFoundError(
        "Cannot locate project root (looked for a 'data/' sibling). "
        "Set TRAJ_PROJECT_ROOT env var or run from the project root."
    )

PROJECT_ROOT = _find_project_root()
os.chdir(PROJECT_ROOT)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR   = PROJECT_ROOT / "data" / "gse242424"
GENES_FILE = DATA_DIR / "GSE242423_scRNA_genes.tsv"
MTX_DIR    = DATA_DIR / "scRNA-seq"
OUT_CSV    = PROJECT_ROOT / "benchmark" / "reports" / "qc" / "gse242424_qc_dryrun_summary.csv"

# ---------------------------------------------------------------------------
# QC thresholds (mirrors existing project defaults from 01_load_and_qc.py)
# ---------------------------------------------------------------------------
MIN_COUNTS = 500
MIN_GENES  = 200
MAX_GENES  = 8000
MAX_PCT_MT = 20.0

# ---------------------------------------------------------------------------
# Sample manifest
# abs_day = -1 for iPSC: terminal / not yet formally assigned (see README)
# ---------------------------------------------------------------------------
SAMPLES: list[tuple[str, str, int]] = [
    ("GSM7763419_D0",   "D0",   0),
    ("GSM7763420_D2",   "D2",   2),
    ("GSM7763421_D4",   "D4",   4),
    ("GSM7763422_D6",   "D6",   6),
    ("GSM7763423_D8",   "D8",   8),
    ("GSM7763424_D10",  "D10",  10),
    ("GSM7763425_D12",  "D12",  12),
    ("GSM7763426_D14",  "D14",  14),
    ("GSM7763427_iPSC", "iPSC", -1),
]


# ===========================================================================
# 1. Gene loading
# ===========================================================================

def load_genes(genes_file: Path) -> tuple[list[str], list[str], np.ndarray]:
    """
    Read a Cell Ranger genes/features TSV (tab-separated, no header).
    Expected columns: ENSEMBL_id  gene_symbol  feature_type

    Returns (gene_ids, gene_symbols, mt_mask)
    where mt_mask is a bool array True for genes whose symbol starts with MT-.
    """
    gene_ids, gene_symbols = [], []
    with open(genes_file, "r") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            gene_ids.append(parts[0])
            gene_symbols.append(parts[1] if len(parts) > 1 else parts[0])
    mt_mask = np.array([s.startswith("MT-") for s in gene_symbols], dtype=bool)
    return gene_ids, gene_symbols, mt_mask


# ===========================================================================
# 2. Barcode counting
# ===========================================================================

def count_barcodes_gz(barcodes_gz: Path) -> int:
    """
    Return line count of a gzipped barcodes TSV.
    Tries zcat | wc -l (C-speed) first; falls back to Python gzip loop.
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


# ===========================================================================
# 3. MTX streaming QC
# ===========================================================================

def stream_mtx_qc(
    mtx_gz: Path,
    n_genes_ref: int,
    mt_mask: np.ndarray,
) -> tuple[int, int, int, np.ndarray, np.ndarray, np.ndarray]:
    """
    Stream a gzipped MTX file and compute per-barcode QC metrics.

    After reading the dimension header the file handle is passed directly to
    pandas.read_csv(chunksize=_PD_CHUNK) which uses the compiled C tokeniser
    (~20-50x faster than a Python line loop).  Each chunk is reduced to
    per-barcode bin-counts via np.bincount and accumulated into O(n_barcodes)
    arrays; the chunk is then discarded.

    Returns
    -------
    n_genes_mtx, n_barcodes_mtx, nnz_actual,
    total_counts, n_genes_by_counts, total_counts_mt
    """
    with gzip.open(mtx_gz, "rt") as fh:
        # skip comment / metadata lines
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()

        # first non-comment line is the dimension header
        parts = line.split()
        n_genes_mtx    = int(parts[0])
        n_barcodes_mtx = int(parts[1])

        if n_genes_mtx != n_genes_ref:
            raise ValueError(
                f"{mtx_gz.name}: gene dimension {n_genes_mtx} "
                f"!= reference {n_genes_ref}"
            )

        # O(n_barcodes) accumulators
        n = n_barcodes_mtx
        total_counts      = np.zeros(n, dtype=np.float64)
        n_genes_by_counts = np.zeros(n, dtype=np.int32)
        total_counts_mt   = np.zeros(n, dtype=np.float64)
        nnz_actual = 0

        # stream via pandas C tokeniser
        chunk_iter = pd.read_csv(
            fh,
            sep=r"\s+",
            header=None,
            names=["r", "c", "v"],
            dtype={"r": np.int32, "c": np.int32, "v": np.int32},
            engine="c",
            chunksize=_PD_CHUNK,
        )

        for chunk in chunk_iter:
            rows = chunk["r"].values - 1   # 0-based gene index
            cols = chunk["c"].values - 1   # 0-based barcode index
            vals = chunk["v"].values.astype(np.float64)

            total_counts      += np.bincount(cols, weights=vals, minlength=n)
            n_genes_by_counts += np.bincount(cols, minlength=n).astype(np.int32)

            mt_sel = mt_mask[rows]
            if mt_sel.any():
                total_counts_mt += np.bincount(
                    cols[mt_sel], weights=vals[mt_sel], minlength=n
                )
            nnz_actual += len(chunk)

    return (n_genes_mtx, n_barcodes_mtx, nnz_actual,
            total_counts, n_genes_by_counts, total_counts_mt)


# ===========================================================================
# 4. Per-sample summary
# ===========================================================================

def summarise_sample(
    sample_id: str,
    time_label: str,
    abs_day: int,
    mtx_gz: Path,
    barcodes_gz: Path,
    n_genes_ref: int,
    mt_mask: np.ndarray,
) -> dict:
    """Run streaming QC for one sample; return a summary dict."""
    t0 = _time.time()

    (n_genes_mtx, n_barcodes_mtx, nnz,
     total_counts, n_genes_by_counts, total_counts_mt) = stream_mtx_qc(
        mtx_gz, n_genes_ref, mt_mask
    )

    n_barcodes_file = count_barcodes_gz(barcodes_gz)
    if n_barcodes_file != n_barcodes_mtx:
        raise ValueError(
            f"{sample_id}: barcodes file has {n_barcodes_file} lines "
            f"but MTX col-dim is {n_barcodes_mtx}"
        )

    elapsed = _time.time() - t0

    # pct_counts_mt - safe against zero-count barcodes
    with np.errstate(invalid="ignore", divide="ignore"):
        pct_counts_mt = np.where(
            total_counts > 0,
            100.0 * total_counts_mt / np.maximum(total_counts, 1.0),
            0.0,
        )

    # QC filter
    pass_mask = (
        (total_counts      >= MIN_COUNTS)
        & (n_genes_by_counts >= MIN_GENES)
        & (n_genes_by_counts <= MAX_GENES)
        & (pct_counts_mt    <= MAX_PCT_MT)
    )
    n_pass = int(pass_mask.sum())

    def _med(arr: np.ndarray, mask=None) -> float:
        x = arr[mask] if mask is not None else arr
        return float(np.median(x)) if len(x) > 0 else float("nan")

    row = {
        "sample_id":                sample_id,
        "time_label":               time_label,
        "abs_day":                  abs_day,
        "matrix_rows":              n_genes_mtx,
        "matrix_cols":              n_barcodes_mtx,
        "nnz":                      nnz,
        "n_barcodes":               n_barcodes_file,
        "n_pass_qc":                n_pass,
        "pass_qc_fraction":         round(n_pass / n_barcodes_mtx, 6) if n_barcodes_mtx else float("nan"),
        "median_total_counts_all":  round(_med(total_counts), 2),
        "median_total_counts_pass": round(_med(total_counts, pass_mask), 2),
        "median_n_genes_all":       round(_med(n_genes_by_counts.astype(float)), 2),
        "median_n_genes_pass":      round(_med(n_genes_by_counts.astype(float), pass_mask), 2),
        "median_pct_mt_all":        round(_med(pct_counts_mt), 4),
        "median_pct_mt_pass":       round(_med(pct_counts_mt, pass_mask), 4),
        "min_total_counts_pass":    float(total_counts[pass_mask].min()) if n_pass else float("nan"),
        "max_total_counts_pass":    float(total_counts[pass_mask].max()) if n_pass else float("nan"),
        # internal diagnostics (popped before CSV)
        "_elapsed_s":               round(elapsed, 1),
        "_gene_dim_ok":             (n_genes_mtx == n_genes_ref),
        "_barcode_dim_ok":          (n_barcodes_file == n_barcodes_mtx),
    }
    return row


# ===========================================================================
# 5. Main
# ===========================================================================

def main() -> None:
    T0 = _time.time()

    if not GENES_FILE.exists():
        sys.exit(f"ERROR: genes file not found: {GENES_FILE}")
    if not MTX_DIR.is_dir():
        sys.exit(f"ERROR: MTX directory not found: {MTX_DIR}")

    print("=" * 70)
    print("GSE242424 -- QC DRY-RUN")
    print("=" * 70)
    gene_ids, gene_symbols, mt_mask = load_genes(GENES_FILE)
    n_genes_ref = len(gene_ids)
    mt_names = [s for s in gene_symbols if s.startswith("MT-")]
    print(f"  Genes file  : {GENES_FILE.name}  ({n_genes_ref:,} genes)")
    print(f"  MT- genes   : {len(mt_names)} ({', '.join(mt_names[:6])} ...)")
    print(f"  QC thresholds: counts>={MIN_COUNTS} | genes in [{MIN_GENES},{MAX_GENES}] | pct_mt<={MAX_PCT_MT}")
    print()

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for sample_id, time_label, abs_day in SAMPLES:
        mtx_gz      = MTX_DIR / f"{sample_id}.matrix.mtx.gz"
        barcodes_gz = MTX_DIR / f"{sample_id}.barcodes.tsv.gz"

        if not mtx_gz.exists():
            print(f"  [SKIP] {sample_id}: MTX not found")
            continue
        if not barcodes_gz.exists():
            print(f"  [SKIP] {sample_id}: barcodes file not found")
            continue

        size_mb = mtx_gz.stat().st_size / 1024 / 1024
        print(f"  Processing {sample_id:25s} ({size_mb:6.1f} MB) ... ", end="", flush=True)

        try:
            rec = summarise_sample(
                sample_id, time_label, abs_day,
                mtx_gz, barcodes_gz, n_genes_ref, mt_mask,
            )
        except Exception as exc:
            print(f"ERROR: {exc}")
            continue

        dim_ok  = rec.pop("_gene_dim_ok") and rec.pop("_barcode_dim_ok")
        elapsed = rec.pop("_elapsed_s")
        rows.append(rec)

        print(
            f"done in {elapsed:5.1f}s | dims OK={dim_ok} | "
            f"n_barcodes={rec['n_barcodes']:>9,} | "
            f"n_pass={rec['n_pass_qc']:>7,} ({100*rec['pass_qc_fraction']:.2f}%)"
        )

    if not rows:
        sys.exit("ERROR: no samples processed")

    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    print()
    print(f"  CSV written -> {OUT_CSV}")
    print()

    display_cols = [
        "sample_id", "time_label", "abs_day",
        "n_barcodes", "n_pass_qc", "pass_qc_fraction",
        "median_total_counts_pass", "median_n_genes_pass", "median_pct_mt_pass",
    ]
    print("=" * 70)
    print("SUMMARY TABLE")
    print("=" * 70)
    print(df[display_cols].to_string(index=False))
    print()
    print(f"Total wall time: {_time.time() - T0:.1f}s")


if __name__ == "__main__":
    main()


# ===========================================================================
# Entry-point with optional batch slicing
# ===========================================================================
# Called by batch runner below; main() above handles full run.
