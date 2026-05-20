#!/usr/bin/env python3
"""
Build the GSE242424 post-QC raw/full-gene h5ad for benchmark integration.

Dataset : GSE242424 / GSE242423 (chemical iPSC reprogramming time-course,
          human fibroblast -> iPSC, D0-D14 + iPSC terminal state)
Output  : data/processed/gse242424_human/
              GSE242424_raw_full_gene_benchmark_input.h5ad

Key differences from GSE230659 / GSE178325 builders
-----------------------------------------------------
1. File layout : flat files in data/gse242424/scRNA-seq/ named
     {GSM}_{label}.matrix.mtx.gz / {GSM}_{label}.barcodes.tsv.gz
   (no per-sample subdirectories; single shared genes TSV)
2. Raw barcode matrices : ~2.4-2.7 M raw barcodes per sample -- streamed
   and QC-filtered via benchmark.shared.dataset.builders.mtx_streaming.
3. iPSC abs_day = 16 : PROVISIONAL terminal time.  Review before any
   formal benchmark run that uses abs_day for trajectory inference.

Usage
-----
    conda activate traj_env
    cd C:\\Users\\37620\\trajectory

    # Full build (all 9 samples, ~4-5 min):
    python scripts/build_gse242424_raw_full_gene_input.py

    # Smoke test (single sample):
    python scripts/build_gse242424_raw_full_gene_input.py --samples GSM7763420_D2

Primary output
    data/processed/gse242424_human/GSE242424_raw_full_gene_benchmark_input.h5ad

Secondary outputs
    benchmark/reports/qc/{TIMESTAMP}_gse242424_qc_summary.csv

Revision history
    v0  2026-05-14  Skeleton -- paths, manifest, helpers, TODO stubs.
    v1  2026-05-14  Implemented via mtx_streaming reusable module.
                    Single streaming call per sample; gene stats from CSR.
                    Added --samples CLI argument for smoke testing.
"""

from __future__ import annotations

import argparse
import gc
import os
import sys
import time as _time
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    import psutil
except ImportError:
    psutil = None


TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")
T0 = _time.time()

DATASET_ID = "GSE242424"


# ===========================================================================
# Project-root detection
# ===========================================================================

def _find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
        raise FileNotFoundError(f"TRAJ_PROJECT_ROOT={env!r} does not exist.")
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "data").exists():
            return candidate
    return here.parent


PROJECT_ROOT = _find_project_root()

# Load mtx_streaming directly from its file path to avoid triggering
# benchmark/shared/dataset/__init__.py's pkgutil auto-import (which would
# require scanpy and other heavy deps not needed here).
import importlib.util as _ilu

_MTX_MOD_PATH = (
    PROJECT_ROOT / "benchmark" / "shared" / "dataset" / "builders" / "mtx_streaming.py"
)
_spec = _ilu.spec_from_file_location("mtx_streaming", _MTX_MOD_PATH)
_mtx_mod = _ilu.module_from_spec(_spec)
sys.modules["mtx_streaming"] = _mtx_mod  # required so @dataclass can resolve __module__
_spec.loader.exec_module(_mtx_mod)
StreamResult = _mtx_mod.StreamResult
stream_mtx_and_filter = _mtx_mod.stream_mtx_and_filter


# ===========================================================================
# Paths
# ===========================================================================

DATA_DIR      = PROJECT_ROOT / "data" / "gse242424"
MTX_DIR       = DATA_DIR / "scRNA-seq"
GENES_FILE    = DATA_DIR / "GSE242423_scRNA_genes.tsv"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed" / "gse242424_human"
OUTPUT_H5AD   = PROCESSED_DIR / "GSE242424_raw_full_gene_benchmark_input.h5ad"
REPORT_DIR    = PROJECT_ROOT / "benchmark" / "reports" / "qc"


# ===========================================================================
# QC constants (project-wide defaults)
# ===========================================================================

QC_PARAMS: dict = {
    "min_genes_per_cell":  200,
    "max_genes_per_cell":  8000,
    "min_cells_per_gene":  3,
    "max_pct_mito":        20.0,
    "min_counts_per_cell": 500,
}
NORM_TARGET_SUM: float = 1e4


# ===========================================================================
# Sample manifest
#
# (sample_id, time_label, abs_day)
#
# NOTE: iPSC abs_day=16 is PROVISIONAL.
# The terminal iPSC state has no precise experimental day in the GEO
# submission.  16 is a placeholder (14-day protocol + small offset).
# Revise before any formal benchmark run using abs_day for trajectory models.
# ===========================================================================

SAMPLES: list[tuple[str, str, float]] = [
    ("GSM7763419_D0",   "D0",    0.0),
    ("GSM7763420_D2",   "D2",    2.0),
    ("GSM7763421_D4",   "D4",    4.0),
    ("GSM7763422_D6",   "D6",    6.0),
    ("GSM7763423_D8",   "D8",    8.0),
    ("GSM7763424_D10",  "D10",  10.0),
    ("GSM7763425_D12",  "D12",  12.0),
    ("GSM7763426_D14",  "D14",  14.0),
    # PROVISIONAL terminal time -- abs_day=16 must be reviewed before formal runs
    ("GSM7763427_iPSC", "iPSC", 16.0),
]

TIMELABEL_TO_CELLTYPE_PROXY: dict[str, str] = {
    "D0":   "somatic",
    "D2":   "early_reprogramming",
    "D4":   "early_reprogramming",
    "D6":   "intermediate",
    "D8":   "intermediate",
    "D10":  "intermediate",
    "D12":  "late_reprogramming",
    "D14":  "late_reprogramming",
    "iPSC": "pluripotent",
}


# ===========================================================================
# Utilities
# ===========================================================================

def rss_gb() -> float:
    if psutil is None:
        return float("nan")
    return psutil.Process().memory_info().rss / 1e9


def mtx_path(sample_id: str) -> Path:
    return MTX_DIR / f"{sample_id}.matrix.mtx.gz"


def barcodes_path(sample_id: str) -> Path:
    return MTX_DIR / f"{sample_id}.barcodes.tsv.gz"


def check_input_files(samples: list[tuple[str, str, float]]) -> None:
    missing: list[str] = []
    if not GENES_FILE.exists():
        missing.append(str(GENES_FILE))
    for sample_id, _, _ in samples:
        for p in (mtx_path(sample_id), barcodes_path(sample_id)):
            if not p.exists():
                missing.append(str(p))
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} required file(s) not found:\n"
            + "\n".join(f"  {p}" for p in missing)
        )


def load_shared_genes(
    genes_file: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Read the shared genes TSV (no header, columns: ENSEMBL_id gene_symbol feature_type).
    Returns (gene_ids, gene_names, feature_types, mt_mask, ribo_mask).
    """
    if not genes_file.exists():
        raise FileNotFoundError(f"Genes file not found: {genes_file}")
    df = pd.read_csv(
        genes_file, sep="\t", header=None,
        names=["gene_id", "gene_name", "feature_type"], dtype=str,
    )
    gene_ids      = df["gene_id"].to_numpy(dtype=object)
    gene_names    = df["gene_name"].to_numpy(dtype=object)
    feature_types = df["feature_type"].to_numpy(dtype=object)
    mt_mask   = np.array([str(g).startswith("MT-") for g in gene_names], dtype=bool)
    ribo_mask = np.array(
        [str(g).startswith("RPS") or str(g).startswith("RPL") for g in gene_names],
        dtype=bool,
    )
    return gene_ids, gene_names, feature_types, mt_mask, ribo_mask


# ===========================================================================
# obs / var builders
# ===========================================================================

def build_obs(
    sample_id: str,
    time_label: str,
    abs_day: float,
    result: StreamResult,
) -> pd.DataFrame:
    """
    Build the per-cell obs DataFrame for one sample from a StreamResult.
    All cells in result are already QC-passing.
    """
    cell_type_proxy = TIMELABEL_TO_CELLTYPE_PROXY.get(time_label, "unknown")
    index = [f"{sample_id}_{bc}" for bc in result.barcodes]
    return pd.DataFrame(
        {
            "cell_id":               index,
            "barcode":               result.barcodes,
            "sample_id":             sample_id,
            "time_label":            time_label,
            "abs_day":               abs_day,
            "stage":                 "iPSC" if time_label == "iPSC" else "reprogramming",
            "stage_day_label":       time_label,
            "scTimeBench_timepoint": abs_day,
            "scTimeBench_cell_type": cell_type_proxy,
            "dataset_id":            DATASET_ID,
            "n_genes_by_counts":     result.qc_df["n_genes_by_counts"].to_numpy(),
            "total_counts":          result.qc_df["total_counts"].to_numpy(),
            "total_counts_mt":       result.qc_df["total_counts_mt"].to_numpy(),
            "pct_counts_mt":         result.qc_df["pct_counts_mt"].to_numpy(),
            "total_counts_ribo":     result.qc_df["total_counts_ribo"].to_numpy(),
            "pct_counts_ribo":       result.qc_df["pct_counts_ribo"].to_numpy(),
        },
        index=index,
    )


def build_var(
    gene_ids: np.ndarray,
    gene_names: np.ndarray,
    feature_types: np.ndarray,
    mt_mask: np.ndarray,
    ribo_mask: np.ndarray,
    gene_passes: np.ndarray,
    gene_n_cells: np.ndarray,
    gene_sum: np.ndarray,
    total_cells: int,
) -> pd.DataFrame:
    """Build the per-gene var DataFrame after the min-cells-per-gene filter."""
    idx = np.where(gene_passes)[0]
    var = pd.DataFrame(
        {
            "gene_id":               gene_ids[idx].astype(str),
            "gene_symbol":           gene_names[idx].astype(str),
            "feature_type":          feature_types[idx].astype(str),
            "gene_ids":              gene_ids[idx].astype(str),
            "gene_name":             gene_names[idx].astype(str),
            "feature_types":         feature_types[idx].astype(str),
            "mt":                    mt_mask[idx],
            "ribo":                  ribo_mask[idx],
            "n_cells":               gene_n_cells[idx].astype(np.int32),
            "n_cells_by_counts":     gene_n_cells[idx].astype(np.int32),
            "mean_counts":           (gene_sum[idx] / max(total_cells, 1)).astype(np.float32),
            "pct_dropout_by_counts": (
                100.0 * (1.0 - gene_n_cells[idx] / max(total_cells, 1))
            ).astype(np.float32),
            "total_counts":          gene_sum[idx].astype(np.float32),
            "passes_filter":         True,
            "highly_variable":       False,
            "is_hvg":                False,
        },
        index=gene_names[idx].astype(str),
    )
    # Deduplicate gene-symbol index (6 genes in the TSV appear twice).
    # Append -1/-2 suffixes so AnnData index is unique without warning.
    seen: dict[str, int] = {}
    deduped = []
    for name in var.index:
        if name in seen:
            seen[name] += 1
            deduped.append(f"{name}-{seen[name]}")
        else:
            seen[name] = 0
            deduped.append(name)
    var.index = deduped
    var.index.name = "var_names"
    return var


# ===========================================================================
# Thin wrapper around the reusable streaming module
# ===========================================================================

def _stream_sample(
    sample_id: str,
    n_genes_ref: int,
    mt_mask: np.ndarray,
    ribo_mask: np.ndarray,
) -> StreamResult:
    """
    Call stream_mtx_and_filter() for one GSE242424 sample using the
    project-wide QC thresholds.  Paths are resolved from the module-level
    MTX_DIR constant.
    """
    return stream_mtx_and_filter(
        mtx_gz=mtx_path(sample_id),
        barcodes_gz=barcodes_path(sample_id),
        n_genes=n_genes_ref,
        mt_mask=mt_mask,
        ribo_mask=ribo_mask,
        min_counts=float(QC_PARAMS["min_counts_per_cell"]),
        min_genes=int(QC_PARAMS["min_genes_per_cell"]),
        max_genes=int(QC_PARAMS["max_genes_per_cell"]),
        max_pct_mt=float(QC_PARAMS["max_pct_mito"]),
    )


# ===========================================================================
# Main build pipeline
# ===========================================================================

def main(samples_subset: list[str] | None = None) -> None:
    """
    Single-pass build pipeline for GSE242424.

    For each sample:
      - Call _stream_sample() which internally does Pass 0 (QC) + Pass 1
        (filtered COO collection) and returns a StreamResult with X_csr.
      - Accumulate per-gene cell counts and sums from X_csr (fast CSR ops).
      - Cache the StreamResult for the normalization step.

    After all samples are streamed:
      - Apply min_cells_per_gene filter.
      - For each cached X_csr, apply the gene filter + total-count
        normalization + log1p, writing COO triples into accumulator lists.
      - Assemble final CSR, construct AnnData, write h5ad.
      - Write QC summary CSV.
    """
    print(f"[{TIMESTAMP}] {DATASET_ID} raw/full-gene builder  v1")
    print(f"  Project root : {PROJECT_ROOT}")
    print(f"  Output h5ad  : {OUTPUT_H5AD}")
    if samples_subset:
        print(f"  Sample limit : {samples_subset}")
    print()

    # Select samples
    samples_to_run = [
        (sid, tl, ad_day) for sid, tl, ad_day in SAMPLES
        if (samples_subset is None or sid in samples_subset)
    ]
    if not samples_to_run:
        sys.exit(f"ERROR: none of {samples_subset} matched SAMPLES manifest.")

    check_input_files(samples_to_run)
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    # Load gene annotations
    gene_ids, gene_names, feature_types, mt_mask, ribo_mask = load_shared_genes(
        GENES_FILE
    )
    N_GENES = len(gene_ids)
    print(f"  Genes : {N_GENES:,}  MT-genes : {mt_mask.sum()}  "
          f"Ribo-genes : {ribo_mask.sum()}")
    print()

    # Per-gene accumulators (int64 / float64 for precision)
    gene_n_cells = np.zeros(N_GENES, dtype=np.int64)
    gene_sum     = np.zeros(N_GENES, dtype=np.float64)

    obs_parts:      list[pd.DataFrame]  = []
    cached_results: list[StreamResult]  = []
    qc_rows:        list[dict]          = []
    total_cells_passed = 0

    # ------------------------------------------------------------------
    # Stream + filter all samples; cache StreamResults for normalization
    # ------------------------------------------------------------------
    print("Streaming and QC-filtering samples ...")
    for sample_id, time_label, abs_day in samples_to_run:
        t0s = _time.time()
        print(f"  {sample_id:25s}  ... ", end="", flush=True)

        result = _stream_sample(sample_id, N_GENES, mt_mask, ribo_mask)

        # Accumulate per-gene stats from the filtered CSR (fast sparse ops)
        gene_n_cells += np.asarray(result.X_csr.getnnz(axis=0)).ravel().astype(np.int64)
        gene_sum     += np.asarray(result.X_csr.sum(axis=0)).ravel()

        obs_parts.append(build_obs(sample_id, time_label, abs_day, result))
        cached_results.append(result)
        total_cells_passed += result.n_pass

        qc_rows.append({
            "sample_id":                 sample_id,
            "time_label":                time_label,
            "abs_day":                   abs_day,
            "n_barcodes_raw":            result.n_barcodes_raw,
            "nnz_raw":                   result.nnz_raw,
            "n_pass_qc":                 result.n_pass,
            "pass_qc_fraction":          round(result.pass_fraction, 6),
            "median_total_counts_pass":  round(float(np.median(
                result.qc_df["total_counts"])), 2) if result.n_pass else float("nan"),
            "median_n_genes_pass":       round(float(np.median(
                result.qc_df["n_genes_by_counts"])), 2) if result.n_pass else float("nan"),
            "median_pct_mt_pass":        round(float(np.median(
                result.qc_df["pct_counts_mt"])), 4) if result.n_pass else float("nan"),
        })

        elapsed = _time.time() - t0s
        print(
            f"done {elapsed:5.1f}s | raw={result.n_barcodes_raw:>9,} | "
            f"pass={result.n_pass:>7,} ({100*result.pass_fraction:.2f}%) | "
            f"RSS={rss_gb():.1f}GB"
        )

    obs_adata = pd.concat(obs_parts)
    obs_adata.index.name = "obs_names"
    del obs_parts
    print(f"\n  Total pass-QC cells : {total_cells_passed:,}")

    # ------------------------------------------------------------------
    # Gene filter
    # ------------------------------------------------------------------
    gene_passes   = gene_n_cells >= QC_PARAMS["min_cells_per_gene"]
    gene_indices  = np.where(gene_passes)[0]
    gene_rank     = np.full(N_GENES, -1, dtype=np.int32)
    gene_rank[gene_indices] = np.arange(len(gene_indices), dtype=np.int32)
    n_genes_out   = len(gene_indices)

    var = build_var(
        gene_ids, gene_names, feature_types,
        mt_mask, ribo_mask,
        gene_passes, gene_n_cells, gene_sum, total_cells_passed,
    )
    print(
        f"  Genes (min_cells>={QC_PARAMS['min_cells_per_gene']}): "
        f"{n_genes_out:,} / {N_GENES:,}"
    )

    # ------------------------------------------------------------------
    # Build log-normalized sparse matrix from cached StreamResults
    # ------------------------------------------------------------------
    print("\nBuilding log-normalized sparse matrix ...")
    rows_all: list[np.ndarray] = []
    cols_all: list[np.ndarray] = []
    data_all: list[np.ndarray] = []
    cell_offset = 0
    cell_totals = obs_adata["total_counts"].to_numpy(dtype=np.float64)

    for i, (sample_id, time_label, abs_day) in enumerate(samples_to_run):
        result = cached_results[i]
        if result.n_pass == 0:
            continue

        # Convert to COO: X_csr rows=cells, cols=genes
        coo = result.X_csr.tocoo()
        cell_idx = coo.row   # 0-based cell index within this sample
        gene_idx = coo.col   # 0-based gene index

        # Apply gene filter
        g_ok = gene_rank[gene_idx] >= 0
        if not g_ok.any():
            cell_offset += result.n_pass
            del coo, cell_idx, gene_idx, g_ok
            gc.collect()
            continue

        global_cell_rows = cell_offset + cell_idx[g_ok]
        target_gene_cols = gene_rank[gene_idx[g_ok]]

        # Total-count normalize then log1p
        scale  = NORM_TARGET_SUM / cell_totals[global_cell_rows]
        values = np.log1p(
            coo.data[g_ok].astype(np.float64) * scale
        ).astype(np.float32)

        rows_all.append(global_cell_rows.astype(np.int32))
        cols_all.append(target_gene_cols.astype(np.int32))
        data_all.append(values)

        nnz_kept = int(g_ok.sum())
        print(f"  {sample_id}: nnz_filtered={nnz_kept:,}  RSS={rss_gb():.1f}GB")

        cell_offset += result.n_pass
        del coo, cell_idx, gene_idx, g_ok, global_cell_rows, target_gene_cols
        del scale, values
        gc.collect()

    if cell_offset != total_cells_passed:
        raise AssertionError(
            f"Cell count mismatch: accumulated {cell_offset} vs "
            f"expected {total_cells_passed}"
        )

    # Free cached StreamResults before allocating the final matrix
    del cached_results
    gc.collect()

    X = sp.coo_matrix(
        (
            np.concatenate(data_all),
            (np.concatenate(rows_all), np.concatenate(cols_all)),
        ),
        shape=(total_cells_passed, n_genes_out),
        dtype=np.float32,
    ).tocsr()
    del rows_all, cols_all, data_all
    gc.collect()

    # ------------------------------------------------------------------
    # Assemble AnnData and write h5ad
    # ------------------------------------------------------------------
    adata = ad.AnnData(X=X, obs=obs_adata, var=var)
    # 6 gene symbols in the genes TSV are duplicated (e.g. TBCE, MATR3).
    # Deduplicate by appending -1/-2 suffixes so anndata index is unique.
    adata.var_names_make_unique()
    adata.uns["dataset_id"]  = DATASET_ID
    adata.uns["time_axis"]   = "observed_time"
    adata.uns["preprocessing_summary"] = {
        "script":               "scripts/build_gse242424_raw_full_gene_input.py",
        "streaming_module":     "benchmark.shared.dataset.builders.mtx_streaming",
        "dataset":              DATASET_ID,
        "timestamp":            TIMESTAMP,
        "samples_built":        [s[0] for s in samples_to_run],
        "qc_params":            QC_PARAMS,
        "norm_target_sum":      NORM_TARGET_SUM,
        "n_reference_genes":    int(N_GENES),
        "n_genes_after_filter": int(n_genes_out),
        "x_transform":          "log1p(total-count normalized to 1e4)",
        "intended_use":         "full-gene marker annotation; not model training",
        "ipsc_abs_day_note":    "iPSC abs_day=16 is PROVISIONAL",
    }
    adata.write_h5ad(OUTPUT_H5AD)

    elapsed_total = _time.time() - T0
    print(f"\nDONE")
    print(f"  h5ad   : {OUTPUT_H5AD}")
    print(f"  Shape  : {adata.n_obs:,} x {adata.n_vars:,}")
    print(f"  Sparse : {sp.issparse(adata.X)}")
    print(f"  Time   : {elapsed_total/60:.1f} min  RSS={rss_gb():.2f}GB")

    # ------------------------------------------------------------------
    # QC summary CSV
    # ------------------------------------------------------------------
    qc_csv = REPORT_DIR / f"{TIMESTAMP}_gse242424_qc_summary.csv"
    pd.DataFrame(qc_rows).to_csv(qc_csv, index=False)
    print(f"  QC CSV : {qc_csv}")
    print()

    # Print summary table
    qc_df_out = pd.DataFrame(qc_rows)
    cols_show = ["sample_id", "time_label", "abs_day",
                 "n_pass_qc", "pass_qc_fraction",
                 "median_total_counts_pass", "median_n_genes_pass"]
    print(qc_df_out[cols_show].to_string(index=False))


# ===========================================================================
# Entry point
# ===========================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build GSE242424 raw/full-gene benchmark h5ad."
    )
    parser.add_argument(
        "--samples",
        metavar="GSM_ID[,GSM_ID,...]",
        default=None,
        help=(
            "Comma-separated sample IDs to process (e.g. GSM7763420_D2). "
            "Useful for smoke-testing without running all 9 samples. "
            "Omit to process the full manifest."
        ),
    )
    args = parser.parse_args()
    subset = [s.strip() for s in args.samples.split(",")] if args.samples else None
    main(samples_subset=subset)
