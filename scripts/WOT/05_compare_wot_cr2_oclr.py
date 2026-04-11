#!/usr/bin/env python3
"""
05_compare_wot_cr2_oclr.py  [v3 — OCLR dual-endpoint benchmark, margin metric]
==========================================================================================
Project : Comparative Study of Trajectory Inference Models for Chemical iPSC
          Reprogramming
Dataset : GSE230659 (human, Liuyang et al. 2023 Cell Stem Cell)
Step    : 5 — Compare WOT and CellRank2 fate predictions against OCLR dual-endpoint

Purpose
-------
This standalone post-processing script evaluates the predictive accuracy of:
  (A) WOT       : p_oclr_margin  = p_oclr_success − p_oclr_failure
                  (OCLR dual-endpoint design, script 02c v14+)
  (B) CellRank2 : cr2_oclr_margin = cr2_p_oclr_success − cr2_p_oclr_failure
                  (OCLR dual-endpoint annotation, script 04 v12+)

against the OCLR-defined dual-endpoint ground truth (script 02b2):
  SUCCESS = top-10% OCLR within hCiPSC (label = 1 in PRIMARY binary)
  FAILURE = bottom-10% OCLR within hCiPSC (label = 0 in PRIMARY binary)
  AMBIGUOUS = middle 80% — excluded from PRIMARY evaluation

Independence guarantee
----------------------
The OCLR endpoints are computed in 02b2 BEFORE WOT (02c) and CellRank2 (04) run.
The margin score is defined post-hoc using OCLR labels, not inputs to the models.
Evaluation is circularity-free.

NOTE: This script requires the hpc_full CellRank2 result.  A FileNotFoundError is
raised if the hpc_full h5ad is absent — there is no silent local_debug fallback.

Inputs  (auto-detected as most-recent file matching pattern)
-------
  results/metrics/*_GSE230659_wot_fates.h5ad          — WOT fate probs (full dataset)
  results/metrics/*_hpcfull_GSE230659_cellrank2.h5ad   — CellRank2 fate probs (hpc_full)
  data/processed/*_oclr_endpoint_labels_hcipsc.tsv     — OCLR dual-endpoint labels

Evaluation sets
---------------
  PRIMARY    : shared hCiPSC cells with non-ambiguous OCLR label
               (success=1 OR failure=0, ambiguous excluded)
               MAIN METRIC: AUROC/AUPRC/Spearman of margin score vs binary label
  SUPPLEMENTARY_success : same cells, p_oclr_success vs success=1 label
  SUPPLEMENTARY_all_wot  : all WOT cells
  SUPPLEMENTARY_all_cr2  : all CR2 cells

Metric types
------------
  binary_OCLR_endpoint   : AUROC, AUPRC, precision@k, recall@k
                           ground truth = success(1) vs failure(0)
                           score = p_oclr_margin or cr2_oclr_margin
  continuous_OCLR_score  : Spearman r vs continuous oclr_score

Outputs
-------
  results/cellrank2/{TIMESTAMP}_primary_endpoint_metrics.csv   — PRIMARY metrics
  results/cellrank2/{TIMESTAMP}_oclr_comparison_summary.tsv    — all eval sets
  results/cellrank2/{TIMESTAMP}_primary_endpoint_eval_cells.tsv — PRIMARY cell list
  results/cellrank2/{TIMESTAMP}_oclr_comparison_per_cell.csv.gz — full per-cell table
  results/cellrank2/{TIMESTAMP}_oclr_comparison_wot_vs_cr2.tsv  — margin correlation
  results/figures/{TIMESTAMP}_oclr_comparison_scatter.png

Changelog v3
-------------
- PRIMARY set: shared hCiPSC non-ambiguous cells (not shared non-hCiPSC).
- Main benchmark scores: p_oclr_margin (WOT) and cr2_oclr_margin (CR2).
- Add supplementary: p_oclr_success vs success-only binary label.
- New output files: *_primary_endpoint_metrics.csv, *_primary_endpoint_eval_cells.tsv.
- Load OCLR endpoint labels from *_oclr_endpoint_labels_hcipsc.tsv.

Changelog v2
-------------
- Remove PREFER_HPC flag; hard-error if hpc_full CellRank2 h5ad is absent.
- Column updates: p_iPSC → p_hcipsc, cr2_fate_ips → cr2_ips_fate_prob.
- Restructure metrics: PRIMARY (shared + non-hCiPSC) vs SUPPLEMENTARY (all cells).
- Add continuous OCLR score correlation (Spearman vs oclr_score).
- Output table has metric_type and eval_set columns.

Usage
-----
  conda activate traj_env
  python scripts/WOT/05_compare_wot_cr2_oclr.py
"""

# =============================================================================
# Imports
# =============================================================================
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# =============================================================================
# 0.  CONFIGURATION
# =============================================================================
PROJECT_ROOT   = Path(__file__).resolve().parents[2]
PROCESSED_DIR  = PROJECT_ROOT / "data" / "processed"
METRICS_DIR    = PROJECT_ROOT / "results" / "metrics"
CR2_RESULTS    = PROJECT_ROOT / "results" / "cellrank2"
FIGURES_DIR    = PROJECT_ROOT / "results" / "figures"
CR2_RESULTS.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")

# OCLR terminal sample (must match 02b2 / 02c)
OCLR_ENDPOINT_SAMPLE_ID = "GSM7230012_hCiPSCs-0618"

# Minimum precision-at-k threshold fractions to report
TOPK_FRACS: list = [0.05, 0.10, 0.20]

# =============================================================================
# 1.  HELPER FUNCTIONS
# =============================================================================

def find_latest(directory: Path, pattern: str) -> Path | None:
    """Return the most recently modified file matching *pattern* in *directory*."""
    candidates = sorted(directory.glob(pattern), key=lambda p: p.stat().st_mtime)
    return candidates[-1] if candidates else None


def load_oclr_ground_truth(adata_obs: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    """
    Derive OCLR ground-truth binary labels for cells in adata_obs.

    Priority order:
      1. oclr_ips_subset_mask column already in adata_obs  (new pipeline)
      2. *_oclr_endpoint_barcodes.txt file in PROCESSED_DIR  (02b2 output)
      3. marker_ips_subset_mask column in adata_obs  (legacy — warn user)

    Returns
    -------
    gt_binary  : pd.Series[int]   (0/1, index = adata_obs.index)
    oclr_score : pd.Series[float] (NaN where unavailable)
    """
    barcodes = adata_obs.index
    gt_binary = pd.Series(np.zeros(len(barcodes), dtype=int), index=barcodes)

    # --- Option 1: pre-computed mask column --------------------------------
    if "oclr_ips_subset_mask" in adata_obs.columns:
        print("  Ground truth source: oclr_ips_subset_mask column (new pipeline)")
        gt_binary  = adata_obs["oclr_ips_subset_mask"].astype(int)
        oclr_score = adata_obs.get("oclr_score", pd.Series(np.nan, index=barcodes))
        return gt_binary, oclr_score

    # --- Option 2: barcode file from 02b2 ----------------------------------
    bc_files = sorted(PROCESSED_DIR.glob("*_oclr_endpoint_barcodes.txt"),
                      key=lambda p: p.stat().st_mtime)
    if bc_files:
        bc_file = bc_files[-1]
        print(f"  Ground truth source: {bc_file.name}")

        try:
            bc_df = pd.read_csv(bc_file, sep="\t", index_col=0, dtype=str)
            bc_df.index = bc_df.index.str.strip()
            if bc_df.index[0].lower() == "barcode":
                bc_df = pd.read_csv(bc_file, sep="\t", header=None,
                                    names=["barcode", "oclr_score"], dtype=str,
                                    skiprows=0)
                bc_df = bc_df[bc_df["barcode"].str.lower() != "barcode"]
                bc_df = bc_df.set_index("barcode")
                bc_df.index = bc_df.index.str.strip()
        except Exception as _parse_err:
            print(f"    [WARN] TSV parse failed ({_parse_err}); "
                  f"falling back to single-column read")
            bc_df = pd.read_csv(bc_file, sep="\t", header=None, dtype=str,
                                 usecols=[0], names=["barcode"])
            bc_df = bc_df[bc_df["barcode"].str.lower() != "barcode"]
            bc_df = bc_df.set_index("barcode")
            bc_df.index = bc_df.index.str.strip()

        endpoint_barcodes = set(bc_df.index.tolist())
        _n_total = len(endpoint_barcodes)
        print(f"    File contains {_n_total:,} terminal barcodes")

        oclr_score = pd.Series(np.nan, index=barcodes, dtype=float)
        if "oclr_score" in bc_df.columns:
            _scores_from_bc = (bc_df["oclr_score"]
                               .apply(lambda x: float(x) if x not in ("nan", "", "None") else np.nan))
            _scores_from_bc = _scores_from_bc.reindex(barcodes)
            oclr_score.update(_scores_from_bc.dropna())
            print(f"    oclr_score column loaded from barcode file "
                  f"({_scores_from_bc.notna().sum():,} values)")
        else:
            score_files = sorted(PROCESSED_DIR.glob("*_oclr_score_all_hcipsc.csv"),
                                 key=lambda p: p.stat().st_mtime)
            if score_files:
                score_df = pd.read_csv(score_files[-1], index_col=0)
                score_col_name = next(
                    (c for c in score_df.columns if "oclr_score" in c.lower()), None
                ) or next(
                    (c for c in score_df.columns if "score" in c.lower()), None
                )
                if score_col_name:
                    _scores = score_df[score_col_name].reindex(barcodes)
                    oclr_score.update(_scores.dropna())
                    print(f"    oclr_score loaded from {score_files[-1].name} "
                          f"({_scores.notna().sum():,} values)")

        gt_binary = pd.Series(
            barcodes.isin(endpoint_barcodes).astype(int), index=barcodes
        )
        _n_match = int(gt_binary.sum())
        print(f"    Matched {_n_match:,} / {_n_total:,} terminal barcodes in h5ad "
              f"({_n_match/_n_total*100:.1f}%)")
        if _n_match == 0:
            print(f"    [WARN] Zero barcodes matched. Sample endpoint barcodes:")
            for _b in list(endpoint_barcodes)[:3]:
                print(f"      {_b!r}")
            print(f"    Sample h5ad barcodes:")
            for _b in list(barcodes)[:3]:
                print(f"      {_b!r}")
        elif _n_match < _n_total * 0.50:
            print(f"    [WARN] Less than 50% of OCLR barcodes found in h5ad. "
                  f"Ensure the full-dataset WOT h5ad is being used.")
        return gt_binary, oclr_score

    # --- Option 3: legacy marker mask (warn) --------------------------------
    if "marker_ips_subset_mask" in adata_obs.columns:
        print("  [WARN] Ground truth source: marker_ips_subset_mask (LEGACY — "
              "run 02b2 + 02c to generate OCLR endpoint)")
        gt_binary  = adata_obs["marker_ips_subset_mask"].astype(int)
        oclr_score = pd.Series(np.nan, index=barcodes, dtype=float)
        return gt_binary, oclr_score

    raise RuntimeError(
        "Cannot determine OCLR ground truth: no oclr_ips_subset_mask column, "
        "no *_oclr_endpoint_barcodes.txt in data/processed/, "
        "and no legacy marker_ips_subset_mask column. "
        "Run 02b2_define_oclr_endpoint.py first."
    )


def compute_metrics(
    y_true: np.ndarray,
    y_score: np.ndarray,
    label: str,
    topk_fracs: list = TOPK_FRACS,
) -> dict:
    """
    Compute classification + ranking metrics for a continuous score against binary
    OCLR terminal ground truth (metric_type = binary_OCLR_terminal).

    Metrics
    -------
    auroc                    : Area under ROC curve
    auprc                    : Area under precision-recall curve
    pearson_r / p            : Pearson correlation vs binary label
    spearman_r / p           : Spearman correlation vs binary label
    precision_at_top{X}pct  : Precision among top k = round(X/100 * n_cells) cells
    recall_at_top{X}pct     : Recall at the same cutoff k.
    """
    from scipy.stats import pearsonr, spearmanr
    from sklearn.metrics import average_precision_score, roc_auc_score

    finite_mask = np.isfinite(y_score)
    yt = y_true[finite_mask]
    ys = y_score[finite_mask]
    n_pos = int(yt.sum())

    row: dict = {
        "model": label,
        "metric_type": "binary_OCLR_terminal",
        "n_cells": int(finite_mask.sum()),
        "n_pos": n_pos,
    }

    if n_pos == 0 or n_pos == len(yt):
        print(f"  [WARN] {label}: degenerate ground truth (n_pos={n_pos}/{len(yt)})")
        return row

    try:
        row["auroc"]  = float(roc_auc_score(yt, ys))
    except Exception as e:
        row["auroc"] = np.nan
        print(f"  [WARN] {label} AUROC failed: {e}")

    try:
        row["auprc"] = float(average_precision_score(yt, ys))
    except Exception as e:
        row["auprc"] = np.nan
        print(f"  [WARN] {label} AUPRC failed: {e}")

    try:
        pr, pp = pearsonr(yt.astype(float), ys)
        row["pearson_r"] = float(pr)
        row["pearson_p"] = float(pp)
    except Exception:
        row["pearson_r"] = row["pearson_p"] = np.nan

    try:
        sr, sp = spearmanr(yt.astype(float), ys)
        row["spearman_r"] = float(sr)
        row["spearman_p"] = float(sp)
    except Exception:
        row["spearman_r"] = row["spearman_p"] = np.nan

    sorted_idx = np.argsort(ys)[::-1]
    for frac in topk_fracs:
        k = max(1, int(round(frac * len(yt))))
        topk_true = yt[sorted_idx[:k]]
        row[f"precision_at_top{int(frac*100)}pct"] = float(topk_true.mean())
        row[f"recall_at_top{int(frac*100)}pct"]    = float(topk_true.sum() / n_pos)

    return row


def compute_oclr_score_correlation(
    fate_prob: np.ndarray,
    oclr_score: np.ndarray,
    label: str,
) -> dict:
    """
    Compute Spearman correlation of fate_prob vs continuous OCLR score
    (metric_type = continuous_OCLR_score).

    Only cells where both fate_prob and oclr_score are finite are used.
    """
    from scipy.stats import spearmanr

    finite = np.isfinite(fate_prob) & np.isfinite(oclr_score)
    row: dict = {
        "model": label,
        "metric_type": "continuous_OCLR_score",
        "n_cells": int(finite.sum()),
        "n_pos": np.nan,
    }
    if finite.sum() < 10:
        print(f"  [WARN] {label}: too few finite (fate_prob, oclr_score) pairs "
              f"({finite.sum()}) — Spearman skipped")
        return row

    try:
        sr, sp = spearmanr(fate_prob[finite], oclr_score[finite])
        row["spearman_r"] = float(sr)
        row["spearman_p"] = float(sp)
    except Exception as e:
        row["spearman_r"] = row["spearman_p"] = np.nan
        print(f"  [WARN] {label} continuous Spearman failed: {e}")

    return row


# =============================================================================
# 2.  MAIN
# =============================================================================
def main():
    print("=" * 70)
    print("05_compare_wot_cr2_oclr.py  v3  —  WOT vs CellRank2 OCLR dual-endpoint benchmark")
    print(f"  Timestamp : {TIMESTAMP}")
    print(f"  Project   : {PROJECT_ROOT}")
    print("=" * 70)

    # ── 2.1  Locate input files ────────────────────────────────────────────────
    print("\n[1] Locating input files ...")

    wot_h5ad_path = find_latest(METRICS_DIR, "*_GSE230659_wot_fates.h5ad")
    if wot_h5ad_path is None:
        raise FileNotFoundError(
            f"No *_GSE230659_wot_fates.h5ad found in {METRICS_DIR}. "
            "Run 02c_compute_transport.py first."
        )
    print(f"  WOT h5ad    : {wot_h5ad_path.name}")

    # CellRank2: hpc_full ONLY — no silent local_debug fallback.
    # Local-debug results are a stratified subset for pipeline smoke-testing;
    # they cannot be used for the final WOT vs CellRank2 comparison.
    cr2_h5ad_path = find_latest(METRICS_DIR, "*_hpcfull_GSE230659_cellrank2.h5ad")
    if cr2_h5ad_path is None:
        raise FileNotFoundError(
            f"No *_hpcfull_GSE230659_cellrank2.h5ad found in {METRICS_DIR}.\n"
            "Script 05 requires the full-dataset CellRank2 result.\n"
            "Run 04_cellrank2.py with RUN_MODE=hpc_full on Shirokane, then "
            "re-run this script.\n"
            "Local-debug (stratified subsample) results are intentionally NOT "
            "accepted here — use them only for pipeline smoke-testing."
        )

    print(f"  CR2 h5ad    : {cr2_h5ad_path.name}")

    # ── 2.2  Load h5ad files ──────────────────────────────────────────────────
    print("\n[2] Loading h5ad files ...")
    import anndata as ad

    adata_wot = ad.read_h5ad(str(wot_h5ad_path))
    print(f"  WOT  : {adata_wot.shape[0]:,} cells × {adata_wot.shape[1]:,} genes")
    print(f"         obs cols: {sorted(adata_wot.obs.columns.tolist())}")

    adata_cr2 = ad.read_h5ad(str(cr2_h5ad_path))
    print(f"  CR2  : {adata_cr2.shape[0]:,} cells × {adata_cr2.shape[1]:,} genes")
    print(f"         obs cols: {sorted(adata_cr2.obs.columns.tolist())}")

    # ── 2.3  Validate required columns ───────────────────────────────────────
    # PRIMARY columns: dual-endpoint margins (v14 / v12 pipeline)
    # Fallback: legacy single-endpoint columns (v13 / v11 pipeline)
    _wot_primary_col   = ("p_oclr_margin"   if "p_oclr_margin"   in adata_wot.obs.columns
                          else "p_hcipsc")
    _cr2_primary_col   = ("cr2_oclr_margin" if "cr2_oclr_margin" in adata_cr2.obs.columns
                          else "cr2_ips_fate_prob")
    _wot_success_col   = ("p_oclr_success"  if "p_oclr_success"  in adata_wot.obs.columns
                          else "p_hcipsc")
    _cr2_success_col   = ("cr2_p_oclr_success" if "cr2_p_oclr_success" in adata_cr2.obs.columns
                          else "cr2_ips_fate_prob")

    # Fatal if neither margin nor legacy columns exist
    _wot_req  = [_wot_primary_col]
    _cr2_req  = [_cr2_primary_col]
    _wot_miss = [c for c in _wot_req if c not in adata_wot.obs.columns]
    _cr2_miss = [c for c in _cr2_req if c not in adata_cr2.obs.columns]
    if _wot_miss:
        raise ValueError(
            f"WOT h5ad missing required columns: {_wot_miss}\n"
            "Ensure 02c_compute_transport.py v14+ (OCLR dual-endpoint) has run.\n"
            "Fallback p_hcipsc also not found."
        )
    if _cr2_miss:
        raise ValueError(
            f"CR2 h5ad missing required columns: {_cr2_miss}\n"
            "Ensure 04_cellrank2.py v12+ (OCLR dual-endpoint annotation) has run.\n"
            "Fallback cr2_ips_fate_prob also not found."
        )
    print(f"  WOT primary column  : {_wot_primary_col}")
    print(f"  CR2 primary column  : {_cr2_primary_col}")
    print(f"  WOT success column  : {_wot_success_col}")
    print(f"  CR2 success column  : {_cr2_success_col}")

    # ── 2.4  Load OCLR endpoint labels ────────────────────────────────────────
    print("\n[3] Loading OCLR dual-endpoint labels ...")
    # Prefer authoritative TSV from 02b2 v14+
    _label_files = sorted(PROCESSED_DIR.glob("*_oclr_endpoint_labels_hcipsc.tsv"),
                          key=lambda p: p.stat().st_mtime)
    if _label_files:
        _lf = _label_files[-1]
        print(f"  OCLR endpoint labels : {_lf.name}")
        _ldf = pd.read_csv(str(_lf), sep="\t", index_col=0)
        _ldf.index = _ldf.index.astype(str)
        _ep_label_global = _ldf.get("oclr_endpoint_label", pd.Series(dtype=str))
        _suc_set_global  = set(_ldf.index[
            _ldf.get("is_oclr_success_endpoint",
                     _ldf.get("is_oclr_terminal", pd.Series(0, index=_ldf.index))) == 1
        ].tolist())
        _fail_set_global = set(_ldf.index[
            _ldf["is_oclr_failure_endpoint"] == 1
        ].tolist()) if "is_oclr_failure_endpoint" in _ldf.columns else set()
        _oclr_score_map  = _ldf["oclr_score"].to_dict() if "oclr_score" in _ldf.columns else {}
        print(f"    SUCCESS: {len(_suc_set_global):,}  FAILURE: {len(_fail_set_global):,}")
    else:
        # Fall back to legacy barcode files
        print("  [FALLBACK] No *_oclr_endpoint_labels_hcipsc.tsv found — "
              "loading legacy barcode files")
        _suc_files  = sorted(PROCESSED_DIR.glob("*_oclr_success_endpoint_barcodes.txt"),
                             key=lambda p: p.stat().st_mtime)
        _fail_files = sorted(PROCESSED_DIR.glob("*_oclr_failure_endpoint_barcodes.txt"),
                             key=lambda p: p.stat().st_mtime)
        _suc_set_global  = set()
        _fail_set_global = set()
        _oclr_score_map  = {}
        if _suc_files:
            _s = pd.read_csv(str(_suc_files[-1]), sep="\t", index_col=0)
            _suc_set_global = set(_s.index.astype(str).tolist())
        if _fail_files:
            _f = pd.read_csv(str(_fail_files[-1]), sep="\t", index_col=0)
            _fail_set_global = set(_f.index.astype(str).tolist())
        _score_sf = sorted(PROCESSED_DIR.glob("*_oclr_score_all_hcipsc.csv"),
                           key=lambda p: p.stat().st_mtime)
        if _score_sf:
            _sf_df = pd.read_csv(str(_score_sf[-1]), index_col=0)
            if "oclr_score" in _sf_df.columns:
                _oclr_score_map = _sf_df["oclr_score"].to_dict()
        print(f"    SUCCESS: {len(_suc_set_global):,}  FAILURE: {len(_fail_set_global):,}")

    # ── Build per-cell endpoint label series for WOT and CR2 ────────────────
    def _build_labels(obs_names):
        bcs = obs_names.astype(str)
        suc  = np.array([bc in _suc_set_global  for bc in bcs], dtype=bool)
        fail = np.array([bc in _fail_set_global for bc in bcs], dtype=bool)
        # Binary for PRIMARY: success=1, failure=0, ambiguous=NaN
        binary = np.where(suc, 1.0, np.where(fail, 0.0, np.nan))
        oclr_sc = np.array([_oclr_score_map.get(bc, np.nan) for bc in bcs], dtype=float)
        return (
            pd.Series(binary, index=obs_names, name="oclr_endpoint_binary"),
            pd.Series(suc.astype(int), index=obs_names, name="is_oclr_success"),
            pd.Series(fail.astype(int), index=obs_names, name="is_oclr_failure"),
            pd.Series(oclr_sc, index=obs_names, name="oclr_score"),
        )

    (gt_binary_wot, gt_suc_wot, gt_fail_wot, oclr_score_wot) = _build_labels(
        adata_wot.obs_names)
    (gt_binary_cr2, gt_suc_cr2, gt_fail_cr2, oclr_score_cr2) = _build_labels(
        adata_cr2.obs_names)

    # Legacy aliases (binary ground truth = success=1 vs rest for supplementary)
    gt_wot = gt_suc_wot   # success vs all (for supplementary)
    gt_cr2 = gt_suc_cr2

    n_suc_wot  = int(gt_suc_wot.sum())
    n_fail_wot = int(gt_fail_wot.sum())
    print(f"  WOT: success={n_suc_wot:,}  failure={n_fail_wot:,}  "
          f"ambiguous/unscored={len(adata_wot.obs_names)-n_suc_wot-n_fail_wot:,}")

    # ── 2.5  Identify shared barcodes and PRIMARY eval set ────────────────────
    shared_bc = adata_wot.obs_names.intersection(adata_cr2.obs_names)
    n_shared  = len(shared_bc)
    print(f"\n[4] Shared barcodes: {n_shared:,}")

    if n_shared < 50:
        raise ValueError(
            f"Only {n_shared} shared barcodes between WOT and CellRank2 h5ads. "
            "Cannot compute the PRIMARY (shared-cells) evaluation. "
            "Ensure both h5ads come from the full dataset run."
        )

    # PRIMARY set: shared cells that are hCiPSC AND have non-ambiguous OCLR label
    # (i.e., gt_binary is 0 or 1, not NaN)
    _gt_binary_shared = gt_binary_wot.reindex(shared_bc)
    _non_ambig_shared = _gt_binary_shared.notna().values
    primary_bc        = shared_bc[_non_ambig_shared]
    n_primary         = len(primary_bc)
    _n_suc_prim  = int((_gt_binary_shared.reindex(primary_bc) == 1).sum())
    _n_fail_prim = int((_gt_binary_shared.reindex(primary_bc) == 0).sum())
    print(f"  PRIMARY eval set: {n_primary:,} shared hCiPSC non-ambiguous cells")
    print(f"    success (label=1): {_n_suc_prim:,}  failure (label=0): {_n_fail_prim:,}")

    if n_primary < 50:
        raise ValueError(
            f"PRIMARY eval set has only {n_primary} cells "
            f"(success={_n_suc_prim}, failure={_n_fail_prim}). "
            "Ensure 02b2 generated dual endpoints and the hCiPSC barcodes are present in both h5ads."
        )
    if _n_suc_prim == 0 or _n_fail_prim == 0:
        print(f"  [WARN] PRIMARY is degenerate: success={_n_suc_prim}, failure={_n_fail_prim}. "
              "AUROC/AUPRC will be undefined.")

    # ── 2.6  Compute evaluation metrics ──────────────────────────────────────
    print("\n[5] Computing evaluation metrics ...")
    records = []

    # ── PRIMARY: shared non-ambiguous hCiPSC cells, margin score ──────────────
    # Ground truth: success=1, failure=0 (ambiguous cells already excluded)
    _gt_primary       = _gt_binary_shared.reindex(primary_bc).values.astype(int)
    _wot_margin_prim  = adata_wot.obs.loc[primary_bc, _wot_primary_col].values.astype(float)
    _cr2_margin_prim  = adata_cr2.obs.loc[primary_bc, _cr2_primary_col].values.astype(float)
    _oclr_sc_primary  = oclr_score_wot.reindex(primary_bc).values.astype(float)

    print(f"\n  -- PRIMARY (shared non-ambiguous hCiPSC, n={n_primary:,}) --")
    _r_wot_prim = compute_metrics(
        _gt_primary, _wot_margin_prim, "WOT_primary_margin")
    _r_wot_prim["eval_set"]  = "primary_shared_hcipsc"
    _r_wot_prim["score_col"] = _wot_primary_col
    records.append(_r_wot_prim)

    _r_cr2_prim = compute_metrics(
        _gt_primary, _cr2_margin_prim, "CR2_primary_margin")
    _r_cr2_prim["eval_set"]  = "primary_shared_hcipsc"
    _r_cr2_prim["score_col"] = _cr2_primary_col
    records.append(_r_cr2_prim)

    # Continuous OCLR score correlation — PRIMARY (margin vs raw OCLR score)
    _rc_wot_prim = compute_oclr_score_correlation(
        _wot_margin_prim, _oclr_sc_primary, "WOT_primary_margin")
    _rc_wot_prim["eval_set"] = "primary_shared_hcipsc"
    records.append(_rc_wot_prim)

    _rc_cr2_prim = compute_oclr_score_correlation(
        _cr2_margin_prim, _oclr_sc_primary, "CR2_primary_margin")
    _rc_cr2_prim["eval_set"] = "primary_shared_hcipsc"
    records.append(_rc_cr2_prim)

    # ── SUPPLEMENTARY_success: same PRIMARY cells, p_oclr_success vs success ──
    # Uses success-probability column vs success=1 / failure=0 binary label
    _wot_suc_prim = adata_wot.obs.loc[primary_bc, _wot_success_col].values.astype(float)
    _cr2_suc_prim = adata_cr2.obs.loc[primary_bc, _cr2_success_col].values.astype(float)

    print(f"\n  -- SUPPLEMENTARY_success (shared non-ambiguous hCiPSC, n={n_primary:,}) --")
    _r_wot_supp_s = compute_metrics(
        _gt_primary, _wot_suc_prim, "WOT_supplementary_success")
    _r_wot_supp_s["eval_set"]  = "supplementary_success"
    _r_wot_supp_s["score_col"] = _wot_success_col
    records.append(_r_wot_supp_s)

    _r_cr2_supp_s = compute_metrics(
        _gt_primary, _cr2_suc_prim, "CR2_supplementary_success")
    _r_cr2_supp_s["eval_set"]  = "supplementary_success"
    _r_cr2_supp_s["score_col"] = _cr2_success_col
    records.append(_r_cr2_supp_s)

    # ── SUPPLEMENTARY: all WOT cells, margin vs success-binary ────────────────
    print(f"\n  -- SUPPLEMENTARY (all WOT cells, n={adata_wot.n_obs:,}) --")
    _wot_all_margin  = adata_wot.obs[_wot_primary_col].values.astype(float)
    _oclr_sc_wot_all = oclr_score_wot.values.astype(float)

    _r_wot_all = compute_metrics(
        gt_suc_wot.values.astype(int), _wot_all_margin, "WOT_all")
    _r_wot_all["eval_set"]  = "supplementary_all_wot"
    _r_wot_all["score_col"] = _wot_primary_col
    records.append(_r_wot_all)

    _rc_wot_all = compute_oclr_score_correlation(
        _wot_all_margin, _oclr_sc_wot_all, "WOT_all")
    _rc_wot_all["eval_set"] = "supplementary_all_wot"
    records.append(_rc_wot_all)

    # ── SUPPLEMENTARY: all CR2 cells, margin vs success-binary ────────────────
    print(f"\n  -- SUPPLEMENTARY (all CR2 cells, n={adata_cr2.n_obs:,}) --")
    _cr2_all_margin = adata_cr2.obs[_cr2_primary_col].values.astype(float)
    _oclr_sc_cr2    = oclr_score_cr2.values.astype(float)

    _r_cr2_all = compute_metrics(
        gt_suc_cr2.values.astype(int), _cr2_all_margin, "CR2_all")
    _r_cr2_all["eval_set"]  = "supplementary_all_cr2"
    _r_cr2_all["score_col"] = _cr2_primary_col
    records.append(_r_cr2_all)

    _rc_cr2_all = compute_oclr_score_correlation(
        _cr2_all_margin, _oclr_sc_cr2, "CR2_all")
    _rc_cr2_all["eval_set"] = "supplementary_all_cr2"
    records.append(_rc_cr2_all)

    # Print summary
    metrics_df = pd.DataFrame(records)
    _print_cols = ["model", "metric_type", "eval_set", "score_col",
                   "n_cells", "auroc", "auprc", "spearman_r", "pearson_r"]
    _print_cols = [c for c in _print_cols if c in metrics_df.columns]
    print("\n  Metrics summary vs OCLR ground truth:")
    print(metrics_df[_print_cols].to_string(index=False))

    # ── 2.7  WOT vs CellRank2 margin agreement (shared cells) ───────────────────
    # Correlate WOT margin (p_oclr_margin) vs CR2 margin (cr2_oclr_margin)
    # on ALL shared cells and on the PRIMARY subset.
    print("\n[6] WOT vs CellRank2 margin agreement (shared cells) ...")
    from scipy.stats import pearsonr, spearmanr

    # Build margin vectors for all shared cells
    _wot_margin_shared = adata_wot.obs.loc[shared_bc, _wot_primary_col].values.astype(float)
    _cr2_margin_shared = adata_cr2.obs.loc[shared_bc, _cr2_primary_col].values.astype(float)

    wot_vs_cr2_rows = []
    _fin_shared = np.isfinite(_wot_margin_shared) & np.isfinite(_cr2_margin_shared)
    if _fin_shared.sum() < 50:
        print("  [WARN] Fewer than 50 finite shared-cell margin pairs — correlation skipped")
    else:
        pr_all, pp_all = pearsonr(_wot_margin_shared[_fin_shared],
                                  _cr2_margin_shared[_fin_shared])
        sr_all, sp_all = spearmanr(_wot_margin_shared[_fin_shared],
                                   _cr2_margin_shared[_fin_shared])
        wot_vs_cr2_rows += [
            {"metric": "score_col_wot",  "value": _wot_primary_col},
            {"metric": "score_col_cr2",  "value": _cr2_primary_col},
            {"metric": "pearson_r_all",  "value": pr_all},
            {"metric": "pearson_p_all",  "value": pp_all},
            {"metric": "spearman_r_all", "value": sr_all},
            {"metric": "spearman_p_all", "value": sp_all},
            {"metric": "n_cells_all",    "value": int(_fin_shared.sum())},
        ]
        print(f"  WOT col: {_wot_primary_col}   CR2 col: {_cr2_primary_col}")
        print(f"  Pearson  r (all shared)    = {pr_all:.4f}  (p={pp_all:.2e})")
        print(f"  Spearman r (all shared)    = {sr_all:.4f}  (p={sp_all:.2e})")

        # PRIMARY set correlation (non-ambiguous hCiPSC)
        _fin_primary = np.isfinite(_wot_margin_prim) & np.isfinite(_cr2_margin_prim)
        if _fin_primary.sum() > 50:
            pr_p, pp_p = pearsonr(_wot_margin_prim[_fin_primary],
                                  _cr2_margin_prim[_fin_primary])
            sr_p, sp_p = spearmanr(_wot_margin_prim[_fin_primary],
                                   _cr2_margin_prim[_fin_primary])
            wot_vs_cr2_rows += [
                {"metric": "pearson_r_primary_hcipsc",  "value": pr_p},
                {"metric": "pearson_p_primary_hcipsc",  "value": pp_p},
                {"metric": "spearman_r_primary_hcipsc", "value": sr_p},
                {"metric": "spearman_p_primary_hcipsc", "value": sp_p},
                {"metric": "n_cells_primary_hcipsc",    "value": int(_fin_primary.sum())},
            ]
            print(f"  Pearson  r (PRIMARY hCiPSC) = {pr_p:.4f}  (p={pp_p:.2e})")
            print(f"  Spearman r (PRIMARY hCiPSC) = {sr_p:.4f}  (p={sp_p:.2e})")

    # ── 2.8  Per-cell predictions table ──────────────────────────────────────
    # Primary score columns: p_oclr_margin (WOT) and cr2_oclr_margin (CR2).
    # Legacy aliases (p_hcipsc, cr2_ips_fate_prob) included if present.
    print("\n[7] Building per-cell predictions table ...")

    # WOT columns — prefer new dual-endpoint columns, fall back to legacy
    _wot_obs = adata_wot.obs
    per_cell_df = pd.DataFrame({
        "p_oclr_margin":       _wot_obs.get(_wot_primary_col,
                               pd.Series(np.nan, index=adata_wot.obs_names)),
        "p_oclr_success":      _wot_obs.get("p_oclr_success",
                               pd.Series(np.nan, index=adata_wot.obs_names)),
        "p_oclr_failure":      _wot_obs.get("p_oclr_failure",
                               pd.Series(np.nan, index=adata_wot.obs_names)),
        "oclr_endpoint_binary": gt_binary_wot,   # success=1, failure=0, ambiguous=NaN
        "is_oclr_success":     gt_suc_wot,
        "is_oclr_failure":     gt_fail_wot,
        "oclr_score":          oclr_score_wot,
    }, index=adata_wot.obs_names)

    # Legacy alias p_hcipsc
    if "p_hcipsc" in _wot_obs.columns:
        per_cell_df["p_hcipsc"] = _wot_obs["p_hcipsc"].values

    # CR2 columns — join by barcode (left join keeps all WOT cells)
    if n_shared > 0:
        for _cr2_col_name in ["cr2_oclr_margin", "cr2_p_oclr_success",
                               "cr2_p_oclr_failure", "cr2_ips_fate_prob",
                               "cr2_ips_lineage_name", "cr2_success_lineage_name",
                               "cr2_failure_lineage_name"]:
            if _cr2_col_name in adata_cr2.obs.columns:
                _s = pd.Series(adata_cr2.obs[_cr2_col_name].values,
                               index=adata_cr2.obs_names,
                               name=_cr2_col_name)
                per_cell_df = per_cell_df.join(_s, how="left")

    # Metadata columns from WOT h5ad
    for _meta in ["sample_id", "stage", "stage_std", "abs_day", "stage_day_label"]:
        if _meta in _wot_obs.columns:
            per_cell_df[_meta] = _wot_obs[_meta].values

    print(f"  Per-cell table shape : {per_cell_df.shape}")
    print(f"  Columns              : {per_cell_df.columns.tolist()}")

    # ── 2.9  Save outputs ─────────────────────────────────────────────────────
    print("\n[8] Saving outputs ...")

    # a) PRIMARY endpoint metrics (main deliverable)
    _prim_mask = metrics_df["eval_set"] == "primary_shared_hcipsc"
    out_prim_metrics = CR2_RESULTS / f"{TIMESTAMP}_primary_endpoint_metrics.csv"
    metrics_df[_prim_mask].to_csv(out_prim_metrics, index=False)
    print(f"  PRIMARY metrics  → {out_prim_metrics.name}")

    # b) PRIMARY eval cell list (barcode + scores + ground truth)
    _prim_cell_cols = [c for c in ["p_oclr_margin", "p_oclr_success", "p_oclr_failure",
                                    "cr2_oclr_margin", "cr2_p_oclr_success",
                                    "cr2_p_oclr_failure", "oclr_score",
                                    "oclr_endpoint_binary", "is_oclr_success",
                                    "is_oclr_failure", "sample_id", "stage"]
                       if c in per_cell_df.columns]
    out_prim_cells = CR2_RESULTS / f"{TIMESTAMP}_primary_endpoint_eval_cells.tsv"
    per_cell_df.loc[primary_bc, _prim_cell_cols].to_csv(out_prim_cells, sep="\t")
    print(f"  PRIMARY cells    → {out_prim_cells.name}  (n={len(primary_bc):,})")

    # c) Full metrics summary (all eval sets)
    out_summary = CR2_RESULTS / f"{TIMESTAMP}_oclr_comparison_summary.tsv"
    _summary_df = metrics_df.copy()
    _summary_df.insert(0, "wot_n_cells_total", adata_wot.shape[0])
    _summary_df.insert(1, "cr2_n_cells_total", adata_cr2.shape[0])
    _summary_df.insert(2, "n_shared_cells",     n_shared)
    _summary_df.insert(3, "n_primary_cells",    n_primary)
    _summary_df.to_csv(out_summary, sep="\t", index=False)
    print(f"  Metrics summary  → {out_summary.name}")

    # d) WOT vs CR2 margin correlation
    out_corr = CR2_RESULTS / f"{TIMESTAMP}_oclr_comparison_wot_vs_cr2.tsv"
    if wot_vs_cr2_rows:
        _corr_df = pd.DataFrame(wot_vs_cr2_rows)
        _corr_df.insert(0, "wot_n_cells",  adata_wot.shape[0])
        _corr_df.insert(1, "cr2_n_cells",  adata_cr2.shape[0])
        _corr_df.to_csv(out_corr, sep="\t", index=False)
        print(f"  WOT vs CR2 corr  → {out_corr.name}")

    # e) Full per-cell predictions table (gzipped)
    out_percell = CR2_RESULTS / f"{TIMESTAMP}_oclr_comparison_per_cell.csv.gz"
    per_cell_df.to_csv(out_percell, compression="gzip")
    print(f"  Per-cell table   → {out_percell.name}")

    # ── 2.10  Diagnostic figure ───────────────────────────────────────────────
    print("\n[9] Generating comparison figure ...")

    # Determine which CR2 column is available for the figure
    _cr2_fig_col = ("cr2_oclr_margin" if "cr2_oclr_margin" in per_cell_df.columns
                    else ("cr2_ips_fate_prob" if "cr2_ips_fate_prob" in per_cell_df.columns
                          else None))

    _metrics_indexed = metrics_df.set_index("model")
    _make_comparison_figure(
        per_cell_df        = per_cell_df,
        metrics_df         = _metrics_indexed,
        gt_col             = "oclr_endpoint_binary",   # success=1, failure=0, ambiguous=NaN
        wot_col            = "p_oclr_margin",          # primary WOT score
        cr2_col            = _cr2_fig_col,             # primary CR2 score (margin or legacy)
        score_col          = "oclr_score",
        wot_metrics_label  = "WOT_primary_margin",
        cr2_metrics_label  = "CR2_primary_margin",
    )

    print("\n" + "=" * 70)
    print("DONE")
    print(f"  wot_n_cells={adata_wot.n_obs:,}  cr2_n_cells={adata_cr2.n_obs:,}")
    print(f"  n_shared={n_shared:,}  n_primary={n_primary:,}")
    print(f"  PRIMARY (non-ambiguous hCiPSC): {n_primary:,} cells  "
          f"(success={_n_suc_prim:,}, failure={_n_fail_prim:,})")
    print("=" * 70)


# =============================================================================
# 3.  FIGURE
# =============================================================================

def _make_comparison_figure(
    per_cell_df: pd.DataFrame,
    metrics_df: pd.DataFrame,
    gt_col: str,
    wot_col: str,
    cr2_col: str | None,
    score_col: str,
    wot_metrics_label: str = "",
    cr2_metrics_label: str = "",
) -> None:
    """Generate a 2×3 comparison figure for the OCLR dual-endpoint benchmark.

    Uses PRIMARY rows (success/failure non-ambiguous hCiPSC) for ROC/PR curves.
    gt_col must contain 1=success, 0=failure, NaN=ambiguous; NaN rows are excluded
    from ROC/PR.

    wot_metrics_label / cr2_metrics_label : index into metrics_df for AUROC/AUPRC.
    """
    from sklearn.metrics import precision_recall_curve, roc_curve

    # Resolve metrics index entries
    _wot_idx = (wot_metrics_label if wot_metrics_label and wot_metrics_label in metrics_df.index
                else next((i for i in metrics_df.index if i.startswith("WOT_primary")), None))
    _cr2_idx = (cr2_metrics_label if cr2_metrics_label and cr2_metrics_label in metrics_df.index
                else next((i for i in metrics_df.index if i.startswith("CR2_primary")), None))

    has_cr2 = (cr2_col is not None) and (cr2_col in per_cell_df.columns)

    # For ROC/PR: restrict to non-ambiguous rows (gt finite = success or failure)
    _gt_raw   = per_cell_df[gt_col].values.astype(float)
    _gt_valid = np.isfinite(_gt_raw)
    gt_bin    = _gt_raw[_gt_valid].astype(int)   # 0/1 for ROC/PR

    wot_p_all = per_cell_df[wot_col].values.astype(float) if wot_col in per_cell_df.columns \
                else np.full(len(per_cell_df), np.nan)
    wot_p_prim = wot_p_all[_gt_valid]
    fin_wot    = np.isfinite(wot_p_prim)

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle(
        f"WOT vs CellRank2 — OCLR dual-endpoint benchmark  ({TIMESTAMP})\n"
        f"PRIMARY: non-ambiguous hCiPSC (success=1 / failure=0)  |  "
        f"scores: {wot_col} vs {cr2_col or 'N/A'}",
        fontsize=11, fontweight="bold"
    )

    # ── Row 0, col 0: ROC curves ──────────────────────────────────────────────
    ax_roc = axes[0, 0]
    ax_roc.set_title(f"ROC — PRIMARY  ({wot_col})")
    ax_roc.plot([0, 1], [0, 1], "k--", lw=0.8, label="chance")
    try:
        fpr_w, tpr_w, _ = roc_curve(gt_bin[fin_wot], wot_p_prim[fin_wot])
        auroc_w = (float(metrics_df.loc[_wot_idx, "auroc"])
                   if _wot_idx and "auroc" in metrics_df.columns else np.nan)
        ax_roc.plot(fpr_w, tpr_w, lw=2, color="#1f77b4",
                    label=f"WOT  AUROC={auroc_w:.3f}")
    except Exception as e:
        ax_roc.text(0.5, 0.5, f"WOT ROC error:\n{e}", ha="center", va="center",
                    transform=ax_roc.transAxes)
    if has_cr2:
        cr2_p_all  = per_cell_df[cr2_col].values.astype(float)
        cr2_p_prim = cr2_p_all[_gt_valid]
        fin_cr2    = np.isfinite(cr2_p_prim)
        try:
            fpr_c, tpr_c, _ = roc_curve(gt_bin[fin_cr2], cr2_p_prim[fin_cr2])
            auroc_c = (float(metrics_df.loc[_cr2_idx, "auroc"])
                       if _cr2_idx and "auroc" in metrics_df.columns else np.nan)
            ax_roc.plot(fpr_c, tpr_c, lw=2, color="#ff7f0e",
                        label=f"CR2  AUROC={auroc_c:.3f}")
        except Exception as e:
            ax_roc.text(0.5, 0.4, f"CR2 ROC error:\n{e}", ha="center", va="center",
                        transform=ax_roc.transAxes)
    ax_roc.set_xlabel("FPR"); ax_roc.set_ylabel("TPR")
    ax_roc.legend(fontsize=8)

    # ── Row 0, col 1: PR curves ───────────────────────────────────────────────
    ax_pr = axes[0, 1]
    ax_pr.set_title(f"Precision-Recall — PRIMARY  ({wot_col})")
    try:
        prec_w, rec_w, _ = precision_recall_curve(gt_bin[fin_wot], wot_p_prim[fin_wot])
        auprc_w = (float(metrics_df.loc[_wot_idx, "auprc"])
                   if _wot_idx and "auprc" in metrics_df.columns else np.nan)
        ax_pr.plot(rec_w, prec_w, lw=2, color="#1f77b4",
                   label=f"WOT  AUPRC={auprc_w:.3f}")
    except Exception as e:
        ax_pr.text(0.5, 0.5, f"WOT PR error:\n{e}", ha="center", va="center",
                   transform=ax_pr.transAxes)
    if has_cr2:
        try:
            prec_c, rec_c, _ = precision_recall_curve(gt_bin[fin_cr2], cr2_p_prim[fin_cr2])
            auprc_c = (float(metrics_df.loc[_cr2_idx, "auprc"])
                       if _cr2_idx and "auprc" in metrics_df.columns else np.nan)
            ax_pr.plot(rec_c, prec_c, lw=2, color="#ff7f0e",
                       label=f"CR2  AUPRC={auprc_c:.3f}")
        except Exception:
            pass
    _prev = gt_bin.mean() if len(gt_bin) > 0 else np.nan
    if np.isfinite(_prev):
        ax_pr.axhline(_prev, color="k", ls="--", lw=0.8,
                      label=f"chance ({_prev:.3f})")
    ax_pr.set_xlabel("Recall"); ax_pr.set_ylabel("Precision")
    ax_pr.legend(fontsize=8)

    # ── Row 0, col 2: WOT vs CR2 margin scatter ───────────────────────────────
    ax_sc = axes[0, 2]
    _cr2_label = cr2_col or "CR2"
    ax_sc.set_title(f"{wot_col} vs {_cr2_label}\n(non-ambiguous hCiPSC)")
    if has_cr2:
        # Use non-ambiguous primary cells only (gt finite)
        _both_fin = _gt_valid & np.isfinite(wot_p_all) & np.isfinite(cr2_p_all)
        _c = _gt_raw[_both_fin]   # 0.0 or 1.0 — colour by endpoint label
        ax_sc.scatter(wot_p_all[_both_fin], cr2_p_all[_both_fin],
                      c=_c, cmap="RdBu_r", vmin=0, vmax=1,
                      s=3, alpha=0.5, rasterized=True)
        from scipy.stats import spearmanr as _sr
        if _both_fin.sum() > 10:
            _r, _ = _sr(wot_p_all[_both_fin], cr2_p_all[_both_fin])
            ax_sc.text(0.05, 0.95, f"Spearman r={_r:.3f}", transform=ax_sc.transAxes,
                       fontsize=9, va="top")
        ax_sc.set_xlabel(wot_col); ax_sc.set_ylabel(_cr2_label)
    else:
        ax_sc.text(0.5, 0.5, "CellRank2 data\nnot available",
                   ha="center", va="center", transform=ax_sc.transAxes)
        ax_sc.axis("off")

    # ── Row 1, col 0: WOT score boxplot by OCLR endpoint label ──────────────────
    # success vs failure (non-ambiguous primary cells)
    ax_hist_wot = axes[1, 0]
    ax_hist_wot.set_title(f"WOT {wot_col}\nby OCLR endpoint (primary)")
    _wot_suc = wot_p_all[_gt_valid & (np.round(_gt_raw) == 1) & np.isfinite(wot_p_all)]
    _wot_fail = wot_p_all[_gt_valid & (np.round(_gt_raw) == 0) & np.isfinite(wot_p_all)]
    if len(_wot_suc) > 0:
        ax_hist_wot.hist(_wot_suc, bins=40, density=True,
                         alpha=0.7, color="#2a9d8f",
                         label=f"SUCCESS n={len(_wot_suc):,}")
    if len(_wot_fail) > 0:
        ax_hist_wot.hist(_wot_fail, bins=40, density=True,
                         alpha=0.7, color="#e63946",
                         label=f"FAILURE n={len(_wot_fail):,}")
    ax_hist_wot.set_xlabel(wot_col); ax_hist_wot.set_ylabel("Density")
    ax_hist_wot.legend(fontsize=8)

    # ── Row 1, col 1: CR2 score histogram by OCLR endpoint label ─────────────
    ax_hist_cr2 = axes[1, 1]
    ax_hist_cr2.set_title(f"CR2 {_cr2_label}\nby OCLR endpoint (primary)")
    if has_cr2:
        _cr2_suc  = cr2_p_all[_gt_valid & (np.round(_gt_raw) == 1) & np.isfinite(cr2_p_all)]
        _cr2_fail = cr2_p_all[_gt_valid & (np.round(_gt_raw) == 0) & np.isfinite(cr2_p_all)]
        if len(_cr2_suc) > 0:
            ax_hist_cr2.hist(_cr2_suc, bins=40, density=True,
                             alpha=0.7, color="#2a9d8f",
                             label=f"SUCCESS n={len(_cr2_suc):,}")
        if len(_cr2_fail) > 0:
            ax_hist_cr2.hist(_cr2_fail, bins=40, density=True,
                             alpha=0.7, color="#e63946",
                             label=f"FAILURE n={len(_cr2_fail):,}")
        ax_hist_cr2.set_xlabel(_cr2_label); ax_hist_cr2.set_ylabel("Density")
        ax_hist_cr2.legend(fontsize=8)
    else:
        ax_hist_cr2.text(0.5, 0.5, "No CR2 data", ha="center", va="center",
                         transform=ax_hist_cr2.transAxes)
        ax_hist_cr2.axis("off")

    # ── Row 1, col 2: OCLR score distribution ────────────────────────────────
    ax_oclr = axes[1, 2]
    ax_oclr.set_title(f"OCLR score distribution\n({OCLR_ENDPOINT_SAMPLE_ID})")
    _scores = per_cell_df[score_col].dropna()
    if len(_scores) > 0 and _scores.std() > 0:
        ax_oclr.hist(_scores.values, bins=60, color="#2ca02c", alpha=0.8,
                     label=f"n={len(_scores):,}")
        _threshold_files = sorted(PROCESSED_DIR.glob("*_oclr_threshold_summary.csv"),
                                  key=lambda p: p.stat().st_mtime)
        if _threshold_files:
            thr_df = pd.read_csv(_threshold_files[-1])
            _key_col = "threshold_key" if "threshold_key" in thr_df.columns else "key"
            _val_col = "cutoff_value"  if "cutoff_value"  in thr_df.columns else "threshold"
            _final_row = thr_df[thr_df.get(_key_col, pd.Series()).astype(str) == "top_10pct"]
            if len(_final_row) > 0 and _val_col in thr_df.columns:
                _thr = float(_final_row[_val_col].iloc[0])
                ax_oclr.axvline(_thr, color="red", lw=1.5, ls="--",
                                label=f"top-10% threshold={_thr:.3f}")
        ax_oclr.set_xlabel("OCLR score"); ax_oclr.set_ylabel("Count")
        ax_oclr.legend(fontsize=8)
    else:
        ax_oclr.text(0.5, 0.5, "OCLR score\nnot available",
                     ha="center", va="center", transform=ax_oclr.transAxes)
        ax_oclr.axis("off")

    plt.tight_layout()
    out_fig = FIGURES_DIR / f"{TIMESTAMP}_oclr_comparison_scatter.png"
    fig.savefig(str(out_fig), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Comparison figure → {out_fig.name}")


# =============================================================================
# Entry point
# =============================================================================
if __name__ == "__main__":
    main()
