#!/usr/bin/env python3
"""
05_compare_wot_cr2_oclr.py  [v2 — multi-fate design, hpc_full only]
==========================================================================================
Project : Comparative Study of Trajectory Inference Models for Chemical iPSC
          Reprogramming
Dataset : GSE230659 (human, Liuyang et al. 2023 Cell Stem Cell)
Step    : 5 — Compare WOT and CellRank2 fate predictions against OCLR ground truth

Purpose
-------
This standalone post-processing script evaluates the predictive accuracy of:
  (A) WOT          : p_hcipsc          — fate probability toward all hCiPSC cells
                                          (multi-fate design, script 02c v13+)
  (B) CellRank2    : cr2_ips_fate_prob  — fate probability toward the OCLR-annotated
                                          iPSC terminal lineage (script 04 v11+)

against the OCLR-defined high-confidence iPSC terminal set (script 02b2):
  Ground truth = top 10% of OCLR stemness scores within GSM7230012_hCiPSCs-0618

Independence guarantee
----------------------
The OCLR endpoint is computed in 02b2 BEFORE WOT (02c) and CellRank2 (04) are run.
The endpoint definition does not use any WOT or CellRank2 outputs.  Therefore this
evaluation is circularity-free.

NOTE: This script requires the hpc_full CellRank2 result.  The local_debug subset
(from 04_cellrank2.py in local_debug mode) covers only a stratified sample and cannot
be used for final biological conclusions.  A FileNotFoundError is raised if the
hpc_full h5ad is absent — there is no silent local_debug fallback.

Inputs  (auto-detected as most-recent file matching pattern)
-------
  results/metrics/*_GSE230659_wot_fates.h5ad         — WOT fate probs (full dataset)
  results/metrics/*_hpcfull_GSE230659_cellrank2.h5ad  — CellRank2 fate probs (hpc_full)
  data/processed/*_oclr_endpoint_barcodes.txt         — OCLR terminal barcode list

Evaluation sets
---------------
  PRIMARY    : shared barcodes (WOT ∩ CellRank2) + non-hCiPSC cells
               Rationale: (a) identical cell selection removes dataset bias; (b)
               excluding hCiPSC cells avoids trivially high fate probs at the
               terminal stage, giving a more discriminative signal for comparing
               trajectory quality.
  SUPPLEMENTARY : all cells (WOT h5ad) and all CellRank2 cells separately

Metric types
------------
  binary_OCLR_terminal  : AUROC, AUPRC, precision@k, recall@k
                          ground truth = oclr_ips_subset_mask (0/1)
  continuous_OCLR_score : Spearman r vs oclr_score (continuous)
                          Among cells where oclr_score is not NaN

Outputs
-------
  results/cellrank2/{TIMESTAMP}_oclr_comparison_summary.tsv
      — AUROC, AUPRC, Spearman r vs binary label and continuous OCLR score
        (metric_type column distinguishes binary vs continuous)
  results/cellrank2/{TIMESTAMP}_oclr_comparison_per_cell.csv.gz
      — Per-cell: p_hcipsc, cr2_ips_fate_prob, oclr_ground_truth, oclr_score
  results/cellrank2/{TIMESTAMP}_oclr_comparison_wot_vs_cr2.tsv
      — WOT vs CellRank2 correlation on shared cells
  results/figures/{TIMESTAMP}_oclr_comparison_scatter.png

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
    print("05_compare_wot_cr2_oclr.py  v2  —  WOT vs CellRank2 under OCLR endpoint")
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
    _wot_req  = ["p_hcipsc"]
    _cr2_req  = ["cr2_ips_fate_prob"]
    _wot_miss = [c for c in _wot_req if c not in adata_wot.obs.columns]
    _cr2_miss = [c for c in _cr2_req if c not in adata_cr2.obs.columns]
    if _wot_miss:
        raise ValueError(
            f"WOT h5ad missing required columns: {_wot_miss}\n"
            "Ensure 02c_compute_transport.py v13+ (multi-fate design) has run."
        )
    if _cr2_miss:
        raise ValueError(
            f"CR2 h5ad missing required columns: {_cr2_miss}\n"
            "Ensure 04_cellrank2.py v11+ (standard GPCCA + OCLR annotation) has run."
        )

    # ── 2.4  Ground truth on WOT (full dataset) ───────────────────────────────
    print("\n[3] Building OCLR ground truth ...")
    gt_wot, oclr_score_wot = load_oclr_ground_truth(adata_wot.obs)
    n_terminal = int(gt_wot.sum())
    n_total    = len(gt_wot)
    print(f"  OCLR terminal cells (WOT): {n_terminal:,} / {n_total:,} "
          f"({n_terminal/n_total*100:.2f}%)")

    gt_cr2, oclr_score_cr2 = load_oclr_ground_truth(adata_cr2.obs)
    n_terminal_cr2 = int(gt_cr2.sum())
    print(f"  OCLR terminal cells (CR2): {n_terminal_cr2:,} / {len(gt_cr2):,} "
          f"({n_terminal_cr2/max(len(gt_cr2),1)*100:.2f}%)")

    # ── 2.5  Identify shared barcodes and non-hCiPSC mask ────────────────────
    shared_bc = adata_wot.obs_names.intersection(adata_cr2.obs_names)
    n_shared  = len(shared_bc)
    print(f"\n[4] Shared barcodes: {n_shared:,}")

    if n_shared < 50:
        raise ValueError(
            f"Only {n_shared} shared barcodes between WOT and CellRank2 h5ads. "
            "Cannot compute the PRIMARY (shared-cells) evaluation. "
            "Ensure both h5ads come from the full dataset run."
        )

    # Identify hCiPSC cells to exclude from the PRIMARY evaluation.
    # Prefer stage_std column (set by 02c v13); fall back to stage column.
    def _hcipsc_mask(obs: pd.DataFrame, bc_index, label: str) -> np.ndarray:
        for _col in ["stage_std", "stage"]:
            if _col in obs.columns:
                _m = obs.loc[bc_index, _col].astype(str).str.contains(
                    "hCiPSC|iPSC|iPsc|ipsc|iPS", case=False, na=False
                )
                n_hcipsc = int(_m.sum())
                print(f"  hCiPSC mask ({label}, col='{_col}'): {n_hcipsc:,} / {len(bc_index):,} cells")
                return _m.values.astype(bool)
        print(f"  [WARN] No stage column found for hCiPSC mask ({label}) — "
              f"PRIMARY evaluation will use ALL shared cells.")
        return np.zeros(len(bc_index), dtype=bool)

    _hcipsc_mask_wot_shared = _hcipsc_mask(adata_wot.obs, shared_bc, "WOT-shared")
    _non_hcipsc_shared      = ~_hcipsc_mask_wot_shared
    primary_bc              = shared_bc[_non_hcipsc_shared]
    n_primary               = len(primary_bc)
    print(f"  PRIMARY eval set: {n_primary:,} shared non-hCiPSC cells")

    if n_primary < 50:
        print(f"  [WARN] PRIMARY set has only {n_primary} cells — "
              f"falling back to all shared cells for PRIMARY evaluation.")
        primary_bc = shared_bc
        n_primary  = len(primary_bc)

    # ── 2.6  Compute evaluation metrics ──────────────────────────────────────
    print("\n[5] Computing evaluation metrics ...")

    records = []

    # ── PRIMARY: shared + non-hCiPSC ──────────────────────────────────────────
    _gt_primary        = gt_wot.reindex(primary_bc).values.astype(int)
    _wot_primary       = adata_wot.obs.loc[primary_bc, "p_hcipsc"].values.astype(float)
    _cr2_primary       = adata_cr2.obs.loc[primary_bc, "cr2_ips_fate_prob"].values.astype(float)
    _oclr_sc_primary   = oclr_score_wot.reindex(primary_bc).values.astype(float)

    print(f"\n  -- PRIMARY (shared, non-hCiPSC, n={n_primary:,}) --")
    _r_wot = compute_metrics(_gt_primary, _wot_primary,
                             "WOT_primary_shared_nonterminal")
    _r_wot["eval_set"] = "primary_shared_nonterminal"
    records.append(_r_wot)

    _r_cr2 = compute_metrics(_gt_primary, _cr2_primary,
                             "CR2_primary_shared_nonterminal")
    _r_cr2["eval_set"] = "primary_shared_nonterminal"
    records.append(_r_cr2)

    # Continuous OCLR score correlation — PRIMARY
    _rc_wot = compute_oclr_score_correlation(
        _wot_primary, _oclr_sc_primary, "WOT_primary_shared_nonterminal")
    _rc_wot["eval_set"] = "primary_shared_nonterminal"
    records.append(_rc_wot)

    _rc_cr2 = compute_oclr_score_correlation(
        _cr2_primary, _oclr_sc_primary, "CR2_primary_shared_nonterminal")
    _rc_cr2["eval_set"] = "primary_shared_nonterminal"
    records.append(_rc_cr2)

    # ── SUPPLEMENTARY: all WOT cells ──────────────────────────────────────────
    print(f"\n  -- SUPPLEMENTARY (all WOT cells, n={adata_wot.n_obs:,}) --")
    _wot_all = adata_wot.obs["p_hcipsc"].values.astype(float)
    _oclr_sc_wot_all = oclr_score_wot.values.astype(float)

    _r_wot_all = compute_metrics(gt_wot.values.astype(int), _wot_all, "WOT_all")
    _r_wot_all["eval_set"] = "supplementary_all_wot"
    records.append(_r_wot_all)

    _rc_wot_all = compute_oclr_score_correlation(
        _wot_all, _oclr_sc_wot_all, "WOT_all")
    _rc_wot_all["eval_set"] = "supplementary_all_wot"
    records.append(_rc_wot_all)

    # ── SUPPLEMENTARY: all CR2 cells ──────────────────────────────────────────
    print(f"\n  -- SUPPLEMENTARY (all CR2 cells, n={adata_cr2.n_obs:,}) --")
    _cr2_all     = adata_cr2.obs["cr2_ips_fate_prob"].values.astype(float)
    _oclr_sc_cr2 = oclr_score_cr2.values.astype(float)

    _r_cr2_all = compute_metrics(gt_cr2.values.astype(int), _cr2_all, "CR2_all")
    _r_cr2_all["eval_set"] = "supplementary_all_cr2"
    records.append(_r_cr2_all)

    _rc_cr2_all = compute_oclr_score_correlation(
        _cr2_all, _oclr_sc_cr2, "CR2_all")
    _rc_cr2_all["eval_set"] = "supplementary_all_cr2"
    records.append(_rc_cr2_all)

    # ── SUPPLEMENTARY: shared cells (all, including hCiPSC) ───────────────────
    print(f"\n  -- SUPPLEMENTARY (all shared cells, n={n_shared:,}) --")
    _gt_shared       = gt_wot.reindex(shared_bc).values.astype(int)
    _wot_shared      = adata_wot.obs.loc[shared_bc, "p_hcipsc"].values.astype(float)
    _cr2_shared      = adata_cr2.obs.loc[shared_bc, "cr2_ips_fate_prob"].values.astype(float)
    _oclr_sc_shared  = oclr_score_wot.reindex(shared_bc).values.astype(float)

    _r_wot_sh = compute_metrics(_gt_shared, _wot_shared, "WOT_supplementary_shared_all")
    _r_wot_sh["eval_set"] = "supplementary_shared_all"
    records.append(_r_wot_sh)

    _r_cr2_sh = compute_metrics(_gt_shared, _cr2_shared, "CR2_supplementary_shared_all")
    _r_cr2_sh["eval_set"] = "supplementary_shared_all"
    records.append(_r_cr2_sh)

    _rc_wot_sh = compute_oclr_score_correlation(
        _wot_shared, _oclr_sc_shared, "WOT_supplementary_shared_all")
    _rc_wot_sh["eval_set"] = "supplementary_shared_all"
    records.append(_rc_wot_sh)

    _rc_cr2_sh = compute_oclr_score_correlation(
        _cr2_shared, _oclr_sc_shared, "CR2_supplementary_shared_all")
    _rc_cr2_sh["eval_set"] = "supplementary_shared_all"
    records.append(_rc_cr2_sh)

    # Print summary
    metrics_df = pd.DataFrame(records)
    _print_cols = ["model", "metric_type", "eval_set", "n_cells",
                   "auroc", "auprc", "spearman_r", "pearson_r"]
    _print_cols = [c for c in _print_cols if c in metrics_df.columns]
    print("\n  Metrics summary vs OCLR ground truth:")
    print(metrics_df[_print_cols].to_string(index=False))

    # ── 2.7  WOT vs CellRank2 agreement (shared cells) ───────────────────────
    print("\n[6] WOT vs CellRank2 agreement (shared cells) ...")
    from scipy.stats import pearsonr, spearmanr

    wot_vs_cr2_rows = []
    _fin_shared = np.isfinite(_wot_shared) & np.isfinite(_cr2_shared)
    if _fin_shared.sum() < 50:
        print("  [WARN] Fewer than 50 finite shared-cell pairs — correlation skipped")
    else:
        pr_all, pp_all = pearsonr(_wot_shared[_fin_shared], _cr2_shared[_fin_shared])
        sr_all, sp_all = spearmanr(_wot_shared[_fin_shared], _cr2_shared[_fin_shared])
        wot_vs_cr2_rows += [
            {"metric": "pearson_r_all",  "value": pr_all},
            {"metric": "pearson_p_all",  "value": pp_all},
            {"metric": "spearman_r_all", "value": sr_all},
            {"metric": "spearman_p_all", "value": sp_all},
            {"metric": "n_cells_all",    "value": int(_fin_shared.sum())},
        ]
        print(f"  Pearson  r (all shared) = {pr_all:.4f}  (p={pp_all:.2e})")
        print(f"  Spearman r (all shared) = {sr_all:.4f}  (p={sp_all:.2e})")

        # Non-hCiPSC only (PRIMARY set)
        _fin_primary = np.isfinite(_wot_primary) & np.isfinite(_cr2_primary)
        if _fin_primary.sum() > 50:
            pr_nt, pp_nt = pearsonr(_wot_primary[_fin_primary], _cr2_primary[_fin_primary])
            sr_nt, sp_nt = spearmanr(_wot_primary[_fin_primary], _cr2_primary[_fin_primary])
            wot_vs_cr2_rows += [
                {"metric": "pearson_r_nonterminal",  "value": pr_nt},
                {"metric": "pearson_p_nonterminal",  "value": pp_nt},
                {"metric": "spearman_r_nonterminal", "value": sr_nt},
                {"metric": "spearman_p_nonterminal", "value": sp_nt},
                {"metric": "n_cells_nonterminal",    "value": int(_fin_primary.sum())},
            ]
            print(f"  Pearson  r (non-hCiPSC) = {pr_nt:.4f}  (p={pp_nt:.2e})")
            print(f"  Spearman r (non-hCiPSC) = {sr_nt:.4f}  (p={sp_nt:.2e})")

    # ── 2.8  Per-cell predictions table ──────────────────────────────────────
    print("\n[7] Building per-cell predictions table ...")
    per_cell_df = pd.DataFrame({
        "p_hcipsc_wot":       adata_wot.obs["p_hcipsc"],
        "oclr_ground_truth":  gt_wot,
        "oclr_score":         oclr_score_wot,
    }, index=adata_wot.obs_names)

    if n_shared > 0:
        cr2_col_series = pd.Series(
            adata_cr2.obs["cr2_ips_fate_prob"].values,
            index=adata_cr2.obs_names,
            name="cr2_ips_fate_prob",
        )
        per_cell_df = per_cell_df.join(cr2_col_series, how="left")
        # Also join cr2_ips_lineage_name if present
        if "cr2_ips_lineage_name" in adata_cr2.obs.columns:
            _lin_series = pd.Series(
                adata_cr2.obs["cr2_ips_lineage_name"].values,
                index=adata_cr2.obs_names,
                name="cr2_ips_lineage_name",
            )
            per_cell_df = per_cell_df.join(_lin_series, how="left")

    for _meta in ["sample_id", "stage", "stage_std", "abs_day", "stage_day_label"]:
        if _meta in adata_wot.obs.columns:
            per_cell_df[_meta] = adata_wot.obs[_meta].values

    print(f"  Per-cell table shape: {per_cell_df.shape}")

    # ── 2.9  Save outputs ─────────────────────────────────────────────────────
    print("\n[8] Saving outputs ...")

    # a) Metrics summary
    out_summary = CR2_RESULTS / f"{TIMESTAMP}_oclr_comparison_summary.tsv"
    _summary_df = metrics_df.copy()
    _summary_df.insert(0, "wot_n_cells_total",  adata_wot.shape[0])
    _summary_df.insert(1, "cr2_n_cells_total",  adata_cr2.shape[0])
    _summary_df.insert(2, "n_shared_cells",      n_shared)
    _summary_df.insert(3, "n_primary_cells",     n_primary)
    _summary_df.to_csv(out_summary, sep="\t", index=False)
    print(f"  Metrics summary  → {out_summary.name}")

    # b) WOT vs CR2 agreement table
    out_corr = CR2_RESULTS / f"{TIMESTAMP}_oclr_comparison_wot_vs_cr2.tsv"
    if wot_vs_cr2_rows:
        _corr_df = pd.DataFrame(wot_vs_cr2_rows)
        _corr_df.insert(0, "wot_n_cells",  adata_wot.shape[0])
        _corr_df.insert(1, "cr2_n_cells",  adata_cr2.shape[0])
        _corr_df.to_csv(out_corr, sep="\t", index=False)
        print(f"  WOT vs CR2 corr  → {out_corr.name}")

    # c) Per-cell predictions
    out_percell = CR2_RESULTS / f"{TIMESTAMP}_oclr_comparison_per_cell.csv.gz"
    per_cell_df.to_csv(out_percell, compression="gzip")
    print(f"  Per-cell table   → {out_percell.name}")

    # ── 2.10  Diagnostic figure ───────────────────────────────────────────────
    print("\n[9] Generating comparison figure ...")

    # Build a metrics_df subset indexed by model name for the figure function
    _metrics_indexed = metrics_df.set_index("model")

    _make_comparison_figure(
        per_cell_df        = per_cell_df,
        metrics_df         = _metrics_indexed,
        gt_col             = "oclr_ground_truth",
        wot_col            = "p_hcipsc_wot",
        cr2_col            = "cr2_ips_fate_prob" if "cr2_ips_fate_prob" in per_cell_df.columns else None,
        score_col          = "oclr_score",
        cr2_metrics_label  = "CR2_primary_shared_nonterminal",
    )

    print("\n" + "=" * 70)
    print("DONE")
    print(f"  wot_n_cells={adata_wot.n_obs:,}  cr2_n_cells={adata_cr2.n_obs:,}")
    print(f"  n_shared={n_shared:,}  n_primary={n_primary:,}")
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
    cr2_metrics_label: str = "",
) -> None:
    """Generate a 3×2 comparison figure and save to FIGURES_DIR.

    cr2_metrics_label : the metrics_df index row to use for CR2 AUROC/AUPRC
                        (e.g. 'CR2_primary_shared_nonterminal').
                        If empty or not found, the first 'CR2_.*_all' row is used.
    """
    if cr2_metrics_label and cr2_metrics_label in metrics_df.index:
        _cr2_idx = cr2_metrics_label
    else:
        _cr2_idx = next(
            (idx for idx in metrics_df.index
             if idx.startswith("CR2_") and idx.endswith("_all")),
            None
        )
    from sklearn.metrics import precision_recall_curve, roc_curve

    has_cr2 = (cr2_col is not None) and (cr2_col in per_cell_df.columns)
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle(
        f"WOT vs CellRank2 — OCLR endpoint evaluation\n{TIMESTAMP}",
        fontsize=13, fontweight="bold"
    )

    gt      = per_cell_df[gt_col].values.astype(int)
    wot_p   = per_cell_df[wot_col].values.astype(float)
    fin_wot = np.isfinite(wot_p)

    # ── Row 0: ROC curves ─────────────────────────────────────────────────────
    ax_roc = axes[0, 0]
    ax_roc.set_title("ROC — WOT p_hcipsc vs OCLR binary label")
    ax_roc.plot([0, 1], [0, 1], "k--", lw=0.8, label="chance")
    try:
        fpr_w, tpr_w, _ = roc_curve(gt[fin_wot], wot_p[fin_wot])
        _wot_idx = next((i for i in metrics_df.index
                         if i.startswith("WOT_") and "all" in i), None)
        auroc_w = float(metrics_df.loc[_wot_idx, "auroc"]) if _wot_idx else np.nan
        ax_roc.plot(fpr_w, tpr_w, lw=2, color="#1f77b4",
                    label=f"WOT  AUROC={auroc_w:.3f}")
    except Exception as e:
        ax_roc.text(0.5, 0.5, f"ROC error:\n{e}", ha="center", va="center",
                    transform=ax_roc.transAxes)
    if has_cr2:
        cr2_p   = per_cell_df[cr2_col].values.astype(float)
        gt_cr2  = gt.copy()
        fin_cr2 = np.isfinite(cr2_p)
        try:
            fpr_c, tpr_c, _ = roc_curve(gt_cr2[fin_cr2], cr2_p[fin_cr2])
            auroc_c = (float(metrics_df.loc[_cr2_idx, "auroc"])
                       if _cr2_idx is not None and "auroc" in metrics_df.columns
                       else np.nan)
            _cr2_short = _cr2_idx or "CR2"
            ax_roc.plot(fpr_c, tpr_c, lw=2, color="#ff7f0e",
                        label=f"{_cr2_short}  AUROC={auroc_c:.3f}")
        except Exception as e:
            ax_roc.text(0.5, 0.4, f"CR2 ROC error:\n{e}", ha="center", va="center",
                        transform=ax_roc.transAxes)
    ax_roc.set_xlabel("FPR"); ax_roc.set_ylabel("TPR")
    ax_roc.legend(fontsize=8)

    # ── Row 0, col 1: PR curves ───────────────────────────────────────────────
    ax_pr = axes[0, 1]
    ax_pr.set_title("Precision-Recall — WOT vs CellRank2")
    try:
        prec_w, rec_w, _ = precision_recall_curve(gt[fin_wot], wot_p[fin_wot])
        _wot_idx_pr = next((i for i in metrics_df.index
                            if i.startswith("WOT_") and "all" in i), None)
        auprc_w = float(metrics_df.loc[_wot_idx_pr, "auprc"]) if _wot_idx_pr else np.nan
        ax_pr.plot(rec_w, prec_w, lw=2, color="#1f77b4",
                   label=f"WOT  AUPRC={auprc_w:.3f}")
    except Exception as e:
        ax_pr.text(0.5, 0.5, f"PR error:\n{e}", ha="center", va="center",
                   transform=ax_pr.transAxes)
    if has_cr2:
        try:
            prec_c, rec_c, _ = precision_recall_curve(gt_cr2[fin_cr2], cr2_p[fin_cr2])
            auprc_c = (float(metrics_df.loc[_cr2_idx, "auprc"])
                       if _cr2_idx is not None and "auprc" in metrics_df.columns
                       else np.nan)
            ax_pr.plot(rec_c, prec_c, lw=2, color="#ff7f0e",
                       label=f"{_cr2_short}  AUPRC={auprc_c:.3f}")
        except Exception:
            pass
    _prev = gt[fin_wot].mean() if fin_wot.sum() > 0 else np.nan
    if np.isfinite(_prev):
        ax_pr.axhline(_prev, color="k", ls="--", lw=0.8, label=f"chance ({_prev:.3f})")
    ax_pr.set_xlabel("Recall"); ax_pr.set_ylabel("Precision")
    ax_pr.legend(fontsize=8)

    # ── Row 0, col 2: WOT vs CR2 scatter ─────────────────────────────────────
    ax_sc = axes[0, 2]
    ax_sc.set_title("WOT p_hcipsc vs CR2 cr2_ips_fate_prob\n(shared cells)")
    if has_cr2:
        shared_mask = fin_wot & np.isfinite(cr2_p)
        _c = gt[shared_mask].astype(float)
        ax_sc.scatter(wot_p[shared_mask], cr2_p[shared_mask],
                      c=_c, cmap="RdBu_r", s=2, alpha=0.4, rasterized=True)
        from scipy.stats import spearmanr as _sr
        _r, _ = _sr(wot_p[shared_mask], cr2_p[shared_mask])
        ax_sc.set_xlabel("WOT p_hcipsc")
        ax_sc.set_ylabel("CR2 cr2_ips_fate_prob")
        ax_sc.text(0.05, 0.95, f"Spearman r={_r:.3f}", transform=ax_sc.transAxes,
                   fontsize=9, va="top")
    else:
        ax_sc.text(0.5, 0.5, "CellRank2 data\nnot available",
                   ha="center", va="center", transform=ax_sc.transAxes)
        ax_sc.axis("off")

    # ── Row 1: WOT histogram by OCLR label ────────────────────────────────────
    ax_hist_wot = axes[1, 0]
    ax_hist_wot.set_title("WOT p_hcipsc distribution by OCLR label")
    _pos_mask = (gt == 1) & fin_wot
    _neg_mask = (gt == 0) & fin_wot
    if _pos_mask.sum() > 0:
        ax_hist_wot.hist(wot_p[_pos_mask], bins=50, density=True,
                         alpha=0.7, color="#d62728", label=f"OCLR+ n={_pos_mask.sum():,}")
    if _neg_mask.sum() > 0:
        ax_hist_wot.hist(wot_p[_neg_mask], bins=50, density=True,
                         alpha=0.5, color="#aec7e8", label=f"OCLR− n={_neg_mask.sum():,}")
    ax_hist_wot.set_xlabel("p_hcipsc"); ax_hist_wot.set_ylabel("Density")
    ax_hist_wot.legend(fontsize=8)

    # ── Row 1, col 1: CR2 histogram by OCLR label ────────────────────────────
    ax_hist_cr2 = axes[1, 1]
    ax_hist_cr2.set_title("CR2 cr2_ips_fate_prob distribution by OCLR label")
    if has_cr2:
        _pos_cr2 = (gt_cr2 == 1) & fin_cr2
        _neg_cr2 = (gt_cr2 == 0) & fin_cr2
        if _pos_cr2.sum() > 0:
            ax_hist_cr2.hist(cr2_p[_pos_cr2], bins=50, density=True,
                             alpha=0.7, color="#d62728",
                             label=f"OCLR+ n={_pos_cr2.sum():,}")
        if _neg_cr2.sum() > 0:
            ax_hist_cr2.hist(cr2_p[_neg_cr2], bins=50, density=True,
                             alpha=0.5, color="#ffbb78",
                             label=f"OCLR− n={_neg_cr2.sum():,}")
        ax_hist_cr2.set_xlabel("cr2_ips_fate_prob"); ax_hist_cr2.set_ylabel("Density")
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
