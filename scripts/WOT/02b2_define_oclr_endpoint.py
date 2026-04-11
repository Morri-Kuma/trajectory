#!/usr/bin/env python3
"""
02b2_define_oclr_endpoint.py  —  OCLR-based iPS terminal-cell endpoint definition.

Project : Comparative Study of Trajectory Inference Models for Chemical iPSC Reprogramming
Dataset : GSE230659 (human, Liuyang et al. 2023 Cell Stem Cell)
Step    : 2b-2 — Load / generate OCLR stemness scores for GSM7230012_hCiPSCs-0618,
          threshold to a high-confidence iPS terminal set, and save the barcode list.

Method reference
----------------
  Mäkinen V-P et al. "A machine learning one-class logistic regression model to predict
  stemness for single cell transcriptomics and spatial omics." (OCLR method)

  Per-cell OCLR stemness score:
    1. Identify overlap genes between the OCLR signature and the cell expression matrix.
    2. Compute the Spearman correlation between each cell's log-normalised expression
       vector (restricted to overlap genes) and the OCLR signature weights.
    3. Linearly rescale the resulting scores to [0, 1].
    4. Higher score = more stem-like / more PSC-like.

  The score is NOT a trajectory or clustering result.  It is used here as an
  independent, model-based criterion to identify high-confidence iPS cells in the
  hCiPSC endpoint sample BEFORE running WOT or CellRank2.

Input
-----
  OCLR score CSV  (primary):
      data/processed/GSE230659hCiPSCs0618_oclr_score.csv
      Expected columns: a cell-ID column (barcode) + an OCLR score column.
      The cell-ID column may be:
          • Full project barcodes: GSM7230012_hCiPSCs-0618_AAACCCAAGTCCCAGC-1
          • Short barcodes:        AAACCCAAGTCCCAGC-1          (prefixed automatically)
          • Bare barcodes:         AAACCCAAGTCCCAGC             (-1 appended automatically)
      Barcode mismatches are reported explicitly; unmatched cells are not silently dropped.

  WOT-ready AnnData (for barcode cross-reference and expression matrix):
      data/processed/*_GSE230659_wot_ready.h5ad   or
      data/processed/*_GSE230659_wot_gr.h5ad

Outputs (all timestamped, in data/processed/)
---------------------------------------------------------------------------
  YYYYMMDD_HHMM_oclr_endpoint_barcodes.txt       — final selected barcode list
  YYYYMMDD_HHMM_oclr_score_all_hcipsc.csv        — OCLR scores for all 9,142 hCiPSC cells
  YYYYMMDD_HHMM_oclr_threshold_summary.csv        — distribution + threshold candidates
  YYYYMMDD_HHMM_oclr_score_distribution.png       — diagnostic histogram + threshold lines

Usage
-----
    conda activate traj_env
    cd C:\\Users\\37620\\trajectory
    python scripts/WOT/02b2_define_oclr_endpoint.py

Prerequisites
-------------
    Step 02a output  : *_GSE230659_wot_ready.h5ad  (or *_wot_gr.h5ad)
    OCLR score CSV   : data/processed/GSE230659hCiPSCs0618_oclr_score.csv
                       (place it there before running this script)

Next step
---------
    Run 02c_compute_transport.py  (reads the barcodes file automatically)
"""

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
from datetime import datetime
from pathlib import Path

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")

# ── Project root (auto-detected) ──────────────────────────────────────────────
_candidates = [
    Path(r"C:\Users\37620\trajectory"),
    Path("/sessions/dazzling-peaceful-brahmagupta/mnt/trajectory"),
    Path("/sessions/cool-admiring-hawking/mnt/trajectory"),
    Path(__file__).resolve().parents[2],
]
PROJECT_ROOT  = next((p for p in _candidates if p.exists()), _candidates[-1])
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
FIGURES_DIR   = PROJECT_ROOT / "results" / "figures" / "oclr_endpoint"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ── OCLR file search locations (in priority order) ────────────────────────────
# The primary location is data/processed/; we also check a few fallbacks.
OCLR_SEARCH_PATHS = [
    PROCESSED_DIR / "GSE230659hCiPSCs0618_oclr_score.csv",
    PROJECT_ROOT  / "data" / "GSE230659hCiPSCs0618_oclr_score.csv",
    PROJECT_ROOT  / "GSE230659hCiPSCs0618_oclr_score.csv",
    # Add more candidate paths here if needed
]

# ── Sample identity ───────────────────────────────────────────────────────────
ENDPOINT_SAMPLE_ID: str = "GSM7230012_hCiPSCs-0618"
# Prefix used in the full-project barcode format
BARCODE_PREFIX: str = "GSM7230012_hCiPSCs-0618_"

# ── Thresholding candidates to evaluate ──────────────────────────────────────
# Each entry: (label, quantile_of_hcipsc_oclr_scores)
# We will evaluate all of them and choose the final threshold explicitly.
THRESHOLD_CANDIDATES = [
    ("top_5pct",  0.95),   # top 5% — high specificity
    ("top_10pct", 0.90),   # top 10% — primary candidate
    ("top_15pct", 0.85),   # top 15%
    ("top_20pct", 0.80),   # top 20%
]

# ── Final threshold choice ────────────────────────────────────────────────────
# PROJECT-SPECIFIC CUTOFF — NOT a threshold defined by the OCLR paper.
#
# The original OCLR method (Mäkinen V-P et al.) provides per-cell stemness
# scores but does NOT define a universal binary threshold for "terminal iPSC".
# The top-10% cutoff below is a project decision made here to create a
# high-confidence evaluation set for trajectory benchmarking.
#
# Rationale: top-10% balances specificity (high-pluripotency requirement)
# against sufficient cell count (≥ MIN_TERMINAL_CELLS).  This threshold is
# independent of WOT and CellRank2 outputs — no circularity.  If the score
# distribution shows a clear bimodal upper tail, switch to "top_5pct" and
# re-run.
FINAL_THRESHOLD_KEY: str = "top_10pct"

# Minimum cells in the final terminal set (safety check)
MIN_TERMINAL_CELLS: int = 50

# ── OCLR score column name(s) to try ─────────────────────────────────────────
# We try these column names in the CSV, in order.
OCLR_SCORE_COLS = ["oclr_score", "score", "stemness_score", "OCLR_score",
                   "oclr", "stemness", "Score"]

# NOTE: No proxy / fallback computation is provided.
# The real OCLR CSV (GSE230659hCiPSCs0618_oclr_score.csv) is required.
# If the CSV is missing the script stops with a clear error message.
# Place the CSV at data/processed/ before running this script.

# =============================================================================
# 1. IMPORTS
# =============================================================================
import gc
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

print(f"[{TIMESTAMP}]  OCLR Endpoint Definition  —  Step 02b2")
print(f"  Project root     : {PROJECT_ROOT}")
print(f"  Endpoint sample  : {ENDPOINT_SAMPLE_ID}")
print(f"  OCLR CSV search  :")
for _p in OCLR_SEARCH_PATHS:
    _status = "✓ FOUND" if _p.exists() else "  missing"
    print(f"    {_status}  {_p}")

# =============================================================================
# 2. LOCATE / LOAD AnnData (for barcodes + fallback expression)
# =============================================================================
def find_latest(pattern: str) -> Path:
    matches = sorted(PROCESSED_DIR.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No files matching '{pattern}' in {PROCESSED_DIR}."
        )
    return matches[-1]

print("\n--- Loading AnnData ---")
# Prefer wot_gr (has growth rates); fall back to wot_ready
try:
    H5AD_PATH = find_latest("*_GSE230659_wot_gr.h5ad")
except FileNotFoundError:
    H5AD_PATH = find_latest("*_GSE230659_wot_ready.h5ad")

print(f"  Input AnnData : {H5AD_PATH.name}")
adata_full = ad.read_h5ad(str(H5AD_PATH))
print(f"  Shape         : {adata_full.n_obs:,} × {adata_full.n_vars:,}")

# Identify hCiPSC cells for OCLR thresholding.
#
# Priority:
#   1. Exact sample_id match: sample_id == ENDPOINT_SAMPLE_ID
#      → guarantees only the one intended endpoint sample is included
#   2. Fallback: stage column contains 'cipsc' / 'ipsc'
#      → used only when sample_id is unavailable; printed explicitly
#
# IMPORTANT: do NOT silently merge multiple iPSC-like samples into the
# thresholding pool.  If sample_id is present, always use exact matching.

_POOL_METHOD: str = ""

if "sample_id" in adata_full.obs.columns:
    _hcipsc_mask = adata_full.obs["sample_id"].astype(str) == ENDPOINT_SAMPLE_ID
    _n_by_sample_id = int(_hcipsc_mask.sum())
    if _n_by_sample_id > 0:
        _POOL_METHOD = f"exact sample_id == '{ENDPOINT_SAMPLE_ID}'"
        print(f"  Pool selection  : {_POOL_METHOD}")
        print(f"  hCiPSC cells    : {_n_by_sample_id:,}")
    else:
        # sample_id column exists but no cells match — warn and fall through
        print(f"  [WARN] sample_id column found but 0 cells matched "
              f"'{ENDPOINT_SAMPLE_ID}'.  Unique sample_id values present:")
        for _sid in adata_full.obs["sample_id"].unique()[:10]:
            print(f"    {_sid!r}")
        _hcipsc_mask = None
else:
    _hcipsc_mask = None

if _hcipsc_mask is None or int(_hcipsc_mask.sum()) == 0:
    # Fallback: stage column string match
    _stage_col = "stage"
    assert _stage_col in adata_full.obs.columns, (
        f"Neither 'sample_id' nor '{_stage_col}' found in adata.obs. "
        "Run 02a first and ensure obs metadata is present."
    )
    _hcipsc_mask = (
        adata_full.obs[_stage_col].astype(str)
        .str.lower()
        .str.contains("cipsc|hcipsc|ipsc", na=False)
    )
    _POOL_METHOD = "stage fallback (contains 'cipsc'/'ipsc')"
    print(f"  Pool selection  : {_POOL_METHOD}  [FALLBACK — sample_id not available]")
    print(f"  [NOTE] This fallback may include cells from multiple iPSC-like samples.")
    print(f"         Provide a dataset with a 'sample_id' column for exact selection.")

_hcipsc_barcodes = adata_full.obs_names[_hcipsc_mask].tolist()
n_hcipsc = len(_hcipsc_barcodes)
print(f"  hCiPSC cells selected : {n_hcipsc:,}  (pool method: {_POOL_METHOD})")
assert n_hcipsc > 0, (
    f"No hCiPSC cells found using pool method '{_POOL_METHOD}'. "
    "Check obs metadata (sample_id / stage columns)."
)

# =============================================================================
# 3. LOAD OCLR SCORE CSV  (real CSV only — no fallback)
# =============================================================================
_oclr_found  = None
for _candidate in OCLR_SEARCH_PATHS:
    if _candidate.exists():
        _oclr_found = _candidate
        break

if _oclr_found is None:
    raise FileNotFoundError(
        "OCLR score CSV not found at any of the following locations:\n"
        + "\n".join(f"  {p}" for p in OCLR_SEARCH_PATHS)
        + "\n\nPlace the real OCLR CSV (GSE230659hCiPSCs0618_oclr_score.csv) at:\n"
        f"  {OCLR_SEARCH_PATHS[0]}\nand re-run this script.\n"
        "The proxy / fallback computation has been removed — only the real\n"
        "published OCLR scores (Mäkinen V-P et al.) are accepted."
    )

# ─── Load the provided OCLR CSV ─────────────────────────────────────────────
print(f"\n--- Loading OCLR score CSV ---")
print(f"  File : {_oclr_found}")
_raw = pd.read_csv(_oclr_found)
print(f"  Raw shape  : {_raw.shape}")
print(f"  Raw columns: {list(_raw.columns)}")
print(f"  Preview:\n{_raw.head(3).to_string()}")

# Identify the score column
_score_col = None
for _cname in OCLR_SCORE_COLS:
    if _cname in _raw.columns:
        _score_col = _cname
        break
if _score_col is None:
    # Heuristic: the score column is the last numeric column
    _num_cols = _raw.select_dtypes(include=[np.number]).columns.tolist()
    if not _num_cols:
        raise ValueError(
            f"Cannot find an OCLR score column in {_oclr_found}.\n"
            f"Tried: {OCLR_SCORE_COLS}\nColumns present: {list(_raw.columns)}"
        )
    _score_col = _num_cols[-1]
    print(f"  [WARN] Score column not in standard list; using '{_score_col}'")
else:
    print(f"  Score column : '{_score_col}'")

# Identify the barcode / cell-ID column (first non-numeric column, or index)
_id_col = None
if _raw.index.name is not None and _raw.index.name != "":
    # scores indexed by barcode
    _raw_barcodes = _raw.index.astype(str).tolist()
else:
    _str_cols = _raw.select_dtypes(exclude=[np.number]).columns.tolist()
    if _str_cols:
        _id_col      = _str_cols[0]
        _raw_barcodes = _raw[_id_col].astype(str).tolist()
    else:
        # All columns numeric → use row 0 index
        _raw_barcodes = _raw.index.astype(str).tolist()
_raw_scores = _raw[_score_col].astype(float).values

print(f"  Barcode column source: "
      f"{'index' if _id_col is None else repr(_id_col)}")
print(f"  N barcodes in CSV    : {len(_raw_barcodes):,}")

# ── Barcode alignment ─────────────────────────────────────────────────────
# Three possible CSV formats:
#   1. Full project barcodes: GSM7230012_hCiPSCs-0618_AAACCCAAGTCCCAGC-1  (exact match)
#   2. Short barcodes with -1: AAACCCAAGTCCCAGC-1  (needs prefix)
#   3. Bare barcodes without -1: AAACCCAAGTCCCAGC   (needs prefix + suffix)
#
# We attempt exact match first, then progressive correction.

_proj_set = set(_hcipsc_barcodes)   # full project barcodes for GSM7230012

# Attempt 1: exact match
_exact_match = sum(1 for b in _raw_barcodes if b in _proj_set)
if _exact_match > 0:
    print(f"\n  Barcode format: EXACT match ({_exact_match}/{len(_raw_barcodes)} "
          f"of CSV rows found in project)")
    _bc_to_score = dict(zip(_raw_barcodes, _raw_scores))
else:
    # Attempt 2: short barcode with -1  (AAACCCAAGTCCCAGC-1)
    _prefixed = [BARCODE_PREFIX + b for b in _raw_barcodes]
    _match_short = sum(1 for b in _prefixed if b in _proj_set)
    if _match_short > 0:
        print(f"\n  Barcode format: SHORT with -1  "
              f"(prefixed {_match_short}/{len(_prefixed)} found)")
        _bc_to_score = dict(zip(_prefixed, _raw_scores))
    else:
        # Attempt 3: bare barcode (no -1 suffix)
        _with_suffix = [BARCODE_PREFIX + b + "-1" for b in _raw_barcodes]
        _match_bare = sum(1 for b in _with_suffix if b in _proj_set)
        if _match_bare > 0:
            print(f"\n  Barcode format: BARE (no -1 suffix)  "
                  f"(prefixed+suffix: {_match_bare}/{len(_with_suffix)} found)")
            _bc_to_score = dict(zip(_with_suffix, _raw_scores))
        else:
            print("\n  [WARN] No barcodes matched in any format.")
            print("  Sample CSV barcodes (first 5):")
            for _b in _raw_barcodes[:5]:
                print(f"    {_b!r}")
            print("  Sample project barcodes (first 5):")
            for _b in _hcipsc_barcodes[:5]:
                print(f"    {_b!r}")
            raise ValueError(
                "No OCLR barcodes matched the project's hCiPSC barcodes.\n"
                "Check the CSV barcode format and BARCODE_PREFIX constant."
            )

# Build aligned score series (NaN for project cells not in CSV)
_aligned_scores = np.array(
    [_bc_to_score.get(bc, np.nan) for bc in _hcipsc_barcodes],
    dtype=np.float64,
)

n_matched   = int(np.isfinite(_aligned_scores).sum())
n_unmatched = int(np.isnan(_aligned_scores).sum())
n_csv_extra = sum(1 for b in _bc_to_score if b not in _proj_set)

print(f"\n  --- Barcode alignment report ---")
print(f"  Project hCiPSC cells      : {n_hcipsc:,}")
print(f"  CSV entries               : {len(_raw_barcodes):,}")
print(f"  Matched (CSV ∩ project)   : {n_matched:,}")
print(f"  Unmatched in project      : {n_unmatched:,}  "
      f"(project cells absent from CSV — will be EXCLUDED from terminal set)")
print(f"  Extra in CSV (not in proj): {n_csv_extra:,}  (ignored)")

if n_unmatched > 0:
    print(f"\n  [INFO] {n_unmatched} project hCiPSC cells have no OCLR score.")
    print(f"  These cells CANNOT be selected as terminal-set members.")
    print(f"  They will receive fate_ips=0, fate_other=1 (complementary fate).")

# Handle NaN: cells without a score receive score=0 (excluded from terminal set)
# This is safe because the terminal set is defined as "top-X% of scored cells".
_aligned_scores_filled = np.where(
    np.isfinite(_aligned_scores), _aligned_scores, 0.0
)

# Linearly rescale to [0, 1]  (the paper specifies this explicitly)
_s_min = _aligned_scores_filled[_aligned_scores_filled > 0].min() if n_matched > 0 else 0
_s_max = _aligned_scores_filled.max()
if _s_max > _s_min:
    _aligned_scores_01 = np.where(
        np.isfinite(_aligned_scores),
        (_aligned_scores - _s_min) / (_s_max - _s_min),
        0.0,
    )
    _aligned_scores_01 = np.clip(_aligned_scores_01, 0.0, 1.0)
else:
    print("  [WARN] Score range is zero — all scores are identical.")
    _aligned_scores_01 = np.zeros_like(_aligned_scores_filled)

oclr_df = pd.DataFrame(
    {
        "oclr_score":          _aligned_scores_01.astype(np.float32),
        "oclr_score_raw":      _aligned_scores.astype(np.float32),
        "matched_in_csv":      np.isfinite(_aligned_scores).astype(int),
        "source":              "oclr_csv",
    },
    index=_hcipsc_barcodes,
)
oclr_df.index.name = "barcode"
_oclr_scores = _aligned_scores_01   # use for all downstream thresholding

# =============================================================================
# 4.  OCLR SCORE DISTRIBUTION ANALYSIS
# =============================================================================
# _scored_mask selects cells that have a real OCLR score from the CSV
# (matched_in_csv == 1).  Cells not present in the CSV (NaN score, filled to
# 0.0) are excluded from all distribution statistics and quantile thresholds.
# They are never selected as terminal-set members.
_scored_mask = oclr_df["matched_in_csv"].values.astype(bool)

n_scored   = int(_scored_mask.sum())
n_unscored = len(oclr_df) - n_scored

print("\n--- OCLR score distribution (scored hCiPSC cells only) ---")
print(f"  Total hCiPSC cells in pool : {len(oclr_df):,}")
print(f"  Scored cells (used for threshold) : {n_scored:,}")
if n_unscored > 0:
    print(f"  Unmatched / unscored cells : {n_unscored:,}  "
          f"(excluded from threshold — will be non-terminal by definition)")
    print(f"  NOTE: threshold candidates are computed on scored cells only.")

_scores_all   = oclr_df["oclr_score"].values.astype(float)
_scores_scored = _scores_all[_scored_mask]   # the only array used for statistics

_stats = {
    "n_cells":   n_scored,
    "mean":      float(np.mean(_scores_scored)),
    "std":       float(np.std(_scores_scored)),
    "min":       float(np.min(_scores_scored)),
    "q25":       float(np.quantile(_scores_scored, 0.25)),
    "median":    float(np.quantile(_scores_scored, 0.50)),
    "q75":       float(np.quantile(_scores_scored, 0.75)),
    "q90":       float(np.quantile(_scores_scored, 0.90)),
    "q95":       float(np.quantile(_scores_scored, 0.95)),
    "max":       float(np.max(_scores_scored)),
}
print(f"  mean±std : {_stats['mean']:.4f} ± {_stats['std']:.4f}")
print(f"  [min, Q25, median, Q75, Q90, Q95, max] =")
print(f"    [{_stats['min']:.4f}, {_stats['q25']:.4f}, {_stats['median']:.4f}, "
      f"{_stats['q75']:.4f}, {_stats['q90']:.4f}, {_stats['q95']:.4f}, {_stats['max']:.4f}]")

# =============================================================================
# 5.  THRESHOLD CANDIDATES
# =============================================================================
print("\n--- Threshold candidates (computed on scored cells only) ---")
print(f"  {'Label':<15}  {'Quantile':>9}  {'Cutoff':>10}  "
      f"{'N scored (≥)':>13}  {'% of scored':>11}")

_threshold_rows = []
for _label, _q in THRESHOLD_CANDIDATES:
    _cutoff  = float(np.quantile(_scores_scored, _q))       # scored cells only
    _n_above = int((_scores_scored >= _cutoff).sum())        # scored cells only
    _pct     = 100.0 * _n_above / max(n_scored, 1)
    _marker  = " ← FINAL" if _label == FINAL_THRESHOLD_KEY else ""
    print(f"  {_label:<15}  {_q:>9.3f}  {_cutoff:>10.4f}  {_n_above:>12,}  {_pct:>9.1f}%{_marker}")
    _threshold_rows.append({
        "threshold_key":   _label,
        "quantile":        _q,
        "cutoff_value":    _cutoff,
        "n_selected":      _n_above,
        "pct_selected":    _pct,
        "is_final":        (_label == FINAL_THRESHOLD_KEY),
        "note":            (
            "FINAL SELECTION: high-specificity endpoint for trajectory benchmarking"
            if _label == FINAL_THRESHOLD_KEY
            else ""
        ),
    })

# Extract the final threshold
_final_thr_row = next(r for r in _threshold_rows if r["is_final"])
FINAL_CUTOFF        = _final_thr_row["cutoff_value"]
FINAL_N_SELECTED    = _final_thr_row["n_selected"]
FINAL_PCT_SELECTED  = _final_thr_row["pct_selected"]

print(f"\n  Final threshold selected  : {FINAL_THRESHOLD_KEY}")
print(f"  Final OCLR score cutoff   : {FINAL_CUTOFF:.4f}")
print(f"  Final terminal-set size   : {FINAL_N_SELECTED:,} cells "
      f"({FINAL_PCT_SELECTED:.1f}% of {n_scored:,} scored cells"
      f"{f', out of {n_hcipsc} total hCiPSC in pool' if n_unscored > 0 else ''})")

if FINAL_N_SELECTED < MIN_TERMINAL_CELLS:
    raise RuntimeError(
        f"Final terminal set has only {FINAL_N_SELECTED} cells "
        f"(minimum required: {MIN_TERMINAL_CELLS}).\n"
        f"Consider using a less restrictive FINAL_THRESHOLD_KEY."
    )

# =============================================================================
# 6.  DEFINE TERMINAL SET
# =============================================================================
# is_oclr_terminal:  True  ↔  cell has a real OCLR score (scored_mask)
#                              AND oclr_score ≥ FINAL_CUTOFF
#
# Unmatched cells are explicitly excluded even if their filled score of 0.0
# happened to satisfy the cutoff (which cannot happen for a q90 cutoff, but
# we guard against it explicitly for correctness).
_above_cutoff = (oclr_df["oclr_score"].values >= FINAL_CUTOFF)
oclr_df["is_oclr_terminal"] = (_above_cutoff & _scored_mask).astype(int)

# Barcode list for the terminal set
_terminal_barcodes = oclr_df.index[oclr_df["is_oclr_terminal"] == 1].tolist()

print(f"\n  Sanity check:")
print(f"    is_oclr_terminal == 1 : {int(oclr_df['is_oclr_terminal'].sum()):,}")
print(f"    is_oclr_terminal == 0 : {int((oclr_df['is_oclr_terminal'] == 0).sum()):,}")
print(f"    Sum = total hCiPSC    : {len(oclr_df):,}  ({'✓' if len(oclr_df) == n_hcipsc else '!!'})")

# =============================================================================
# 7.  DIAGNOSTIC PLOT
# =============================================================================
print("\n--- Generating diagnostic plot ---")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# ── Left: histogram of scored cells only ────────────────────────────────────
# _scores_scored contains only cells with a real OCLR score (matched_in_csv==1
# for the CSV path, or all cells for the proxy path).  This is exactly the
# distribution used for quantile thresholding, so the plot is consistent with
# the actual terminal-set definition.
ax0 = axes[0]
_bins = np.linspace(0, 1, 60)
ax0.hist(_scores_scored, bins=_bins,
         color="#7fbbe3", edgecolor="white", linewidth=0.3, alpha=0.9)
ax0.set_xlabel("OCLR stemness score (rescaled [0, 1])")
ax0.set_ylabel("Cell count")
ax0.set_title(f"OCLR score distribution — scored cells\n"
              f"{ENDPOINT_SAMPLE_ID}  (n={n_scored:,})")

# Add threshold candidate lines — cutoffs re-computed from _scores_scored
# to guarantee exact consistency with the thresholding section above.
_colors = ["#e63946", "#f77f00", "#457b9d", "#2a9d8f"]
for (_lbl, _q), _color in zip(THRESHOLD_CANDIDATES, _colors):
    _cut = float(np.quantile(_scores_scored, _q))
    _lw  = 2.5 if _lbl == FINAL_THRESHOLD_KEY else 1.0
    _ls  = "-"  if _lbl == FINAL_THRESHOLD_KEY else "--"
    ax0.axvline(_cut, color=_color, linewidth=_lw, linestyle=_ls,
                label=f"{_lbl}  (≥{_cut:.3f})")

ax0.legend(fontsize=9, loc="upper left")

# ── Right: empirical CDF of scored cells only ────────────────────────────────
ax1 = axes[1]
_sorted_scores = np.sort(_scores_scored)
_cdf = np.arange(1, len(_sorted_scores) + 1) / len(_sorted_scores)
ax1.plot(_sorted_scores, 1 - _cdf, color="#1d3557", linewidth=1.5)
ax1.set_xlabel("OCLR stemness score (rescaled [0, 1])")
ax1.set_ylabel("Fraction of scored cells with score ≥ x  (1 − CDF)")
ax1.set_title(f"OCLR score complementary CDF — scored cells\n{ENDPOINT_SAMPLE_ID}")

for (_lbl, _q), _color in zip(THRESHOLD_CANDIDATES, _colors):
    _cut  = float(np.quantile(_scores_scored, _q))
    _frac = 1 - _q
    _lw   = 2.5 if _lbl == FINAL_THRESHOLD_KEY else 1.0
    _ls   = "-"  if _lbl == FINAL_THRESHOLD_KEY else "--"
    ax1.axvline(_cut, color=_color, linewidth=_lw, linestyle=_ls,
                label=f"{_lbl}  ({100*_frac:.0f}% of scored)")

ax1.set_ylim(0, 0.50)
ax1.legend(fontsize=9)

_pool_note = (f"  [{n_unscored:,} unscored / unmatched cells excluded from plot]"
              if n_unscored > 0 else "")
fig.suptitle(
    f"OCLR-based iPS terminal endpoint — {ENDPOINT_SAMPLE_ID}\n"
    f"Final selection: {FINAL_THRESHOLD_KEY}  "
    f"(cutoff={FINAL_CUTOFF:.4f},  n={FINAL_N_SELECTED:,} / {n_scored:,} scored cells, "
    f"{FINAL_PCT_SELECTED:.1f}%){_pool_note}",
    fontsize=11,
)
fig.tight_layout()
_plot_path = FIGURES_DIR / f"{TIMESTAMP}_oclr_score_distribution.png"
fig.savefig(str(_plot_path), dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {_plot_path}")

# =============================================================================
# 8.  SAVE OUTPUT FILES
# =============================================================================
print("\n--- Saving outputs ---")

# 8a: All hCiPSC cells with OCLR scores
_scores_out = PROCESSED_DIR / f"{TIMESTAMP}_oclr_score_all_hcipsc.csv"
oclr_df.to_csv(str(_scores_out))
print(f"  OCLR scores (all hCiPSC) → {_scores_out.name}")

# 8b: Threshold summary table
_thr_df = pd.DataFrame(_threshold_rows)
_thr_df["method_note"] = (
    "OCLR stemness score thresholding within GSM7230012_hCiPSCs-0618.\n"
    "Top-N% defined by quantile of the score distribution in this sample."
)
_thr_df["barcode_prefix"]   = BARCODE_PREFIX
_thr_df["n_hcipsc_total"]   = n_hcipsc
_thr_df["n_csv_matched"]    = n_matched
_thr_df["n_csv_unmatched"]  = n_unmatched
_thr_df["score_min"]        = _stats["min"]
_thr_df["score_max"]        = _stats["max"]
_thr_df["score_mean"]       = _stats["mean"]
_thr_df["score_median"]     = _stats["median"]
_thr_df["oclr_score_source"] = "real_csv"
_thr_df["pool_method"]       = _POOL_METHOD

_thr_out = PROCESSED_DIR / f"{TIMESTAMP}_oclr_threshold_summary.csv"
_thr_df.to_csv(str(_thr_out), index=False)
print(f"  Threshold summary         → {_thr_out.name}")

# 8c: Final terminal barcode list
# Format: one barcode per line (full project barcode)
# Column header: "barcode"
# Extra column: oclr_score  (for traceability)
_bc_df = oclr_df[oclr_df["is_oclr_terminal"] == 1][["oclr_score"]].copy()
_bc_df.index.name = "barcode"
_bc_out = PROCESSED_DIR / f"{TIMESTAMP}_oclr_endpoint_barcodes.txt"
_bc_df.to_csv(str(_bc_out), sep="\t")
print(f"  Terminal barcode list     → {_bc_out.name}")
print(f"    N terminal cells : {len(_bc_df):,}")
print(f"    Sample entries   :")
for _bc in _bc_df.head(5).itertuples():
    print(f"      {_bc.Index}  oclr={_bc.oclr_score:.4f}")

# =============================================================================
# 9.  SUMMARY
# =============================================================================
print(f"\n{'='*65}")
print(f"[{TIMESTAMP}]  Step 02b2 COMPLETE")
print(f"\n  OCLR Endpoint Definition Summary")
print(f"  ─────────────────────────────────────────────────────────────")
print(f"  Source dataset          : {ENDPOINT_SAMPLE_ID}")
print(f"  OCLR score source       : PROVIDED CSV (real published scores)")
print(f"  CSV file                : {_oclr_found.name}")
print(f"  CSV barcodes matched    : {n_matched:,} / {n_hcipsc:,}")
print(f"  Unmatched in project    : {n_unmatched:,}")
print(f"\n  Score distribution      : [{_stats['min']:.3f}, {_stats['max']:.3f}]")
print(f"  Median OCLR score       : {_stats['median']:.4f}")
print(f"\n  Threshold candidates evaluated:")
for _r in _threshold_rows:
    _marker = "← FINAL" if _r["is_final"] else ""
    print(f"    {_r['threshold_key']:<15}  q={_r['quantile']:.2f}  "
          f"cut={_r['cutoff_value']:.4f}  n={_r['n_selected']:>5,}  "
          f"({_r['pct_selected']:.1f}%)  {_marker}")
print(f"\n  FINAL terminal set:")
print(f"    Key         : {FINAL_THRESHOLD_KEY}")
print(f"    Cutoff      : {FINAL_CUTOFF:.4f}")
print(f"    N cells     : {FINAL_N_SELECTED:,} / {n_scored:,} scored cells  "
      f"({FINAL_PCT_SELECTED:.1f}% of scored)"
      + (f"\n               out of {n_hcipsc:,} total hCiPSC cells in pool"
         if n_unscored > 0 else ""))
print(f"    Barcodes    : {_bc_out.name}")
print(f"\n  Project-specific cutoff note")
print(f"    The top-{int((1-THRESHOLD_CANDIDATES[1][1])*100)}% threshold is a PROJECT DECISION, not a threshold")
print(f"    specified by the original OCLR paper (Mäkinen V-P et al.).")
print(f"    It defines a high-confidence binary evaluation set for")
print(f"    trajectory benchmarking, independent of WOT / CellRank2.")
print(f"    The continuous OCLR score is the primary biological signal.")
print(f"    If the distribution is clearly bimodal in the upper tail,")
print(f"    switch to FINAL_THRESHOLD_KEY='top_5pct' and re-run.")
print(f"\n  Output files:")
print(f"    {_scores_out.name}")
print(f"    {_thr_out.name}")
print(f"    {_bc_out.name}")
print(f"    {_plot_path.name}")
print(f"\n  Next: run 02c_compute_transport.py  (reads the barcode file above)")
print(f"{'='*65}")
