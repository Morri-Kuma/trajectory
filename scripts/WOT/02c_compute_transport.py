#!/usr/bin/env python3
"""
02c_compute_transport.py  [v13 — multi-fate biological terminal design]
==========================================================================================
Project : Comparative Study of Trajectory Inference Models for Chemical iPSC
          Reprogramming
Dataset : GSE230659 (human, Liuyang et al. 2023 Cell Stem Cell)
Step    : 2c — WOT transport-map computation, fate probabilities, trajectory trends
Env     : conda activate traj_env

Multi-fate terminal design (v13)
----------------------------------
  This version follows the original WOT reprogramming paper logic (Schiebinger et al.
  2019 Cell) in which multiple biologically distinct terminal populations compete as
  fate targets, making the backward propagation non-degenerate without relying on any
  scoring model to split the terminal pool.

  Terminal fate sets (biological, independent of OCLR / WOT outputs):
    fate_hcipsc        — ALL hCiPSC cells (stage_std == "hCiPSC")
    fate_stageIII_late — StageIII cells at their latest observed timepoint

  Because the two fate sets partition (or at least cover) the terminal cells, the
  backward-propagated fate probabilities p_hcipsc and p_stageIII_late are non-degenerate
  and sum to ≤ 1 for all earlier cells.

  OCLR scores (Mäkinen V-P et al.) are loaded from the pre-computed file produced by
  02b2_define_oclr_endpoint.py and attached as obs annotations ONLY — they serve as an
  external biological anchor for downstream evaluation (script 05) and are NOT used to
  define or split any terminal fate set in this script.

  Fate computation uses the WOT Python API:
    wot.tmap.TransportMapModel.from_directory()  +  .get_fate_probabilities()

  This is the canonical WOT approach and correctly handles fate sets defined at different
  timepoints (e.g. fate_hcipsc at day30, fate_stageIII_late at day21).

  Run 02b2_define_oclr_endpoint.py BEFORE this script (for OCLR annotation).

Version history
---------------
  v9 : Two-fate design introduced (partition final cells into ips/other).
  v10: 6-class UMAP label.
  v11: CLI OT pipeline unchanged.
  v12: Endpoint replaced from marker genes → OCLR stemness scores (OCLR split).
  v13: OCLR no longer defines the fate split.  Multi-fate biological design using
       hCiPSC + StageIII_late as competing fates.  WOT Python API for fate computation.
       OCLR scores attached as obs annotation only.

Usage
-----
    conda activate traj_env
    cd C:\\Users\\37620\\trajectory
    python scripts/WOT/02c_compute_transport.py
"""

# =============================================================================
# 0.  CONFIGURATION  ←  EDIT THESE FLAGS TO CONTROL BEHAVIOUR
# =============================================================================
from datetime import datetime
from pathlib import Path
import time as _time

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")
T0        = _time.time()

# ── Project root ──────────────────────────────────────────────────────────────
_candidates = [
    Path(r"C:\Users\37620\trajectory"),
    Path("/sessions/cool-admiring-hawking/mnt/trajectory"),
    Path(__file__).resolve().parents[2],
]
PROJECT_ROOT  = next((p for p in _candidates if p.exists()), _candidates[-1])
PROCESSED_DIR = PROJECT_ROOT / "data"    / "processed"
FIGURES_DIR   = PROJECT_ROOT / "results" / "figures"
METRICS_DIR   = PROJECT_ROOT / "results" / "metrics"
TMAPS_DIR     = PROJECT_ROOT / "results" / "tmaps" / f"{TIMESTAMP}_tmaps"
WOT_INPUTS    = TMAPS_DIR / "wot_inputs"
for _d in [FIGURES_DIR, METRICS_DIR, TMAPS_DIR, WOT_INPUTS]:
    _d.mkdir(parents=True, exist_ok=True)

def find_latest(pattern: str) -> Path:
    matches = sorted(PROCESSED_DIR.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No files matching '{pattern}' in {PROCESSED_DIR}.\n"
            f"Run steps 02a and 02b first."
        )
    return matches[-1]

H5AD_PATH = find_latest("*_GSE230659_wot_gr.h5ad")

# ══════════════════════════════════════════════════════════════════════════════
# WOT CLI PARAMETERS
# ══════════════════════════════════════════════════════════════════════════════
EPSILON:   float = 0.05
LAMBDA1:   float = 1.0
LAMBDA2:   float = 50.0
LOCAL_PCA: int   = 30

USE_GROWTH_RATES: bool = True

# ══════════════════════════════════════════════════════════════════════════════
# DEBUG CONTROLS
# ══════════════════════════════════════════════════════════════════════════════
DEBUG_SUBSAMPLE_PER_TIMEPOINT: int | None = None
DEBUG_ONLY_FIRST_N_PAIRS:      int | None = None
FORCE_DENSE_FLOAT64:           bool       = False

# ══════════════════════════════════════════════════════════════════════════════
# FATE COMPUTATION PARAMETERS  (multi-fate biological design)
# ══════════════════════════════════════════════════════════════════════════════
#
# TERMINAL FATE SETS (biological, no OCLR split)
# -----------------------------------------------
# Two biologically distinct terminal populations are defined using the stage
# annotation in adata.obs["stage_std"].  This follows the WOT reprogramming
# paper design (Schiebinger et al. 2019 Cell) in which multiple cell types at
# the endpoint serve as competing fate targets.
#
#   fate_hcipsc        — ALL cells with stage_std == "hCiPSC"
#                        (fully reprogrammed iPSCs; typically at day30)
#   fate_stageIII_late — StageIII cells at their LATEST observed timepoint
#                        (the most advanced non-iPSC intermediate; may be
#                         at a different day than hCiPSC)
#
# WHY THIS DESIGN IS NON-DEGENERATE
# ----------------------------------
# If the entire terminal pool were one fate (e.g., all day30 cells), the
# backward propagation collapses to p=1 for all cells.  Having two distinct
# competing fates makes p_hcipsc non-trivial and biologically interpretable:
#   p_hcipsc(i)        = probability that cell i eventually becomes hCiPSC
#   p_stageIII_late(i) = probability that cell i terminates as StageIII
#
# ROLE OF OCLR SCORES
# --------------------
# OCLR scores (Mäkinen V-P et al.) are loaded as obs annotations for external
# evaluation only.  They are NOT used to define or split any fate set here.
# Use script 05 (compare_wot_cr2_oclr.py) to correlate p_hcipsc with OCLR.

# Minimum number of cells required in each fate set.
# If either falls below this the script stops with a clear error.
MIN_FATE_CELLS: int = 20

# ── OCLR endpoint sample ID (for annotation only) ────────────────────────────
# OCLR scores are loaded by load_oclr_endpoint() and attached to adata.obs.
# This constant is used only for validation / logging.
OCLR_ENDPOINT_SAMPLE_ID: str = "GSM7230012_hCiPSCs-0618"

# ══════════════════════════════════════════════════════════════════════════════
# DOWNSTREAM STEP TOGGLES
# ══════════════════════════════════════════════════════════════════════════════
RUN_FATE_COMPUTATION: bool = True
RUN_PLOTTING:         bool = True

GR_CLIP = 1.0

# =============================================================================
# PLOTTING PALETTE  (stage-based, no OCLR split)
# =============================================================================
# hCiPSC   — all fully reprogrammed cells (fate_hcipsc target)
# StageIII — StageIII cells (fate_stageIII_late target at latest day)
# StageI/II — early-intermediate cells
# Other    — any other stage label
STAGE_PALETTE = {
    "hCiPSC":   "#d62728",   # strong red
    "StageIII": "#2ca02c",   # green
    "StageII":  "#ff7f0e",   # orange
    "StageI":   "#1f77b4",   # blue
    "Other":    "#BDBDBD",   # grey
}

# =============================================================================
# 1.  IMPORTS
# =============================================================================
import sys
import re
import subprocess
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
from scipy.spatial.distance import cdist as _pdist
from sklearn.decomposition import PCA as _SklearnPCA
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import wot
except ImportError:
    raise ImportError("WOT not installed.\n  conda activate traj_env && pip install wot")

try:
    import importlib.metadata as _imm
    _wot_ver = _imm.version("wot")
except Exception:
    _wot_ver = getattr(wot, "__version__", "unknown")

_SEP  = "─" * 70
_SEP2 = "· " * 35

print(f"[{TIMESTAMP}]  WOT Transport Map Computation  (v13 — multi-fate biological design)")
print(f"  wot version     : {_wot_ver}")
print(f"  Input AnnData   : {H5AD_PATH.name}")
print(f"  Output dir      : {TMAPS_DIR}")
print(f"\n  WOT CLI parameters")
print(f"    EPSILON          : {EPSILON}")
print(f"    LAMBDA1/LAMBDA2  : {LAMBDA1} / {LAMBDA2}")
print(f"    LOCAL_PCA        : {LOCAL_PCA}")
print(f"    USE_GROWTH_RATES : {USE_GROWTH_RATES}")
print(f"\n  Fate computation design (multi-fate biological, WOT Python API)")
print(f"    fate_hcipsc        : ALL hCiPSC cells  (stage_std == 'hCiPSC')")
print(f"    fate_stageIII_late : StageIII cells at their latest timepoint")
print(f"    p_hcipsc           : probability of joining the hCiPSC fate")
print(f"    p_stageIII_late    : probability of joining the StageIII_late fate")
print(f"    OCLR scores        : loaded as obs annotation only (external anchor)")
print(f"\n  Debug controls")
print(f"    DEBUG_SUBSAMPLE_PER_TIMEPOINT : {DEBUG_SUBSAMPLE_PER_TIMEPOINT}")
print(f"    DEBUG_ONLY_FIRST_N_PAIRS      : {DEBUG_ONLY_FIRST_N_PAIRS}")
print(f"    FORCE_DENSE_FLOAT64           : {FORCE_DENSE_FLOAT64}")

# =============================================================================
# HELPER: stage-label normalizer
# =============================================================================
_STAGE_ROMAN = {1: "StageI", 2: "StageII", 3: "StageIII", 4: "StageIV"}

def normalize_stage(s) -> str:
    """
    Map any stage-string variant to one of:
      'StageI', 'StageII', 'StageIII', 'hCiPSC', or the original string.

    Examples
    --------
    'hCiPSCs' → 'hCiPSC'
    'hcipsc'  → 'hCiPSC'
    'iPSC'    → 'hCiPSC'
    'StageII' → 'StageII'
    'stage2'  → 'StageII'
    'STAGEIII'→ 'StageIII'
    """
    if s is None or (isinstance(s, float) and np.isnan(s)):
        return "Unknown"
    sl = str(s).strip().lower()
    # hCiPSC variants
    if any(kw in sl for kw in ("cipsc", "cipsc", "hcipsc", "ipsc")):
        return "hCiPSC"
    # Stage with Roman or Arabic numerals
    m = re.search(r"stage\s*([ivxIVX\d]+)", sl)
    if m:
        raw = m.group(1)
        # Arabic
        if raw.isdigit():
            n = int(raw)
        else:
            # Roman numeral: map i=1, ii=2, iii=3, iv=4
            _rom = {"i": 1, "ii": 2, "iii": 3, "iv": 4}
            n = _rom.get(raw.lower(), 0)
        return _STAGE_ROMAN.get(n, f"Stage{raw.upper()}")
    return str(s).strip()   # fallback: return as-is


# =============================================================================
# HELPER: OCLR-defined iPSC fate set
# =============================================================================
def load_oclr_endpoint(adata) -> "tuple[pd.Series, pd.Series]":
    """
    Load the OCLR-based iPS terminal barcode list produced by
    02b2_define_oclr_endpoint.py and return (oclr_score_series, oclr_subset_mask).

    The endpoint is defined by:
      - One-class logistic regression (OCLR) stemness scores for every cell in
        GSM7230012_hCiPSCs-0618 (pre-computed, independent of WOT / CellRank2).
      - High-confidence iPS cells = top 10% by OCLR score within that sample
        (default in 02b2; see FINAL_THRESHOLD_KEY there for the exact quantile).

    This function searches PROCESSED_DIR for the most recent
    *_oclr_endpoint_barcodes.txt file and aligns it to adata.obs_names.

    Barcode format compatibility
    ----------------------------
    The barcode file contains full project barcodes:
        GSM7230012_hCiPSCs-0618_AAACCCAAGTCCCAGC-1
    which must match adata.obs_names exactly.  Any barcodes in the file that
    are not in adata.obs_names are reported as mismatches (never silently dropped).

    Parameters
    ----------
    adata : AnnData — the full project AnnData (all cells, all stages)

    Returns
    -------
    oclr_score_series : pd.Series[float]  — OCLR score per cell (NaN for non-hCiPSC)
    subset_mask       : pd.Series[bool]   — True for cells in the OCLR top-10% set
    """
    # ── Find the most recent barcode list ─────────────────────────────────────
    _bc_files = sorted(PROCESSED_DIR.glob("*_oclr_endpoint_barcodes.txt"))
    if not _bc_files:
        raise FileNotFoundError(
            f"No OCLR endpoint barcode file found in {PROCESSED_DIR}.\n"
            f"Run 02b2_define_oclr_endpoint.py first to generate it.\n"
            f"Expected pattern: *_oclr_endpoint_barcodes.txt"
        )
    _bc_file = _bc_files[-1]   # most recent (lexicographic sort on timestamp prefix)
    print(f"  OCLR barcode file : {_bc_file.name}")

    # ── Load barcode list ─────────────────────────────────────────────────────
    _bc_df = pd.read_csv(str(_bc_file), sep="\t", index_col=0)
    _oclr_terminal_barcodes = set(_bc_df.index.astype(str).tolist())
    _oclr_score_map  = {}
    if "oclr_score" in _bc_df.columns:
        _oclr_score_map = _bc_df["oclr_score"].to_dict()
    print(f"  Terminal barcodes in file : {len(_oclr_terminal_barcodes):,}")

    # ── Also load the per-cell OCLR score file if available ───────────────────
    _score_files = sorted(PROCESSED_DIR.glob("*_oclr_score_all_hcipsc.csv"))
    if _score_files:
        _score_file = _score_files[-1]
        print(f"  OCLR score file (all hCiPSC) : {_score_file.name}")
        _all_scores_df = pd.read_csv(str(_score_file), index_col=0)
        if "oclr_score" in _all_scores_df.columns:
            _oclr_score_map.update(
                _all_scores_df["oclr_score"].to_dict()
            )

    # ── Align to adata.obs_names ──────────────────────────────────────────────
    _proj_barcodes = adata.obs_names.astype(str).tolist()
    _proj_set      = set(_proj_barcodes)

    # Barcodes in file that are NOT in adata
    _extra_in_file = _oclr_terminal_barcodes - _proj_set
    # Barcodes in adata that ARE in the terminal set
    _matched_terminal = _oclr_terminal_barcodes & _proj_set

    if _extra_in_file:
        print(f"\n  [WARN] {len(_extra_in_file)} terminal barcodes in file "
              f"not found in adata.obs_names:")
        for _b in sorted(_extra_in_file)[:5]:
            print(f"    {_b!r}  (not in project AnnData)")
        if len(_extra_in_file) > 5:
            print(f"    ... ({len(_extra_in_file) - 5} more not shown)")
    else:
        print(f"  Barcode alignment: all terminal barcodes found in adata  ✓")

    print(f"\n  Matched terminal barcodes  : {len(_matched_terminal):,}")
    print(f"  Unmatched (file − adata)   : {len(_extra_in_file):,}")

    # ── Build per-cell OCLR score series (NaN for non-hCiPSC / unscored cells)─
    _oclr_scores = pd.Series(
        [_oclr_score_map.get(bc, np.nan) for bc in _proj_barcodes],
        index=_proj_barcodes,
        name="oclr_score",
        dtype=float,
    )

    # ── Build subset mask (True = in OCLR terminal set) ───────────────────────
    subset_mask = pd.Series(
        [bc in _matched_terminal for bc in _proj_barcodes],
        index=_proj_barcodes,
        name="oclr_ips_subset_mask",
        dtype=bool,
    )

    return _oclr_scores, subset_mask


# =============================================================================
# 2.  LOAD ANNDATA AND VALIDATE INPUTS
# =============================================================================
print(f"\n{_SEP}")
print("STEP 2  —  Load AnnData and validate inputs")
print(_SEP)

adata_full = sc.read_h5ad(str(H5AD_PATH))
_is_sparse = sp.issparse(adata_full.X)
print(f"  Loaded : {adata_full.n_obs:,} cells × {adata_full.n_vars:,} genes")
print(f"  obsm   : {list(adata_full.obsm.keys())}")
print(f"  obs    : {list(adata_full.obs.columns)}")
print(f"  X type : {'sparse' if _is_sparse else 'dense'}  "
      f"dtype={adata_full.X.data.dtype if _is_sparse else adata_full.X.dtype}")

# ── Ensure 'day' column ───────────────────────────────────────────────────────
if "day" not in adata_full.obs.columns:
    if "abs_day" in adata_full.obs.columns:
        adata_full.obs["day"] = adata_full.obs["abs_day"].astype(float)
        print("  Created obs['day'] from 'abs_day'")
    else:
        raise ValueError("Neither 'day' nor 'abs_day' in obs. Re-run 02a/02b.")

REQUIRED_OBS = ["day", "growth_rate", "stage", "stage_day_label"]
_missing = [c for c in REQUIRED_OBS if c not in adata_full.obs.columns]
if _missing:
    raise ValueError(f"Required obs columns missing: {_missing}")
print(f"  Required obs columns : OK  {REQUIRED_OBS}")

# ── NaN/Inf quick check ───────────────────────────────────────────────────────
_samp = (adata_full.X[:500].toarray() if _is_sparse else np.asarray(adata_full.X[:500]))
if np.isnan(_samp).any() or np.isinf(_samp).any():
    raise ValueError("adata.X contains NaN/Inf in first 500 cells. Re-run preprocessing.")
print(f"  NaN/Inf check (500-cell sample) : PASSED")

# ── Growth-rate zeroing / clipping ────────────────────────────────────────────
_gr_raw     = adata_full.obs["growth_rate"].values.astype(float)
_gr_clipped = np.clip(_gr_raw, -GR_CLIP, GR_CLIP)
adata_full.obs["growth_rate"] = _gr_clipped
if not USE_GROWTH_RATES:
    adata_full.obs["growth_rate"] = 0.0
    print("  Growth rates : set to 0.0  (USE_GROWTH_RATES=False)")
else:
    print(f"  Growth rates : [{_gr_clipped.min():.4f}, {_gr_clipped.max():.4f}]  (clipped)")

_sorted_days = sorted(adata_full.obs["day"].unique())
print(f"\n  Timepoints ({len(_sorted_days)}):")
for _d in _sorted_days:
    _n   = int((adata_full.obs["day"] == _d).sum())
    _lbl = adata_full.obs.loc[adata_full.obs["day"] == _d, "stage_day_label"].iloc[0]
    print(f"    Day {_d:7.3f}  |  {_n:>6,} cells  |  {_lbl}")

# =============================================================================
# 3a.  adata.X CONTENT AUDIT
# =============================================================================
print(f"\n{_SEP}")
print("STEP 3a —  adata.X content audit")
print(_SEP)

_rng_audit = np.random.default_rng(99)
_audit_n   = min(2_000, adata_full.n_obs)
_audit_idx = _rng_audit.choice(adata_full.n_obs, _audit_n, replace=False)
_audit_raw = adata_full.X[_audit_idx]
_audit_mat = (_audit_raw.toarray().astype(np.float64) if sp.issparse(_audit_raw)
              else np.asarray(_audit_raw, dtype=np.float64))

_flat      = _audit_mat.ravel()
_nz        = _flat[_flat != 0.0]
frac_zero  = float((_flat == 0.0).sum()) / _flat.size
frac_neg   = float((_flat  < 0.0).sum()) / _flat.size

print(f"  Sample: {_audit_n} cells  |  dtype: "
      f"{adata_full.X.data.dtype if _is_sparse else adata_full.X.dtype}")
print(f"  Global min/max  : {_flat.min():.4f} / {_flat.max():.4f}")
print(f"  Fraction zero   : {frac_zero:.3f}")
print(f"  Fraction neg    : {frac_neg:.6f}")

_nz_max   = float(_nz.max()) if len(_nz) > 0 else 0.0
_frac_int = float((np.abs(_nz - np.round(_nz)) < 0.01).sum()) / max(len(_nz), 1)
_row_sums_audit = _audit_mat.sum(axis=1)
_row_cv         = _row_sums_audit.std() / max(_row_sums_audit.mean(), 1e-9)

if len(_nz) > 0:
    _uniq_nz = np.unique(_nz)[:20]
    print(f"  Nonzero min/max : {_nz.min():.4f} / {_nz_max:.4f}")
    print(f"  Frac near-int   : {_frac_int:.3f}  "
          f"({'count-like' if _frac_int > 0.8 else 'transformed'})")
    print(f"  Sample vals     : {' '.join(f'{v:.4f}' for v in _uniq_nz)}")
print(f"  Row-sum mean/std: {_row_sums_audit.mean():.2f} / {_row_sums_audit.std():.2f}  "
      f"CV={_row_cv:.3f}")

_has_neg = frac_neg > 1e-6
if _has_neg:
    _verdict_short = "scaled/centered"
    _verdict = "SCALED/CENTERED — negative values; adata.X looks like z-scored data"
elif _frac_int > 0.80 and _nz_max >= 500:
    _verdict_short = "raw counts"
    _verdict = "RAW COUNTS — integer-like values, large max"
elif _row_sums_audit.mean() > 5_000 and not (_nz_max < 30):
    _verdict_short = "CPM/library-normalised"
    _verdict = "CPM/LIBRARY-NORMALISED — large row sums, not log-transformed"
elif _nz_max < 30 and _frac_int < 0.80 and not _has_neg:
    _verdict_short = "log1p-normalised"
    _verdict = "LOG1P-NORMALISED — bounded, non-integer, non-negative  ✓"
else:
    _verdict_short = "uncertain"
    _verdict = "UNCERTAIN — inspect sampled values above"

print(f"\n  Verdict : {_verdict}")
if _verdict_short == "log1p-normalised":
    print("  ✓  Correct input type for WOT.")
elif _has_neg:
    print("  !! PROBLEM: negative values — WOT requires non-negative input.")
elif _verdict_short in ("raw counts", "CPM/library-normalised"):
    print(f"  !! WARNING: {_verdict_short} data will cause huge local-PCA costs.")
else:
    print("  ?  Verify adata.X is log1p-normalised before trusting OT results.")

# =============================================================================
# 3b.  LOCAL-PCA COST DIAGNOSTIC  (first pair, pre-flight check)
# =============================================================================
print(f"\n{_SEP}")
print("STEP 3b —  Local-PCA cost diagnostic  (first pair sample)")
print(_SEP)

_t0d, _t1d = _sorted_days[0], _sorted_days[1]
_rng_d     = np.random.default_rng(42)
_i0 = np.where(adata_full.obs["day"] == _t0d)[0]
_i1 = np.where(adata_full.obs["day"] == _t1d)[0]
_i0 = _rng_d.choice(_i0, min(500, len(_i0)), replace=False)
_i1 = _rng_d.choice(_i1, min(500, len(_i1)), replace=False)

_X0 = (adata_full.X[_i0].toarray() if _is_sparse else np.asarray(adata_full.X[_i0])).astype(np.float64)
_X1 = (adata_full.X[_i1].toarray() if _is_sparse else np.asarray(adata_full.X[_i1])).astype(np.float64)

_n_pca_d = min(LOCAL_PCA, _X0.shape[1] - 1, len(_i0) + len(_i1) - 1) if LOCAL_PCA > 0 else 0
if _n_pca_d > 0:
    _pca_d  = _SklearnPCA(n_components=_n_pca_d, random_state=0)
    _Xp     = _pca_d.fit_transform(np.vstack([_X0, _X1]))
    _Xp0, _Xp1 = _Xp[:len(_i0)], _Xp[len(_i0):]
    _varexp = _pca_d.explained_variance_ratio_.sum() * 100
    print(f"  PCA {_n_pca_d} components, var. explained {_varexp:.1f}%")
else:
    _Xp0, _Xp1 = _X0, _X1
    print(f"  LOCAL_PCA=0: raw gene-expression distances ({_X0.shape[1]} dims)")

_ns = min(300, len(_Xp0))
_nt = min(300, len(_Xp1))
_s0 = _rng_d.choice(len(_Xp0), _ns, replace=False)
_s1 = _rng_d.choice(len(_Xp1), _nt, replace=False)
_C  = _pdist(_Xp0[_s0], _Xp1[_s1], metric="sqeuclidean")
_c_min, _c_med, _c_p99, _c_max = (float(_C.min()), float(np.median(_C)),
                                   float(np.percentile(_C, 99)), float(_C.max()))

print(f"\n  Cost distribution (Day {_t0d:.3f}→{_t1d:.3f}, {_ns}×{_nt} sample):")
print(f"    min={_c_min:.1f}  median={_c_med:.1f}  p99={_c_p99:.1f}  max={_c_max:.1f}")
_eps_rec_lo = _c_med / 10.0
_eps_rec_hi = _c_med / 2.0
print(f"\n  Recommended epsilon range : {_eps_rec_lo:.1f} – {_eps_rec_hi:.1f}")
print(f"  Current EPSILON           : {EPSILON}")
if EPSILON < _eps_rec_lo:
    print(f"  !! WARNING: EPSILON={EPSILON} may be too small.  "
          f"Set ≥ {_eps_rec_lo:.0f} to avoid overflow.")
else:
    print(f"  ✓  EPSILON looks acceptable for this cost distribution.")

_check_eps = sorted(set([EPSILON, _eps_rec_lo, _eps_rec_hi, 0.05, 0.5]))
print(f"\n  K-table:")
for _e in _check_eps:
    _km = np.exp(-_c_med / _e)
    _kn = np.exp(-_c_max / _e)
    _st = ("OK" if _kn > 1e-15 and _km > 1e-4
           else "MARGINAL" if _kn > 1e-30 else "WILL OVERFLOW")
    _tag = "  ◄ current" if _e == EPSILON else ""
    print(f"    ε={_e:<8.2f}  median_cost/ε={_c_med/_e:7.2f}  "
          f"K_median={_km:.2e}  K_min={_kn:.2e}  [{_st}]{_tag}")

# =============================================================================
# 4.  OPTIONAL SUBSAMPLE
# =============================================================================
if DEBUG_SUBSAMPLE_PER_TIMEPOINT is not None:
    print(f"\n{_SEP}")
    print(f"STEP 4  —  Subsample to {DEBUG_SUBSAMPLE_PER_TIMEPOINT} cells/timepoint")
    print(_SEP)
    _rng_sub = np.random.default_rng(0)
    _keep    = []
    for _d in _sorted_days:
        _idx = np.where(adata_full.obs["day"].values == _d)[0]
        _n   = min(DEBUG_SUBSAMPLE_PER_TIMEPOINT, len(_idx))
        _keep.extend(_rng_sub.choice(_idx, _n, replace=False).tolist())
    _keep.sort()
    adata = adata_full[_keep].copy()
    print(f"  After subsampling : {adata.n_obs:,} cells")
else:
    adata = adata_full

_days_ot   = sorted(adata.obs["day"].unique())
_all_pairs = [(float(_days_ot[i]), float(_days_ot[i+1]))
              for i in range(len(_days_ot)-1)]
if DEBUG_ONLY_FIRST_N_PAIRS is not None:
    _days_keep = sorted(set(
        d for t1, t2 in _all_pairs[:DEBUG_ONLY_FIRST_N_PAIRS] for d in (t1, t2)
    ))
    adata      = adata[adata.obs["day"].isin(_days_keep)].copy()
    _all_pairs = _all_pairs[:DEBUG_ONLY_FIRST_N_PAIRS]
    print(f"  Restricted to {len(_days_keep)} timepoints / {len(_all_pairs)} pair(s)")

print(f"\n  Cells for OT  : {adata.n_obs:,}")
print(f"  Pairs to run  : {len(_all_pairs)}")

# =============================================================================
# 5.  PREPARE WOT CLI INPUT FILES
# =============================================================================
print(f"\n{_SEP}")
print("STEP 5  —  Prepare WOT CLI input files")
print(_SEP)

# ── a) ExprMatrix.h5ad ────────────────────────────────────────────────────────
adata_wot = adata.copy()
if "day" in adata_wot.obs.columns:
    adata_wot.obs = adata_wot.obs.rename(columns={"day": "day_in_stage"})
    print("  Renamed obs['day'] → 'day_in_stage'  (avoids WOT internal collision)")

if FORCE_DENSE_FLOAT64:
    _Xd = (adata_wot.X.toarray().astype(np.float64) if sp.issparse(adata_wot.X)
           else np.asarray(adata_wot.X, dtype=np.float64))
    adata_wot = ad.AnnData(X=_Xd, obs=adata_wot.obs.copy(), var=adata_wot.var.copy())
    print(f"  X converted to dense float64  (FORCE_DENSE_FLOAT64=True)")
else:
    _x_type = "sparse" if sp.issparse(adata_wot.X) else "dense"
    _x_dt   = adata_wot.X.data.dtype if sp.issparse(adata_wot.X) else adata_wot.X.dtype
    print(f"  X kept as {_x_type} {_x_dt}")

MATRIX_FILE_PATH = WOT_INPUTS / "ExprMatrix.h5ad"
adata_wot.write_h5ad(str(MATRIX_FILE_PATH))
print(f"  ExprMatrix.h5ad  → {MATRIX_FILE_PATH.name}  "
      f"({adata_wot.n_obs:,} cells × {adata_wot.n_vars:,} genes)")

# ── b) cell_days.txt ──────────────────────────────────────────────────────────
CELL_DAYS_PATH = WOT_INPUTS / "cell_days.txt"
_cd = pd.DataFrame({
    "id":  adata.obs_names.astype(str),
    "day": adata.obs["day"].values.astype(float),
})
_cd.to_csv(str(CELL_DAYS_PATH), sep="\t", index=False)
print(f"  cell_days.txt    → {CELL_DAYS_PATH.name}  "
      f"(timepoints: {sorted(_cd['day'].unique())})")

# ── c) growth_rates.txt ───────────────────────────────────────────────────────
GROWTH_RATES_PATH = WOT_INPUTS / "growth_rates.txt"
_gr = pd.DataFrame({
    "id":          adata.obs_names.astype(str),
    "growth_rate": adata.obs["growth_rate"].values.astype(float),
})
_gr.to_csv(str(GROWTH_RATES_PATH), sep="\t", index=False)
print(f"  growth_rates.txt → {GROWTH_RATES_PATH.name}  "
      f"(all zero: {(_gr['growth_rate'] == 0).all()})")

# =============================================================================
# 6.  RUN WOT CLI  via subprocess
# =============================================================================
print(f"\n{_SEP}")
print("STEP 6  —  Run WOT CLI  (wot optimal_transport)")
print(_SEP)

TMAPS_PREFIX = str(TMAPS_DIR / "tmaps")

def _find_wot_exe() -> list[str]:
    """Return the command prefix to invoke the WOT CLI."""
    import shutil, sys as _sys
    if shutil.which("wot"):
        return ["wot"]
    try:
        _r = subprocess.run(
            [_sys.executable, "-m", "wot", "--help"],
            capture_output=True, timeout=10
        )
        if _r.returncode == 0 or b"optimal_transport" in _r.stdout + _r.stderr:
            return [_sys.executable, "-m", "wot"]
    except Exception:
        pass
    try:
        import wot as _wot_pkg
        _pkg_dir = Path(_wot_pkg.__file__).parent
        for _c in [_pkg_dir.parent / "Scripts" / "wot",
                   _pkg_dir.parent / "Scripts" / "wot.exe",
                   _pkg_dir.parent / "bin" / "wot"]:
            if _c.exists():
                return [str(_c)]
    except Exception:
        pass
    raise FileNotFoundError(
        "Cannot find the 'wot' CLI executable.\n"
        "  conda activate traj_env && pip install wot\n"
        "Then verify with:  wot --help"
    )

_wot_cmd = _find_wot_exe()
print(f"  WOT executable : {' '.join(_wot_cmd)}")

_cmd = (
    _wot_cmd + [
        "optimal_transport",
        "--matrix",       str(MATRIX_FILE_PATH),
        "--cell_days",    str(CELL_DAYS_PATH),
        "--epsilon",      str(EPSILON),
        "--lambda1",      str(LAMBDA1),
        "--lambda2",      str(LAMBDA2),
        "--local_pca",    str(LOCAL_PCA),
        "--growth_iters", "1",
        "--out",          TMAPS_PREFIX,
        "--verbose",
    ]
)
if USE_GROWTH_RATES:
    _cmd += ["--cell_growth_rates", str(GROWTH_RATES_PATH)]

print(f"\n  Command:\n    {' '.join(_cmd)}")
print(f"\n  Running … (may take several minutes for large datasets)")

_t_cli_start = _time.time()
try:
    _proc = subprocess.run(
        _cmd, capture_output=True, text=True,
        cwd=str(TMAPS_DIR), timeout=3600,
    )
    _cli_elapsed = _time.time() - _t_cli_start
    print(f"\n  Return code : {_proc.returncode}  ({_cli_elapsed:.1f} s)")

    if _proc.stdout.strip():
        print(f"\n  --- WOT stdout ---")
        for _ln in _proc.stdout.strip().splitlines():
            print(f"  {_ln}")
    if _proc.stderr.strip():
        print(f"\n  --- WOT stderr ---")
        for _ln in _proc.stderr.strip().splitlines():
            print(f"  {_ln}")
    if _proc.returncode != 0:
        print(f"\n  !! WOT CLI exited with non-zero return code {_proc.returncode}.")

except subprocess.TimeoutExpired:
    print("  !! Timeout: WOT CLI did not complete within 1 hour.")
    sys.exit(1)
except FileNotFoundError as _e:
    print(f"  !! {_e}")
    sys.exit(1)

# =============================================================================
# 7.  DISCOVER AND VALIDATE OUTPUT TRANSPORT MAPS
# =============================================================================
print(f"\n{_SEP}")
print("STEP 7  —  Discover output transport maps")
print(_SEP)

def _parse_tmap_pair(fp: Path):
    """Extract (t0, t1) from tmaps_<t0>_<t1>.h5ad filename."""
    m = re.search(r"tmaps_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)\.h5ad$", fp.name)
    return (float(m.group(1)), float(m.group(2))) if m else None

_tmap_files = sorted(TMAPS_DIR.glob("tmaps_*.h5ad"))
print(f"  Found {len(_tmap_files)} tmap file(s)")

successful_tmaps: dict = {}   # (t0, t1) → Path
for _fp in _tmap_files:
    _pr = _parse_tmap_pair(_fp)
    if _pr is None:
        print(f"  [SKIP] Cannot parse pair from: {_fp.name}")
        continue
    _t0p, _t1p = _pr
    successful_tmaps[(_t0p, _t1p)] = _fp
    _tm = sc.read_h5ad(str(_fp))
    print(f"  ✓  T({_t0p:.3f}→{_t1p:.3f})  obs={_tm.n_obs}  var={_tm.n_vars}  "
          f"→ {_fp.name}")

if not successful_tmaps:
    print(f"\n  !! No tmap files found.  Possible causes:")
    print(f"     • WOT CLI returned non-zero (see Step 6 stderr)")
    print(f"     • EPSILON={EPSILON} caused overflow — try ≥ {_eps_rec_lo:.0f}")
    print(f"     • Output written elsewhere — check {TMAPS_DIR}")
    RUN_FATE_COMPUTATION = False
    RUN_PLOTTING         = False
else:
    print(f"\n  {len(successful_tmaps)} transport map(s) ready for fate computation.")

# =============================================================================
# 8.  FATE COMPUTATION  (multi-fate biological design, WOT Python API)
# =============================================================================
#
# Framework (original WOT paper design, Schiebinger et al. 2019)
# ---------------------------------------------------------------
# Define K ≥ 2 biologically distinct terminal fate sets; compute per-cell
# fate probabilities using the WOT TransportMapModel Python API, which
# handles fate sets at different timepoints (e.g. hCiPSC at day30,
# StageIII_late at day21) correctly via multi-step backward propagation.
#
# Fate sets (stage annotation, independent of OCLR / trajectory models):
#   fate_hcipsc        — ALL cells with stage_std == "hCiPSC"
#   fate_stageIII_late — StageIII cells at their latest observed timepoint
#
# Output columns on adata_full.obs:
#   p_hcipsc        — fate probability toward hCiPSC
#   p_stageIII_late — fate probability toward StageIII_late
#   oclr_score      — OCLR stemness score (annotation only, NaN for non-hCiPSC)
#   oclr_ips_subset_mask — OCLR top-10% binary flag (annotation only)
#
print(f"\n{_SEP}")
print("STEP 8  —  Fate computation  (multi-fate biological design, WOT Python API)")
print(_SEP)

# Tracking variables initialised here; updated inside the if-block below.
_n_hcipsc              = 0
_n_stageIII_late       = 0
_T_stageIII_latest     = np.nan
_mask_cipsc            = None
_p_hcipsc_stats        = {}
_p_stageIII_stats      = {}

if not RUN_FATE_COMPUTATION or not successful_tmaps:
    print("  Skipped  (no transport maps available)")
else:

    # ─────────────────────────────────────────────────────────────────────────
    # 8a.  Standardize stage labels → obs["stage_std"]
    # ─────────────────────────────────────────────────────────────────────────
    print(f"\n{_SEP2}")
    print("  Step 8a — Standardize stage labels → obs['stage_std']")

    adata_full.obs["stage_std"] = (
        adata_full.obs["stage"].apply(normalize_stage).astype(str)
    )
    _stage_counts = adata_full.obs["stage_std"].value_counts().sort_index()
    print(f"  Standardized stage distribution:")
    for _st, _cnt in _stage_counts.items():
        print(f"    {_st:<15} : {_cnt:>6,} cells")

    _mask_cipsc = adata_full.obs["stage_std"] == "hCiPSC"
    _n_day30_total = int(_mask_cipsc.sum())
    print(f"\n  day30/hCiPSC cells total : {_n_day30_total:,}")

    if _n_day30_total == 0:
        print("  !! FATAL: No cells with stage_std == 'hCiPSC'.")
        print(f"     Observed values: {sorted(adata_full.obs['stage_std'].unique())}")
        RUN_FATE_COMPUTATION = False

    # ─────────────────────────────────────────────────────────────────────────
    # 8b.  Define multi-fate biological terminal sets
    # ─────────────────────────────────────────────────────────────────────────
    # fate_hcipsc        — ALL hCiPSC cells (fully reprogrammed iPSCs)
    # fate_stageIII_late — StageIII cells at their latest observed timepoint
    #
    # OCLR scores are attached as obs annotations for external evaluation only.
    if RUN_FATE_COMPUTATION:
        print(f"\n{_SEP2}")
        print("  Step 8b — Define multi-fate biological terminal sets")

        # ── fate_hcipsc: all hCiPSC cells ────────────────────────────────────
        _mask_cipsc_bool   = _mask_cipsc.values if hasattr(_mask_cipsc, "values") else _mask_cipsc
        _n_hcipsc          = int(_mask_cipsc_bool.sum())
        adata_full.obs["fate_hcipsc"] = _mask_cipsc_bool.astype(int)

        print(f"\n  fate_hcipsc (stage_std == 'hCiPSC') : {_n_hcipsc:,} cells")
        if _n_hcipsc < MIN_FATE_CELLS:
            print(f"  !! FATAL: Only {_n_hcipsc} hCiPSC cells (need ≥ {MIN_FATE_CELLS}).")
            RUN_FATE_COMPUTATION = False

    if RUN_FATE_COMPUTATION:
        # ── fate_stageIII_late: StageIII at latest observed timepoint ─────────
        _mask_stageIII_all = (adata_full.obs["stage_std"] == "StageIII").values
        if _mask_stageIII_all.any():
            _T_stageIII_latest  = float(
                adata_full.obs.loc[_mask_stageIII_all, "day"].max()
            )
            _mask_stageIII_late = (
                _mask_stageIII_all &
                (adata_full.obs["day"].values == _T_stageIII_latest)
            )
            _n_stageIII_late = int(_mask_stageIII_late.sum())
        else:
            _mask_stageIII_late = np.zeros(adata_full.n_obs, dtype=bool)
            _n_stageIII_late    = 0
            _T_stageIII_latest  = np.nan

        adata_full.obs["fate_stageIII_late"] = _mask_stageIII_late.astype(int)

        print(f"  fate_stageIII_late (StageIII @ day {_T_stageIII_latest}) : "
              f"{_n_stageIII_late:,} cells")
        if _n_stageIII_late < MIN_FATE_CELLS:
            print(f"  !! FATAL: Only {_n_stageIII_late} StageIII_late cells "
                  f"(need ≥ {MIN_FATE_CELLS}).")
            print(f"     Check that 'StageIII' cells are present in obs['stage_std'].")
            print(f"     Observed values: {sorted(adata_full.obs['stage_std'].unique())}")
            RUN_FATE_COMPUTATION = False

    if RUN_FATE_COMPUTATION:
        # ── OCLR scores: load as obs annotation only ──────────────────────────
        print(f"\n{_SEP2}")
        print("  Step 8b (cont.) — Attach OCLR scores as obs annotation")
        print(f"  (Annotation only — NOT used to define fate sets)")

        try:
            _oclr_score_series, _oclr_subset_mask = load_oclr_endpoint(adata_full)
            adata_full.obs["oclr_score"]           = _oclr_score_series.values
            adata_full.obs["oclr_ips_subset_mask"] = _oclr_subset_mask.astype(int).values
            _n_oclr_annotated = int(np.isfinite(_oclr_score_series.values).sum())
            print(f"  OCLR scores attached: {_n_oclr_annotated:,} cells have finite scores")
        except FileNotFoundError as _e:
            print(f"  [WARN] OCLR score file not found — skipping annotation.")
            print(f"    {_e}")
            adata_full.obs["oclr_score"]           = np.nan
            adata_full.obs["oclr_ips_subset_mask"] = 0

    # ─────────────────────────────────────────────────────────────────────────
    # 8c.  Build cell_set_matrix for WOT Python API
    # ─────────────────────────────────────────────────────────────────────────
    if RUN_FATE_COMPUTATION:
        print(f"\n{_SEP2}")
        print("  Step 8c — Build cell_set_matrix (cells × fates, binary)")

        _fate_names = ["fate_hcipsc", "fate_stageIII_late"]
        _cell_bcs   = adata_full.obs_names.astype(str).tolist()

        _cell_set_df = pd.DataFrame(
            0, index=_cell_bcs, columns=_fate_names, dtype=np.float32
        )
        _cell_set_df.loc[
            adata_full.obs_names[_mask_cipsc_bool].astype(str), "fate_hcipsc"
        ] = 1.0
        if _n_stageIII_late > 0:
            _cell_set_df.loc[
                adata_full.obs_names[_mask_stageIII_late].astype(str),
                "fate_stageIII_late"
            ] = 1.0

        _n_hcipsc_in_csm = int((_cell_set_df["fate_hcipsc"] > 0).sum())
        _n_s3_in_csm     = int((_cell_set_df["fate_stageIII_late"] > 0).sum())
        print(f"  cell_set_matrix : {len(_cell_bcs):,} cells × {len(_fate_names)} fates")
        print(f"    fate_hcipsc entries        : {_n_hcipsc_in_csm:,}")
        print(f"    fate_stageIII_late entries : {_n_s3_in_csm:,}")

        _cell_set_adata = ad.AnnData(
            X=_cell_set_df.values,
            obs=pd.DataFrame(index=_cell_set_df.index),
            var=pd.DataFrame(index=_cell_set_df.columns),
        )

    # ─────────────────────────────────────────────────────────────────────────
    # 8d.  Fate computation — WOT Python API
    # ─────────────────────────────────────────────────────────────────────────
    #
    # wot.tmap.TransportMapModel.from_directory() loads the transport maps
    # produced by the WOT CLI (step 6) as a linked chain.
    #
    # .get_fate_probabilities(cell_set_matrix) propagates fate indicators
    # backward through the chain, handling fate sets at different timepoints
    # correctly (e.g., fate_hcipsc at day30, fate_stageIII_late at day21).
    #
    # Returns an AnnData: obs = all project cells, var = fate names,
    # X[i, k] = fate probability of cell i toward fate k.
    # ─────────────────────────────────────────────────────────────────────────
    if RUN_FATE_COMPUTATION:
        print(f"\n{_SEP2}")
        print("  Step 8d — Fate computation  (WOT Python API)")

        # Load transport map model from directory
        try:
            _tmap_model = wot.tmap.TransportMapModel.from_directory(
                str(TMAPS_DIR), prefix="tmaps"
            )
            print(f"  TransportMapModel loaded from : {TMAPS_DIR.name}/")
        except AttributeError:
            # Older WOT versions may have a different module path
            import wot.tmap as _wot_tmap
            _tmap_model = _wot_tmap.TransportMapModel.from_directory(
                str(TMAPS_DIR), prefix="tmaps"
            )
            print(f"  TransportMapModel (fallback import) loaded from : {TMAPS_DIR.name}/")

        # Compute fate probabilities
        print(f"  Computing get_fate_probabilities() ...")
        _fate_ds = _tmap_model.get_fate_probabilities(_cell_set_adata)
        # _fate_ds: AnnData with obs=cells, var=fates, X=fate_prob matrix

        print(f"  Fate probability AnnData : {_fate_ds.n_obs:,} cells × {_fate_ds.n_vars} fates")
        print(f"  Fate names               : {list(_fate_ds.var_names)}")

        # ── Align fate probabilities back to adata_full ───────────────────────
        _fate_df = pd.DataFrame(
            _fate_ds.X if not sp.issparse(_fate_ds.X) else _fate_ds.X.toarray(),
            index=_fate_ds.obs_names.astype(str),
            columns=_fate_ds.var_names.astype(str),
        )

        _p_hcipsc_vec     = np.full(adata_full.n_obs, np.nan, dtype=np.float64)
        _p_stageIII_vec   = np.full(adata_full.n_obs, np.nan, dtype=np.float64)

        for _i, _bc in enumerate(adata_full.obs_names.astype(str)):
            if _bc in _fate_df.index:
                if "fate_hcipsc" in _fate_df.columns:
                    _p_hcipsc_vec[_i] = float(_fate_df.loc[_bc, "fate_hcipsc"])
                if "fate_stageIII_late" in _fate_df.columns:
                    _p_stageIII_vec[_i] = float(_fate_df.loc[_bc, "fate_stageIII_late"])

        adata_full.obs["p_hcipsc"]        = _p_hcipsc_vec
        adata_full.obs["p_stageIII_late"] = _p_stageIII_vec

        # ── Sanity checks ─────────────────────────────────────────────────────
        _n_finite_hcipsc  = int(np.isfinite(_p_hcipsc_vec).sum())
        _n_finite_s3      = int(np.isfinite(_p_stageIII_vec).sum())
        _early_mask_bool  = ~_mask_cipsc_bool & np.isfinite(_p_hcipsc_vec)
        _p_early_hcipsc   = _p_hcipsc_vec[_early_mask_bool]

        print(f"\n  p_hcipsc finite values        : {_n_finite_hcipsc:,}  "
              f"({_n_finite_hcipsc / adata_full.n_obs:.1%})")
        print(f"  p_stageIII_late finite values : {_n_finite_s3:,}  "
              f"({_n_finite_s3 / adata_full.n_obs:.1%})")

        if len(_p_early_hcipsc) > 0:
            _p_hcipsc_stats = {
                "min":    float(np.nanmin(_p_early_hcipsc)),
                "median": float(np.nanmedian(_p_early_hcipsc)),
                "max":    float(np.nanmax(_p_early_hcipsc)),
            }
            print(f"  p_hcipsc (non-terminal cells): "
                  f"min={_p_hcipsc_stats['min']:.4f}  "
                  f"median={_p_hcipsc_stats['median']:.4f}  "
                  f"max={_p_hcipsc_stats['max']:.4f}")

            if _p_early_hcipsc.min() > 0.99:
                print("  !! WARNING: p_hcipsc ≈ 1.0 for all non-terminal cells — "
                      "degenerate!  Check that fate_stageIII_late is non-empty.")
            elif _p_early_hcipsc.max() < 0.01:
                print("  !! WARNING: p_hcipsc ≈ 0.0 for all non-terminal cells — "
                      "inverse collapse.  Check barcode matching.")
            else:
                print("  ✓  p_hcipsc shows variation across non-terminal cells — "
                      "non-degenerate.")

    # ── Save adata_full with all new columns ──────────────────────────────────
    if RUN_FATE_COMPUTATION:
        _out_h5ad = METRICS_DIR / f"{TIMESTAMP}_GSE230659_wot_fates.h5ad"
        adata_full.write_h5ad(str(_out_h5ad))
        print(f"\n  Saved adata_full with fate columns → {_out_h5ad.name}")
        print(f"    obs columns added: p_hcipsc, p_stageIII_late, "
              f"fate_hcipsc, fate_stageIII_late, oclr_score, oclr_ips_subset_mask")

# =============================================================================
# 9.  PLOTTING
# =============================================================================
print(f"\n{_SEP}")
print("STEP 9  —  Plotting")
print(_SEP)

_umap_recomputed = False

if not RUN_PLOTTING:
    print("  Skipped  (RUN_PLOTTING=False)")
elif "p_hcipsc" not in adata_full.obs.columns:
    print("  Skipped  (no fate probabilities available)")
else:

    # ── Ensure UMAP embedding exists ──────────────────────────────────────────
    try:
        if "X_umap" in adata_full.obsm:
            print("  UMAP : reusing existing X_umap from adata_full.obsm")
            _umap_recomputed = False
        else:
            if "X_pca" in adata_full.obsm:
                print("  UMAP : computing from X_pca (n_neighbors=15) …")
                sc.pp.neighbors(adata_full, use_rep="X_pca", n_neighbors=15)
            else:
                print("  UMAP : no X_pca available; running PCA on X first …")
                sc.pp.pca(adata_full, n_comps=30)
                sc.pp.neighbors(adata_full, use_rep="X_pca", n_neighbors=15)
            sc.tl.umap(adata_full)
            _umap_recomputed = True
            print("  UMAP computation complete.")
    except Exception as _ue:
        print(f"  UMAP computation failed: {_ue}")

    # ── Plot A: continuous day vs p_hcipsc scatter ────────────────────────────
    try:
        _fig, _ax = plt.subplots(figsize=(8, 5))
        _sc_data  = adata_full.obs[["day", "p_hcipsc", "stage_std"]].dropna(
            subset=["p_hcipsc"]
        )
        for _st, _col in STAGE_PALETTE.items():
            _mask_st = _sc_data["stage_std"].apply(
                lambda x: normalize_stage(x) == _st
            )
            _sub = _sc_data[_mask_st]
            if len(_sub):
                _ax.scatter(_sub["day"], _sub["p_hcipsc"],
                            c=_col, s=2, alpha=0.3, label=_st, rasterized=True)
        _ax.set_xlabel("Day")
        _ax.set_ylabel("p(hCiPSC fate)")
        _ax.set_title("WOT fate probability — transition toward hCiPSC  (multi-fate)")
        _ax.legend(markerscale=4, fontsize=8)
        _fig.tight_layout()
        _fig_path_A = FIGURES_DIR / f"{TIMESTAMP}_fate_prob_day.png"
        _fig.savefig(str(_fig_path_A), dpi=150, bbox_inches="tight")
        plt.close(_fig)
        print(f"  Plot A (day vs p_hcipsc) → {_fig_path_A.name}")
    except Exception as _pe:
        print(f"  Plot A failed: {_pe}")

    # ── Plot B: stage-colored UMAP with p_hcipsc continuous ──────────────────
    if "X_umap" in adata_full.obsm:
        try:
            _u     = adata_full.obsm["X_umap"]
            _stage = adata_full.obs["stage_std"].apply(normalize_stage).values

            _fig_B, _ax_B = plt.subplots(figsize=(8, 7))
            _draw_order = ["Other", "StageI", "StageII", "StageIII", "hCiPSC"]
            for _lbl in _draw_order:
                _sel = _stage == _lbl
                if not _sel.any():
                    continue
                _col = STAGE_PALETTE.get(_lbl, "#BDBDBD")
                _alpha = 0.15 if _lbl == "Other" else 0.55 if "Stage" in _lbl else 0.80
                _size  = 3    if _lbl == "Other" else 5    if "Stage" in _lbl else 7
                _ax_B.scatter(
                    _u[_sel, 0], _u[_sel, 1],
                    c=_col, s=_size, alpha=_alpha,
                    label=_lbl, rasterized=True, linewidths=0,
                )
            _ax_B.set_xlabel("UMAP 1"); _ax_B.set_ylabel("UMAP 2")
            _ax_B.set_title("WOT — stage annotation  (hCiPSC = fate target)")
            _ax_B.legend(markerscale=3, fontsize=9, framealpha=0.8, loc="best")
            _ax_B.set_xticks([]); _ax_B.set_yticks([])
            _fig_B.tight_layout()
            _fig_path_B = FIGURES_DIR / f"{TIMESTAMP}_umap_stage.png"
            _fig_B.savefig(str(_fig_path_B), dpi=150, bbox_inches="tight")
            plt.close(_fig_B)
            print(f"  Plot B (UMAP stage) → {_fig_path_B.name}")
        except Exception as _pe:
            print(f"  Plot B failed: {_pe}")

    # ── Plot C: continuous p_hcipsc on UMAP ──────────────────────────────────
    if "X_umap" in adata_full.obsm:
        try:
            _fig_C, _ax_C = plt.subplots(figsize=(7, 6))
            _p_all  = adata_full.obs["p_hcipsc"].values
            _finite = np.isfinite(_p_all)
            _sc_C   = _ax_C.scatter(
                adata_full.obsm["X_umap"][_finite, 0],
                adata_full.obsm["X_umap"][_finite, 1],
                c=_p_all[_finite], cmap="RdYlBu_r",
                s=3, alpha=0.5, vmin=0, vmax=1, rasterized=True, linewidths=0,
            )
            plt.colorbar(_sc_C, ax=_ax_C, label="p(hCiPSC fate)")
            _ax_C.set_xlabel("UMAP 1"); _ax_C.set_ylabel("UMAP 2")
            _ax_C.set_title("WOT p(hCiPSC) — continuous  (RdYlBu_r)")
            _ax_C.set_xticks([]); _ax_C.set_yticks([])
            _fig_C.tight_layout()
            _fig_path_C = FIGURES_DIR / f"{TIMESTAMP}_umap_p_hcipsc.png"
            _fig_C.savefig(str(_fig_path_C), dpi=150, bbox_inches="tight")
            plt.close(_fig_C)
            print(f"  Plot C (UMAP p_hcipsc) → {_fig_path_C.name}")
        except Exception as _pe:
            print(f"  Plot C failed: {_pe}")

# =============================================================================
# 10.  FINAL SUMMARY
# =============================================================================
_elapsed = _time.time() - T0
print(f"\n{_SEP}")
print("STEP 10 —  Final summary")
print(_SEP)
print(f"  Total elapsed              : {_elapsed:.1f} s")
print(f"  adata.X audit              : {_verdict_short}")
print(f"  EPSILON used               : {EPSILON}  "
      f"(recommended {_eps_rec_lo:.1f}–{_eps_rec_hi:.1f})")
print(f"  Transport maps found       : {len(successful_tmaps)}")
print(f"  USE_GROWTH_RATES           : {USE_GROWTH_RATES}")
print()
print(f"  ── Fate computation  (multi-fate biological design, WOT Python API) ───────")
print(f"  fate_hcipsc                : {_n_hcipsc:,} cells  "
      f"(ALL hCiPSC; biological fate target)")
print(f"  fate_stageIII_late         : {_n_stageIII_late:,} cells  "
      f"(StageIII @ day {_T_stageIII_latest}; competing biological fate)")
if _p_hcipsc_stats:
    print(f"  p_hcipsc (non-terminal)    : "
          f"min={_p_hcipsc_stats['min']:.4f}  "
          f"median={_p_hcipsc_stats['median']:.4f}  "
          f"max={_p_hcipsc_stats['max']:.4f}")
print(f"  OCLR annotation            : obs['oclr_score'] + obs['oclr_ips_subset_mask']")
print(f"    (external biological anchor for evaluation in script 05 only)")
print()
print(f"  UMAP                       : "
      + ("newly computed" if _umap_recomputed else "reused from obsm"))

if not successful_tmaps:
    print(f"\n  !! All OT pairs failed.  Suggested next steps:")
    print(f"     1. Check Step 6 stderr for specific WOT error messages.")
    print(f"     2. EPSILON={EPSILON} — if overflow, set ≥ {_eps_rec_lo:.0f}.")
    print(f"     3. If adata.X is not log1p-normalised, fix upstream.")
    print(f"     4. Try LOCAL_PCA=5 or FORCE_DENSE_FLOAT64=True.")
elif not RUN_FATE_COMPUTATION:
    print(f"\n  !! WOT OT succeeded but fate computation was skipped.")
    print(f"     Check warnings above  "
          f"(fate_hcipsc: {_n_hcipsc}, fate_stageIII_late: {_n_stageIII_late}).")
else:
    print(f"\n  ✓  WOT CLI + multi-fate computation complete.")
    print(f"     Key outputs (adata_full.obs):")
    print(f"       p_hcipsc            : fate probability toward hCiPSC")
    print(f"       p_stageIII_late     : fate probability toward StageIII_late")
    print(f"       fate_hcipsc         : binary fate membership (0/1)")
    print(f"       fate_stageIII_late  : binary fate membership (0/1)")
    print(f"       oclr_score          : OCLR stemness score (annotation, NaN for non-hCiPSC)")
    print(f"       oclr_ips_subset_mask: OCLR top-10% flag (annotation)")
    _h5ad_out = METRICS_DIR / f"{TIMESTAMP}_GSE230659_wot_fates.h5ad"
    print(f"     h5ad saved to : {_h5ad_out.name}")
