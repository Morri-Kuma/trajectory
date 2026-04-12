#!/usr/bin/env python3
"""
02c_compute_transport.py  [v14 — OCLR-anchored dual-endpoint terminal design]
==========================================================================================
Project : Comparative Study of Trajectory Inference Models for Chemical iPSC
          Reprogramming
Dataset : GSE230659 (human, Liuyang et al. 2023 Cell Stem Cell)
Step    : 2c — WOT transport-map computation, fate probabilities, trajectory trends
Env     : conda activate traj_env

OCLR-anchored dual-endpoint design (v14)
-----------------------------------------
  Two competing fate sets are defined within the hCiPSC endpoint sample using the
  OCLR stemness labels produced by 02b2_define_oclr_endpoint.py:

    fate_oclr_success  — top-10% OCLR score within hCiPSC  (high-stemness)
    fate_oclr_failure  — bottom-10% OCLR score within hCiPSC  (low-stemness)

  Both sets are at the same timepoint (day30 hCiPSC sample), satisfying the WOT
  fates() single-timepoint constraint.  The middle 80% (ambiguous) cells are NOT
  used as a fate target.

  Backward-propagated fate probabilities:
    p_oclr_success(i)  = P(cell i trajectory → high-stemness hCiPSC)
    p_oclr_failure(i)  = P(cell i trajectory → low-stemness hCiPSC)
    p_oclr_margin(i)   = p_oclr_success(i) − p_oclr_failure(i)  [main benchmark metric]

  Legacy alias: p_hcipsc = p_oclr_success  (kept for backward compatibility).

  Fate computation uses the WOT Python API:
    wot.tmap.TransportMapModel.from_directory()  +  .fates(populations)
  where populations are built with population_from_cell_sets(cell_sets, at_time).

  Run 02b2_define_oclr_endpoint.py BEFORE this script.

Version history
---------------
  v9 : Two-fate design introduced.
  v10: 6-class UMAP label.
  v11: CLI OT pipeline unchanged.
  v12: Endpoint replaced from marker genes → OCLR stemness scores.
  v13: Multi-fate biological design (hCiPSC + StageIII_late).
  v14: OCLR-anchored dual-endpoint design.  fate_oclr_success / fate_oclr_failure
       replace hCiPSC + StageIII_late.  Uses fates() API.
       p_oclr_margin added as primary benchmark metric.

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
EPSILON:   float = 80.0
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
# FATE COMPUTATION PARAMETERS  (OCLR-anchored dual-endpoint design)
# ══════════════════════════════════════════════════════════════════════════════
#
# TERMINAL FATE SETS (OCLR-anchored, loaded from 02b2 output)
# ------------------------------------------------------------
# Both fate sets are within the hCiPSC endpoint sample at day30.
# This satisfies WOT fates() single-timepoint constraint.
#
#   fate_oclr_success  — top-10% OCLR score within hCiPSC (is_oclr_success_endpoint == 1)
#   fate_oclr_failure  — bottom-10% OCLR score within hCiPSC (is_oclr_failure_endpoint == 1)
#
# WHY THIS IS NON-DEGENERATE
# ---------------------------
# Both fate sets are subsets of the same sample, so they partition a biologically
# meaningful within-sample axis (stemness level).  The backward-propagated
# probabilities are:
#   p_oclr_success(i)  = P(cell i trajectory → high-stemness hCiPSC endpoint)
#   p_oclr_failure(i)  = P(cell i trajectory → low-stemness hCiPSC endpoint)
#   p_oclr_margin(i)   = p_oclr_success − p_oclr_failure  ← primary benchmark metric
#
# Minimum number of cells required in each fate set.
# If either falls below this the script stops with a clear error.
MIN_FATE_CELLS: int = 20

# ── OCLR endpoint sample ID ──────────────────────────────────────────────────
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
# hCiPSC   — fully reprogrammed cells (contains OCLR success/failure endpoints)
# StageIII — intermediate reprogramming stage
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

print(f"[{TIMESTAMP}]  WOT Transport Map Computation  (v14 — OCLR dual-endpoint design)")
print(f"  wot version     : {_wot_ver}")
print(f"  Input AnnData   : {H5AD_PATH.name}")
print(f"  Output dir      : {TMAPS_DIR}")
print(f"\n  WOT CLI parameters")
print(f"    EPSILON          : {EPSILON}")
print(f"    LAMBDA1/LAMBDA2  : {LAMBDA1} / {LAMBDA2}")
print(f"    LOCAL_PCA        : {LOCAL_PCA}")
print(f"    USE_GROWTH_RATES : {USE_GROWTH_RATES}")
print(f"\n  Fate computation design (OCLR dual-endpoint, WOT fates() API)")
print(f"    fate_oclr_success  : top-10% OCLR within GSM7230012_hCiPSCs-0618")
print(f"    fate_oclr_failure  : bottom-10% OCLR within GSM7230012_hCiPSCs-0618")
print(f"    p_oclr_success     : P(cell → high-stemness hCiPSC endpoint)")
print(f"    p_oclr_failure     : P(cell → low-stemness hCiPSC endpoint)")
print(f"    p_oclr_margin      : p_oclr_success − p_oclr_failure  [PRIMARY metric]")
print(f"    p_hcipsc           : legacy alias = p_oclr_success")
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
# HELPER: load OCLR dual endpoints from 02b2 output
# =============================================================================
def load_oclr_endpoints(adata):
    """
    Load the OCLR dual-endpoint labels produced by 02b2_define_oclr_endpoint.py.

    Searches PROCESSED_DIR for the most recent *_oclr_endpoint_labels_hcipsc.tsv
    (the authoritative per-cell file from 02b2 v14+).  Falls back to loading
    separate *_oclr_success_endpoint_barcodes.txt / *_oclr_failure_endpoint_barcodes.txt
    if the tsv is not found.

    Returns
    -------
    oclr_score_series    : pd.Series[float] — OCLR score per cell (NaN for non-hCiPSC)
    success_mask         : pd.Series[bool]  — True = fate_oclr_success endpoint
    failure_mask         : pd.Series[bool]  — True = fate_oclr_failure endpoint
    endpoint_label_series: pd.Series[str]   — 'success'|'failure'|'ambiguous'|'unscored'|''
    """
    _proj_barcodes = adata.obs_names.astype(str).tolist()
    _proj_set      = set(_proj_barcodes)

    # ── Try authoritative TSV first ───────────────────────────────────────────
    _label_files = sorted(PROCESSED_DIR.glob("*_oclr_endpoint_labels_hcipsc.tsv"))
    if _label_files:
        _lf = _label_files[-1]
        print(f"  OCLR endpoint labels file : {_lf.name}")
        _ldf = pd.read_csv(str(_lf), sep="\t", index_col=0)
        _ldf.index = _ldf.index.astype(str)

        # Build per-cell series aligned to adata
        _score_map   = _ldf["oclr_score"].to_dict()    if "oclr_score"   in _ldf.columns else {}
        _label_map   = _ldf["oclr_endpoint_label"].to_dict() if "oclr_endpoint_label" in _ldf.columns else {}
        _suc_set     = set(_ldf.index[_ldf.get("is_oclr_success_endpoint",
                           _ldf.get("is_oclr_terminal", pd.Series(0))) == 1].tolist())
        _fail_set    = set(_ldf.index[_ldf.get("is_oclr_failure_endpoint",
                           pd.Series(0)) == 1].tolist()) if "is_oclr_failure_endpoint" in _ldf.columns else set()

    else:
        # Fallback: load separate barcode files
        _suc_files  = sorted(PROCESSED_DIR.glob("*_oclr_success_endpoint_barcodes.txt"))
        _fail_files = sorted(PROCESSED_DIR.glob("*_oclr_failure_endpoint_barcodes.txt"))
        if not _suc_files:
            raise FileNotFoundError(
                f"No OCLR endpoint label file found in {PROCESSED_DIR}.\n"
                "Run 02b2_define_oclr_endpoint.py first.\n"
                "Expected: *_oclr_endpoint_labels_hcipsc.tsv"
            )
        print(f"  [FALLBACK] Loading separate barcode files:")
        _suc_df  = pd.read_csv(str(_suc_files[-1]),  sep="\t", index_col=0)
        _suc_set = set(_suc_df.index.astype(str).tolist())
        _fail_set = set()
        if _fail_files:
            _fail_df  = pd.read_csv(str(_fail_files[-1]), sep="\t", index_col=0)
            _fail_set = set(_fail_df.index.astype(str).tolist())
        _score_map = {}
        if "oclr_score" in _suc_df.columns:
            _score_map.update(_suc_df["oclr_score"].to_dict())
        # Also try all-hcipsc score file
        _sf = sorted(PROCESSED_DIR.glob("*_oclr_score_all_hcipsc.csv"))
        if _sf:
            _sdf = pd.read_csv(str(_sf[-1]), index_col=0)
            if "oclr_score" in _sdf.columns:
                _score_map.update(_sdf["oclr_score"].to_dict())
        _label_map = {}
        print(f"    SUCCESS: {len(_suc_set):,} barcodes  FAILURE: {len(_fail_set):,} barcodes")

    # Report alignment
    _extra_suc  = _suc_set  - _proj_set
    _extra_fail = _fail_set - _proj_set
    if _extra_suc or _extra_fail:
        print(f"  [WARN] {len(_extra_suc)} success / {len(_extra_fail)} failure barcodes "
              f"not in adata.obs_names")
    print(f"  Matched SUCCESS: {len(_suc_set & _proj_set):,}  "
          f"FAILURE: {len(_fail_set & _proj_set):,}")

    # Build aligned series
    oclr_score_series = pd.Series(
        [_score_map.get(bc, np.nan) for bc in _proj_barcodes],
        index=_proj_barcodes, name="oclr_score", dtype=float,
    )
    success_mask = pd.Series(
        [bc in _suc_set for bc in _proj_barcodes],
        index=_proj_barcodes, name="oclr_success_endpoint_mask", dtype=bool,
    )
    failure_mask = pd.Series(
        [bc in _fail_set for bc in _proj_barcodes],
        index=_proj_barcodes, name="oclr_failure_endpoint_mask", dtype=bool,
    )
    endpoint_label_series = pd.Series(
        [_label_map.get(bc, "") for bc in _proj_barcodes],
        index=_proj_barcodes, name="oclr_endpoint_label", dtype=str,
    )

    return oclr_score_series, success_mask, failure_mask, endpoint_label_series


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

REQUIRED_OBS = ["day", "cell_growth_rate", "stage", "stage_day_label"]
_missing = [c for c in REQUIRED_OBS if c not in adata_full.obs.columns]
if _missing:
    raise ValueError(f"Required obs columns missing: {_missing}")
print(f"  Required obs columns : OK  {REQUIRED_OBS}")

# ── NaN/Inf quick check ───────────────────────────────────────────────────────
_samp = (adata_full.X[:500].toarray() if _is_sparse else np.asarray(adata_full.X[:500]))
if np.isnan(_samp).any() or np.isinf(_samp).any():
    raise ValueError("adata.X contains NaN/Inf in first 500 cells. Re-run preprocessing.")
print(f"  NaN/Inf check (500-cell sample) : PASSED")

# ── Growth prior: cell_growth_rate is WOT's multiplicative growth factor ───────
# Must be strictly positive.  log_growth_rate is used for reporting only.
# USE_GROWTH_RATES=False  →  neutral prior: 1.0 for every cell (no mass gain/loss)
# USE_GROWTH_RATES=True   →  use cell_growth_rate = exp(log_growth_rate) from 02b
if not USE_GROWTH_RATES:
    adata_full.obs["cell_growth_rate"] = 1.0
    print("  Growth prior : cell_growth_rate = 1.0 (neutral multiplicative; "
          "USE_GROWTH_RATES=False)")
else:
    _cgr = adata_full.obs["cell_growth_rate"].values.astype(float)
    _n_nan    = int((~np.isfinite(_cgr)).sum())
    _n_nonpos = int((_cgr <= 0).sum())
    if _n_nan > 0:
        raise ValueError(
            f"cell_growth_rate contains {_n_nan} non-finite values. "
            "Re-run 02b_estimate_growth_rates.py."
        )
    if _n_nonpos > 0:
        raise ValueError(
            f"cell_growth_rate contains {_n_nonpos} non-positive values. "
            "WOT requires strictly positive multiplicative growth factors. "
            "Re-run 02b_estimate_growth_rates.py."
        )
    _lgr_col = ("log_growth_rate"
                if "log_growth_rate" in adata_full.obs.columns else None)
    if _lgr_col:
        _lgr = adata_full.obs[_lgr_col].values.astype(float)
        print(f"  Growth prior : USE_GROWTH_RATES=True")
        print(f"    log_growth_rate  : [{_lgr.min():+.6f}, {_lgr.max():+.6f}]  (reporting only)")
    print(f"    cell_growth_rate : [{_cgr.min():.6f}, {_cgr.max():.6f}]  "
          f"(all finite: ✓  all >0: ✓)")

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

_gr_obs_cols = [c for c in ["growth_rate", "log_growth_rate", "cell_growth_rate"]
                if c in adata_wot.obs.columns]
if _gr_obs_cols:
    adata_wot.obs = adata_wot.obs.drop(columns=_gr_obs_cols)
    print(f"  Removed growth columns from ExprMatrix.h5ad obs: {_gr_obs_cols}")
    print(f"  — cell_growth_rate provided via external growth-rate file")

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
    "id":               adata.obs_names.astype(str),
    "cell_growth_rate": adata.obs["cell_growth_rate"].values.astype(float),
})
_gr.to_csv(str(GROWTH_RATES_PATH), sep="\t", index=False)
print(f"  growth_rates.txt → {GROWTH_RATES_PATH.name}  "
      f"(column: 'cell_growth_rate';  all zero: {(_gr['cell_growth_rate'] == 0).all()})")

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
    _cmd += ["--cell_growth_rates",       str(GROWTH_RATES_PATH),
             "--cell_growth_rates_field", "cell_growth_rate"]

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
# 8.  FATE COMPUTATION  (OCLR dual-endpoint design, WOT fates() API)
# =============================================================================
#
# Two competing fate sets are defined WITHIN the hCiPSC terminal sample using
# OCLR stemness labels from 02b2_define_oclr_endpoint.py:
#
#   fate_oclr_success  — top-10% OCLR within GSM7230012_hCiPSCs-0618
#   fate_oclr_failure  — bottom-10% OCLR within GSM7230012_hCiPSCs-0618
#
# Both sets are at the same timepoint (day30), satisfying the WOT fates()
# single-timepoint constraint.  If they span different days, a RuntimeError
# is raised immediately (not a warning).
#
# Output columns on adata_full.obs (primary):
#   p_oclr_success  — P(cell → high-stemness hCiPSC endpoint)
#   p_oclr_failure  — P(cell → low-stemness hCiPSC endpoint)
#   p_oclr_margin   — p_oclr_success − p_oclr_failure  [PRIMARY benchmark metric]
# Legacy alias:
#   p_hcipsc        = p_oclr_success  (kept for backward compatibility)
#
print(f"\n{_SEP}")
print("STEP 8  —  Fate computation  (OCLR dual-endpoint, WOT fates() API)")
print(_SEP)

# Tracking variables initialised here; updated inside the if-block below.
_n_hcipsc       = 0
_mask_cipsc     = None
_p_hcipsc_stats = {}

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
    # 8b.  Load OCLR dual endpoints and define fate sets
    # ─────────────────────────────────────────────────────────────────────────
    # fate_oclr_success  — top-10% OCLR within hCiPSC  (loaded from 02b2)
    # fate_oclr_failure  — bottom-10% OCLR within hCiPSC  (loaded from 02b2)
    #
    # Both sets must be at the same timepoint (day30 hCiPSC sample).
    # This satisfies the WOT fates() single-timepoint constraint.
    if RUN_FATE_COMPUTATION:
        print(f"\n{_SEP2}")
        print("  Step 8b — Load OCLR dual endpoints (from 02b2 output)")

        _mask_cipsc_bool = (
            _mask_cipsc.values if hasattr(_mask_cipsc, "values") else _mask_cipsc
        )
        _n_hcipsc = int(_mask_cipsc_bool.sum())
        print(f"  hCiPSC cells in adata (stage_std == 'hCiPSC') : {_n_hcipsc:,}")

        try:
            (
                _oclr_score_series,
                _oclr_success_mask,
                _oclr_failure_mask,
                _oclr_label_series,
            ) = load_oclr_endpoints(adata_full)

            adata_full.obs["oclr_score"]                = _oclr_score_series.values
            adata_full.obs["oclr_success_endpoint_mask"]= _oclr_success_mask.astype(int).values
            adata_full.obs["oclr_failure_endpoint_mask"]= _oclr_failure_mask.astype(int).values
            adata_full.obs["oclr_endpoint_label"]       = _oclr_label_series.values
            # Legacy alias
            adata_full.obs["oclr_ips_subset_mask"]      = _oclr_success_mask.astype(int).values

            _n_suc = int(_oclr_success_mask.sum())
            _n_fail= int(_oclr_failure_mask.sum())
            print(f"  fate_oclr_success : {_n_suc:,} cells")
            print(f"  fate_oclr_failure : {_n_fail:,} cells")

            if _n_suc < MIN_FATE_CELLS:
                print(f"  !! FATAL: Only {_n_suc} success endpoint cells "
                      f"(need ≥ {MIN_FATE_CELLS}).")
                RUN_FATE_COMPUTATION = False
            if _n_fail < MIN_FATE_CELLS:
                print(f"  !! FATAL: Only {_n_fail} failure endpoint cells "
                      f"(need ≥ {MIN_FATE_CELLS}).")
                RUN_FATE_COMPUTATION = False

            # Verify both endpoints are at the same day
            _suc_days  = sorted(adata_full.obs.loc[
                _oclr_success_mask.values.astype(bool), "day"
            ].unique())
            _fail_days = sorted(adata_full.obs.loc[
                _oclr_failure_mask.values.astype(bool), "day"
            ].unique())
            if len(set(_suc_days) | set(_fail_days)) > 1:
                raise RuntimeError(
                    "OCLR success and failure endpoints span different timepoints.\n"
                    f"  Success endpoint days : {_suc_days}\n"
                    f"  Failure endpoint days : {_fail_days}\n"
                    "Both endpoints must be within the same hCiPSC terminal sample "
                    "(same day).\nThe WOT fates() API requires all population objects "
                    "at a single timepoint.\n"
                    "Verify that 02b2_define_oclr_endpoint.py is restricting scoring "
                    "to a single sample (ENDPOINT_SAMPLE_ID = "
                    f"'{OCLR_ENDPOINT_SAMPLE_ID}')."
                )
            _endpoint_day = float(_suc_days[0]) if _suc_days else np.nan
            print(f"  Endpoint timepoint (day) : {_endpoint_day}")

        except FileNotFoundError as _e:
            print(f"  !! FATAL: OCLR endpoint file not found:\n    {_e}")
            print(f"  Run 02b2_define_oclr_endpoint.py first.")
            RUN_FATE_COMPUTATION = False
            _oclr_success_mask = pd.Series(
                np.zeros(adata_full.n_obs, dtype=bool),
                index=adata_full.obs_names,
            )
            _oclr_failure_mask = _oclr_success_mask.copy()
            _n_suc = _n_fail = 0
            _endpoint_day = np.nan

        adata_full.obs["fate_oclr_success"] = _oclr_success_mask.astype(int).values
        adata_full.obs["fate_oclr_failure"] = _oclr_failure_mask.astype(int).values
        # Legacy alias
        adata_full.obs["fate_hcipsc"]       = _oclr_success_mask.astype(int).values

    # ─────────────────────────────────────────────────────────────────────────
    # 8c.  Build cell_sets dict for WOT fates() API
    # ─────────────────────────────────────────────────────────────────────────
    # wot.tmap.TransportMapModel.fates(populations) takes a list of Population
    # objects all at the same timepoint.
    # population_from_cell_sets(cell_sets_dict, at_time) converts a dict of
    # {name: [cell_ids]} to the required Population list.
    if RUN_FATE_COMPUTATION:
        print(f"\n{_SEP2}")
        print("  Step 8c — Build cell_sets dict for WOT fates() API")

        _cell_bcs = adata_full.obs_names.astype(str).tolist()

        _success_barcodes = [
            bc for bc in _cell_bcs
            if bc in _oclr_success_mask.index and bool(_oclr_success_mask.loc[bc])
        ]
        _failure_barcodes = [
            bc for bc in _cell_bcs
            if bc in _oclr_failure_mask.index and bool(_oclr_failure_mask.loc[bc])
        ]

        _cell_sets_dict = {
            "fate_oclr_success": _success_barcodes,
            "fate_oclr_failure": _failure_barcodes,
        }
        print(f"  fate_oclr_success : {len(_success_barcodes):,} cells  "
              f"(at day {_endpoint_day})")
        print(f"  fate_oclr_failure : {len(_failure_barcodes):,} cells  "
              f"(at day {_endpoint_day})")

    # ─────────────────────────────────────────────────────────────────────────
    # 8d.  Fate computation — WOT Python API  (fates() method)
    # ─────────────────────────────────────────────────────────────────────────
    #
    # wot.tmap.TransportMapModel.from_directory() loads the transport maps.
    # .population_from_cell_sets(cell_sets_dict, at_time) builds Population list.
    # .fates(populations) propagates backward and returns an AnnData with
    #   obs = all project cells, var = fate names, X[i,k] = fate probability.
    #
    # Both endpoints are at the same timepoint (day30) — fates() constraint met.
    # ─────────────────────────────────────────────────────────────────────────
    if RUN_FATE_COMPUTATION:
        print(f"\n{_SEP2}")
        print("  Step 8d — Fate computation  (WOT fates() API)")

        # Load transport map model
        # WOT 1.0.8.post2: from_directory() takes a single positional argument
        # which must be the full tmap prefix (same string passed to --out in
        # the CLI step).  The prefix= keyword was removed and must not be used.
        try:
            _tmap_model = wot.tmap.TransportMapModel.from_directory(
                str(TMAPS_DIR / "tmaps")
            )
            print(f"  TransportMapModel loaded from prefix : {TMAPS_DIR.name}/tmaps")
        except AttributeError:
            import wot.tmap as _wot_tmap
            _tmap_model = _wot_tmap.TransportMapModel.from_directory(
                str(TMAPS_DIR / "tmaps")
            )
            print(f"  TransportMapModel (fallback import) loaded from prefix : {TMAPS_DIR.name}/tmaps")

        # Build Population objects from cell_sets dict
        print(f"  Building Population objects (at_time={_endpoint_day}) …")
        _populations = _tmap_model.population_from_cell_sets(
            _cell_sets_dict, at_time=_endpoint_day
        )
        print(f"  Populations built : {[p.name for p in _populations]}")

        # Compute fate probabilities using fates()
        print(f"  Computing fates() …")
        _fate_ds = _tmap_model.fates(_populations)
        # _fate_ds: AnnData — obs = all cells, var = fate names, X = fate probs

        print(f"  Fate probability AnnData : {_fate_ds.n_obs:,} cells × {_fate_ds.n_vars} fates")
        print(f"  Fate names               : {list(_fate_ds.var_names)}")

        # ── Align fate probabilities back to adata_full ───────────────────────
        _fate_df = pd.DataFrame(
            _fate_ds.X if not sp.issparse(_fate_ds.X) else _fate_ds.X.toarray(),
            index=_fate_ds.obs_names.astype(str),
            columns=_fate_ds.var_names.astype(str),
        )

        _p_success_vec = np.full(adata_full.n_obs, np.nan, dtype=np.float64)
        _p_failure_vec = np.full(adata_full.n_obs, np.nan, dtype=np.float64)

        for _i, _bc in enumerate(adata_full.obs_names.astype(str)):
            if _bc in _fate_df.index:
                if "fate_oclr_success" in _fate_df.columns:
                    _p_success_vec[_i] = float(_fate_df.loc[_bc, "fate_oclr_success"])
                if "fate_oclr_failure" in _fate_df.columns:
                    _p_failure_vec[_i] = float(_fate_df.loc[_bc, "fate_oclr_failure"])

        # Primary columns
        adata_full.obs["p_oclr_success"] = _p_success_vec
        adata_full.obs["p_oclr_failure"] = _p_failure_vec
        adata_full.obs["p_oclr_margin"]  = _p_success_vec - _p_failure_vec
        # Legacy alias
        adata_full.obs["p_hcipsc"]       = _p_success_vec

        # ── Sanity checks ─────────────────────────────────────────────────────
        _n_fin_suc  = int(np.isfinite(_p_success_vec).sum())
        _n_fin_fail = int(np.isfinite(_p_failure_vec).sum())
        _non_terminal_mask = (
            ~_mask_cipsc_bool
            & np.isfinite(_p_success_vec)
            & np.isfinite(_p_failure_vec)
        )
        _p_suc_nt  = _p_success_vec[_non_terminal_mask]
        _p_fail_nt = _p_failure_vec[_non_terminal_mask]
        _p_marg_nt = _p_suc_nt - _p_fail_nt

        print(f"\n  p_oclr_success finite : {_n_fin_suc:,}  "
              f"({_n_fin_suc / adata_full.n_obs:.1%})")
        print(f"  p_oclr_failure finite : {_n_fin_fail:,}  "
              f"({_n_fin_fail / adata_full.n_obs:.1%})")

        def _fate_stats_line(arr: np.ndarray, label: str) -> str:
            q25, med, q75 = np.nanpercentile(arr, [25, 50, 75])
            rng = float(np.nanmax(arr) - np.nanmin(arr))
            iqr = float(q75 - q25)
            return (
                f"  {label} (non-terminal, n={len(arr):,}):\n"
                f"    min={np.nanmin(arr):.4f}  q25={q25:.4f}  "
                f"median={med:.4f}  q75={q75:.4f}  max={np.nanmax(arr):.4f}\n"
                f"    std={np.nanstd(arr):.4f}  range={rng:.4f}  iqr={iqr:.4f}"
            )

        if len(_p_suc_nt) > 0:
            print(_fate_stats_line(_p_suc_nt,  "p_oclr_success"))
            print(_fate_stats_line(_p_fail_nt, "p_oclr_failure"))
            print(_fate_stats_line(_p_marg_nt, "p_oclr_margin "))

            _suc_range  = float(np.nanmax(_p_suc_nt)  - np.nanmin(_p_suc_nt))
            _marg_range = float(np.nanmax(_p_marg_nt) - np.nanmin(_p_marg_nt))
            _marg_std   = float(np.nanstd(_p_marg_nt))

            if _suc_range < 0.01 and _marg_range < 0.01:
                print("  !! WARNING: near-degenerate fate signal — "
                      "success-range and margin-range both < 0.01.")
                print("     Fate probabilities are near-constant.  "
                      "Consider enabling USE_GROWTH_RATES=True.")
            elif _marg_std < 0.005:
                print("  !! WARNING: near-degenerate fate signal — "
                      f"p_oclr_margin std={_marg_std:.5f} (< 0.005).")
                print("     Margin has almost no variance across non-terminal cells.")
            else:
                print(f"  ✓  Fate signal non-degenerate  "
                      f"(suc_range={_suc_range:.4f}  "
                      f"marg_range={_marg_range:.4f}  "
                      f"marg_std={_marg_std:.4f})")

            # Keep slim dict for backward-compatible final-summary reference
            _p_hcipsc_stats = {
                "min":    float(np.nanmin(_p_suc_nt)),
                "median": float(np.nanmedian(_p_suc_nt)),
                "max":    float(np.nanmax(_p_suc_nt)),
                "range":  _suc_range,
                "marg_std": _marg_std,
            }

    # ── Save adata_full with all new columns ──────────────────────────────────
    if RUN_FATE_COMPUTATION:
        _out_h5ad = METRICS_DIR / f"{TIMESTAMP}_GSE230659_wot_fates.h5ad"
        adata_full.write_h5ad(str(_out_h5ad))
        print(f"\n  Saved adata_full with fate columns → {_out_h5ad.name}")
        print(f"    obs columns added (primary): "
              f"p_oclr_success, p_oclr_failure, p_oclr_margin, "
              f"fate_oclr_success, fate_oclr_failure, "
              f"oclr_score, oclr_success_endpoint_mask, oclr_failure_endpoint_mask, "
              f"oclr_endpoint_label")
        print(f"    obs columns added (legacy): p_hcipsc, fate_hcipsc, oclr_ips_subset_mask")

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

    # ── Plot A: continuous day vs p_oclr_success scatter ─────────────────────
    try:
        _fig, _ax = plt.subplots(figsize=(8, 5))
        _plot_col = "p_oclr_success" if "p_oclr_success" in adata_full.obs.columns else "p_hcipsc"
        _sc_data  = adata_full.obs[["day", _plot_col, "stage_std"]].dropna(
            subset=[_plot_col]
        )
        for _st, _col in STAGE_PALETTE.items():
            _mask_st = _sc_data["stage_std"].apply(
                lambda x: normalize_stage(x) == _st
            )
            _sub = _sc_data[_mask_st]
            if len(_sub):
                _ax.scatter(_sub["day"], _sub[_plot_col],
                            c=_col, s=2, alpha=0.3, label=_st, rasterized=True)
        _ax.set_xlabel("Day")
        _ax.set_ylabel(_plot_col)
        _ax.set_title("WOT p_oclr_success — OCLR dual-endpoint fate probability")
        _ax.legend(markerscale=4, fontsize=8)
        _fig.tight_layout()
        _fig_path_A = FIGURES_DIR / f"{TIMESTAMP}_fate_prob_day.png"
        _fig.savefig(str(_fig_path_A), dpi=150, bbox_inches="tight")
        plt.close(_fig)
        print(f"  Plot A (day vs {_plot_col}) → {_fig_path_A.name}")
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

    # ── Plot C: p_oclr_success on UMAP ───────────────────────────────────────
    if "X_umap" in adata_full.obsm:
        try:
            _fig_C, _axes_C = plt.subplots(1, 2, figsize=(14, 6))
            for _ci, (_pcol, _clabel) in enumerate([
                ("p_oclr_success", "p_oclr_success"),
                ("p_oclr_margin",  "p_oclr_margin (success − failure)"),
            ]):
                if _pcol not in adata_full.obs.columns:
                    continue
                _p_all  = adata_full.obs[_pcol].values.astype(float)
                _finite = np.isfinite(_p_all)
                _cmap   = "RdYlBu_r" if _pcol == "p_oclr_success" else "RdBu_r"
                _vcent  = None if _pcol == "p_oclr_success" else 0.0
                _vabs   = np.nanmax(np.abs(_p_all[_finite])) if _vcent == 0.0 else 1.0
                _sc = _axes_C[_ci].scatter(
                    adata_full.obsm["X_umap"][_finite, 0],
                    adata_full.obsm["X_umap"][_finite, 1],
                    c=_p_all[_finite], cmap=_cmap,
                    s=3, alpha=0.5,
                    vmin=-_vabs if _vcent == 0.0 else 0,
                    vmax= _vabs if _vcent == 0.0 else 1,
                    rasterized=True, linewidths=0,
                )
                plt.colorbar(_sc, ax=_axes_C[_ci], label=_clabel)
                _axes_C[_ci].set_xlabel("UMAP 1"); _axes_C[_ci].set_ylabel("UMAP 2")
                _axes_C[_ci].set_title(f"WOT {_pcol}")
                _axes_C[_ci].set_xticks([]); _axes_C[_ci].set_yticks([])
            _fig_C.tight_layout()
            _fig_path_C = FIGURES_DIR / f"{TIMESTAMP}_umap_p_oclr.png"
            _fig_C.savefig(str(_fig_path_C), dpi=150, bbox_inches="tight")
            plt.close(_fig_C)
            print(f"  Plot C (UMAP p_oclr) → {_fig_path_C.name}")
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
print(f"  ── Fate computation  (OCLR dual-endpoint, WOT fates() API) ─────────────")
print(f"  fate_oclr_success  : "
      + (f"{_n_suc:,} cells  (top-10% OCLR within hCiPSC)" if RUN_FATE_COMPUTATION else "skipped"))
print(f"  fate_oclr_failure  : "
      + (f"{_n_fail:,} cells  (bottom-10% OCLR within hCiPSC)" if RUN_FATE_COMPUTATION else "skipped"))
if RUN_FATE_COMPUTATION and _p_hcipsc_stats:
    print(f"  p_oclr_success (non-terminal): "
          f"min={_p_hcipsc_stats['min']:.4f}  "
          f"median={_p_hcipsc_stats['median']:.4f}  "
          f"max={_p_hcipsc_stats['max']:.4f}  "
          f"range={_p_hcipsc_stats['range']:.4f}  "
          f"marg_std={_p_hcipsc_stats['marg_std']:.4f}")
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
    print(f"     Check warnings above (OCLR endpoint file missing or too few cells).")
else:
    print(f"\n  ✓  WOT CLI + OCLR dual-endpoint fate computation complete.")
    print(f"     Key outputs (adata_full.obs):")
    print(f"       p_oclr_success       : P(cell → high-stemness hCiPSC fate)")
    print(f"       p_oclr_failure       : P(cell → low-stemness hCiPSC fate)")
    print(f"       p_oclr_margin        : p_oclr_success − p_oclr_failure  [PRIMARY metric]")
    print(f"       fate_oclr_success    : binary endpoint membership")
    print(f"       fate_oclr_failure    : binary endpoint membership")
    print(f"       oclr_score           : continuous OCLR stemness score")
    print(f"       oclr_success_endpoint_mask / oclr_failure_endpoint_mask: binary flags")
    print(f"       oclr_endpoint_label  : 'success'|'failure'|'ambiguous'|'unscored'")
    print(f"     Legacy aliases: p_hcipsc = p_oclr_success, fate_hcipsc = fate_oclr_success")
    _h5ad_out = METRICS_DIR / f"{TIMESTAMP}_GSE230659_wot_fates.h5ad"
    print(f"     h5ad saved to : {_h5ad_out.name}")