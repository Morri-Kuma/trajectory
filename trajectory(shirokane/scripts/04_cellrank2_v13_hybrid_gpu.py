#!/usr/bin/env python3
"""
04_cellrank2_v12_gpt.py  [v12 — dynamic-vs-biological iPSC endpoint comparison]
==========================================================================================
Project : Comparative Study of Trajectory Inference Models for Chemical iPSC
          Reprogramming
Dataset : GSE230659 (human, Liuyang et al. 2023 Cell Stem Cell)
Step    : 4 — CellRank2 trajectory analysis (second method, comparative study)
Env     : conda activate traj_env
Analysis design
---------------
This script implements CellRank2 RealTimeKernel as a WOT-based trajectory model.
  - WOT               = baseline fate-probability model (Problem A, script 02c)
  - CellRank2 RTK     = WOT-informed CellRank2 model that adds within-timepoint
                        connectivity on top of the same WOT transport maps
  - The goal is to compare classical WOT fate mapping with CellRank2 fate mapping
    while keeping the cross-time transport information fixed.
IMPORTANT CHANGE FROM v9
------------------------
This version further tightens terminal-state handling to avoid late-stage
collapse when multiple macrostates are all majority-labeled as hCiPSC.
  1. Compute macrostates with GPCCA.
  2. Use GPCCA.predict_terminal_states(method="stability") to automatically
     identify terminal macrostates from coarse_T.
  3. Summarize each macrostate using biological evidence available in this
     dataset (hCiPSC stage composition, fate_ips, p_iPSC, late-time enrichment).
  4. If the automatic call misses a strongly iPSC-like late macrostate, add it
     manually via GPCCA.set_terminal_states(...), mirroring the official tutorial
     pattern of “automatic prediction first, then biologically informed manual
     confirmation/supplement”.
  5. Compute fate probabilities on the final terminal-state set, then define
     cr2_fate_ips as the summed fate probability toward the iPSC terminal
     lineage/lineages only.
This means:
  - We NO LONGER set all macrostates as terminal states.
  - We NO LONGER select iPSC lineages from the full macrostate pool.
  - Instead, we first define biologically plausible terminal macrostates, then
    read out iPSC-specific fate probabilities from that terminal set.
Changelog v10
-------------------------------------------------------------
Rewrite — terminal-state workflow aligned with official CellRank tutorial
  Fix A — use predict_terminal_states(method="stability") as the default first step
  Fix B — summarize macrostates and manually supplement ONLY if auto-prediction
           misses a strongly iPSC-like late macrostate
  Fix C — compute cr2_fate_ips only from confirmed iPSC terminal lineages
  Fix D — save terminal-state summary table instead of Scheme-B enrichment table
  Fix E — driver-gene computation preferentially restricts to fate-relevant stages
Changelog v12
-------------------------------------------------------------
Dual endpoint definition — compare CellRank dynamics vs biological markers
  Fix K — define dynamic iPSC terminal cells from CellRank-selected iPSC terminal lineages
  Fix L — define biological iPSC terminal cells using mutually exclusive marker classes:
           POU5F1-only, single_high, bi_high, tri_high, tri_very_high
  Fix M — compare dynamic vs biological iPSC terminal cells within the terminal-cell universe
           and report overlap counts, coverage rates and Jaccard index
  Fix N — final iPSC detection is BIOLOGY-FIRST (method 2), stored in ips_final_mask
  Fix O — manual driver-gene fallback uses final biological iPSC labels when available

Changelog v11
-------------------------------------------------------------
Tightening — finer macrostate annotation + stricter iPSC-terminal selection
  Fix F — compute Leiden clusters and use cluster_key='leiden' for macrostate naming
  Fix G — remove name-based iPSC terminal selection ('contains ips')
  Fix H — manual supplement is allowed ONLY when auto terminal states do not
           already contain a convincing iPSC endpoint
  Fix I — require late-time enrichment + multi-metric evidence to define iPSC terminals
  Fix J — cap the number of selected iPSC terminal lineages to avoid cr2_fate_ips ≈ 1
Changelog v9
-------------------------------------------------------------
Fix A — terminal-states obs key ('term_states_fwd', not 'terminal_states_fwd')
Fix B — eigenspectrum: real_only=True, no ax=, capture via plt.gcf()
Fix C — lineage selection: macrostate-fraction-based 4-tier cascade
Fix D — driver-gene trend: .values on boolean Series before sparse indexing
Fix E — suppress 'more than two extensions' CellRank INFO messages
Changelog v6
-------------------------------------------------------------
Fix 8 — Windows multiprocessing guard (main() + freeze_support)
Changelog v5
-------------------------------------------------------------
Fix 7 — local_debug tmap cropping
Changelog v4
-------------------------------------------------------------
New: dual execution mode (RUN_MODE = "local_debug" | "hpc_full")
Changelog v2 (source-verified against scverse/cellrank main)
-------------------------------------------------------------
Fix 1  — adata.obs["day"] must be pd.Categorical.
Fix 2  — compute_transition_matrix() correct parameters.
Fix 3  — ConnectivityKernel mix removed (redundant).
Fix 4  — compute_lineage_drivers() column-name convention.
Fix 5  — fate_probabilities Lineage access pattern.
Usage
-----
  # Local debug (default, Windows-safe):
  conda activate traj_env
  cd C:\\Users\\37620\\trajectory
  python scripts/WOT/04_cellrank2.py
  # HPC full run (two ways):
  RUN_MODE=hpc_full python scripts/WOT/04_cellrank2.py
  # — or change RUN_MODE in the CONFIG section below —
Reference
---------
  Weiler P, et al. (2024) CellRank 2: unified fate mapping in multiview
  single-cell data. Nature Methods 21, 1196–1205.
  https://doi.org/10.1038/s41592-024-02303-9
"""
# =============================================================================
# MODULE-LEVEL: imports, config constants, helper functions
#
# IMPORTANT — Windows multiprocessing safety
# ------------------------------------------
# On Windows, Python's multiprocessing uses the "spawn" start method.
# Spawned child processes import this file as a module.  Code at module level
# (outside any function or if __name__ == "__main__":) therefore runs again
# inside each child process.
#
# Safe at module level: imports, simple constant assignments, function defs.
# NOT safe at module level: I/O, data loading, CellRank API calls, prints.
#
# All pipeline execution logic lives in main() below.
# =============================================================================
import os
import re
import shutil
import sys
import time as _time
import warnings
from datetime import datetime
from multiprocessing import freeze_support
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ── Optional GPU acceleration (hybrid mode) ──────────────────────────────────
# This script can accelerate Scanpy-like preprocessing / embedding steps on
# NVIDIA GPUs via rapids-singlecell, but CellRank/GPCCA remains on CPU unless
# deep library internals are changed.
USE_GPU: bool = os.environ.get("USE_GPU", "1").lower() not in {"0", "false", "no"}
GPU_FORCE_RECOMPUTE_PCA: bool = os.environ.get("GPU_FORCE_RECOMPUTE_PCA", "1").lower() not in {"0", "false", "no"}
GPU_PCA_COMPONENTS: int = int(os.environ.get("GPU_PCA_COMPONENTS", "50"))
# =============================================================================
# 0.  CONFIGURATION  (module level — safe to read on child-process import)
# =============================================================================
# NOTE: TERMINAL_STAGE constant has been removed (v3).  The iPSC stage label is
# detected dynamically from adata.obs["stage"] using case-insensitive "ips"
# matching.  _ips_stage_label is set inside main() and used throughout.
# ── Execution mode ────────────────────────────────────────────────────────────
# "local_debug" : stratified subsample + tmap cropping + n_jobs=1 fate probs
#                 → Windows-safe pipeline smoke test on a 32 GB workstation
# "hpc_full"    : full dataset, original tmaps, n_jobs=None (all cores)
#                 → final biological-conclusion run on Shirokane
#
# Override via environment variable:
#   RUN_MODE=hpc_full python scripts/WOT/04_cellrank2.py
RUN_MODE: str = os.environ.get("RUN_MODE", "local_debug")   # ← edit here OR env var
assert RUN_MODE in ("local_debug", "hpc_full"), (
    f"RUN_MODE must be 'local_debug' or 'hpc_full', got '{RUN_MODE}'"
)
# ── Local-debug subsampling parameters ───────────────────────────────────────
# Target total cells after stratified (stage × day) sampling.
# 8 000 cells → transition matrix ≈ 0.5 GB sparse; Schur decomp ≈ 3–4 GB peak.
LOCAL_TARGET_CELLS:    int = 8_000
# Minimum cells kept per (stage × day) stratum.
# Prevents rare terminal populations (e.g. hCiPSC) from being dropped.
LOCAL_MIN_PER_STRATUM: int = 20
# Random seed for reproducible subsampling.
LOCAL_SEED:            int = 42
# ── Project root ──────────────────────────────────────────────────────────────
_candidates = [
    Path("/home/xzy0723/projects/trajectory"),
    Path(r"C:\Users\37620\trajectory"),
    Path("/sessions/cool-admiring-hawking/mnt/trajectory"),
    Path(__file__).resolve().parents[1],
]
PROJECT_ROOT = next((p for p in _candidates if p.exists()), Path("/home/xzy0723/projects/trajectory"))
PROCESSED_DIR = PROJECT_ROOT / "data"    / "processed"
FIGURES_DIR   = PROJECT_ROOT / "results" / "figures"
METRICS_DIR   = PROJECT_ROOT / "results" / "metrics"
TMAPS_BASE    = PROJECT_ROOT / "results" / "tmaps"
CR2_OUT_DIR   = PROJECT_ROOT / "results" / "cellrank2"
# ── RealTimeKernel parameters ─────────────────────────────────────────────────
# self_transitions="all" → every timepoint block gets:
#   diag(t)          = conn_weight × connectivity
#   off-diag(t→t+1)  = (1 - conn_weight) × WOT transport
CONN_WEIGHT: float = 0.2    # connectivity weight at each timepoint
N_NEIGHBORS: int   = 30     # kNN for sc.pp.neighbors (used in "all" mode)
THRESHOLD:   str   = "auto" # "auto" | "auto_local" | float ∈ [0, 100]
# ── GPCCA parameters ──────────────────────────────────────────────────────────
N_SCHUR_COMPONENTS: int = 20
N_MACROSTATES:      int = 6
# ── Fate probability parallelism ──────────────────────────────────────────────
# compute_fate_probabilities() uses joblib.Parallel internally to solve the
# linear systems for each terminal state.
#
# n_jobs=1  → joblib runs serially; NO child processes are spawned.
#             Required in local_debug on Windows to avoid recursive re-execution.
#             Verified against cellrank source:
#               _parallelize.py _get_n_cores() returns 1 when n_jobs==1
#               joblib.Parallel(n_jobs=1) always executes inline (no spawn).
#
# n_jobs=None → joblib default: uses all available cores via loky backend.
#               Appropriate for hpc_full on a many-core server.
FATE_PROB_N_JOBS: int | None = 1 if RUN_MODE == "local_debug" else None
# ── Leiden clustering + stricter terminal-state workflow ─────────────────────
# CellRank's `cluster_key` should point to a fine-grained annotation, not a very
# coarse stage label.  We therefore compute Leiden clusters (if absent) and use
# them for GPCCA macrostate naming/composition.
LEIDEN_KEY: str           = "leiden"
LEIDEN_RESOLUTION: float  = 0.6
LEIDEN_RANDOM_STATE: int  = 0

# Step 1: auto-predict terminal macrostates from coarse_T using stability.
# Step 2: inspect macrostate composition using late-time/iPSC evidence.
# Step 3: ONLY if the auto-predicted terminal set does not already contain a
#         convincing iPSC endpoint, manually add at most a small number of
#         missed late macrostates with strong evidence.
# Step 4: among the final terminal-state set, keep only the strongest
#         evidence-supported iPSC terminal lineage(s).
TERMINAL_METHOD: str                = "stability"
TERMINAL_STABILITY_THRESHOLD: float = 0.96
TERMINAL_N_CELLS: int               = 30
MIN_MACROSTATE_CELLS: int           = 5

# Manual supplement is intentionally strict to avoid collapsing all late
# macrostates into iPSC terminal states.
MANUAL_ADD_IPS_STAGE_FRAC: float    = 0.80
MANUAL_ADD_FATE_IPS_FRAC: float     = 0.35
MANUAL_ADD_PIPSC_THRESHOLD: float   = 0.45
MANUAL_ADD_LAST_DAY_FRAC: float     = 0.80
MANUAL_ADD_MIN_EVIDENCE_COUNT: int  = 2
MANUAL_ADD_MAX_STATES: int          = 1

# iPSC terminal lineage identification inside the final terminal-state set.
IPS_TERMINAL_STAGE_FRAC: float      = 0.80
IPS_TERMINAL_FATE_IPS_FRAC: float   = 0.35
IPS_TERMINAL_PIPSC_THRESHOLD: float = 0.45
IPS_TERMINAL_LAST_DAY_FRAC: float   = 0.80
IPS_TERMINAL_MIN_EVIDENCE_COUNT: int = 2
IPS_TERMINAL_MAX_STATES: int        = 2

# ── Biological iPSC marker-based endpoint definition (method 2) ─────────────
# Final iPSC detection will use these mutually exclusive classes among terminal
# cells: POU5F1-only, single_high, bi_high, tri_high, tri_very_high.
# Thresholds are estimated from the log-normalized expression distribution of
# all cells with non-zero expression for each marker gene.
BIO_MARKER_ALIASES: dict[str, tuple[str, ...]] = {
    "POU5F1": ("POU5F1", "OCT4", "OCT3/4"),
    "SOX2": ("SOX2",),
    "NANOG": ("NANOG",),
}
BIO_HIGH_QUANTILE: float = 0.80
BIO_VERY_HIGH_QUANTILE: float = 0.95
BIO_MIN_HIGH_EXPR: float = 1.00
BIO_MIN_VERY_HIGH_EXPR: float = 1.50
BIO_MIN_VERY_HIGH_DELTA: float = 0.25

# ── Lineage driver genes ──────────────────────────────────────────────────────
N_DRIVER_GENES: int = 300
DRIVER_TOP_QUANTILE: float = 0.75
DRIVER_MIN_STAGE_SHARE: float = 0.05
# ── Plotting ──────────────────────────────────────────────────────────────────
UMAP_SIZE:  float = 2.0
UMAP_ALPHA: float = 0.5
SAVE_FMT:   str   = "png"
# ── Output tag (injected into every output file name) ─────────────────────────
_MODE_TAG = "localdebug" if RUN_MODE == "local_debug" else "hpcfull"
# =============================================================================
# HELPER FUNCTIONS  (module level — safe to import in child processes)
# =============================================================================
def find_latest(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern), key=lambda p: p.stat().st_mtime)
    if not matches:
        raise FileNotFoundError(f"No files matching '{pattern}' in {directory}.")
    return matches[-1]
def extract_ts(filename: str) -> str | None:
    m = re.match(r"(\d{8}_\d{4})", filename)
    return m.group(1) if m else None
def find_tmaps_dir(base: Path, ts: str | None) -> Path:
    """Return the tmaps directory matching *ts*, or the most-recently modified one."""
    if ts is not None:
        cand = base / (ts + "_tmaps")
        if cand.exists() and list(cand.glob("tmaps_*.h5ad")):
            return cand
    dirs = sorted(
        [d for d in base.iterdir()
         if d.is_dir() and d.name.endswith("_tmaps")
            and list(d.glob("tmaps_*.h5ad"))],
        key=lambda d: d.stat().st_mtime,
    )
    if not dirs:
        raise FileNotFoundError(f"No transport-map directories found under {base}.")
    return dirs[-1]
def savefig(fig: plt.Figure, name: str, figures_dir: Path, dpi: int = 150) -> None:
    p = figures_dir / (name + "." + SAVE_FMT)
    fig.savefig(str(p), dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  [FIG] {p.name}")

def _bool_to_int(x: bool) -> int:
    return int(bool(x))


def _safe_mean(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    return float(np.nanmean(values)) if values.size else np.nan


def _safe_frac(mask: np.ndarray) -> float:
    mask = np.asarray(mask, dtype=bool)
    return float(mask.mean()) if mask.size else np.nan


def stratified_subsample(
    adata_full:      sc.AnnData,
    target_n:        int,
    min_per_stratum: int,
    strat_cols:      list[str],
    seed:            int,
) -> sc.AnnData:
    """
    Return a biologically stratified subsample of *adata_full*.
    Strategy
    --------
    1. Form strata from unique combinations of *strat_cols*.
    2. Allocate cells to each stratum proportionally to its size.
    3. Clip each allocation to [min_per_stratum, stratum_size].
    4. Sample without replacement using *seed*.
    5. Return adata[selected_barcodes].copy().
    """
    rng    = np.random.default_rng(seed)
    obs    = adata_full.obs.copy()
    n_full = len(obs)
    _key = obs[strat_cols].astype(str).agg("__".join, axis=1)
    obs["_stratum"] = _key
    strata   = _key.unique()
    n_strata = len(strata)
    selected_idx: list[str] = []
    for _s in strata:
        _mask  = obs["_stratum"] == _s
        _idx   = obs.index[_mask].tolist()
        _n     = len(_idx)
        _alloc = int(round(target_n * _n / n_full))
        _alloc = max(min_per_stratum, min(_alloc, _n))
        selected_idx.extend(rng.choice(_idx, size=_alloc, replace=False).tolist())
    selected_idx = list(dict.fromkeys(selected_idx))   # deduplicate, preserve order
    adata_sub    = adata_full[selected_idx].copy()
    print(f"\n  [SUBSAMPLE] stratified by {strat_cols}")
    print(f"    Full dataset  : {n_full:,} cells")
    print(f"    Target        : {target_n:,}  (min_per_stratum={min_per_stratum})")
    print(f"    Actual sample : {adata_sub.n_obs:,} cells  "
          f"({100 * adata_sub.n_obs / n_full:.1f}% of full)")
    print(f"    N strata      : {n_strata}")
    print("\n    Per-stage counts (full → sampled):")
    for _stg in sorted(adata_full.obs["stage"].unique()):
        _nf = int((adata_full.obs["stage"] == _stg).sum())
        _ns = int((adata_sub.obs["stage"]  == _stg).sum())
        print(f"      {_stg:<22s}  {_nf:>6,}  →  {_ns:>5,}")
    print("\n    Per-day counts (full → sampled):")
    _days_sorted = sorted(
        adata_full.obs["day"].cat.categories.tolist()
        if hasattr(adata_full.obs["day"], "cat")
        else adata_full.obs["day"].unique().tolist()
    )
    for _d in _days_sorted:
        _nf = int((adata_full.obs["day"].astype(str) == str(_d)).sum())
        _ns = int((adata_sub.obs["day"].astype(str)  == str(_d)).sum())
        print(f"      day {str(_d):<8s}  {_nf:>6,}  →  {_ns:>5,}")
    return adata_sub
def build_cropped_local_tmaps(
    adata_work:    sc.AnnData,
    src_tmaps_dir: Path,
    out_dir:       Path,
) -> Path:
    """
    Crop every tmap file in *src_tmaps_dir* to exactly the cells present in
    *adata_work* at the corresponding timepoints, and write results to *out_dir*.
    WHY THIS IS NECESSARY
    ---------------------
    RealTimeKernel.from_wot() does NOT silently filter coupling matrices.
    It requires that for each (src_day, tgt_day) pair:
      - tmap.obs_names == adata.obs_names at src_day   (exact match, same order)
      - tmap.var_names == adata.obs_names at tgt_day   (exact match, same order)
    Passing original full-dataset tmaps after subsampling causes:
      IndexError: Source observations for (...) don't match with adata.obs_names
    MATHEMATICAL NOTE
    -----------------
    Cropping a WOT transport plan to a row/column subset is a valid restriction
    of the original optimal transport plan to the sampled marginals.  No
    biological assumption of the WOT extension design is altered.
    ORDERING
    --------
    obs / var in each cropped file follow the order cells appear in adata_work,
    ensuring index-based alignment with from_wot() is consistent throughout.
    VALIDATION
    ----------
    After writing, each file is read back (backed="r") and obs_names / var_names
    are verified to exactly match the expected barcode lists.  Raises on mismatch.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    day_col = adata_work.obs["day"].astype(str)
    # day-string → ordered barcode list for cells in adata_work at that day
    day_to_barcodes: dict[str, list[str]] = {}
    for _d in adata_work.obs["day"].cat.categories:
        _mask = day_col == str(_d)
        day_to_barcodes[str(_d)] = adata_work.obs_names[_mask].tolist()
    tmap_files = sorted(src_tmaps_dir.glob("tmaps_*.h5ad"))
    if not tmap_files:
        raise FileNotFoundError(f"No tmaps_*.h5ad files in {src_tmaps_dir}")
    print(f"\n  [CROP TMAPS] {len(tmap_files)} files  →  {out_dir.name}/")
    for _fp in tmap_files:
        _m = re.search(r"tmaps_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)\.h5ad$", _fp.name)
        if not _m:
            print(f"    [SKIP] Unrecognised filename pattern: {_fp.name}")
            continue
        src_day_str = _m.group(1)
        tgt_day_str = _m.group(2)
        src_bcs = day_to_barcodes.get(src_day_str, [])
        tgt_bcs = day_to_barcodes.get(tgt_day_str, [])
        if len(src_bcs) == 0:
            raise ValueError(
                f"Tmap pair ({src_day_str}, {tgt_day_str}): "
                f"no sampled cells at src day '{src_day_str}'.\n"
                f"Days present in adata_work: {list(day_to_barcodes.keys())}"
            )
        if len(tgt_bcs) == 0:
            raise ValueError(
                f"Tmap pair ({src_day_str}, {tgt_day_str}): "
                f"no sampled cells at tgt day '{tgt_day_str}'.\n"
                f"Days present in adata_work: {list(day_to_barcodes.keys())}"
            )
        _tm_full      = sc.read_h5ad(str(_fp))
        _full_src_set = set(_tm_full.obs_names)
        _full_tgt_set = set(_tm_full.var_names)
        _miss_src = [b for b in src_bcs if b not in _full_src_set]
        _miss_tgt = [b for b in tgt_bcs if b not in _full_tgt_set]
        if _miss_src:
            raise ValueError(
                f"Tmap {_fp.name}: {len(_miss_src)} sampled src-barcodes absent "
                f"from full tmap — subsample may contain cells outside the WOT run.\n"
                f"First 5: {_miss_src[:5]}"
            )
        if _miss_tgt:
            raise ValueError(
                f"Tmap {_fp.name}: {len(_miss_tgt)} sampled tgt-barcodes absent "
                f"from full tmap — subsample may contain cells outside the WOT run.\n"
                f"First 5: {_miss_tgt[:5]}"
            )
        _tm_crop = _tm_full[src_bcs, tgt_bcs].copy()
        if _tm_crop.n_obs == 0 or _tm_crop.n_vars == 0:
            raise ValueError(
                f"Cropped tmap for ({src_day_str}, {tgt_day_str}) is empty "
                f"(shape {_tm_crop.shape}).  Cannot continue."
            )
        _out_fp = out_dir / _fp.name
        _tm_crop.write_h5ad(str(_out_fp))
        # Exact-match post-write validation (backed="r" reads only index)
        _check = sc.read_h5ad(str(_out_fp), backed="r")
        _check_obs = _check.obs_names.tolist()
        _check_var = _check.var_names.tolist()
        _check.file.close()
        if _check_obs != src_bcs:
            raise RuntimeError(
                f"Post-write validation FAILED for {_fp.name}: obs_names mismatch.\n"
                f"Expected {len(src_bcs)} barcodes; got {len(_check_obs)}."
            )
        if _check_var != tgt_bcs:
            raise RuntimeError(
                f"Post-write validation FAILED for {_fp.name}: var_names mismatch.\n"
                f"Expected {len(tgt_bcs)} barcodes; got {len(_check_var)}."
            )
        print(
            f"    {_fp.name}  "
            f"{_tm_full.shape[0]:>6,}×{_tm_full.shape[1]:>6,}"
            f"  →  {_tm_crop.shape[0]:>5,}×{_tm_crop.shape[1]:>5,}  ✓"
        )
    print(f"  Cropped tmaps written to: {out_dir}")
    return out_dir


def resolve_var_name(var_names: pd.Index, aliases: tuple[str, ...]) -> str | None:
    upper_map = {str(v).upper(): str(v) for v in var_names}
    for alias in aliases:
        if alias.upper() in upper_map:
            return upper_map[alias.upper()]
    return None


def extract_gene_expression(adata: sc.AnnData, gene_name: str) -> np.ndarray:
    idx = adata.var_names.get_loc(gene_name)
    x = adata.X[:, idx]
    if sp.issparse(x):
        x = np.asarray(x.toarray()).ravel()
    else:
        x = np.asarray(x).ravel()
    return x.astype(np.float64)


def calc_expression_threshold(
    values: np.ndarray,
    quantile: float,
    minimum: float,
    floor: float | None = None,
) -> float:
    vals = np.asarray(values, dtype=np.float64)
    vals = vals[np.isfinite(vals)]
    vals = vals[vals > 0]
    if vals.size == 0:
        return float(minimum if floor is None else max(minimum, floor))
    thr = float(np.quantile(vals, quantile))
    if floor is not None:
        thr = max(thr, floor)
    return max(float(minimum), thr)


def classify_biological_ips_terminal_cells(
    adata: sc.AnnData,
    terminal_mask: np.ndarray,
    gene_aliases: dict[str, tuple[str, ...]],
    high_quantile: float,
    very_high_quantile: float,
    min_high_expr: float,
    min_very_high_expr: float,
    min_very_high_delta: float,
) -> tuple[pd.DataFrame, dict[str, str], dict[str, float], dict[str, np.ndarray], np.ndarray, pd.Series]:
    terminal_mask = np.asarray(terminal_mask, dtype=bool)
    if terminal_mask.sum() == 0:
        raise RuntimeError("No terminal cells available for biological iPSC classification.")

    resolved: dict[str, str] = {}
    expr: dict[str, np.ndarray] = {}
    thr_high: dict[str, float] = {}
    thr_very: dict[str, float] = {}

    for logical_name, aliases in gene_aliases.items():
        gene_name = resolve_var_name(adata.var_names, aliases)
        if gene_name is None:
            raise KeyError(
                f"Required biological marker gene for {logical_name} not found in var_names. "
                f"Tried aliases={aliases}"
            )
        resolved[logical_name] = gene_name
        expr_vec = extract_gene_expression(adata, gene_name)
        expr[logical_name] = expr_vec
        thr_high[logical_name] = calc_expression_threshold(
            expr_vec, quantile=high_quantile, minimum=min_high_expr
        )
        thr_very[logical_name] = calc_expression_threshold(
            expr_vec, quantile=very_high_quantile, minimum=min_very_high_expr,
            floor=thr_high[logical_name] + min_very_high_delta,
        )

    high_masks = {k: (expr[k] >= thr_high[k]) for k in expr}
    very_masks = {k: (expr[k] >= thr_very[k]) for k in expr}

    n_high = np.zeros(adata.n_obs, dtype=int)
    n_very = np.zeros(adata.n_obs, dtype=int)
    for k in expr:
        n_high += high_masks[k].astype(int)
        n_very += very_masks[k].astype(int)

    category = np.full(adata.n_obs, "non_terminal", dtype=object)
    category[terminal_mask] = "terminal_non_ips_bio"

    tri_very = terminal_mask & (n_very == 3)
    tri_high = terminal_mask & (n_high == 3) & (~tri_very)
    bi_high = terminal_mask & (n_high == 2)
    pou_only = terminal_mask & high_masks["POU5F1"] & (~high_masks["SOX2"]) & (~high_masks["NANOG"])
    single_high = terminal_mask & (n_high == 1) & (~pou_only)

    # Assign in descending priority to keep classes mutually exclusive.
    category[single_high] = "single_high"
    category[pou_only] = "POU5F1-only"
    category[bi_high] = "bi_high"
    category[tri_high] = "tri_high"
    category[tri_very] = "tri_very_high"

    bio_positive = np.isin(category, ["POU5F1-only", "single_high", "bi_high", "tri_high", "tri_very_high"])

    out_df = pd.DataFrame({
        "bio_ips_category": category,
        "bio_ips_terminal_cell": bio_positive.astype(bool),
    }, index=adata.obs_names)
    for logical_name in ["POU5F1", "SOX2", "NANOG"]:
        out_df[f"expr_{logical_name}"] = expr[logical_name]
        out_df[f"high_{logical_name}"] = high_masks[logical_name].astype(bool)
        out_df[f"very_high_{logical_name}"] = very_masks[logical_name].astype(bool)
        out_df[f"thr_high_{logical_name}"] = thr_high[logical_name]
        out_df[f"thr_very_high_{logical_name}"] = thr_very[logical_name]

    category_series = pd.Series(category, index=adata.obs_names, dtype="object")
    return out_df, resolved, {**{f"{k}_high": v for k, v in thr_high.items()}, **{f"{k}_very_high": v for k, v in thr_very.items()}}, expr, bio_positive, category_series


def _init_gpu_backend() -> tuple[object | None, bool]:
    """
    Try to initialize rapids-singlecell for GPU-accelerated preprocessing.

    Returns
    -------
    (rsc_module_or_none, gpu_enabled)
    """
    if not USE_GPU:
        print("  [GPU] USE_GPU=0 -> GPU acceleration disabled.")
        return None, False
    try:
        import rapids_singlecell as rsc  # type: ignore
        import cupy  # noqa: F401
        print("  [GPU] rapids-singlecell backend available.")
        return rsc, True
    except Exception as e:
        print(f"  [GPU] rapids-singlecell unavailable -> falling back to CPU. Details: {e}")
        return None, False


def ensure_graph_and_embedding(
    adata_work: sc.AnnData,
    *,
    n_neighbors: int,
    leiden_key: str,
    leiden_resolution: float,
    leiden_random_state: int,
    use_gpu: bool,
    rsc=None,
) -> bool:
    """
    Ensure PCA / neighbors / Leiden / UMAP exist.

    Returns
    -------
    used_gpu : bool
        Whether GPU acceleration was actually used.
    """
    used_gpu = False

    if use_gpu and rsc is not None:
        print("  [GPU] Moving AnnData to GPU memory for PCA / neighbors / Leiden / UMAP ...")
        rsc.get.anndata_to_GPU(adata_work)
        used_gpu = True
        try:
            need_pca = ("X_pca" not in adata_work.obsm) or GPU_FORCE_RECOMPUTE_PCA
            if need_pca:
                print(f"  [GPU] rsc.pp.pca(n_comps={GPU_PCA_COMPONENTS}) ...")
                rsc.pp.pca(adata_work, n_comps=GPU_PCA_COMPONENTS)
            else:
                print("  [GPU] Existing X_pca detected; reusing it.")

            if "connectivities" not in adata_work.obsp:
                print(f"  [GPU] rsc.pp.neighbors(n_neighbors={n_neighbors}) ...")
                if "X_pca" in adata_work.obsm:
                    _n_pcs = min(
                        GPU_PCA_COMPONENTS,
                        int(adata_work.obsm["X_pca"].shape[1]),
                    )
                    rsc.pp.neighbors(adata_work, n_neighbors=n_neighbors, n_pcs=_n_pcs)
                else:
                    rsc.pp.neighbors(adata_work, n_neighbors=n_neighbors)
            else:
                print("  [GPU] Existing connectivities detected; reusing graph.")

            if leiden_key not in adata_work.obs.columns:
                print(
                    f"  [GPU] rsc.tl.leiden(key_added='{leiden_key}', "
                    f"resolution={leiden_resolution}, random_state={leiden_random_state}) ..."
                )
                rsc.tl.leiden(
                    adata_work,
                    resolution=leiden_resolution,
                    key_added=leiden_key,
                    random_state=leiden_random_state,
                )
            else:
                print(f"  [GPU] Existing `{leiden_key}` annotation found; reusing it.")

            if "X_umap" not in adata_work.obsm:
                print("  [GPU] rsc.tl.umap() ...")
                rsc.tl.umap(adata_work)
            else:
                print("  [GPU] Existing X_umap detected; reusing it.")
        finally:
            print("  [GPU] Moving AnnData back to CPU memory for CellRank / matplotlib ...")
            rsc.get.anndata_to_CPU(adata_work)

        return used_gpu

    # CPU fallback
    if "connectivities" in adata_work.obsp:
        print("  kNN graph already present in adata_work.obsp.")
    else:
        if "X_pca" not in adata_work.obsm:
            print("  X_pca not found — computing PCA (50 components) …")
            sc.pp.pca(adata_work, n_comps=50)
        else:
            print(f"  X_pca found: shape {adata_work.obsm['X_pca'].shape}")
        print(f"  sc.pp.neighbors(n_neighbors={n_neighbors}) …")
        sc.pp.neighbors(adata_work, n_neighbors=n_neighbors, use_rep="X_pca")
        print("  Done.")

    if leiden_key in adata_work.obs.columns:
        print(f"  Existing `{leiden_key}` annotation found.")
    else:
        print(
            f"  sc.tl.leiden(key_added='{leiden_key}', resolution={leiden_resolution}, "
            f"random_state={leiden_random_state}) …"
        )
        sc.tl.leiden(
            adata_work,
            resolution=leiden_resolution,
            key_added=leiden_key,
            random_state=leiden_random_state,
        )
        print("  Leiden done.")

    if "X_umap" not in adata_work.obsm:
        print("  [INFO] X_umap not found — computing UMAP …")
        sc.tl.umap(adata_work)

    return used_gpu



# =============================================================================
# PIPELINE  (inside main — only runs when __name__ == "__main__")
# =============================================================================
def main() -> None:
    """
    Full CellRank2 WOT-extension pipeline.
    local_debug mode  → pipeline smoke test on a 32 GB Windows workstation.
                        Results are NOT intended as final biological conclusions.
    hpc_full mode     → final analysis on Shirokane server with full 75 k-cell
                        dataset and PETSc/SLEPc-backed GPCCA.
    """
    # ── Per-run globals ────────────────────────────────────────────────────────
    TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")
    T0        = _time.time()
    _SEP      = "=" * 78
    for _d in [FIGURES_DIR, METRICS_DIR, CR2_OUT_DIR]:
        _d.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir    = str(FIGURES_DIR)
    sc.settings.autoshow  = False
    sc.settings.verbosity = 2
    # ── CellRank import ────────────────────────────────────────────────────────
    # Imported here (inside main) so child processes that import this file as
    # a module do not trigger CellRank's module-level side effects.
    try:
        import cellrank as cr
        print(f"[INFO] CellRank version : {cr.__version__}")
    except ImportError as _e:
        print("[ERROR] CellRank2 not installed.")
        print("  conda activate traj_env && pip install cellrank")
        print(f"  Details: {_e}")
        sys.exit(1)

    rsc, _gpu_backend_ok = _init_gpu_backend()
    # ── Suppress noisy CellRank log messages ──────────────────────────────────
    # from_wot() emits repeated INFO-level messages like:
    #   "File 'tmaps_0.5_2.0.h5ad' has more than two extensions"
    # for every tmap whose name contains multiple dots (e.g. 0.5, 2.0).
    # These are harmless — the files load correctly.  We silence them here
    # without touching real ERROR/WARNING diagnostics.
    import logging as _logging
    _logging.getLogger("cellrank").setLevel(_logging.WARNING)
    # Also suppress via the warnings module for any variant emitted that way.
    warnings.filterwarnings("ignore", message=".*more than.*extension.*")
    # =========================================================================
    # STEP 1  —  Environment / backend check
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 1  —  Environment / backend check")
    print(_SEP)
    print(f"  RUN_MODE          : {RUN_MODE}  (tag={_MODE_TAG})")
    print(f"  FATE_PROB_N_JOBS  : {FATE_PROB_N_JOBS}  "
          f"({'serial — no child processes' if FATE_PROB_N_JOBS == 1 else 'parallel'})")
    try:
        import petsc4py  # noqa: F401
        _petsc_ok = True
    except ImportError:
        _petsc_ok = False
    try:
        import slepc4py  # noqa: F401
        _slepc_ok = True
    except ImportError:
        _slepc_ok = False
    print(f"  petsc4py          : {'available' if _petsc_ok else 'NOT available'}")
    print(f"  slepc4py          : {'available' if _slepc_ok else 'NOT available'}")
    if not (_petsc_ok and _slepc_ok):
        if RUN_MODE == "hpc_full":
            print(
                "\n  [WARNING] PETSc/SLEPc NOT available in hpc_full mode.\n"
                "  GPCCA.compute_schur() will fall back to method='brandts',\n"
                "  which densifies the full (~75 k × ~75 k) transition matrix.\n"
                "  Estimated peak RAM: ~42 GiB (float64).  Likely OOM on <64 GB.\n"
                "  Install before running hpc_full:\n"
                "    conda install -c conda-forge petsc4py slepc4py"
            )
        else:
            print(
                "\n  [INFO] PETSc/SLEPc not available — acceptable in local_debug.\n"
                "  Brandts dense fallback on ~8 k-cell subset needs only ~0.5 GB."
            )
    # =========================================================================
    # STEP 2  —  Load AnnData  (WOT fate results from script 02c)
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 2  —  Load AnnData  (WOT fate results)")
    print(_SEP)
    fates_h5ad = find_latest(METRICS_DIR, "*_GSE230659_wot_fates.h5ad")
    print(f"  Input : {fates_h5ad.name}")
    adata = sc.read_h5ad(str(fates_h5ad))
    print(f"  Shape : {adata.n_obs:,} cells  ×  {adata.n_vars:,} genes")
    print(f"  obs   : {list(adata.obs.columns[:20])}")
    print(f"  obsm  : {list(adata.obsm.keys())}")
    # ── Validate required columns ─────────────────────────────────────────────
    _req  = ["day", "stage", "p_iPSC"]
    _miss = [c for c in _req if c not in adata.obs.columns]
    if _miss:
        raise ValueError(f"Required obs columns missing: {_miss}")
    # ── Convert day to ordered Categorical  (Fix 1) ───────────────────────────
    # RealTimeKernel._read_from_adata() raises TypeError if day is not CategoricalDtype.
    adata.obs["day"] = (
        pd.to_numeric(adata.obs["day"], errors="coerce")
        .astype(float)
        .astype("category")
    )
    _day_cats = sorted(adata.obs["day"].cat.categories.tolist())
    print(f"\n  Timepoints (day categories) : {_day_cats}")
    _stage_cts = adata.obs["stage"].value_counts().to_dict()
    print(f"  Stage distribution          : {_stage_cts}")
    # ── Dynamically detect the iPSC stage label ───────────────────────────────
    # Accepts any label variant: "hCiPSC", "hCiPSCs", "iPSC", "hcipsc", etc.
    # Picks the most-frequent matching label for stability across batches.
    _all_stages     = sorted(adata.obs["stage"].unique())
    _ips_candidates = [s for s in _all_stages if "ips" in str(s).lower()]
    if _ips_candidates:
        _ips_stage_label = max(_ips_candidates, key=lambda s: _stage_cts.get(s, 0))
        print(f"  Detected iPSC stage label   : '{_ips_stage_label}'")
        print(f"  All iPSC-like candidates    : {_ips_candidates}")
    else:
        _ips_stage_label = None
        print(f"  [WARN] No stage label containing 'ips' found.")
        print(f"  Available stages: {_all_stages}")
        print("  Terminal-state selection will rely on GPCCA macrostate names only.")
    # =========================================================================
    # STEP 3  —  Locate WOT transport maps  (original full-dataset)
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 3  —  Locate WOT transport maps  (original full-dataset)")
    print(_SEP)
    _ts        = extract_ts(fates_h5ad.name)
    tmaps_dir  = find_tmaps_dir(TMAPS_BASE, _ts)
    _tmap_files = sorted(tmaps_dir.glob("tmaps_*.h5ad"))
    print(f"  Transport-map dir   : {tmaps_dir}")
    print(f"  Files found         : {len(_tmap_files)}")
    _tmap_pairs: list[tuple[float, float]] = []
    for _fp in _tmap_files:
        _m = re.search(r"tmaps_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)\.h5ad$", _fp.name)
        if _m:
            _tmap_pairs.append((float(_m.group(1)), float(_m.group(2))))
    print(f"  Pairs in dir        : {sorted(_tmap_pairs)}")
    _expected = [(float(_day_cats[i]), float(_day_cats[i + 1]))
                 for i in range(len(_day_cats) - 1)]
    _missing  = [p for p in _expected if p not in _tmap_pairs]
    if _missing:
        print(f"  [WARN] Missing tmap pairs: {_missing}")
    # Guard: every day in adata must appear in at least one tmap
    _tmap_days  = set()
    for _src, _tgt in _tmap_pairs:
        _tmap_days.add(str(_src))
        _tmap_days.add(str(_tgt))
    _adata_days         = {str(float(d)) for d in _day_cats}
    _days_not_in_tmaps  = _adata_days - _tmap_days
    if _days_not_in_tmaps:
        raise ValueError(
            f"Days in adata not represented in any tmap: {_days_not_in_tmaps}\n"
            f"Tmap days: {sorted(_tmap_days)}  |  Adata days: {sorted(_adata_days)}"
        )
    # =========================================================================
    # STEP 4  —  Prepare working AnnData  (subsample in local_debug)
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 4  —  Prepare working AnnData")
    print(_SEP)
    print(f"  Mode : {RUN_MODE}")
    if RUN_MODE == "local_debug":
        adata.obs["_day_str"] = adata.obs["day"].astype(str)
        adata_work = stratified_subsample(
            adata_full       = adata,
            target_n         = LOCAL_TARGET_CELLS,
            min_per_stratum  = LOCAL_MIN_PER_STRATUM,
            strat_cols       = ["stage", "_day_str"],
            seed             = LOCAL_SEED,
        )
        # Re-apply ordered categorical day in the copied subset
        adata_work.obs["day"] = pd.Categorical(
            pd.to_numeric(adata_work.obs["day"], errors="coerce").astype(float),
            categories=sorted(adata_work.obs["day"].astype(float).unique()),
            ordered=True,
        )
        _day_cats_work = sorted(adata_work.obs["day"].cat.categories.tolist())
        print(f"\n  Working AnnData timepoints : {_day_cats_work}")
        print(f"  Working AnnData shape      : {adata_work.n_obs:,} × {adata_work.n_vars:,}")
        # Guard: all days in adata_work must be covered by original tmaps
        _work_days = {str(float(d)) for d in _day_cats_work}
        _uncovered = _work_days - _tmap_days
        if _uncovered:
            raise ValueError(
                f"Sampled AnnData contains days not in any tmap: {_uncovered}\n"
                "This should not happen — the subsample was drawn from outside the "
                "WOT run?"
            )
    else:   # hpc_full
        adata_work     = adata   # alias — no copy, no extra memory
        _day_cats_work = _day_cats
        print(f"  Full dataset : {adata_work.n_obs:,} cells (no subsampling)")
    # =========================================================================
    # STEP 5  —  Prepare working tmaps directory
    # =========================================================================
    # local_debug: crop tmap files to exactly the sampled cells, write to a new
    #              directory, then point from_wot() at that directory.
    # hpc_full:    use original tmaps directory unchanged.
    #
    # WHY CROPPING IS DONE HERE (not after from_wot):
    # from_wot() requires exact barcode alignment.  We crop before calling it.
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 5  —  Prepare tmaps directory for working AnnData")
    print(_SEP)
    if RUN_MODE == "local_debug":
        _cropped_dir_name = f"{TIMESTAMP}_tmaps_localdebug_{adata_work.n_obs}"
        tmaps_dir_work    = TMAPS_BASE / _cropped_dir_name
        print(f"  local_debug: cropping tmaps to {adata_work.n_obs:,} sampled cells")
        print(f"  Cropped dir : {tmaps_dir_work}")
        build_cropped_local_tmaps(
            adata_work    = adata_work,
            src_tmaps_dir = tmaps_dir,
            out_dir       = tmaps_dir_work,
        )
    else:
        tmaps_dir_work = tmaps_dir
        print(f"  hpc_full: using original tmaps dir")
        print(f"  Tmaps dir   : {tmaps_dir_work}")
    _tmap_files_work = sorted(tmaps_dir_work.glob("tmaps_*.h5ad"))
    print(f"  Files in working tmaps dir : {len(_tmap_files_work)}")
    # =========================================================================
    # STEP 6  —  Ensure PCA / kNN graph / Leiden / UMAP
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 6  —  Ensure PCA / kNN graph / Leiden / UMAP")
    print(_SEP)

    _used_gpu_preproc = ensure_graph_and_embedding(
        adata_work,
        n_neighbors=N_NEIGHBORS,
        leiden_key=LEIDEN_KEY,
        leiden_resolution=LEIDEN_RESOLUTION,
        leiden_random_state=LEIDEN_RANDOM_STATE,
        use_gpu=_gpu_backend_ok,
        rsc=rsc,
    )

    _leiden_counts = adata_work.obs[LEIDEN_KEY].value_counts().sort_index()
    print(f"  Leiden clusters ({len(_leiden_counts)}):")
    print("    " + ", ".join([f"{k}:{v}" for k, v in _leiden_counts.items()]))
    print(f"  Preprocessing backend used: {'GPU' if _used_gpu_preproc else 'CPU'}")

    # =========================================================================
    # STEP 7  —  Build RealTimeKernel  (WOT extension design)
    # =========================================================================
    # Analysis design note:
    # RTK is built from the WOT transport maps (same as script 02c).
    # This is deliberate: we compare classical WOT fate probabilities against
    # CellRank2-extended fate probabilities that add within-timepoint kNN
    # connectivity on top of the same WOT plans.
    # RTK = WOT extension, NOT an independent CellRank2 model.
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 7  —  Build RealTimeKernel  (from_wot)")
    print(_SEP)
    print(f"  tmaps dir  : {tmaps_dir_work}")
    print(f"  time_key   : 'day'  (dtype={adata_work.obs['day'].dtype})")
    print(f"  n_cells    : {adata_work.n_obs:,}")
    try:
        rtk = cr.kernels.RealTimeKernel.from_wot(
            adata_work,
            path=str(tmaps_dir_work),
            time_key="day",
        )
        print(f"  RealTimeKernel created: {type(rtk).__name__}")
        print(f"  Couplings loaded: {list(rtk.couplings.keys())}")
    except Exception as _e:
        print(f"  [FATAL] from_wot failed: {_e}")
        sys.exit(1)
    # ── compute_transition_matrix  (Fix 2) ────────────────────────────────────
    # Verified parameters from cellrank._real_time_kernel.py:
    #   threshold        = "auto" | "auto_local" | float ∈ [0, 100]
    #   self_transitions = "all" → 80/20 WOT/connectivity at every timepoint
    #   conn_weight      = float in (0, 1)
    print(f"\n  compute_transition_matrix("
          f"threshold='{THRESHOLD}', self_transitions='all', "
          f"conn_weight={CONN_WEIGHT}) …")
    _t0 = _time.time()
    try:
        rtk.compute_transition_matrix(
            threshold=THRESHOLD,
            self_transitions="all",
            conn_weight=CONN_WEIGHT,
            n_neighbors=N_NEIGHBORS,
        )
    except TypeError as _e:
        print(f"  [WARN] Retrying with conn_kwargs dict: {_e}")
        rtk.compute_transition_matrix(
            threshold=THRESHOLD,
            self_transitions="all",
            conn_weight=CONN_WEIGHT,
            conn_kwargs={"n_neighbors": N_NEIGHBORS},
        )
    print(f"  Done in {_time.time() - _t0:.1f} s")
    print(f"  Transition matrix  : {rtk.transition_matrix.shape}")
    # =========================================================================
    # STEP 8  —  GPCCA: Schur decomposition
    # =========================================================================
    # MEMORY NOTE:
    #   With PETSc/SLEPc  → sparse solver; memory ∝ n_cells × n_components
    #   Without PETSc/SLEPc → brandts dense fallback; memory ∝ n_cells²
    #     local_debug (~8 k): ~0.5 GB → fine on 32 GB
    #     hpc_full  (~75 k): ~42 GiB → requires ≥64 GB RAM or PETSc/SLEPc
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 8  —  GPCCA: Schur decomposition")
    print(_SEP)
    print(f"  Matrix : {rtk.transition_matrix.shape[0]:,} × {rtk.transition_matrix.shape[1]:,}")
    print(f"  Backend: {'PETSc/SLEPc' if (_petsc_ok and _slepc_ok) else 'brandts (dense fallback)'}")
    if RUN_MODE == "hpc_full" and not (_petsc_ok and _slepc_ok):
        print(
            f"  [WARNING] hpc_full without PETSc/SLEPc: dense Schur on "
            f"{rtk.transition_matrix.shape[0]:,} cells may exhaust RAM.  Proceeding."
        )
    gpcca = cr.estimators.GPCCA(rtk)
    _t0   = _time.time()
    try:
        gpcca.compute_schur(n_components=N_SCHUR_COMPONENTS)
        print(f"  compute_schur(n_components={N_SCHUR_COMPONENTS}) done in "
              f"{_time.time() - _t0:.1f} s")
    except Exception as _e:
        print(f"  [WARN] n_components={N_SCHUR_COMPONENTS} failed: {_e}")
        print("  Retrying with n_components=10 …")
        gpcca.compute_schur(n_components=10)
        print(f"  Done in {_time.time() - _t0:.1f} s")
    # Eigenspectrum plot — confirmed pattern from CellRank2 reproducibility repo:
    #   estimator.plot_spectrum(real_only=True)   ← always real_only=True, no ax=
    # Passing ax= to the complex-eigenvalue scatter variant triggers:
    #   "Axes.scatter() got multiple values for argument 'ax'"
    # because CellRank passes ax positionally AND as a keyword internally.
    # Fix: use real_only=True (bar-chart; no scatter kwargs conflict), no ax=.
    # Capture the figure via plt.gcf() after the call.
    try:
        gpcca.plot_spectrum(real_only=True)
        _fig_eig = plt.gcf()
        _fig_eig.suptitle("GPCCA eigenspectrum (real part)", fontsize=11,
                           y=1.02)
        savefig(_fig_eig, f"{TIMESTAMP}_{_MODE_TAG}_cr2_eigenspectrum", FIGURES_DIR)
    except Exception as _e:
        print(f"  [WARN] Eigenspectrum plot skipped: {_e}")
        plt.close("all")
    # =========================================================================
    # STEP 9  —  GPCCA: compute macrostates
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 9  —  GPCCA: compute macrostates")
    print(_SEP)
    print(f"  n_states={N_MACROSTATES},  cluster_key='{LEIDEN_KEY}'")
    _t0 = _time.time()
    try:
        gpcca.compute_macrostates(
            n_states=N_MACROSTATES,
            cluster_key=LEIDEN_KEY,
            n_cells=TERMINAL_N_CELLS,
        )
        print(f"  Done in {_time.time() - _t0:.1f} s")
    except Exception as _e:
        print(f"  [WARN] n_states={N_MACROSTATES} failed: {_e}. Trying 4 …")
        gpcca.compute_macrostates(n_states=4, cluster_key=LEIDEN_KEY)
        print(f"  Done in {_time.time() - _t0:.1f} s")
    _ms_obs_key = "macrostates_fwd"
    _ms_names   = (
        list(gpcca.macrostates.cat.categories)
        if gpcca.macrostates is not None else []
    )
    print(f"  Macrostates found : {_ms_names}")
    _has_umap = "X_umap" in adata_work.obsm
    try:
        if not _has_umap:
            print("  [INFO] X_umap not found — computing UMAP …")
            sc.tl.umap(adata_work)
        _fig_ms, _ax_ms = plt.subplots(figsize=(7, 6))
        sc.pl.umap(adata_work, color=_ms_obs_key,
                   ax=_ax_ms, show=False, size=UMAP_SIZE, alpha=UMAP_ALPHA,
                   title=f"CellRank2 macrostates (cluster_key={LEIDEN_KEY})")
        savefig(_fig_ms, f"{TIMESTAMP}_{_MODE_TAG}_cr2_macrostates_umap", FIGURES_DIR)
        try:
            gpcca.plot_macrostate_composition(key="stage", show=False)
            _fig_comp = plt.gcf()
            _fig_comp.suptitle("Macrostate composition by stage", fontsize=11, y=1.02)
            savefig(_fig_comp, f"{TIMESTAMP}_{_MODE_TAG}_cr2_macrostate_composition_stage", FIGURES_DIR)
        except Exception as _e_comp:
            print(f"  [WARN] Macrostate composition plot skipped: {_e_comp}")
    except Exception as _e:
        print(f"  [WARN] Macrostates UMAP skipped: {_e}")

    # =========================================================================
    # STEP 10  —  GPCCA: predict + confirm terminal states
    # =========================================================================
    # v12 logic:
    #   1. auto-predict terminal states from coarse_T
    #   2. summarise each macrostate using stage / fate_ips / p_iPSC / late-time evidence
    #   3. treat auto terminal states as default terminal set
    #   4. ONLY if the auto set does not already contain a convincing iPSC endpoint,
    #      manually add at most a very small number of strong late candidates
    #   5. Method 1 (dynamics): define dynamic iPSC terminal CELLS from CellRank-selected
    #      iPSC terminal lineages among the final terminal-state set
    #   6. Method 2 (biology): among terminal cells, define mutually exclusive classes
    #      POU5F1-only / single_high / bi_high / tri_high / tri_very_high using expression
    #      thresholds for POU5F1, SOX2 and NANOG. Final iPSC detection uses method 2.
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 10  —  GPCCA: predict + confirm terminal states")
    print(_SEP)
    if _ms_obs_key not in adata_work.obs.columns and gpcca.macrostates is not None:
        adata_work.obs[_ms_obs_key] = gpcca.macrostates

    _auto_terminal_states: list[str] = []
    print(f"  Auto-predict terminal states: method='{TERMINAL_METHOD}', "
          f"stability_threshold={TERMINAL_STABILITY_THRESHOLD}")
    try:
        gpcca.predict_terminal_states(
            method=TERMINAL_METHOD,
            n_cells=TERMINAL_N_CELLS,
            stability_threshold=TERMINAL_STABILITY_THRESHOLD,
        )
        _auto_terminal_states = [str(x) for x in list(gpcca.terminal_states.cat.categories)]
        print(f"  Auto terminal states ({len(_auto_terminal_states)}): {_auto_terminal_states}")
    except Exception as _e:
        print(f"  [WARN] predict_terminal_states() failed: {_e}")

    if not _auto_terminal_states:
        print("  [WARN] Auto terminal-state detection returned nothing.")
        print("  Falling back to top_n terminal-state prediction (n_states=2) …")
        gpcca.predict_terminal_states(
            method="top_n",
            n_states=min(2, max(1, len(_ms_names))),
            n_cells=TERMINAL_N_CELLS,
        )
        _auto_terminal_states = [str(x) for x in list(gpcca.terminal_states.cat.categories)]
        print(f"  Fallback terminal states: {_auto_terminal_states}")

    _day_numeric = adata_work.obs["day"].astype(float).values
    _last_day    = float(np.max(_day_numeric))

    _term_rows: list[dict] = []
    for _ms in _ms_names:
        _ms_str = str(_ms)
        _mask_ms = (adata_work.obs[_ms_obs_key] == _ms).values
        _n_ms    = int(_mask_ms.sum())
        if _n_ms == 0:
            continue

        _stage_vals = adata_work.obs.loc[_mask_ms, "stage"].astype(str)
        _dominant_stage = _stage_vals.mode().iloc[0] if len(_stage_vals) else "NA"
        _dominant_stage_frac = float((_stage_vals == _dominant_stage).mean()) if len(_stage_vals) else np.nan

        _leiden_vals = adata_work.obs.loc[_mask_ms, LEIDEN_KEY].astype(str)
        _dominant_leiden = _leiden_vals.mode().iloc[0] if len(_leiden_vals) else "NA"
        _dominant_leiden_frac = float((_leiden_vals == _dominant_leiden).mean()) if len(_leiden_vals) else np.nan

        _frac_ips_stage = (
            float(np.mean((adata_work.obs.loc[_mask_ms, "stage"] == _ips_stage_label).values))
            if _ips_stage_label is not None else np.nan
        )
        _mean_fate_ips = (
            float(np.nanmean(adata_work.obs.loc[_mask_ms, "fate_ips"].astype(float)))
            if "fate_ips" in adata_work.obs.columns else np.nan
        )
        _mean_p_iPSC = float(np.nanmean(adata_work.obs.loc[_mask_ms, "p_iPSC"].astype(float)))
        _frac_last_day = float(np.mean(_day_numeric[_mask_ms] == _last_day))
        _mean_day = float(np.nanmean(_day_numeric[_mask_ms]))

        _ev_stage = bool(np.isfinite(_frac_ips_stage) and _frac_ips_stage >= IPS_TERMINAL_STAGE_FRAC)
        _ev_fate  = bool(np.isfinite(_mean_fate_ips) and _mean_fate_ips >= IPS_TERMINAL_FATE_IPS_FRAC)
        _ev_pips  = bool(np.isfinite(_mean_p_iPSC) and _mean_p_iPSC >= IPS_TERMINAL_PIPSC_THRESHOLD)
        _ev_late  = bool(np.isfinite(_frac_last_day) and _frac_last_day >= IPS_TERMINAL_LAST_DAY_FRAC)
        _evidence_count = int(_ev_stage) + int(_ev_fate) + int(_ev_pips)
        _evidence_score = (
            (0.35 * (0.0 if not np.isfinite(_frac_ips_stage) else _frac_ips_stage)) +
            (0.35 * (0.0 if not np.isfinite(_mean_fate_ips) else _mean_fate_ips)) +
            (0.25 * (0.0 if not np.isfinite(_mean_p_iPSC) else _mean_p_iPSC)) +
            (0.05 * (0.0 if not np.isfinite(_frac_last_day) else _frac_last_day))
        )

        _term_rows.append({
            "macrostate": _ms_str,
            "n_assigned_cells": _n_ms,
            "is_auto_terminal": _ms_str in _auto_terminal_states,
            "mean_day": _mean_day,
            "frac_last_day": _frac_last_day,
            "mean_p_iPSC": _mean_p_iPSC,
            "frac_ips_stage": _frac_ips_stage,
            "mean_fate_ips": _mean_fate_ips,
            "dominant_stage": _dominant_stage,
            "dominant_stage_frac": _dominant_stage_frac,
            "dominant_leiden": _dominant_leiden,
            "dominant_leiden_frac": _dominant_leiden_frac,
            "evidence_stage": _ev_stage,
            "evidence_fate_ips": _ev_fate,
            "evidence_p_iPSC": _ev_pips,
            "evidence_late": _ev_late,
            "evidence_count": _evidence_count,
            "evidence_score": _evidence_score,
        })

    _term_summary_df = pd.DataFrame(_term_rows)
    if _term_summary_df.empty:
        raise RuntimeError("Macrostate summary table is empty — cannot continue.")

    _term_summary_df["manual_added_terminal"] = False
    _term_summary_df["is_final_terminal"]     = _term_summary_df["is_auto_terminal"]
    _term_summary_df["is_ips_terminal"]       = False
    _term_summary_df["ips_selection_rule"]    = ""

    print("\nMacrostate summary:")
    print(
        _term_summary_df.sort_values(
            ["is_auto_terminal", "evidence_score", "frac_last_day", "mean_p_iPSC"],
            ascending=[False, False, False, False],
        ).to_string(index=False)
    )

    _auto_ips_mask = (
        _term_summary_df["is_auto_terminal"] &
        (_term_summary_df["n_assigned_cells"] >= MIN_MACROSTATE_CELLS) &
        _term_summary_df["evidence_late"] &
        (_term_summary_df["evidence_count"] >= IPS_TERMINAL_MIN_EVIDENCE_COUNT)
    )
    _n_auto_ips_like = int(_auto_ips_mask.sum())
    print(f"\n  Auto terminal states with convincing iPSC evidence: {_n_auto_ips_like}")

    _manual_add_states: list[str] = []
    if _n_auto_ips_like == 0:
        _manual_mask = (
            (~_term_summary_df["is_auto_terminal"]) &
            (_term_summary_df["n_assigned_cells"] >= MIN_MACROSTATE_CELLS) &
            (_term_summary_df["frac_last_day"] >= MANUAL_ADD_LAST_DAY_FRAC) &
            (
                (_term_summary_df["frac_ips_stage"] >= MANUAL_ADD_IPS_STAGE_FRAC).fillna(False).astype(bool) +
                (_term_summary_df["mean_fate_ips"] >= MANUAL_ADD_FATE_IPS_FRAC).fillna(False).astype(bool) +
                (_term_summary_df["mean_p_iPSC"] >= MANUAL_ADD_PIPSC_THRESHOLD).fillna(False).astype(bool)
            >= MANUAL_ADD_MIN_EVIDENCE_COUNT)
        )
        if _manual_mask.any():
            _manual_candidates = (
                _term_summary_df.loc[_manual_mask]
                .sort_values(["evidence_score", "mean_p_iPSC", "frac_last_day"], ascending=False)
                .head(MANUAL_ADD_MAX_STATES)
            )
            _manual_add_states = _manual_candidates["macrostate"].tolist()
            _term_summary_df.loc[
                _term_summary_df["macrostate"].isin(_manual_add_states),
                "manual_added_terminal"
            ] = True
            print("\n  Manual supplement — auto set lacks convincing iPSC endpoint; adding:")
            for _ms in _manual_add_states:
                print(f"    • {_ms}")
        else:
            print("\n  Manual supplement — no strong missed iPSC-like late macrostate found.")
    else:
        print("\n  Manual supplement skipped — auto terminal set already contains convincing iPSC endpoint.")

    _final_terminal_states = list(dict.fromkeys(_auto_terminal_states + _manual_add_states))
    if not _final_terminal_states:
        raise RuntimeError("No terminal states available after auto + manual selection.")

    gpcca.set_terminal_states(states=_final_terminal_states)
    _ts_names = [str(x) for x in list(gpcca.terminal_states.cat.categories)]
    _term_summary_df["is_final_terminal"] = _term_summary_df["macrostate"].isin(_ts_names)
    print(f"\n  Final terminal states ({len(_ts_names)}): {_ts_names}")

    _ips_candidates = _term_summary_df[
        _term_summary_df["is_final_terminal"] &
        (_term_summary_df["n_assigned_cells"] >= MIN_MACROSTATE_CELLS) &
        _term_summary_df["evidence_late"] &
        (_term_summary_df["evidence_count"] >= IPS_TERMINAL_MIN_EVIDENCE_COUNT)
    ].copy()

    if _ips_candidates.empty:
        _cand = _term_summary_df[_term_summary_df["is_final_terminal"]].copy()
        if _cand.empty:
            raise RuntimeError("Final terminal-state set is empty after set_terminal_states().")
        _cand = _cand.sort_values(["evidence_score", "mean_p_iPSC", "frac_last_day"], ascending=False)
        _selected_lineages = _cand.head(1)["macrostate"].tolist()
        _term_summary_df.loc[
            _term_summary_df["macrostate"].isin(_selected_lineages),
            "is_ips_terminal"
        ] = True
        _term_summary_df.loc[
            _term_summary_df["macrostate"].isin(_selected_lineages),
            "ips_selection_rule"
        ] = "fallback: top-1 final terminal by evidence_score"
    else:
        _ips_candidates = _ips_candidates.sort_values(
            ["evidence_score", "mean_p_iPSC", "frac_last_day"],
            ascending=False
        ).head(IPS_TERMINAL_MAX_STATES)
        _selected_lineages = _ips_candidates["macrostate"].tolist()
        _term_summary_df.loc[
            _term_summary_df["macrostate"].isin(_selected_lineages),
            "is_ips_terminal"
        ] = True
        _term_summary_df.loc[
            _term_summary_df["macrostate"].isin(_selected_lineages),
            "ips_selection_rule"
        ] = (
            f"strict evidence: final terminal + last_day>={IPS_TERMINAL_LAST_DAY_FRAC} "
            f"+ at least {IPS_TERMINAL_MIN_EVIDENCE_COUNT}/3 of "
            f"(frac_ips_stage, mean_fate_ips, mean_p_iPSC); top-{IPS_TERMINAL_MAX_STATES} by evidence_score"
        )

    print(f"\n  Dynamic method — iPSC terminal lineage(s) ({len(_selected_lineages)}): {_selected_lineages}")
    for _, _row in _term_summary_df[_term_summary_df["is_ips_terminal"]].iterrows():
        print(
            f"    • {_row['macrostate']:<25s}  "
            f"late={_row['frac_last_day']:.3f}  "
            f"frac_ips_stage={_row['frac_ips_stage']:.3f}  "
            f"mean_fate_ips={_row['mean_fate_ips']:.3f}  "
            f"mean_p_iPSC={_row['mean_p_iPSC']:.3f}  "
            f"evidence_count={int(_row['evidence_count'])}  "
            f"score={_row['evidence_score']:.3f}  "
            f"rule={_row['ips_selection_rule']}"
        )

    _term_summary_path = (
        CR2_OUT_DIR / (TIMESTAMP + f"_{_MODE_TAG}_cr2_terminal_state_summary.tsv")
    )

    _ts_obs_key: str | None = None
    for _k_cand in ("term_states_fwd", "terminal_states_fwd"):
        if _k_cand in adata_work.obs.columns:
            _ts_obs_key = _k_cand
            break
    if _ts_obs_key is None:
        _fwd_cols = [
            c for c in adata_work.obs.columns
            if (c.endswith("_states_fwd") and "initial" not in c.lower())
        ]
        _ts_obs_key = _fwd_cols[0] if _fwd_cols else None
    print(f"  Terminal-states obs key detected: '{_ts_obs_key}'")

    if _ts_obs_key is not None:
        try:
            _fig_ts, _ax_ts = plt.subplots(figsize=(7, 6))
            sc.pl.umap(
                adata_work,
                color=_ts_obs_key,
                ax=_ax_ts,
                show=False,
                size=UMAP_SIZE,
                alpha=UMAP_ALPHA,
                title=f"CellRank2 terminal states [{_ts_obs_key}]",
            )
            savefig(_fig_ts, f"{TIMESTAMP}_{_MODE_TAG}_cr2_terminal_states_umap", FIGURES_DIR)
        except Exception as _e:
            print(f"  [WARN] Terminal states UMAP skipped: {_e}")
    else:
        print("  [WARN] Terminal-states obs key not found in adata_work.obs — UMAP skipped.")

    # ── v12: cell-level dynamic-vs-biology iPSC endpoint comparison ──────────
    if _ts_obs_key is None:
        raise RuntimeError("Terminal-states obs key not found — cannot define terminal-cell universe.")

    _ts_series = adata_work.obs[_ts_obs_key].astype(object)
    _terminal_cell_mask = _ts_series.isin(_ts_names).values
    _dyn_ips_terminal_mask = _ts_series.isin(_selected_lineages).values
    adata_work.obs["dyn_ips_terminal_cell"] = _dyn_ips_terminal_mask
    adata_work.obs["terminal_cell"] = _terminal_cell_mask
    adata_work.obs["dyn_ips_terminal_label"] = np.where(
        _dyn_ips_terminal_mask,
        "dynamic_ips_terminal",
        np.where(_terminal_cell_mask, "dynamic_terminal_non_ips", "non_terminal")
    )

    print(f"\n  Terminal-cell universe: {_terminal_cell_mask.sum():,} cells")
    print(f"  Dynamic method positive cells: {_dyn_ips_terminal_mask.sum():,}")

    _bio_df, _bio_resolved_genes, _bio_thresholds, _bio_expr, _bio_positive_mask, _bio_category = (
        classify_biological_ips_terminal_cells(
            adata=adata_work,
            terminal_mask=_terminal_cell_mask,
            gene_aliases=BIO_MARKER_ALIASES,
            high_quantile=BIO_HIGH_QUANTILE,
            very_high_quantile=BIO_VERY_HIGH_QUANTILE,
            min_high_expr=BIO_MIN_HIGH_EXPR,
            min_very_high_expr=BIO_MIN_VERY_HIGH_EXPR,
            min_very_high_delta=BIO_MIN_VERY_HIGH_DELTA,
        )
    )
    for _col in _bio_df.columns:
        adata_work.obs[_col] = _bio_df[_col]

    # Final detection uses the biological method.
    adata_work.obs["ips_final_mask"] = adata_work.obs["bio_ips_terminal_cell"].astype(bool)
    adata_work.obs["ips_final_method"] = np.where(
        adata_work.obs["bio_ips_terminal_cell"].astype(bool),
        adata_work.obs["bio_ips_category"].astype(str),
        np.where(_terminal_cell_mask, "terminal_non_ips_bio", "non_terminal")
    )

    _bio_counts = (
        adata_work.obs.loc[_terminal_cell_mask, "bio_ips_category"]
        .value_counts()
        .reindex(["POU5F1-only", "single_high", "bi_high", "tri_high", "tri_very_high", "terminal_non_ips_bio"], fill_value=0)
    )
    print("\n  Biological method — terminal-cell categories (mutually exclusive):")
    for _cat, _ncat in _bio_counts.items():
        print(f"    • {_cat:<20s} {_ncat:>6,}")

    _overlap_mask = _dyn_ips_terminal_mask & _bio_positive_mask
    _dyn_only_mask = _dyn_ips_terminal_mask & (~_bio_positive_mask)
    _bio_only_mask = (~_dyn_ips_terminal_mask) & _bio_positive_mask & _terminal_cell_mask
    _neither_mask = _terminal_cell_mask & (~_dyn_ips_terminal_mask) & (~_bio_positive_mask)

    _n_dyn = int(_dyn_ips_terminal_mask.sum())
    _n_bio = int((_bio_positive_mask & _terminal_cell_mask).sum())
    _n_overlap = int(_overlap_mask.sum())
    _n_union = int((_dyn_ips_terminal_mask | (_bio_positive_mask & _terminal_cell_mask)).sum())
    _jaccard = (_n_overlap / _n_union) if _n_union > 0 else np.nan
    _dyn_in_bio = (_n_overlap / _n_dyn) if _n_dyn > 0 else np.nan
    _bio_in_dyn = (_n_overlap / _n_bio) if _n_bio > 0 else np.nan

    adata_work.obs["dyn_bio_overlap_label"] = np.where(
        _overlap_mask, "overlap_both",
        np.where(_dyn_only_mask, "dynamic_only",
                 np.where(_bio_only_mask, "biology_only",
                          np.where(_terminal_cell_mask, "terminal_neither", "non_terminal")))
    )

    print("\n  Dynamic vs biological iPSC terminal-cell overlap:")
    print(f"    Dynamic positives : {_n_dyn:,}")
    print(f"    Biology positives : {_n_bio:,}")
    print(f"    Overlap           : {_n_overlap:,}")
    print(f"    Jaccard           : {_jaccard:.4f}" if np.isfinite(_jaccard) else "    Jaccard           : NA")
    print(f"    Overlap / Dynamic : {_dyn_in_bio:.4f}" if np.isfinite(_dyn_in_bio) else "    Overlap / Dynamic : NA")
    print(f"    Overlap / Biology : {_bio_in_dyn:.4f}" if np.isfinite(_bio_in_dyn) else "    Overlap / Biology : NA")
    print("  Final iPSC detection method: biology-based marker classes (method 2)")

    _bio_threshold_rows = []
    for _logical in ["POU5F1", "SOX2", "NANOG"]:
        _bio_threshold_rows.append({
            "marker": _logical,
            "resolved_gene": _bio_resolved_genes[_logical],
            "high_threshold": _bio_thresholds[f"{_logical}_high"],
            "very_high_threshold": _bio_thresholds[f"{_logical}_very_high"],
        })
    _bio_threshold_df = pd.DataFrame(_bio_threshold_rows)
    _bio_threshold_path = CR2_OUT_DIR / (TIMESTAMP + f"_{_MODE_TAG}_cr2_bio_ips_thresholds.tsv")
    _bio_threshold_df.to_csv(str(_bio_threshold_path), sep="	", index=False)

    _overlap_df = pd.DataFrame([
        {"metric": "n_terminal_cells", "value": int(_terminal_cell_mask.sum())},
        {"metric": "n_dynamic_positive", "value": _n_dyn},
        {"metric": "n_biology_positive", "value": _n_bio},
        {"metric": "n_overlap", "value": _n_overlap},
        {"metric": "n_dynamic_only", "value": int(_dyn_only_mask.sum())},
        {"metric": "n_biology_only", "value": int(_bio_only_mask.sum())},
        {"metric": "n_terminal_neither", "value": int(_neither_mask.sum())},
        {"metric": "jaccard", "value": _jaccard},
        {"metric": "overlap_over_dynamic", "value": _dyn_in_bio},
        {"metric": "overlap_over_biology", "value": _bio_in_dyn},
    ])
    _overlap_path = CR2_OUT_DIR / (TIMESTAMP + f"_{_MODE_TAG}_cr2_dynamic_vs_bio_overlap.tsv")
    _overlap_df.to_csv(str(_overlap_path), sep="	", index=False)

    _term_summary_df.to_csv(str(_term_summary_path), sep="	", index=False)
    print(f"\n  Terminal-state summary: {_term_summary_path.name}")
    print(f"  Biological thresholds : {_bio_threshold_path.name}")
    print(f"  Dynamic-vs-bio overlap: {_overlap_path.name}")
    # =========================================================================
    # STEP 11  —  GPCCA: compute fate probabilities
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 11  —  GPCCA: compute fate probabilities")
    print(_SEP)
    print(f"  n_jobs = {FATE_PROB_N_JOBS}  "
          f"({'serial' if FATE_PROB_N_JOBS == 1 else 'parallel (all cores)'})")
    _t0 = _time.time()
    gpcca.compute_fate_probabilities(
        n_jobs=FATE_PROB_N_JOBS,
        show_progress_bar=True,
    )
    print(f"  Done in {_time.time() - _t0:.1f} s")
    _fp_names = list(gpcca.fate_probabilities.names)
    print(f"  Fate probability lineages ({len(_fp_names)}): {_fp_names}")
    if not _selected_lineages:
        raise RuntimeError("No iPSC terminal lineage selected — cannot build cr2_fate_ips.")
    _all_fp_arrays: dict[str, np.ndarray] = {}
    for _lin in _fp_names:
        try:
            _fp_arr = np.asarray(gpcca.fate_probabilities[:, _lin]).flatten()
            _fp_arr = np.where(np.isfinite(_fp_arr), _fp_arr, 0.0).astype(np.float64)
        except Exception as _e_fp:
            print(f"  [WARN] Cannot extract fate_prob['{_lin}']: {_e_fp}")
            _fp_arr = np.zeros(adata_work.n_obs, dtype=np.float64)
        _all_fp_arrays[_lin] = _fp_arr
        _safe = _lin.replace(" ", "_").replace(",", "").replace("/", "_")
        adata_work.obs[f"cr2_fp_{_safe}"] = _fp_arr
    _missing_ips_lineages = [l for l in _selected_lineages if l not in _all_fp_arrays]
    if _missing_ips_lineages:
        raise RuntimeError(
            f"Selected iPSC terminal lineage(s) absent from fate_probabilities: {_missing_ips_lineages}"
        )
    _fp_ips_sum = np.zeros(adata_work.n_obs, dtype=np.float64)
    for _lin in _selected_lineages:
        _fp_ips_sum += _all_fp_arrays[_lin]
    _fp_ips_sum = np.clip(_fp_ips_sum, 0.0, 1.0)
    adata_work.obs["cr2_fate_ips"] = _fp_ips_sum
    _fin = np.isfinite(_fp_ips_sum)
    print(f"\n  cr2_fate_ips  "
          f"min={_fp_ips_sum[_fin].min():.4f}  "
          f"median={np.median(_fp_ips_sum[_fin]):.4f}  "
          f"max={_fp_ips_sum[_fin].max():.4f}")
    print(f"  (sum over {len(_selected_lineages)} confirmed iPSC terminal lineage(s))")
    if "marker_ips_subset_mask" in adata_work.obs.columns:
        try:
            from scipy.stats import pointbiserialr
            _marker_vec  = adata_work.obs["marker_ips_subset_mask"].astype(int).values
            _fin_mask    = np.isfinite(_fp_ips_sum)
            _r_pb, _p_pb = pointbiserialr(
                _marker_vec[_fin_mask], _fp_ips_sum[_fin_mask]
            )
            print(f"\n  Marker-subset sensitivity:")
            print(f"    Point-biserial r(cr2_fate_ips, marker_ips_subset_mask) = "
                  f"{_r_pb:.4f}  (p={_p_pb:.2e})")
        except Exception as _e:
            print(f"  [WARN] Marker-subset correlation failed: {_e}")
    try:
        _per_lin_cols = [
            f"cr2_fp_{_lin.replace(' ', '_').replace(',', '').replace('/', '_')}"
            for _lin in _selected_lineages
            if f"cr2_fp_{_lin.replace(' ', '_').replace(',', '').replace('/', '_')}"
               in adata_work.obs.columns
        ]
        _ncols      = 2 + len(_per_lin_cols)
        _fig_fp, _axes_fp = plt.subplots(1, _ncols, figsize=(7 * _ncols, 6))
        _axf = list(_axes_fp) if _ncols > 1 else [_axes_fp]
        sc.pl.umap(adata_work, color="cr2_fate_ips", ax=_axf[0],
                   show=False, size=UMAP_SIZE, alpha=UMAP_ALPHA,
                   cmap="viridis",
                   title=f"CR2 cr2_fate_ips\n({len(_selected_lineages)} iPSC terminal lineage(s))")
        sc.pl.umap(adata_work, color="p_iPSC", ax=_axf[1],
                   show=False, size=UMAP_SIZE, alpha=UMAP_ALPHA,
                   cmap="viridis", title="WOT p_iPSC")
        for _li, (_lin, _col) in enumerate(
            zip(_selected_lineages, _per_lin_cols)
        ):
            sc.pl.umap(adata_work, color=_col, ax=_axf[2 + _li],
                       show=False, size=UMAP_SIZE, alpha=UMAP_ALPHA,
                       cmap="viridis", title=f"CR2 fate: {_lin}")
        plt.suptitle(f"Fate probabilities: CellRank2  vs  WOT  [{_MODE_TAG}]",
                     fontsize=13)
        plt.tight_layout()
        savefig(_fig_fp, f"{TIMESTAMP}_{_MODE_TAG}_cr2_vs_wot_fate_umap", FIGURES_DIR)
    except Exception as _e:
        print(f"  [WARN] Fate UMAP skipped: {_e}")
    try:
        gpcca.plot_fate_probabilities(
            same_plot=False, show=False, basis="umap", size=UMAP_SIZE,
            save=f"{TIMESTAMP}_{_MODE_TAG}_cr2_fate_probs_all.{SAVE_FMT}",
        )
        print(f"  [FIG] {TIMESTAMP}_{_MODE_TAG}_cr2_fate_probs_all.{SAVE_FMT}")
    except Exception as _e:
        print(f"  [WARN] plot_fate_probabilities skipped: {_e}")
    # =========================================================================
    # STEP 12  —  Compare CR2 fate prob vs WOT p_iPSC
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 12  —  Compare CR2 fate prob  vs  WOT p_iPSC")
    print(_SEP)
    _pr = _sr = np.nan
    try:
        from scipy.stats import pearsonr, spearmanr
        _cr2_p   = adata_work.obs["cr2_fate_ips"].values.astype(float)
        _wot_p   = adata_work.obs["p_iPSC"].values.astype(float)
        # Compare on ALL cells (including iPSC-stage) — official terminal-state workflow cr2_fate_ips is
        # no longer trivially 1.0 for terminal cells, so full-dataset correlation
        # is meaningful.  We also report the non-iPSC-stage subset separately.
        _all_fin = np.isfinite(_cr2_p) & np.isfinite(_wot_p)
        _n_all   = int(_all_fin.sum())
        _non_ips = (
            (adata_work.obs["stage"] != _ips_stage_label)
            if _ips_stage_label is not None
            else pd.Series(True, index=adata_work.obs.index)
        )
        _mask    = (_non_ips & _all_fin)
        _n_cmp   = int(_mask.sum())
        print(f"  All finite pairs: {_n_all:,}  |  Non-iPSC-stage: {_n_cmp:,}")
        _pr, _pp = pearsonr(_cr2_p[_all_fin], _wot_p[_all_fin])
        _sr, _sp = spearmanr(_cr2_p[_all_fin], _wot_p[_all_fin])
        print(f"  All cells   — Pearson r={_pr:.4f}  Spearman r={_sr:.4f}")
        if _n_cmp >= 10:
            _pr2, _pp2 = pearsonr(_cr2_p[_mask], _wot_p[_mask])
            _sr2, _sp2 = spearmanr(_cr2_p[_mask], _wot_p[_mask])
            print(f"  Non-iPSC    — Pearson r={_pr2:.4f}  Spearman r={_sr2:.4f}")
        else:
            _pr2, _pp2, _sr2, _sp2 = np.nan, np.nan, np.nan, np.nan
        _corr_df = pd.DataFrame({
            "metric": [
                "pearson_r_all", "pearson_p_all",
                "spearman_r_all", "spearman_p_all",
                "pearson_r_nonips", "pearson_p_nonips",
                "spearman_r_nonips", "spearman_p_nonips",
                "n_cells_all", "n_cells_nonips",
                "n_ips_terminal_lineages",
            ],
            "value": [
                _pr, _pp, _sr, _sp,
                _pr2, _pp2, _sr2, _sp2,
                _n_all, _n_cmp,
                len(_selected_lineages),
            ],
        })
        _corr_path = (
            CR2_OUT_DIR / (TIMESTAMP + f"_{_MODE_TAG}_cr2_vs_wot_correlation.tsv")
        )
        _corr_df.to_csv(str(_corr_path), sep="\t", index=False)
        print(f"  Saved: {_corr_path.name}")
        _fig_sc, _ax_sc = plt.subplots(figsize=(6, 5))
        _stages   = adata_work.obs["stage"].values
        _stg_uniq = sorted(set(_stages[_all_fin]))
        _cmap_s   = plt.get_cmap("tab20", max(len(_stg_uniq), 1))
        for _si, _stg in enumerate(_stg_uniq):
            _m2 = _all_fin & (_stages == _stg)
            _ax_sc.scatter(_wot_p[_m2], _cr2_p[_m2],
                           s=1.5, alpha=0.3, color=_cmap_s(_si),
                           label=_stg, rasterized=True)
        _ax_sc.set_xlabel("WOT p_iPSC", fontsize=11)
        _ax_sc.set_ylabel("CR2 cr2_fate_ips (official terminal workflow)", fontsize=11)
        _ax_sc.set_title(
            f"CR2 vs WOT [{_MODE_TAG}]\n"
            f"Pearson r={_pr:.3f}  Spearman r={_sr:.3f}  "
            f"({len(_selected_lineages)} iPSC terminal lineage(s))",
            fontsize=9,
        )
        _ax_sc.legend(markerscale=5, fontsize=7, loc="lower right",
                      ncol=max(1, len(_stg_uniq) // 5))
        plt.tight_layout()
        savefig(_fig_sc, f"{TIMESTAMP}_{_MODE_TAG}_cr2_vs_wot_scatter", FIGURES_DIR)
    except Exception as _e:
        print(f"  [WARN] Comparison failed: {_e}")
    # =========================================================================
    # STEP 13  —  Lineage driver genes
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 13  —  Lineage driver genes")
    print(_SEP)
    _drivers  = None
    _corr_col = "corr"
    _driver_cluster_subset: list[str] | None = None
    if LEIDEN_KEY in adata_work.obs.columns and "cr2_fate_ips" in adata_work.obs.columns:
        try:
            _cr2_tmp = adata_work.obs["cr2_fate_ips"].values.astype(float)
            _fin_tmp = np.isfinite(_cr2_tmp)
            if _fin_tmp.any():
                _q = float(np.nanquantile(_cr2_tmp[_fin_tmp], DRIVER_TOP_QUANTILE))
                _hi = _fin_tmp & (_cr2_tmp >= _q)
                if _hi.any():
                    _cluster_freq = adata_work.obs.loc[_hi, LEIDEN_KEY].astype(str).value_counts(normalize=True)
                    _driver_cluster_subset = _cluster_freq[
                        _cluster_freq >= DRIVER_MIN_STAGE_SHARE
                    ].index.tolist()
            if _driver_cluster_subset:
                print(f"  Driver-cluster restriction ({LEIDEN_KEY} labels): {_driver_cluster_subset}")
            else:
                print("  Driver-cluster restriction unavailable — using all cells.")
        except Exception as _e:
            print(f"  [WARN] Could not derive driver-cluster subset: {_e}")
            _driver_cluster_subset = None
    def _run_lineage_drivers(lineage_name: str) -> pd.DataFrame | None:
        try:
            _kwargs = {
                "lineages": lineage_name,
                "method": "fisher",
                "use_raw": False,
                "return_drivers": True,
            }
            if _driver_cluster_subset:
                _kwargs.update({"cluster_key": LEIDEN_KEY, "clusters": _driver_cluster_subset})
            _df = gpcca.compute_lineage_drivers(**_kwargs)
            if not isinstance(_df, pd.DataFrame):
                raise TypeError(f"Expected DataFrame, got {type(_df)}")
            _cc = lineage_name + "_corr"
            if _cc not in _df.columns:
                _cc_candidates = [c for c in _df.columns if c.endswith("_corr")]
                if not _cc_candidates:
                    raise KeyError(f"No *_corr column in {list(_df.columns)}")
                _cc = _cc_candidates[0]
            _df = _df.sort_values(_cc, ascending=False).head(N_DRIVER_GENES)
            _df["_corr_col_used"] = _cc
            return _df
        except Exception as _e:
            print(f"  [WARN] compute_lineage_drivers('{lineage_name}') failed: {_e}")
            return None
    if len(_selected_lineages) == 1:
        _lin_primary = _selected_lineages[0]
        print(f"  compute_lineage_drivers(lineages='{_lin_primary}') …")
        _drivers = _run_lineage_drivers(_lin_primary)
        if _drivers is not None:
            _corr_col = _drivers["_corr_col_used"].iloc[0]
            _drivers  = _drivers.drop(columns=["_corr_col_used"])
            print(f"  Drivers for '{_lin_primary}': {len(_drivers)} genes  "
                  f"(corr col='{_corr_col}')")
    elif len(_selected_lineages) > 1:
        print(f"  Multiple selected iPSC lineages → computing consensus drivers …")
        _per_lin_drivers: dict[str, pd.DataFrame] = {}
        for _lin_d in _selected_lineages:
            print(f"    compute_lineage_drivers('{_lin_d}') …")
            _df_d = _run_lineage_drivers(_lin_d)
            if _df_d is not None:
                _cc_used = _df_d["_corr_col_used"].iloc[0]
                _df_d    = _df_d.drop(columns=["_corr_col_used"])
                _per_lin_drivers[_lin_d] = _df_d
                print(f"      {len(_df_d)} genes  (corr col='{_cc_used}')")
        if _per_lin_drivers:
            _top_sets = [set(_df.head(N_DRIVER_GENES).index)
                         for _df in _per_lin_drivers.values()]
            _consensus_genes = set.intersection(*_top_sets) if _top_sets else set()
            print(f"  Consensus (intersection of top-{N_DRIVER_GENES} across "
                  f"{len(_per_lin_drivers)} lineages): {len(_consensus_genes)} genes")
            if not _consensus_genes:
                print("  [INFO] Empty intersection — relaxing to union + mean corr.")
                _all_genes = set.union(*(set(_df.index) for _df in _per_lin_drivers.values()))
                _consensus_genes = _all_genes
            _corr_arrays = {}
            for _lin_d, _df_d in _per_lin_drivers.items():
                _cc = [c for c in _df_d.columns if c.endswith("_corr")]
                if _cc:
                    _corr_arrays[_lin_d] = _df_d[_cc[0]].rename(_lin_d)
            if _corr_arrays:
                _corr_merged = pd.concat(list(_corr_arrays.values()), axis=1)
                _corr_merged = _corr_merged.loc[
                    _corr_merged.index.isin(_consensus_genes)
                ]
                _corr_merged["mean_corr"] = _corr_merged.mean(axis=1)
                _drivers     = _corr_merged.sort_values("mean_corr", ascending=False)
                _drivers     = _drivers.head(N_DRIVER_GENES)
                _corr_col    = "mean_corr"
                print(f"  Consensus drivers: {len(_drivers)} genes")
                for _lin_d, _df_d in _per_lin_drivers.items():
                    _safe_lin = _lin_d.replace(" ", "_").replace(",", "")
                    _pld_path = CR2_OUT_DIR / (
                        TIMESTAMP + f"_{_MODE_TAG}_cr2_drivers_{_safe_lin}.tsv"
                    )
                    _df_d.to_csv(str(_pld_path), sep="	")
                    print(f"    Saved per-lineage: {_pld_path.name}")
    if _drivers is None:
        print("  Falling back to manual Pearson correlation with final biological iPSC labels if available …")
        if "ips_final_mask" in adata_work.obs.columns and adata_work.obs["ips_final_mask"].nunique() > 1:
            _target_vec = adata_work.obs["ips_final_mask"].astype(int).values.astype(np.float64)
            _corr_col = "manual_corr_bio_final"
            print("  Using `ips_final_mask` (biology-first final labels) as manual-correlation target.")
        elif "cr2_fate_ips" in adata_work.obs.columns:
            _target_vec = adata_work.obs["cr2_fate_ips"].values.astype(np.float64)
            _corr_col = "manual_corr_cr2_fate"
            print("  [WARN] `ips_final_mask` unavailable or constant — using `cr2_fate_ips` instead.")
        else:
            _target_vec = None
            print("  [WARN] Neither `ips_final_mask` nor `cr2_fate_ips` available. Skipping.")

        if _target_vec is not None:
            _ok       = np.isfinite(_target_vec)
            _target_f = _target_vec[_ok]
            _Xsub     = adata_work.X[_ok]
            if sp.issparse(_Xsub):
                _Xsub = _Xsub.toarray()
            _Xsub  = np.asarray(_Xsub, dtype=np.float64)
            _tz    = (_target_f - _target_f.mean()) / (_target_f.std() + 1e-12)
            _n_g   = _Xsub.shape[1]
            _corrs = np.full(_n_g, np.nan, dtype=np.float64)
            for _j in range(_n_g):
                _gv = _Xsub[:, _j]
                _gs = _gv.std()
                if _gs < 1e-12:
                    continue
                _corrs[_j] = float(
                    np.dot(_tz, (_gv - _gv.mean()) / _gs)
                ) / len(_target_f)
            _gene_names = np.asarray(adata_work.var_names, dtype=str)
            _drivers    = (
                pd.DataFrame({_corr_col: _corrs}, index=_gene_names)
                .dropna()
                .sort_values(_corr_col, ascending=False)
                .head(N_DRIVER_GENES)
            )
            print(f"  Manual correlation done: {len(_drivers)} genes")
    if _drivers is not None:
        _drv_path = CR2_OUT_DIR / (TIMESTAMP + f"_{_MODE_TAG}_cr2_driver_genes.tsv")
        _drivers.to_csv(str(_drv_path), sep="	")
        print(f"  Saved: {_drv_path.name}")
        print(f"\n  Top 20 driver genes (by {_corr_col}):")
        print(_drivers.head(20).to_string())
        _wot_files = sorted(
            list(METRICS_DIR.glob("*driver_genes_overlap.tsv")) +
            list((PROJECT_ROOT / "results" / "tables").glob("*driver*genes*.tsv")),
            key=lambda p: p.stat().st_mtime,
        )
        if _wot_files:
            _wot_drv = pd.read_csv(str(_wot_files[-1]), sep="	")
            _wot_col = "gene" if "gene" in _wot_drv.columns else _wot_drv.columns[0]
            _wot_set = set(_wot_drv[_wot_col].astype(str))
            _cr2_set = set(_drivers.index.astype(str))
            _overlap = _wot_set & _cr2_set
            print(f"\n  WOT drivers  : {len(_wot_set)}")
            print(f"  CR2 drivers  : {len(_cr2_set)}")
            print(f"  Overlap (∩)  : {len(_overlap)}")
            if _overlap:
                print(f"  Genes        : {sorted(_overlap)[:30]}")
                _ovlp_df   = _drivers.loc[[g for g in _drivers.index if g in _overlap]]
                _ovlp_path = CR2_OUT_DIR / (
                    TIMESTAMP + f"_{_MODE_TAG}_cr2_wot_overlap_drivers.tsv"
                )
                _ovlp_df.to_csv(str(_ovlp_path), sep="	")
                print(f"  Saved: {_ovlp_path.name}")
        else:
            print("  [INFO] WOT driver table not found — skipping overlap.")
        try:
            _top10    = list(_drivers.head(10).index)
            _in_adata = [g for g in _top10 if g in adata_work.var_names]
            _stage_order = sorted(
                adata_work.obs["stage"].unique(),
                key=lambda s: float(re.search(r"\d+", s).group(0))
                    if re.search(r"\d+", s) else 999.0,
            )
            if _in_adata:
                _nc      = len(_in_adata)
                _nc_cols = max(1, (_nc + 1) // 2)
                _fig_gt, _axes_gt = plt.subplots(
                    2, _nc_cols, figsize=(4 * _nc_cols, 8))
                _axes_flat = np.asarray(_axes_gt).flatten()
                for _gi, _gene in enumerate(_in_adata):
                    _ax     = _axes_flat[_gi]
                    _gi_idx = list(adata_work.var_names).index(_gene)
                    _means, _stds = [], []
                    for _stg in _stage_order:
                        _idx = (adata_work.obs["stage"] == _stg).values
                        if not _idx.any():
                            _means.append(np.nan); _stds.append(0.0); continue
                        _expr = (
                            np.asarray(
                                adata_work.X[_idx, _gi_idx].todense()
                            ).flatten()
                            if sp.issparse(adata_work.X)
                            else np.asarray(
                                adata_work.X[_idx, _gi_idx]
                            ).flatten()
                        )
                        _means.append(float(np.nanmean(_expr)))
                        _stds.append(float(np.nanstd(_expr)))
                    _means = np.asarray(_means, dtype=float)
                    _stds  = np.asarray(_stds,  dtype=float)
                    _ax.plot(_stage_order, _means, "o-", lw=2, ms=5)
                    _ax.fill_between(
                        range(len(_stage_order)),
                        np.where(np.isfinite(_means - _stds), _means - _stds, _means),
                        np.where(np.isfinite(_means + _stds), _means + _stds, _means),
                        alpha=0.2,
                    )
                    _ax.set_xticks(range(len(_stage_order)))
                    _ax.set_xticklabels(
                        _stage_order, rotation=45, ha="right", fontsize=7)
                    _ax.set_title(_gene, fontsize=9)
                    _ax.set_ylabel("mean log1p", fontsize=7)
                for _gi2 in range(len(_in_adata), len(_axes_flat)):
                    _axes_flat[_gi2].set_visible(False)
                plt.suptitle(
                    f"CR2 driver gene expression trends [{_MODE_TAG}]", fontsize=11)
                plt.tight_layout()
                savefig(
                    _fig_gt,
                    f"{TIMESTAMP}_{_MODE_TAG}_cr2_driver_gene_trends",
                    FIGURES_DIR,
                )
        except Exception as _e:
            print(f"  [WARN] Gene trend plots failed: {_e}")
    # =========================================================================
    # STEP 14  —  Save outputs
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 14  —  Save outputs")
    print(_SEP)
    try:
        if _ms_obs_key not in adata_work.obs.columns:
            adata_work.obs[_ms_obs_key] = gpcca.macrostates
            print(f"  Added obs['{_ms_obs_key}'] from GPCCA.")
    except Exception as _e:
        print(f"  [WARN] Could not embed macrostates: {_e}")
    _out_h5ad = METRICS_DIR / (TIMESTAMP + f"_{_MODE_TAG}_GSE230659_cellrank2.h5ad")
    try:
        adata_work.write_h5ad(str(_out_h5ad))
        print(f"  Saved h5ad : {_out_h5ad.name}")
    except Exception as _e:
        print(f"  [WARN] h5ad save failed: {_e}")
    _cr2_fate_ips_vals = adata_work.obs.get("cr2_fate_ips", pd.Series(dtype=float)).values
    _fin_cr2 = np.isfinite(_cr2_fate_ips_vals.astype(float))
    _cr2_median = (
        float(np.median(_cr2_fate_ips_vals[_fin_cr2])) if _fin_cr2.any() else np.nan
    )
    _rows = [
        ("script",                 "04_cellrank2_v12_gpt.py (dynamic-vs-biological iPSC endpoint comparison)"),
        ("run_mode",               RUN_MODE),
        ("mode_tag",               _MODE_TAG),
        ("timestamp",              TIMESTAMP),
        ("input_h5ad",             fates_h5ad.name),
        ("tmaps_dir_original",     str(tmaps_dir)),
        ("tmaps_dir_work",         str(tmaps_dir_work)),
        ("n_tmap_files_work",      len(_tmap_files_work)),
        ("n_cells_full",           adata.n_obs),
        ("n_cells_work",           adata_work.n_obs),
        ("n_genes",                adata_work.n_vars),
        ("n_timepoints",           len(_day_cats_work)),
        ("ips_stage_label",        str(_ips_stage_label)),
        ("conn_weight",            CONN_WEIGHT),
        ("threshold",              THRESHOLD),
        ("n_macrostates",          N_MACROSTATES),
        ("cluster_key",            LEIDEN_KEY),
        ("leiden_resolution",      LEIDEN_RESOLUTION),
        ("all_macrostate_names",   str(_ms_names)),
        ("n_terminal_lineages",    len(_fp_names)),
        ("terminal_method",        TERMINAL_METHOD),
        ("terminal_stability_threshold", TERMINAL_STABILITY_THRESHOLD),
        ("min_macrostate_cells",    MIN_MACROSTATE_CELLS),
        ("ips_terminal_max_states", IPS_TERMINAL_MAX_STATES),
        ("n_auto_terminal_states", len(_auto_terminal_states)),
        ("auto_terminal_states",   str(_auto_terminal_states)),
        ("n_final_terminal_states", len(_ts_names)),
        ("final_terminal_states",  str(_ts_names)),
        ("n_ips_terminal_lineages", len(_selected_lineages)),
        ("ips_terminal_lineages",  str(_selected_lineages)),
        ("n_terminal_cells",       int(_terminal_cell_mask.sum())),
        ("n_dynamic_ips_cells",    int(_n_dyn)),
        ("n_bio_ips_cells",        int(_n_bio)),
        ("n_dyn_bio_overlap",      int(_n_overlap)),
        ("dyn_bio_jaccard",        f"{_jaccard:.4f}"),
        ("overlap_over_dynamic",   f"{_dyn_in_bio:.4f}"),
        ("overlap_over_biology",   f"{_bio_in_dyn:.4f}"),
        ("final_detection_method", "biology_marker_classes"),
        ("bio_positive_categories", "['POU5F1-only','single_high','bi_high','tri_high','tri_very_high']"),
        ("cr2_fate_ips_median",    f"{_cr2_median:.4f}"),
        ("fate_prob_n_jobs",       str(FATE_PROB_N_JOBS)),
        ("n_driver_genes",         len(_drivers) if _drivers is not None else "N/A"),
        ("driver_corr_col",        _corr_col),
        ("pearson_r_all",          f"{_pr:.4f}"),
        ("spearman_r_all",         f"{_sr:.4f}"),
        ("terminal_state_summary_tsv", _term_summary_path.name),
        ("bio_threshold_tsv",      _bio_threshold_path.name),
        ("dynamic_vs_bio_overlap_tsv", _overlap_path.name),
        ("petsc_ok",               str(_petsc_ok)),
        ("slepc_ok",               str(_slepc_ok)),
    ]
    _sum_df   = pd.DataFrame(_rows, columns=["parameter", "value"])
    _sum_path = CR2_OUT_DIR / (TIMESTAMP + f"_{_MODE_TAG}_cr2_summary.tsv")
    _sum_df.to_csv(str(_sum_path), sep="\t", index=False)
    print(f"  Saved summary : {_sum_path.name}")
    print("\n" + _sum_df.to_string(index=False))
    # =========================================================================
    # DONE
    # =========================================================================
    print(f"\n{_SEP}")
    _elapsed = _time.time() - T0
    print(f"04_cellrank2_v12_gpt.py  [{_MODE_TAG}]  COMPLETE  ({_elapsed:.1f} s total)")
    print(_SEP)
    print(f"\n  Outputs in : {CR2_OUT_DIR}")
    print(f"  Figures in : {FIGURES_DIR}")
    print(f"\n  Key outputs:")
    print(f"    {_out_h5ad.name}")
    print(f"    {TIMESTAMP}_{_MODE_TAG}_cr2_driver_genes.tsv")
    print(f"    {TIMESTAMP}_{_MODE_TAG}_cr2_terminal_state_summary.tsv")
    print(f"    {TIMESTAMP}_{_MODE_TAG}_cr2_summary.tsv")
    print(f"\n  Leiden-guided terminal-state workflow: {len(_selected_lineages)} iPSC terminal lineage(s) selected")
    print(f"    → {_selected_lineages}")
    if RUN_MODE == "local_debug":
        print(
            "\n  [LOCAL DEBUG SMOKE TEST COMPLETE]\n"
            f"  Worked on {adata_work.n_obs:,} cells "
            f"(stratified subsample of {adata.n_obs:,}).\n"
            f"  Cropped tmaps : {tmaps_dir_work.name}/\n"
            f"  Fate prob n_jobs = {FATE_PROB_N_JOBS} (serial — no child processes).\n"
            "  This run is for pipeline validation only, not final biological conclusions.\n"
            "  To run the full analysis on Shirokane:\n"
            "    RUN_MODE=hpc_full python scripts/WOT/04_cellrank2.py"
        )
# =============================================================================
# ENTRY POINT
# =============================================================================
# freeze_support() is required when packaging with PyInstaller / cx_Freeze.
# It is harmless when running as a plain script and costs zero overhead.
#
# The if __name__ == "__main__": guard is CRITICAL on Windows.
# Windows multiprocessing uses the "spawn" start method, which imports the
# script as a module in each child process.  Without this guard, the entire
# pipeline would re-execute inside each worker, causing infinite recursion
# and an immediate RuntimeError.
# =============================================================================
if __name__ == "__main__":
    freeze_support()
    main()
