#!/usr/bin/env python3
"""
04_cellrank2.py  [v12 — standard GPCCA flow, OCLR dual-endpoint lineage annotation]
==========================================================================================
Project : Comparative Study of Trajectory Inference Models for Chemical iPSC
          Reprogramming
Dataset : GSE230659 (human, Liuyang et al. 2023 Cell Stem Cell)
Step    : 4 — CellRank2 trajectory analysis (second method, comparative study)
Env     : conda activate traj_env

Analysis design
---------------
This script implements Problem B: CellRank2 RealTimeKernel as a WOT extension.
  - WOT               = baseline fate-probability model (Problem A, script 02c)
  - CellRank2 RTK     = WOT extension: adds within-timepoint connectivity on top
                        of the same WOT transport maps
  - The goal is NOT to build an independent CellRank2 model; it is to test
    whether adding kNN-based connectivity to WOT improves fate mapping.
  Therefore RealTimeKernel.from_wot() is used deliberately, and the same tmaps
  produced by script 02c are the primary input.

Standard GPCCA flow (v12)
--------------------------
  STEP 10: predict_terminal_states() — standard GPCCA method to identify
           biologically terminal macrostates.  NOT all macrostates are set as
           terminal (that was Scheme B — a non-standard rescue heuristic).

  STEP 11: compute_fate_probabilities() — standard CellRank2 API.
           After computing per-lineage fate probabilities, the OCLR dual-endpoint
           labels (from 02b2) are used as EXTERNAL ANNOTATION to identify
           the SUCCESS and FAILURE lineages:
             - For each terminal lineage l:
                 mean_fp_success(l) = mean(fate_prob_l among success endpoint cells)
                 mean_fp_failure(l) = mean(fate_prob_l among failure endpoint cells)
             - cr2_success_lineage = lineage with highest mean_fp_success
             - cr2_failure_lineage = lineage with highest mean_fp_failure
               (if same as success, use second-best)
           Also retained: Spearman(fate_prob_l, oclr_score) for diagnostics.

  Key output columns:
    cr2_p_oclr_success   = P(cell → success lineage)  [PRIMARY metric]
    cr2_p_oclr_failure   = P(cell → failure lineage)
    cr2_oclr_margin      = cr2_p_oclr_success − cr2_p_oclr_failure  [PRIMARY]
    cr2_ips_fate_prob    = alias for cr2_p_oclr_success  (legacy)
    cr2_ips_lineage_name = alias for cr2_success_lineage_name  (legacy)

Changelog v12
-------------------------------------------------------------
Extend OCLR annotation to dual-endpoint: mean(fate_prob in success/failure cells).
Add: cr2_p_oclr_success, cr2_p_oclr_failure, cr2_oclr_margin,
     cr2_success_lineage_name, cr2_failure_lineage_name.
Load OCLR endpoint labels from *_oclr_endpoint_labels_hcipsc.tsv (02b2 v14).
Keep cr2_ips_fate_prob = cr2_p_oclr_success (legacy alias).

Changelog v11
-------------------------------------------------------------
Remove Scheme B entirely.  Restore standard GPCCA flow:
  predict_terminal_states() + compute_fate_probabilities().
Use OCLR scores as post-hoc annotation only to identify cr2_ips_lineage.
Remove: MAX_SELECTED_LINEAGES, DELTA_FP_MIN, delta_fp, custom cascade,
        cr2_fate_ips pooling, fate_ips/fate_other column references.
Add: cr2_ips_fate_prob, cr2_ips_lineage_name (OCLR-annotated single lineage).

Changelog v10 (superseded — Scheme B removed in v11)
-------------------------------------------------------------
Fix F — lineage selection: delta_fp within day30 terminal pool.
Fix G — eigenspectrum: removed show= kwarg.

Changelog v9
-------------------------------------------------------------
Fix A — terminal-states obs key ('term_states_fwd', not 'terminal_states_fwd')
Fix B — eigenspectrum: real_only=True, no ax=, capture via plt.gcf()
Fix C — lineage selection: macrostate-fraction-based 4-tier cascade
Fix D — driver-gene trend: .values on boolean Series before sparse indexing
Fix E — suppress 'more than two extensions' CellRank INFO messages

Changelog v8
-------------------------------------------------------------
Scheme B — per-lineage fate_ips enrichment + selective pooling

  ROOT CAUSE OF v6/v7 DEGENERACY
  --------------------------------
  In v6/v7, STEP 10 set ONLY the iPSC-matched macrostates as terminal states.
  When GPCCA finds multiple macrostates all labeled "hCiPSCs_N", summing ALL
  their fate probabilities yields ≈1.0 for every cell (since the terminal pool
  covers essentially all fate probability mass).  This gives a completely
  uninformative cr2_fate_terminal_pool column.

  SCHEME B FIX
  ------------
  1. STEP 10: Set ALL macrostates as terminal states.
     Rationale: all late-time macrostates (iPSC, partial, fail) are biologically
     terminal in the observed timeframe.  With multiple distinct terminal states,
     fate probabilities vary meaningfully across cells.

  2. STEP 11: Per-lineage iPSC enrichment analysis.
     For each terminal lineage l, compute:
       weighted_frac_fate_ips(l) = Σ(fate_prob_l · fate_ips) / Σ(fate_prob_l)
     where fate_ips = 1 for day30 hCiPSC AND OCLR-defined terminal (02b2/02c definition).
     This detects which macrostates correspond to the high-confidence, fully
     reprogrammed iPSC vs partially reprogrammed / failed cells.

  3. Selection (4-tier cascade): primary = frac_fate_ips_macrostate >= MS_IPS_FRAC_THRESHOLD
     AND n_assigned >= MIN_MACROSTATE_CELLS.  Fallback tiers avoid empty selection.
     Pool selected lineages → cr2_fate_ips (replaces cr2_fate_terminal_pool).
     Save per-lineage enrichment summary table (_cr2_lineage_enrichment.tsv).

  4. Comparison and driver genes use cr2_fate_ips exclusively.
     Multiple selected lineages → consensus driver genes (intersection of top-N).

Changelog v6
-------------------------------------------------------------
Fix 8 — Windows multiprocessing guard (main() + freeze_support)

  ROOT CAUSE OF v5 FAILURE
  -------------------------
  On Windows, Python's multiprocessing module uses the "spawn" start method
  (unlike Linux/macOS which default to "fork").  When CellRank calls
  compute_fate_probabilities(), joblib spawns worker processes to solve the
  linear systems in parallel.  Each worker process imports the script as a
  Python module.  Because v5 had all pipeline logic at module level (top level
  of the script), the entire pipeline re-executed inside each child process.
  The child process then tried to load data, locate tmaps, subsample, etc. and
  eventually crashed with:
    RuntimeError: An attempt has been made to start a new process before the
    current process has finished its bootstrapping phase.
    ...consider using `if __name__ == '__main__':`

  FIX
  ---
  1. All pipeline logic (STEP 1 through save outputs) is moved into main().
  2. Helper functions and config constants remain at module level (safe to import).
  3. The entry point is:
       if __name__ == "__main__":
           freeze_support()   # required for frozen executables; harmless otherwise
           main()
  4. In local_debug mode, compute_fate_probabilities() is called with n_jobs=1,
     which forces joblib into serial execution — no child processes are spawned
     at all.  This is verified against the CellRank source:
       cellrank/estimators/mixins/_fate_probabilities.py  L178-189
       cellrank/_utils/_parallelize.py  L128-150 (_get_n_cores returns 1 if n_jobs==1)
  5. In hpc_full mode, n_jobs=None (joblib default) is used to utilise all
     available cores on the Shirokane server.

Changelog v5
-------------------------------------------------------------
Fix 7 — local_debug tmap cropping

  RealTimeKernel.from_wot() does NOT silently subset coupling matrices.
  In v4 the full-dataset tmaps were passed after subsampling, causing:
    IndexError: Source observations for (...) don't match with adata.obs_names
  Fix: build_cropped_local_tmaps() crops each tmap to exactly the sampled
  barcodes before calling from_wot().

Changelog v4
-------------------------------------------------------------
New: dual execution mode (RUN_MODE = "local_debug" | "hpc_full")

Changelog v3
-------------------------------------------------------------
Fix 6 — dynamic iPSC-stage label detection + three-tier terminal-state selection.

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
    Path(r"C:\Users\37620\trajectory"),
    Path("/sessions/cool-admiring-hawking/mnt/trajectory"),
    Path(__file__).resolve().parents[2],
]
PROJECT_ROOT  = next((p for p in _candidates if p.exists()), _candidates[-1])
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

# ── Lineage driver genes ──────────────────────────────────────────────────────
N_DRIVER_GENES: int = 300

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
    if RUN_MODE == "local_debug":
        print(f"  ┌─────────────────────────────────────────────────────────────┐")
        print(f"  │  LOCAL DEBUG — stratified subsample + n_jobs=1             │")
        print(f"  │  Results are for pipeline smoke-testing only.              │")
        print(f"  │  For final model comparison, run with RUN_MODE=hpc_full   │")
        print(f"  │  on Shirokane and use the resulting hpcfull h5ad.          │")
        print(f"  └─────────────────────────────────────────────────────────────┘")
    else:
        print(f"  *** hpc_full mode: full dataset, parallel fate probabilities ***")
        print(f"  *** Output from this run is suitable for final conclusions.  ***")
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
    _req  = ["day", "stage", "p_hcipsc"]
    _miss = [c for c in _req if c not in adata.obs.columns]
    if _miss:
        raise ValueError(
            f"Required obs columns missing: {_miss}\n"
            f"Ensure script 02c (v13+, multi-fate design) has run and produced "
            f"p_hcipsc in the wot_fates.h5ad."
        )

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
    # STEP 6  —  Ensure kNN neighbor graph
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 6  —  Ensure kNN neighbor graph")
    print(_SEP)

    if "connectivities" in adata_work.obsp:
        print("  kNN graph already present in adata_work.obsp.")
    else:
        if "X_pca" not in adata_work.obsm:
            print("  X_pca not found — computing PCA (50 components) …")
            sc.pp.pca(adata_work, n_comps=50)
        else:
            print(f"  X_pca found: shape {adata_work.obsm['X_pca'].shape}")
        print(f"  sc.pp.neighbors(n_neighbors={N_NEIGHBORS}) …")
        sc.pp.neighbors(adata_work, n_neighbors=N_NEIGHBORS, use_rep="X_pca")
        print("  Done.")

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
        # Fix G: show= is not accepted by all CellRank builds; backend is Agg,
        # so no window appears regardless.  Capture via plt.gcf() after the call.
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

    print(f"  n_states={N_MACROSTATES},  cluster_key='stage'")
    _t0 = _time.time()
    try:
        gpcca.compute_macrostates(
            n_states=N_MACROSTATES,
            cluster_key="stage",
            n_cells=30,
        )
        print(f"  Done in {_time.time() - _t0:.1f} s")
    except Exception as _e:
        print(f"  [WARN] n_states={N_MACROSTATES} failed: {_e}. Trying 4 …")
        gpcca.compute_macrostates(n_states=4, cluster_key="stage")
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
                   title="CellRank2 macrostates")
        savefig(_fig_ms, f"{TIMESTAMP}_{_MODE_TAG}_cr2_macrostates_umap", FIGURES_DIR)
    except Exception as _e:
        print(f"  [WARN] Macrostates UMAP skipped: {_e}")

    # =========================================================================
    # STEP 10  —  GPCCA: predict terminal states  (standard GPCCA)
    # =========================================================================
    # Uses the standard GPCCA method predict_terminal_states() which identifies
    # terminal macrostates based on the Schur decomposition.  This is the
    # canonical CellRank2 / pyGPCCA approach and does NOT require external
    # labels or manual selection.
    #
    # NOTE: "set all macrostates as terminal" (Scheme B, v8–v10) was a
    # non-standard rescue heuristic.  It has been removed in v11.
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 10  —  GPCCA: predict terminal states  (standard GPCCA)")
    print(_SEP)

    if not _ms_names:
        print("  [FATAL] No macrostates found — cannot proceed.")
        sys.exit(1)

    try:
        gpcca.predict_terminal_states()
        print("  predict_terminal_states() completed.")
    except Exception as _e:
        print(f"  [FATAL] predict_terminal_states() failed: {_e}")
        sys.exit(1)

    _ts_names = list(gpcca.terminal_states.cat.categories)
    print(f"\n  Active terminal states ({len(_ts_names)}): {_ts_names}")

    # ── Detect the actual obs-column key written by CellRank ──────────────────
    # CellRank 2 writes 'term_states_fwd' (confirmed in reproducibility repo).
    # Older builds / future refactors may use 'terminal_states_fwd'.
    # We probe adata_work.obs directly to pick whichever key exists.
    _ts_obs_key: str | None = None
    for _k_cand in ("term_states_fwd", "terminal_states_fwd"):
        if _k_cand in adata_work.obs.columns:
            _ts_obs_key = _k_cand
            break
    if _ts_obs_key is None:
        # Last-resort: any column whose name ends with '_states_fwd'
        _fwd_cols = [c for c in adata_work.obs.columns
                     if c.endswith("_states_fwd") and "terminal" not in c.lower()
                     or c.endswith("states_fwd")]
        _ts_obs_key = _fwd_cols[0] if _fwd_cols else None
    print(f"  Terminal-states obs key detected: '{_ts_obs_key}'")

    # ── UMAP: terminal states ──────────────────────────────────────────────────
    if _ts_obs_key is not None:
        try:
            _fig_ts, _ax_ts = plt.subplots(figsize=(7, 6))
            sc.pl.umap(adata_work, color=_ts_obs_key,
                       ax=_ax_ts, show=False, size=UMAP_SIZE, alpha=UMAP_ALPHA,
                       title=f"CellRank2 terminal states [{_ts_obs_key}]")
            savefig(_fig_ts, f"{TIMESTAMP}_{_MODE_TAG}_cr2_terminal_states_umap",
                    FIGURES_DIR)
        except Exception as _e:
            print(f"  [WARN] Terminal states UMAP skipped: {_e}")
    else:
        print("  [WARN] Terminal-states obs key not found in adata_work.obs — "
              "UMAP skipped.")

    # =========================================================================
    # STEP 11  —  GPCCA: compute fate probabilities  (standard API)
    #             + OCLR-based post-hoc annotation of iPSC lineage
    # =========================================================================
    # Standard CellRank2 API: compute_fate_probabilities() using the terminal
    # states identified by predict_terminal_states() in STEP 10.
    #
    # After computing per-lineage fate probabilities, use OCLR scores
    # (from 02b2/02c) as an EXTERNAL ANNOTATION to identify the most
    # iPSC-enriched terminal lineage:
    #
    #   For each terminal lineage l:
    #     oclr_spearman_r(l) = Spearman(fate_prob_l, oclr_score)
    #                          among hCiPSC cells with finite OCLR scores
    #
    #   cr2_ips_lineage = lineage with highest oclr_spearman_r
    #   cr2_ips_fate_prob = fate probability toward cr2_ips_lineage
    #
    # OCLR is used POST-HOC only — NOT to define terminal states.
    # =========================================================================
    print(f"\n{_SEP}")
    print("STEP 11  —  GPCCA: compute fate probabilities  (standard)")
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

    # ── Extract per-lineage fate probability arrays ───────────────────────────
    _all_fp_arrays: dict[str, np.ndarray] = {}
    for _lin in _fp_names:
        try:
            _fp_arr = np.asarray(gpcca.fate_probabilities[:, _lin]).flatten()
            _fp_arr = np.where(np.isfinite(_fp_arr), _fp_arr, 0.0).astype(np.float64)
        except Exception as _e_fp:
            print(f"  [WARN] Cannot extract fate_prob['{_lin}']: {_e_fp}")
            _fp_arr = np.zeros(adata_work.n_obs, dtype=np.float64)
        _all_fp_arrays[_lin] = _fp_arr
        # Also store in adata_work for plotting / downstream use
        _safe = _lin.replace(" ", "_").replace(",", "").replace("/", "_")
        adata_work.obs[f"cr2_fp_{_safe}"] = _fp_arr

    # ── OCLR dual-endpoint post-hoc annotation ───────────────────────────────
    #
    # Load the OCLR endpoint labels produced by 02b2_define_oclr_endpoint.py.
    # For each terminal lineage l:
    #   mean_fp_success(l)       = mean(fate_prob_l among success-endpoint cells)
    #   mean_fp_failure(l)       = mean(fate_prob_l among failure-endpoint cells)
    #   delta_success_failure(l) = mean_fp_success(l) − mean_fp_failure(l)
    #
    # Lineage selection (symmetric around the separation metric):
    #   cr2_success_lineage = lineage with LARGEST  delta_success_failure
    #                         (most enriched in success over failure)
    #   cr2_failure_lineage = lineage with SMALLEST delta_success_failure
    #                         (most enriched in failure over success)
    # If success and failure collapse to the same lineage (only one lineage total),
    # failure is left unresolved: cr2_p_oclr_failure = NaN, cr2_oclr_margin = NaN.
    # Also retain: Spearman(fate_prob_l, oclr_score) for diagnostics.
    _cr2_ips_lineage_name:     str | None = None
    _cr2_success_lineage_name: str | None = None
    _cr2_failure_lineage_name: str | None = None
    _cr2_ips_oclr_spearman_r: float       = np.nan
    _oclr_annotation_rows: list[dict]     = []

    # ── Load OCLR endpoint labels ─────────────────────────────────────────────
    _oclr_label_files = sorted(PROCESSED_DIR.glob("*_oclr_endpoint_labels_hcipsc.tsv"))
    if _oclr_label_files:
        _ldf = pd.read_csv(str(_oclr_label_files[-1]), sep="\t", index_col=0)
        _ldf.index = _ldf.index.astype(str)
        _adata_bc = adata_work.obs_names.astype(str)
        _suc_set  = set(_ldf.index[_ldf.get(
            "is_oclr_success_endpoint",
            _ldf.get("is_oclr_terminal", pd.Series(0, index=_ldf.index))) == 1
        ].tolist())
        _fail_set = set(_ldf.index[
            _ldf["is_oclr_failure_endpoint"] == 1
        ].tolist()) if "is_oclr_failure_endpoint" in _ldf.columns else set()
        _suc_mask_bool  = np.array([bc in _suc_set  for bc in _adata_bc])
        _fail_mask_bool = np.array([bc in _fail_set for bc in _adata_bc])
        print(f"  OCLR endpoint labels loaded: {_ldf.shape[0]:,} cells")
        print(f"    success endpoint in adata : {int(_suc_mask_bool.sum()):,}")
        print(f"    failure endpoint in adata : {int(_fail_mask_bool.sum()):,}")
    else:
        # Fallback: use obs columns if present (loaded from wot_fates.h5ad)
        _suc_mask_bool  = (adata_work.obs.get(
            "oclr_success_endpoint_mask",
            adata_work.obs.get("oclr_ips_subset_mask", pd.Series(0, index=adata_work.obs.index))
        ) == 1).values
        _fail_mask_bool = (adata_work.obs.get(
            "oclr_failure_endpoint_mask", pd.Series(0, index=adata_work.obs.index)
        ) == 1).values
        print(f"  [FALLBACK] OCLR endpoint labels from obs columns  "
              f"(success: {int(_suc_mask_bool.sum()):,}  "
              f"failure: {int(_fail_mask_bool.sum()):,})")

    # Continuous OCLR score for Spearman diagnostics
    _oclr_scores_work = adata_work.obs.get(
        "oclr_score", pd.Series(np.nan, index=adata_work.obs.index)
    ).values.astype(float)

    # Mask for hCiPSC cells (broad — includes all OCLR-labeled cells)
    _hcipsc_mask_bool = _suc_mask_bool | _fail_mask_bool
    if not _hcipsc_mask_bool.any():
        # Fall back to stage label
        _hcipsc_mask_bool = adata_work.obs["stage"].apply(
            lambda s: "ips" in str(s).lower()
        ).values.astype(bool)

    print(f"\n  OCLR dual-endpoint annotation:")
    print(f"    hCiPSC cells (suc|fail)   : {int(_hcipsc_mask_bool.sum()):,}")
    print(f"    success endpoint cells    : {int(_suc_mask_bool.sum()):,}")
    print(f"    failure endpoint cells    : {int(_fail_mask_bool.sum()):,}")

    try:
        from scipy.stats import spearmanr as _spearmanr_fn
        for _lin in _fp_names:
            _fp_arr = _all_fp_arrays[_lin]

            # Mean fate_prob among success / failure endpoint cells
            _suc_valid  = _suc_mask_bool  & np.isfinite(_fp_arr)
            _fail_valid = _fail_mask_bool & np.isfinite(_fp_arr)
            _mean_suc   = float(_fp_arr[_suc_valid].mean())  if _suc_valid.any()  else np.nan
            _mean_fail  = float(_fp_arr[_fail_valid].mean()) if _fail_valid.any() else np.nan
            _delta_sf   = (_mean_suc - _mean_fail) if (
                np.isfinite(_mean_suc) and np.isfinite(_mean_fail)) else np.nan

            # Spearman vs continuous score (diagnostic)
            _sp_valid = _hcipsc_mask_bool & np.isfinite(_oclr_scores_work) & np.isfinite(_fp_arr)
            _n_sp     = int(_sp_valid.sum())
            if _n_sp >= 10:
                _r_sp, _p_sp = _spearmanr_fn(
                    _fp_arr[_sp_valid], _oclr_scores_work[_sp_valid])
                _r_sp, _p_sp = float(_r_sp), float(_p_sp)
            else:
                _r_sp = _p_sp = np.nan

            _oclr_annotation_rows.append({
                "lineage":             _lin,
                "n_success_valid":     int(_suc_valid.sum()),
                "n_failure_valid":     int(_fail_valid.sum()),
                "mean_fp_success":     _mean_suc,
                "mean_fp_failure":     _mean_fail,
                "delta_success_failure": _delta_sf,
                "n_spearman_valid":    _n_sp,
                "oclr_spearman_r":     _r_sp,
                "oclr_spearman_p":     _p_sp,
                "fp_mean_hcipsc":      float(_fp_arr[_hcipsc_mask_bool].mean())
                                       if _hcipsc_mask_bool.any() else np.nan,
                "fp_mean_all":         float(_fp_arr.mean()),
            })

        _oclr_ann_df = pd.DataFrame(_oclr_annotation_rows)

        print(f"\n  Per-lineage dual-endpoint annotation:")
        _disp = ["lineage", "n_success_valid", "mean_fp_success",
                 "n_failure_valid", "mean_fp_failure",
                 "delta_success_failure", "oclr_spearman_r"]
        print(_oclr_ann_df[[c for c in _disp if c in _oclr_ann_df.columns]].to_string(
            index=False))

        # ── Select success lineage: largest delta_success_failure ────────────
        # Uses the separation metric (success_enrichment − failure_enrichment),
        # which is symmetric and directly interpretable as benchmark alignment.
        _valid_delta = _oclr_ann_df["delta_success_failure"].dropna()
        if not _valid_delta.empty:
            _suc_best_idx = _valid_delta.idxmax()
            _cr2_success_lineage_name = str(
                _oclr_ann_df.loc[_suc_best_idx, "lineage"])
            print(f"\n  cr2_success_lineage_name = '{_cr2_success_lineage_name}'  "
                  f"(delta_sf={_valid_delta[_suc_best_idx]:.4f}  "
                  f"mean_fp_success="
                  f"{_oclr_ann_df.loc[_suc_best_idx, 'mean_fp_success']:.4f})")
        else:
            _name_match = [l for l in _fp_names if "ips" in l.lower()]
            _cr2_success_lineage_name = _name_match[0] if _name_match else (
                _fp_names[0] if _fp_names else None)
            print(f"\n  [WARN] delta_success_failure all NaN — name fallback: "
                  f"'{_cr2_success_lineage_name}'")

        # ── Select failure lineage: smallest delta_success_failure ────────────
        # Symmetric to success: most failure-enriched lineage.
        # If it matches success (only one lineage), failure is left unresolved.
        if not _valid_delta.empty:
            _fail_best_idx    = _valid_delta.idxmin()
            _fail_candidate   = str(_oclr_ann_df.loc[_fail_best_idx, "lineage"])
            if _fail_candidate != _cr2_success_lineage_name:
                _cr2_failure_lineage_name = _fail_candidate
                print(f"  cr2_failure_lineage_name = '{_cr2_failure_lineage_name}'  "
                      f"(delta_sf={_valid_delta[_fail_best_idx]:.4f}  "
                      f"mean_fp_failure="
                      f"{_oclr_ann_df.loc[_fail_best_idx, 'mean_fp_failure']:.4f})")
            else:
                # Collision: only one distinct lineage — dual-endpoint unresolvable.
                # Do NOT reuse the success lineage for failure; leave failure as NaN.
                _cr2_failure_lineage_name = None
                print(f"  [WARN] Dual-endpoint separation UNRESOLVED: success and failure "
                      f"both map to '{_cr2_success_lineage_name}'.\n"
                      f"         Likely cause: GPCCA identified only one terminal lineage.\n"
                      f"         cr2_p_oclr_failure and cr2_oclr_margin will be NaN.\n"
                      f"         This dataset may not support a dual-endpoint benchmark.")
        else:
            _cr2_failure_lineage_name = None
            print(f"  [WARN] delta_success_failure all NaN — failure lineage not assigned.")

        # ── Get Spearman r for the success lineage (diagnostic) ──────────────
        _suc_row = _oclr_ann_df[
            _oclr_ann_df["lineage"] == _cr2_success_lineage_name]
        if not _suc_row.empty:
            _cr2_ips_oclr_spearman_r = float(
                _suc_row["oclr_spearman_r"].values[0])

        # Legacy alias
        _cr2_ips_lineage_name = _cr2_success_lineage_name

    except Exception as _e_ann:
        print(f"  [WARN] OCLR dual-endpoint annotation failed: {_e_ann}")
        import traceback; traceback.print_exc()
        _oclr_ann_df = pd.DataFrame()
        _name_match  = [l for l in _fp_names if "ips" in l.lower()]
        _cr2_success_lineage_name = (
            _name_match[0] if _name_match else (_fp_names[0] if _fp_names else None))
        _cr2_failure_lineage_name = None
        _cr2_ips_lineage_name     = _cr2_success_lineage_name

    # ── Attach dual-endpoint columns to adata_work ───────────────────────────
    # cr2_p_oclr_success / cr2_p_oclr_failure / cr2_oclr_margin
    if _cr2_success_lineage_name and _cr2_success_lineage_name in _all_fp_arrays:
        _suc_fp = _all_fp_arrays[_cr2_success_lineage_name]
        adata_work.obs["cr2_p_oclr_success"]    = _suc_fp
        adata_work.obs["cr2_success_lineage_name"] = _cr2_success_lineage_name
    else:
        _suc_fp = np.full(adata_work.n_obs, np.nan)
        adata_work.obs["cr2_p_oclr_success"]    = np.nan
        adata_work.obs["cr2_success_lineage_name"] = ""

    if _cr2_failure_lineage_name and _cr2_failure_lineage_name in _all_fp_arrays:
        _fail_fp = _all_fp_arrays[_cr2_failure_lineage_name]
        adata_work.obs["cr2_p_oclr_failure"]    = _fail_fp
        adata_work.obs["cr2_failure_lineage_name"] = _cr2_failure_lineage_name
    else:
        _fail_fp = np.full(adata_work.n_obs, np.nan)
        adata_work.obs["cr2_p_oclr_failure"]    = np.nan
        adata_work.obs["cr2_failure_lineage_name"] = ""

    _fin_both = np.isfinite(_suc_fp) & np.isfinite(_fail_fp)
    _margin   = np.where(_fin_both, _suc_fp - _fail_fp, np.nan)
    adata_work.obs["cr2_oclr_margin"] = _margin

    # Legacy aliases
    adata_work.obs["cr2_ips_fate_prob"]    = _suc_fp
    adata_work.obs["cr2_ips_lineage_name"] = _cr2_success_lineage_name or ""

    _cr2_ips_fp = adata_work.obs["cr2_ips_fate_prob"].values.astype(float)
    _fin_cr2    = np.isfinite(_cr2_ips_fp)
    if _fin_cr2.any():
        print(f"\n  cr2_p_oclr_success  "
              f"min={_cr2_ips_fp[_fin_cr2].min():.4f}  "
              f"median={np.median(_cr2_ips_fp[_fin_cr2]):.4f}  "
              f"max={_cr2_ips_fp[_fin_cr2].max():.4f}")
        _fin_m = np.isfinite(_margin)
        if _fin_m.any():
            print(f"  cr2_oclr_margin     "
                  f"min={float(np.nanmin(_margin)):.4f}  "
                  f"median={float(np.nanmedian(_margin)):.4f}  "
                  f"max={float(np.nanmax(_margin)):.4f}")
    else:
        print(f"  [WARN] cr2_p_oclr_success not set — no valid lineage identified.")

    # ── Save OCLR annotation table ────────────────────────────────────────────
    _enrich_path = (
        CR2_OUT_DIR / (TIMESTAMP + f"_{_MODE_TAG}_cr2_oclr_annotation.tsv")
    )
    try:
        if not _oclr_ann_df.empty:
            _oclr_ann_df.to_csv(str(_enrich_path), sep="\t", index=False)
            print(f"\n  OCLR annotation table: {_enrich_path.name}")
    except Exception as _e:
        print(f"  [WARN] Cannot save OCLR annotation table: {_e}")

    # ── UMAP: cr2_oclr_margin vs WOT p_oclr_margin (primary benchmark scores) ─
    try:
        _ncols = 2
        _fig_fp, _axes_fp = plt.subplots(1, _ncols, figsize=(7 * _ncols, 6))
        _axf = list(_axes_fp)

        _cr2_umap_col = ("cr2_oclr_margin" if "cr2_oclr_margin" in adata_work.obs.columns
                         else "cr2_ips_fate_prob")
        _wot_umap_col = ("p_oclr_margin" if "p_oclr_margin" in adata_work.obs.columns
                         else "p_hcipsc")
        sc.pl.umap(adata_work, color=_cr2_umap_col, ax=_axf[0],
                   show=False, size=UMAP_SIZE, alpha=UMAP_ALPHA,
                   cmap="RdBu_r",
                   title=f"CR2 {_cr2_umap_col}\n"
                         f"(suc: '{_cr2_success_lineage_name}'  "
                         f"fail: '{_cr2_failure_lineage_name or 'unresolved'}')")
        if _wot_umap_col in adata_work.obs.columns:
            sc.pl.umap(adata_work, color=_wot_umap_col, ax=_axf[1],
                       show=False, size=UMAP_SIZE, alpha=UMAP_ALPHA,
                       cmap="RdBu_r", title=f"WOT {_wot_umap_col}")
        else:
            _axf[1].set_visible(False)

        plt.suptitle(
            f"Primary benchmark scores: CR2 vs WOT  [{_MODE_TAG}]\n"
            f"margin = p_oclr_success − p_oclr_failure",
            fontsize=12,
        )
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
    # STEP 12  —  Primary benchmark: CR2 cr2_oclr_margin  vs  WOT p_oclr_margin
    # =========================================================================
    # PRIMARY comparison uses the margin scores (success − failure), which are
    # the main benchmark metric for the OCLR dual-endpoint terminal benchmark.
    # Legacy success-only columns (cr2_ips_fate_prob, p_hcipsc) are NOT the
    # primary comparison here; they are retained as aliases only.
    print(f"\n{_SEP}")
    print("STEP 12  —  Primary benchmark: CR2 cr2_oclr_margin  vs  WOT p_oclr_margin")
    print(_SEP)

    _pr = _sr = np.nan
    try:
        from scipy.stats import pearsonr, spearmanr

        # Use margin columns if available; fall back to success-only with a clear label
        _cr2_cmp_col = ("cr2_oclr_margin"  if "cr2_oclr_margin"  in adata_work.obs.columns
                        else "cr2_ips_fate_prob")
        _wot_cmp_col = ("p_oclr_margin"    if "p_oclr_margin"    in adata_work.obs.columns
                        else "p_hcipsc")
        _using_legacy = (_cr2_cmp_col == "cr2_ips_fate_prob" or
                         _wot_cmp_col == "p_hcipsc")
        if _using_legacy:
            print(f"  [NOTE] Margin columns not found; comparing legacy success-only scores "
                  f"({_wot_cmp_col} vs {_cr2_cmp_col}). "
                  f"This is a supplementary comparison, not the primary benchmark.")

        _cr2_p   = adata_work.obs[_cr2_cmp_col].values.astype(float)
        _wot_p   = adata_work.obs[_wot_cmp_col].values.astype(float)
        _all_fin = np.isfinite(_cr2_p) & np.isfinite(_wot_p)
        _n_all   = int(_all_fin.sum())

        _non_ips = (
            (adata_work.obs["stage"] != _ips_stage_label)
            if _ips_stage_label is not None
            else pd.Series(True, index=adata_work.obs.index)
        )
        _mask  = (_non_ips & _all_fin)
        _n_cmp = int(_mask.sum())

        print(f"  Comparing: WOT {_wot_cmp_col}  vs  CR2 {_cr2_cmp_col}")
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
                "wot_score_col", "cr2_score_col",
                "pearson_r_all", "pearson_p_all",
                "spearman_r_all", "spearman_p_all",
                "pearson_r_nonips", "pearson_p_nonips",
                "spearman_r_nonips", "spearman_p_nonips",
                "n_cells_all", "n_cells_nonips",
                "cr2_success_lineage_name", "cr2_failure_lineage_name",
            ],
            "value": [
                _wot_cmp_col, _cr2_cmp_col,
                _pr, _pp, _sr, _sp,
                _pr2, _pp2, _sr2, _sp2,
                _n_all, _n_cmp,
                str(_cr2_success_lineage_name or ""),
                str(_cr2_failure_lineage_name or ""),
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
        _ax_sc.set_xlabel(f"WOT {_wot_cmp_col}", fontsize=11)
        _ax_sc.set_ylabel(f"CR2 {_cr2_cmp_col}", fontsize=11)
        _ax_sc.set_title(
            f"CR2 vs WOT  [{_MODE_TAG}]\n"
            f"Pearson r={_pr:.3f}  Spearman r={_sr:.3f}",
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
    _corr_col = "corr"   # will be updated per method

    # ── Helper: compute_lineage_drivers for one lineage, return sorted DF ─────
    def _run_lineage_drivers(lineage_name: str) -> pd.DataFrame | None:
        try:
            _df = gpcca.compute_lineage_drivers(
                lineages=lineage_name,
                method="fisher",
                use_raw=False,
                return_drivers=True,
            )
            if not isinstance(_df, pd.DataFrame):
                raise TypeError(f"Expected DataFrame, got {type(_df)}")
            # Find correlation column: "<lineage>_corr" pattern
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

    # ── Method A: CellRank2 compute_lineage_drivers ───────────────────────────
    # Use the OCLR-annotated iPSC lineage identified in Step 11.
    # This is a single lineage — no multi-lineage pooling (Scheme B removed).
    if _cr2_ips_lineage_name is not None:
        print(f"  compute_lineage_drivers(lineages='{_cr2_ips_lineage_name}') …")
        _drivers = _run_lineage_drivers(_cr2_ips_lineage_name)
        if _drivers is not None:
            _corr_col = _drivers["_corr_col_used"].iloc[0]
            _drivers  = _drivers.drop(columns=["_corr_col_used"])
            print(f"  Drivers for '{_cr2_ips_lineage_name}': {len(_drivers)} genes  "
                  f"(corr col='{_corr_col}')")

    # ── Method B fallback: manual Pearson(gene, cr2_ips_fate_prob) ───────────
    if _drivers is None:
        print("  Falling back to manual Pearson correlation with cr2_ips_fate_prob …")
        if "cr2_ips_fate_prob" not in adata_work.obs.columns:
            print("  [WARN] cr2_ips_fate_prob not in adata_work.obs. Skipping.")
        else:
            _fate_vec = adata_work.obs["cr2_ips_fate_prob"].values.astype(np.float64)
            _ok       = np.isfinite(_fate_vec)
            _fate_f   = _fate_vec[_ok]
            _Xsub     = adata_work.X[_ok]
            if sp.issparse(_Xsub):
                _Xsub = _Xsub.toarray()
            _Xsub  = np.asarray(_Xsub, dtype=np.float64)
            _fz    = (_fate_f - _fate_f.mean()) / (_fate_f.std() + 1e-12)
            _n_g   = _Xsub.shape[1]
            _corrs = np.full(_n_g, np.nan, dtype=np.float64)
            for _j in range(_n_g):
                _gv = _Xsub[:, _j]
                _gs = _gv.std()
                if _gs < 1e-12:
                    continue
                _corrs[_j] = float(
                    np.dot(_fz, (_gv - _gv.mean()) / _gs)
                ) / len(_fate_f)
            _gene_names = np.asarray(adata_work.var_names, dtype=str)
            _corr_col   = "manual_corr"
            _drivers    = (
                pd.DataFrame({"manual_corr": _corrs}, index=_gene_names)
                .dropna()
                .sort_values("manual_corr", ascending=False)
                .head(N_DRIVER_GENES)
            )
            print(f"  Manual correlation done: {len(_drivers)} genes")

    if _drivers is not None:
        _drv_path = CR2_OUT_DIR / (TIMESTAMP + f"_{_MODE_TAG}_cr2_driver_genes.tsv")
        _drivers.to_csv(str(_drv_path), sep="\t")
        print(f"  Saved: {_drv_path.name}")
        print(f"\n  Top 20 driver genes (by {_corr_col}):")
        print(_drivers.head(20).to_string())

        # ── Overlap with WOT driver genes ──────────────────────────────────────
        _wot_files = sorted(
            list(METRICS_DIR.glob("*driver_genes_overlap.tsv")) +
            list((PROJECT_ROOT / "results" / "tables").glob("*driver*genes*.tsv")),
            key=lambda p: p.stat().st_mtime,
        )
        if _wot_files:
            _wot_drv = pd.read_csv(str(_wot_files[-1]), sep="\t")
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
                _ovlp_df.to_csv(str(_ovlp_path), sep="\t")
                print(f"  Saved: {_ovlp_path.name}")
        else:
            print("  [INFO] WOT driver table not found — skipping overlap.")

        # ── Driver gene expression trends (top 10) ─────────────────────────────
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
                        # .values converts pandas Series → numpy bool array;
                        # sparse matrix indexing (.nonzero()) requires numpy,
                        # not a pandas Series — this was the root cause of
                        # "'Series' object has no attribute 'nonzero'".
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

    _cr2_ips_fp_vals = adata_work.obs.get(
        "cr2_p_oclr_success",
        adata_work.obs.get("cr2_ips_fate_prob", pd.Series(dtype=float))
    ).values
    _fin_cr2 = np.isfinite(_cr2_ips_fp_vals.astype(float))
    _cr2_suc_fp_median = (
        float(np.median(_cr2_ips_fp_vals[_fin_cr2])) if _fin_cr2.any() else np.nan
    )
    _cr2_margin_vals = adata_work.obs.get(
        "cr2_oclr_margin", pd.Series(np.nan, index=adata_work.obs.index)
    ).values.astype(float)
    _fin_m = np.isfinite(_cr2_margin_vals)
    _cr2_margin_median = float(np.nanmedian(_cr2_margin_vals[_fin_m])) if _fin_m.any() else np.nan

    _rows = [
        ("script",                    "04_cellrank2.py v12 (standard GPCCA + OCLR dual-endpoint annotation)"),
        ("run_mode",                  RUN_MODE),
        ("mode_tag",                  _MODE_TAG),
        ("is_final_result",           "YES — hpc_full full-dataset"
                                      if RUN_MODE == "hpc_full"
                                      else "NO — local_debug subset (pipeline smoke-test only)"),
        ("timestamp",                 TIMESTAMP),
        ("input_h5ad",                fates_h5ad.name),
        ("tmaps_dir_original",        str(tmaps_dir)),
        ("tmaps_dir_work",            str(tmaps_dir_work)),
        ("n_tmap_files_work",         len(_tmap_files_work)),
        ("n_cells_full",              adata.n_obs),
        ("n_cells_work",              adata_work.n_obs),
        ("n_genes",                   adata_work.n_vars),
        ("n_timepoints",              len(_day_cats_work)),
        ("ips_stage_label",           str(_ips_stage_label)),
        ("conn_weight",               CONN_WEIGHT),
        ("threshold",                 THRESHOLD),
        ("n_macrostates",             N_MACROSTATES),
        ("all_macrostate_names",      str(_ms_names)),
        ("n_terminal_lineages",       len(_fp_names)),
        ("terminal_state_method",     "predict_terminal_states() — standard GPCCA"),
        ("cr2_success_lineage_name",  str(_cr2_success_lineage_name)),
        ("cr2_failure_lineage_name",  str(_cr2_failure_lineage_name)),
        ("cr2_ips_lineage_name",      str(_cr2_ips_lineage_name)),
        ("cr2_ips_oclr_spearman_r",   f"{_cr2_ips_oclr_spearman_r:.4f}"),
        ("cr2_p_oclr_success_median", f"{_cr2_suc_fp_median:.4f}"),
        ("cr2_oclr_margin_median",    f"{_cr2_margin_median:.4f}"),
        ("fate_prob_n_jobs",          str(FATE_PROB_N_JOBS)),
        ("n_driver_genes",            len(_drivers) if _drivers is not None else "N/A"),
        ("driver_corr_col",           _corr_col),
        ("pearson_r_all",             f"{_pr:.4f}"),
        ("spearman_r_all",            f"{_sr:.4f}"),
        ("oclr_annotation_tsv",       _enrich_path.name),
        ("petsc_ok",                  str(_petsc_ok)),
        ("slepc_ok",                  str(_slepc_ok)),
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
    print(f"04_cellrank2.py v12  [{_MODE_TAG}]  COMPLETE  ({_elapsed:.1f} s total)")
    print(_SEP)
    print(f"\n  Outputs in : {CR2_OUT_DIR}")
    print(f"  Figures in : {FIGURES_DIR}")
    print(f"\n  Key outputs:")
    print(f"    {_out_h5ad.name}")
    print(f"    {TIMESTAMP}_{_MODE_TAG}_cr2_driver_genes.tsv")
    print(f"    {TIMESTAMP}_{_MODE_TAG}_cr2_oclr_annotation.tsv")
    print(f"    {TIMESTAMP}_{_MODE_TAG}_cr2_summary.tsv")
    print(f"\n  OCLR dual-endpoint annotation:")
    print(f"    SUCCESS lineage : '{_cr2_success_lineage_name}'")
    print(f"    FAILURE lineage : '{_cr2_failure_lineage_name}'")
    print(f"    OCLR Spearman r (success lineage) = {_cr2_ips_oclr_spearman_r:.4f}")
    print(f"    cr2_oclr_margin median = {_cr2_margin_median:.4f}")

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
