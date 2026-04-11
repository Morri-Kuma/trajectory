#!/usr/bin/env python3
"""
03_driver_genes.py  —  WOT driver gene analysis  [v2]
======================================================
Project : Comparative Study of Trajectory Inference Models for Chemical iPSC
          Reprogramming
Dataset : GSE230659 (Liuyang et al. 2023 Cell Stem Cell)
Step    : 03 — Identify genes that drive cells toward iPSC fate

Changelog v2
------------
Fix 1 — Analysis C now uses the exact same terminal iPSC fate definition
         as 02c_compute_transport.py.
         Old (wrong): marker_ips_subset_mask == 1  (global; includes markers
                      outside day30, inflating the group to ~4,680 cells)
         New (correct): fate_ips == 1  (day30 AND marker-positive; the actual
                        WOT fate target, ~4,252 cells)
         Backward-compat fallback: if fate_ips is absent but
         marker_ips_subset_mask is present, the mask is intersected with
         the terminal-pool mask to recover the same cells before running DE.

Fix 2 — Trajectory label names updated to match 02c v10+:
         Old: IPS_seed, hCiPSC_other
         New: marker_ips_subset, terminal_pool_other
         Backward-compatible: code checks for both sets of labels and uses
         whichever is actually present in the data.

Fix 3 — UMAP fallback: if X_umap is absent from the saved h5ad but X_pca is
         available, neighbors and UMAP are computed locally so that Figure 5
         is never silently skipped.

Fix 4 — Filtered consensus driver gene list: a second overlap table is saved
         that excludes ribosomal genes (RPL* / RPS*).  Figures use the
         filtered list when it has ≥5 genes; the unfiltered list is always
         preserved.

What this script does
---------------------
1. Loads wot_fates.h5ad produced by 02c (p_iPSC, right_traj,
   fate_ips, fate_other, traj_stage_seed_failgrey, X_umap).
2. Runs three complementary DE analyses:
     Analysis A — continuous: top-10% p_iPSC vs bottom-10% p_iPSC
                  (non-terminal cells only)
     Analysis B — binary:     right_traj==1 vs right_traj==0
                  (non-terminal cells only)
     Analysis C — identity:   fate_ips (exact 02c terminal set) vs early cells
                  (sensitivity check: genes specific to the WOT iPSC endpoint)
3. Saves ranked TSV tables for each analysis.
4. Saves filtered consensus table (A∩B minus ribosomal genes).
5. Produces publication-ready figures:
     • Dot plot   — top 30 driver genes across trajectory stages
     • Heatmap    — top 50 genes (high vs low p_iPSC)
     • Violin     — top 10 driver genes by stage
     • Scatter    — p_iPSC vs top driver gene expression
     • UMAP       — coloured by p_iPSC, trajectory class, top driver genes

Usage
-----
    conda activate traj_env
    cd C:\\Users\\37620\\trajectory
    python scripts/WOT/03_driver_genes.py
"""

# =============================================================================
# 0.  CONFIGURATION
# =============================================================================
from datetime import datetime
from pathlib import Path
import time as _time

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")
T0        = _time.time()

_candidates = [
    Path(r"C:\Users\37620\trajectory"),
    Path("/sessions/cool-admiring-hawking/mnt/trajectory"),
    Path(__file__).resolve().parents[2],
]
PROJECT_ROOT = next((p for p in _candidates if p.exists()), _candidates[-1])
METRICS_DIR  = PROJECT_ROOT / "results" / "metrics"
FIGURES_DIR  = PROJECT_ROOT / "results" / "figures"
TABLES_DIR   = PROJECT_ROOT / "results" / "tables"
for _d in [FIGURES_DIR, TABLES_DIR]:
    _d.mkdir(parents=True, exist_ok=True)

def find_latest(pattern: str, directory: Path) -> Path:
    matches = sorted(directory.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No files matching '{pattern}' in {directory}.\n"
            f"Run 02c_compute_transport.py first."
        )
    return matches[-1]

FATES_PATH = find_latest("*_GSE230659_wot_fates.h5ad", METRICS_DIR)

# ── Analysis parameters ───────────────────────────────────────────────────────
# Analysis A: top vs bottom p_iPSC quantile split
Q_LOW:  float = 0.10   # bottom 10% = "failing" cells
Q_HIGH: float = 0.90   # top 10%    = "succeeding" cells

# Number of top driver genes to save / plot
TOP_N_SAVE:   int = 300
TOP_N_PLOT:   int = 30
TOP_N_VIOLIN: int = 10
TOP_N_UMAP:   int = 6

# DE method
DE_METHOD: str = "wilcoxon"

# UMAP recomputation parameters (used only when X_umap is absent)
UMAP_N_NEIGHBORS: int = 15
UMAP_N_PCS:       int = 50

# Ribosomal gene prefixes to exclude from filtered consensus
RIBO_PREFIXES: tuple[str, ...] = ("RPL", "RPS")

# Stage order for plots (current 02c labels)
STAGE_ORDER = ["StageI", "StageII", "StageIII", "hCiPSCs"]

STAGE_COLORS = {
    "StageI":   "#2196F3",
    "StageII":  "#FF9800",
    "StageIII": "#4CAF50",
    "hCiPSCs":  "#E91E63",
}

# Trajectory category palette — matches 02c TRAJ_PALETTE (current label names)
# Backward-compatible: old label names are listed as aliases below.
TRAJ_COLORS = {
    # Current labels (02c v10+)
    "marker_ips_subset":   "#d62728",   # strong red  — fate_ips terminal cells
    "terminal_pool_other": "#f4a3a3",   # light pink  — fate_other terminal cells
    "StageIII":            "#2ca02c",   # green
    "StageII":             "#ff7f0e",   # orange
    "StageI":              "#1f77b4",   # blue
    "Fail":                "#BDBDBD",   # grey
    # Legacy aliases (02c v9 and earlier) — kept for backward compatibility
    "IPS_seed":            "#d62728",
    "hCiPSC_other":        "#f4a3a3",
}

# Preferred draw order: Fail first (background), marker_ips_subset last (foreground)
# Both current and legacy labels included so either set is handled gracefully.
_TRAJ_ORDER_CURRENT = [
    "Fail", "StageI", "StageII", "StageIII",
    "terminal_pool_other", "marker_ips_subset",
]
_TRAJ_ORDER_LEGACY  = [
    "Fail", "StageI", "StageII", "StageIII",
    "hCiPSC_other", "IPS_seed",
]

# =============================================================================
# 1.  IMPORTS
# =============================================================================
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sc.settings.verbosity = 1
sc.settings.figdir    = str(FIGURES_DIR)
sc.settings.autoshow  = False

_SEP  = "─" * 70
_SEP2 = "· " * 35

print(f"[{TIMESTAMP}]  WOT Driver Gene Analysis  (Step 03 v2)")
print(f"  Input  : {FATES_PATH.name}")
print(f"  Tables → {TABLES_DIR}")
print(f"  Figures→ {FIGURES_DIR}")

# =============================================================================
# 2.  LOAD wot_fates.h5ad AND VALIDATE
# =============================================================================
print(f"\n{_SEP}")
print("STEP 2  —  Load wot_fates AnnData")
print(_SEP)

adata = sc.read_h5ad(str(FATES_PATH))
print(f"  Shape  : {adata.n_obs:,} cells × {adata.n_vars:,} genes")
print(f"  obs    : {list(adata.obs.columns)}")
print(f"  obsm   : {list(adata.obsm.keys())}")

# ── Validate required columns ──────────────────────────────────────────────────
REQUIRED = ["p_iPSC", "right_traj", "stage", "day"]
_missing  = [c for c in REQUIRED if c not in adata.obs.columns]
if _missing:
    raise ValueError(
        f"Required obs columns missing: {_missing}\n"
        f"Re-run 02c_compute_transport.py to regenerate wot_fates.h5ad."
    )
print(f"  Required columns OK : {REQUIRED}")

# ── Detect optional fate columns written by 02c ───────────────────────────────
_has_fate_ips           = "fate_ips"              in adata.obs.columns
_has_fate_other         = "fate_other"            in adata.obs.columns
_has_marker_subset_mask = "marker_ips_subset_mask" in adata.obs.columns

print(f"  fate_ips column             : {'present' if _has_fate_ips else 'absent'}")
print(f"  fate_other column           : {'present' if _has_fate_other else 'absent'}")
print(f"  marker_ips_subset_mask col  : {'present' if _has_marker_subset_mask else 'absent'}")

# ── Report fate_ips / fate_other sizes ────────────────────────────────────────
if _has_fate_ips:
    _n_fate_ips   = int((adata.obs["fate_ips"].astype(int) == 1).sum())
    print(f"  fate_ips  cells (day30 + marker+) : {_n_fate_ips:,}")
if _has_fate_other:
    _n_fate_other = int((adata.obs["fate_other"].astype(int) == 1).sum())
    print(f"  fate_other cells (day30 + marker-): {_n_fate_other:,}")

# ── Validate p_iPSC ───────────────────────────────────────────────────────────
_p      = adata.obs["p_iPSC"].values.astype(float)
_finite = np.isfinite(_p)
print(f"\n  p_iPSC : {_finite.sum():,} finite ({_finite.mean():.1%})  "
      f"min={np.nanmin(_p):.4f}  median={np.nanmedian(_p):.4f}  "
      f"max={np.nanmax(_p):.4f}")

if _finite.mean() < 0.5:
    raise ValueError(
        "Fewer than 50% of cells have finite p_iPSC. "
        "Check wot_fates.h5ad — fate computation may have failed."
    )

# ── Build terminal-pool mask (ALL hCiPSC/day30 cells) ─────────────────────────
# Defined as the stage-level group.  Used to exclude terminal cells from
# Analyses A and B (they ARE the endpoint; including them confounds DE).
_terminal_pool_mask = adata.obs["stage"].str.lower().str.contains(
    "cipsc|ipsc", na=False
)
print(f"  Terminal reference pool (ALL hCiPSC/day30) : {_terminal_pool_mask.sum():,} cells")

# If fate_ips / fate_other are present we can also verify the partition
if _has_fate_ips and _has_fate_other:
    _fate_union = (
        (adata.obs["fate_ips"].astype(int) == 1) |
        (adata.obs["fate_other"].astype(int) == 1)
    )
    _diff = int(_fate_union.sum()) - int(_terminal_pool_mask.sum())
    if abs(_diff) > 10:
        print(f"  [WARN] fate_ips + fate_other = {_fate_union.sum():,}  "
              f"vs terminal_pool_mask = {_terminal_pool_mask.sum():,}  "
              f"(diff={_diff}; minor discrepancy expected if stages overlap)")

# ── Build adata_de: non-terminal cells with finite p_iPSC ─────────────────────
adata_de = adata[~_terminal_pool_mask & _finite].copy()
print(f"  Cells available for DE (non-terminal, finite p) : {adata_de.n_obs:,}")

# ── Check X looks like log1p expression ───────────────────────────────────────
_x_sample = (adata_de.X[:200].toarray() if sp.issparse(adata_de.X)
             else np.asarray(adata_de.X[:200]))
if _x_sample.min() < 0 or _x_sample.max() > 50:
    print(f"  [WARN] adata.X values look unexpected for log1p expression: "
          f"min={_x_sample.min():.3f}  max={_x_sample.max():.3f}")
else:
    print(f"  X range check OK : min={_x_sample.min():.3f}  "
          f"max={_x_sample.max():.3f}")

# =============================================================================
# 3.  ANALYSIS A — high p_iPSC (top 10%) vs low p_iPSC (bottom 10%)
# =============================================================================
print(f"\n{_SEP}")
print("STEP 3  —  Analysis A: high vs low p_iPSC  (Wilcoxon DE)")
print(_SEP)

_p_de   = adata_de.obs["p_iPSC"].values.astype(float)
_lo_thr = float(np.nanquantile(_p_de, Q_LOW))
_hi_thr = float(np.nanquantile(_p_de, Q_HIGH))

print(f"  Quantile thresholds : low ≤ {_lo_thr:.4f} (Q{Q_LOW:.0%})  "
      f"| high ≥ {_hi_thr:.4f} (Q{Q_HIGH:.0%})")

_grp = np.full(adata_de.n_obs, "mid", dtype=object)
_grp[_p_de <= _lo_thr] = "low"
_grp[_p_de >= _hi_thr] = "high"
adata_de.obs["driver_group"] = _grp

_n_low  = int((_grp == "low").sum())
_n_high = int((_grp == "high").sum())
_n_mid  = int((_grp == "mid").sum())
print(f"  Group sizes : low={_n_low:,}  mid={_n_mid:,}  high={_n_high:,}")

if _n_low < 50 or _n_high < 50:
    print("  [WARN] Group too small — relaxing to Q0.20 / Q0.80")
    Q_LOW, Q_HIGH = 0.20, 0.80
    _lo_thr = float(np.nanquantile(_p_de, Q_LOW))
    _hi_thr = float(np.nanquantile(_p_de, Q_HIGH))
    _grp[_p_de <= _lo_thr] = "low"
    _grp[_p_de >= _hi_thr] = "high"
    adata_de.obs["driver_group"] = _grp
    _n_low  = int((_grp == "low").sum())
    _n_high = int((_grp == "high").sum())
    print(f"  Adjusted group sizes : low={_n_low:,}  high={_n_high:,}")

_adata_a = adata_de[adata_de.obs["driver_group"].isin(["high", "low"])].copy()

print(f"  Running rank_genes_groups (method={DE_METHOD}) … ", end="", flush=True)
sc.tl.rank_genes_groups(
    _adata_a,
    groupby   = "driver_group",
    groups    = ["high"],
    reference = "low",
    method    = DE_METHOD,
    n_genes   = adata_de.n_vars,
)
print("done")


def _extract_de(adata_r, group: str) -> pd.DataFrame:
    """Extract rank_genes_groups results into a clean DataFrame."""
    res       = adata_r.uns["rank_genes_groups"]
    names     = np.asarray(res["names"][group])
    scores    = np.asarray(res["scores"][group], dtype=float)
    pvals_adj = np.asarray(res["pvals_adj"][group], dtype=float)
    try:
        logfc = np.asarray(res["logfoldchanges"][group], dtype=float)
    except (KeyError, TypeError):
        logfc = np.full(len(names), np.nan)
    df = pd.DataFrame({
        "gene":     names,
        "score":    scores,
        "pvals_adj": pvals_adj,
        "logFC":    logfc,
    })
    return df.dropna(subset=["gene"]).drop_duplicates("gene")


df_a = _extract_de(_adata_a, "high")
print(f"  Genes ranked : {len(df_a):,}")
print(f"\n  Top 15 driver genes (high p_iPSC):")
print(df_a.head(15).to_string(index=False))

_path_a = TABLES_DIR / f"{TIMESTAMP}_driver_genes_pIPSC_highvslow.tsv"
df_a.head(TOP_N_SAVE).to_csv(_path_a, sep="\t", index=False)
print(f"\n  Saved → {_path_a.name}  (top {TOP_N_SAVE} genes)")

TOP_GENES_A = df_a["gene"].head(TOP_N_PLOT).tolist()

# =============================================================================
# 4.  ANALYSIS B — right_traj == 1 vs right_traj == 0
# =============================================================================
print(f"\n{_SEP}")
print("STEP 4  —  Analysis B: right_traj 1 vs 0  (binary DE)")
print(_SEP)

_rt = adata_de.obs["right_traj"].values.astype(int)
_n1 = int((_rt == 1).sum())
_n0 = int((_rt == 0).sum())
print(f"  right_traj=1 : {_n1:,}   right_traj=0 : {_n0:,}")

_adata_b = adata_de.copy()
_adata_b.obs["traj_bin"] = np.where(_rt == 1, "on_traj", "off_traj")

print(f"  Running rank_genes_groups (method={DE_METHOD}) … ", end="", flush=True)
sc.tl.rank_genes_groups(
    _adata_b,
    groupby   = "traj_bin",
    groups    = ["on_traj"],
    reference = "off_traj",
    method    = DE_METHOD,
    n_genes   = adata_de.n_vars,
)
print("done")

df_b = _extract_de(_adata_b, "on_traj")
print(f"  Genes ranked : {len(df_b):,}")
print(f"\n  Top 15 driver genes (right_traj=1):")
print(df_b.head(15).to_string(index=False))

_path_b = TABLES_DIR / f"{TIMESTAMP}_driver_genes_righttraj_1vs0.tsv"
df_b.head(TOP_N_SAVE).to_csv(_path_b, sep="\t", index=False)
print(f"\n  Saved → {_path_b.name}  (top {TOP_N_SAVE} genes)")

TOP_GENES_B = df_b["gene"].head(TOP_N_PLOT).tolist()

# Union of top genes from both analyses (for some plots)
TOP_GENES_UNION = list(dict.fromkeys(TOP_GENES_A[:20] + TOP_GENES_B[:20]))[:TOP_N_PLOT]

# =============================================================================
# 5.  OVERLAP TABLE — genes in both Analyses A and B
# =============================================================================
print(f"\n{_SEP}")
print("STEP 5  —  Overlap between Analysis A and B")
print(_SEP)

_set_a   = set(df_a["gene"].head(TOP_N_SAVE))
_set_b   = set(df_b["gene"].head(TOP_N_SAVE))
_overlap = _set_a & _set_b
print(f"  Top-{TOP_N_SAVE} A only   : {len(_set_a - _set_b)}")
print(f"  Top-{TOP_N_SAVE} B only   : {len(_set_b - _set_a)}")
print(f"  In both (overlap) : {len(_overlap)}")

_df_overlap = (
    df_a[df_a["gene"].isin(_overlap)][["gene", "score", "logFC"]]
    .rename(columns={"score": "score_A", "logFC": "logFC_A"})
    .merge(
        df_b[df_b["gene"].isin(_overlap)][["gene", "score", "logFC"]]
        .rename(columns={"score": "score_B", "logFC": "logFC_B"}),
        on="gene",
    )
    .assign(rank_sum=lambda d: d["score_A"] + d["score_B"])
    .sort_values("rank_sum", ascending=False)
    .reset_index(drop=True)
)

_path_overlap = TABLES_DIR / f"{TIMESTAMP}_driver_genes_overlap.tsv"
_df_overlap.to_csv(_path_overlap, sep="\t", index=False)
print(f"\n  Top overlapping genes (unfiltered):")
print(_df_overlap.head(20).to_string(index=False))
print(f"\n  Saved → {_path_overlap.name}")

# ── Filtered consensus: exclude ribosomal genes (RPL* / RPS*) ─────────────────
# Rationale: the raw A∩B overlap is often dominated by ribosomal proteins
# (RPL41, RPS10, RPL8, etc.), which reflect general translational activity
# rather than specific pluripotency drivers.  A second filtered table removes
# these for biological interpretation; the unfiltered table is always preserved.
def _is_ribo(gene: str) -> bool:
    return any(gene.upper().startswith(p) for p in RIBO_PREFIXES)

_df_overlap_filtered = _df_overlap[
    ~_df_overlap["gene"].apply(_is_ribo)
].reset_index(drop=True)

_n_ribo_removed = len(_df_overlap) - len(_df_overlap_filtered)
print(f"\n  Ribosomal genes removed from filtered consensus : {_n_ribo_removed}")
print(f"  Filtered overlap size                          : {len(_df_overlap_filtered)}")

_path_overlap_filt = TABLES_DIR / f"{TIMESTAMP}_driver_genes_overlap_filtered.tsv"
_df_overlap_filtered.to_csv(_path_overlap_filt, sep="\t", index=False)
print(f"  Saved → {_path_overlap_filt.name}")

if len(_df_overlap_filtered) > 0:
    print(f"\n  Top filtered consensus driver genes:")
    print(_df_overlap_filtered.head(20).to_string(index=False))

# Choose which gene list drives the main figures
# Use filtered list if it has enough coverage; otherwise fall back to unfiltered.
_USE_FILTERED_FOR_PLOTS = len(_df_overlap_filtered) >= 5

TOP_GENES_CONSENSUS = (
    _df_overlap_filtered["gene"].head(TOP_N_PLOT).tolist()
    if _USE_FILTERED_FOR_PLOTS
    else _df_overlap["gene"].head(TOP_N_PLOT).tolist()
)
if len(TOP_GENES_CONSENSUS) < 10:
    # Fall back to Analysis A if overlap is too small
    TOP_GENES_CONSENSUS = TOP_GENES_A[:TOP_N_PLOT]

print(f"\n  Gene list used for figures: "
      f"{'filtered consensus (ribo excluded)' if _USE_FILTERED_FOR_PLOTS else 'unfiltered consensus'}"
      f"  ({len(TOP_GENES_CONSENSUS)} genes)")

# =============================================================================
# 5b.  ANALYSIS C — fate_ips (exact 02c terminal set) vs early cells
# =============================================================================
# DESIGN NOTE
# -----------
# In 02c, the WOT two-fate design defines:
#   fate_ips   = day30/hCiPSC cells AND marker-positive (POU5F1/SOX2/NANOG high)
#   fate_other = day30/hCiPSC cells AND marker-negative
# fate_ips is the actual WOT fate target — the column that WOT backward-
# propagates to give p_iPSC.
#
# Analysis C compares fate_ips cells against early (non-terminal) cells.
# This answers: "which genes distinguish the confirmed iPSC endpoint from
# early reprogramming intermediates?"
#
# IMPORTANT: we must use fate_ips == 1, NOT marker_ips_subset_mask == 1 alone.
# The global marker_ips_subset_mask can include marker-positive cells outside
# the day30 stage, inflating the positive group.  fate_ips is already
# restricted to day30 AND marker-positive, matching the 02c WOT design.
# =============================================================================
print(f"\n{_SEP}")
print("STEP 5b  —  Analysis C: fate_ips vs early cells  (identity / sensitivity)")
print(_SEP)

df_c        = None
TOP_GENES_C = []
_path_c     = None

# ── Resolve the positive group for Analysis C ──────────────────────────────────
if _has_fate_ips:
    # Preferred path: use the exact fate_ips column from 02c
    _mask_c_pos = adata.obs["fate_ips"].astype(int) == 1
    _c_source   = "fate_ips column (exact 02c WOT terminal set)"

elif _has_marker_subset_mask:
    # Backward-compat path: restrict global marker mask to terminal pool
    # so we recover the same set of cells as fate_ips would give.
    _mask_c_pos = (
        (adata.obs["marker_ips_subset_mask"].astype(int) == 1) &
        _terminal_pool_mask
    )
    _c_source = (
        "marker_ips_subset_mask ∩ terminal_pool (backward-compat fallback; "
        "fate_ips column absent — re-run 02c for exact alignment)"
    )

else:
    _mask_c_pos = None
    _c_source   = "unavailable"

_mask_c_ref = ~_terminal_pool_mask & _finite   # early cells = reference

if _mask_c_pos is None:
    print("  Neither fate_ips nor marker_ips_subset_mask found — Analysis C skipped.")
    print("  Re-run 02c_compute_transport.py to generate fate_ips.")

else:
    _n_pos = int(_mask_c_pos.sum())
    _n_ref = int(_mask_c_ref.sum())
    print(f"  Positive group : {_n_pos:,} cells  ({_c_source})")
    print(f"  Reference group: {_n_ref:,} early cells")

    if _n_pos < 20:
        print(f"  [WARN] Fewer than 20 cells in positive group — Analysis C skipped.")

    else:
        _adata_c = adata[_mask_c_pos | _mask_c_ref].copy()

        # Label column: "fate_ips_cell" | "early"
        _grp_c = np.where(
            (_adata_c.obs["fate_ips"].astype(int) == 1
             if _has_fate_ips
             else ((_adata_c.obs["marker_ips_subset_mask"].astype(int) == 1) &
                   _adata_c.obs["stage"].str.lower().str.contains("cipsc|ipsc", na=False))),
            "fate_ips_cell",
            "early",
        )
        _adata_c.obs["c_group"] = _grp_c

        print(f"  Running rank_genes_groups (method={DE_METHOD}) … ",
              end="", flush=True)
        sc.tl.rank_genes_groups(
            _adata_c,
            groupby   = "c_group",
            groups    = ["fate_ips_cell"],
            reference = "early",
            method    = DE_METHOD,
            n_genes   = adata.n_vars,
        )
        print("done")

        df_c = _extract_de(_adata_c, "fate_ips_cell")
        print(f"  Genes ranked : {len(df_c):,}")
        print(f"\n  Top 15 driver genes (fate_ips vs early):")
        print(df_c.head(15).to_string(index=False))

        _path_c = TABLES_DIR / f"{TIMESTAMP}_driver_genes_fateips_vs_early.tsv"
        df_c.head(TOP_N_SAVE).to_csv(_path_c, sep="\t", index=False)
        print(f"\n  Saved → {_path_c.name}  (top {TOP_N_SAVE} genes)")

        TOP_GENES_C = df_c["gene"].head(TOP_N_PLOT).tolist()

        # Cross-overlap: genes robust across trajectory AND identity analysis
        _overlap_ac = (
            set(df_a["gene"].head(TOP_N_SAVE)) & set(df_c["gene"].head(TOP_N_SAVE))
        )
        _overlap_ac_filt = {g for g in _overlap_ac if not _is_ribo(g)}
        print(f"\n  A ∩ C overlap (total): {len(_overlap_ac)} genes")
        print(f"  A ∩ C overlap (non-ribosomal): {len(_overlap_ac_filt)} genes")
        if _overlap_ac_filt:
            print(f"  → {sorted(_overlap_ac_filt)[:30]}")

# =============================================================================
# 6.  UMAP — ensure X_umap is available
# =============================================================================
# Fix 3: if X_umap is absent from the saved h5ad (e.g. because 02c was run
# before UMAP computation was added, or the file was written without obsm),
# compute neighbors + UMAP locally so Figure 5 is not silently skipped.
print(f"\n{_SEP}")
print("STEP 6  —  UMAP availability check")
print(_SEP)

if "X_umap" in adata.obsm:
    print("  X_umap found in obsm — reusing existing UMAP coordinates.")
    _umap_recomputed = False

elif "X_pca" in adata.obsm:
    print(f"  X_umap absent but X_pca present — computing UMAP locally.")
    print(f"  sc.pp.neighbors(n_neighbors={UMAP_N_NEIGHBORS}, "
          f"use_rep='X_pca') …", end="", flush=True)
    sc.pp.neighbors(adata, n_neighbors=UMAP_N_NEIGHBORS, use_rep="X_pca",
                    n_pcs=min(UMAP_N_PCS, adata.obsm["X_pca"].shape[1]))
    print(" done")
    print("  sc.tl.umap() …", end="", flush=True)
    sc.tl.umap(adata)
    print(" done")
    _umap_recomputed = True
    print(f"  UMAP recomputed: {adata.obsm['X_umap'].shape}")

else:
    print("  [WARN] Neither X_umap nor X_pca found in obsm.")
    print("  Figure 5 (UMAP) will be skipped.")
    _umap_recomputed = False

# =============================================================================
# 7.  FIGURES
# =============================================================================
print(f"\n{_SEP}")
print("STEP 7  —  Figures")
print(_SEP)


def _ordered_stage_obs(ad: sc.AnnData) -> sc.AnnData:
    """Add stage_ord as ordered categorical for consistent plot ordering."""
    _present = [s for s in STAGE_ORDER if s in ad.obs["stage"].values]
    ad.obs["stage_ord"] = pd.Categorical(
        ad.obs["stage"].values, categories=_present, ordered=True
    )
    return ad


adata = _ordered_stage_obs(adata)

# ── Resolve trajectory column and category order ────────────────────────────────
# Try current label names first; fall back to legacy names if absent.
_traj_col = (
    "traj_stage_seed_failgrey"
    if "traj_stage_seed_failgrey" in adata.obs.columns
    else "stage"
)

def _resolve_traj_order(col: str) -> list[str]:
    """
    Return ordered category list for the trajectory column.
    Checks for current labels first, then legacy labels, then uses whatever
    categories are actually present.
    """
    _present_vals = set(adata.obs[col].values)
    # Current label set
    _order = [c for c in _TRAJ_ORDER_CURRENT if c in _present_vals]
    if len(_order) >= 2:
        return _order
    # Legacy label set
    _order = [c for c in _TRAJ_ORDER_LEGACY if c in _present_vals]
    if len(_order) >= 2:
        return _order
    # Fallback: alphabetical
    return sorted(_present_vals)

_traj_order = _resolve_traj_order(_traj_col)
_traj_label_source = (
    "current (marker_ips_subset / terminal_pool_other)"
    if "marker_ips_subset" in _traj_order
    else ("legacy (IPS_seed / hCiPSC_other)" if "IPS_seed" in _traj_order
          else "stage column fallback")
)
print(f"  Trajectory column    : {_traj_col}")
print(f"  Trajectory label set : {_traj_label_source}")
print(f"  Category order       : {_traj_order}")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 1 — p_iPSC distribution by stage (violin + median trajectory)
# ─────────────────────────────────────────────────────────────────────────────
_fig1, _axes1 = plt.subplots(1, 2, figsize=(14, 5))

_stage_vals = {
    s: adata.obs.loc[adata.obs["stage"] == s, "p_iPSC"].dropna().values
    for s in STAGE_ORDER if s in adata.obs["stage"].values
}
_vp = _axes1[0].violinplot(
    list(_stage_vals.values()), showmedians=True, showextrema=True
)
for _pc, _s in zip(_vp["bodies"], _stage_vals.keys()):
    _pc.set_facecolor(STAGE_COLORS.get(_s, "#999"))
    _pc.set_alpha(0.7)
_axes1[0].set_xticks(range(1, len(_stage_vals) + 1))
_axes1[0].set_xticklabels(list(_stage_vals.keys()), rotation=30, ha="right")
_axes1[0].set_ylabel("p(iPSC fate)")
_axes1[0].set_title("p_iPSC distribution by stage")

if "day" in adata.obs.columns:
    _day_med = (
        adata.obs[np.isfinite(adata.obs["p_iPSC"].values)]
        .groupby(["day", "stage"], observed=True)["p_iPSC"]
        .median()
        .reset_index()
    )
    for _s, _grp in _day_med.groupby("stage"):
        _axes1[1].scatter(
            _grp["day"], _grp["p_iPSC"],
            color=STAGE_COLORS.get(_s, "#999"),
            s=60, label=_s, zorder=3,
        )
    _axes1[1].set_xlabel("Day")
    _axes1[1].set_ylabel("Median p(iPSC fate)")
    _axes1[1].set_title("p_iPSC median trajectory over time")
    _axes1[1].legend(fontsize=8)

_fig1.suptitle("WOT Fate Probability — iPSC Reprogramming  (GSE230659)",
               fontsize=12, fontweight="bold")
_fig1.tight_layout()
_p1 = FIGURES_DIR / f"{TIMESTAMP}_piPSC_by_stage_day.png"
_fig1.savefig(_p1, dpi=150, bbox_inches="tight")
plt.close(_fig1)
print(f"  Fig 1 (p_iPSC by stage/day)  → {_p1.name}")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 2 — Dot plot: top driver genes across trajectory categories
# ─────────────────────────────────────────────────────────────────────────────
_genes_in_var = [g for g in TOP_GENES_CONSENSUS if g in adata.var_names]
if len(_genes_in_var) < 5:
    _genes_in_var = [g for g in TOP_GENES_A if g in adata.var_names][:TOP_N_PLOT]

if _genes_in_var and len(_traj_order) >= 2:
    try:
        _adata_dot = adata.copy()
        _adata_dot.obs[_traj_col] = pd.Categorical(
            _adata_dot.obs[_traj_col], categories=_traj_order, ordered=True
        )
        sc.pl.dotplot(
            _adata_dot,
            var_names        = _genes_in_var[:TOP_N_PLOT],
            groupby          = _traj_col,
            categories_order = _traj_order,
            standard_scale   = "var",
            show             = False,
            save             = f"_{TIMESTAMP}_driver_dotplot.png",
            title            = ("Top driver genes (consensus A∩B, ribo-excluded)"
                                if _USE_FILTERED_FOR_PLOTS
                                else "Top driver genes (consensus A∩B)"),
            figsize          = (max(12, len(_genes_in_var) * 0.5), 5),
        )
        print(f"  Fig 2 (dot plot)             → {TIMESTAMP}_driver_dotplot.png")
    except Exception as _e:
        print(f"  Fig 2 (dot plot) skipped: {_e}")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 3 — Heatmap: top 50 genes from Analysis A, high vs low p_iPSC
# ─────────────────────────────────────────────────────────────────────────────
_genes_heat = [g for g in df_a["gene"].head(50).tolist() if g in adata.var_names]
_adata_heat = _adata_a[_adata_a.obs["driver_group"].isin(["high", "low"])].copy()

if _genes_heat and _adata_heat.n_obs > 0:
    try:
        _adata_heat.obs["driver_group"] = pd.Categorical(
            _adata_heat.obs["driver_group"],
            categories=["high", "low"],
            ordered=True,
        )
        sc.pl.heatmap(
            _adata_heat,
            var_names       = _genes_heat,
            groupby         = "driver_group",
            standard_scale  = "var",
            show_gene_labels= True,
            show            = False,
            save            = f"_{TIMESTAMP}_driver_heatmap.png",
            figsize         = (14, 6),
        )
        print(f"  Fig 3 (heatmap)              → {TIMESTAMP}_driver_heatmap.png")
    except Exception as _e:
        print(f"  Fig 3 (heatmap) skipped: {_e}")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 4 — Violin: top 10 driver genes across stages
# ─────────────────────────────────────────────────────────────────────────────
_genes_vln = [g for g in df_a["gene"].head(TOP_N_VIOLIN).tolist() if g in adata.var_names]
if _genes_vln:
    try:
        sc.pl.violin(
            adata,
            keys    = _genes_vln,
            groupby = "stage",
            order   = [s for s in STAGE_ORDER if s in adata.obs["stage"].values],
            rotation= 30,
            show    = False,
            save    = f"_{TIMESTAMP}_driver_violin.png",
        )
        print(f"  Fig 4 (violin)               → {TIMESTAMP}_driver_violin.png")
    except Exception as _e:
        print(f"  Fig 4 (violin) skipped: {_e}")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 5 — UMAP panels (p_iPSC + trajectory class + top driver genes)
# ─────────────────────────────────────────────────────────────────────────────
if "X_umap" in adata.obsm:
    _umap_note = "recomputed locally" if _umap_recomputed else "from h5ad"
    print(f"  Generating UMAP panels ({_umap_note}) …")

    _genes_umap = [g for g in TOP_GENES_CONSENSUS[:TOP_N_UMAP]
                   if g in adata.var_names]
    _umap_keys  = ["p_iPSC", _traj_col] + _genes_umap
    _umap_keys  = [k for k in _umap_keys
                   if k in adata.obs.columns or k in adata.var_names]

    _ncols = 4
    _nrows = int(np.ceil(len(_umap_keys) / _ncols))
    _fig5, _axes5 = plt.subplots(
        _nrows, _ncols, figsize=(_ncols * 4, _nrows * 3.5))
    _axes5 = np.asarray(_axes5).ravel()

    _umap_xy = adata.obsm["X_umap"]

    for _ax_i, _key in enumerate(_umap_keys):
        _ax = _axes5[_ax_i]
        if _key in adata.obs.columns:
            _vals = adata.obs[_key].values
            if pd.api.types.is_numeric_dtype(_vals):
                _vals = pd.to_numeric(_vals, errors="coerce")
                _sc   = _ax.scatter(
                    _umap_xy[:, 0], _umap_xy[:, 1],
                    c=_vals, cmap="viridis", s=1, alpha=0.5, rasterized=True,
                )
                plt.colorbar(_sc, ax=_ax, shrink=0.7)
            else:
                _cats = pd.Categorical(_vals)
                _cmap = plt.cm.get_cmap("tab10", len(_cats.categories))
                for _ci, _cat in enumerate(_cats.categories):
                    _cmask = _cats == _cat
                    _ax.scatter(
                        _umap_xy[_cmask, 0], _umap_xy[_cmask, 1],
                        c=[TRAJ_COLORS.get(str(_cat), _cmap(_ci))],
                        s=1, alpha=0.4, label=str(_cat), rasterized=True,
                    )
                _ax.legend(fontsize=5, markerscale=4,
                           loc="upper right", framealpha=0.5)
        elif _key in adata.var_names:
            _gi  = list(adata.var_names).index(_key)
            _col = (adata.X[:, _gi].toarray().ravel()
                    if sp.issparse(adata.X)
                    else np.asarray(adata.X[:, _gi]).ravel())
            _sc  = _ax.scatter(
                _umap_xy[:, 0], _umap_xy[:, 1],
                c=_col, cmap="Reds", s=1, alpha=0.5, rasterized=True,
            )
            plt.colorbar(_sc, ax=_ax, shrink=0.7)
        _ax.set_title(_key, fontsize=9)
        _ax.axis("off")

    for _ax_i in range(len(_umap_keys), len(_axes5)):
        _axes5[_ax_i].axis("off")

    _umap_title = (
        "UMAP — p_iPSC, trajectory class, top driver genes"
        + ("\n(UMAP recomputed locally — X_umap absent from h5ad)"
           if _umap_recomputed else "")
    )
    _fig5.suptitle(_umap_title, fontsize=11, fontweight="bold")
    _fig5.tight_layout()
    _p5 = FIGURES_DIR / f"{TIMESTAMP}_umap_driver_genes.png"
    _fig5.savefig(_p5, dpi=150, bbox_inches="tight")
    plt.close(_fig5)
    print(f"  Fig 5 (UMAP panels)          → {_p5.name}")
else:
    print("  Fig 5 (UMAP) skipped : X_umap not available and X_pca absent.")

# ─────────────────────────────────────────────────────────────────────────────
# Figure 6 — Scatter: top 4 driver genes vs p_iPSC
# ─────────────────────────────────────────────────────────────────────────────
_genes_scatter = [g for g in df_a["gene"].head(4).tolist() if g in adata.var_names]
if _genes_scatter:
    _fig6, _ax6s = plt.subplots(
        1, len(_genes_scatter), figsize=(5 * len(_genes_scatter), 4))
    if len(_genes_scatter) == 1:
        _ax6s = [_ax6s]
    _p_vals = adata.obs["p_iPSC"].values.astype(float)
    for _ax6, _gene in zip(_ax6s, _genes_scatter):
        _gi   = list(adata.var_names).index(_gene)
        _expr = (adata.X[:, _gi].toarray().ravel()
                 if sp.issparse(adata.X)
                 else np.asarray(adata.X[:, _gi]).ravel())
        _mask6 = np.isfinite(_p_vals) & np.isfinite(_expr)
        _sc6   = _ax6.scatter(
            _expr[_mask6], _p_vals[_mask6],
            c=_p_vals[_mask6], cmap="viridis",
            s=2, alpha=0.3, rasterized=True,
        )
        plt.colorbar(_sc6, ax=_ax6, label="p_iPSC")
        _ax6.set_xlabel(f"{_gene} expression")
        _ax6.set_ylabel("p_iPSC")
        _ax6.set_title(_gene)
    _fig6.suptitle("Top driver gene expression vs p_iPSC fate probability",
                   fontsize=10, fontweight="bold")
    _fig6.tight_layout()
    _p6 = FIGURES_DIR / f"{TIMESTAMP}_driver_scatter_piPSC.png"
    _fig6.savefig(_p6, dpi=150, bbox_inches="tight")
    plt.close(_fig6)
    print(f"  Fig 6 (scatter expr vs p)    → {_p6.name}")

# =============================================================================
# 8.  FINAL SUMMARY
# =============================================================================
print(f"\n{_SEP}")
print("STEP 8  —  Final summary")
print(_SEP)

_elapsed = _time.time() - T0
print(f"  Total elapsed : {_elapsed:.1f} s")

print(f"\n  Output tables:")
print(f"    {_path_a.name}             (Analysis A)")
print(f"    {_path_b.name}           (Analysis B)")
print(f"    {_path_overlap.name}            (A∩B overlap, unfiltered)")
print(f"    {_path_overlap_filt.name}  (A∩B overlap, ribo-excluded)")
if _path_c is not None:
    print(f"    {_path_c.name}   (Analysis C — fate_ips vs early)")

print(f"\n  Analysis C positive group    : {_c_source if _mask_c_pos is not None else 'N/A (skipped)'}")
print(f"  Trajectory labels in figures : {_traj_label_source}")
print(f"  UMAP source                  : {'recomputed locally' if _umap_recomputed else 'from h5ad' if 'X_umap' in adata.obsm else 'unavailable'}")
print(f"  Ribosomal filter applied     : {_USE_FILTERED_FOR_PLOTS} ({_n_ribo_removed} genes removed)")

print(f"\n  Top 20 consensus driver genes (A∩B, {'ribo-excluded' if _USE_FILTERED_FOR_PLOTS else 'unfiltered'}):")
_top_df = _df_overlap_filtered if _USE_FILTERED_FOR_PLOTS else _df_overlap
for _i, (_gene_c, _row_c) in enumerate(
    _top_df.head(20).iterrows(), 1
):
    _g = _top_df["gene"].iloc[_i - 1]
    _row_a = df_a[df_a["gene"] == _g]
    _lfc_a = float(_row_a["logFC"].values[0]) if len(_row_a) else float("nan")
    print(f"    {_i:>2}. {_g:<14}  logFC(A)={_lfc_a:+.3f}")

print(f"\n  Next steps:")
print(f"    • Inspect figures in {FIGURES_DIR}")
print(f"    • Top driver genes are candidates for functional follow-up")
print(f"    • Step 04: run CellRank2 RealTimeKernel for WOT-extension comparison")
