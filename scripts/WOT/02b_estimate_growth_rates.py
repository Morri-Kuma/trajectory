#!/usr/bin/env python3
"""
02b_estimate_growth_rates.py — Estimate per-cell growth rates for WOT.

Project : Comparative Study of Trajectory Inference Models for Chemical iPSC Reprogramming
Dataset : GSE230659 (human, Liuyang et al. 2023 Cell Stem Cell)
Step    : 2b — Growth rate estimation via cell cycle + apoptosis gene scoring
Author  : Kuma
Env     : conda activate traj_env

Background:
    WOT models cell populations as distributions evolving under optimal transport.
    By default it assumes mass conservation (no growth/death). For iPSC reprogramming,
    where cells proliferate at different rates and some undergo apoptosis, providing
    per-cell growth rates significantly improves trajectory accuracy.

    Growth rate = birth_rate - death_rate
    We estimate:
        birth_rate  ∝ cell cycle score (S-phase + G2M gene sets, Tirosh et al. 2016)
        death_rate  ∝ apoptosis score (Hallmark_Apoptosis gene set, MSigDB)

    Rates are then rescaled so that the POPULATION-LEVEL growth matches the
    observed fold-change in cell counts between consecutive timepoints.

Prerequisites:
    Step 02a output: *_GSE230659_wot_ready.h5ad

Usage:
    conda activate traj_env
    cd C:\\Users\\37620\\trajectory
    python scripts/WOT/02b_estimate_growth_rates.py

Outputs (all timestamped):
    data/processed/  YYYYMMDD_HHMM_cell_growth_rates.txt   — WOT growth rates file
    data/processed/  YYYYMMDD_HHMM_GSE230659_wot_gr.h5ad  — AnnData with growth rates
    results/figures/ YYYYMMDD_HHMM_growth_rates.png        — Growth rate distributions
    results/figures/ YYYYMMDD_HHMM_cell_cycle_scores.png   — Cell cycle score plots
"""

# =============================================================================
# 0. CONFIGURATION
# =============================================================================
from datetime import datetime
from pathlib import Path
import gc

TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")

_candidates = [
    Path(r"C:\Users\37620\trajectory"),
    Path("/sessions/cool-admiring-hawking/mnt/trajectory"),
    Path(__file__).resolve().parents[2],
]
PROJECT_ROOT = next((p for p in _candidates if p.exists()), _candidates[-1])
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
FIGURES_DIR   = PROJECT_ROOT / "results" / "figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

def find_latest(pattern: str) -> Path:
    matches = sorted(PROCESSED_DIR.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No files matching '{pattern}' in {PROCESSED_DIR}.\n"
            f"Run scripts/WOT/02a_prepare_wot_inputs.py first."
        )
    return matches[-1]

H5AD_PATH = find_latest("*_GSE230659_wot_ready.h5ad")

# =============================================================================
# 1. GENE SIGNATURES
# =============================================================================
# Tirosh et al. 2016 (Science) cell cycle gene sets — widely used in scRNA-seq
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG",
    "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP",
    "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51",
    "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM",
    "CASP8AP2", "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8",
]

G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
    "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF",
    "TACC3", "FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB",
    "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP",
    "CDCA3", "HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1",
    "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
    "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF",
    "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA",
]

# MSigDB Hallmark_Apoptosis (representative subset of high-confidence genes)
APOPTOSIS_GENES = [
    "CASP3", "CASP6", "CASP7", "CASP8", "CASP9", "CASP10",
    "BAX", "BAK1", "BCL2", "BCL2L1", "MCL1", "BBC3", "PMAIP1",
    "BID", "DIABLO", "APAF1", "CYCS",
    "TP53", "MDM2", "CDKN1A",
    "FAS", "FASLG", "TNF", "TNFRSF10A", "TNFRSF10B",
    "DFFA", "DFFB",
    "XIAP", "BIRC2", "BIRC3",
    "PARP1",
    "AKT1", "AKT2", "PIK3CA",
    "LMNA", "VIM",
]

# =============================================================================
# 2. IMPORTS
# =============================================================================
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

print(f"[{TIMESTAMP}] WOT Growth Rate Estimation")
print(f"  Input AnnData: {H5AD_PATH.name}")

# =============================================================================
# 3. LOAD ANNDATA
# =============================================================================
print("\n--- Loading AnnData ---")
adata = sc.read_h5ad(H5AD_PATH)
print(f"  Shape: {adata.n_obs:,} x {adata.n_vars:,}")
print(f"  Genes in HVG set: {adata.n_vars}")

# Store log-normalized values for scoring
# adata.X is already log-normalized from Step 01
adata.layers["log_norm"] = adata.X.copy()

# =============================================================================
# 4. CELL CYCLE SCORING
# =============================================================================
# scanpy.tl.score_genes computes the mean expression of a gene set,
# controlled for expression level using a reference set of equal-expression genes.
print("\n--- Cell cycle scoring (Tirosh 2016 gene sets) ---")

# Filter gene lists to genes present in our HVG set
s_genes_present   = [g for g in S_GENES   if g in adata.var_names]
g2m_genes_present = [g for g in G2M_GENES if g in adata.var_names]

print(f"  S-phase genes   : {len(s_genes_present)}/{len(S_GENES)} present in HVG set")
print(f"  G2M genes       : {len(g2m_genes_present)}/{len(G2M_GENES)} present in HVG set")

# Score S phase
sc.tl.score_genes(adata, s_genes_present, score_name="s_score", use_raw=False)
# Score G2M phase
sc.tl.score_genes(adata, g2m_genes_present, score_name="g2m_score", use_raw=False)

# Assign cell cycle phase based on highest score
adata.obs["cc_phase"] = "G1"
s_higher  = adata.obs["s_score"] > adata.obs["g2m_score"]
g2m_above = adata.obs["g2m_score"] > 0
s_above   = adata.obs["s_score"] > 0
adata.obs.loc[s_higher  & s_above,   "cc_phase"] = "S"
adata.obs.loc[~s_higher & g2m_above, "cc_phase"] = "G2M"

phase_counts = adata.obs["cc_phase"].value_counts()
print(f"\n  Phase distribution:")
for phase, count in phase_counts.items():
    print(f"    {phase}: {count:>7,}  ({count/len(adata)*100:.1f}%)")

# =============================================================================
# 5. APOPTOSIS SCORING
# =============================================================================
print("\n--- Apoptosis scoring ---")
apop_genes_present = [g for g in APOPTOSIS_GENES if g in adata.var_names]
print(f"  Apoptosis genes : {len(apop_genes_present)}/{len(APOPTOSIS_GENES)} present in HVG set")

if len(apop_genes_present) >= 5:
    sc.tl.score_genes(adata, apop_genes_present, score_name="apoptosis_score",
                      use_raw=False)
    print(f"  Apoptosis score range: "
          f"[{adata.obs['apoptosis_score'].min():.3f}, "
          f"{adata.obs['apoptosis_score'].max():.3f}]")
else:
    # Fallback: use mito% as proxy for apoptosis (stressed/dying cells)
    print(f"  WARNING: <5 apoptosis genes in HVG set. Using mito% as proxy.")
    adata.obs["apoptosis_score"] = (
        adata.obs["pct_counts_mt"] / adata.obs["pct_counts_mt"].max()
    ).astype(np.float32)

# =============================================================================
# 6. COMPUTE PROLIFERATION + GROWTH RATES
# =============================================================================
print("\n--- Computing growth rates ---")

# Proliferation rate: combine S and G2M scores
# Shift scores to be non-negative, then normalize to [0, 1]
s_score   = adata.obs["s_score"].values.copy()
g2m_score = adata.obs["g2m_score"].values.copy()
apop      = adata.obs["apoptosis_score"].values.copy()

# Proliferation = max(s_score, g2m_score) since they represent sequential phases
prolif = np.maximum(s_score, g2m_score)

# Clip to non-negative (negative scores = below background)
prolif = np.clip(prolif, 0, None)
apop   = np.clip(apop,   0, None)

# Normalize each to [0, 1] range
def norm01(x):
    rng = x.max() - x.min()
    return (x - x.min()) / rng if rng > 0 else np.zeros_like(x)

prolif_norm = norm01(prolif)
apop_norm   = norm01(apop)

# Raw growth rate = proliferation - apoptosis
# Use weights: proliferation dominates (iPSC reprogramming is primarily a survival
# + proliferation problem; apoptosis is a secondary selection mechanism)
BIRTH_WEIGHT = 0.7
DEATH_WEIGHT = 0.3
raw_growth = BIRTH_WEIGHT * prolif_norm - DEATH_WEIGHT * apop_norm

# =============================================================================
# 7. POPULATION-LEVEL RESCALING
# =============================================================================
# WOT uses growth rates to reweight transport maps.
# Rescale so that the sum of growth rates per timepoint reflects the OBSERVED
# fold-change in cell numbers between consecutive timepoints.
#
# For time point t with n_t cells and next time point t+1 with n_{t+1} cells:
#   Expected total growth factor G_t = n_{t+1} / n_t
#   We set: growth_rate[cell in t] such that mean(exp(growth_rate)) = G_t
#
# WOT internally uses: exp(growth_rate * dt) as the growth factor
# So: growth_rate = log(G_t) / dt

print("  Population-level rescaling (observed cell count fold-changes):")

# Cell counts per timepoint
time_counts = (
    adata.obs.groupby("abs_day", observed=True)
    .size()
    .sort_index()
)
sorted_days = time_counts.index.values
sorted_n    = time_counts.values

# Compute observed fold-changes
fold_changes = {}
for i in range(len(sorted_days) - 1):
    t1, t2 = sorted_days[i], sorted_days[i+1]
    n1, n2 = sorted_n[i], sorted_n[i+1]
    dt = t2 - t1
    fc = n2 / n1
    # Expected log growth rate
    expected_log_gr = np.log(fc) / dt
    fold_changes[t1] = {"fc": fc, "dt": dt, "expected_log_gr": expected_log_gr}
    print(f"    t={t1:.2f}→{t2:.2f}: n={n1:>6,}→{n2:>6,}  "
          f"fc={fc:.3f}  log_gr={expected_log_gr:+.4f}")

# Assign rescaled growth rates per timepoint
growth_rates = np.zeros(adata.n_obs, dtype=np.float32)

for t_day, stats in fold_changes.items():
    mask = adata.obs["abs_day"].values == t_day
    raw_t = raw_growth[mask]

    # Center raw scores around the expected population mean
    expected = stats["expected_log_gr"]
    # Shift so that mean of growth_rates for this timepoint = expected
    current_mean = raw_t.mean()
    shift = expected - current_mean
    growth_rates[mask] = (raw_t + shift).astype(np.float32)

# Last timepoint = terminal reference pool (all hCiPSC/day30 cells).
# Growth rate is set to 0: these cells are the absorbing endpoint of the model,
# not an intermediate state undergoing further dynamics.
# NOTE: growth_rate=0 is assigned to ALL hCiPSC cells (the full day30 pool).
# Within this pool, the high-confidence iPSC terminal subset is later defined
# by OCLR stemness score (top 10%) in 02b2, then loaded into fate_ips/fate_other
# columns by 02c.  The growth-rate assignment here is independent of that split.
last_day = sorted_days[-1]
mask_last = adata.obs["abs_day"].values == last_day
growth_rates[mask_last] = 0.0

adata.obs["growth_rate"] = growth_rates

print(f"\n  Growth rate summary:")
for day in sorted_days:
    mask = adata.obs["abs_day"].values == day
    gr = growth_rates[mask]
    print(f"    Day {day:5.2f}: mean={gr.mean():+.4f}  "
          f"std={gr.std():.4f}  "
          f"range=[{gr.min():+.4f}, {gr.max():+.4f}]")

# =============================================================================
# 8. WRITE cell_growth_rates.txt FOR WOT
# =============================================================================
print("\n--- Writing cell_growth_rates.txt ---")

gr_df = adata.obs[["growth_rate"]].copy()
gr_df.index.name = "id"

gr_path = PROCESSED_DIR / f"{TIMESTAMP}_cell_growth_rates.txt"
gr_df.to_csv(gr_path, sep="\t")
print(f"  Saved → {gr_path}")

# =============================================================================
# 9. SAVE UPDATED ANNDATA
# =============================================================================
print("\n--- Saving updated AnnData ---")

gr_h5ad_path = PROCESSED_DIR / f"{TIMESTAMP}_GSE230659_wot_gr.h5ad"
adata.write_h5ad(gr_h5ad_path, compression="gzip")
print(f"  Saved → {gr_h5ad_path}")

# =============================================================================
# 10. VISUALIZATIONS
# =============================================================================
print("\n--- Generating plots ---")

STAGE_COLORS = {
    "StageI": "#2196F3", "StageII": "#FF9800",
    "StageIII": "#4CAF50", "hCiPSCs": "#E91E63",
}
stage_order = ["StageI", "StageII", "StageIII", "hCiPSCs"]

# --- Cell cycle scores by stage ---
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for ax, score, title in zip(
    axes,
    ["s_score", "g2m_score", "apoptosis_score"],
    ["S-phase Score", "G2M Score", "Apoptosis Score"],
):
    sns.violinplot(data=adata.obs, x="stage", y=score, order=stage_order,
                   palette=STAGE_COLORS, ax=ax, inner="box", cut=0, linewidth=0.5)
    ax.set_title(title)
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=30)

fig.suptitle("Cell Cycle & Apoptosis Scores by Stage", fontsize=13, y=1.01)
fig.tight_layout()
cc_path = FIGURES_DIR / f"{TIMESTAMP}_cell_cycle_scores.png"
fig.savefig(cc_path, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {cc_path}")

# --- Cell cycle phase composition per timepoint ---
phase_comp = (
    adata.obs.groupby(["stage_day_label", "cc_phase"], observed=True)
    .size()
    .unstack(fill_value=0)
)
# Sort by abs_day
sort_order = (
    adata.obs.drop_duplicates("stage_day_label")
    .sort_values("abs_day")["stage_day_label"]
    .values
)
phase_comp = phase_comp.reindex(sort_order)
phase_pct = phase_comp.div(phase_comp.sum(axis=1), axis=0) * 100

fig, ax = plt.subplots(figsize=(15, 5))
bottom = np.zeros(len(phase_pct))
phase_palette = {"G1": "#90CAF9", "S": "#FF8F00", "G2M": "#EF5350"}
for phase in ["G1", "S", "G2M"]:
    if phase in phase_pct.columns:
        ax.bar(range(len(phase_pct)), phase_pct[phase], bottom=bottom,
               label=phase, color=phase_palette[phase])
        bottom += phase_pct[phase].values

ax.set_xticks(range(len(phase_pct)))
ax.set_xticklabels(phase_pct.index, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Cell fraction (%)")
ax.set_title("Cell cycle phase composition per timepoint")
ax.legend(loc="upper right")
fig.tight_layout()
phase_path = FIGURES_DIR / f"{TIMESTAMP}_cell_cycle_phases.png"
fig.savefig(phase_path, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {phase_path}")

# --- Growth rates by timepoint (box + scatter) ---
abs_days_sorted = sorted(adata.obs["abs_day"].unique())
labels = [
    adata.obs.loc[adata.obs["abs_day"] == d, "stage_day_label"].iloc[0]
    for d in abs_days_sorted
]

fig, ax = plt.subplots(figsize=(15, 5))
gr_by_day = [
    adata.obs.loc[adata.obs["abs_day"] == d, "growth_rate"].values
    for d in abs_days_sorted
]
bp = ax.boxplot(gr_by_day, patch_artist=True, medianprops={"color": "black", "lw": 2},
                flierprops={"marker": ".", "markersize": 2, "alpha": 0.3},
                whiskerprops={"linewidth": 0.8})

# Color boxes by stage
for i, d in enumerate(abs_days_sorted):
    stg = adata.obs.loc[adata.obs["abs_day"] == d, "stage"].iloc[0]
    bp["boxes"][i].set_facecolor(STAGE_COLORS.get(stg, "#999"))
    bp["boxes"][i].set_alpha(0.8)

ax.axhline(0, color="black", linestyle="--", alpha=0.4, linewidth=1)
ax.set_xticks(range(1, len(abs_days_sorted) + 1))
ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
ax.set_ylabel("Growth rate (log scale)")
ax.set_title("Per-cell growth rates by timepoint (rescaled to observed population fold-changes)")
fig.tight_layout()
gr_path_fig = FIGURES_DIR / f"{TIMESTAMP}_growth_rates.png"
fig.savefig(gr_path_fig, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"  Saved → {gr_path_fig}")

# =============================================================================
# 11. SUMMARY
# =============================================================================
print(f"\n{'='*60}")
print(f"[{TIMESTAMP}] Step 02b COMPLETE")
print(f"\n  Growth rate files:")
print(f"    {gr_path.name}")
print(f"    {gr_h5ad_path.name}")
print(f"\n  Cell cycle scoring:")
print(f"    S genes used    : {len(s_genes_present)}")
print(f"    G2M genes used  : {len(g2m_genes_present)}")
print(f"    Apoptosis genes : {len(apop_genes_present)}")
print(f"\n  Phase distribution:")
for phase, n in phase_counts.items():
    print(f"    {phase}: {n:,}  ({n/len(adata)*100:.1f}%)")
print(f"\n  Next: run 02c_compute_transport.py")
print(f"{'='*60}")
