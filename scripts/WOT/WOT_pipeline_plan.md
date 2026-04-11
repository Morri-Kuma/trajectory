# WOT Pipeline Plan for GSE230659

## Overview

WaddingtonOT (Schiebinger et al., Science 2019) models cell populations as probability
distributions evolving over time under optimal transport. It assumes **deterministic,
mass-conserving transport** between adjacent time points, optionally with growth/death
corrections.

This document outlines the specific steps to prepare GSE230659 data for WOT and run the
full trajectory inference pipeline.

---

## Prerequisites

- Completed: `scripts/C_traj/01_load_and_qc.py` → produces `*_GSE230659_qc.h5ad`
- Python packages: `wot`, `scanpy`, `anndata`, `numpy`, `pandas`, `matplotlib`
- Install wot: `pip install wot`

---

## Step-by-Step Pipeline

### Step 2A: Prepare WOT Input Files

**Script**: `scripts/WOT/02a_prepare_wot_inputs.py`

WOT requires specific input formats:

1. **Expression matrix** (genes × cells, or AnnData):
   - Use the QC-filtered, normalized, PCA-reduced AnnData from Step 01
   - WOT can operate on the full expression matrix or PCA embeddings
   - **Recommendation**: Use PCA (50 PCs) as the feature space for initial runs

2. **Cell days file** (`cell_days.txt`):
   - Two-column TSV: `id` (cell barcode) and `day` (absolute time)
   - Extracted from `adata.obs['abs_day']`
   - Format: `{cell_barcode}\t{absolute_day}`

3. **Cell growth rates file** (`cell_growth_rates.txt`) — *Critical for WOT*:
   - Two-column TSV: `id` and `cell_growth_rate`
   - Estimated from gene signatures of proliferation and apoptosis
   - See Step 2B below

### Step 2B: Estimate Cell Growth Rates

**Script**: `scripts/WOT/02b_estimate_growth_rates.py`

Growth rate estimation is the most biology-sensitive step in WOT. Options:

#### Option A: Gene-signature-based (Recommended for first pass)
- Score each cell for **proliferation** using S-phase + G2M-phase gene sets
  (from Tirosh et al. 2016 or Regev lab cell cycle genes)
- Score each cell for **apoptosis** using curated apoptosis gene signatures
  (e.g., GO:0006915 or Hallmark_Apoptosis from MSigDB)
- Compute: `growth_rate = birth_rate - death_rate`
  - `birth_rate ∝ proliferation_score`
  - `death_rate ∝ apoptosis_score`
- Normalize so that the **population-level growth matches observed cell count changes**
  between consecutive time points

#### Option B: Population-level only
- Use observed cell counts at each time point to estimate bulk growth factors
- Assign uniform growth rates to all cells within each time point
- Simpler but loses cell-level heterogeneity

#### Option C: WOT built-in (if no growth rate file provided)
- WOT defaults to uniform growth (mass conservation)
- This is the **null assumption** — good baseline, but biologically unrealistic
  for reprogramming where differential proliferation/death is expected

**Recommended approach**: Start with Option C (no growth rates) as baseline,
then implement Option A and compare trajectory results.

### Step 2C: Compute Temporal Couplings (Transport Maps)

**Script**: `scripts/WOT/02c_compute_transport.py`

```python
import wot

# Initialize OTModel
ot_model = wot.ot.OTModel(
    adata,                           # QC AnnData
    day_field='abs_day',             # column in .obs
    growth_rate_field='growth_rate', # column in .obs (or None for uniform)
    # Key hyperparameters:
    epsilon=0.05,                    # entropic regularization (start with 0.05)
    lambda1=1.0,                    # unbalanced penalty source
    lambda2=50.0,                   # unbalanced penalty target
)

# Compute transport maps between all consecutive time pairs
# This produces a TransportMapModel with T_{t→t+1} for each adjacent pair
tmap_model = ot_model.compute_all_transport_maps()
```

**Key hyperparameters to tune**:
- `epsilon` (entropic regularization): Controls transport map smoothness
  - Small ε → sparser, more deterministic couplings
  - Large ε → smoother, more diffuse couplings
  - Start: 0.05, grid search [0.01, 0.05, 0.1, 0.5]
- `lambda1`, `lambda2` (unbalanced OT penalties):
  - Control how much growth/death is allowed
  - Higher values → closer to balanced (mass-conserving) OT
  - Start: lambda1=1, lambda2=50
- `growth_iters`: Number of growth rate estimation iterations (default: 3)

### Step 2D: Trajectory Analysis

**Script**: `scripts/WOT/02d_trajectory_analysis.py`

With transport maps computed, perform:

1. **Ancestor/Descendant distributions**:
   ```python
   # Fate probabilities: what fraction of cells at time t end up in each fate?
   # Define cell sets at the final timepoint (hCiPSCs vs non-iPSC)
   cell_sets = wot.tmap.trajectory_trends(tmap_model, cell_sets_dict)
   ```

2. **Trajectory trends** (gene expression along trajectories):
   - Track how pluripotency markers (POU5F1/OCT4, SOX2, NANOG, KLF4) change
   - Track mesenchymal markers (VIM, CDH2) for EMT/MET transitions
   - Track apoptosis markers

3. **Fate probability at early timepoints**:
   - For each cell at early time points, compute P(reaching hCiPSC state)
   - This enables "early fate prediction" — identifying successful reprogramming cells

4. **Validation metrics**:
   - Temporal consistency of trajectories
   - Smoothness of gene expression trends
   - Interpolation at held-out time points (leave-one-out)

### Step 2E: Visualization & Metrics

**Script**: `scripts/WOT/02e_visualize.py`

- UMAP overlaid with fate probabilities
- Trajectory trend plots (gene expression × absolute time)
- Sankey/flow diagram of cell populations across time
- Heatmap of transport coupling matrices

---

## Time Point Considerations for WOT

The 15 timepoints in GSE230659 have **uneven spacing**:

| abs_day | stage_day_label      | gap to next |
|---------|---------------------|-------------|
| 0.50    | StageI_Day0.5       | 1.50        |
| 2.00    | StageI_Day2         | 2.00        |
| 4.00    | StageI_Day4         | 4.00        |
| 8.00    | StageI_Day8         | 4.00        |
| 12.00   | StageII_Day4        | 4.00        |
| 16.00   | StageII_Day8        | 0.33        |
| 16.33   | StageIII_Day0.33    | 0.34        |
| 16.67   | StageIII_Day0.67    | 0.33        |
| 17.00   | StageIII_Day1       | 1.00        |
| 18.00   | StageIII_Day2       | 2.00        |
| 20.00   | StageIII_Day4       | 2.00        |
| 22.00   | StageIII_Day6       | 2.00        |
| 24.00   | StageIII_Day8       | 4.00        |
| 28.00   | StageIII_Day12      | 2.00        |
| 30.00   | hCiPSCs             | endpoint    |

WOT handles uneven time gaps naturally through its temporal coupling formulation,
but large gaps (e.g., Day 4→8, Day 8→12) mean more uncertainty in the inferred
transport. The dense sampling around Stage III transition (Day 16-17) is valuable
for capturing early reprogramming dynamics.

---

## File Outputs (all timestamped)

| Step  | Output                                    | Location            |
|-------|-------------------------------------------|---------------------|
| 02a   | `cell_days.txt`, `cell_growth_rates.txt`  | `data/processed/`   |
| 02c   | `tmaps/tmap_{t1}_{t2}.h5ad`              | `results/tmaps/`    |
| 02d   | `fate_probabilities.csv`                  | `results/metrics/`  |
| 02d   | `trajectory_trends.csv`                   | `results/metrics/`  |
| 02e   | Various `.png` plots                      | `results/figures/`  |
