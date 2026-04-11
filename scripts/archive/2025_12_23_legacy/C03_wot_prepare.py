# C03_wot_prepare.py
# Build WOT inputs using discrete timepoints t_wot (0..K-1)
# Fixes:
#   1) Rename obs['day'] -> obs['day_in_stage'] in ExprMatrix to avoid WOT collision
#   2) Use non-numeric cell ids (prefix "cell_") to avoid pandas int/str join mismatch
#
# Outputs:
#   wot/input/timepoint_map.tsv
#   wot/input/cell_days.txt                 (id + day, where day == t_wot)
#   wot/input/hvg_genes.txt
#   wot/input/ExprMatrix.var.genes.h5ad
#   wot/input/wot_id_map.tsv                (wot_id <-> orig_id)

import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path

BASE_RESULTS = Path(r"E:\scgpt\results\C_traj")


def get_run_dir() -> Path:
    latest = BASE_RESULTS / "_latest_run.txt"
    if not latest.exists():
        raise FileNotFoundError(f"Missing {latest}. Run C01/C02 first.")
    run_dir = Path(latest.read_text(encoding="utf-8").strip())
    if not run_dir.exists():
        raise FileNotFoundError(f"Run dir not found: {run_dir}")
    return run_dir


def main():
    run_dir = get_run_dir()
    adir = run_dir / "adata"
    wotin = run_dir / "wot" / "input"
    wotin.mkdir(parents=True, exist_ok=True)

    in_h5ad = adir / "C2_pseudotime.h5ad"
    if not in_h5ad.exists():
        raise FileNotFoundError(f"Missing input h5ad: {in_h5ad}")

    adata = sc.read_h5ad(in_h5ad)

    # ------------------------------
    # 1) Build continuous time axis day_wot (R-style offsets)
    # ------------------------------
    if "stage" not in adata.obs or "day" not in adata.obs:
        raise ValueError("adata.obs must contain 'stage' and 'day' columns")

    adata.obs["stage"] = adata.obs["stage"].astype(str)
    adata.obs["day"] = adata.obs["day"].astype(float)

    # 1) 只对 StageI/II/III 计算 offset（更稳健：从数据里自动取 max）
    stI = adata.obs["stage"].eq("StageI")
    stII = adata.obs["stage"].eq("StageII")
    stIII = adata.obs["stage"].eq("StageIII")

    if not (stI.any() and stII.any() and stIII.any()):
        raise ValueError("Need StageI/StageII/StageIII to build continuous day_wot")

    max_I = float(np.nanmax(adata.obs.loc[stI, "day"].values))   # 期望 8.0
    max_II = float(np.nanmax(adata.obs.loc[stII, "day"].values)) # 期望 8.0

    stage_offset = {
        "StageI": 0.0,
        "StageII": max_I,
        "StageIII": max_I + max_II,
    }

    adata.obs["day_wot"] = np.nan
    mask_base = adata.obs["stage"].isin(stage_offset)
    adata.obs.loc[mask_base, "day_wot"] = (
        adata.obs.loc[mask_base, "stage"].map(stage_offset).astype(float).values
        + adata.obs.loc[mask_base, "day"].values
    )

    # 2) 终点 hCiPSC：统一放到最后（避免 B01 里 +9 的占位影响）
    mask_h = adata.obs["stage"].str.lower().str.contains("hcipsc", na=False)
    adata.obs.loc[mask_h, "stage"] = "hCiPSC"

    max_non = float(np.nanmax(adata.obs.loc[mask_base, "day_wot"].values))
    adata.obs.loc[mask_h, "day_wot"] = max_non + 1.0

    if adata.obs["day_wot"].isna().any():
        bad = adata.obs.loc[adata.obs["day_wot"].isna(), ["stage", "day"]].drop_duplicates()
        raise ValueError(f"Some cells have missing day_wot:\n{bad}")

    # 写 timepoint_map：保留 stage/day_in_stage 对应的 day_wot
    u = (
        adata.obs.groupby(["stage", "day"], observed=True)
        .agg(n_cells=("stage", "size"), day_wot=("day_wot", "median"))
        .reset_index()
        .sort_values(["day_wot", "stage", "day"])
        .reset_index(drop=True)
    )
    u.to_csv(wotin / "timepoint_map.tsv", sep="\t", index=False)

    # ------------------------------
    # 3) Write cell_days.txt (WOT expects columns: id, day)
    # Here day == day_wot (continuous)
    # ------------------------------
    cell_days = pd.DataFrame({"id": adata.obs_names.astype(str), "day": adata.obs["day_wot"].values})
    cell_days.to_csv(wotin / "cell_days.txt", sep="\t", index=False)

    # ------------------------------
    # 4) Write HVG genes list
    # ------------------------------
    adata.var_names.to_series().to_csv(
        wotin / "hvg_genes.txt", index=False, header=False
    )

    # ------------------------------
    # 5) Make sure X is log-normalized
    # ------------------------------
    if "log1p" not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    # ------------------------------
    # 6) Avoid column name collision with WOT day
    # ------------------------------
    if "day" in adata.obs.columns:
        adata.obs = adata.obs.rename(columns={"day": "day_in_stage"})

    # Save matrix for WOT
    adata.write_h5ad(wotin / "ExprMatrix.var.genes.h5ad")

    print("[OK] wrote WOT input to:", wotin)
    print("[OK] wrote:", wotin / "timepoint_map.tsv")
    print("[OK] wrote:", wotin / "cell_days.txt")
    print("[OK] wrote:", wotin / "hvg_genes.txt")
    print("[OK] wrote:", wotin / "ExprMatrix.var.genes.h5ad")
    print("[OK] wrote:", wotin / "wot_id_map.tsv")
    print("[NOTE] ExprMatrix.obs column 'day' renamed to 'day_in_stage' to avoid WOT collision.")
    print("[NOTE] Using WOT-safe ids with prefix 'cell_' to avoid int/str join issues.")


if __name__ == "__main__":
    main()
