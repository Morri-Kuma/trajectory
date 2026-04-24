"""Recompute remaining Forecast Accuracy metrics for scNODE reduced runs."""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import anndata as ad

from benchmark.evaluation.eval_forecast import run_forecast_evaluation


RUNS = [
    "benchmark/results/scnode/scenario_A_scgpt_v1_full_hvg2000_reduced",
    "benchmark/results/scnode/scenario_B_scgpt_v1_hvg2000_reduced",
    "benchmark/results/scnode/scenario_C_scgpt_v1_hvg2000_reduced",
]
HVG_ADATA = "benchmark/results/scnode/hvg_inputs/adata_scnode_full_A_hvg2000.h5ad"


def main() -> None:
    adata = ad.read_h5ad(ROOT / HVG_ADATA, backed="r")
    for run_dir_str in RUNS:
        run_dir = ROOT / run_dir_str
        print(f"recomputing forecast metrics: {run_dir}")
        metrics = run_forecast_evaluation(
            projected_expression_path=str(run_dir / "projected_expression.npy"),
            adata=adata,
            output_dir=str(run_dir),
            time_key="abs_day",
            max_cells_per_timepoint=1000,
        )
        print(metrics)


if __name__ == "__main__":
    main()
