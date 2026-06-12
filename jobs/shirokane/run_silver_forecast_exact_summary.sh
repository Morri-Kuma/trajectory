#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_silver_forecast_exact_summary.$JOB_ID.log
#$ -l s_vmem=16G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
RUN_TESTS="${RUN_TESTS:-1}"
RUN_OFFICIAL_SUMMARY="${RUN_OFFICIAL_SUMMARY:-1}"
REQUIRE_ALL_EXACT="${REQUIRE_ALL_EXACT:-1}"
LEIDEN_N_NEIGHBORS="${LEIDEN_N_NEIGHBORS:-15}"
LEIDEN_RESOLUTION="${LEIDEN_RESOLUTION:-0.5}"

cd "${PROJECT_ROOT}"
mkdir -p logs benchmark/reports/official_silver

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export RUN_TESTS
export RUN_OFFICIAL_SUMMARY
export REQUIRE_ALL_EXACT
export LEIDEN_N_NEIGHBORS
export LEIDEN_RESOLUTION
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Forecast exact summary started at : $(date)"
echo "Project root                      : ${PROJECT_ROOT}"
echo "Run tests                         : ${RUN_TESTS}"
echo "Run official summary              : ${RUN_OFFICIAL_SUMMARY}"
echo "Require all exact                 : ${REQUIRE_ALL_EXACT}"
echo "Leiden n_neighbors                : ${LEIDEN_N_NEIGHBORS}"
echo "Leiden resolution                 : ${LEIDEN_RESOLUTION}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

if [ "${RUN_TESTS}" = "1" ]; then
  python -m pytest benchmark/evaluation/tests/test_eval_forecast.py -q
fi

python - <<'PY'
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


ROOT = Path(os.environ["TRAJ_PROJECT_ROOT"]).resolve()
REQUIRE_ALL_EXACT = os.environ.get("REQUIRE_ALL_EXACT", "1") == "1"

CONFIGS = [
    "benchmark/configs/runtime/scnode_gse175634_cardiac_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse175634_cardiac_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse175634_cardiac_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse175634_cardiac_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse175634_cardiac_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse175634_cardiac_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse175634_cardiac_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse175634_cardiac_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse175634_cardiac_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse178325_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse178325_marker_fm_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse178325_marker_fm_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse178325_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse178325_marker_fm_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse178325_marker_fm_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse178325_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse178325_marker_fm_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse178325_marker_fm_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse230659_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse230659_marker_fm_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse230659_marker_fm_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/scnode_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml",
    "benchmark/configs/scnode_gse242424_oskm_ground_truth_B_hvg2000_formal.yaml",
    "benchmark/configs/scnode_gse242424_oskm_ground_truth_C_hvg2000_formal.yaml",
    "benchmark/configs/mioflow_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml",
    "benchmark/configs/mioflow_gse242424_oskm_ground_truth_B_hvg2000_formal.yaml",
    "benchmark/configs/mioflow_gse242424_oskm_ground_truth_C_hvg2000_formal.yaml",
    "benchmark/configs/prescient_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml",
    "benchmark/configs/prescient_gse242424_oskm_ground_truth_B_hvg2000_formal.yaml",
    "benchmark/configs/prescient_gse242424_oskm_ground_truth_C_hvg2000_formal.yaml",
]


def read_yaml(path: Path) -> dict:
    with open(path, encoding="utf-8-sig") as f:
        return yaml.safe_load(f) or {}


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        return json.load(f) or {}


def resolve_path(value):
    if value in (None, ""):
        return None
    path = Path(str(value))
    return path if path.is_absolute() else ROOT / path


def is_finite_number(value) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except Exception:
        return False


rows = []
errors = []
for rel_config in CONFIGS:
    config_path = ROOT / rel_config
    cfg = read_yaml(config_path)
    output_dir = resolve_path((cfg.get("output") or {}).get("base_dir"))
    metrics_path = output_dir / "forecast_metrics.json" if output_dir else None
    projected_expression = output_dir / "projected_expression.npy" if output_dir else None
    fm = read_json(metrics_path) if metrics_path else {}
    metric_keys = [
        "wasserstein_distance",
        "gaussian_mmd",
        "energy_distance_mmd",
        "hausdorff_loss",
    ]
    exact = (
        fm.get("metric_backend") == "scTimeBench_exact"
        and fm.get("metric_protocol") == "sctimebench_gex_prediction_otloss"
        and fm.get("exact_cell_usage") is True
        and fm.get("lognorm") is False
        and fm.get("aggregate") == "mean"
        and fm.get("status") == "completed"
        and all(is_finite_number(fm.get(key)) for key in metric_keys)
    )
    row = {
        "config": rel_config,
        "dataset_id": (cfg.get("dataset") or {}).get("id"),
        "method": cfg.get("method"),
        "scenario": cfg.get("scenario"),
        "output_dir": str(output_dir.relative_to(ROOT)) if output_dir else None,
        "projected_expression_exists": projected_expression.exists() if projected_expression else False,
        "forecast_metrics_exists": metrics_path.exists() if metrics_path else False,
        "metric_backend": fm.get("metric_backend"),
        "metric_protocol": fm.get("metric_protocol"),
        "exact_cell_usage": fm.get("exact_cell_usage"),
        "lognorm": fm.get("lognorm"),
        "aggregate": fm.get("aggregate"),
        "status": fm.get("status"),
        "n_eval_timepoints": fm.get("n_eval_timepoints"),
        "wasserstein_distance": fm.get("wasserstein_distance"),
        "gaussian_mmd": fm.get("gaussian_mmd"),
        "energy_distance_mmd": fm.get("energy_distance_mmd"),
        "hausdorff_loss": fm.get("hausdorff_loss"),
        "official_forecast_exact": exact,
    }
    rows.append(row)
    if not exact:
        errors.append(f"{rel_config}: forecast exact validation failed")

inventory = pd.DataFrame(rows)
report_dir = ROOT / "benchmark" / "reports" / "official_silver"
report_dir.mkdir(parents=True, exist_ok=True)
inventory_path = report_dir / "official_silver_forecast_exact_rerun_inventory.csv"
inventory.to_csv(inventory_path, index=False)

print(json.dumps({
    "inventory": str(inventory_path.relative_to(ROOT)),
    "n_configs": len(CONFIGS),
    "n_exact": int(inventory["official_forecast_exact"].sum()),
    "n_failed": len(errors),
}, indent=2))

if errors:
    print("Forecast exact validation failures:")
    for error in errors:
        print(f"- {error}")
    if REQUIRE_ALL_EXACT:
        raise SystemExit("One or more forecast exact reruns are missing or invalid.")
PY

if [ "${RUN_OFFICIAL_SUMMARY}" = "1" ]; then
  python -m benchmark.evaluation.summarize_official_silver \
    --leiden-n-neighbors "${LEIDEN_N_NEIGHBORS}" \
    --leiden-resolution "${LEIDEN_RESOLUTION}"
fi

echo "============================================================"
echo "Forecast exact summary finished at: $(date)"
echo "Inventory: benchmark/reports/official_silver/official_silver_forecast_exact_rerun_inventory.csv"
echo "============================================================"
