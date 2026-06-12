#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_silver_forecast_exact.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=96G
#$ -t 1-36

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"
RUN_TAG="${RUN_TAG:-forecast_exact_$(date +%Y%m%d_%H%M%S)}"
BACKUP_OLD="${BACKUP_OLD:-1}"

CONFIGS=(
  "benchmark/configs/runtime/scnode_gse175634_cardiac_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse175634_cardiac_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse175634_cardiac_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse175634_cardiac_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse175634_cardiac_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse175634_cardiac_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse175634_cardiac_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse175634_cardiac_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse175634_cardiac_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse178325_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse178325_marker_fm_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse178325_marker_fm_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse178325_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse178325_marker_fm_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse178325_marker_fm_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse178325_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse178325_marker_fm_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse178325_marker_fm_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse230659_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse230659_marker_fm_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse230659_marker_fm_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/scnode_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml"
  "benchmark/configs/scnode_gse242424_oskm_ground_truth_B_hvg2000_formal.yaml"
  "benchmark/configs/scnode_gse242424_oskm_ground_truth_C_hvg2000_formal.yaml"
  "benchmark/configs/mioflow_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml"
  "benchmark/configs/mioflow_gse242424_oskm_ground_truth_B_hvg2000_formal.yaml"
  "benchmark/configs/mioflow_gse242424_oskm_ground_truth_C_hvg2000_formal.yaml"
  "benchmark/configs/prescient_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml"
  "benchmark/configs/prescient_gse242424_oskm_ground_truth_B_hvg2000_formal.yaml"
  "benchmark/configs/prescient_gse242424_oskm_ground_truth_C_hvg2000_formal.yaml"
)

if [ "${TASK_ID}" -lt 1 ] || [ "${TASK_ID}" -gt "${#CONFIGS[@]}" ]; then
  echo "ERROR: unsupported task id ${TASK_ID}; expected 1-${#CONFIGS[@]}" >&2
  exit 2
fi

CONFIG="${CONFIGS[$((TASK_ID - 1))]}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export CONFIG
export RUN_TAG
export BACKUP_OLD
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Forecast exact rerun started at : $(date)"
echo "Task ID                         : ${TASK_ID}/${#CONFIGS[@]}"
echo "Config                          : ${CONFIG}"
echo "Project root                    : ${PROJECT_ROOT}"
echo "Run tag                         : ${RUN_TAG}"
echo "Backup old forecast outputs     : ${BACKUP_OLD}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

python - <<'PY'
import json
import os
import shutil
from pathlib import Path

import anndata as ad
import numpy as np
import yaml

from benchmark.evaluation.eval_forecast import run_forecast_evaluation


ROOT = Path(os.environ["TRAJ_PROJECT_ROOT"]).resolve()
CONFIG_PATH = Path(os.environ["CONFIG"])
if not CONFIG_PATH.is_absolute():
    CONFIG_PATH = ROOT / CONFIG_PATH
RUN_TAG = os.environ.get("RUN_TAG", "forecast_exact")
BACKUP_OLD = os.environ.get("BACKUP_OLD", "1") == "1"


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


def backup_file(path: Path, suffix: str) -> None:
    if not BACKUP_OLD or not path.exists():
        return
    backup = path.with_name(f"{path.stem}.{suffix}_{RUN_TAG}{path.suffix}")
    if not backup.exists():
        shutil.copy2(path, backup)
        print(f"[forecast_exact] Backed up {path} -> {backup}")


cfg = read_yaml(CONFIG_PATH)
method = str(cfg.get("method", "")).lower()
if method not in {"scnode", "mioflow", "prescient"}:
    raise SystemExit(f"Forecast exact rerun only supports projection methods, got {method!r}")

dataset = cfg.get("dataset") or {}
scenario_params = cfg.get("scenario_params") or {}
output_dir = resolve_path((cfg.get("output") or {}).get("base_dir"))
h5ad_path = resolve_path(
    dataset.get("h5ad_path")
    or (cfg.get("state_system") or {}).get("source_h5ad")
)
time_key = dataset.get("time_key") or cfg.get("time_key") or "abs_day"

if output_dir is None:
    raise SystemExit(f"output.base_dir missing in {CONFIG_PATH}")
projected_expression_path = output_dir / "projected_expression.npy"
forecast_metrics_path = output_dir / "forecast_metrics.json"
per_tp_path = output_dir / "per_timepoint_forecast_metrics.csv"

for required in [CONFIG_PATH, h5ad_path, output_dir, projected_expression_path]:
    if required is None or not Path(required).exists():
        raise SystemExit(f"required path missing: {required}")

existing_metrics = read_json(forecast_metrics_path)
eval_timepoints = scenario_params.get("heldout_times") or existing_metrics.get("eval_timepoints")

print(json.dumps({
    "config": str(CONFIG_PATH.relative_to(ROOT)),
    "method": method,
    "scenario": cfg.get("scenario"),
    "dataset": dataset.get("id"),
    "output_dir": str(output_dir.relative_to(ROOT)),
    "projected_expression": str(projected_expression_path.relative_to(ROOT)),
    "h5ad_path": str(h5ad_path.relative_to(ROOT)),
    "time_key": time_key,
    "eval_timepoints_source": (
        "scenario_params.heldout_times"
        if scenario_params.get("heldout_times")
        else "existing forecast_metrics.json"
        if existing_metrics.get("eval_timepoints")
        else "adata.obs unique timepoints"
    ),
}, indent=2))

adata = ad.read_h5ad(str(h5ad_path))
if eval_timepoints is None:
    if time_key not in adata.obs.columns:
        raise KeyError(f"time_key {time_key!r} not found in {h5ad_path}")
    eval_timepoints = sorted(
        float(t) for t in np.unique(adata.obs[time_key].astype(float).to_numpy())
    )
eval_timepoints = [float(t) for t in eval_timepoints]

proj = np.load(projected_expression_path, mmap_mode="r")
if proj.ndim == 3 and proj.shape[0] != len(eval_timepoints):
    raise SystemExit(
        f"projected_expression has {proj.shape[0]} timepoint blocks but "
        f"{len(eval_timepoints)} eval_timepoints were resolved: {eval_timepoints}"
    )
if proj.ndim == 2 and proj.shape[0] % len(eval_timepoints) != 0:
    raise SystemExit(
        f"2D projected_expression rows={proj.shape[0]} is not divisible by "
        f"n_eval_timepoints={len(eval_timepoints)}"
    )
del proj

backup_file(forecast_metrics_path, "pre_exact")
backup_file(per_tp_path, "pre_exact")

metrics = run_forecast_evaluation(
    projected_expression_path=str(projected_expression_path),
    adata=adata,
    output_dir=str(output_dir),
    eval_timepoints=eval_timepoints,
    time_key=str(time_key),
    max_cells_per_timepoint=None,
    lognorm=False,
    aggregate="mean",
)

required_meta = {
    "metric_backend": "scTimeBench_exact",
    "metric_protocol": "sctimebench_gex_prediction_otloss",
    "exact_cell_usage": True,
    "lognorm": False,
    "aggregate": "mean",
    "status": "completed",
}
for key, expected in required_meta.items():
    if metrics.get(key) != expected:
        raise SystemExit(f"unexpected {key}: {metrics.get(key)!r}, expected {expected!r}")
for key in ("wasserstein_distance", "gaussian_mmd", "energy_distance_mmd", "hausdorff_loss"):
    if metrics.get(key) is None:
        raise SystemExit(f"missing forecast metric {key}")

print(json.dumps({
    "status": metrics.get("status"),
    "metric_backend": metrics.get("metric_backend"),
    "metric_protocol": metrics.get("metric_protocol"),
    "exact_cell_usage": metrics.get("exact_cell_usage"),
    "n_eval_timepoints": metrics.get("n_eval_timepoints"),
    "wasserstein_distance": metrics.get("wasserstein_distance"),
    "gaussian_mmd": metrics.get("gaussian_mmd"),
    "energy_distance_mmd": metrics.get("energy_distance_mmd"),
    "hausdorff_loss": metrics.get("hausdorff_loss"),
}, indent=2))
PY

echo "============================================================"
echo "Forecast exact rerun finished at: $(date)"
echo "Config: ${CONFIG}"
echo "============================================================"
