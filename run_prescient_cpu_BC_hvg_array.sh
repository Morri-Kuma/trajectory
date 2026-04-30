#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_prescient_cpu_BC_hvg.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=48G
#$ -t 1-2

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
PRESCIENT_REPO="${PRESCIENT_REPO:-${PROJECT_ROOT}/benchmark/methods/PRESCIENT/prescient_module}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
ADATA="benchmark/inputs/gse230659_scgpt_hvg2000/GSE230659_scGPT_annotated_HVG2000_benchmark_input.h5ad"

case "${SGE_TASK_ID:-${TASK_ID:-1}}" in
  1)
    SCENARIO="B"
    BASE_CONFIG="benchmark/configs/prescient_gse230659_observed_scgpt_v1_scenarioB_hvg2000_formal.yaml"
    CPU_CONFIG="benchmark/configs/prescient_gse230659_observed_scgpt_v1_scenarioB_hvg2000_cpu.yaml"
    OUTPUT_DIR="benchmark/results/prescient/scenario_B_scgpt_v1_hvg2000_cpu"
    ;;
  2)
    SCENARIO="C"
    BASE_CONFIG="benchmark/configs/prescient_gse230659_observed_scgpt_v1_scenarioC_hvg2000_formal.yaml"
    CPU_CONFIG="benchmark/configs/prescient_gse230659_observed_scgpt_v1_scenarioC_hvg2000_cpu.yaml"
    OUTPUT_DIR="benchmark/results/prescient/scenario_C_scgpt_v1_hvg2000_cpu"
    ;;
  *)
    echo "Unsupported SGE_TASK_ID=${SGE_TASK_ID:-${TASK_ID:-unset}}; expected 1 or 2" >&2
    exit 2
    ;;
esac

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PRESCIENT_REPO="${PRESCIENT_REPO}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Job started at: $(date)"
echo "Host          : $(hostname)"
echo "Project root  : ${TRAJ_PROJECT_ROOT}"
echo "PRESCIENT repo: ${PRESCIENT_REPO}"
echo "Scenario      : ${SCENARIO} HVG2000 CPU low-memory"
echo "Base config   : ${BASE_CONFIG}"
echo "CPU config    : ${CPU_CONFIG}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Task ID       : ${SGE_TASK_ID:-${TASK_ID:-N/A}}"
echo "============================================================"

if [ ! -f "${PRESCIENT_REPO}/prescient/train/model.py" ]; then
  echo "ERROR: PRESCIENT_REPO does not contain prescient/train/model.py: ${PRESCIENT_REPO}" >&2
  echo "Sync benchmark/methods/PRESCIENT/prescient_module or set PRESCIENT_REPO to the original PRESCIENT checkout." >&2
  exit 2
fi

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

python - <<PY
from pathlib import Path
import yaml

base = Path("${BASE_CONFIG}")
out = Path("${CPU_CONFIG}")
with open(base, encoding="utf-8-sig") as f:
    cfg = yaml.safe_load(f)
cfg["run_id"] = f"prescient_gse230659_${SCENARIO}_scgpt_v1_hvg2000_cpu"
cfg["result_class"] = "reduced_validation"
cfg["formal_benchmark"] = False
cfg["scenario_name"] = f"Scenario ${SCENARIO} PRESCIENT CPU low-memory run (HVG2000, scGPT-v1)"
cfg["output"]["base_dir"] = "${OUTPUT_DIR}"
params = cfg.setdefault("prescient_params", {})
params.update({
    "use_cuda": False,
    "num_pcs": 20,
    "k_dim": 64,
    "pretrain_epochs": 10,
    "train_epochs": 50,
    "train_batch": 0.02,
    "n_sim_cells": 500,
    "metric_sample_cells": 300,
    "subsample_per_timepoint": 500,
    "embedding_reference_cells_per_timepoint": 500,
})
params["memory_note"] = (
    "CPU low-memory run: capped at 500 training cells per timepoint, "
    "20 PCs, k_dim=64, n_sim_cells=500."
)
out.parent.mkdir(parents=True, exist_ok=True)
with open(out, "w", encoding="utf-8") as f:
    yaml.safe_dump(cfg, f, sort_keys=False)
print(f"wrote {out}")
print("train_times:", cfg.get("scenario_params", {}).get("train_times"))
print("heldout_times:", cfg.get("scenario_params", {}).get("heldout_times"))
PY

python - <<'PY'
import torch
print("torch:", torch.__version__)
print("cuda_available:", torch.cuda.is_available())
PY

python benchmark/evaluation/eval_dispatch.py \
  --method prescient \
  --scenario "${SCENARIO}" \
  --adata "${ADATA}" \
  --method-config "${CPU_CONFIG}" \
  --output-dir "${OUTPUT_DIR}"

python - <<PY
from pathlib import Path
out = Path("${OUTPUT_DIR}")
required = [
    "run_metadata.json",
    "forecast_metrics.json",
    "embedding_metrics.json",
    "lineage_metrics.json",
    "projected_expression.npy",
    "projected_embedding.npy",
    "projected_cluster_labels.csv",
    "state_transition_matrix.csv",
    "lineage_graph_edges.csv",
]
missing = [name for name in required if not (out / name).exists()]
print("missing:", missing)
if missing:
    raise SystemExit(1)
PY

echo "Job finished at: $(date)"
