#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_mioflow_formal_BC_hvg.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=128G
#$ -t 1-2
## To run on GPU nodes, enable the appropriate Shirokane GPU resource line
## for your partition, for example one of:
##$ -l gpu=1
##$ -l cuda=1
## Keep the line commented until you confirm the resource name with qstat/qconf.

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
ADATA="benchmark/inputs/gse230659_scgpt_hvg2000/GSE230659_scGPT_annotated_HVG2000_benchmark_input.h5ad"

case "${SGE_TASK_ID:-${TASK_ID:-1}}" in
  1)
    SCENARIO="B"
    CONFIG="benchmark/configs/mioflow_gse230659_observed_scgpt_v1_scenarioB_hvg2000_formal.yaml"
    OUTPUT_DIR="benchmark/results/mioflow/scenario_B_scgpt_v1_hvg2000_formal"
    ;;
  2)
    SCENARIO="C"
    CONFIG="benchmark/configs/mioflow_gse230659_observed_scgpt_v1_scenarioC_hvg2000_formal.yaml"
    OUTPUT_DIR="benchmark/results/mioflow/scenario_C_scgpt_v1_hvg2000_formal"
    ;;
  *)
    echo "Unsupported SGE_TASK_ID=${SGE_TASK_ID:-${TASK_ID:-unset}}; expected 1 or 2" >&2
    exit 2
    ;;
esac

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Job started at: $(date)"
echo "Host          : $(hostname)"
echo "Project root  : ${TRAJ_PROJECT_ROOT}"
echo "Method        : MIOFlow"
echo "Scenario      : ${SCENARIO} HVG2000 formal"
echo "Config        : ${CONFIG}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Task ID       : ${SGE_TASK_ID:-${TASK_ID:-N/A}}"
echo "Memory request: s_vmem=128G"
echo "============================================================"

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

python - <<'PY'
import torch
import torchdiffeq
import torchsde
import ot
print("MIOFlow dependencies OK")
print("torch:", torch.__version__)
print("cuda_available:", torch.cuda.is_available())
print("cuda_device_count:", torch.cuda.device_count())
if torch.cuda.is_available():
    print("cuda_device_name:", torch.cuda.get_device_name(0))
PY

python benchmark/evaluation/eval_dispatch.py \
  --method mioflow \
  --scenario "${SCENARIO}" \
  --adata "${ADATA}" \
  --method-config "${CONFIG}" \
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
