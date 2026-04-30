#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_mioflow_cpu_A_hvg.$JOB_ID.log
#$ -l s_vmem=48G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"

CONFIG="benchmark/configs/mioflow_gse230659_observed_scgpt_v1_A_hvg2000_cpu.yaml"
ADATA="benchmark/inputs/gse230659_scgpt_hvg2000/GSE230659_scGPT_annotated_HVG2000_benchmark_input.h5ad"
OUTPUT_DIR="benchmark/results/mioflow/scenario_A_scgpt_v1_hvg2000_cpu"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export CUDA_VISIBLE_DEVICES=""
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Job started at: $(date)"
echo "Host          : $(hostname)"
echo "Project root  : ${TRAJ_PROJECT_ROOT}"
echo "Method        : MIOFlow"
echo "Scenario      : A HVG2000 CPU reduced validation"
echo "Config        : ${CONFIG}"
echo "Output dir    : ${OUTPUT_DIR}"
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
PY

python benchmark/evaluation/eval_dispatch.py \
  --method mioflow \
  --scenario A \
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
