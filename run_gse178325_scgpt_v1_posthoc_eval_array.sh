#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_gse178325_scgpt_v1_posthoc_eval_array.$JOB_ID.$TASK_ID.log
#$ -t 1-9
#$ -l s_vmem=128G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

case "${SGE_TASK_ID}" in
  1) METHOD="mioflow";   SCENARIO="A"; MODULE="benchmark.methods.MIOFlow.run" ;;
  2) METHOD="mioflow";   SCENARIO="B"; MODULE="benchmark.methods.MIOFlow.run" ;;
  3) METHOD="mioflow";   SCENARIO="C"; MODULE="benchmark.methods.MIOFlow.run" ;;
  4) METHOD="prescient"; SCENARIO="A"; MODULE="benchmark.methods.PRESCIENT.run" ;;
  5) METHOD="prescient"; SCENARIO="B"; MODULE="benchmark.methods.PRESCIENT.run" ;;
  6) METHOD="prescient"; SCENARIO="C"; MODULE="benchmark.methods.PRESCIENT.run" ;;
  7) METHOD="scnode";    SCENARIO="A"; MODULE="benchmark.methods.scNODE.run" ;;
  8) METHOD="scnode";    SCENARIO="B"; MODULE="benchmark.methods.scNODE.run" ;;
  9) METHOD="scnode";    SCENARIO="C"; MODULE="benchmark.methods.scNODE.run" ;;
  *) echo "Invalid SGE_TASK_ID=${SGE_TASK_ID}" >&2; exit 2 ;;
esac

CONFIG="benchmark/configs/${METHOD}_gse178325_observed_0618_hvg2000_${SCENARIO}_scgpt_v1.yaml"
INPUT_H5AD="benchmark/inputs/gse178325_human_hvg2000/GSE178325_human_HVG2000_scGPTv1_HVG2000_benchmark_input.h5ad"
PROVIDER_JSON="benchmark/ground_truth/providers/gse178325_0618_scgpt_v1/reference_graph.json"

echo "============================================================"
echo "Job started at : $(date)"
echo "Host           : $(hostname)"
echo "Project root   : ${PROJECT_ROOT}"
echo "Task ID        : ${SGE_TASK_ID}"
echo "Method         : ${METHOD}"
echo "Scenario       : ${SCENARIO}"
echo "Config         : ${CONFIG}"
echo "Input h5ad     : ${INPUT_H5AD}"
echo "Provider graph : ${PROVIDER_JSON}"
echo "Memory request : s_vmem=128G"
echo "============================================================"

for path in "${CONFIG}" "${INPUT_H5AD}" "${PROVIDER_JSON}"; do
  if [ ! -f "${path}" ]; then
    echo "ERROR: required file not found: ${path}" >&2
    exit 2
  fi
done

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

echo "------------------------------------------------------------"
echo "Environment check"
echo "------------------------------------------------------------"
echo "Python: $(which python)"
python -V
python - <<'PY'
import anndata
import numpy
import pandas
import torch
print("imports OK")
print("anndata:", anndata.__version__)
print("numpy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("torch:", torch.__version__)
print("cuda_available:", torch.cuda.is_available())
PY

echo "------------------------------------------------------------"
echo "Running ${METHOD} ${SCENARIO} posthoc scGPT-v1 rerun"
echo "------------------------------------------------------------"
time python -m "${MODULE}" --config "${CONFIG}"

echo "------------------------------------------------------------"
echo "Output validation"
echo "------------------------------------------------------------"
python - <<PY
import json
from pathlib import Path
import yaml

cfg = yaml.safe_load(Path("${CONFIG}").read_text(encoding="utf-8-sig"))
out = Path(cfg["output"]["base_dir"])
required = [
    "embedding_metrics.json",
    "per_timepoint_embedding_metrics.csv",
    "projected_cluster_labels.csv",
    "state_transition_matrix.csv",
    "lineage_graph_edges.csv",
    "lineage_metrics.json",
]
missing = [name for name in required if not (out / name).exists()]
print("output_dir:", out)
print("missing:", missing)
if missing:
    raise SystemExit(1)
lineage = json.loads((out / "lineage_metrics.json").read_text(encoding="utf-8"))
embedding = json.loads((out / "embedding_metrics.json").read_text(encoding="utf-8"))
print("lineage_status:", lineage.get("status"))
print("lineage_provider:", (lineage.get("ground_truth") or {}).get("provider_id"))
print("embedding_status:", embedding.get("status"))
if lineage.get("status") != "completed":
    raise SystemExit(1)
PY

echo "============================================================"
echo "Job finished at : $(date)"
echo "============================================================"
