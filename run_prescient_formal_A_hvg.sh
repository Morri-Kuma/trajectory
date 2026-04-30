#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_prescient_formal_A_hvg.$JOB_ID.log
#$ -l s_vmem=128G
## To run on a GPU node, enable the appropriate Shirokane GPU resource line
## for your partition, for example one of:
##$ -l gpu=1
##$ -l cuda=1
## Keep the line commented until you confirm the resource name with qstat/qconf.

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
PRESCIENT_REPO="${PRESCIENT_REPO:-${PROJECT_ROOT}/benchmark/methods/PRESCIENT/prescient_module}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"

CONFIG="benchmark/configs/prescient_gse230659_observed_scgpt_v1_A_hvg2000_formal.yaml"
ADATA="benchmark/inputs/gse230659_scgpt_hvg2000/GSE230659_scGPT_annotated_HVG2000_benchmark_input.h5ad"
OUTPUT_DIR="benchmark/results/prescient/scenario_A_scgpt_v1_hvg2000_formal"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PRESCIENT_REPO="${PRESCIENT_REPO}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"

echo "============================================================"
echo "Job started at: $(date)"
echo "Host          : $(hostname)"
echo "Project root  : ${TRAJ_PROJECT_ROOT}"
echo "PRESCIENT repo: ${PRESCIENT_REPO}"
echo "Scenario      : A HVG2000 formal"
echo "Config        : ${CONFIG}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Memory request: s_vmem=128G"
echo "============================================================"

if [ ! -f "${PRESCIENT_REPO}/prescient/train/model.py" ]; then
  echo "ERROR: PRESCIENT_REPO does not contain prescient/train/model.py: ${PRESCIENT_REPO}" >&2
  echo "Sync benchmark/methods/PRESCIENT/prescient_module or set PRESCIENT_REPO to the original PRESCIENT checkout." >&2
  exit 2
fi

if [ "${PRESCIENT_CLEAN_INTERNAL:-1}" = "1" ]; then
  echo "Cleaning PRESCIENT internal cache: ${OUTPUT_DIR}/prescient_internal"
  rm -rf "${OUTPUT_DIR}/prescient_internal"
fi

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

python - <<'PY'
import torch
print("torch:", torch.__version__)
print("cuda_available:", torch.cuda.is_available())
print("cuda_device_count:", torch.cuda.device_count())
if torch.cuda.is_available():
    print("cuda_device_name:", torch.cuda.get_device_name(0))
PY

python benchmark/evaluation/eval_dispatch.py \
  --method prescient \
  --scenario A \
  --adata "${ADATA}" \
  --method-config "${CONFIG}" \
  --output-dir "${OUTPUT_DIR}"

python - <<'PY'
from pathlib import Path
out = Path("benchmark/results/prescient/scenario_A_scgpt_v1_hvg2000_formal")
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
