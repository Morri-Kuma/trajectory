#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_prescient_gse178325_0618_hvg.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=128G
#$ -t 1-3
## To run on GPU nodes, enable the appropriate Shirokane GPU resource line
## for your partition, for example one of:
##$ -l gpu=1
##$ -l cuda=1

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
PRESCIENT_REPO="${PRESCIENT_REPO:-${PROJECT_ROOT}/benchmark/methods/PRESCIENT/prescient_module}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
METHOD="prescient"
ADATA="benchmark/inputs/gse178325_scgpt_hvg2000/GSE178325_scGPT_annotated_HVG2000_benchmark_input.h5ad"

case "${SGE_TASK_ID:-${TASK_ID:-1}}" in
  1)
    SCENARIO="A"
    CONFIG="benchmark/configs/prescient_gse178325_observed_0618_hvg2000_A_scgpt_v1.yaml"
    OUTPUT_DIR="benchmark/results/prescient/gse178325_0618_scenario_A_hvg2000"
    ;;
  2)
    SCENARIO="B"
    CONFIG="benchmark/configs/prescient_gse178325_observed_0618_hvg2000_B_scgpt_v1.yaml"
    OUTPUT_DIR="benchmark/results/prescient/gse178325_0618_scenario_B_hvg2000"
    ;;
  3)
    SCENARIO="C"
    CONFIG="benchmark/configs/prescient_gse178325_observed_0618_hvg2000_C_scgpt_v1.yaml"
    OUTPUT_DIR="benchmark/results/prescient/gse178325_0618_scenario_C_hvg2000"
    ;;
  *)
    echo "Unsupported SGE_TASK_ID=${SGE_TASK_ID:-${TASK_ID:-unset}}; expected 1, 2, or 3" >&2
    exit 2
    ;;
esac

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PRESCIENT_REPO="${PRESCIENT_REPO}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Job started at: $(date)"
echo "Host          : $(hostname)"
echo "Project root  : ${TRAJ_PROJECT_ROOT}"
echo "PRESCIENT repo: ${PRESCIENT_REPO}"
echo "Method        : ${METHOD}"
echo "Dataset       : GSE178325 batch 0618"
echo "Scenario      : ${SCENARIO}"
echo "Adata         : ${ADATA}"
echo "Config        : ${CONFIG}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Task ID       : ${SGE_TASK_ID:-${TASK_ID:-N/A}}"
echo "Memory request: s_vmem=128G"
echo "============================================================"

if [ ! -f "${PRESCIENT_REPO}/prescient/train/model.py" ]; then
  echo "ERROR: PRESCIENT_REPO does not contain prescient/train/model.py: ${PRESCIENT_REPO}" >&2
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
print("PRESCIENT dependencies OK")
print("torch:", torch.__version__)
print("cuda_available:", torch.cuda.is_available())
print("cuda_device_count:", torch.cuda.device_count())
if torch.cuda.is_available():
    print("cuda_device_name:", torch.cuda.get_device_name(0))
PY

python - <<PY
import anndata as ad
adata = ad.read_h5ad("${ADATA}", backed="r")
print("adata_shape:", adata.shape)
print("timepoints:", sorted(adata.obs["abs_day"].astype(float).unique()))
print("state_key:", "scgpt_pseudostate_provisional")
print("state_key_values:", sorted(set(adata.obs["scgpt_pseudostate_provisional"].astype(str))))
print("benchmark_input_id:", adata.uns.get("benchmark_input_id"))
adata.file.close()
PY

if [ "${CLEAN_OUTPUT_DIR:-1}" = "1" ]; then
  echo "Cleaning output dir before aligned rerun: ${OUTPUT_DIR}"
  rm -rf "${OUTPUT_DIR}"
  mkdir -p "${OUTPUT_DIR}"
fi

python -m benchmark.evaluation.eval_dispatch \
  --method "${METHOD}" \
  --scenario "${SCENARIO}" \
  --adata "${ADATA}" \
  --method-config "${CONFIG}" \
  --output-dir "${OUTPUT_DIR}"

python - <<PY
from pathlib import Path
import json

out = Path("${OUTPUT_DIR}")
required = [
    "run_metadata.json",
    "forecast_metrics.json",
    "per_timepoint_forecast_metrics.csv",
    "embedding_metrics.json",
    "per_timepoint_embedding_metrics.csv",
    "projected_expression.npy",
    "projected_embedding.npy",
    "projected_cluster_labels.csv",
    "state_transition_matrix.csv",
    "lineage_graph_edges.csv",
    "lineage_metrics.json",
]
missing = [name for name in required if not (out / name).exists()]
print("missing:", missing)
if missing:
    raise SystemExit(1)

with open(out / "lineage_metrics.json", encoding="utf-8") as f:
    lineage = json.load(f)
print("lineage_status:", lineage.get("status"))
PY

echo "Job finished at: $(date)"
