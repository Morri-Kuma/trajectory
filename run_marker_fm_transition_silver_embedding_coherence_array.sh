#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_marker_fm_transition_silver_embedding.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G
#$ -t 1-18

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"

case "${TASK_ID}" in
  1)  DATASET_ID="GSE178325"; METHOD="scnode";    SCENARIO="A"; OUTPUT_DIR="benchmark/results/scnode/gse178325_marker_fm_silver_A_hvg2000_formal" ;;
  2)  DATASET_ID="GSE178325"; METHOD="scnode";    SCENARIO="B"; OUTPUT_DIR="benchmark/results/scnode/gse178325_marker_fm_silver_B_hvg2000_formal" ;;
  3)  DATASET_ID="GSE178325"; METHOD="scnode";    SCENARIO="C"; OUTPUT_DIR="benchmark/results/scnode/gse178325_marker_fm_silver_C_hvg2000_formal" ;;
  4)  DATASET_ID="GSE178325"; METHOD="mioflow";   SCENARIO="A"; OUTPUT_DIR="benchmark/results/mioflow/gse178325_marker_fm_silver_A_hvg2000_formal" ;;
  5)  DATASET_ID="GSE178325"; METHOD="mioflow";   SCENARIO="B"; OUTPUT_DIR="benchmark/results/mioflow/gse178325_marker_fm_silver_B_hvg2000_formal" ;;
  6)  DATASET_ID="GSE178325"; METHOD="mioflow";   SCENARIO="C"; OUTPUT_DIR="benchmark/results/mioflow/gse178325_marker_fm_silver_C_hvg2000_formal" ;;
  7)  DATASET_ID="GSE178325"; METHOD="prescient"; SCENARIO="A"; OUTPUT_DIR="benchmark/results/prescient/gse178325_marker_fm_silver_A_hvg2000_formal" ;;
  8)  DATASET_ID="GSE178325"; METHOD="prescient"; SCENARIO="B"; OUTPUT_DIR="benchmark/results/prescient/gse178325_marker_fm_silver_B_hvg2000_formal" ;;
  9)  DATASET_ID="GSE178325"; METHOD="prescient"; SCENARIO="C"; OUTPUT_DIR="benchmark/results/prescient/gse178325_marker_fm_silver_C_hvg2000_formal" ;;
  10) DATASET_ID="GSE230659"; METHOD="scnode";    SCENARIO="A"; OUTPUT_DIR="benchmark/results/scnode/gse230659_marker_fm_silver_A_hvg2000_formal" ;;
  11) DATASET_ID="GSE230659"; METHOD="scnode";    SCENARIO="B"; OUTPUT_DIR="benchmark/results/scnode/gse230659_marker_fm_silver_B_hvg2000_formal" ;;
  12) DATASET_ID="GSE230659"; METHOD="scnode";    SCENARIO="C"; OUTPUT_DIR="benchmark/results/scnode/gse230659_marker_fm_silver_C_hvg2000_formal" ;;
  13) DATASET_ID="GSE230659"; METHOD="mioflow";   SCENARIO="A"; OUTPUT_DIR="benchmark/results/mioflow/gse230659_marker_fm_silver_A_hvg2000_formal" ;;
  14) DATASET_ID="GSE230659"; METHOD="mioflow";   SCENARIO="B"; OUTPUT_DIR="benchmark/results/mioflow/gse230659_marker_fm_silver_B_hvg2000_formal" ;;
  15) DATASET_ID="GSE230659"; METHOD="mioflow";   SCENARIO="C"; OUTPUT_DIR="benchmark/results/mioflow/gse230659_marker_fm_silver_C_hvg2000_formal" ;;
  16) DATASET_ID="GSE230659"; METHOD="prescient"; SCENARIO="A"; OUTPUT_DIR="benchmark/results/prescient/gse230659_marker_fm_silver_A_hvg2000_formal" ;;
  17) DATASET_ID="GSE230659"; METHOD="prescient"; SCENARIO="B"; OUTPUT_DIR="benchmark/results/prescient/gse230659_marker_fm_silver_B_hvg2000_formal" ;;
  18) DATASET_ID="GSE230659"; METHOD="prescient"; SCENARIO="C"; OUTPUT_DIR="benchmark/results/prescient/gse230659_marker_fm_silver_C_hvg2000_formal" ;;
  *) echo "ERROR: unsupported task id ${TASK_ID}; expected 1-18" >&2; exit 2 ;;
esac

if [ "${DATASET_ID}" = "GSE178325" ]; then
  REFERENCE_H5AD="benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"
  PROVIDER_ID="gse178325_marker_fm_transition_silver_v1"
else
  REFERENCE_H5AD="benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"
  PROVIDER_ID="gse230659_marker_fm_transition_silver_v1"
fi

EMBEDDING_OUT="${OUTPUT_DIR}/embedding_milestone_eval"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Job started at : $(date)"
echo "Task ID        : ${TASK_ID}"
echo "Dataset        : ${DATASET_ID}"
echo "Method         : ${METHOD}"
echo "Scenario       : ${SCENARIO}"
echo "Provider       : ${PROVIDER_ID}"
echo "Reference h5ad : ${REFERENCE_H5AD}"
echo "Output dir     : ${OUTPUT_DIR}"
echo "Embedding out  : ${EMBEDDING_OUT}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

if [ ! -f "${REFERENCE_H5AD}" ]; then
  echo "ERROR: reference h5ad not found: ${REFERENCE_H5AD}" >&2
  exit 2
fi
if [ ! -d "${OUTPUT_DIR}" ]; then
  echo "ERROR: formal run output dir not found: ${OUTPUT_DIR}" >&2
  exit 2
fi
if [ ! -f "${OUTPUT_DIR}/projected_cluster_labels.csv" ]; then
  echo "ERROR: projected_cluster_labels.csv missing in ${OUTPUT_DIR}; run primary array first." >&2
  exit 2
fi
if [ ! -f "${OUTPUT_DIR}/projected_milestone_labels.csv" ]; then
  echo "ERROR: projected_milestone_labels.csv missing in ${OUTPUT_DIR}; rerun primary task to create it." >&2
  exit 2
fi

mkdir -p "${EMBEDDING_OUT}"

echo "------------------------------------------------------------"
echo "Running official_silver milestone embedding coherence"
echo "------------------------------------------------------------"
python -m benchmark.evaluation.eval_embedding_milestone \
  --input-h5ad "${REFERENCE_H5AD}" \
  --run-output-dir "${OUTPUT_DIR}" \
  --output-dir "${EMBEDDING_OUT}" \
  --dataset-id "${DATASET_ID}" \
  --mode run-output \
  --label-mode official_silver \
  --exclude-label ambiguous \
  --exclude-label unknown_or_ood

echo "------------------------------------------------------------"
echo "Validating official_silver embedding output"
echo "------------------------------------------------------------"
python - <<PY
import json
from pathlib import Path

out = Path("${EMBEDDING_OUT}")
metrics_path = out / "embedding_metrics_official_silver.json"
if not metrics_path.exists():
    raise SystemExit(f"missing {metrics_path}")
m = json.load(open(metrics_path, encoding="utf-8"))
print(json.dumps({
    "status": m.get("status"),
    "label_mode": m.get("label_mode"),
    "provider_id": m.get("provider_id"),
    "state_key": m.get("state_key"),
    "n_cells_evaluated": m.get("n_cells_evaluated"),
    "adjusted_rand_index": m.get("adjusted_rand_index"),
    "mean_normalized_entropy": m.get("mean_normalized_entropy"),
}, indent=2))
if m.get("status") != "completed":
    raise SystemExit(f"embedding metrics not completed: {m.get('status')}")
if m.get("label_mode") != "official_silver":
    raise SystemExit(f"unexpected label_mode: {m.get('label_mode')}")
if m.get("provider_id") != "${PROVIDER_ID}":
    raise SystemExit(f"unexpected provider_id: {m.get('provider_id')}")
if m.get("state_key") != "final_milestone_label_expanded":
    raise SystemExit(f"unexpected state_key: {m.get('state_key')}")
if not isinstance(m.get("n_cells_evaluated"), int) or m["n_cells_evaluated"] <= 0:
    raise SystemExit(f"invalid n_cells_evaluated: {m.get('n_cells_evaluated')}")
PY

echo "============================================================"
echo "Embedding coherence task finished at: $(date)"
echo "Output dir: ${OUTPUT_DIR}"
echo "============================================================"
