#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N v21_embed
#$ -o logs/run_v21_new_datasets_embedding.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G
#$ -t 1-18
# ── Official-silver milestone embedding coherence for GSE298212 + GSE218855 ───
#
# Same role as run_pseudotime_embedding_coherence_array.sh: produce the official
# embedding metric <run_dir>/embedding_milestone_eval/embedding_metrics_official_silver.json
# that summarize_official_silver.py requires for projection methods. Run AFTER
# run_v21_new_datasets_benchmark_array.sh, then add the result dirs to
# benchmark/results/result_manifest.yaml and re-run the lineage graph-sim summary.
#
# 18 tasks = 2 datasets x 3 projection methods x A/B/C. WOT/CellRank2 are
# lineage-only and need no embedding metric.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"

case "${TASK_ID}" in
  1)  DATASET_ID="GSE298212"; METHOD="scnode";    SCENARIO="A" ;;
  2)  DATASET_ID="GSE298212"; METHOD="scnode";    SCENARIO="B" ;;
  3)  DATASET_ID="GSE298212"; METHOD="scnode";    SCENARIO="C" ;;
  4)  DATASET_ID="GSE298212"; METHOD="mioflow";   SCENARIO="A" ;;
  5)  DATASET_ID="GSE298212"; METHOD="mioflow";   SCENARIO="B" ;;
  6)  DATASET_ID="GSE298212"; METHOD="mioflow";   SCENARIO="C" ;;
  7)  DATASET_ID="GSE298212"; METHOD="prescient"; SCENARIO="A" ;;
  8)  DATASET_ID="GSE298212"; METHOD="prescient"; SCENARIO="B" ;;
  9)  DATASET_ID="GSE298212"; METHOD="prescient"; SCENARIO="C" ;;
  10) DATASET_ID="GSE218855"; METHOD="scnode";    SCENARIO="A" ;;
  11) DATASET_ID="GSE218855"; METHOD="scnode";    SCENARIO="B" ;;
  12) DATASET_ID="GSE218855"; METHOD="scnode";    SCENARIO="C" ;;
  13) DATASET_ID="GSE218855"; METHOD="mioflow";   SCENARIO="A" ;;
  14) DATASET_ID="GSE218855"; METHOD="mioflow";   SCENARIO="B" ;;
  15) DATASET_ID="GSE218855"; METHOD="mioflow";   SCENARIO="C" ;;
  16) DATASET_ID="GSE218855"; METHOD="prescient"; SCENARIO="A" ;;
  17) DATASET_ID="GSE218855"; METHOD="prescient"; SCENARIO="B" ;;
  18) DATASET_ID="GSE218855"; METHOD="prescient"; SCENARIO="C" ;;
  *) echo "ERROR: unsupported task id ${TASK_ID}; expected 1-18" >&2; exit 2 ;;
esac

DL=$(echo "${DATASET_ID}" | tr '[:upper:]' '[:lower:]')
OUTPUT_DIR="benchmark/results/${METHOD}/${DL}_marker_fm_silver_${SCENARIO}_hvg2000_formal"
if [ "${DATASET_ID}" = "GSE298212" ]; then
  REFERENCE_H5AD="benchmark/inputs/gse298212_marker_fm_transition_silver_hvg2000/GSE298212_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"
  PROVIDER_ID="gse298212_marker_fm_transition_silver_v1"
else
  REFERENCE_H5AD="benchmark/inputs/gse218855_marker_fm_transition_silver_hvg2000/GSE218855_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"
  PROVIDER_ID="gse218855_marker_fm_transition_silver_v1"
fi
EMBEDDING_OUT="${OUTPUT_DIR}/embedding_milestone_eval"

mkdir -p "${PROJECT_ROOT}/logs"; cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"; export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"; export MKL_NUM_THREADS="${NSLOTS:-1}"; export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "=== v21 embedding coherence: ${DATASET_ID} ${METHOD} ${SCENARIO} (task ${TASK_ID}) at $(date) ==="
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"

if [ ! -f "${REFERENCE_H5AD}" ]; then echo "ERROR: reference h5ad not found: ${REFERENCE_H5AD}" >&2; exit 2; fi
if [ ! -d "${OUTPUT_DIR}" ]; then echo "ERROR: run output dir not found: ${OUTPUT_DIR} (run the benchmark array first)" >&2; exit 2; fi
if [ ! -f "${OUTPUT_DIR}/embedding.npy" ]; then echo "ERROR: embedding.npy missing in ${OUTPUT_DIR}" >&2; exit 2; fi
if [ ! -f "${OUTPUT_DIR}/next_timepoint_embedding.npy" ] && [ ! -f "${OUTPUT_DIR}/projected_embedding.npy" ]; then
  echo "ERROR: projected/next_timepoint embedding missing in ${OUTPUT_DIR}" >&2; exit 2; fi

mkdir -p "${EMBEDDING_OUT}"
python -m benchmark.evaluation.eval_embedding_milestone \
  --input-h5ad "${REFERENCE_H5AD}" \
  --run-output-dir "${OUTPUT_DIR}" \
  --output-dir "${EMBEDDING_OUT}" \
  --dataset-id "${DATASET_ID}" \
  --mode run-output \
  --label-mode official_silver \
  --exclude-label ambiguous \
  --exclude-label unknown_or_ood \
  --leiden-n-neighbors 15 \
  --leiden-resolution 0.5

python - <<PY
import json
from pathlib import Path
m = json.load(open(Path("${EMBEDDING_OUT}") / "embedding_metrics_official_silver.json", encoding="utf-8"))
print(json.dumps({k: m.get(k) for k in
      ("status","label_mode","provider_id","state_key","n_cells_evaluated",
       "adjusted_rand_index","pred_tp_avg_normalized_entropy")}, indent=2))
assert m.get("status") == "completed", m.get("status")
assert m.get("provider_id") == "${PROVIDER_ID}", m.get("provider_id")
assert m.get("state_key") == "final_milestone_label_coarse", m.get("state_key")
PY
echo "=== task ${TASK_ID} DONE at $(date); wrote ${EMBEDDING_OUT}/embedding_metrics_official_silver.json ==="
