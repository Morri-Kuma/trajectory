#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N pt_embed
#$ -o logs/run_pseudotime_embedding.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G
#$ -t 1-18
# ── Official-silver milestone embedding coherence for pseudotime D/E/F ────────
#
# The official embedding metric the summarizer requires is
#   <run_dir>/embedding_milestone_eval/embedding_metrics_official_silver.json
# produced by benchmark.evaluation.eval_embedding_milestone — NOT the base
# embedding_metrics.json the dispatcher writes. A/B/C went through this step via
# run_marker_fm_transition_silver_embedding_coherence_array.sh; the pseudotime
# D/E/F projection runs need the same step before the summary will include them.
#
# This mirrors that job exactly (same eval, same Leiden 15/0.5, same provider/
# state_key) but points at the pseudotime result dirs and the pseudotime input
# h5ad (so embedding.npy rows match the cell universe each run trained on).
#
# 18 tasks = 2 datasets x 3 projection methods x D/E/F. WOT/CellRank2 are
# lineage-only and need no embedding metric.
#
# Run AFTER run_pseudotime_benchmark_array.sh and BEFORE the summary:
#   qsub jobs/shirokane/run_pseudotime_embedding_coherence_array.sh
#   qsub -hold_jid pt_embed jobs/shirokane/run_silver_lineage_graphsim_summary.sh
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"

case "${TASK_ID}" in
  1)  DATASET_ID="GSE178325"; METHOD="scnode";    SCENARIO="D" ;;
  2)  DATASET_ID="GSE178325"; METHOD="scnode";    SCENARIO="E" ;;
  3)  DATASET_ID="GSE178325"; METHOD="scnode";    SCENARIO="F" ;;
  4)  DATASET_ID="GSE178325"; METHOD="mioflow";   SCENARIO="D" ;;
  5)  DATASET_ID="GSE178325"; METHOD="mioflow";   SCENARIO="E" ;;
  6)  DATASET_ID="GSE178325"; METHOD="mioflow";   SCENARIO="F" ;;
  7)  DATASET_ID="GSE178325"; METHOD="prescient"; SCENARIO="D" ;;
  8)  DATASET_ID="GSE178325"; METHOD="prescient"; SCENARIO="E" ;;
  9)  DATASET_ID="GSE178325"; METHOD="prescient"; SCENARIO="F" ;;
  10) DATASET_ID="GSE230659"; METHOD="scnode";    SCENARIO="D" ;;
  11) DATASET_ID="GSE230659"; METHOD="scnode";    SCENARIO="E" ;;
  12) DATASET_ID="GSE230659"; METHOD="scnode";    SCENARIO="F" ;;
  13) DATASET_ID="GSE230659"; METHOD="mioflow";   SCENARIO="D" ;;
  14) DATASET_ID="GSE230659"; METHOD="mioflow";   SCENARIO="E" ;;
  15) DATASET_ID="GSE230659"; METHOD="mioflow";   SCENARIO="F" ;;
  16) DATASET_ID="GSE230659"; METHOD="prescient"; SCENARIO="D" ;;
  17) DATASET_ID="GSE230659"; METHOD="prescient"; SCENARIO="E" ;;
  18) DATASET_ID="GSE230659"; METHOD="prescient"; SCENARIO="F" ;;
  *) echo "ERROR: unsupported task id ${TASK_ID}; expected 1-18" >&2; exit 2 ;;
esac

DL=$(echo "${DATASET_ID}" | tr '[:upper:]' '[:lower:]')
OUTPUT_DIR="benchmark/results/${METHOD}/${DL}_marker_fm_silver_${SCENARIO}_pseudotime_formal"

if [ "${DATASET_ID}" = "GSE178325" ]; then
  REFERENCE_H5AD="benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_marker_fm_transition_silver_HVG2000_pseudotime_benchmark_input.h5ad"
  PROVIDER_ID="gse178325_marker_fm_transition_silver_v1"
else
  REFERENCE_H5AD="benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_marker_fm_transition_silver_HVG2000_pseudotime_benchmark_input.h5ad"
  PROVIDER_ID="gse230659_marker_fm_transition_silver_v1"
fi
EMBEDDING_OUT="${OUTPUT_DIR}/embedding_milestone_eval"

mkdir -p "${PROJECT_ROOT}/logs"; cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"; export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"; export MKL_NUM_THREADS="${NSLOTS:-1}"; export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "=== pseudotime embedding coherence: ${DATASET_ID} ${METHOD} ${SCENARIO} (task ${TASK_ID}) at $(date) ==="
echo "  run dir   : ${OUTPUT_DIR}"
echo "  reference : ${REFERENCE_H5AD}"

source "${CONDA_SH}"; conda activate "${CONDA_ENV}"

if [ ! -f "${REFERENCE_H5AD}" ]; then echo "ERROR: pseudotime reference h5ad not found: ${REFERENCE_H5AD} (run run_pseudotime_axis_build.sh)" >&2; exit 2; fi
if [ ! -d "${OUTPUT_DIR}" ]; then echo "ERROR: run output dir not found: ${OUTPUT_DIR} (run run_pseudotime_benchmark_array.sh)" >&2; exit 2; fi
if [ ! -f "${OUTPUT_DIR}/embedding.npy" ]; then echo "ERROR: embedding.npy missing in ${OUTPUT_DIR}" >&2; exit 2; fi
if [ ! -f "${OUTPUT_DIR}/next_timepoint_embedding.npy" ] && [ ! -f "${OUTPUT_DIR}/projected_embedding.npy" ]; then
  echo "ERROR: next_timepoint_embedding.npy/projected_embedding.npy missing in ${OUTPUT_DIR}" >&2; exit 2; fi

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
