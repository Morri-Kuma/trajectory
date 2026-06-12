#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_representation_metric_backfill.$JOB_ID.log
#$ -l s_vmem=96G
# Backfill representation-space metrics for every completed representation result
# directory by pairing it with its config and calling
# benchmark.evaluation.eval_representation_run_from_config.
#
# Scores both Geneformer and scGPT arms through the same evaluator path. Skips
# result directories that do not yet have projected_expression.npy, so it is safe
# to run while some trajectory tasks are still pending.
#
# Usage:
#   bash jobs/shirokane/run_representation_metric_backfill.sh
#   REPS="rep_geneformer_cls_pca50" bash jobs/shirokane/run_representation_metric_backfill.sh
#   BACKEND=numpy bash jobs/shirokane/run_representation_metric_backfill.sh

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
REPS="${REPS:-rep_geneformer_cls_pca50 rep_scgpt_cls_pca50}"
METHODS="${METHODS:-mioflow scnode prescient}"
SCENARIOS="${SCENARIOS:-A B C}"
BACKEND="${BACKEND:-auto}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

echo "============================================================"
echo "Representation metric backfill started at: $(date)"
echo "Project root : ${PROJECT_ROOT}"
echo "Conda env    : ${CONDA_ENV}"
echo "Reps         : ${REPS}"
echo "Methods      : ${METHODS}"
echo "Scenarios    : ${SCENARIOS}"
echo "Backend      : ${BACKEND}"
echo "============================================================"

CONFIG_DIR="benchmark/configs/representation"
RESULTS_ROOT="benchmark/results/representation_dynamics"

n_done=0
n_skip=0
for rep in ${REPS}; do
  for method in ${METHODS}; do
    for scenario in ${SCENARIOS}; do
      config="${CONFIG_DIR}/${method}_gse230659_${rep}_${scenario}.yaml"
      result_dir="${RESULTS_ROOT}/${method}/gse230659_${rep}_${scenario}"
      if [ ! -f "${config}" ]; then
        echo "SKIP (no config): ${config}"
        n_skip=$((n_skip + 1))
        continue
      fi
      if [ ! -f "${result_dir}/projected_expression.npy" ]; then
        echo "SKIP (no projected_expression.npy): ${result_dir}"
        n_skip=$((n_skip + 1))
        continue
      fi
      echo "=== backfilling ${method}/${rep}/${scenario} ==="
      python -m benchmark.evaluation.eval_representation_run_from_config \
        --config "${config}" \
        --result-dir "${result_dir}" \
        --backend "${BACKEND}"
      n_done=$((n_done + 1))
    done
  done
done

echo "============================================================"
echo "Representation metric backfill done: ${n_done} processed, ${n_skip} skipped"
echo "============================================================"
