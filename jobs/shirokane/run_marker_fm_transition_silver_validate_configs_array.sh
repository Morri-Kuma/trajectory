#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -t 1-20
#$ -o logs/run_marker_fm_transition_silver_validate_configs.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=16G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"

echo "============================================================"
echo "Config validation job started at: $(date)"
echo "Host                         : $(hostname)"
echo "Project root                 : ${PROJECT_ROOT}"
echo "Job ID / task                : ${JOB_ID:-N/A} / ${SGE_TASK_ID:-N/A}"
echo "BLAS threads                 : OPENBLAS=${OPENBLAS_NUM_THREADS}, OMP=${OMP_NUM_THREADS}, MKL=${MKL_NUM_THREADS}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

mapfile -t CONFIGS < <(find benchmark/configs/runtime -maxdepth 1 -type f -name '*marker_fm_silver*_formal.yaml' | sort)

N_CONFIGS="${#CONFIGS[@]}"
TASK_ID="${SGE_TASK_ID:-1}"

echo "Discovered runtime configs: ${N_CONFIGS}"

if [ "${N_CONFIGS}" -eq 0 ]; then
  echo "ERROR: no marker_fm_silver runtime configs found." >&2
  exit 3
fi

if [ "${TASK_ID}" -lt 1 ] || [ "${TASK_ID}" -gt "${N_CONFIGS}" ]; then
  echo "Task ${TASK_ID} is outside config range 1-${N_CONFIGS}; nothing to do."
  exit 0
fi

CONFIG="${CONFIGS[$((TASK_ID - 1))]}"
echo "Validating config ${TASK_ID}/${N_CONFIGS}: ${CONFIG}"

python -m benchmark.annotation.validate_milestone_configs "${CONFIG}"

echo "============================================================"
echo "Config validation job finished at: $(date)"
echo "============================================================"
