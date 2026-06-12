#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scanvi_score_ood.$JOB_ID.log
#$ -l s_vmem=96G
#$ -l h100=1
# Re-score the COMPLETED server run with a geometric OOD flag on the scANVI latent.
# No retraining; loads models/scanvi_reference_gse242424/. CONFIG/PROFILE default config.yaml/server.
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${SCVI_CONDA_ENV:-scvi}"; CONFIG="${CONFIG:-config.yaml}"; PROFILE="${PROFILE:-server}"
mkdir -p "${PROJECT_ROOT}/logs"; cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}" PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
# shellcheck disable=SC1090
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
python scripts/score_ood.py --config "${CONFIG}" --profile "${PROFILE}"
echo "score_ood DONE $(date)"
