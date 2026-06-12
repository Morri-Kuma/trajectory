#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scanvi_positive_control.$JOB_ID.log
#$ -l s_vmem=128G
#$ -l h100=1
# Positive control: scANVI self-transfer on GSE175634 (labels known). Decisive.
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"; CONDA_ENV="${SCVI_CONDA_ENV:-scvi}"
mkdir -p "${PROJECT_ROOT}/logs"; cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}" PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
# shellcheck disable=SC1090
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
python -c "import torch,sys; sys.exit(0 if torch.cuda.is_available() else 3)"
python scripts/positive_control_gse175634.py
echo "positive control DONE $(date)"
