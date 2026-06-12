#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scanvi_train_reference.$JOB_ID.log
#$ -l s_vmem=128G
#$ -l h100=1
# ── scANVI annotation branch — STEP 2/3: train reference scANVI (GPU) ─────────
# Trains scVI -> scANVI on the prepped GSE242424 reference and saves the model.
# Hard-fails if the reference matrix is not raw integer counts, and enforces the
# reference label-recovery gate before any query mapping.
# GPU flag: adjust `-l h100=1` to your queue (a100=1 / v100=1).
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${SCVI_CONDA_ENV:-scvi}"
CONFIG="${CONFIG:-config.yaml}"
PROFILE="${PROFILE:-server}"
REQUIRE_CUDA="${REQUIRE_CUDA:-1}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}" MKL_NUM_THREADS="${NSLOTS:-1}" OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "scANVI train reference | started $(date) | host $(hostname)"
echo "Project root : ${PROJECT_ROOT} | env ${CONDA_ENV} | config ${CONFIG} | profile ${PROFILE}"
echo "============================================================"
[ -f "${CONDA_SH}" ] || { echo "ERROR: conda init not found: ${CONDA_SH}" >&2; exit 2; }
# shellcheck disable=SC1090
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

if [ "${REQUIRE_CUDA}" = "1" ]; then
  python - <<'PYX'
import torch, sys
print("CUDA available:", torch.cuda.is_available())
sys.exit(0 if torch.cuda.is_available() else 3)
PYX
fi

python scripts/train_scanvi_reference.py --config "${CONFIG}" --profile "${PROFILE}"
echo "scANVI train reference DONE $(date)"
