#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scanvi_map.$JOB_ID.log
#$ -l s_vmem=128G
#$ -l h100=1
# Generic scArches map + comparison + figures via src.pipeline.
# Honors CONFIG + PROFILE env (default config.yaml / server). Requires the
# matching trained model (run_scanvi_train_reference.sh with the same CONFIG/PROFILE).
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${SCVI_CONDA_ENV:-scvi}"
CONFIG="${CONFIG:-config.yaml}"; PROFILE="${PROFILE:-server}"; REQUIRE_CUDA="${REQUIRE_CUDA:-1}"
mkdir -p "${PROJECT_ROOT}/logs"; cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}" PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}" MPLBACKEND=Agg
export OMP_NUM_THREADS="${NSLOTS:-1}" MKL_NUM_THREADS="${NSLOTS:-1}" OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
echo "map | $(date) | host $(hostname) | config ${CONFIG} | profile ${PROFILE}"
[ -f "${CONDA_SH}" ] || { echo "ERROR: conda init not found" >&2; exit 2; }
# shellcheck disable=SC1090
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
if [ "${REQUIRE_CUDA}" = "1" ]; then
  python -c "import torch,sys; print('cuda',torch.cuda.is_available()); sys.exit(0 if torch.cuda.is_available() else 3)"
fi
TMP_CFG="results/_run_config_${PROFILE}.yaml"
python - "${CONFIG}" "${TMP_CFG}" "${PROFILE}" <<'PYX'
import sys, os, yaml
src, dst, prof = sys.argv[1], sys.argv[2], sys.argv[3]
c = yaml.safe_load(open(src)); c["mode"] = prof
os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
open(dst, "w").write(yaml.safe_dump(c, sort_keys=False)); print("wrote", dst)
PYX
python -m src.pipeline --config "${TMP_CFG}"
echo "map DONE $(date)"
