#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scanvi_map_compare.$JOB_ID.log
#$ -l s_vmem=128G
#$ -l h100=1
# ── scANVI annotation branch — STEP 3/3: map queries + compare (GPU) ──────────
# scArches-maps GSE178325/GSE230659 onto the frozen reference scANVI model and
# runs the full annotation-comparison + trajectory + figures via src.pipeline in
# `server` mode. Requires STEP 2 to have produced the reference model.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${SCVI_CONDA_ENV:-scvi}"
CONFIG="${CONFIG:-config.yaml}"
REQUIRE_CUDA="${REQUIRE_CUDA:-1}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}" MKL_NUM_THREADS="${NSLOTS:-1}" OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
export MPLBACKEND=Agg

echo "============================================================"
echo "scANVI map + compare | started $(date) | host $(hostname)"
echo "Project root : ${PROJECT_ROOT} | env ${CONDA_ENV} | config ${CONFIG}"
echo "============================================================"
[ -f "${CONDA_SH}" ] || { echo "ERROR: conda init not found: ${CONDA_SH}" >&2; exit 2; }
# shellcheck disable=SC1090
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

if [ "${REQUIRE_CUDA}" = "1" ]; then
  python - <<'PYX'
import torch, sys
print("CUDA available:", torch.cuda.is_available()); sys.exit(0 if torch.cuda.is_available() else 3)
PYX
fi

# Run src.pipeline in server mode (generate a temp config with mode: server).
TMP_CFG="results/annotation_branch/_server_config.yaml"
python - "${CONFIG}" "${TMP_CFG}" <<'PYX'
import sys, os, yaml
src, dst = sys.argv[1], sys.argv[2]
c = yaml.safe_load(open(src)); c["mode"] = "server"
os.makedirs(os.path.dirname(dst), exist_ok=True)
open(dst, "w").write(yaml.safe_dump(c, sort_keys=False))
print("wrote", dst)
PYX

python -m src.pipeline --config "${TMP_CFG}"
echo "scANVI map + compare DONE $(date)"
