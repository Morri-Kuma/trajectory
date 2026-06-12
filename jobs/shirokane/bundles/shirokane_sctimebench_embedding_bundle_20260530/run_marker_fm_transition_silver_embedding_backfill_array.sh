#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_marker_fm_transition_silver_embedding_backfill.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=128G
#$ -t 1-18

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-1}"
export VECLIB_MAXIMUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Marker-FM scTimeBench embedding backfill"
echo "Started at   : $(date)"
echo "Host         : $(hostname)"
echo "Project root : ${PROJECT_ROOT}"
echo "Task ID      : ${TASK_ID}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

mapfile -t CONFIGS < <(python - <<'PY'
import glob
from pathlib import Path
import yaml

for p in sorted(glob.glob("benchmark/configs/runtime/*marker_fm_silver*_formal.yaml")):
    cfg = yaml.safe_load(Path(p).read_text(encoding="utf-8-sig")) or {}
    if str(cfg.get("method", "")).lower() in {"scnode", "mioflow", "prescient"}:
        print(p)
PY
)

N_CONFIGS="${#CONFIGS[@]}"
if [ "${N_CONFIGS}" -eq 0 ]; then
  echo "ERROR: no projection-capable marker_fm_silver configs found." >&2
  exit 3
fi
if [ "${TASK_ID}" -lt 1 ] || [ "${TASK_ID}" -gt "${N_CONFIGS}" ]; then
  echo "Task ${TASK_ID} is outside config range 1-${N_CONFIGS}; nothing to do."
  exit 0
fi

CONFIG="${CONFIGS[$((TASK_ID - 1))]}"
echo "Config: ${CONFIG}"

python -m benchmark.evaluation.backfill_sctimebench_embeddings \
  --config "${CONFIG}"

echo "============================================================"
echo "Embedding backfill task finished at: $(date)"
echo "============================================================"
