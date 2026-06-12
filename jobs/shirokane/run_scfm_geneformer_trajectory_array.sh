#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scfm_geneformer_trajectory.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=96G
#$ -t 1-9

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"

case "${TASK_ID}" in
  1) METHOD="mioflow";   SCENARIO="A"; CONFIG="benchmark/configs/representation/mioflow_gse230659_rep_geneformer_cls_pca50_A.yaml" ;;
  2) METHOD="mioflow";   SCENARIO="B"; CONFIG="benchmark/configs/representation/mioflow_gse230659_rep_geneformer_cls_pca50_B.yaml" ;;
  3) METHOD="mioflow";   SCENARIO="C"; CONFIG="benchmark/configs/representation/mioflow_gse230659_rep_geneformer_cls_pca50_C.yaml" ;;
  4) METHOD="scnode";    SCENARIO="A"; CONFIG="benchmark/configs/representation/scnode_gse230659_rep_geneformer_cls_pca50_A.yaml" ;;
  5) METHOD="scnode";    SCENARIO="B"; CONFIG="benchmark/configs/representation/scnode_gse230659_rep_geneformer_cls_pca50_B.yaml" ;;
  6) METHOD="scnode";    SCENARIO="C"; CONFIG="benchmark/configs/representation/scnode_gse230659_rep_geneformer_cls_pca50_C.yaml" ;;
  7) METHOD="prescient"; SCENARIO="A"; CONFIG="benchmark/configs/representation/prescient_gse230659_rep_geneformer_cls_pca50_A.yaml" ;;
  8) METHOD="prescient"; SCENARIO="B"; CONFIG="benchmark/configs/representation/prescient_gse230659_rep_geneformer_cls_pca50_B.yaml" ;;
  9) METHOD="prescient"; SCENARIO="C"; CONFIG="benchmark/configs/representation/prescient_gse230659_rep_geneformer_cls_pca50_C.yaml" ;;
  *)
    echo "ERROR: unsupported task id ${TASK_ID}; expected 1-9" >&2
    exit 2
    ;;
esac

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export CONFIG
export METHOD
export SCENARIO
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Geneformer trajectory task started at: $(date)"
echo "Task ID      : ${TASK_ID}/9"
echo "Method       : ${METHOD}"
echo "Scenario     : ${SCENARIO}"
echo "Config       : ${CONFIG}"
echo "Project root : ${PROJECT_ROOT}"
echo "Conda env    : ${CONDA_ENV}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

if [ ! -f "${CONFIG}" ]; then
  echo "ERROR: config not found: ${CONFIG}" >&2
  exit 2
fi

python - <<'PY'
import os
from pathlib import Path

import yaml

root = Path(os.environ["TRAJ_PROJECT_ROOT"])
cfg_path = root / os.environ["CONFIG"]
with open(cfg_path, encoding="utf-8-sig") as f:
    cfg = yaml.safe_load(f) or {}
h5ad = root / cfg["dataset"]["h5ad_path"]
if not h5ad.exists():
    raise SystemExit(f"config input h5ad missing: {h5ad}")
print("validated config input:", h5ad)
PY

case "${METHOD}" in
  mioflow)
    python -m benchmark.methods.MIOFlow.run --config "${CONFIG}"
    ;;
  scnode)
    python -m benchmark.methods.scNODE.run --config "${CONFIG}"
    ;;
  prescient)
    python -m benchmark.methods.PRESCIENT.run --config "${CONFIG}"
    ;;
  *)
    echo "ERROR: unsupported method ${METHOD}" >&2
    exit 2
    ;;
esac

echo "============================================================"
echo "Geneformer trajectory task finished at: $(date)"
echo "Method: ${METHOD}"
echo "Scenario: ${SCENARIO}"
echo "Config: ${CONFIG}"
echo "============================================================"
