#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_build_gse178325_raw_full_gene_input.$JOB_ID.log
#$ -l s_vmem=128G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
SCRIPT="scripts/build_gse178325_raw_full_gene_input.py"
OUTPUT_H5AD="data/processed/gse178325_human/GSE178325_0618_raw_full_gene_benchmark_input.h5ad"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Job started at : $(date)"
echo "Host           : $(hostname)"
echo "Project root   : ${PROJECT_ROOT}"
echo "Script         : ${SCRIPT}"
echo "Output h5ad    : ${OUTPUT_H5AD}"
echo "Memory request : s_vmem=128G"
echo "Job ID         : ${JOB_ID:-N/A}"
echo "============================================================"

if [ ! -f "${SCRIPT}" ]; then
  echo "ERROR: script not found: ${SCRIPT}" >&2
  exit 2
fi
if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

echo "------------------------------------------------------------"
echo "Environment check"
echo "------------------------------------------------------------"
echo "Python: $(which python)"
python -V
python - <<'PY'
import anndata
import numpy
import pandas
import scipy
print("imports OK")
PY

echo "------------------------------------------------------------"
echo "Building GSE178325 0618 raw/full-gene scGPT source h5ad"
echo "------------------------------------------------------------"
time python "${SCRIPT}"

echo "------------------------------------------------------------"
echo "Output validation"
echo "------------------------------------------------------------"
python - <<'PY'
import anndata as ad
import scipy.sparse as sp
from pathlib import Path

p = Path("data/processed/gse178325_human/GSE178325_0618_raw_full_gene_benchmark_input.h5ad")
if not p.exists():
    raise SystemExit(f"missing output: {p}")
a = ad.read_h5ad(p, backed="r")
print("shape:", a.shape)
print("obs columns:", list(a.obs.columns))
print("var columns:", list(a.var.columns))
print("times:", sorted(map(float, a.obs["abs_day"].unique().tolist())))
print("X dtype:", a.X.dtype)
a.file.close()
PY

echo "============================================================"
echo "Job finished at : $(date)"
echo "============================================================"
