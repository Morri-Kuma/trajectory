#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scfm_geneformer_build_xrep.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=96G
#$ -t 1-3

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"

INPUT_H5AD="${INPUT_H5AD:-benchmark/inputs/representation/gse230659/source/GSE230659_geneformer_cls_full_gene.h5ad}"
REPRESENTATION_ID="${REPRESENTATION_ID:-rep_geneformer_cls_pca50}"
TIME_KEY="${TIME_KEY:-abs_day}"
SEED="${SEED:-42}"

case "${TASK_ID}" in
  1)
    SCENARIO="A"
    TRAIN_TIMES=(0.5 2.0 4.0 8.0 12.0 16.0 16.33 16.67 17.0 18.0 20.0 22.0 24.0 28.0 30.0)
    ;;
  2)
    SCENARIO="B"
    TRAIN_TIMES=(0.5 2.0 4.0 8.0 12.0 16.0)
    ;;
  3)
    SCENARIO="C"
    TRAIN_TIMES=(0.5 2.0 8.0 16.0 16.33 16.67 17.0 18.0 22.0 24.0)
    ;;
  *)
    echo "ERROR: unsupported task id ${TASK_ID}; expected 1-3" >&2
    exit 2
    ;;
esac

OUTPUT_H5AD="benchmark/inputs/representation/gse230659/${REPRESENTATION_ID}/${SCENARIO}/GSE230659_${REPRESENTATION_ID}_${SCENARIO}_X_rep.h5ad"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export INPUT_H5AD
export OUTPUT_H5AD
export REPRESENTATION_ID
export TIME_KEY
export SEED
export SCENARIO
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Build Geneformer X_rep started at: $(date)"
echo "Task ID           : ${TASK_ID}/3"
echo "Scenario          : ${SCENARIO}"
echo "Input h5ad        : ${INPUT_H5AD}"
echo "Output h5ad       : ${OUTPUT_H5AD}"
echo "Representation ID : ${REPRESENTATION_ID}"
echo "Time key          : ${TIME_KEY}"
echo "Train times       : ${TRAIN_TIMES[*]}"
echo "Seed              : ${SEED}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

if [ ! -f "${INPUT_H5AD}" ]; then
  echo "ERROR: input h5ad not found: ${INPUT_H5AD}" >&2
  exit 2
fi

python -m benchmark.representations.build_representation_inputs \
  --input-h5ad "${INPUT_H5AD}" \
  --output-h5ad "${OUTPUT_H5AD}" \
  --representation-id "${REPRESENTATION_ID}" \
  --time-key "${TIME_KEY}" \
  --fit-scope train_only \
  --train-times "${TRAIN_TIMES[@]}" \
  --seed "${SEED}"

python - <<'PY'
import os
from pathlib import Path

import anndata as ad

path = Path(os.environ["OUTPUT_H5AD"])
adata = ad.read_h5ad(path, backed="r")
print("output:", path)
print("shape:", adata.shape)
print("obsm keys:", list(adata.obsm.keys()))
if "X_rep" not in adata.obsm:
    raise SystemExit("X_rep missing from output h5ad")
print("X_rep shape:", adata.obsm["X_rep"].shape)
required = "final_milestone_label_coarse"
if required not in adata.obs.columns:
    raise SystemExit(f"{required} missing from obs")
adata.file.close()
PY

echo "============================================================"
echo "Build Geneformer X_rep finished at: $(date)"
echo "Scenario: ${SCENARIO}"
echo "Output h5ad: ${OUTPUT_H5AD}"
echo "============================================================"
