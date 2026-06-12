#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scfm_scgpt_extract_embedding.$JOB_ID.log
#$ -l s_vmem=128G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${SCGPT_CONDA_ENV:-scgpt}"

SCGPT_MODEL_DIR="${SCGPT_MODEL_DIR:-/home/xzy0723/projects/trajectory/models/scgpt_whole_human}"
INPUT_H5AD="${INPUT_H5AD:-benchmark/inputs/representation/gse230659/source/GSE230659_full_gene_with_final_labels.h5ad}"
OUTPUT_H5AD="${OUTPUT_H5AD:-benchmark/inputs/representation/gse230659/source/GSE230659_scgpt_cls_full_gene.h5ad}"
GENE_COL="${GENE_COL:-gene_name}"
BATCH_SIZE="${BATCH_SIZE:-64}"
DEVICE="${DEVICE:-cuda}"
REQUIRE_CUDA="${REQUIRE_CUDA:-1}"
UPSTREAM_COMMIT_OR_RELEASE="${UPSTREAM_COMMIT_OR_RELEASE:-official_scgpt_local_checkout}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export SCGPT_MODEL_DIR
export INPUT_H5AD
export OUTPUT_H5AD
export GENE_COL
export BATCH_SIZE
export DEVICE
export REQUIRE_CUDA
export UPSTREAM_COMMIT_OR_RELEASE
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Extract scGPT embedding started at: $(date)"
echo "Project root      : ${PROJECT_ROOT}"
echo "Conda env         : ${CONDA_ENV}"
echo "Model dir         : ${SCGPT_MODEL_DIR}"
echo "Input h5ad        : ${INPUT_H5AD}"
echo "Output h5ad       : ${OUTPUT_H5AD}"
echo "Gene col          : ${GENE_COL}"
echo "Batch size        : ${BATCH_SIZE}"
echo "Device            : ${DEVICE}"
echo "Require CUDA      : ${REQUIRE_CUDA}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

for f in vocab.json args.json best_model.pt; do
  if [ ! -f "${SCGPT_MODEL_DIR}/${f}" ]; then
    echo "ERROR: missing ${SCGPT_MODEL_DIR}/${f}" >&2
    exit 2
  fi
done
if [ ! -f "${INPUT_H5AD}" ]; then
  echo "ERROR: input h5ad not found: ${INPUT_H5AD}" >&2
  exit 2
fi

python - <<'PY'
import scgpt
import os
import torch
from scgpt.tasks.cell_emb import embed_data
print("scgpt:", scgpt.__file__)
print("embed_data OK:", embed_data)
print("torch:", torch.__version__)
print("cuda available:", torch.cuda.is_available())
if os.environ.get("DEVICE") == "cuda" and os.environ.get("REQUIRE_CUDA", "1") == "1":
    if not torch.cuda.is_available():
        raise SystemExit(
            "DEVICE=cuda but torch.cuda.is_available() is false. "
            "The job did not receive a CUDA-visible GPU. Fix qsub GPU resources "
            "or rerun with DEVICE=cpu REQUIRE_CUDA=0 for a slow CPU extraction."
        )
PY

python -m benchmark.representations.extract_scfm_embeddings \
  --model scgpt \
  --input-h5ad "${INPUT_H5AD}" \
  --output-h5ad "${OUTPUT_H5AD}" \
  --model-path "${SCGPT_MODEL_DIR}" \
  --gene-col "${GENE_COL}" \
  --batch-size "${BATCH_SIZE}" \
  --device "${DEVICE}" \
  --upstream-commit-or-release "${UPSTREAM_COMMIT_OR_RELEASE}"

python - <<'PY'
import os
from pathlib import Path

import anndata as ad

path = Path(os.environ["OUTPUT_H5AD"])
adata = ad.read_h5ad(path, backed="r")
print("output:", path)
print("shape:", adata.shape)
print("obsm keys:", list(adata.obsm.keys()))
if "X_scgpt_cls" not in adata.obsm:
    raise SystemExit("X_scgpt_cls missing from output h5ad")
print("X_scgpt_cls shape:", adata.obsm["X_scgpt_cls"].shape)
adata.file.close()
PY

echo "============================================================"
echo "Extract scGPT embedding finished at: $(date)"
echo "Output h5ad: ${OUTPUT_H5AD}"
echo "============================================================"
