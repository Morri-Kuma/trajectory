#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scfm_geneformer_extract_embedding.$JOB_ID.log
#$ -l s_vmem=128G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${GENEFORMER_CONDA_ENV:-geneformer}"

GENEFORMER_MODEL_DIR="${GENEFORMER_MODEL_DIR:-/home/xzy0723/projects/trajectory/models/Geneformer-V2-104M}"
INPUT_H5AD="${INPUT_H5AD:-benchmark/inputs/representation/gse230659/source/GSE230659_full_gene_with_final_labels.h5ad}"
OUTPUT_H5AD="${OUTPUT_H5AD:-benchmark/inputs/representation/gse230659/source/GSE230659_geneformer_cls_full_gene.h5ad}"
ENSEMBL_COL="${ENSEMBL_COL:-gene_ids}"
COUNTS_COL="${COUNTS_COL:-total_counts}"
MODEL_VERSION="${MODEL_VERSION:-V2}"
EMB_MODE="${EMB_MODE:-cls}"
EMB_LAYER="${EMB_LAYER:--1}"
BATCH_SIZE="${BATCH_SIZE:-64}"
NPROC="${NPROC:-4}"
UPSTREAM_COMMIT_OR_RELEASE="${UPSTREAM_COMMIT_OR_RELEASE:-official_geneformer_v2_104m}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export GENEFORMER_MODEL_DIR
export INPUT_H5AD
export OUTPUT_H5AD
export ENSEMBL_COL
export COUNTS_COL
export MODEL_VERSION
export EMB_MODE
export EMB_LAYER
export BATCH_SIZE
export NPROC
export UPSTREAM_COMMIT_OR_RELEASE
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Extract Geneformer embedding started at: $(date)"
echo "Project root      : ${PROJECT_ROOT}"
echo "Conda env         : ${CONDA_ENV}"
echo "Model dir         : ${GENEFORMER_MODEL_DIR}"
echo "Input h5ad        : ${INPUT_H5AD}"
echo "Output h5ad       : ${OUTPUT_H5AD}"
echo "Ensembl col       : ${ENSEMBL_COL}"
echo "Counts col        : ${COUNTS_COL}"
echo "Model version     : ${MODEL_VERSION}"
echo "Emb mode          : ${EMB_MODE}"
echo "Emb layer         : ${EMB_LAYER}"
echo "Batch size        : ${BATCH_SIZE}"
echo "Nproc             : ${NPROC}"
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
if [ ! -d "${GENEFORMER_MODEL_DIR}" ]; then
  echo "ERROR: Geneformer model directory not found: ${GENEFORMER_MODEL_DIR}" >&2
  exit 2
fi

python - <<'PY'
import os
from pathlib import Path

import geneformer
from geneformer import TranscriptomeTokenizer, EmbExtractor  # noqa: F401
import anndata as ad

print("geneformer:", geneformer.__file__)
print("TranscriptomeTokenizer OK:", TranscriptomeTokenizer)
print("EmbExtractor OK:", EmbExtractor)

model_dir = Path(os.environ["GENEFORMER_MODEL_DIR"])
print("model path:", model_dir, "exists:", model_dir.exists())

in_path = Path(os.environ["INPUT_H5AD"])
adata = ad.read_h5ad(in_path, backed="r")
print("input h5ad:", in_path)
print("input shape:", adata.shape)

ensembl_col = os.environ["ENSEMBL_COL"]
counts_col = os.environ["COUNTS_COL"]
has_ensembl = ensembl_col in adata.var.columns
has_counts = counts_col in adata.obs.columns
print(f"var['{ensembl_col}'] present:", has_ensembl)
print(f"obs['{counts_col}'] present:", has_counts)
if has_ensembl:
    print(f"sample {ensembl_col}:", list(adata.var[ensembl_col].astype(str).values[:5]))
if has_counts:
    print(f"sample {counts_col}:", list(adata.obs[counts_col].values[:5]))
adata.file.close()
if not has_ensembl:
    raise SystemExit(f"input h5ad var is missing required Ensembl column '{ensembl_col}'")
if not has_counts:
    raise SystemExit(f"input h5ad obs is missing required counts column '{counts_col}'")
PY

python -m benchmark.representations.extract_scfm_embeddings \
  --model geneformer \
  --input-h5ad "${INPUT_H5AD}" \
  --output-h5ad "${OUTPUT_H5AD}" \
  --model-path "${GENEFORMER_MODEL_DIR}" \
  --ensembl-col "${ENSEMBL_COL}" \
  --counts-col "${COUNTS_COL}" \
  --emb-mode "${EMB_MODE}" \
  --emb-layer "${EMB_LAYER}" \
  --model-version "${MODEL_VERSION}" \
  --batch-size "${BATCH_SIZE}" \
  --nproc "${NPROC}" \
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
if "X_geneformer_cls" not in adata.obsm:
    raise SystemExit("X_geneformer_cls missing from output h5ad")
print("X_geneformer_cls shape:", adata.obsm["X_geneformer_cls"].shape)
adata.file.close()
PY

echo "============================================================"
echo "Extract Geneformer embedding finished at: $(date)"
echo "Output h5ad: ${OUTPUT_H5AD}"
echo "============================================================"
