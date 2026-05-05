#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_build_gse178325_scgpt_provider.$JOB_ID.log
#$ -q gpuh.q
#$ -l h100=1
#$ -l s_vmem=128G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-scgpt}"
SCGPT_REPO="${SCGPT_REPO:-/home/xzy0723/projects/scGPT}"

SCRIPT="benchmark/scgpt/build_gse178325_0618_provider.py"
INPUT_H5AD="data/processed/gse178325_human/GSE178325_0618_raw_full_gene_benchmark_input.h5ad"
PROVIDER_DIR="benchmark/ground_truth/providers/scgpt_v1_gse178325_0618"
OUTPUT_DIR="benchmark/results/scgpt/gse178325_0618/full"
HVG_SCRIPT="benchmark/scgpt/build_gse178325_hvg2000_from_scgpt_full.py"
HVG_OUTPUT="benchmark/inputs/gse178325_scgpt_hvg2000/GSE178325_scGPT_annotated_HVG2000_benchmark_input.h5ad"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export SCGPT_REPO="${SCGPT_REPO}"
export PYTHONPATH="${PROJECT_ROOT}:${SCGPT_REPO}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Job started at : $(date)"
echo "Host           : $(hostname)"
echo "Project root   : ${PROJECT_ROOT}"
echo "scGPT repo     : ${SCGPT_REPO}"
echo "Script         : ${SCRIPT}"
echo "Input h5ad     : ${INPUT_H5AD}"
echo "Provider dir   : ${PROVIDER_DIR}"
echo "Output dir     : ${OUTPUT_DIR}"
echo "HVG output     : ${HVG_OUTPUT}"
echo "Memory request : s_vmem=128G"
echo "Queue          : gpuh.q"
echo "GPU request    : h100=1"
echo "Job ID         : ${JOB_ID:-N/A}"
echo "============================================================"

if [ ! -f "${SCRIPT}" ]; then
  echo "ERROR: script not found: ${SCRIPT}" >&2
  exit 2
fi

if [ ! -f "${INPUT_H5AD}" ]; then
  echo "ERROR: input h5ad not found: ${INPUT_H5AD}" >&2
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
import scanpy
import torch
import scgpt
print("imports OK")
print("anndata:", anndata.__version__)
print("scanpy:", scanpy.__version__)
print("torch:", torch.__version__)
print("cuda_available:", torch.cuda.is_available())
print("cuda_device_count:", torch.cuda.device_count())
if torch.cuda.is_available():
    print("cuda_device_name:", torch.cuda.get_device_name(0))
print("scgpt:", getattr(scgpt, "__version__", "unknown"))
if not torch.cuda.is_available():
    raise SystemExit("ERROR: CUDA unavailable; this script must run on gpuh.q/H100.")
PY

echo "------------------------------------------------------------"
echo "Building GSE178325 0618 scGPT-derived provider"
echo "------------------------------------------------------------"
time python "${SCRIPT}" \
  --input-h5ad "${INPUT_H5AD}" \
  --output-dir "${OUTPUT_DIR}" \
  --provider-dir "${PROVIDER_DIR}" \
  --require-cuda \
  --num-workers 0 \
  --skip-embedding-if-exists

echo "------------------------------------------------------------"
echo "Exporting GSE178325 scGPT-annotated HVG2000 benchmark input"
echo "------------------------------------------------------------"
time python "${HVG_SCRIPT}" \
  --input-h5ad "${OUTPUT_DIR}/adata_scgpt_annotated.h5ad" \
  --output-h5ad "${HVG_OUTPUT}" \
  --n-top-genes 2000

echo "------------------------------------------------------------"
echo "Provider validation"
echo "------------------------------------------------------------"
python - <<'PY'
from pathlib import Path
from benchmark.ground_truth import load_ground_truth
from benchmark.evaluation.eval_lineage import load_reference_graph

spec = load_ground_truth(provider_id="scgpt_v1_gse178325_0618")
print(spec)
required = [
    spec.label_path,
    spec.metadata_path,
    spec.reference_graph_path,
    spec.reference_edges_path,
]
missing = [p for p in required if not p or not Path(p).exists()]
print("missing:", missing)
if missing:
    raise SystemExit(1)
mat, edges, nodes = load_reference_graph(
    spec.reference_graph_path,
    edge_confidence_mode=spec.confidence_mode,
    exclude_uncertain_states=spec.exclude_uncertain_states,
)
print("nodes:", len(nodes))
print("edges used:", len(edges))
if len(nodes) == 0 or len(edges) == 0:
    raise SystemExit("empty provider graph")
PY

echo "============================================================"
echo "Job finished at : $(date)"
echo "============================================================"
