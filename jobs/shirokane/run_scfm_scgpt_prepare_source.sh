#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scfm_scgpt_prepare_source.$JOB_ID.log
#$ -l s_vmem=96G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"

FULL_H5AD="${FULL_H5AD:-benchmark/inputs/gse230659_marker_fm_transition_silver_full_gene/GSE230659_stage1_marker_seed_full_gene.h5ad}"
HVG_LABEL_H5AD="${HVG_LABEL_H5AD:-benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad}"
OUTPUT_H5AD="${OUTPUT_H5AD:-benchmark/inputs/representation/gse230659/source/GSE230659_full_gene_with_final_labels.h5ad}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export FULL_H5AD
export HVG_LABEL_H5AD
export OUTPUT_H5AD
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Prepare scGPT source h5ad started at: $(date)"
echo "Project root       : ${PROJECT_ROOT}"
echo "Conda env          : ${CONDA_ENV}"
echo "Full-gene h5ad     : ${FULL_H5AD}"
echo "Label source h5ad  : ${HVG_LABEL_H5AD}"
echo "Output h5ad        : ${OUTPUT_H5AD}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

python - <<'PY'
import os
from pathlib import Path

import anndata as ad

root = Path(os.environ["TRAJ_PROJECT_ROOT"])
full_path = root / os.environ["FULL_H5AD"]
hvg_path = root / os.environ["HVG_LABEL_H5AD"]
out_path = root / os.environ["OUTPUT_H5AD"]

for path in (full_path, hvg_path):
    if not path.exists():
        raise SystemExit(f"missing required h5ad: {path}")

out_path.parent.mkdir(parents=True, exist_ok=True)

full = ad.read_h5ad(full_path)
hvg = ad.read_h5ad(hvg_path, backed="r")

if list(full.obs_names) != list(hvg.obs_names):
    raise SystemExit("obs_names differ between full-gene and label h5ad")

cols = [
    "trajectory_membership_label",
    "trajectory_membership_score",
    "trajectory_membership_margin",
    "stage2_label_coarse",
    "stage2_label_expanded",
    "stage2_confidence",
    "stage2_source",
    "final_milestone_label_coarse",
    "final_milestone_label_expanded",
    "final_milestone_confidence",
    "final_milestone_source",
]

missing = [c for c in cols if c not in hvg.obs.columns]
if missing:
    raise SystemExit(f"label source h5ad is missing columns: {missing}")

for col in cols:
    full.obs[col] = hvg.obs[col].astype(str).values

hvg.file.close()
full.write_h5ad(out_path)

print(f"wrote {out_path} shape={full.shape}")
print("has final_milestone_label_coarse:", "final_milestone_label_coarse" in full.obs.columns)
print("obs_names_first:", list(full.obs_names[:3]))
PY

echo "============================================================"
echo "Prepare scGPT source h5ad finished at: $(date)"
echo "Output h5ad: ${OUTPUT_H5AD}"
echo "============================================================"
