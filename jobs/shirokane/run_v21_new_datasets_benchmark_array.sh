#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N v21_bench
#$ -o logs/run_v21_new_datasets_benchmark.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G
#$ -t 1-22
# ── v2.1 new-dataset benchmark array: GSE298212 + GSE218855 (A/B/C) ───────────
#
# Requires the HVG2000 inputs and official_silver providers from
# run_v21_new_datasets_preprocess.sh.
# 22 tasks = 2 datasets x (scnode/prescient/mioflow x A/B/C + wot/cellrank2 x A).
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"
RT="benchmark/configs/runtime"

CONFIGS=(
  # --- GSE298212 (tasks 1-11) ---
  "${RT}/scnode_gse298212_marker_fm_silver_A_hvg2000_formal.yaml"
  "${RT}/scnode_gse298212_marker_fm_silver_B_hvg2000_formal.yaml"
  "${RT}/scnode_gse298212_marker_fm_silver_C_hvg2000_formal.yaml"
  "${RT}/prescient_gse298212_marker_fm_silver_A_hvg2000_formal.yaml"
  "${RT}/prescient_gse298212_marker_fm_silver_B_hvg2000_formal.yaml"
  "${RT}/prescient_gse298212_marker_fm_silver_C_hvg2000_formal.yaml"
  "${RT}/mioflow_gse298212_marker_fm_silver_A_hvg2000_formal.yaml"
  "${RT}/mioflow_gse298212_marker_fm_silver_B_hvg2000_formal.yaml"
  "${RT}/mioflow_gse298212_marker_fm_silver_C_hvg2000_formal.yaml"
  "${RT}/wot_gse298212_marker_fm_silver_A_hvg2000_formal.yaml"
  "${RT}/cellrank2_gse298212_marker_fm_silver_A_hvg2000_formal.yaml"
  # --- GSE218855 (tasks 12-22) ---
  "${RT}/scnode_gse218855_marker_fm_silver_A_hvg2000_formal.yaml"
  "${RT}/scnode_gse218855_marker_fm_silver_B_hvg2000_formal.yaml"
  "${RT}/scnode_gse218855_marker_fm_silver_C_hvg2000_formal.yaml"
  "${RT}/prescient_gse218855_marker_fm_silver_A_hvg2000_formal.yaml"
  "${RT}/prescient_gse218855_marker_fm_silver_B_hvg2000_formal.yaml"
  "${RT}/prescient_gse218855_marker_fm_silver_C_hvg2000_formal.yaml"
  "${RT}/mioflow_gse218855_marker_fm_silver_A_hvg2000_formal.yaml"
  "${RT}/mioflow_gse218855_marker_fm_silver_B_hvg2000_formal.yaml"
  "${RT}/mioflow_gse218855_marker_fm_silver_C_hvg2000_formal.yaml"
  "${RT}/wot_gse218855_marker_fm_silver_A_hvg2000_formal.yaml"
  "${RT}/cellrank2_gse218855_marker_fm_silver_A_hvg2000_formal.yaml"
)

CONFIG="${CONFIGS[$((TASK_ID-1))]}"
[ -z "${CONFIG:-}" ] && { echo "ERROR: no config for task ${TASK_ID}" >&2; exit 2; }

mkdir -p "${PROJECT_ROOT}/logs"; cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"; export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"; export MKL_NUM_THREADS="${NSLOTS:-1}"; export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"

read -r METHOD SCENARIO H5AD OUTDIR PROV < <(python - "${CONFIG}" <<'PY'
import sys, yaml
from pathlib import Path
c = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8-sig")) or {}
print(str(c.get("method","")).lower(), str(c.get("scenario","")),
      str((c.get("dataset") or {}).get("h5ad_path","")), str((c.get("output") or {}).get("base_dir","")),
      str((c.get("ground_truth") or {}).get("provider_id","")))
PY
)
echo "=== task ${TASK_ID}: ${METHOD} ${SCENARIO}  ${CONFIG} at $(date) ==="
if [ ! -f "${H5AD}" ]; then
  echo "ERROR: benchmark h5ad not found: ${H5AD}" >&2
  echo "Run jobs/shirokane/run_v21_new_datasets_preprocess.sh first." >&2; exit 4
fi
GRAPH="benchmark/ground_truth/providers/${PROV}/reference_graph.json"
if [ ! -f "${GRAPH}" ]; then
  echo "ERROR: provider reference_graph not found: ${GRAPH}" >&2
  echo "Run jobs/shirokane/run_v21_new_datasets_preprocess.sh (provider step) first." >&2; exit 4
fi
mkdir -p "${OUTDIR}"
python benchmark/evaluation/eval_dispatch.py \
  --method "${METHOD}" --scenario "${SCENARIO}" \
  --adata "${H5AD}" --output-dir "${OUTDIR}" --method-config "${CONFIG}"
echo "=== task ${TASK_ID} DONE at $(date); out=${OUTDIR} ==="
