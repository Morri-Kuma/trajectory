#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N pt_bench
#$ -o logs/run_pseudotime_benchmark.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G
#$ -t 1-22
# ── Pseudotime Scenarios D/E/F benchmark array (GSE230659 + GSE178325) ────────
#
# Requires the pseudotime inputs from run_pseudotime_axis_build.sh.
# 22 tasks = 2 datasets x (scnode/prescient/mioflow x D/E/F + wot/cellrank2 x D).
# WOT/CellRank2 are lineage-only (Scenario D only); projection methods run D/E/F
# across all three core dimensions.  GPU is used by tasks where use_cuda: true;
# the runners fall back to CPU automatically.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"
RT="benchmark/configs/runtime"

CONFIGS=(
  # --- GSE230659 (tasks 1-11) ---
  "${RT}/scnode_gse230659_marker_fm_silver_D_pseudotime_formal.yaml"
  "${RT}/scnode_gse230659_marker_fm_silver_E_pseudotime_formal.yaml"
  "${RT}/scnode_gse230659_marker_fm_silver_F_pseudotime_formal.yaml"
  "${RT}/prescient_gse230659_marker_fm_silver_D_pseudotime_formal.yaml"
  "${RT}/prescient_gse230659_marker_fm_silver_E_pseudotime_formal.yaml"
  "${RT}/prescient_gse230659_marker_fm_silver_F_pseudotime_formal.yaml"
  "${RT}/mioflow_gse230659_marker_fm_silver_D_pseudotime_formal.yaml"
  "${RT}/mioflow_gse230659_marker_fm_silver_E_pseudotime_formal.yaml"
  "${RT}/mioflow_gse230659_marker_fm_silver_F_pseudotime_formal.yaml"
  "${RT}/wot_gse230659_marker_fm_silver_D_pseudotime_formal.yaml"
  "${RT}/cellrank2_gse230659_marker_fm_silver_D_pseudotime_formal.yaml"
  # --- GSE178325 (tasks 12-22) ---
  "${RT}/scnode_gse178325_marker_fm_silver_D_pseudotime_formal.yaml"
  "${RT}/scnode_gse178325_marker_fm_silver_E_pseudotime_formal.yaml"
  "${RT}/scnode_gse178325_marker_fm_silver_F_pseudotime_formal.yaml"
  "${RT}/prescient_gse178325_marker_fm_silver_D_pseudotime_formal.yaml"
  "${RT}/prescient_gse178325_marker_fm_silver_E_pseudotime_formal.yaml"
  "${RT}/prescient_gse178325_marker_fm_silver_F_pseudotime_formal.yaml"
  "${RT}/mioflow_gse178325_marker_fm_silver_D_pseudotime_formal.yaml"
  "${RT}/mioflow_gse178325_marker_fm_silver_E_pseudotime_formal.yaml"
  "${RT}/mioflow_gse178325_marker_fm_silver_F_pseudotime_formal.yaml"
  "${RT}/wot_gse178325_marker_fm_silver_D_pseudotime_formal.yaml"
  "${RT}/cellrank2_gse178325_marker_fm_silver_D_pseudotime_formal.yaml"
)

CONFIG="${CONFIGS[$((TASK_ID-1))]}"
[ -z "${CONFIG:-}" ] && { echo "ERROR: no config for task ${TASK_ID}" >&2; exit 2; }

mkdir -p "${PROJECT_ROOT}/logs"; cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"; export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"; export MKL_NUM_THREADS="${NSLOTS:-1}"; export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"

read -r METHOD SCENARIO H5AD OUTDIR < <(python - "${CONFIG}" <<'PY'
import sys, yaml
from pathlib import Path
c = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8-sig")) or {}
print(str(c.get("method","")).lower(), str(c.get("scenario","")),
      str((c.get("dataset") or {}).get("h5ad_path","")), str((c.get("output") or {}).get("base_dir","")))
PY
)
echo "=== task ${TASK_ID}: ${METHOD} ${SCENARIO}  ${CONFIG} at $(date) ==="
if [ ! -f "${H5AD}" ]; then
  echo "ERROR: pseudotime h5ad not found: ${H5AD}" >&2
  echo "Run jobs/shirokane/run_pseudotime_axis_build.sh first." >&2; exit 4
fi
mkdir -p "${OUTDIR}"
python benchmark/evaluation/eval_dispatch.py \
  --method "${METHOD}" --scenario "${SCENARIO}" \
  --adata "${H5AD}" --output-dir "${OUTDIR}" --method-config "${CONFIG}"
echo "=== task ${TASK_ID} DONE at $(date); out=${OUTDIR} ==="
