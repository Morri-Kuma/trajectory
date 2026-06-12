#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N multiseed
#$ -o logs/run_multiseed_benchmark.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G
#$ -t 1-90
# ── Multi-seed reruns for per-method confidence intervals ─────────────────────
#
# 90 tasks = 3 methods (scnode, prescient, mioflow) x 2 primaries (GSE178325,
# GSE230659) x 3 observed-time scenarios (A/B/C) x 5 seeds (101/202/303/404/505).
# Each run writes to benchmark/results/<method>/<ds>_marker_fm_silver_<sc>_hvg2000_seed<NN>.
#
# PREREQUISITE: build the seed configs first (idempotent):
#   python scripts/generate_multiseed_configs.py
# AFTER the array completes, aggregate to mean+-CI:
#   python scripts/aggregate_multiseed_ci.py
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"
RT="benchmark/configs/runtime/multiseed"

# Deterministic config list (method-major, dataset, scenario, seed) — must match
# the ordering implied by scripts/generate_multiseed_configs.py.
CONFIGS=()
for METHOD in scnode prescient mioflow; do
  for DS in gse178325 gse230659; do
    for SCEN in A B C; do
      for SEED in 101 202 303 404 505; do
        CONFIGS+=("${RT}/${METHOD}_${DS}_marker_fm_silver_${SCEN}_hvg2000_seed${SEED}.yaml")
      done
    done
  done
done

CONFIG="${CONFIGS[$((TASK_ID-1))]}"
[ -z "${CONFIG:-}" ] && { echo "ERROR: no config for task ${TASK_ID}" >&2; exit 2; }

mkdir -p "${PROJECT_ROOT}/logs"; cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"; export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"; export MKL_NUM_THREADS="${NSLOTS:-1}"; export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"

if [ ! -f "${CONFIG}" ]; then
  echo "ERROR: config not found: ${CONFIG}. Run: python scripts/generate_multiseed_configs.py" >&2; exit 4; fi

read -r METHOD SCENARIO H5AD OUTDIR < <(python - "${CONFIG}" <<'PY'
import sys, yaml
from pathlib import Path
c = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8-sig")) or {}
print(str(c.get("method","")).lower(), str(c.get("scenario","")),
      str((c.get("dataset") or {}).get("h5ad_path","")), str((c.get("output") or {}).get("base_dir","")))
PY
)
echo "=== multiseed task ${TASK_ID}: ${METHOD} ${SCENARIO}  ${CONFIG} at $(date) ==="
if [ ! -f "${H5AD}" ]; then echo "ERROR: input h5ad not found: ${H5AD}" >&2; exit 4; fi
mkdir -p "${OUTDIR}"
python benchmark/evaluation/eval_dispatch.py \
  --method "${METHOD}" --scenario "${SCENARIO}" \
  --adata "${H5AD}" --output-dir "${OUTDIR}" --method-config "${CONFIG}"
echo "=== multiseed task ${TASK_ID} DONE at $(date); out=${OUTDIR} ==="
