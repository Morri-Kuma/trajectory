#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N method_expand
#$ -o logs/run_method_expansion_benchmark.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G
#$ -t 1-36
# ── New projection-capable methods benchmark array (scIMF, PI-SDE, Squidiff) ──
#
# 36 tasks = 3 methods x 2 chemical primaries (GSE178325, GSE230659) x 6 scenarios
# (observed A/B/C on the HVG2000 input; pseudotime D/E/F on the DPT-bin input).
#
# PREREQUISITE: each method's source must be vendored under
# benchmark/methods/{scIMF,PISDE,Squidiff}/<...>_module/ and its run.py wired (see the
# run.py docstrings). Until then a task fails fast with a clear VendoringRequired
# message and writes status=failed run_metadata.json — no partial/garbage outputs.
#
# Build configs first (idempotent): python scripts/generate_method_expansion_configs.py
# Pseudotime tasks also require run_pseudotime_axis_build.sh to have produced the
# *_pseudotime_benchmark_input.h5ad.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"
RT="benchmark/configs/runtime"

# Build the deterministic config list (method-major, dataset, scenario).
CONFIGS=()
for METHOD in scimf pisde squiddiff; do
  for DS in gse178325 gse230659; do
    for SCEN in A B C D E F; do
      if [[ "$SCEN" == [DEF] ]]; then TAG="pseudotime"; else TAG="hvg2000"; fi
      CONFIGS+=("${RT}/${METHOD}_${DS}_marker_fm_silver_${SCEN}_${TAG}_formal.yaml")
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
  echo "ERROR: config not found: ${CONFIG}. Run: python scripts/generate_method_expansion_configs.py" >&2; exit 4; fi

read -r METHOD SCENARIO H5AD OUTDIR < <(python - "${CONFIG}" <<'PY'
import sys, yaml
from pathlib import Path
c = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8-sig")) or {}
print(str(c.get("method","")).lower(), str(c.get("scenario","")),
      str((c.get("dataset") or {}).get("h5ad_path","")), str((c.get("output") or {}).get("base_dir","")))
PY
)
echo "=== task ${TASK_ID}: ${METHOD} ${SCENARIO}  ${CONFIG} at $(date) ==="
if [ ! -f "${H5AD}" ]; then echo "ERROR: input h5ad not found: ${H5AD}" >&2; exit 4; fi
mkdir -p "${OUTDIR}"
python benchmark/evaluation/eval_dispatch.py \
  --method "${METHOD}" --scenario "${SCENARIO}" \
  --adata "${H5AD}" --output-dir "${OUTDIR}" --method-config "${CONFIG}"
echo "=== task ${TASK_ID} DONE at $(date); out=${OUTDIR} ==="
