#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_wot_validate_C.$JOB_ID.log
#$ -l s_vmem=32G

set -euo pipefail

# ===== user-configurable section =====
PROJECT_ROOT="/home/xzy0723/projects/trajectory"
SCGPT_REPO="/home/xzy0723/projects/scGPT"
CONDA_SH="/home/xzy0723/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="traj_env"

CONFIG="benchmark/configs/wot_gse230659_observed_scgpt_v1_scenarioB.yaml"
OUTPUT_DIR="benchmark/results/wot/scenario_B_scgpt_v1_shirokane_test"

# For a quick smoke test, set this to a number like 500.
# For full validation, leave it empty.
SUBSAMPLE=""

# Reuse existing tmaps on reruns
SKIP_TMAP_IF_EXISTS="1"
# =====================================

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export SCGPT_REPO="${SCGPT_REPO}"

echo "============================================================"
echo "Job started at: $(date)"
echo "Host          : $(hostname)"
echo "Project root  : ${TRAJ_PROJECT_ROOT}"
echo "SCGPT repo    : ${SCGPT_REPO}"
echo "Config        : ${CONFIG}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Job ID        : ${JOB_ID:-N/A}"
echo "============================================================"

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

echo "Python        : $(which python)"
python -V

python - <<'PY'
import anndata, scanpy, cellrank, scvelo, wot
print("imports OK")
print("anndata :", anndata.__version__)
print("scanpy  :", scanpy.__version__)
print("cellrank:", cellrank.__version__)
print("scvelo  :", scvelo.__version__)
PY

CMD=(
  python benchmark/methods/WOT/run.py
  --config "${CONFIG}"
  --output-dir "${OUTPUT_DIR}"
)

if [[ -n "${SUBSAMPLE}" ]]; then
  CMD+=(--subsample "${SUBSAMPLE}")
fi

if [[ "${SKIP_TMAP_IF_EXISTS}" == "1" ]]; then
  CMD+=(--skip-tmap-if-exists)
fi

echo "------------------------------------------------------------"
echo "Running command:"
printf ' %q' "${CMD[@]}"
echo
echo "------------------------------------------------------------"

time "${CMD[@]}"

echo "============================================================"
echo "Job finished at: $(date)"
echo "============================================================"