#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_cellrank_array.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G
#$ -t 1-2

set -euo pipefail

PROJECT_ROOT="/home/xzy0723/projects/trajectory"
SCGPT_REPO="/home/xzy0723/projects/scGPT"
CONDA_SH="/home/xzy0723/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="traj_env"

ADATA="benchmark/results/scgpt/full/adata_scgpt_annotated.h5ad"
METHOD="cellrank2"

case "${SGE_TASK_ID}" in
  1)
    SCENARIO="B"
    METHOD_CONFIG="benchmark/configs/cellrank2_gse230659_observed_scgpt_v1_scenarioB.yaml"
    OUTPUT_DIR="benchmark/results/cellrank2/scenario_B_scgpt_v1_shirokane_test"
    ;;
  2)
    SCENARIO="C"
    METHOD_CONFIG="benchmark/configs/cellrank2_gse230659_observed_scgpt_v1_scenarioC.yaml"
    OUTPUT_DIR="benchmark/results/cellrank2/scenario_C_scgpt_v1_shirokane_test"
    ;;
  *)
    echo "Unsupported SGE_TASK_ID: ${SGE_TASK_ID}"
    exit 1
    ;;
esac

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export SCGPT_REPO="${SCGPT_REPO}"

echo "============================================================"
echo "Job started at: $(date)"
echo "Host          : $(hostname)"
echo "Job ID        : ${JOB_ID:-N/A}"
echo "Task ID       : ${SGE_TASK_ID:-N/A}"
echo "Project root  : ${TRAJ_PROJECT_ROOT}"
echo "SCGPT repo    : ${SCGPT_REPO}"
echo "Method        : ${METHOD}"
echo "Scenario      : ${SCENARIO}"
echo "ADATA         : ${ADATA}"
echo "Method config : ${METHOD_CONFIG}"
echo "Output dir    : ${OUTPUT_DIR}"
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
  python -m benchmark.evaluation.eval_dispatch
  --method "${METHOD}"
  --scenario "${SCENARIO}"
  --adata "${ADATA}"
  --output-dir "${OUTPUT_DIR}"
  --method-config "${METHOD_CONFIG}"
)

echo "------------------------------------------------------------"
echo "Running command:"
printf ' %q' "${CMD[@]}"
echo
echo "------------------------------------------------------------"

time "${CMD[@]}"

echo "============================================================"
echo "Job finished at: $(date)"
echo "============================================================"