#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_reeval_sctimebench_baselines.$JOB_ID.log
#$ -l s_vmem=128G

set -euo pipefail

PROJECT_ROOT="/home/xzy0723/projects/trajectory"
CONDA_SH="/home/xzy0723/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="traj_env"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"

echo "============================================================"
echo "Job started at: $(date)"
echo "Host          : $(hostname)"
echo "Project root  : ${TRAJ_PROJECT_ROOT}"
echo "Job ID        : ${JOB_ID:-N/A}"
echo "Task          : re-evaluate official WOT/CellRank2 baselines only"
echo "============================================================"

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

python -V

python scripts/reeval_lineage_sctimebench_baseline.py \
  --run benchmark/configs/wot_gse230659_observed_scgpt_v1.yaml=benchmark/results/wot/scenario_A_scgpt_v1 \
  --run benchmark/configs/wot_gse230659_observed_scgpt_v1_scenarioB.yaml=benchmark/results/wot/scenario_B_scgpt_v1 \
  --run benchmark/configs/wot_gse230659_observed_scgpt_v1_scenarioC.yaml=benchmark/results/wot/scenario_C_scgpt_v1 \
  --run benchmark/configs/cellrank2_gse230659_observed_scgpt_v1.yaml=benchmark/results/cellrank2/scenario_A_scgpt_v1 \
  --run benchmark/configs/cellrank2_gse230659_observed_scgpt_v1_scenarioB.yaml=benchmark/results/cellrank2/scenario_B_scgpt_v1 \
  --run benchmark/configs/cellrank2_gse230659_observed_scgpt_v1_scenarioC.yaml=benchmark/results/cellrank2/scenario_C_scgpt_v1

python benchmark/reports/summarize_lineage.py

echo "============================================================"
echo "Job finished at: $(date)"
echo "============================================================"
