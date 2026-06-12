#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_prep_reconstruct.$JOB_ID.log
#$ -l s_vmem=96G
set -euo pipefail
source /home/xzy0723/miniconda3/etc/profile.d/conda.sh
conda activate traj_env
mkdir -p logs results/annotation_branch/inputs
python scripts/prep_scanvi_inputs_reconstruct.py
