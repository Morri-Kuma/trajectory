#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N pt_axis_build
#$ -o logs/run_pseudotime_axis_build.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G
#$ -t 1-2
# ── Pseudotime axis builder for Scenarios D/E/F ───────────────────────────────
#
# Builds the DPT pseudotime benchmark input for the two chemical primaries by
# computing a Scanpy DPT axis on the existing HVG2000 benchmark h5ad and
# discretizing it into 15 equal-cell-count bins. Writes:
#   *_HVG2000_pseudotime_benchmark_input.h5ad   (obs['dpt_pseudotime_bin_numeric'])
# which the D/E/F runtime configs consume (time_key: dpt_pseudotime_bin_numeric).
#
# Task 1 : GSE230659    Task 2 : GSE178325
# Run this BEFORE run_pseudotime_benchmark_array.sh.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"

case "${TASK_ID}" in
  1) IN="benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"
     OUT="benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_marker_fm_transition_silver_HVG2000_pseudotime_benchmark_input.h5ad" ;;
  2) IN="benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"
     OUT="benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_marker_fm_transition_silver_HVG2000_pseudotime_benchmark_input.h5ad" ;;
  *) echo "ERROR: task id ${TASK_ID} not in 1-2" >&2; exit 2 ;;
esac

mkdir -p "${PROJECT_ROOT}/logs"; cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"; export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"

echo "=== pseudotime axis build (task ${TASK_ID}) at $(date) ==="
if [ ! -f "${IN}" ]; then echo "ERROR: input not found: ${IN}" >&2; exit 4; fi
python scripts/build_pseudotime_axis.py --input "${IN}" --output "${OUT}" --time-key abs_day --num-bins 15
if [ ! -f "${OUT}" ]; then echo "ERROR: output not produced: ${OUT}" >&2; exit 3; fi
echo "=== DONE: ${OUT} at $(date) ==="
