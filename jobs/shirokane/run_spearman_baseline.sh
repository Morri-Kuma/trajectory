#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N spearman_base
#$ -o logs/run_spearman_baseline.$JOB_ID.log
#$ -l s_vmem=64G
# ── scTimeBench Spearman correlation baseline for Lineage Fidelity ────────────
#
# Computes the per-dataset Spearman maximum-vote correlation baseline (the null a
# method must beat) for the four marker-silver datasets, using the frozen reference
# graphs. CPU-only; no method retraining. Writes
#   benchmark/reports/official_silver/spearman_baseline.{csv,json}
# Run on a compute node (Spearman over all cells is memory-heavy on the login node).
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
cd "${PROJECT_ROOT}"; mkdir -p logs
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"; export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"
echo "=== spearman baseline at $(date) ==="
python scripts/compute_spearman_baseline.py
echo "=== done at $(date) ==="
