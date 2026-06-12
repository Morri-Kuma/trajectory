#!/bin/bash
# ── GSE175634 Shirokane submission script ─────────────────────────────────────
#
# DATASET ID NOTE:
#   Canonical GEO accession: GSE175634 (not GSE174534 — that was a typo).
#   All source files are named GSE175634_*.
#
# This script submits the full GSE175634 benchmark workflow as two sequential
# SGE jobs:
#
#   Job 1 — Preprocessing (run_gse175634_preprocess.sh):
#     - Builds the HVG2000 benchmark input h5ad from raw MTX files.
#       Silver-standard annotation step is SKIPPED (author labels used directly).
#     - Builds the frozen ground-truth provider (gse175634_cardiac_silver_v1).
#
#   Job 2 — Benchmark array (run_gse175634_benchmark_array.sh):
#     - 11 array tasks: scnode/prescient/mioflow × scenarios A/B/C, plus
#       wot and cellrank2 for scenario A.
#     - Held back with -hold_jid until Job 1 completes successfully.
#
# Usage on Shirokane:
#   cd /home/xzy0723/projects/trajectory
#   bash jobs/shirokane/submit_gse175634_benchmark.sh
#
# To adjust paths, set environment variables before running:
#   export TRAJ_PROJECT_ROOT=/your/path/to/trajectory
#   export CONDA_SH=/your/path/to/miniconda3/etc/profile.d/conda.sh
#   export CONDA_ENV=traj_env
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
JOB_SCRIPT_DIR="${PROJECT_ROOT}/jobs/shirokane"
RUN_TAG="${RUN_TAG:-gse175634_$(date +%Y%m%d_%H%M%S)}"

cd "${PROJECT_ROOT}"
mkdir -p logs

# ── helper: extract numeric job id from qsub output ──────────────────────────
submit_and_get_job_id() {
  local output
  output="$("$@" 2>&1)"
  echo "${output}" >&2
  echo "${output}" | awk '
    /Your job/ || /Your job-array/ {
      for (i = 1; i <= NF; i++) {
        token = $i
        if (token ~ /^[0-9]+(\.[0-9:,-]+)?$/) {
          sub(/\..*/, "", token)
          print token
          exit
        }
      }
    }
  '
}

echo "============================================================"
echo "Submitting GSE175634 benchmark workflow"
echo "Project root : ${PROJECT_ROOT}"
echo "Run tag      : ${RUN_TAG}"
echo "Silver step  : SKIPPED (author 'type' labels used directly)"
echo "============================================================"

# ── Job 1: preprocessing ──────────────────────────────────────────────────────
echo ""
echo "Submitting preprocessing job …"
PREPROCESS_JOB="$(
  submit_and_get_job_id \
    qsub -N gse175634_preprocess_${RUN_TAG} \
         -l s_vmem=64G \
         "${JOB_SCRIPT_DIR}/run_gse175634_preprocess.sh"
)"
if [ -z "${PREPROCESS_JOB}" ]; then
  echo "ERROR: failed to parse preprocessing job id" >&2
  exit 2
fi
echo "Preprocessing job id: ${PREPROCESS_JOB}"

# ── Job 2: benchmark array (held until preprocessing completes) ───────────────
echo ""
echo "Submitting benchmark array job (held on ${PREPROCESS_JOB}) …"
BENCHMARK_JOB="$(
  submit_and_get_job_id \
    qsub -N gse175634_benchmark_${RUN_TAG} \
         -hold_jid "${PREPROCESS_JOB}" \
         "${JOB_SCRIPT_DIR}/run_gse175634_benchmark_array.sh"
)"
if [ -z "${BENCHMARK_JOB}" ]; then
  echo "ERROR: failed to parse benchmark array job id" >&2
  exit 2
fi
echo "Benchmark array job id: ${BENCHMARK_JOB}"

echo ""
echo "============================================================"
echo "Submitted GSE175634 workflow:"
echo "  Preprocessing : ${PREPROCESS_JOB}"
echo "  Benchmark     : ${BENCHMARK_JOB}  (11 tasks; held on ${PREPROCESS_JOB})"
echo ""
echo "Useful commands:"
echo "  qstat -u ${USER}"
echo "  qstat -j ${PREPROCESS_JOB}"
echo "  qstat -j ${BENCHMARK_JOB}"
echo "  ls logs/run_gse175634_preprocess.${PREPROCESS_JOB}.log"
echo "  ls logs/run_gse175634_benchmark.${BENCHMARK_JOB}.*.log"
echo ""
echo "Results will be written to:"
echo "  benchmark/results/{scnode,prescient,mioflow,wot,cellrank2}/gse175634_cardiac_silver_*"
echo "============================================================"
