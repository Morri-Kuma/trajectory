#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
JOB_SCRIPT_DIR="${PROJECT_ROOT}/jobs/shirokane"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
RUN_TAG="${RUN_TAG:-forecast_exact_$(date +%Y%m%d_%H%M%S)}"
TASK_RANGE="${TASK_RANGE:-1-36}"
BACKUP_OLD="${BACKUP_OLD:-1}"
RUN_SUMMARY="${RUN_SUMMARY:-1}"
RUN_TESTS="${RUN_TESTS:-1}"
RUN_OFFICIAL_SUMMARY="${RUN_OFFICIAL_SUMMARY:-1}"
REQUIRE_ALL_EXACT="${REQUIRE_ALL_EXACT:-1}"
LEIDEN_N_NEIGHBORS="${LEIDEN_N_NEIGHBORS:-15}"
LEIDEN_RESOLUTION="${LEIDEN_RESOLUTION:-0.5}"

cd "${PROJECT_ROOT}"
mkdir -p logs

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
echo "Submitting silver Forecast Accuracy exact rerun"
echo "Project root          : ${PROJECT_ROOT}"
echo "Run tag               : ${RUN_TAG}"
echo "Task range            : ${TASK_RANGE}"
echo "Tasks available       : 36 projection-capable formal outputs"
echo "Backup old forecast   : ${BACKUP_OLD}"
echo "Run summary job       : ${RUN_SUMMARY}"
echo "Run tests in summary  : ${RUN_TESTS}"
echo "Official summary      : ${RUN_OFFICIAL_SUMMARY}"
echo "Require all exact     : ${REQUIRE_ALL_EXACT}"
echo "============================================================"

FORECAST_JOB="$(
  submit_and_get_job_id \
    qsub -N "silver_fcex_${RUN_TAG}" \
      -t "${TASK_RANGE}" \
      -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}" \
      -v "CONDA_SH=${CONDA_SH}" \
      -v "CONDA_ENV=${CONDA_ENV}" \
      -v "RUN_TAG=${RUN_TAG}" \
      -v "BACKUP_OLD=${BACKUP_OLD}" \
      "${JOB_SCRIPT_DIR}/run_silver_forecast_exact_array.sh"
)"
if [ -z "${FORECAST_JOB}" ]; then
  echo "ERROR: failed to parse forecast exact array job id" >&2
  exit 2
fi
echo "Forecast exact array job id: ${FORECAST_JOB}"

if [ "${RUN_SUMMARY}" = "1" ]; then
  SUMMARY_JOB="$(
    submit_and_get_job_id \
      qsub -N "silver_fcsum_${RUN_TAG}" \
        -hold_jid "${FORECAST_JOB}" \
        -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}" \
        -v "CONDA_SH=${CONDA_SH}" \
        -v "CONDA_ENV=${CONDA_ENV}" \
        -v "RUN_TESTS=${RUN_TESTS}" \
        -v "RUN_OFFICIAL_SUMMARY=${RUN_OFFICIAL_SUMMARY}" \
        -v "REQUIRE_ALL_EXACT=${REQUIRE_ALL_EXACT}" \
        -v "LEIDEN_N_NEIGHBORS=${LEIDEN_N_NEIGHBORS}" \
        -v "LEIDEN_RESOLUTION=${LEIDEN_RESOLUTION}" \
        "${JOB_SCRIPT_DIR}/run_silver_forecast_exact_summary.sh"
  )"
  if [ -z "${SUMMARY_JOB}" ]; then
    echo "ERROR: failed to parse forecast summary job id" >&2
    exit 2
  fi
  echo "Summary job id: ${SUMMARY_JOB}"
fi

echo "============================================================"
echo "Submitted silver Forecast Accuracy exact rerun:"
echo "  ${FORECAST_JOB} run tasks ${TASK_RANGE}"
if [ "${RUN_SUMMARY}" = "1" ]; then
  echo "  ${SUMMARY_JOB} summarize after ${FORECAST_JOB}"
fi
echo
echo "Useful commands:"
echo "  qstat -u ${USER:-unknown}"
echo "  qstat -j ${FORECAST_JOB}"
echo "  ls logs/run_silver_forecast_exact.${FORECAST_JOB}.*.log"
if [ "${RUN_SUMMARY}" = "1" ]; then
  echo "  qstat -j ${SUMMARY_JOB}"
  echo "  ls logs/run_silver_forecast_exact_summary.${SUMMARY_JOB}.log"
fi
echo
echo "Override examples:"
echo "  TASK_RANGE=10-18 bash jobs/shirokane/submit_silver_forecast_exact.sh"
echo "  RUN_OFFICIAL_SUMMARY=0 bash jobs/shirokane/submit_silver_forecast_exact.sh"
echo "  BACKUP_OLD=0 REQUIRE_ALL_EXACT=0 bash jobs/shirokane/submit_silver_forecast_exact.sh"
echo "============================================================"
