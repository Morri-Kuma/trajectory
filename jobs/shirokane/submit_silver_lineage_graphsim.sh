#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
JOB_SCRIPT_DIR="${PROJECT_ROOT}/jobs/shirokane"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
RUN_TAG="${RUN_TAG:-lineage_graphsim_$(date +%Y%m%d_%H%M%S)}"
TASK_RANGE="${TASK_RANGE:-1-40}"
RUN_BASELINE="${RUN_BASELINE:-1}"
BACKUP_OLD="${BACKUP_OLD:-1}"
RUN_SUMMARY="${RUN_SUMMARY:-1}"
RUN_TESTS="${RUN_TESTS:-1}"
RUN_OFFICIAL_SUMMARY="${RUN_OFFICIAL_SUMMARY:-1}"
REQUIRE_ALL_GRAPH_SIM="${REQUIRE_ALL_GRAPH_SIM:-1}"
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
echo "Submitting silver lineage scTimeBench graph-sim rerun"
echo "Project root          : ${PROJECT_ROOT}"
echo "Run tag               : ${RUN_TAG}"
echo "Task range            : ${TASK_RANGE}"
echo "Tasks available       : 40 formal silver lineage outputs"
echo "Run correlation base  : ${RUN_BASELINE}"
echo "Backup old metrics    : ${BACKUP_OLD}"
echo "Run summary job       : ${RUN_SUMMARY}"
echo "Run tests in summary  : ${RUN_TESTS}"
echo "Official summary      : ${RUN_OFFICIAL_SUMMARY}"
echo "Require all graph-sim : ${REQUIRE_ALL_GRAPH_SIM}"
echo "============================================================"

LINEAGE_JOB="$(
  submit_and_get_job_id \
    qsub -N "silver_lingraph_${RUN_TAG}" \
      -t "${TASK_RANGE}" \
      -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}" \
      -v "CONDA_SH=${CONDA_SH}" \
      -v "CONDA_ENV=${CONDA_ENV}" \
      -v "RUN_TAG=${RUN_TAG}" \
      -v "RUN_BASELINE=${RUN_BASELINE}" \
      -v "BACKUP_OLD=${BACKUP_OLD}" \
      "${JOB_SCRIPT_DIR}/run_silver_lineage_graphsim_array.sh"
)"
if [ -z "${LINEAGE_JOB}" ]; then
  echo "ERROR: failed to parse lineage graph-sim array job id" >&2
  exit 2
fi
echo "Lineage graph-sim array job id: ${LINEAGE_JOB}"

if [ "${RUN_SUMMARY}" = "1" ]; then
  SUMMARY_JOB="$(
    submit_and_get_job_id \
      qsub -N "silver_linsum_${RUN_TAG}" \
        -hold_jid "${LINEAGE_JOB}" \
        -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}" \
        -v "CONDA_SH=${CONDA_SH}" \
        -v "CONDA_ENV=${CONDA_ENV}" \
        -v "RUN_TESTS=${RUN_TESTS}" \
        -v "RUN_OFFICIAL_SUMMARY=${RUN_OFFICIAL_SUMMARY}" \
        -v "REQUIRE_ALL_GRAPH_SIM=${REQUIRE_ALL_GRAPH_SIM}" \
        -v "LEIDEN_N_NEIGHBORS=${LEIDEN_N_NEIGHBORS}" \
        -v "LEIDEN_RESOLUTION=${LEIDEN_RESOLUTION}" \
        "${JOB_SCRIPT_DIR}/run_silver_lineage_graphsim_summary.sh"
  )"
  if [ -z "${SUMMARY_JOB}" ]; then
    echo "ERROR: failed to parse summary job id" >&2
    exit 2
  fi
  echo "Summary job id: ${SUMMARY_JOB}"
fi

echo "============================================================"
echo "Submitted silver lineage graph-sim rerun:"
echo "  ${LINEAGE_JOB} run tasks ${TASK_RANGE}"
if [ "${RUN_SUMMARY}" = "1" ]; then
  echo "  ${SUMMARY_JOB} summarize after ${LINEAGE_JOB}"
fi
echo
echo "Useful commands:"
echo "  qstat -u ${USER:-unknown}"
echo "  qstat -j ${LINEAGE_JOB}"
echo "  ls logs/run_silver_lineage_graphsim.${LINEAGE_JOB}.*.log"
if [ "${RUN_SUMMARY}" = "1" ]; then
  echo "  qstat -j ${SUMMARY_JOB}"
  echo "  ls logs/run_silver_lineage_graphsim_summary.${SUMMARY_JOB}.log"
fi
echo
echo "Override examples:"
echo "  RUN_BASELINE=0 bash jobs/shirokane/submit_silver_lineage_graphsim.sh"
echo "  TASK_RANGE=12-31 REQUIRE_ALL_GRAPH_SIM=0 bash jobs/shirokane/submit_silver_lineage_graphsim.sh"
echo "============================================================"
