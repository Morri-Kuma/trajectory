#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
JOB_SCRIPT_DIR="${PROJECT_ROOT}/jobs/shirokane"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
RUN_TAG="${RUN_TAG:-marker_fm_coarse_$(date +%Y%m%d_%H%M%S)}"
TASK_RANGE="${TASK_RANGE:-1-20}"
RUN_PROVIDER="${RUN_PROVIDER:-1}"
RUN_SUMMARY="${RUN_SUMMARY:-1}"
RUN_OFFICIAL_SUMMARY="${RUN_OFFICIAL_SUMMARY:-1}"
BACKUP_PROVIDER_OLD="${BACKUP_PROVIDER_OLD:-1}"
BACKUP_OLD_LINEAGE="${BACKUP_OLD_LINEAGE:-1}"
PURGE_METHOD_CACHE="${PURGE_METHOD_CACHE:-0}"
PROVIDER_S_VMEM="${PROVIDER_S_VMEM:-32G}"
METHOD_S_VMEM="${METHOD_S_VMEM:-128G}"
SUMMARY_S_VMEM="${SUMMARY_S_VMEM:-32G}"
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
echo "Submitting marker FM coarse-label rerun"
echo "Project root          : ${PROJECT_ROOT}"
echo "Run tag               : ${RUN_TAG}"
echo "Task range            : ${TASK_RANGE}"
echo "Run provider refresh  : ${RUN_PROVIDER}"
echo "Run summary           : ${RUN_SUMMARY}"
echo "Official summary      : ${RUN_OFFICIAL_SUMMARY}"
echo "Backup old provider   : ${BACKUP_PROVIDER_OLD}"
echo "Backup old lineage    : ${BACKUP_OLD_LINEAGE}"
echo "Purge method cache    : ${PURGE_METHOD_CACHE}"
echo "Provider s_vmem       : ${PROVIDER_S_VMEM}"
echo "Method s_vmem         : ${METHOD_S_VMEM}"
echo "Summary s_vmem        : ${SUMMARY_S_VMEM}"
echo "============================================================"

HOLD_JID=""
if [ "${RUN_PROVIDER}" = "1" ]; then
  PROVIDER_JOB="$(
    submit_and_get_job_id \
      qsub -N "mf_coarse_provider_${RUN_TAG}" \
        -l "s_vmem=${PROVIDER_S_VMEM}" \
        -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}" \
        -v "CONDA_SH=${CONDA_SH}" \
        -v "CONDA_ENV=${CONDA_ENV}" \
        -v "RUN_TAG=${RUN_TAG}" \
        -v "BACKUP_PROVIDER_OLD=${BACKUP_PROVIDER_OLD}" \
        "${JOB_SCRIPT_DIR}/run_marker_fm_coarse_provider_refresh.sh"
  )"
  if [ -z "${PROVIDER_JOB}" ]; then
    echo "ERROR: failed to parse provider refresh job id" >&2
    exit 2
  fi
  echo "Provider refresh job id: ${PROVIDER_JOB}"
  HOLD_JID="${PROVIDER_JOB}"
fi

METHOD_QSUB=(
  qsub
  -N "mf_coarse_methods_${RUN_TAG}"
  -t "${TASK_RANGE}"
  -l "s_vmem=${METHOD_S_VMEM}"
  -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}"
  -v "CONDA_SH=${CONDA_SH}"
  -v "CONDA_ENV=${CONDA_ENV}"
  -v "RUN_TAG=${RUN_TAG}"
  -v "BACKUP_OLD_LINEAGE=${BACKUP_OLD_LINEAGE}"
  -v "PURGE_METHOD_CACHE=${PURGE_METHOD_CACHE}"
)
if [ -n "${HOLD_JID}" ]; then
  METHOD_QSUB+=(-hold_jid "${HOLD_JID}")
fi
METHOD_QSUB+=("${JOB_SCRIPT_DIR}/run_marker_fm_coarse_method_rerun_array.sh")

METHOD_JOB="$(submit_and_get_job_id "${METHOD_QSUB[@]}")"
if [ -z "${METHOD_JOB}" ]; then
  echo "ERROR: failed to parse method array job id" >&2
  exit 2
fi
echo "Method rerun array job id: ${METHOD_JOB}"

if [ "${RUN_SUMMARY}" = "1" ]; then
  SUMMARY_JOB="$(
    submit_and_get_job_id \
      qsub -N "mf_coarse_summary_${RUN_TAG}" \
        -hold_jid "${METHOD_JOB}" \
        -l "s_vmem=${SUMMARY_S_VMEM}" \
        -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}" \
        -v "CONDA_SH=${CONDA_SH}" \
        -v "CONDA_ENV=${CONDA_ENV}" \
        -v "RUN_OFFICIAL_SUMMARY=${RUN_OFFICIAL_SUMMARY}" \
        -v "LEIDEN_N_NEIGHBORS=${LEIDEN_N_NEIGHBORS}" \
        -v "LEIDEN_RESOLUTION=${LEIDEN_RESOLUTION}" \
        "${JOB_SCRIPT_DIR}/run_marker_fm_coarse_rerun_summary.sh"
  )"
  if [ -z "${SUMMARY_JOB}" ]; then
    echo "ERROR: failed to parse summary job id" >&2
    exit 2
  fi
  echo "Summary job id: ${SUMMARY_JOB}"
fi

echo "============================================================"
echo "Submitted marker FM coarse-label rerun."
if [ "${RUN_PROVIDER}" = "1" ]; then
  echo "  provider: ${PROVIDER_JOB}"
fi
echo "  methods : ${METHOD_JOB} tasks ${TASK_RANGE}"
if [ "${RUN_SUMMARY}" = "1" ]; then
  echo "  summary : ${SUMMARY_JOB}"
fi
echo
echo "Useful commands:"
echo "  qstat -u ${USER:-unknown}"
echo "  ls logs/run_marker_fm_coarse_method_rerun.${METHOD_JOB}.*.log"
if [ "${RUN_SUMMARY}" = "1" ]; then
  echo "  ls logs/run_marker_fm_coarse_rerun_summary.${SUMMARY_JOB}.log"
fi
echo
echo "Override examples:"
echo "  TASK_RANGE=1-9 bash jobs/shirokane/submit_marker_fm_coarse_rerun.sh"
echo "  RUN_PROVIDER=0 RUN_OFFICIAL_SUMMARY=0 bash jobs/shirokane/submit_marker_fm_coarse_rerun.sh"
echo "  METHOD_S_VMEM=192G PURGE_METHOD_CACHE=1 bash jobs/shirokane/submit_marker_fm_coarse_rerun.sh"
echo "============================================================"
