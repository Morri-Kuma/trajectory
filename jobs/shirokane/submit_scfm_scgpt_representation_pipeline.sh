#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
JOB_SCRIPT_DIR="${PROJECT_ROOT}/jobs/shirokane"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
TRAJ_CONDA_ENV="${TRAJ_CONDA_ENV:-traj_env}"
SCGPT_CONDA_ENV="${SCGPT_CONDA_ENV:-scgpt}"
SCGPT_MODEL_DIR="${SCGPT_MODEL_DIR:-/home/xzy0723/projects/trajectory/models/scgpt_whole_human}"
RUN_TAG="${RUN_TAG:-scgpt_rep_$(date +%Y%m%d_%H%M%S)}"

# Default to smoke: MIOFlow + Scenario B. Use TRAJ_TASK_RANGE=1-9 for all
# MIOFlow/scNODE/PRESCIENT x A/B/C tasks after smoke succeeds.
RUN_TRAJECTORY="${RUN_TRAJECTORY:-1}"
TRAJ_TASK_RANGE="${TRAJ_TASK_RANGE:-2}"

# Shirokane resource flags. The Shirokane user guide selects GPU queues by GPU
# type, e.g. -l h100=1 for gpuh.q, -l a100=1 for gpua.q, or -l v100=1 for
# gpuv.q. Keep this configurable because queue policies can differ by course.
PREP_RESOURCES="${PREP_RESOURCES:--l s_vmem=96G}"
EXTRACT_RESOURCES="${EXTRACT_RESOURCES:--l s_vmem=128G -l h100=1}"
XREP_RESOURCES="${XREP_RESOURCES:--l s_vmem=96G}"
TRAJ_RESOURCES="${TRAJ_RESOURCES:--l s_vmem=96G}"

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
echo "Submitting scGPT representation pipeline"
echo "Project root       : ${PROJECT_ROOT}"
echo "Run tag            : ${RUN_TAG}"
echo "Trajectory env     : ${TRAJ_CONDA_ENV}"
echo "scGPT env          : ${SCGPT_CONDA_ENV}"
echo "scGPT model dir    : ${SCGPT_MODEL_DIR}"
echo "Run trajectory     : ${RUN_TRAJECTORY}"
echo "Trajectory tasks   : ${TRAJ_TASK_RANGE}"
echo "Prep resources     : ${PREP_RESOURCES}"
echo "Extract resources  : ${EXTRACT_RESOURCES}"
echo "X_rep resources    : ${XREP_RESOURCES}"
echo "Trajectory resources: ${TRAJ_RESOURCES}"
echo "============================================================"

if [ ! -f "${SCGPT_MODEL_DIR}/vocab.json" ] || \
   [ ! -f "${SCGPT_MODEL_DIR}/args.json" ] || \
   [ ! -f "${SCGPT_MODEL_DIR}/best_model.pt" ]; then
  echo "ERROR: scGPT checkpoint directory is missing vocab.json/args.json/best_model.pt: ${SCGPT_MODEL_DIR}" >&2
  exit 2
fi

# shellcheck disable=SC2206
PREP_RESOURCE_ARGS=(${PREP_RESOURCES})
# shellcheck disable=SC2206
EXTRACT_RESOURCE_ARGS=(${EXTRACT_RESOURCES})
# shellcheck disable=SC2206
XREP_RESOURCE_ARGS=(${XREP_RESOURCES})
# shellcheck disable=SC2206
TRAJ_RESOURCE_ARGS=(${TRAJ_RESOURCES})

PREP_JOB="$(
  submit_and_get_job_id \
    qsub -N "scgpt_prep_${RUN_TAG}" \
      "${PREP_RESOURCE_ARGS[@]}" \
      -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}" \
      -v "CONDA_SH=${CONDA_SH}" \
      -v "CONDA_ENV=${TRAJ_CONDA_ENV}" \
      "${JOB_SCRIPT_DIR}/run_scfm_scgpt_prepare_source.sh"
)"
if [ -z "${PREP_JOB}" ]; then
  echo "ERROR: failed to parse prepare-source job id" >&2
  exit 2
fi
echo "Prepare-source job id: ${PREP_JOB}"

EXTRACT_JOB="$(
  submit_and_get_job_id \
    qsub -N "scgpt_emb_${RUN_TAG}" \
      -hold_jid "${PREP_JOB}" \
      "${EXTRACT_RESOURCE_ARGS[@]}" \
      -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}" \
      -v "CONDA_SH=${CONDA_SH}" \
      -v "SCGPT_CONDA_ENV=${SCGPT_CONDA_ENV}" \
      -v "SCGPT_MODEL_DIR=${SCGPT_MODEL_DIR}" \
      "${JOB_SCRIPT_DIR}/run_scfm_scgpt_extract_embedding.sh"
)"
if [ -z "${EXTRACT_JOB}" ]; then
  echo "ERROR: failed to parse scGPT embedding job id" >&2
  exit 2
fi
echo "scGPT embedding job id: ${EXTRACT_JOB}"

XREP_JOB="$(
  submit_and_get_job_id \
    qsub -N "scgpt_xrep_${RUN_TAG}" \
      -hold_jid "${EXTRACT_JOB}" \
      -t 1-3 \
      "${XREP_RESOURCE_ARGS[@]}" \
      -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}" \
      -v "CONDA_SH=${CONDA_SH}" \
      -v "CONDA_ENV=${TRAJ_CONDA_ENV}" \
      "${JOB_SCRIPT_DIR}/run_scfm_scgpt_build_xrep_array.sh"
)"
if [ -z "${XREP_JOB}" ]; then
  echo "ERROR: failed to parse X_rep array job id" >&2
  exit 2
fi
echo "X_rep array job id: ${XREP_JOB}"

if [ "${RUN_TRAJECTORY}" = "1" ]; then
  TRAJ_JOB="$(
    submit_and_get_job_id \
      qsub -N "scgpt_traj_${RUN_TAG}" \
        -hold_jid "${XREP_JOB}" \
        -t "${TRAJ_TASK_RANGE}" \
        "${TRAJ_RESOURCE_ARGS[@]}" \
        -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}" \
        -v "CONDA_SH=${CONDA_SH}" \
        -v "CONDA_ENV=${TRAJ_CONDA_ENV}" \
        "${JOB_SCRIPT_DIR}/run_scfm_scgpt_trajectory_array.sh"
  )"
  if [ -z "${TRAJ_JOB}" ]; then
    echo "ERROR: failed to parse trajectory array job id" >&2
    exit 2
  fi
  echo "Trajectory array job id: ${TRAJ_JOB}"
fi

echo "============================================================"
echo "Submitted scGPT representation pipeline"
echo "  ${PREP_JOB} prepare labeled full-gene h5ad"
echo "  ${EXTRACT_JOB} extract scGPT embedding after ${PREP_JOB}"
echo "  ${XREP_JOB} build A/B/C X_rep after ${EXTRACT_JOB}"
if [ "${RUN_TRAJECTORY}" = "1" ]; then
  echo "  ${TRAJ_JOB} run trajectory tasks ${TRAJ_TASK_RANGE} after ${XREP_JOB}"
fi
echo
echo "Useful commands:"
echo "  qstat -u ${USER:-unknown}"
echo "  ls logs/run_scfm_scgpt_prepare_source.*.log"
echo "  ls logs/run_scfm_scgpt_extract_embedding.*.log"
echo "  ls logs/run_scfm_scgpt_build_xrep.*.*.log"
echo "  ls logs/run_scfm_scgpt_trajectory.*.*.log"
echo
echo "After smoke succeeds, run all trajectory tasks with:"
echo "  TRAJ_TASK_RANGE=1-9 bash jobs/shirokane/submit_scfm_scgpt_representation_pipeline.sh"
echo
echo "If H100 is unavailable, retry extraction with an available GPU resource, e.g.:"
echo "  EXTRACT_RESOURCES='-l s_vmem=128G -l a100=1' bash jobs/shirokane/submit_scfm_scgpt_representation_pipeline.sh"
echo "  EXTRACT_RESOURCES='-l s_vmem=128G -l v100=1' bash jobs/shirokane/submit_scfm_scgpt_representation_pipeline.sh"
echo "============================================================"
