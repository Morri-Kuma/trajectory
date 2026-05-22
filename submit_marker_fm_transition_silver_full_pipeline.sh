#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"

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
echo "Submitting marker-FM transition silver full rerun"
echo "Project root: ${PROJECT_ROOT}"
echo "Run tag     : ${RUN_TAG}"
echo "============================================================"

LABEL_JOB="$(submit_and_get_job_id qsub -N mf_rebuild_${RUN_TAG} -v RUN_TAG="${RUN_TAG}" run_marker_fm_transition_silver_rebuild_labels_and_providers.sh)"
if [ -z "${LABEL_JOB}" ]; then
  echo "ERROR: failed to parse label/provider job id" >&2
  exit 2
fi
echo "Label/provider job id: ${LABEL_JOB}"

METHOD_JOB="$(submit_and_get_job_id qsub -N mf_methods_${RUN_TAG} -hold_jid "${LABEL_JOB}" -v RUN_TAG="${RUN_TAG}" run_marker_fm_transition_silver_methods_array.sh)"
if [ -z "${METHOD_JOB}" ]; then
  echo "ERROR: failed to parse method array job id" >&2
  exit 2
fi
echo "Method array job id: ${METHOD_JOB}"

POST_JOB="$(submit_and_get_job_id qsub -N mf_embed_${RUN_TAG} -hold_jid "${METHOD_JOB}" run_marker_fm_transition_silver_post_embedding_array.sh)"
if [ -z "${POST_JOB}" ]; then
  echo "ERROR: failed to parse post-embedding array job id" >&2
  exit 2
fi
echo "Post-embedding array job id: ${POST_JOB}"

echo "============================================================"
echo "Submitted full pipeline:"
echo "  1. ${LABEL_JOB} rebuild labels/providers"
echo "  2. ${METHOD_JOB} run method array, held on ${LABEL_JOB}"
echo "  3. ${POST_JOB} rebuild projected labels + embedding metrics, held on ${METHOD_JOB}"
echo
echo "Useful commands:"
echo "  qstat -j ${LABEL_JOB}"
echo "  qstat -j ${METHOD_JOB}"
echo "  qstat -j ${POST_JOB}"
echo "  ls logs/*${RUN_TAG}*"
echo "============================================================"
