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
echo "Submitting marker-FM expanded-official rerun"
echo "Project root: ${PROJECT_ROOT}"
echo "Run tag     : ${RUN_TAG}"
echo "Steps       : provider refresh -> methods -> embedding metrics"
echo "============================================================"

PROVIDER_JOB="$(
  submit_and_get_job_id \
    qsub -N mf_exp_provider_${RUN_TAG} \
    -v RUN_TAG="${RUN_TAG}" \
    run_marker_fm_expanded_official_refresh_providers.sh
)"
if [ -z "${PROVIDER_JOB}" ]; then
  echo "ERROR: failed to parse provider-refresh job id" >&2
  exit 2
fi
echo "Provider-refresh job id: ${PROVIDER_JOB}"

METHOD_JOB="$(
  submit_and_get_job_id \
    qsub -N mf_exp_methods_${RUN_TAG} \
    -hold_jid "${PROVIDER_JOB}" \
    -v RUN_TAG="${RUN_TAG}" \
    run_marker_fm_transition_silver_methods_array.sh
)"
if [ -z "${METHOD_JOB}" ]; then
  echo "ERROR: failed to parse method-array job id" >&2
  exit 2
fi
echo "Method-array job id: ${METHOD_JOB}"

POST_JOB="$(
  submit_and_get_job_id \
    qsub -N mf_exp_embed_${RUN_TAG} \
    -hold_jid "${METHOD_JOB}" \
    -v RUN_TAG="${RUN_TAG}" \
    run_marker_fm_transition_silver_post_embedding_array.sh
)"
if [ -z "${POST_JOB}" ]; then
  echo "ERROR: failed to parse post-embedding job id" >&2
  exit 2
fi
echo "Post-embedding job id: ${POST_JOB}"

echo "============================================================"
echo "Submitted expanded-official rerun:"
echo "  1. ${PROVIDER_JOB} refresh GSE178325/GSE230659 expanded providers"
echo "  2. ${METHOD_JOB} run all marker-FM method configs, held on ${PROVIDER_JOB}"
echo "  3. ${POST_JOB} rebuild projected labels + embedding metrics, held on ${METHOD_JOB}"
echo
echo "Useful commands:"
echo "  qstat -u ${USER}"
echo "  qstat -j ${PROVIDER_JOB} | grep -i hold"
echo "  qstat -j ${METHOD_JOB} | grep -i hold"
echo "  qstat -j ${POST_JOB} | grep -i hold"
echo "  ls logs/*${RUN_TAG}*"
echo "============================================================"
