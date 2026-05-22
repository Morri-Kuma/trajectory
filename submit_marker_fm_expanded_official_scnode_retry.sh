#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
RUN_TAG="${RUN_TAG:-expanded_official_scnode_retry_$(date +%Y%m%d_%H%M%S)}"

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
echo "Submitting marker-FM expanded-official scNODE retry"
echo "Project root: ${PROJECT_ROOT}"
echo "Run tag     : ${RUN_TAG}"
echo "Steps       : provider sanity refresh -> scNODE methods -> scNODE embedding metrics"
echo "============================================================"

PROVIDER_JOB="$(
  submit_and_get_job_id \
    qsub -N mf_exp_provider_retry_${RUN_TAG} \
    -v RUN_TAG="${RUN_TAG}" \
    run_marker_fm_expanded_official_refresh_providers.sh
)"
if [ -z "${PROVIDER_JOB}" ]; then
  echo "ERROR: failed to parse provider-refresh job id" >&2
  exit 2
fi
echo "Provider-refresh job id: ${PROVIDER_JOB}"

SCNODE_JOB="$(
  submit_and_get_job_id \
    qsub -N mf_exp_scnode_${RUN_TAG} \
    -hold_jid "${PROVIDER_JOB}" \
    -t 14-19 \
    -v RUN_TAG="${RUN_TAG}" \
    run_marker_fm_transition_silver_methods_array.sh
)"
if [ -z "${SCNODE_JOB}" ]; then
  echo "ERROR: failed to parse scNODE method-array job id" >&2
  exit 2
fi
echo "scNODE method-array job id: ${SCNODE_JOB}"

POST_JOB="$(
  submit_and_get_job_id \
    qsub -N mf_exp_scnode_embed_${RUN_TAG} \
    -hold_jid "${SCNODE_JOB}" \
    -t 13-18 \
    -v RUN_TAG="${RUN_TAG}" \
    run_marker_fm_transition_silver_post_embedding_array.sh
)"
if [ -z "${POST_JOB}" ]; then
  echo "ERROR: failed to parse scNODE post-embedding job id" >&2
  exit 2
fi
echo "scNODE post-embedding job id: ${POST_JOB}"

echo "============================================================"
echo "Submitted scNODE retry:"
echo "  1. ${PROVIDER_JOB} refresh/check expanded providers"
echo "  2. ${SCNODE_JOB} rerun methods array tasks 14-19, held on ${PROVIDER_JOB}"
echo "  3. ${POST_JOB} rerun post-embedding tasks 13-18, held on ${SCNODE_JOB}"
echo
echo "Useful commands:"
echo "  qstat -u ${USER}"
echo "  qstat -j ${PROVIDER_JOB} | grep -i hold"
echo "  qstat -j ${SCNODE_JOB} | grep -i hold"
echo "  qstat -j ${POST_JOB} | grep -i hold"
echo "  ls logs/*${RUN_TAG}*"
echo "============================================================"
