#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
RUN_TAG="${RUN_TAG:-reference_ari_$(date +%Y%m%d_%H%M%S)}"
RUN_BACKFILL="${RUN_BACKFILL:-1}"

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
echo "Submitting marker-FM reference-ARI embedding rerun"
echo "Project root: ${PROJECT_ROOT}"
echo "Run tag     : ${RUN_TAG}"
echo "Resolution  : ${LEIDEN_RESOLUTION:-0.5}"
echo "Steps       : rerun 18 projection-capable embedding post-metrics only"
echo "Backfill    : ${RUN_BACKFILL}"
echo "============================================================"

HOLD_ARGS=()
if [ "${RUN_BACKFILL}" = "1" ]; then
  BACKFILL_JOB="$(
    submit_and_get_job_id \
      qsub -N mf_refari_bfill_${RUN_TAG} \
      run_marker_fm_transition_silver_embedding_backfill_array.sh
  )"
  if [ -z "${BACKFILL_JOB}" ]; then
    echo "ERROR: failed to parse embedding backfill array job id" >&2
    exit 2
  fi
  echo "Embedding backfill array job id: ${BACKFILL_JOB}"
  HOLD_ARGS=(-hold_jid "${BACKFILL_JOB}")
fi

EMBED_JOB="$(
  submit_and_get_job_id \
    qsub -N mf_refari_embed_${RUN_TAG} \
    "${HOLD_ARGS[@]}" \
    run_marker_fm_reference_ari_embedding_array.sh
)"
if [ -z "${EMBED_JOB}" ]; then
  echo "ERROR: failed to parse embedding array job id" >&2
  exit 2
fi
echo "Embedding array job id: ${EMBED_JOB}"

echo "============================================================"
echo "Submitted reference-ARI embedding rerun:"
if [ "${RUN_BACKFILL}" = "1" ]; then
  echo "  ${BACKFILL_JOB} backfill embedding.npy / next_timepoint_embedding.npy"
fi
echo "  ${EMBED_JOB} run tasks 1-18"
echo
echo "Useful commands:"
echo "  qstat -u ${USER}"
if [ "${RUN_BACKFILL}" = "1" ]; then
  echo "  qstat -j ${BACKFILL_JOB}"
  echo "  ls logs/run_marker_fm_transition_silver_embedding_backfill.${BACKFILL_JOB}.*.log"
fi
echo "  qstat -j ${EMBED_JOB}"
echo "  ls logs/run_marker_fm_reference_ari_embedding.${EMBED_JOB}.*.log"
echo "============================================================"
