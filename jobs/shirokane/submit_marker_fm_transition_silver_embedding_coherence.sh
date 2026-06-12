#!/bin/bash

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
JOB_SCRIPT_DIR="${PROJECT_ROOT}/jobs/shirokane"
RUN_TAG="${RUN_TAG:-sctimebench_embedding_$(date +%Y%m%d_%H%M%S)}"
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
echo "Submitting marker-FM scTimeBench Embedding Coherence"
echo "Project root : ${PROJECT_ROOT}"
echo "Run tag      : ${RUN_TAG}"
echo "Tasks        : 18 projection-capable official-silver runs"
echo "Protocol     : embedding.npy + next_timepoint_embedding.npy"
echo "Leiden       : n_neighbors=15, resolution=0.5"
echo "Entropy      : RandomForest classifier probability entropy"
echo "Backfill     : ${RUN_BACKFILL}"
echo "============================================================"

HOLD_ARGS=()
if [ "${RUN_BACKFILL}" = "1" ]; then
  BACKFILL_JOB="$(
    submit_and_get_job_id \
      qsub -N mf_embed_bfill_${RUN_TAG} \
      "${JOB_SCRIPT_DIR}/run_marker_fm_transition_silver_embedding_backfill_array.sh"
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
    qsub -N mf_embed_sctb_${RUN_TAG} \
    "${HOLD_ARGS[@]}" \
    "${JOB_SCRIPT_DIR}/run_marker_fm_transition_silver_embedding_coherence_array.sh"
)"
if [ -z "${EMBED_JOB}" ]; then
  echo "ERROR: failed to parse embedding array job id" >&2
  exit 2
fi
echo "Embedding array job id: ${EMBED_JOB}"

echo "============================================================"
echo "Submitted scTimeBench Embedding Coherence rerun:"
if [ "${RUN_BACKFILL}" = "1" ]; then
  echo "  ${BACKFILL_JOB} backfill embedding.npy / next_timepoint_embedding.npy"
fi
echo "  ${EMBED_JOB} run tasks 1-18"
echo
echo "Useful commands:"
echo "  qstat -u ${USER:-unknown}"
if [ "${RUN_BACKFILL}" = "1" ]; then
  echo "  qstat -j ${BACKFILL_JOB}"
  echo "  ls logs/run_marker_fm_transition_silver_embedding_backfill.${BACKFILL_JOB}.*.log"
fi
echo "  qstat -j ${EMBED_JOB}"
echo "  ls logs/run_marker_fm_transition_silver_embedding.${EMBED_JOB}.*.log"
echo
echo "After all tasks complete, regenerate official-silver tables with:"
echo "  python -m benchmark.evaluation.summarize_official_silver --leiden-resolution 0.5"
echo "============================================================"
