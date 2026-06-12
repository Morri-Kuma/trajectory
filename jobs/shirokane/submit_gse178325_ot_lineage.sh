#!/bin/bash
# Submit GSE178325 WOT/CellRank2 scenario-A official-silver lineage jobs.
#
# Typical use on Shirokane:
#   cd /home/xzy0723/projects/trajectory
#   bash jobs/shirokane/submit_gse178325_ot_lineage.sh

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
JOB_SCRIPT_DIR="${PROJECT_ROOT}/jobs/shirokane"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
RUN_TAG="${RUN_TAG:-gse178325_ot_lineage_$(date +%Y%m%d_%H%M%S)}"
TASK_RANGE="${TASK_RANGE:-1-2}"
S_VMEM="${S_VMEM:-192G}"
SUMMARY_S_VMEM="${SUMMARY_S_VMEM:-32G}"
RUN_SUMMARY="${RUN_SUMMARY:-1}"
RUN_OFFICIAL_SUMMARY="${RUN_OFFICIAL_SUMMARY:-1}"
BACKUP_OLD_LINEAGE="${BACKUP_OLD_LINEAGE:-1}"
PURGE_METHOD_CACHE="${PURGE_METHOD_CACHE:-0}"
APPLY_MANIFEST="${APPLY_MANIFEST:-1}"
LEIDEN_N_NEIGHBORS="${LEIDEN_N_NEIGHBORS:-15}"
LEIDEN_RESOLUTION="${LEIDEN_RESOLUTION:-0.5}"

cd "${PROJECT_ROOT}"
mkdir -p logs

REQUIRED_FILES=(
  "benchmark/configs/runtime/wot_gse178325_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/cellrank2_gse178325_marker_fm_silver_A_hvg2000_formal.yaml"
  "jobs/shirokane/run_gse178325_ot_lineage_array.sh"
  "jobs/shirokane/run_gse178325_ot_lineage_summary.sh"
  "scripts/update_gse178325_ot_manifest.py"
  "benchmark/evaluation/eval_dispatch.py"
  "benchmark/evaluation/eval_lineage.py"
  "benchmark/evaluation/lineage_graphsim_sctimebench.py"
  "benchmark/evaluation/summarize_official_silver.py"
  "benchmark/adapters/wot_adapter.py"
  "benchmark/adapters/cellrank2_adapter.py"
)

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
echo "Submitting GSE178325 WOT/CellRank2 official-silver lineage"
echo "Project root      : ${PROJECT_ROOT}"
echo "Run tag           : ${RUN_TAG}"
echo "Task range        : ${TASK_RANGE}"
echo "Method s_vmem     : ${S_VMEM}"
echo "Summary s_vmem    : ${SUMMARY_S_VMEM}"
echo "Run summary       : ${RUN_SUMMARY}"
echo "Official summary  : ${RUN_OFFICIAL_SUMMARY}"
echo "Backup old lineage: ${BACKUP_OLD_LINEAGE}"
echo "Purge method cache: ${PURGE_METHOD_CACHE}"
echo "Apply manifest    : ${APPLY_MANIFEST}"
echo "============================================================"

for required in "${REQUIRED_FILES[@]}"; do
  if [ ! -f "${required}" ]; then
    echo "ERROR: required project file is missing: ${required}" >&2
    echo "Sync the current trajectory project code to Shirokane before submitting." >&2
    exit 3
  fi
done

if [ "${APPLY_MANIFEST}" = "1" ]; then
  if [ ! -f "${CONDA_SH}" ]; then
    echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
    exit 2
  fi
  source "${CONDA_SH}"
  conda activate "${CONDA_ENV}"
  python scripts/update_gse178325_ot_manifest.py
fi

METHOD_JOB="$(
  submit_and_get_job_id \
    qsub -N "gse178325_ot_${RUN_TAG}" \
      -t "${TASK_RANGE}" \
      -l "s_vmem=${S_VMEM}" \
      -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}" \
      -v "CONDA_SH=${CONDA_SH}" \
      -v "CONDA_ENV=${CONDA_ENV}" \
      -v "RUN_TAG=${RUN_TAG}" \
      -v "BACKUP_OLD_LINEAGE=${BACKUP_OLD_LINEAGE}" \
      -v "PURGE_METHOD_CACHE=${PURGE_METHOD_CACHE}" \
      "${JOB_SCRIPT_DIR}/run_gse178325_ot_lineage_array.sh"
)"
if [ -z "${METHOD_JOB}" ]; then
  echo "ERROR: failed to parse method array job id" >&2
  exit 2
fi
echo "Method array job id: ${METHOD_JOB}"

if [ "${RUN_SUMMARY}" = "1" ]; then
  SUMMARY_JOB="$(
    submit_and_get_job_id \
      qsub -N "gse178325_ot_sum_${RUN_TAG}" \
        -hold_jid "${METHOD_JOB}" \
        -l "s_vmem=${SUMMARY_S_VMEM}" \
        -v "TRAJ_PROJECT_ROOT=${PROJECT_ROOT}" \
        -v "CONDA_SH=${CONDA_SH}" \
        -v "CONDA_ENV=${CONDA_ENV}" \
        -v "RUN_OFFICIAL_SUMMARY=${RUN_OFFICIAL_SUMMARY}" \
        -v "LEIDEN_N_NEIGHBORS=${LEIDEN_N_NEIGHBORS}" \
        -v "LEIDEN_RESOLUTION=${LEIDEN_RESOLUTION}" \
        "${JOB_SCRIPT_DIR}/run_gse178325_ot_lineage_summary.sh"
  )"
  if [ -z "${SUMMARY_JOB}" ]; then
    echo "ERROR: failed to parse summary job id" >&2
    exit 2
  fi
  echo "Summary job id: ${SUMMARY_JOB}"
fi

echo "============================================================"
echo "Submitted GSE178325 OT lineage jobs."
echo "  methods: ${METHOD_JOB} tasks ${TASK_RANGE}"
if [ "${RUN_SUMMARY}" = "1" ]; then
  echo "  summary: ${SUMMARY_JOB}"
fi
echo
echo "Useful commands:"
echo "  qstat -u ${USER:-unknown}"
echo "  qstat -j ${METHOD_JOB}"
echo "  ls logs/run_gse178325_ot_lineage.${METHOD_JOB}.*.log"
if [ "${RUN_SUMMARY}" = "1" ]; then
  echo "  qstat -j ${SUMMARY_JOB}"
  echo "  ls logs/run_gse178325_ot_lineage_summary.${SUMMARY_JOB}.log"
fi
echo
echo "Task mapping:"
echo "  1 = WOT GSE178325-A"
echo "  2 = CellRank2 GSE178325-A"
echo
echo "Override examples:"
echo "  S_VMEM=256G bash jobs/shirokane/submit_gse178325_ot_lineage.sh"
echo "  TASK_RANGE=1 RUN_SUMMARY=0 bash jobs/shirokane/submit_gse178325_ot_lineage.sh"
echo "  PURGE_METHOD_CACHE=1 bash jobs/shirokane/submit_gse178325_ot_lineage.sh"
echo "============================================================"
