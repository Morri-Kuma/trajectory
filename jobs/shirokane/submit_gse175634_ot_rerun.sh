#!/bin/bash
# Rerun only the GSE175634 OT-based tasks that need per-timepoint subsampling:
#   task 10: WOT
#   task 11: CellRank2
#
# Upload this script together with the updated adapter/config files, then run:
#   cd /home/xzy0723/projects/trajectory
#   bash jobs/shirokane/submit_gse175634_ot_rerun.sh

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
JOB_SCRIPT_DIR="${PROJECT_ROOT}/jobs/shirokane"
RUN_TAG="${RUN_TAG:-gse175634_ot_rerun_$(date +%Y%m%d_%H%M%S)}"
S_VMEM="${S_VMEM:-96G}"

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
echo "Submitting GSE175634 WOT/CellRank2 rerun"
echo "Project root : ${PROJECT_ROOT}"
echo "Run tag      : ${RUN_TAG}"
echo "s_vmem       : ${S_VMEM}"
echo "Tasks        : 10-11"
echo "============================================================"

BACKUP_DIR="benchmark/results/_rerun_backup/${RUN_TAG}"
mkdir -p "${BACKUP_DIR}"

for result_dir in \
  "benchmark/results/wot/gse175634_cardiac_silver_A_hvg2000_formal" \
  "benchmark/results/cellrank2/gse175634_cardiac_silver_A_hvg2000_formal"
do
  if [ -e "${result_dir}" ]; then
    mkdir -p "${BACKUP_DIR}/$(dirname "${result_dir}")"
    echo "Moving existing ${result_dir} to ${BACKUP_DIR}/${result_dir}"
    mv "${result_dir}" "${BACKUP_DIR}/${result_dir}"
  fi
done

OT_JOB="$(
  submit_and_get_job_id \
    qsub -N "gse175634_ot_${RUN_TAG}" \
         -t 10-11 \
         -l "s_vmem=${S_VMEM}" \
         "${JOB_SCRIPT_DIR}/run_gse175634_benchmark_array.sh"
)"

if [ -z "${OT_JOB}" ]; then
  echo "ERROR: failed to parse OT rerun job id" >&2
  exit 2
fi

echo ""
echo "Submitted OT rerun job: ${OT_JOB}"
echo "Useful commands:"
echo "  qstat -u ${USER}"
echo "  qstat -j ${OT_JOB}"
echo "  ls logs/run_gse175634_benchmark.${OT_JOB}.*.log"
echo "============================================================"
