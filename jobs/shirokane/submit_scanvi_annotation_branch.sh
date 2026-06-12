#!/bin/bash
# Submit the scANVI annotation branch as a 3-stage dependent chain on Shirokane:
#   1) prep inputs (CPU)  ->  2) train reference (GPU)  ->  3) map + compare (GPU)
# Each stage waits for the previous via -hold_jid. Override resources/paths via env.
set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
JOB_DIR="${PROJECT_ROOT}/jobs/shirokane"
PREP_RESOURCES="${PREP_RESOURCES:--l s_vmem=96G}"
TRAIN_RESOURCES="${TRAIN_RESOURCES:--l s_vmem=128G -l h100=1}"
MAP_RESOURCES="${MAP_RESOURCES:--l s_vmem=128G -l h100=1}"

cd "${PROJECT_ROOT}"; mkdir -p logs

get_jid() {  # parse "Your job NNN ..." from qsub output
  awk '/Your job/ {for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/){print $i; exit}}'
}

echo "Submitting scANVI annotation branch from ${PROJECT_ROOT}"

PREP_JID=$(qsub ${PREP_RESOURCES} "${JOB_DIR}/run_scanvi_prep_inputs.sh" | tee /dev/stderr | get_jid)
echo "  prep   job id: ${PREP_JID}"

TRAIN_JID=$(qsub ${TRAIN_RESOURCES} -hold_jid "${PREP_JID}" "${JOB_DIR}/run_scanvi_train_reference.sh" | tee /dev/stderr | get_jid)
echo "  train  job id: ${TRAIN_JID} (holds on ${PREP_JID})"

MAP_JID=$(qsub ${MAP_RESOURCES} -hold_jid "${TRAIN_JID}" "${JOB_DIR}/run_scanvi_map_compare.sh" | tee /dev/stderr | get_jid)
echo "  map    job id: ${MAP_JID} (holds on ${TRAIN_JID})"

echo "Submitted chain: ${PREP_JID} -> ${TRAIN_JID} -> ${MAP_JID}"
echo "Outputs will appear under results/annotation_branch/ and models/scanvi_reference_gse242424/"
