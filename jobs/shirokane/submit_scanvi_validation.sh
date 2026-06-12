#!/bin/bash
# Submit the validation experiments that complete the manuscript.
# A: geometric OOD on the completed run | B: positive control (GSE175634)
# C/D: same-protocol transfers (train->map each). No -hold_jid between A/B/C/D.
set -euo pipefail
cd "${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"; mkdir -p logs
JD="jobs/shirokane"
gid(){ awk '/Your job/{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+$/){print $i;exit}}'; }

echo "== A: geometric OOD on the completed scANVI run =="
A=$(qsub "$JD/run_scanvi_score_ood.sh" | tee /dev/stderr | gid)
echo "== B: positive control (GSE175634, decisive) =="
B=$(qsub "$JD/run_scanvi_positive_control.sh" | tee /dev/stderr | gid)
echo "== C: same-protocol GSE230659 -> GSE178325 =="
T1=$(qsub -v CONFIG=config_same_protocol_230to178.yaml,PROFILE=server "$JD/run_scanvi_train_reference.sh" | tee /dev/stderr | gid)
M1=$(qsub -hold_jid "$T1" -v CONFIG=config_same_protocol_230to178.yaml,PROFILE=server "$JD/run_scanvi_map.sh" | tee /dev/stderr | gid)
echo "== D: same-protocol GSE178325 -> GSE230659 =="
T2=$(qsub -v CONFIG=config_same_protocol_178to230.yaml,PROFILE=server "$JD/run_scanvi_train_reference.sh" | tee /dev/stderr | gid)
M2=$(qsub -hold_jid "$T2" -v CONFIG=config_same_protocol_178to230.yaml,PROFILE=server "$JD/run_scanvi_map.sh" | tee /dev/stderr | gid)
echo ""
echo "Submitted: OOD=$A  positive_control=$B  same230to178=$T1->$M1  same178to230=$T2->$M2"
echo "Outputs: results/annotation_branch/*/ood_scanvi_latent.csv | results/positive_control_gse175634/ | results/same_protocol_*/"
