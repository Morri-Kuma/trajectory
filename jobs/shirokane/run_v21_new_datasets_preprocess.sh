#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N v21_prep
#$ -o logs/run_v21_new_datasets_preprocess.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=128G
#$ -t 1-2
# ── v2.1 preprocessing: build HVG2000 inputs + official_silver providers ──────
#   for GSE298212 (human blood chemical reprogramming) and
#       GSE218855 (mouse MEF fast chemical reprogramming).
#
# Task 1 : GSE298212    Task 2 : GSE218855
#
# THREE sub-steps per dataset:
#   Step 1  Build the HVG2000 benchmark input h5ad from raw data.
#           >>> dataset-specific: edit BUILD_INPUT_CMD below to match the raw
#           layout on Shirokane (mirror scripts/build_gse175634_cardiac_author_input.py
#           or scripts/build_gse242424_*; both datasets ship 10x MTX). The output
#           path must equal H5AD below (the runtime configs already point there).
#   Step 2  Annotate official_silver milestones (final_milestone_label_coarse)
#           with the marker-FM workflow:
#             benchmark/annotation/build_marker_seed_labels.py      (Stage 1)
#             benchmark/annotation/build_trajectory_aware_labels.py (Stage 2)
#           GSE298212 reuses the human chemical milestone markers
#           (benchmark/annotation/milestone_markers.yaml; same hADSCs ->
#           epithelial_like -> intermediate_plastic -> hCiPS schema as
#           GSE178325/GSE230659). GSE218855 (mouse) needs the mouse marker set
#           benchmark/annotation/milestone_markers_mouse.yaml (MEF -> ... -> iPSC).
#   Step 3  Freeze the provider with build_milestone_providers.py ->
#             benchmark/ground_truth/providers/<provider>/reference_graph.json
#
# After this completes, submit run_v21_new_datasets_benchmark_array.sh.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"

case "${TASK_ID}" in
  1) DS="GSE298212"; SPECIES="human"
     RAW="data/gse298212"
     H5AD="benchmark/inputs/gse298212_marker_fm_transition_silver_hvg2000/GSE298212_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"
     PROV="gse298212_marker_fm_transition_silver_v1"
     MARKERS="benchmark/annotation/milestone_markers.yaml"
     BUILD_INPUT_CMD="python scripts/build_gse298212_chemical_input.py"
     START_LABEL="blood_start"; START_SAMPLES="pbmc_epc,PBMC_EPC,EPC"
     ;;
  2) DS="GSE218855"; SPECIES="mouse"
     RAW="data/gse218855"
     H5AD="benchmark/inputs/gse218855_marker_fm_transition_silver_hvg2000/GSE218855_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"
     PROV="gse218855_marker_fm_transition_silver_v1"
     MARKERS="benchmark/annotation/milestone_markers_mouse.yaml"
     BUILD_INPUT_CMD="python scripts/build_gse218855_fcr_input.py"
     START_LABEL="MEF"; START_SAMPLES="D0,d0"
     ;;
  *) echo "ERROR: task id ${TASK_ID} not in 1-2" >&2; exit 2 ;;
esac

mkdir -p "${PROJECT_ROOT}/logs"; cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"; export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"
echo "=== v2.1 preprocess ${DS} (${SPECIES}) at $(date) ==="

# Step 1 — build HVG2000 input
if [ -f "${H5AD}" ]; then
  echo "[skip] HVG2000 input already present: ${H5AD}"
else
  echo "[Step 1] building HVG2000 input -> ${H5AD}"
  echo "         (raw dir: ${RAW}; builder: ${BUILD_INPUT_CMD})"
  if [ ! -d "${RAW}" ]; then
    echo "ERROR: raw data dir not found: ${RAW}." >&2
    echo "  Acquire ${DS} (see datasets/dataset_inventory.md) and point BUILD_INPUT_CMD at it." >&2
    exit 4
  fi
  ${BUILD_INPUT_CMD}
fi
[ -f "${H5AD}" ] || { echo "ERROR: input not produced: ${H5AD}" >&2; exit 3; }

# Step 2 — marker-FM milestone annotation (official_silver)
echo "[Step 2] marker-FM annotation using ${MARKERS}"
if [ ! -f "${MARKERS}" ]; then
  echo "ERROR: marker file missing: ${MARKERS}" >&2
  if [ "${SPECIES}" = "mouse" ]; then
    echo "  Create benchmark/annotation/milestone_markers_mouse.yaml (mouse MEF->iPSC markers)." >&2
  fi
  exit 4
fi
ANN_DIR="benchmark/annotation_runs/${DS,,}_v21"
SEED_H5AD="${H5AD%.h5ad}_stage1_seed.h5ad"
mkdir -p "${ANN_DIR}"
# Stage 1: marker seed labels. var_names are gene symbols (no gene_symbol var col),
# so fall back to the var index. The start milestone is assigned by abs_day==0.
python benchmark/annotation/build_marker_seed_labels.py \
  --input-h5ad "${H5AD}" --output-h5ad "${SEED_H5AD}" \
  --dataset-id "${DS}" --markers-yaml "${MARKERS}" --use-var-index-if-needed \
  --sample-key sample_id --time-key abs_day \
  --hadsc-label "${START_LABEL}" --hadsc-time-values 0,0.0 --hadsc-sample-values "${START_SAMPLES}" \
  --output-dir "${ANN_DIR}"
# Stage 2: trajectory-aware final labels, written BACK into the benchmark input h5ad
# (the dispatcher reads that file and needs final_milestone_label_coarse in it).
python benchmark/annotation/build_trajectory_aware_labels.py \
  --input-h5ad "${SEED_H5AD}" --output-h5ad "${H5AD}" \
  --dataset-id "${DS}" --markers-yaml "${MARKERS}" --output-dir "${ANN_DIR}"

# Step 2.5 — inject a reference UMAP into .obsm so the official embedding-coherence
# evaluator (eval_embedding_milestone.py, detection priority [X_pca, X_umap]) has a
# reference embedding, matching the primary datasets which expose X_umap. Idempotent.
echo "[Step 2.5] injecting reference UMAP into ${H5AD}"
python scripts/inject_reference_umap.py --input-h5ad "${H5AD}"

# Step 3 — freeze official_silver provider from the labelled h5ad
echo "[Step 3] freezing provider ${PROV}"
python benchmark/annotation/build_milestone_providers.py \
  --dataset-id "${DS}" --label-mode official_silver --input-h5ad "${H5AD}" \
  --markers-yaml "${MARKERS}"

GRAPH="benchmark/ground_truth/providers/${PROV}/reference_graph.json"
[ -f "${GRAPH}" ] || { echo "ERROR: provider not produced: ${GRAPH}" >&2; exit 3; }
echo "=== ${DS} preprocess DONE at $(date); provider ${GRAPH} ==="
echo "Next: qsub jobs/shirokane/run_v21_new_datasets_benchmark_array.sh"
