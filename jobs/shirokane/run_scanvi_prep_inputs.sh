#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scanvi_prep_inputs.$JOB_ID.log
#$ -l s_vmem=96G
# ── scANVI annotation branch — STEP 1/3: build shared-gene-space inputs ───────
#
# Joins FULL-GENE count matrices with the per-cell labels (which live only in the
# HVG2000 inputs) by cell barcode, via src/data/build_scanvi_inputs.py. This is
# required because the per-dataset HVG2000 inputs are gene-disjoint
# (ref∩query ≈ 181–238 genes); the full-gene matrices share ~17,829 genes incl.
# all canonical markers. See results/test_outputs/SERVER_DRYRUN_NOTES.md.
#
# REQUIREMENTS (verify before running):
#   * each --full-gene matrix must contain RAW integer counts (scVI input);
#   * the --full-gene and --labels files must share obs_names (cell barcodes).
# Override any path/env below to match your Shirokane checkout.
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${TRAJ_CONDA_ENV:-traj_env}"

OUT_DIR="${OUT_DIR:-results/annotation_branch/inputs}"
N_HVG_REF="${N_HVG_REF:-2000}"

# Reference (GSE242424, author_cluster_label) — full-gene counts + HVG2000 labels
REF_FULLGENE="${REF_FULLGENE:-data/processed/gse242424_human/GSE242424_raw_full_gene_benchmark_input.h5ad}"
REF_LABELS="${REF_LABELS:-benchmark/inputs/gse242424_author_cluster_matched/GSE242424_author_cluster_matched_HVG2000_benchmark_input.h5ad}"
REF_LABEL_KEY="${REF_LABEL_KEY:-author_cluster_label}"

# Queries (final_milestone_label_coarse) — full-gene counts + HVG2000 silver labels
Q178_FULLGENE="${Q178_FULLGENE:-data/processed/gse178325_human/GSE178325_0618_raw_full_gene_benchmark_input.h5ad}"
Q178_LABELS="${Q178_LABELS:-benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad}"
Q230_FULLGENE="${Q230_FULLGENE:-data/processed/gse230659_human/GSE230659_0618_raw_full_gene_benchmark_input.h5ad}"
Q230_LABELS="${Q230_LABELS:-benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad}"
Q_LABEL_KEY="${Q_LABEL_KEY:-final_milestone_label_coarse}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}" MKL_NUM_THREADS="${NSLOTS:-1}" OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "scANVI prep inputs | started $(date) | host $(hostname)"
echo "Project root : ${PROJECT_ROOT} | env ${CONDA_ENV} | out ${OUT_DIR}"
echo "============================================================"
[ -f "${CONDA_SH}" ] || { echo "ERROR: conda init not found: ${CONDA_SH}" >&2; exit 2; }
# shellcheck disable=SC1090
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"

mkdir -p "${OUT_DIR}"
python -m src.data.build_scanvi_inputs --full-gene "${REF_FULLGENE}"  --labels "${REF_LABELS}" \
    --label-key "${REF_LABEL_KEY}" --n-hvg "${N_HVG_REF}" \
    --out "${OUT_DIR}/gse242424_ref_fullgene_labelled.h5ad"
python -m src.data.build_scanvi_inputs --full-gene "${Q178_FULLGENE}" --labels "${Q178_LABELS}" \
    --label-key "${Q_LABEL_KEY}" --n-hvg 0 \
    --out "${OUT_DIR}/gse178325_query_fullgene_labelled.h5ad"
python -m src.data.build_scanvi_inputs --full-gene "${Q230_FULLGENE}" --labels "${Q230_LABELS}" \
    --label-key "${Q_LABEL_KEY}" --n-hvg 0 \
    --out "${OUT_DIR}/gse230659_query_fullgene_labelled.h5ad"

echo "scANVI prep inputs DONE $(date)"
