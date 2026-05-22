#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_marker_fm_transition_silver_rebuild_labels.$JOB_ID.log
#$ -l s_vmem=128G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
REBUILD_RAW_FULL_GENE="${REBUILD_RAW_FULL_GENE:-0}"
STAGE1_MARGIN="${STAGE1_MARGIN:-0.10}"
TRANSITION_MARGIN="${TRANSITION_MARGIN:-0.10}"
BACKUP_EXISTING="${BACKUP_EXISTING:-1}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"

BACKUP_ROOT="${PROJECT_ROOT}/benchmark/_rerun_backup/${RUN_TAG}"
TMP_ROOT="${PROJECT_ROOT}/benchmark/_rerun_tmp/${RUN_TAG}"
mkdir -p "${TMP_ROOT}"

echo "============================================================"
echo "Marker-FM label/provider rebuild started at: $(date)"
echo "Host                  : $(hostname)"
echo "Project root          : ${PROJECT_ROOT}"
echo "Run tag               : ${RUN_TAG}"
echo "Stage1 margin         : ${STAGE1_MARGIN}"
echo "Stage2 transition     : ${TRANSITION_MARGIN}"
echo "Rebuild raw full gene : ${REBUILD_RAW_FULL_GENE}"
echo "Backup existing       : ${BACKUP_EXISTING}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

backup_path() {
  local path="$1"
  if [ "${BACKUP_EXISTING}" != "1" ] || [ ! -e "${path}" ]; then
    return 0
  fi
  local backup="${BACKUP_ROOT}/${path}"
  mkdir -p "$(dirname "${backup}")"
  if [ -e "${backup}" ]; then
    backup="${backup}.$(date +%s)"
  fi
  echo "[backup] ${path} -> ${backup#${PROJECT_ROOT}/}"
  mv "${path}" "${backup}"
}

backup_copy() {
  local path="$1"
  if [ "${BACKUP_EXISTING}" != "1" ] || [ ! -f "${path}" ]; then
    return 0
  fi
  local backup="${BACKUP_ROOT}/${path}"
  mkdir -p "$(dirname "${backup}")"
  if [ -e "${backup}" ]; then
    backup="${backup}.$(date +%s)"
  fi
  echo "[backup-copy] ${path} -> ${backup#${PROJECT_ROOT}/}"
  cp "${path}" "${backup}"
}

promote_file() {
  local tmp="$1"
  local final="$2"
  if [ ! -f "${tmp}" ]; then
    echo "ERROR: temporary file missing: ${tmp}" >&2
    exit 3
  fi
  if [ "${BACKUP_EXISTING}" != "1" ] && [ -e "${final}" ]; then
    echo "ERROR: ${final} already exists and BACKUP_EXISTING is not 1." >&2
    exit 3
  fi
  backup_path "${final}"
  mkdir -p "$(dirname "${final}")"
  mv "${tmp}" "${final}"
  echo "[promote] ${final}"
}

promote_dir() {
  local tmp="$1"
  local final="$2"
  if [ ! -d "${tmp}" ]; then
    echo "ERROR: temporary dir missing: ${tmp}" >&2
    exit 3
  fi
  if [ "${BACKUP_EXISTING}" != "1" ] && [ -e "${final}" ]; then
    echo "ERROR: ${final} already exists and BACKUP_EXISTING is not 1." >&2
    exit 3
  fi
  backup_path "${final}"
  mkdir -p "$(dirname "${final}")"
  mv "${tmp}" "${final}"
  echo "[promote] ${final}"
}

ensure_raw_full_gene() {
  local dataset="$1"
  local raw_path="$2"
  local builder="$3"
  if [ "${REBUILD_RAW_FULL_GENE}" = "1" ] || [ ! -f "${raw_path}" ]; then
    echo "------------------------------------------------------------"
    echo "[${dataset}] Building raw/full-gene input"
    echo "------------------------------------------------------------"
    python "${builder}"
  else
    echo "[${dataset}] Raw/full-gene input exists; skipping rebuild: ${raw_path}"
  fi
  if [ ! -f "${raw_path}" ]; then
    echo "ERROR: raw/full-gene input still missing: ${raw_path}" >&2
    exit 4
  fi
}

transfer_stage1_to_hvg() {
  local dataset="$1"
  local stage1_full="$2"
  local hvg_base="$3"
  local stage1_hvg_tmp="$4"
  python - "${dataset}" "${stage1_full}" "${hvg_base}" "${stage1_hvg_tmp}" "${RUN_TAG}" <<'PY'
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad

dataset, stage1_full, hvg_base, stage1_hvg_tmp, run_tag = sys.argv[1:6]
stage1_full = Path(stage1_full)
hvg_base = Path(hvg_base)
stage1_hvg_tmp = Path(stage1_hvg_tmp)

print(f"[transfer:{dataset}] source stage1 full-gene: {stage1_full}")
print(f"[transfer:{dataset}] HVG base              : {hvg_base}")
src = ad.read_h5ad(stage1_full)
base = ad.read_h5ad(hvg_base)

if not base.obs_names.equals(src.obs_names):
    missing = base.obs_names.difference(src.obs_names)
    extra = src.obs_names.difference(base.obs_names)
    if len(missing) or len(extra):
        raise SystemExit(
            f"obs_names mismatch for {dataset}: missing_in_stage1={len(missing)}, "
            f"extra_in_stage1={len(extra)}"
        )
    src = src[base.obs_names].copy()

stage1_cols = [
    "milestone_marker_label",
    "milestone_marker_score",
    "milestone_marker_margin",
    "stage1_label",
    "stage1_confidence",
    "stage1_margin",
    "stage1_status",
    "stage1_source",
]
score_cols = [c for c in src.obs.columns if c.startswith("milestone_score_")]
cols = score_cols + [c for c in stage1_cols if c in src.obs.columns]
missing_required = [c for c in stage1_cols if c not in src.obs.columns]
if missing_required:
    raise SystemExit(f"missing required Stage 1 columns in {stage1_full}: {missing_required}")

for col in cols:
    base.obs[col] = src.obs[col].to_numpy()

base.uns["milestone_annotation_step3"] = src.uns.get("milestone_annotation_step3", {})
base.uns["stage1_annotation_transfer"] = {
    "dataset_id": dataset,
    "run_tag": run_tag,
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "source_stage1_full_gene_h5ad": str(stage1_full),
    "hvg_base_h5ad": str(hvg_base),
    "output_stage1_hvg_h5ad": str(stage1_hvg_tmp),
    "transferred_obs_columns": cols,
    "transfer_key": "obs_names",
    "note": "Stage 1 marker labels were computed on full-gene expression and transferred to the HVG2000 benchmark matrix.",
}

stage1_hvg_tmp.parent.mkdir(parents=True, exist_ok=True)
base.write_h5ad(stage1_hvg_tmp)

counts = base.obs["stage1_label"].astype(str).value_counts().to_dict()
meta_path = stage1_hvg_tmp.with_suffix(".transfer_metadata.json")
meta_path.write_text(json.dumps({
    "dataset_id": dataset,
    "n_cells": int(base.n_obs),
    "n_vars": int(base.n_vars),
    "stage1_label_counts": {str(k): int(v) for k, v in counts.items()},
    "transferred_obs_columns": cols,
}, indent=2), encoding="utf-8")
print(f"[transfer:{dataset}] wrote {stage1_hvg_tmp}")
print(f"[transfer:{dataset}] stage1_label_counts={counts}")
PY
}

rebuild_dataset() {
  local dataset="$1"
  local raw_path="$2"
  local raw_builder="$3"
  local full_dir="$4"
  local hvg_dir="$5"
  local stage1_full_final="$6"
  local stage1_hvg_final="$7"
  local final_hvg="$8"

  ensure_raw_full_gene "${dataset}" "${raw_path}" "${raw_builder}"

  local dataset_tmp="${TMP_ROOT}/${dataset}"
  mkdir -p "${dataset_tmp}"
  local stage1_full_tmp="${dataset_tmp}/$(basename "${stage1_full_final}")"
  local stage1_qc_tmp="${dataset_tmp}/stage1_qc"
  local stage1_hvg_tmp="${dataset_tmp}/$(basename "${stage1_hvg_final}")"
  local stage2_hvg_tmp="${dataset_tmp}/$(basename "${final_hvg}")"
  local stage2_qc_tmp="${dataset_tmp}/stage2_qc"

  echo "------------------------------------------------------------"
  echo "[${dataset}] Stage 1 marker seed labels"
  echo "------------------------------------------------------------"
  python -m benchmark.annotation.build_marker_seed_labels \
    --dataset-id "${dataset}" \
    --input-h5ad "${raw_path}" \
    --output-h5ad "${stage1_full_tmp}" \
    --output-dir "${stage1_qc_tmp}" \
    --gene-symbol-column gene_symbol \
    --use-var-index-if-needed \
    --sample-key sample_id \
    --time-key abs_day \
    --min-primary-margin "${STAGE1_MARGIN}"

  echo "------------------------------------------------------------"
  echo "[${dataset}] Transfer Stage 1 labels to HVG2000 benchmark matrix"
  echo "------------------------------------------------------------"
  if [ ! -f "${final_hvg}" ]; then
    echo "ERROR: HVG2000 base h5ad missing: ${final_hvg}" >&2
    exit 5
  fi
  transfer_stage1_to_hvg "${dataset}" "${stage1_full_tmp}" "${final_hvg}" "${stage1_hvg_tmp}"

  echo "------------------------------------------------------------"
  echo "[${dataset}] Stage 2 trajectory-aware labels"
  echo "------------------------------------------------------------"
  python -m benchmark.annotation.build_trajectory_aware_labels \
    --dataset-id "${dataset}" \
    --input-h5ad "${stage1_hvg_tmp}" \
    --output-h5ad "${stage2_hvg_tmp}" \
    --output-dir "${stage2_qc_tmp}" \
    --transition-margin "${TRANSITION_MARGIN}"

  echo "------------------------------------------------------------"
  echo "[${dataset}] Promote rebuilt h5ad/QC artifacts"
  echo "------------------------------------------------------------"
  promote_file "${stage1_full_tmp}" "${stage1_full_final}"
  promote_dir "${stage1_qc_tmp}" "${full_dir}/stage1_qc"
  promote_file "${stage1_hvg_tmp}" "${stage1_hvg_final}"
  promote_file "${stage2_hvg_tmp}" "${final_hvg}"
  promote_dir "${stage2_qc_tmp}" "${hvg_dir}/stage2_qc"

  local transfer_meta_tmp="${stage1_hvg_tmp%.h5ad}.transfer_metadata.json"
  if [ -f "${transfer_meta_tmp}" ]; then
    promote_file "${transfer_meta_tmp}" "${hvg_dir}/stage1_annotation_transfer_metadata.json"
  fi

  echo "------------------------------------------------------------"
  echo "[${dataset}] Build official_silver provider"
  echo "------------------------------------------------------------"
  local provider_dir="benchmark/ground_truth/providers/$(echo "${dataset}" | tr '[:upper:]' '[:lower:]')_marker_fm_transition_silver_v1"
  python -m benchmark.annotation.build_milestone_providers \
    --dataset-id "${dataset}" \
    --label-mode official_silver \
    --input-h5ad "${final_hvg}" \
    --dry-run
  backup_path "${provider_dir}"
  python -m benchmark.annotation.build_milestone_providers \
    --dataset-id "${dataset}" \
    --label-mode official_silver \
    --input-h5ad "${final_hvg}"
}

backup_copy "benchmark/ground_truth/registry.yaml"

rebuild_dataset \
  "GSE178325" \
  "data/processed/gse178325_human/GSE178325_0618_raw_full_gene_benchmark_input.h5ad" \
  "scripts/build_gse178325_raw_full_gene_input.py" \
  "benchmark/inputs/gse178325_marker_fm_transition_silver_full_gene" \
  "benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000" \
  "benchmark/inputs/gse178325_marker_fm_transition_silver_full_gene/GSE178325_stage1_marker_seed_full_gene.h5ad" \
  "benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_stage1_marker_seed_HVG2000.h5ad" \
  "benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"

rebuild_dataset \
  "GSE230659" \
  "data/processed/gse230659_human/GSE230659_0618_raw_full_gene_benchmark_input.h5ad" \
  "scripts/build_gse230659_raw_full_gene_input.py" \
  "benchmark/inputs/gse230659_marker_fm_transition_silver_full_gene" \
  "benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000" \
  "benchmark/inputs/gse230659_marker_fm_transition_silver_full_gene/GSE230659_stage1_marker_seed_full_gene.h5ad" \
  "benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_stage1_marker_seed_HVG2000.h5ad" \
  "benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"

echo "------------------------------------------------------------"
echo "Validate all marker_fm_silver runtime configs"
echo "------------------------------------------------------------"
mapfile -t CONFIGS < <(find benchmark/configs/runtime -maxdepth 1 -type f -name '*marker_fm_silver*_formal.yaml' | sort)
python -m benchmark.annotation.validate_milestone_configs "${CONFIGS[@]}"

echo "============================================================"
echo "Marker-FM label/provider rebuild finished at: $(date)"
echo "Backups, if any: ${BACKUP_ROOT}"
echo "Temp root      : ${TMP_ROOT}"
echo "============================================================"
