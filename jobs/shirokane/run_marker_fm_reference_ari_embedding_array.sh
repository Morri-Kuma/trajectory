#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -t 1-18
#$ -o logs/run_marker_fm_reference_ari_embedding.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"
LEIDEN_RESOLUTION="${LEIDEN_RESOLUTION:-0.5}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-1}"
export VECLIB_MAXIMUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Marker-FM reference-ARI embedding post-evaluation"
echo "Started at   : $(date)"
echo "Host         : $(hostname)"
echo "Project root : ${PROJECT_ROOT}"
echo "Task ID      : ${TASK_ID}"
echo "Resolution   : ${LEIDEN_RESOLUTION}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

mapfile -t CONFIGS < <(python - <<'PY'
import glob
from pathlib import Path
import yaml

for p in sorted(glob.glob("benchmark/configs/runtime/*marker_fm_silver*_formal.yaml")):
    cfg = yaml.safe_load(Path(p).read_text(encoding="utf-8-sig")) or {}
    if str(cfg.get("method", "")).lower() in {"scnode", "mioflow", "prescient"}:
        print(p)
PY
)

N_CONFIGS="${#CONFIGS[@]}"
if [ "${N_CONFIGS}" -eq 0 ]; then
  echo "ERROR: no projection-capable marker_fm_silver configs found." >&2
  exit 3
fi
if [ "${TASK_ID}" -lt 1 ] || [ "${TASK_ID}" -gt "${N_CONFIGS}" ]; then
  echo "Task ${TASK_ID} is outside config range 1-${N_CONFIGS}; nothing to do."
  exit 0
fi

CONFIG="${CONFIGS[$((TASK_ID - 1))]}"
read -r METHOD SCENARIO DATASET H5AD_PATH OUTPUT_DIR RUN_ID < <(python - "${CONFIG}" <<'PY'
import sys
from pathlib import Path
import yaml

cfg = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8-sig")) or {}
method = str(cfg.get("method", "")).lower()
scenario = str(cfg.get("scenario", ""))
dataset = str((cfg.get("dataset") or {}).get("id", ""))
h5ad = str((cfg.get("dataset") or {}).get("h5ad_path", ""))
out = str((cfg.get("output") or {}).get("base_dir", ""))
run_id = str(cfg.get("run_id", Path(sys.argv[1]).stem))
print(method, scenario, dataset, h5ad, out, run_id)
PY
)

EMBEDDING_OUT="${OUTPUT_DIR}/embedding_milestone_eval"

echo "Config        : ${CONFIG}"
echo "Run ID        : ${RUN_ID}"
echo "Dataset       : ${DATASET}"
echo "Method        : ${METHOD}"
echo "Scenario      : ${SCENARIO}"
echo "Input h5ad    : ${H5AD_PATH}"
echo "Run output    : ${OUTPUT_DIR}"
echo "Embedding out : ${EMBEDDING_OUT}"

if [ ! -f "${H5AD_PATH}" ]; then
  echo "ERROR: reference h5ad missing: ${H5AD_PATH}" >&2
  exit 4
fi
if [ ! -d "${OUTPUT_DIR}" ]; then
  echo "ERROR: method output dir missing: ${OUTPUT_DIR}" >&2
  exit 4
fi
if [ ! -f "${OUTPUT_DIR}/embedding.npy" ]; then
  echo "ERROR: embedding.npy missing in ${OUTPUT_DIR}; rerun the method task with scTimeBench embedding outputs." >&2
  exit 4
fi
if [ ! -f "${OUTPUT_DIR}/next_timepoint_embedding.npy" ] && [ ! -f "${OUTPUT_DIR}/projected_embedding.npy" ]; then
  echo "ERROR: next_timepoint_embedding.npy/projected_embedding.npy missing in ${OUTPUT_DIR}" >&2
  exit 4
fi

mkdir -p "${EMBEDDING_OUT}"

echo "------------------------------------------------------------"
echo "Run official_silver Embedding Coherence with observed reference ARI"
echo "------------------------------------------------------------"
python -m benchmark.evaluation.eval_embedding_milestone \
  --input-h5ad "${H5AD_PATH}" \
  --run-output-dir "${OUTPUT_DIR}" \
  --output-dir "${EMBEDDING_OUT}" \
  --dataset-id "${DATASET}" \
  --mode run-output \
  --label-mode official_silver \
  --exclude-label ambiguous \
  --exclude-label unknown_or_ood \
  --leiden-n-neighbors 15 \
  --leiden-resolution "${LEIDEN_RESOLUTION}"

echo "------------------------------------------------------------"
echo "Validate reference-ARI fields"
echo "------------------------------------------------------------"
python - "${EMBEDDING_OUT}" "${LEIDEN_RESOLUTION}" <<'PY'
import json
import sys
from pathlib import Path

out = Path(sys.argv[1])
resolution = sys.argv[2]
metrics_path = out / "embedding_metrics_official_silver.json"
if not metrics_path.exists():
    raise SystemExit(f"missing {metrics_path}")

m = json.loads(metrics_path.read_text(encoding="utf-8"))
required = [
    "adjusted_rand_index",
    "reference_adjusted_rand_index",
    "reference_ari_status",
    "ari_retention_fraction",
    "pred_tp_avg_normalized_entropy",
    "entropy_basis",
    "state_key",
]
missing = [k for k in required if k not in m]
if missing:
    raise SystemExit(f"missing reference-ARI field(s): {missing}")
if m.get("state_key") != "final_milestone_label_coarse":
    raise SystemExit(f"unexpected state_key: {m.get('state_key')}")
if m.get("reference_ari_status") != "ok":
    raise SystemExit(f"reference ARI not ok: {m.get('reference_ari_status')}")
if m.get("reference_adjusted_rand_index") is None:
    raise SystemExit("reference_adjusted_rand_index is null")
if m.get("embedding_metric_protocol") != "sctimebench":
    raise SystemExit(f"unexpected embedding_metric_protocol: {m.get('embedding_metric_protocol')}")
if m.get("entropy_basis") != "classifier_probability_vector":
    raise SystemExit(f"unexpected entropy_basis: {m.get('entropy_basis')}")
if m.get("label_source") != "knn_transfer_from_observed_embedding":
    raise SystemExit(f"unexpected label_source: {m.get('label_source')}")
for key in ("cluster_source", "reference_cluster_source"):
    source = str(m.get(key, ""))
    expected = f"resolution={resolution}"
    if expected not in source:
        raise SystemExit(f"{key} does not include {expected!r}: {source}")

print(json.dumps({
    "status": m.get("status"),
    "state_key": m.get("state_key"),
    "cluster_source": m.get("cluster_source"),
    "reference_cluster_source": m.get("reference_cluster_source"),
    "projected_ari": m.get("adjusted_rand_index"),
    "reference_ari": m.get("reference_adjusted_rand_index"),
    "ari_retention_fraction": m.get("ari_retention_fraction"),
    "pred_tp_avg_normalized_entropy": m.get("pred_tp_avg_normalized_entropy"),
    "n_cells_evaluated": m.get("n_cells_evaluated"),
    "reference_n_cells_evaluated": m.get("reference_n_cells_evaluated"),
}, indent=2))
PY

echo "============================================================"
echo "Reference-ARI embedding post-evaluation finished at: $(date)"
echo "Output: ${EMBEDDING_OUT}/embedding_metrics_official_silver.json"
echo "============================================================"
