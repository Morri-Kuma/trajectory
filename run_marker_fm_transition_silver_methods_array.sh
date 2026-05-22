#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -t 1-20
#$ -o logs/run_marker_fm_transition_silver_methods.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=128G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"
RESET_RESULT_DIR="${RESET_RESULT_DIR:-1}"
BACKUP_EXISTING_RESULTS="${BACKUP_EXISTING_RESULTS:-1}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-1}"
export VECLIB_MAXIMUM_THREADS="${NSLOTS:-1}"

RESULT_BACKUP_ROOT="${PROJECT_ROOT}/benchmark/results/_rerun_backup/${RUN_TAG}"

echo "============================================================"
echo "Marker-FM method task started at: $(date)"
echo "Host              : $(hostname)"
echo "Project root      : ${PROJECT_ROOT}"
echo "Task ID           : ${TASK_ID}"
echo "Run tag           : ${RUN_TAG}"
echo "Reset result dir  : ${RESET_RESULT_DIR}"
echo "Backup results    : ${BACKUP_EXISTING_RESULTS}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

mapfile -t CONFIGS < <(find benchmark/configs/runtime -maxdepth 1 -type f -name '*marker_fm_silver*_formal.yaml' | sort)
N_CONFIGS="${#CONFIGS[@]}"

if [ "${N_CONFIGS}" -eq 0 ]; then
  echo "ERROR: no marker_fm_silver runtime configs found." >&2
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

cfg_path = Path(sys.argv[1])
cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8-sig")) or {}
method = str(cfg.get("method", "")).lower()
scenario = str(cfg.get("scenario", ""))
dataset = str((cfg.get("dataset") or {}).get("id", ""))
h5ad = str((cfg.get("dataset") or {}).get("h5ad_path", ""))
out = str((cfg.get("output") or {}).get("base_dir", ""))
run_id = str(cfg.get("run_id", cfg_path.stem))
print(method, scenario, dataset, h5ad, out, run_id)
PY
)

echo "Config     : ${CONFIG}"
echo "Run ID     : ${RUN_ID}"
echo "Method     : ${METHOD}"
echo "Dataset    : ${DATASET}"
echo "Scenario   : ${SCENARIO}"
echo "Input h5ad : ${H5AD_PATH}"
echo "Output dir : ${OUTPUT_DIR}"

if [ ! -f "${H5AD_PATH}" ]; then
  echo "ERROR: input h5ad not found: ${H5AD_PATH}" >&2
  exit 4
fi

if [ "${RESET_RESULT_DIR}" = "1" ] && [ -d "${OUTPUT_DIR}" ]; then
  if [ "${BACKUP_EXISTING_RESULTS}" = "1" ]; then
    BACKUP_DIR="${RESULT_BACKUP_ROOT}/${OUTPUT_DIR}"
    mkdir -p "$(dirname "${BACKUP_DIR}")"
    if [ -e "${BACKUP_DIR}" ]; then
      BACKUP_DIR="${BACKUP_DIR}.$(date +%s)"
    fi
    echo "[backup] ${OUTPUT_DIR} -> ${BACKUP_DIR#${PROJECT_ROOT}/}"
    mv "${OUTPUT_DIR}" "${BACKUP_DIR}"
  else
    echo "ERROR: ${OUTPUT_DIR} exists and BACKUP_EXISTING_RESULTS is not 1." >&2
    exit 4
  fi
fi
mkdir -p "${OUTPUT_DIR}"

case "${METHOD}" in
  scnode)
    python benchmark/methods/scNODE/run.py --config "${CONFIG}"
    ;;
  mioflow)
    python benchmark/methods/MIOFlow/run.py --config "${CONFIG}"
    ;;
  prescient)
    python benchmark/methods/PRESCIENT/run.py --config "${CONFIG}"
    ;;
  wot)
    python benchmark/methods/WOT/run.py --config "${CONFIG}"
    ;;
  cellrank2)
    python -m benchmark.evaluation.eval_dispatch \
      --method cellrank2 \
      --scenario "${SCENARIO}" \
      --adata "${H5AD_PATH}" \
      --output-dir "${OUTPUT_DIR}" \
      --method-config "${CONFIG}"
    ;;
  *)
    echo "ERROR: unsupported method in ${CONFIG}: ${METHOD}" >&2
    exit 5
    ;;
esac

echo "------------------------------------------------------------"
echo "Posthoc Lineage Fidelity metrics if runner did not write them"
echo "------------------------------------------------------------"
python - "${CONFIG}" <<'PY'
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import yaml

from benchmark.evaluation.eval_lineage import run_lineage_evaluation

cfg_path = Path(sys.argv[1])
cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8-sig")) or {}
out = Path((cfg.get("output") or {}).get("base_dir", ""))
stm = out / "state_transition_matrix.csv"
edges = out / "lineage_graph_edges.csv"
metrics = out / "lineage_metrics.json"

if metrics.exists():
    print(f"lineage_metrics.json already present: {metrics}")
    raise SystemExit(0)
if not stm.exists() or not edges.exists():
    print("lineage transition outputs not present; skipping posthoc lineage metrics")
    raise SystemExit(0)

dataset_cfg = cfg.get("dataset") or {}
lineage_cfg = cfg.get("lineage") or {}
gt = cfg.get("ground_truth") or {}
h5ad_path = Path(dataset_cfg.get("h5ad_path", ""))
time_key = dataset_cfg.get("time_key") or cfg.get("time_key") or "abs_day"
cell_state_key = (
    gt.get("state_key")
    or lineage_cfg.get("cell_state_key")
    or cfg.get("cell_state_key")
    or "final_milestone_label_expanded"
)
ref = Path(lineage_cfg.get("reference_graph_path", ""))
edge_mode = gt.get("confidence_mode") or lineage_cfg.get("edge_confidence_mode") or "all"
exclude_uncertain = bool(
    gt.get("exclude_uncertain_states", lineage_cfg.get("exclude_uncertain_states", False))
)
sp = cfg.get("scenario_params") or {}
train_times = sp.get("train_times")

print(f"computing posthoc lineage metrics for {cfg.get('run_id', cfg_path.stem)}")
adata = ad.read_h5ad(h5ad_path, backed="r")
try:
    if train_times:
        train_f = {float(x) for x in train_times}
        t = adata.obs[time_key].astype(float).to_numpy()
        mask = np.array([float(x) in train_f for x in t], dtype=bool)
        print(f"scenario filter: {int(mask.sum())}/{adata.n_obs} cells")
        eval_adata = adata[mask].to_memory()
        try:
            adata.file.close()
        except Exception:
            pass
    else:
        eval_adata = adata

    out_metrics = run_lineage_evaluation(
        state_transition_matrix_path=str(stm),
        lineage_graph_edges_path=str(edges),
        output_dir=str(out),
        reference_graph_path=str(ref),
        edge_confidence_mode=edge_mode,
        exclude_uncertain_states=exclude_uncertain,
        cell_state_key=cell_state_key,
        time_key=time_key,
        adata=eval_adata,
        ground_truth=gt,
    )
    print(
        "posthoc lineage metrics written:",
        out_metrics.get("status"),
        "auroc=",
        out_metrics.get("auroc"),
    )
finally:
    try:
        if hasattr(adata, "file"):
            adata.file.close()
    except Exception:
        pass
PY

echo "------------------------------------------------------------"
echo "Validate task outputs"
echo "------------------------------------------------------------"
python - "${METHOD}" "${OUTPUT_DIR}" <<'PY'
import json
import sys
from pathlib import Path

method = sys.argv[1]
out = Path(sys.argv[2])
missing = []
for name in ["lineage_metrics.json", "run_metadata.json"]:
    if not (out / name).exists():
        missing.append(name)

if method in {"scnode", "mioflow", "prescient"}:
    for name in [
        "forecast_metrics.json",
        "embedding_metrics.json",
        "projected_expression.npy",
        "projected_embedding.npy",
        "projected_cluster_labels.csv",
    ]:
        if not (out / name).exists():
            missing.append(name)

if missing:
    raise SystemExit(f"missing expected output(s) in {out}: {missing}")

meta_path = out / "run_metadata.json"
if meta_path.exists():
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    print(json.dumps({
        "method": meta.get("method"),
        "dataset": meta.get("dataset"),
        "scenario": meta.get("scenario"),
        "status": meta.get("status"),
        "runtime_seconds": meta.get("runtime_seconds"),
    }, indent=2))
print(f"PASS: {out}")
PY

echo "============================================================"
echo "Marker-FM method task finished at: $(date)"
echo "Output dir: ${OUTPUT_DIR}"
echo "============================================================"
