#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_gse242424_ot_lineage.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=192G
#$ -t 1-2

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"
RUN_TAG="${RUN_TAG:-gse242424_ot_lineage_$(date +%Y%m%d_%H%M%S)}"
BACKUP_OLD_LINEAGE="${BACKUP_OLD_LINEAGE:-1}"
PURGE_METHOD_CACHE="${PURGE_METHOD_CACHE:-0}"

CONFIGS=(
  "benchmark/configs/wot_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml"
  "benchmark/configs/cellrank2_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml"
)

if [ "${TASK_ID}" -lt 1 ] || [ "${TASK_ID}" -gt "${#CONFIGS[@]}" ]; then
  echo "ERROR: unsupported task id ${TASK_ID}; expected 1-${#CONFIGS[@]}" >&2
  exit 2
fi

CONFIG="${CONFIGS[$((TASK_ID - 1))]}"

cd "${PROJECT_ROOT}"
mkdir -p logs

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export CONFIG
export RUN_TAG
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-1}"
export VECLIB_MAXIMUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "GSE242424 OT lineage run"
echo "Started at        : $(date)"
echo "Host              : $(hostname)"
echo "Task ID           : ${TASK_ID}/${#CONFIGS[@]}"
echo "Config            : ${CONFIG}"
echo "Project root      : ${PROJECT_ROOT}"
echo "Run tag           : ${RUN_TAG}"
echo "Backup old lineage: ${BACKUP_OLD_LINEAGE}"
echo "Purge method cache: ${PURGE_METHOD_CACHE}"
echo "============================================================"

if [ ! -f "${CONFIG}" ]; then
  echo "ERROR: config not found: ${CONFIG}" >&2
  exit 3
fi
if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

read -r METHOD SCENARIO DATASET_ID H5AD_PATH OUTPUT_DIR CELL_STATE_KEY REF_GRAPH < <(
  python - "${CONFIG}" <<'PY'
from pathlib import Path
import sys

import yaml

cfg = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8-sig")) or {}
method = str(cfg.get("method", "")).lower()
scenario = str(cfg.get("scenario", ""))
dataset = cfg.get("dataset") or {}
lineage = cfg.get("lineage") or {}
ground_truth = cfg.get("ground_truth") or {}
state_system = cfg.get("state_system") or {}
h5ad = dataset.get("h5ad_path") or state_system.get("source_h5ad") or ""
output = (cfg.get("output") or {}).get("base_dir") or ""
cell_state_key = (
    lineage.get("cell_state_key")
    or ground_truth.get("state_key")
    or cfg.get("cell_state_key")
    or "final_milestone_label_coarse"
)
ref_graph = lineage.get("reference_graph_path") or ground_truth.get("reference_graph_path") or ""
print(method, scenario, dataset.get("id", ""), h5ad, output, cell_state_key, ref_graph)
PY
)

echo "Method        : ${METHOD}"
echo "Scenario      : ${SCENARIO}"
echo "Dataset       : ${DATASET_ID}"
echo "Input h5ad    : ${H5AD_PATH}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Cell state key: ${CELL_STATE_KEY}"
echo "Reference     : ${REF_GRAPH}"

if [ "${METHOD}" != "wot" ] && [ "${METHOD}" != "cellrank2" ]; then
  echo "ERROR: this array only supports wot/cellrank2, got ${METHOD}" >&2
  exit 3
fi
if [ "${DATASET_ID}" != "GSE242424" ] || [ "${SCENARIO}" != "A" ]; then
  echo "ERROR: this array is only for GSE242424 scenario A, got ${DATASET_ID}/${SCENARIO}" >&2
  exit 3
fi
if [ "${CELL_STATE_KEY}" != "final_milestone_label_coarse" ]; then
  echo "ERROR: config does not use final_milestone_label_coarse: ${CONFIG}" >&2
  exit 3
fi
if [ ! -f "${H5AD_PATH}" ]; then
  echo "ERROR: benchmark h5ad not found: ${H5AD_PATH}" >&2
  exit 3
fi
if [ ! -f "${REF_GRAPH}" ]; then
  echo "ERROR: reference graph not found: ${REF_GRAPH}" >&2
  exit 3
fi

python - "${CONFIG}" <<'PY'
from pathlib import Path
import json
import sys

import anndata as ad
import yaml

config = Path(sys.argv[1])
cfg = yaml.safe_load(config.read_text(encoding="utf-8-sig")) or {}
dataset = cfg.get("dataset") or {}
lineage = cfg.get("lineage") or {}
h5ad_path = Path(dataset.get("h5ad_path"))
ref_graph = Path(lineage.get("reference_graph_path"))
cell_state_key = lineage.get("cell_state_key")

adata = ad.read_h5ad(str(h5ad_path), backed="r")
try:
    if cell_state_key != "final_milestone_label_coarse":
        raise SystemExit(f"cell_state_key={cell_state_key!r}")
    if cell_state_key not in adata.obs.columns:
        raise SystemExit(f"{cell_state_key!r} not found in h5ad obs")
    labels = {str(x) for x in adata.obs[cell_state_key].dropna().unique()}
    if not labels:
        raise SystemExit(f"{cell_state_key!r} has no labels")
finally:
    adata.file.close()

graph = json.loads(ref_graph.read_text(encoding="utf-8"))
nodes = {str(n["id"]) for n in graph.get("nodes", [])}
meta = graph.get("_meta") or {}
if meta.get("state_key") != "final_milestone_label_coarse":
    raise SystemExit(f"reference graph state_key={meta.get('state_key')!r}")
missing = sorted(nodes - labels)
print(
    f"[OK] preflight labels: {config.name}, "
    f"n_h5ad_labels={len(labels)}, n_graph_nodes={len(nodes)}, "
    f"graph_nodes_missing_from_h5ad={missing}"
)
PY

if [ "${BACKUP_OLD_LINEAGE}" = "1" ]; then
  python - "${OUTPUT_DIR}" "${RUN_TAG}" "${PURGE_METHOD_CACHE}" <<'PY'
from pathlib import Path
import shutil
import sys

out_dir = Path(sys.argv[1])
run_tag = sys.argv[2]
purge_cache = sys.argv[3] == "1"
backup_dir = out_dir / "_rerun_backup" / run_tag
backup_dir.mkdir(parents=True, exist_ok=True)

patterns = [
    "state_transition_matrix.csv",
    "lineage_graph_edges.csv",
    "lineage_metrics.json",
    "lineage_diagnostics.json",
    "baseline_state_transition_matrix.csv",
    "baseline_lineage_graph_edges.csv",
    "lineage_metrics.pre*.json",
    "run_metadata.json",
]
if purge_cache:
    patterns.extend(["tmaps", "tmaps_cr2", "tmaps*", "*.pt", "*.pth"])

for pattern in patterns:
    for path in out_dir.glob(pattern):
        if not path.exists():
            continue
        dest = backup_dir / path.name
        if dest.exists():
            continue
        if path.is_dir():
            shutil.move(str(path), str(dest))
        else:
            shutil.move(str(path), str(dest))
        print(f"[backup] {path} -> {dest}")
PY
fi

python benchmark/evaluation/eval_dispatch.py \
  --method "${METHOD}" \
  --scenario "${SCENARIO}" \
  --adata "${H5AD_PATH}" \
  --output-dir "${OUTPUT_DIR}" \
  --method-config "${CONFIG}"

python - "${CONFIG}" <<'PY'
from pathlib import Path
import json
import sys

import pandas as pd
import yaml

cfg = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8-sig")) or {}
out_dir = Path((cfg.get("output") or {}).get("base_dir"))
ref_graph = Path((cfg.get("lineage") or {}).get("reference_graph_path"))
graph = json.loads(ref_graph.read_text(encoding="utf-8"))
ref_nodes = [str(n["id"]) for n in graph.get("nodes", [])]

stm_path = out_dir / "state_transition_matrix.csv"
edges_path = out_dir / "lineage_graph_edges.csv"
metrics_path = out_dir / "lineage_metrics.json"
run_meta_path = out_dir / "run_metadata.json"
for path in (stm_path, edges_path, metrics_path, run_meta_path):
    if not path.exists():
        raise SystemExit(f"missing required output: {path}")

stm = pd.read_csv(stm_path, index_col=0)
rows = [str(x) for x in stm.index]
cols = [str(x) for x in stm.columns]
missing_rows = [x for x in ref_nodes if x not in set(rows)]
missing_cols = [x for x in ref_nodes if x not in set(cols)]

metrics = json.loads(metrics_path.read_text(encoding="utf-8-sig"))
if metrics.get("status") != "completed":
    raise SystemExit(f"lineage_metrics status={metrics.get('status')!r}")
if metrics.get("cell_state_key") != "final_milestone_label_coarse":
    raise SystemExit(f"cell_state_key={metrics.get('cell_state_key')!r}")
if metrics.get("prediction_state_key") != "final_milestone_label_coarse":
    raise SystemExit(f"prediction_state_key={metrics.get('prediction_state_key')!r}")
expected_label_source = "reference_graph_nodes_sctimebench_zero_filled"
if metrics.get("prediction_label_source") != expected_label_source:
    raise SystemExit(f"prediction_label_source={metrics.get('prediction_label_source')!r}")

report = metrics.get("prediction_label_report") or {}
reported_missing = set(report.get("missing_reference_nodes") or [])
expected_missing = set(missing_rows) | set(missing_cols)
if expected_missing != reported_missing:
    raise SystemExit(
        "prediction_label_report missing_reference_nodes mismatch: "
        f"expected={sorted(expected_missing)}, reported={sorted(reported_missing)}"
    )

run_meta = json.loads(run_meta_path.read_text(encoding="utf-8-sig"))
if run_meta.get("status") != "completed":
    raise SystemExit(f"run_metadata status={run_meta.get('status')!r}")

single = (metrics.get("graph_metrics") or {}).get("single_step") or {}
print(
    "[OK] GSE242424 OT lineage completed: "
    f"{out_dir}; AUROC={single.get('auc_roc')}; AUPRC={single.get('auc_prc')}; "
    f"missing_reference_nodes={sorted(reported_missing)}"
)
PY

echo "============================================================"
echo "GSE242424 OT lineage run finished at: $(date)"
echo "Config: ${CONFIG}"
echo "============================================================"
