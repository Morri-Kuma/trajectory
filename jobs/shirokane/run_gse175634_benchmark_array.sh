#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_gse175634_benchmark.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G
#$ -t 1-11
# NOTE: GPU is NOT requested here — Shirokane does not expose gpu as a hard
# resource flag in the same way as the other clusters. The benchmark configs
# have use_cuda: true but each method runner falls back to CPU automatically
# if CUDA is unavailable, so jobs run correctly without a gpu resource request.
# ── GSE175634 benchmark array job ────────────────────────────────────────────
#
# DATASET ID NOTE:
#   Canonical GEO accession: GSE175634 (not GSE174534, which is a typo).
#   All configs in benchmark/configs/runtime/ use GSE175634_.
#
# Array task mapping:
#   Tasks 1-3  : scnode   scenarios A, B, C
#   Tasks 4-6  : mioflow  scenarios A, B, C
#   Tasks 7-9  : prescient scenarios A, B, C
#   Task  10   : wot      scenario A   (lineage-fidelity only; GPU not used)
#   Task  11   : cellrank2 scenario A  (lineage-fidelity only; GPU not used)
#
# Silver-standard annotation step: SKIPPED.
#   GSE175634 uses author-provided cell-state labels (type column).
#   The h5ad produced by build_gse175634_cardiac_author_input.py already has
#   final_milestone_label_coarse populated from those labels.
#
# Memory / GPU notes:
#   -l s_vmem=64G covers scnode/prescient/mioflow comfortably for HVG2000.
#   -l gpu=1 enables CUDA for tasks 1-9. Tasks 10-11 ignore the GPU.
#   Adjust these resource flags for your Shirokane queue policies.
#
# Adjust the three variables below for your Shirokane environment:
#   TRAJ_PROJECT_ROOT — full path to the cloned repository on Shirokane
#   CONDA_SH          — path to conda.sh for your miniconda/anaconda install
#   CONDA_ENV         — name of the conda environment with the benchmark deps
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"

# ── task → config mapping ─────────────────────────────────────────────────────
case "${TASK_ID}" in
  1)  CONFIG="benchmark/configs/runtime/scnode_gse175634_cardiac_silver_A_hvg2000_formal.yaml" ;;
  2)  CONFIG="benchmark/configs/runtime/scnode_gse175634_cardiac_silver_B_hvg2000_formal.yaml" ;;
  3)  CONFIG="benchmark/configs/runtime/scnode_gse175634_cardiac_silver_C_hvg2000_formal.yaml" ;;
  4)  CONFIG="benchmark/configs/runtime/mioflow_gse175634_cardiac_silver_A_hvg2000_formal.yaml" ;;
  5)  CONFIG="benchmark/configs/runtime/mioflow_gse175634_cardiac_silver_B_hvg2000_formal.yaml" ;;
  6)  CONFIG="benchmark/configs/runtime/mioflow_gse175634_cardiac_silver_C_hvg2000_formal.yaml" ;;
  7)  CONFIG="benchmark/configs/runtime/prescient_gse175634_cardiac_silver_A_hvg2000_formal.yaml" ;;
  8)  CONFIG="benchmark/configs/runtime/prescient_gse175634_cardiac_silver_B_hvg2000_formal.yaml" ;;
  9)  CONFIG="benchmark/configs/runtime/prescient_gse175634_cardiac_silver_C_hvg2000_formal.yaml" ;;
  10) CONFIG="benchmark/configs/runtime/wot_gse175634_cardiac_silver_A_hvg2000_formal.yaml" ;;
  11) CONFIG="benchmark/configs/runtime/cellrank2_gse175634_cardiac_silver_A_hvg2000_formal.yaml" ;;
  *)
    echo "ERROR: unsupported task id ${TASK_ID}; expected 1-11" >&2
    exit 2
    ;;
esac

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-1}"
export VECLIB_MAXIMUM_THREADS="${NSLOTS:-1}"

# Parse config fields for logging
read -r METHOD SCENARIO DATASET H5AD_PATH OUTPUT_DIR RUN_ID < <(python - "${CONFIG}" <<'PY'
import sys
from pathlib import Path
import yaml

cfg = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8-sig")) or {}
method   = str(cfg.get("method", "")).lower()
scenario = str(cfg.get("scenario", ""))
dataset  = str((cfg.get("dataset") or {}).get("id", ""))
h5ad     = str((cfg.get("dataset") or {}).get("h5ad_path", ""))
out      = str((cfg.get("output") or {}).get("base_dir", ""))
run_id   = str(cfg.get("run_id", Path(sys.argv[1]).stem))
print(method, scenario, dataset, h5ad, out, run_id)
PY
)

echo "============================================================"
echo "GSE175634 benchmark array"
echo "Started at   : $(date)"
echo "Host         : $(hostname)"
echo "Task ID      : ${TASK_ID}"
echo "Config       : ${CONFIG}"
echo "Run ID       : ${RUN_ID}"
echo "Dataset      : ${DATASET}"
echo "Method       : ${METHOD}"
echo "Scenario     : ${SCENARIO}"
echo "Input h5ad   : ${H5AD_PATH}"
echo "Output dir   : ${OUTPUT_DIR}"
echo "Silver step  : SKIPPED (author-provided labels used directly)"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

# ── pre-flight checks ─────────────────────────────────────────────────────────
if [ ! -f "${H5AD_PATH}" ]; then
  echo "ERROR: benchmark h5ad not found: ${H5AD_PATH}" >&2
  echo "Run jobs/shirokane/run_gse175634_preprocess.sh first." >&2
  exit 4
fi

PROVIDER_GRAPH="benchmark/ground_truth/providers/gse175634_cardiac_silver_v1/reference_graph.json"
if [ ! -f "${PROVIDER_GRAPH}" ]; then
  echo "ERROR: provider reference_graph not found: ${PROVIDER_GRAPH}" >&2
  echo "Run jobs/shirokane/run_gse175634_preprocess.sh (Step 2) first." >&2
  exit 4
fi

mkdir -p "${OUTPUT_DIR}"

# ── dispatch to method runner ─────────────────────────────────────────────────
echo ""
echo "------------------------------------------------------------"
echo "Running method: ${METHOD}  scenario: ${SCENARIO}"
echo "------------------------------------------------------------"

case "${METHOD}" in
  scnode|prescient|mioflow|wot|cellrank2)
    python benchmark/evaluation/eval_dispatch.py \
      --method "${METHOD}" \
      --scenario "${SCENARIO}" \
      --adata "${H5AD_PATH}" \
      --output-dir "${OUTPUT_DIR}" \
      --method-config "${CONFIG}"
    ;;
  *)
    echo "ERROR: unknown method '${METHOD}' in config ${CONFIG}" >&2
    exit 5
    ;;
esac

# ── verify outputs ────────────────────────────────────────────────────────────
echo ""
echo "------------------------------------------------------------"
echo "Verifying outputs …"
echo "------------------------------------------------------------"
python - "${OUTPUT_DIR}" "${METHOD}" "${SCENARIO}" "${RUN_ID}" <<'PY'
import json, sys
from pathlib import Path

out_dir  = Path(sys.argv[1])
method   = sys.argv[2]
scenario = sys.argv[3]
run_id   = sys.argv[4]

run_meta = out_dir / "run_metadata.json"
if not run_meta.exists():
    raise SystemExit(f"run_metadata.json missing in {out_dir}")
meta = json.loads(run_meta.read_text(encoding="utf-8"))
if meta.get("status") not in {"completed", "ok"}:
    raise SystemExit(f"run_metadata status not completed: {meta.get('status')}")
print(f"[OK] run_metadata.json: status={meta.get('status')!r}")

# Lineage output (all methods produce this)
stm = out_dir / "state_transition_matrix.csv"
if not stm.exists():
    raise SystemExit(f"state_transition_matrix.csv missing in {out_dir}")
print(f"[OK] state_transition_matrix.csv present")

metrics = out_dir / "lineage_metrics.json"
if not metrics.exists():
    raise SystemExit(f"lineage_metrics.json missing in {out_dir}")
lm = json.loads(metrics.read_text(encoding="utf-8"))
print(f"[OK] lineage_metrics.json: {list(lm.keys())}")
PY

echo "============================================================"
echo "Task ${TASK_ID} DONE at: $(date)"
echo "Output: ${OUTPUT_DIR}"
echo "============================================================"
