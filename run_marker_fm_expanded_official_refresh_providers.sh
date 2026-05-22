#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_marker_fm_expanded_official_refresh_providers.$JOB_ID.log
#$ -l s_vmem=64G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
RUN_TAG="${RUN_TAG:-$(date +%Y%m%d_%H%M%S)}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-1}"
export VECLIB_MAXIMUM_THREADS="${NSLOTS:-1}"

BACKUP_ROOT="${PROJECT_ROOT}/benchmark/ground_truth/_expanded_official_backup/${RUN_TAG}"

echo "============================================================"
echo "Refresh expanded official marker-FM providers"
echo "Started at    : $(date)"
echo "Host          : $(hostname)"
echo "Project root  : ${PROJECT_ROOT}"
echo "Run tag       : ${RUN_TAG}"
echo "Backup root   : ${BACKUP_ROOT}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

backup_path() {
  local src="$1"
  local dst="${BACKUP_ROOT}/${src}"
  if [ -e "${src}" ]; then
    mkdir -p "$(dirname "${dst}")"
    cp -a "${src}" "${dst}"
    echo "[backup] ${src} -> ${dst#${PROJECT_ROOT}/}"
  fi
}

refresh_provider() {
  local dataset_id="$1"
  local input_h5ad="$2"
  local provider_id="$3"
  local expected_states="$4"
  local expected_edges="$5"

  echo "------------------------------------------------------------"
  echo "Refresh provider: ${provider_id}"
  echo "Dataset         : ${dataset_id}"
  echo "Input h5ad      : ${input_h5ad}"
  echo "Expected graph  : ${expected_states} states, ${expected_edges} edges"
  echo "------------------------------------------------------------"

  if [ ! -f "${input_h5ad}" ]; then
    echo "ERROR: input h5ad not found: ${input_h5ad}" >&2
    exit 3
  fi

  backup_path "benchmark/ground_truth/providers/${provider_id}"

  python -m benchmark.annotation.build_milestone_providers \
    --dataset-id "${dataset_id}" \
    --label-mode official_silver \
    --input-h5ad "${input_h5ad}" \
    --dry-run

  python -m benchmark.annotation.build_milestone_providers \
    --dataset-id "${dataset_id}" \
    --label-mode official_silver \
    --input-h5ad "${input_h5ad}"

  python - "${provider_id}" "${expected_states}" "${expected_edges}" <<'PY'
import csv
import json
import sys
from pathlib import Path

provider_id = sys.argv[1]
expected_states = int(sys.argv[2])
expected_edges = int(sys.argv[3])
root = Path("benchmark/ground_truth/providers") / provider_id
meta = json.loads((root / "ground_truth_metadata.json").read_text(encoding="utf-8"))
edges = list(csv.DictReader((root / "reference_graph_edges.csv").open(encoding="utf-8")))
states = (root / "reference_graph.json").read_text(encoding="utf-8")
graph = json.loads(states)

summary = {
    "provider_id": provider_id,
    "state_key": meta.get("state_key"),
    "coarse_state_key": meta.get("coarse_state_key"),
    "expanded_state_key": meta.get("expanded_state_key"),
    "n_states": meta.get("n_states"),
    "n_graph_edges": meta.get("n_graph_edges"),
    "graph_nodes": len(graph.get("nodes", [])),
    "graph_edges": len(graph.get("edges", [])),
    "graph_edges_csv": len(edges),
}
print(json.dumps(summary, indent=2))

if meta.get("state_key") != "final_milestone_label_expanded":
    raise SystemExit(f"provider {provider_id} is not expanded official")
if int(meta.get("n_states", -1)) != expected_states:
    raise SystemExit(f"unexpected n_states for {provider_id}: {meta.get('n_states')}")
if int(meta.get("n_graph_edges", -1)) != expected_edges:
    raise SystemExit(
        f"unexpected n_graph_edges for {provider_id}: {meta.get('n_graph_edges')}"
    )
if len(graph.get("edges", [])) != expected_edges:
    raise SystemExit(f"unexpected reference_graph edge count for {provider_id}: {len(graph.get('edges', []))}")
if len(edges) != expected_edges:
    raise SystemExit(f"unexpected edge table length for {provider_id}: {len(edges)}")
PY
}

backup_path "benchmark/ground_truth/registry.yaml"

refresh_provider \
  "GSE178325" \
  "benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad" \
  "gse178325_marker_fm_transition_silver_v1" \
  "9" \
  "8"

refresh_provider \
  "GSE230659" \
  "benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad" \
  "gse230659_marker_fm_transition_silver_v1" \
  "7" \
  "6"

echo "------------------------------------------------------------"
echo "Validate marker-FM runtime configs"
echo "------------------------------------------------------------"
mapfile -t CONFIGS < <(find benchmark/configs/runtime -maxdepth 1 -type f -name '*marker_fm_silver*_formal.yaml' | sort)
if [ "${#CONFIGS[@]}" -eq 0 ]; then
  echo "ERROR: no marker_fm_silver runtime configs found." >&2
  exit 4
fi
python -m benchmark.annotation.validate_milestone_configs "${CONFIGS[@]}"

echo "------------------------------------------------------------"
echo "Provider sanity summary"
echo "------------------------------------------------------------"
python - <<'PY'
import json
from pathlib import Path

provider_ids = [
    "gse178325_marker_fm_transition_silver_v1",
    "gse230659_marker_fm_transition_silver_v1",
    "gse242424_annotation_v1",
]

for provider_id in provider_ids:
    root = Path("benchmark/ground_truth/providers") / provider_id
    meta = json.loads((root / "ground_truth_metadata.json").read_text(encoding="utf-8"))
    graph = json.loads((root / "reference_graph.json").read_text(encoding="utf-8"))
    print(json.dumps({
        "provider_id": provider_id,
        "state_key": meta.get("state_key"),
        "n_states": meta.get("n_states"),
        "n_reference_edges": meta.get("n_reference_edges"),
        "graph_nodes": len(graph.get("nodes", [])),
        "graph_edges": len(graph.get("edges", [])),
    }, indent=2))
PY

echo "============================================================"
echo "Expanded official provider refresh finished at: $(date)"
echo "============================================================"
