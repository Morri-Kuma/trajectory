#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_marker_fm_coarse_provider_refresh.$JOB_ID.log
#$ -l s_vmem=32G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
RUN_TAG="${RUN_TAG:-marker_fm_coarse_$(date +%Y%m%d_%H%M%S)}"
BACKUP_PROVIDER_OLD="${BACKUP_PROVIDER_OLD:-1}"

cd "${PROJECT_ROOT}"
mkdir -p logs

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"

echo "============================================================"
echo "Marker FM coarse provider refresh"
echo "Started at       : $(date)"
echo "Project root     : ${PROJECT_ROOT}"
echo "Run tag          : ${RUN_TAG}"
echo "Backup old files : ${BACKUP_PROVIDER_OLD}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

if [ "${BACKUP_PROVIDER_OLD}" = "1" ]; then
  BACKUP_ROOT="benchmark/ground_truth/_provider_backup/${RUN_TAG}"
  mkdir -p "${BACKUP_ROOT}"
  for provider in \
    gse178325_marker_fm_transition_silver_v1 \
    gse230659_marker_fm_transition_silver_v1
  do
    if [ -d "benchmark/ground_truth/providers/${provider}" ]; then
      cp -a "benchmark/ground_truth/providers/${provider}" "${BACKUP_ROOT}/${provider}"
      echo "Backed up ${provider} -> ${BACKUP_ROOT}/${provider}"
    fi
  done
fi

python benchmark/annotation/build_milestone_providers.py \
  --dataset-id GSE230659 \
  --label-mode official_silver \
  --input-h5ad benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad

python benchmark/annotation/build_milestone_providers.py \
  --dataset-id GSE178325 \
  --label-mode official_silver \
  --input-h5ad benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad

python - <<'PY'
from pathlib import Path
import json

import pandas as pd
import yaml

root = Path(".").resolve()
providers = {
    "gse230659_marker_fm_transition_silver_v1": {
        "nodes": {"hADSCs", "epithelial_like", "intermediate_plastic", "hCiPS"},
        "labels": {"hADSCs", "epithelial_like", "intermediate_plastic", "hCiPS", "ambiguous"},
    },
    "gse178325_marker_fm_transition_silver_v1": {
        "nodes": {"hADSCs", "epithelial_like", "intermediate_plastic", "xen_like", "hCiPS"},
        "labels": {
            "hADSCs", "epithelial_like", "intermediate_plastic",
            "xen_like", "hCiPS", "ambiguous",
        },
    },
}

registry = yaml.safe_load(
    (root / "benchmark/ground_truth/registry.yaml").read_text(encoding="utf-8-sig")
) or {}
registry_providers = registry.get("providers") or {}

for provider_id, expected in providers.items():
    provider_dir = root / "benchmark/ground_truth/providers" / provider_id
    graph = json.loads((provider_dir / "reference_graph.json").read_text(encoding="ascii"))
    meta = graph.get("_meta") or {}
    nodes = {str(n["id"]) for n in graph.get("nodes", [])}
    labels = set(
        pd.read_csv(provider_dir / "state_labels.tsv", sep="\t")["state_id"].astype(str)
    )
    state_meta = set(
        pd.read_csv(provider_dir / "state_metadata.tsv", sep="\t")["state_id"].astype(str)
    )
    stage = sorted(x for x in nodes | labels | state_meta if x.startswith("stage_"))
    if stage:
        raise SystemExit(f"{provider_id}: found expanded stage labels: {stage[:20]}")
    if meta.get("state_key") != "final_milestone_label_coarse":
        raise SystemExit(f"{provider_id}: reference_graph state_key={meta.get('state_key')!r}")
    if nodes != expected["nodes"]:
        raise SystemExit(f"{provider_id}: nodes={sorted(nodes)}, expected={sorted(expected['nodes'])}")
    if state_meta != expected["nodes"]:
        raise SystemExit(f"{provider_id}: state_metadata={sorted(state_meta)}")
    if labels != expected["labels"]:
        raise SystemExit(f"{provider_id}: state_labels={sorted(labels)}")
    reg = registry_providers.get(provider_id) or {}
    if reg.get("state_key") != "final_milestone_label_coarse":
        raise SystemExit(f"{provider_id}: registry state_key={reg.get('state_key')!r}")
    print(f"[OK] {provider_id}: coarse provider labels validated.")
PY

echo "============================================================"
echo "Marker FM coarse provider refresh finished at: $(date)"
echo "============================================================"
