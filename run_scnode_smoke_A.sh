#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scnode_smoke_A.$JOB_ID.log
#$ -l s_vmem=64G

set -euo pipefail

# ===== user-configurable section =====
PROJECT_ROOT="/home/xzy0723/projects/trajectory"
SCGPT_REPO="/home/xzy0723/projects/scGPT"
CONDA_SH="/home/xzy0723/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="traj_env"

METHOD="scnode"
SCENARIO="A"
BASE_CONFIG="benchmark/configs/scnode_gse230659_observed_scgpt_v1.yaml"
SMOKE_CONFIG="benchmark/configs/scnode_gse230659_observed_scgpt_v1_smoke.yaml"
SOURCE_ADATA="benchmark/results/scgpt/full/adata_scgpt_annotated.h5ad"
SMOKE_ADATA="benchmark/results/scnode/smoke_inputs/adata_scnode_smoke_A.h5ad"
OUTPUT_DIR="benchmark/results/scnode/scenario_A_scgpt_v1_smoke"

# Smoke-test size and model settings. Increase CELLS_PER_TIMEPOINT only after
# the adapter path has been validated.
CELLS_PER_TIMEPOINT=100
N_SIM_CELLS=100
PRETRAIN_ITERS=5
EPOCHS=1
ITERS=2
BATCH_SIZE=32
LATENT_DIM=10
# =====================================

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export SCGPT_REPO="${SCGPT_REPO}"

echo "============================================================"
echo "Job started at: $(date)"
echo "Host          : $(hostname)"
echo "Project root  : ${TRAJ_PROJECT_ROOT}"
echo "SCGPT repo    : ${SCGPT_REPO}"
echo "Method        : ${METHOD}"
echo "Scenario      : ${SCENARIO} smoke"
echo "Source adata  : ${SOURCE_ADATA}"
echo "Smoke adata   : ${SMOKE_ADATA}"
echo "Base config   : ${BASE_CONFIG}"
echo "Smoke config  : ${SMOKE_CONFIG}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Job ID        : ${JOB_ID:-N/A}"
echo "============================================================"

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

echo "Python        : $(which python)"
python -V

python - <<'PY'
import anndata
import torch
import torchdiffeq
print("imports OK")
print("anndata  :", anndata.__version__)
print("torch    :", torch.__version__)
print("torchdiffeq import OK")
PY

echo "------------------------------------------------------------"
echo "Building smoke-test h5ad and config"
echo "------------------------------------------------------------"

python - <<PY
from pathlib import Path
import yaml
import anndata as ad

root = Path("${PROJECT_ROOT}")
source_adata = root / "${SOURCE_ADATA}"
smoke_adata = root / "${SMOKE_ADATA}"
base_config = root / "${BASE_CONFIG}"
smoke_config = root / "${SMOKE_CONFIG}"

cells_per_timepoint = int("${CELLS_PER_TIMEPOINT}")
n_sim_cells = int("${N_SIM_CELLS}")

smoke_adata.parent.mkdir(parents=True, exist_ok=True)
smoke_config.parent.mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(source_adata, backed="r")
selected = []
for _, obs in adata.obs.groupby("abs_day", observed=True):
    selected.extend(obs.index[:cells_per_timepoint])

subset = adata[selected].to_memory()
subset.write_h5ad(smoke_adata)
print(f"wrote smoke adata: {smoke_adata} shape={subset.shape}")

with open(base_config, encoding="utf-8-sig") as f:
    cfg = yaml.safe_load(f)

cfg["run_id"] = "scnode_gse230659_A_scgpt_v1_smoke"
cfg["scenario"] = "A"
cfg["dataset"]["id"] = "GSE230659"
cfg["dataset"]["h5ad_path"] = str(smoke_adata.relative_to(root))
cfg["dataset"]["time_key"] = "abs_day"
cfg["output"]["base_dir"] = "${OUTPUT_DIR}"
cfg["scnode_params"].update({
    "latent_dim": int("${LATENT_DIM}"),
    "drift_latent_size": [16, 16],
    "enc_latent_list": [32],
    "dec_latent_list": [32],
    "pretrain_iters": int("${PRETRAIN_ITERS}"),
    "epochs": int("${EPOCHS}"),
    "iters": int("${ITERS}"),
    "batch_size": int("${BATCH_SIZE}"),
    "n_sim_cells": n_sim_cells,
    "n_sim_cells_cap": n_sim_cells,
    "seed": 42,
})

# Scenario A smoke uses all time points present in the smoke h5ad.
cfg.pop("scenario_params", None)

with open(smoke_config, "w", encoding="utf-8") as f:
    yaml.safe_dump(cfg, f, sort_keys=False)
print(f"wrote smoke config: {smoke_config}")
PY

CMD=(
  python -m benchmark.evaluation.eval_dispatch
  --method "${METHOD}"
  --scenario "${SCENARIO}"
  --adata "${SMOKE_ADATA}"
  --output-dir "${OUTPUT_DIR}"
  --method-config "${SMOKE_CONFIG}"
)

echo "------------------------------------------------------------"
echo "Running command:"
printf ' %q' "${CMD[@]}"
echo
echo "------------------------------------------------------------"

time "${CMD[@]}"

echo "------------------------------------------------------------"
echo "Output check"
echo "------------------------------------------------------------"

python - <<PY
from pathlib import Path
import json

out = Path("${PROJECT_ROOT}") / "${OUTPUT_DIR}"
required = [
    "projected_expression.npy",
    "forecast_metrics.json",
    "per_timepoint_forecast_metrics.csv",
    "projected_embedding.npy",
    "embedding_metrics.json",
    "projected_cluster_labels.csv",
    "state_transition_matrix.csv",
    "lineage_graph_edges.csv",
    "lineage_metrics.json",
    "run_metadata.json",
]
missing = [name for name in required if not (out / name).exists()]
print(f"output_dir: {out}")
print("missing:", missing)
if missing:
    raise SystemExit(1)

for name in ["run_metadata.json", "forecast_metrics.json", "embedding_metrics.json", "lineage_metrics.json"]:
    print(f"\\n== {name} ==")
    with open(out / name, encoding="utf-8") as f:
        data = json.load(f)
    print(json.dumps(data, indent=2)[:2000])
PY

echo "============================================================"
echo "Job finished at: $(date)"
echo "============================================================"
