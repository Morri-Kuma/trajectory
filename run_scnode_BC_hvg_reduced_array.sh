#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scnode_BC_hvg_reduced.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=128G
#$ -t 1-2

set -euo pipefail

# ===== user-configurable section =====
PROJECT_ROOT="/home/xzy0723/projects/trajectory"
SCGPT_REPO="/home/xzy0723/projects/scGPT"
CONDA_SH="/home/xzy0723/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="traj_env"

METHOD="scnode"
HVG_ADATA="benchmark/inputs/gse230659_scgpt_hvg2000/GSE230659_scGPT_annotated_HVG2000_benchmark_input.h5ad"

N_SIM_CELLS=1000
PRETRAIN_ITERS=50
EPOCHS=3
ITERS=20
BATCH_SIZE=32
LATENT_DIM=50
# =====================================

case "${SGE_TASK_ID:-${TASK_ID:-1}}" in
  1)
    SCENARIO="B"
    BASE_CONFIG="benchmark/configs/scnode_gse230659_observed_scgpt_v1_scenarioB.yaml"
    RUN_ID="scnode_gse230659_B_scgpt_v1_hvg2000_reduced"
    REDUCED_CONFIG="benchmark/configs/scnode_gse230659_observed_scgpt_v1_scenarioB_hvg2000_reduced.yaml"
    OUTPUT_DIR="benchmark/results/scnode/scenario_B_scgpt_v1_hvg2000_reduced"
    ;;
  2)
    SCENARIO="C"
    BASE_CONFIG="benchmark/configs/scnode_gse230659_observed_scgpt_v1_scenarioC.yaml"
    RUN_ID="scnode_gse230659_C_scgpt_v1_hvg2000_reduced"
    REDUCED_CONFIG="benchmark/configs/scnode_gse230659_observed_scgpt_v1_scenarioC_hvg2000_reduced.yaml"
    OUTPUT_DIR="benchmark/results/scnode/scenario_C_scgpt_v1_hvg2000_reduced"
    ;;
  *)
    echo "Unsupported SGE_TASK_ID=${SGE_TASK_ID:-${TASK_ID:-unset}}; expected 1 or 2" >&2
    exit 2
    ;;
esac

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
echo "Scenario      : ${SCENARIO} HVG2000 reduced-training"
echo "HVG adata     : ${HVG_ADATA}"
echo "Base config   : ${BASE_CONFIG}"
echo "Reduced config: ${REDUCED_CONFIG}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Job ID        : ${JOB_ID:-N/A}"
echo "Task ID       : ${SGE_TASK_ID:-${TASK_ID:-N/A}}"
echo "============================================================"

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

echo "Python        : $(which python)"
python -V

python - <<'PY'
import anndata
import scanpy
import scipy
import torch
import torchdiffeq
print("imports OK")
print("anndata  :", anndata.__version__)
print("scanpy   :", scanpy.__version__)
print("scipy    :", scipy.__version__)
print("torch    :", torch.__version__)
print("torchdiffeq import OK")
PY

echo "------------------------------------------------------------"
echo "Preparing Scenario ${SCENARIO} reduced config"
echo "------------------------------------------------------------"

python - <<PY
from pathlib import Path
import yaml
import anndata as ad

root = Path("${PROJECT_ROOT}")
hvg_adata = root / "${HVG_ADATA}"
base_config = root / "${BASE_CONFIG}"
reduced_config = root / "${REDUCED_CONFIG}"

if not hvg_adata.exists():
    raise FileNotFoundError(
        f"Benchmark HVG adata not found: {hvg_adata}. "
        "Run/sync run_scnode_full_A_hvg.sh first."
    )
if not base_config.exists():
    raise FileNotFoundError(f"Base config not found: {base_config}")

adata_meta = ad.read_h5ad(hvg_adata, backed="r")
print(f"HVG adata shape={adata_meta.shape}")
print("timepoints:", sorted(adata_meta.obs["abs_day"].astype(float).unique()))
print("states:", len(set(adata_meta.obs["scgpt_pseudostate_provisional"].astype(str))))

with open(base_config, encoding="utf-8-sig") as f:
    cfg = yaml.safe_load(f)

cfg["run_id"] = "${RUN_ID}"
cfg["scenario"] = "${SCENARIO}"
cfg["dataset"]["id"] = "GSE230659"
cfg["dataset"]["h5ad_path"] = "${HVG_ADATA}"
cfg["dataset"]["time_key"] = "abs_day"
cfg["output"]["base_dir"] = "${OUTPUT_DIR}"
cfg["ground_truth"] = {
    "provider_id": "scgpt_v1",
    "state_key": "scgpt_pseudostate_provisional",
    "confidence_mode": "medium_and_above",
    "exclude_uncertain_states": False,
}
cfg["scnode_params"].update({
    "latent_dim": int("${LATENT_DIM}"),
    "drift_latent_size": [50, 50],
    "enc_latent_list": [64, 64],
    "dec_latent_list": [64, 64],
    "pretrain_iters": int("${PRETRAIN_ITERS}"),
    "epochs": int("${EPOCHS}"),
    "iters": int("${ITERS}"),
    "batch_size": int("${BATCH_SIZE}"),
    "n_sim_cells": int("${N_SIM_CELLS}"),
    "n_sim_cells_cap": int("${N_SIM_CELLS}"),
    "seed": 42,
})

reduced_config.parent.mkdir(parents=True, exist_ok=True)
with open(reduced_config, "w", encoding="utf-8") as f:
    yaml.safe_dump(cfg, f, sort_keys=False)
print(f"wrote reduced config: {reduced_config}")
print("train_times:", cfg.get("scenario_params", {}).get("train_times"))
print("heldout_times:", cfg.get("scenario_params", {}).get("heldout_times"))
PY

CMD=(
  python -m benchmark.evaluation.eval_dispatch
  --method "${METHOD}"
  --scenario "${SCENARIO}"
  --adata "${HVG_ADATA}"
  --output-dir "${OUTPUT_DIR}"
  --method-config "${REDUCED_CONFIG}"
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
import numpy as np
import pandas as pd

out = Path("${PROJECT_ROOT}") / "${OUTPUT_DIR}"
required = [
    "trained_scnode_model.pth",
    "projected_expression.npy",
    "forecast_metrics.json",
    "per_timepoint_forecast_metrics.csv",
    "projected_embedding.npy",
    "embedding_metrics.json",
    "projected_cluster_labels.csv",
    "state_transition_matrix.csv",
    "lineage_graph_edges.csv",
    "lineage_diagnostics.json",
    "lineage_metrics.json",
    "run_metadata.json",
]
missing = [name for name in required if not (out / name).exists()]
print(f"output_dir: {out}")
print("missing:", missing)
if missing:
    raise SystemExit(1)

stm = pd.read_csv(out / "state_transition_matrix.csv", index_col=0)
vals = stm.values.astype(float)
print("stm_shape:", stm.shape)
print("stm_nonfinite:", int((~np.isfinite(vals)).sum()))
print("stm_zero_rows:", stm.index[np.isclose(vals.sum(axis=1), 0)].tolist())
if not np.isfinite(vals).all():
    raise SystemExit(1)

for name in [
    "run_metadata.json",
    "forecast_metrics.json",
    "embedding_metrics.json",
    "lineage_diagnostics.json",
    "lineage_metrics.json",
]:
    print(f"\\n== {name} ==")
    with open(out / name, encoding="utf-8") as f:
        data = json.load(f)
    print(json.dumps(data, indent=2)[:2400])
PY

echo "============================================================"
echo "Job finished at: $(date)"
echo "============================================================"
