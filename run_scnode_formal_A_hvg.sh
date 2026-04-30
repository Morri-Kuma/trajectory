#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scnode_formal_A_hvg.$JOB_ID.log
#$ -l s_vmem=128G

set -euo pipefail

PROJECT_ROOT="/home/xzy0723/projects/trajectory"
SCGPT_REPO="/home/xzy0723/projects/scGPT"
CONDA_SH="/home/xzy0723/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="traj_env"

METHOD="scnode"
SCENARIO="A"
BASE_CONFIG="benchmark/configs/scnode_gse230659_observed_scgpt_v1.yaml"
FORMAL_CONFIG="benchmark/configs/scnode_gse230659_observed_scgpt_v1_A_hvg2000_formal.yaml"
HVG_ADATA="benchmark/inputs/gse230659_scgpt_hvg2000/GSE230659_scGPT_annotated_HVG2000_benchmark_input.h5ad"
OUTPUT_DIR="benchmark/results/scnode/scenario_A_scgpt_v1_hvg2000_formal"

N_HVG=2000
N_SIM_CELLS=2000
PRETRAIN_ITERS=200
EPOCHS=10
ITERS=100
BATCH_SIZE=32
LATENT_DIM=50

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export SCGPT_REPO="${SCGPT_REPO}"

echo "============================================================"
echo "Job started at: $(date)"
echo "Host          : $(hostname)"
echo "Project root  : ${TRAJ_PROJECT_ROOT}"
echo "Method        : ${METHOD}"
echo "Scenario      : ${SCENARIO} HVG${N_HVG} formal"
echo "HVG adata     : ${HVG_ADATA}"
echo "Formal config : ${FORMAL_CONFIG}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Job ID        : ${JOB_ID:-N/A}"
echo "============================================================"

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

python - <<'PY'
import anndata, scanpy, scipy, torch, torchdiffeq
print("imports OK")
print("anndata:", anndata.__version__)
print("scanpy:", scanpy.__version__)
print("scipy:", scipy.__version__)
print("torch:", torch.__version__)
PY

echo "------------------------------------------------------------"
echo "Preparing Scenario A formal config"
echo "------------------------------------------------------------"

python - <<PY
from pathlib import Path
import yaml
import anndata as ad

root = Path("${PROJECT_ROOT}")
hvg_adata = root / "${HVG_ADATA}"
base_config = root / "${BASE_CONFIG}"
formal_config = root / "${FORMAL_CONFIG}"

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

cfg["run_id"] = "scnode_gse230659_A_scgpt_v1_hvg2000_formal"
cfg["scenario"] = "A"
cfg["scenario_name"] = "Observed-time formal scNODE run (HVG2000, scGPT-v1)"
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
cfg.pop("scenario_params", None)

formal_config.parent.mkdir(parents=True, exist_ok=True)
with open(formal_config, "w", encoding="utf-8") as f:
    yaml.safe_dump(cfg, f, sort_keys=False)
print(f"wrote formal config: {formal_config}")
PY

CMD=(
  python -m benchmark.evaluation.eval_dispatch
  --method "${METHOD}"
  --scenario "${SCENARIO}"
  --adata "${HVG_ADATA}"
  --output-dir "${OUTPUT_DIR}"
  --method-config "${FORMAL_CONFIG}"
)

echo "------------------------------------------------------------"
echo "Running command:"
printf ' %q' "${CMD[@]}"
echo
echo "------------------------------------------------------------"
time "${CMD[@]}"

echo "------------------------------------------------------------"
echo "Completing forecast metrics and formal metadata"
echo "------------------------------------------------------------"

python - <<PY
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad

from benchmark.evaluation.eval_forecast import run_forecast_evaluation

root = Path("${PROJECT_ROOT}")
out = root / "${OUTPUT_DIR}"
hvg_adata = root / "${HVG_ADATA}"

adata = ad.read_h5ad(hvg_adata, backed="r")
run_forecast_evaluation(
    projected_expression_path=str(out / "projected_expression.npy"),
    adata=adata,
    output_dir=str(out),
    time_key="abs_day",
    max_cells_per_timepoint=1000,
)

meta_path = out / "run_metadata.json"
with open(meta_path, encoding="utf-8") as f:
    meta = json.load(f)
meta["result_class"] = "official"
meta["formal_benchmark"] = True
meta["feature_space"] = {
    "type": "HVG",
    "n_genes": int("${N_HVG}"),
    "source": "${HVG_ADATA}",
}
meta["training_setting"] = {
    "label": "hvg2000_formal",
    "pretrain_iters": int("${PRETRAIN_ITERS}"),
    "epochs": int("${EPOCHS}"),
    "iters": int("${ITERS}"),
    "batch_size": int("${BATCH_SIZE}"),
    "n_sim_cells": int("${N_SIM_CELLS}"),
}
with open(meta_path, "w", encoding="utf-8") as f:
    json.dump(meta, f, indent=2)

stm = pd.read_csv(out / "state_transition_matrix.csv", index_col=0)
vals = stm.values.astype(float)
print("stm_shape:", stm.shape)
print("stm_nonfinite:", int((~np.isfinite(vals)).sum()))
print("stm_zero_rows:", stm.index[np.isclose(vals.sum(axis=1), 0)].tolist())
if not np.isfinite(vals).all():
    raise SystemExit(1)
PY

echo "------------------------------------------------------------"
echo "Output check"
echo "------------------------------------------------------------"

python - <<PY
from pathlib import Path
import json

out = Path("${PROJECT_ROOT}") / "${OUTPUT_DIR}"
required = [
    "trained_scnode_model.pth",
    "projected_expression.npy",
    "forecast_metrics.json",
    "per_timepoint_forecast_metrics.csv",
    "projected_embedding.npy",
    "embedding_metrics.json",
    "per_timepoint_embedding_metrics.csv",
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

for name in ["run_metadata.json", "forecast_metrics.json", "embedding_metrics.json", "lineage_metrics.json"]:
    print(f"\\n== {name} ==")
    with open(out / name, encoding="utf-8") as f:
        data = json.load(f)
    print(json.dumps(data, indent=2)[:2400])
PY

echo "============================================================"
echo "Job finished at: $(date)"
echo "============================================================"
