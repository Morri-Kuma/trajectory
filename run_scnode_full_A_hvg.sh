#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_scnode_full_A_hvg.$JOB_ID.log
#$ -l s_vmem=128G

set -euo pipefail

# ===== user-configurable section =====
PROJECT_ROOT="/home/xzy0723/projects/trajectory"
SCGPT_REPO="/home/xzy0723/projects/scGPT"
CONDA_SH="/home/xzy0723/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="traj_env"

METHOD="scnode"
SCENARIO="A"
BASE_CONFIG="benchmark/configs/scnode_gse230659_observed_scgpt_v1.yaml"
HVG_CONFIG="benchmark/configs/scnode_gse230659_observed_scgpt_v1_full_hvg2000_reduced.yaml"
SOURCE_ADATA="benchmark/results/scgpt/full/adata_scgpt_annotated.h5ad"
HVG_ADATA="benchmark/inputs/gse230659_scgpt_hvg2000/GSE230659_scGPT_annotated_HVG2000_benchmark_input.h5ad"
OUTPUT_DIR="benchmark/results/scnode/scenario_A_scgpt_v1_full_hvg2000_reduced"

N_HVG=2000

# Conservative first full-cell run on HVG2000. If this completes, increase
# toward formal values: pretrain_iters=200, epochs=10, iters=100, n_sim_cells=2000.
N_SIM_CELLS=1000
PRETRAIN_ITERS=50
EPOCHS=3
ITERS=20
BATCH_SIZE=32
LATENT_DIM=50
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
echo "Scenario      : ${SCENARIO} full-cell HVG${N_HVG} reduced-training"
echo "Source adata  : ${SOURCE_ADATA}"
echo "HVG adata     : ${HVG_ADATA}"
echo "Base config   : ${BASE_CONFIG}"
echo "HVG config    : ${HVG_CONFIG}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Job ID        : ${JOB_ID:-N/A}"
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
echo "Building/reusing full-cell HVG${N_HVG} adata"
echo "------------------------------------------------------------"

python - <<PY
from pathlib import Path
import numpy as np
import yaml
import anndata as ad
import scanpy as sc
import scipy.sparse as sp

root = Path("${PROJECT_ROOT}")
source_adata = root / "${SOURCE_ADATA}"
hvg_adata = root / "${HVG_ADATA}"
base_config = root / "${BASE_CONFIG}"
hvg_config = root / "${HVG_CONFIG}"
n_hvg = int("${N_HVG}")

if not source_adata.exists():
    raise FileNotFoundError(f"Source adata not found: {source_adata}")

hvg_adata.parent.mkdir(parents=True, exist_ok=True)
hvg_config.parent.mkdir(parents=True, exist_ok=True)

if hvg_adata.exists():
    print(f"reuse existing HVG adata: {hvg_adata}")
    adata_check = ad.read_h5ad(hvg_adata, backed="r")
    print(f"HVG adata shape={adata_check.shape}")
else:
    print(f"loading source adata into memory for HVG selection: {source_adata}")
    adata = ad.read_h5ad(source_adata)
    print(f"source shape={adata.shape}")

    required_obs = ["abs_day", "scgpt_pseudostate_provisional"]
    missing = [c for c in required_obs if c not in adata.obs.columns]
    if missing:
        raise KeyError(f"Missing required obs columns: {missing}")

    if "highly_variable" in adata.var.columns and int(adata.var["highly_variable"].sum()) >= n_hvg:
        print("using existing adata.var['highly_variable']")
        hv = adata.var["highly_variable"].astype(bool).values
        if hv.sum() > n_hvg:
            if "highly_variable_rank" in adata.var.columns:
                ranks = adata.var["highly_variable_rank"].to_numpy()
                order = np.argsort(np.where(np.isnan(ranks), np.inf, ranks))
                keep = np.zeros(adata.n_vars, dtype=bool)
                keep[order[:n_hvg]] = True
                hv = keep
            elif "dispersions_norm" in adata.var.columns:
                score = adata.var["dispersions_norm"].to_numpy()
                order = np.argsort(np.nan_to_num(score, nan=-np.inf))[::-1]
                keep = np.zeros(adata.n_vars, dtype=bool)
                keep[order[:n_hvg]] = True
                hv = keep
            else:
                keep = np.zeros(adata.n_vars, dtype=bool)
                keep[np.flatnonzero(hv)[:n_hvg]] = True
                hv = keep
    else:
    print(f"computing top {n_hvg} benchmark HVGs with finite sparse mean/variance dispersion")
        X = adata.X
        if sp.issparse(X):
            X = X.tocsr(copy=True)
            if X.data.size:
                bad = ~np.isfinite(X.data)
                if bad.any():
                    print(f"WARNING: replacing {int(bad.sum())} non-finite sparse values with 0 before HVG scoring")
                    X.data[bad] = 0.0
            means = np.asarray(X.mean(axis=0)).ravel()
            mean_sq = np.asarray(X.power(2).mean(axis=0)).ravel()
        else:
            X = np.asarray(X)
            bad = ~np.isfinite(X)
            if bad.any():
                print(f"WARNING: replacing {int(bad.sum())} non-finite dense values with 0 before HVG scoring")
                X = np.where(bad, 0.0, X)
            means = X.mean(axis=0)
            mean_sq = np.square(X).mean(axis=0)

        variances = np.maximum(mean_sq - np.square(means), 0.0)
        dispersion = variances / np.maximum(means, 1e-12)
        valid = np.isfinite(dispersion) & np.isfinite(means) & (means > 0)
        if int(valid.sum()) < n_hvg:
            raise ValueError(f"Only {int(valid.sum())} valid genes available for HVG selection; requested {n_hvg}")

        scores = np.where(valid, dispersion, -np.inf)
        top_idx = np.argpartition(scores, -n_hvg)[-n_hvg:]
        top_idx = top_idx[np.argsort(scores[top_idx])[::-1]]
        hv = np.zeros(adata.n_vars, dtype=bool)
        hv[top_idx] = True
        adata.var["highly_variable"] = hv
        adata.var["benchmark_hvg_score"] = scores

    if int(hv.sum()) != n_hvg:
        print(f"WARNING: selected {int(hv.sum())} HVGs, expected {n_hvg}")

    adata_hvg = adata[:, hv].copy()
    adata_hvg.uns["benchmark_input_id"] = "GSE230659_scGPT_annotated_HVG2000"
    adata_hvg.uns["benchmark_input_label"] = (
        "GSE230659 scGPT-annotated HVG2000 benchmark input"
    )
    adata_hvg.uns["benchmark_hvg_source"] = str(source_adata)
    adata_hvg.uns["benchmark_hvg_n_top_genes"] = n_hvg
    adata_hvg.uns["benchmark_hvg_selection_method"] = (
        "Top genes by finite sparse mean/variance dispersion from "
        "scGPT-annotated full-gene AnnData."
    )

    # Keep expression in float32 to reduce downstream dense memory pressure.
    if sp.issparse(adata_hvg.X):
        adata_hvg.X = adata_hvg.X.astype(np.float32)
    else:
        adata_hvg.X = np.asarray(adata_hvg.X, dtype=np.float32)

    adata_hvg.write_h5ad(hvg_adata)
    print(f"wrote HVG adata: {hvg_adata} shape={adata_hvg.shape}")

adata_meta = ad.read_h5ad(hvg_adata, backed="r")
print(f"final HVG adata shape={adata_meta.shape}")
print("timepoints:", sorted(adata_meta.obs["abs_day"].astype(float).unique()))
print("states:", len(set(adata_meta.obs["scgpt_pseudostate_provisional"].astype(str))))

with open(base_config, encoding="utf-8-sig") as f:
    cfg = yaml.safe_load(f)

cfg["run_id"] = "scnode_gse230659_A_scgpt_v1_full_hvg2000_reduced"
cfg["scenario"] = "A"
cfg["dataset"]["id"] = "GSE230659"
cfg["dataset"]["h5ad_path"] = "${HVG_ADATA}"
cfg["dataset"]["time_key"] = "abs_day"
cfg["output"]["base_dir"] = "${OUTPUT_DIR}"
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

# Scenario A uses all observed time points.
cfg.pop("scenario_params", None)

with open(hvg_config, "w", encoding="utf-8") as f:
    yaml.safe_dump(cfg, f, sort_keys=False)
print(f"wrote HVG reduced config: {hvg_config}")
PY

CMD=(
  python -m benchmark.evaluation.eval_dispatch
  --method "${METHOD}"
  --scenario "${SCENARIO}"
  --adata "${HVG_ADATA}"
  --output-dir "${OUTPUT_DIR}"
  --method-config "${HVG_CONFIG}"
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
    "trained_scnode_model.pth",
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
