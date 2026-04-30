#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_cellrank2_formal_C.$JOB_ID.log
#$ -l s_vmem=128G

set -euo pipefail

PROJECT_ROOT="/home/xzy0723/projects/trajectory"
SCGPT_REPO="/home/xzy0723/projects/scGPT"
CONDA_SH="/home/xzy0723/miniconda3/etc/profile.d/conda.sh"
CONDA_ENV="traj_env"

METHOD="cellrank2"
SCENARIO="C"
ADATA="benchmark/results/scgpt/full/adata_scgpt_annotated.h5ad"
METHOD_CONFIG="benchmark/configs/cellrank2_gse230659_observed_scgpt_v1_scenarioC.yaml"
OUTPUT_DIR="benchmark/results/cellrank2/scenario_C_scgpt_v1"

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
echo "Scenario      : ${SCENARIO} formal full-data"
echo "ADATA         : ${ADATA}"
echo "Method config : ${METHOD_CONFIG}"
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
import cellrank
import scvelo
import wot
print("imports OK")
print("anndata :", anndata.__version__)
print("scanpy  :", scanpy.__version__)
print("cellrank:", cellrank.__version__)
print("scvelo  :", scvelo.__version__)
PY

echo "------------------------------------------------------------"
echo "Preflight and pilot-output backup"
echo "------------------------------------------------------------"

python - <<PY
from pathlib import Path
import json
import shutil
import time
import yaml
import anndata as ad

root = Path("${PROJECT_ROOT}")
adata_path = root / "${ADATA}"
config_path = root / "${METHOD_CONFIG}"
out = root / "${OUTPUT_DIR}"

if not adata_path.exists():
    raise FileNotFoundError(f"Missing adata: {adata_path}")
if not config_path.exists():
    raise FileNotFoundError(f"Missing method config: {config_path}")

with open(config_path, encoding="utf-8-sig") as f:
    cfg = yaml.safe_load(f)
expected = cfg.get("scenario_params", {}).get("n_train_cells_expected")
print("expected_train_cells:", expected)

adata = ad.read_h5ad(adata_path, backed="r")
time_key = cfg.get("dataset", {}).get("time_key", "abs_day")
train_times = {float(x) for x in cfg.get("scenario_params", {}).get("train_times", [])}
obs_times = adata.obs[time_key].astype(float)
n_train = int(obs_times.isin(train_times).sum())
print("actual_train_cells:", n_train)
if expected is not None and n_train != int(expected):
    raise RuntimeError(f"Scenario-C train cell count mismatch: expected {expected}, got {n_train}")

meta_path = out / "run_metadata.json"
if meta_path.exists():
    with open(meta_path, encoding="utf-8") as f:
        meta = json.load(f)
    result_class = meta.get("result_class")
    sampling = meta.get("sampling")
    if result_class == "official" and not sampling:
        raise RuntimeError(f"Official output already exists at {out}; refusing to overwrite.")
    stamp = time.strftime("%Y%m%d_%H%M%S")
    backup = out.with_name(out.name + f"_pilot_backup_{stamp}")
    print(f"Backing up existing non-official output: {out} -> {backup}")
    shutil.move(str(out), str(backup))
elif out.exists():
    stamp = time.strftime("%Y%m%d_%H%M%S")
    backup = out.with_name(out.name + f"_backup_{stamp}")
    print(f"Backing up existing output without metadata: {out} -> {backup}")
    shutil.move(str(out), str(backup))
else:
    print("No existing output directory to back up.")
PY

CMD=(
  python -m benchmark.evaluation.eval_dispatch
  --method "${METHOD}"
  --scenario "${SCENARIO}"
  --adata "${ADATA}"
  --output-dir "${OUTPUT_DIR}"
  --method-config "${METHOD_CONFIG}"
)

echo "------------------------------------------------------------"
echo "Running command:"
printf ' %q' "${CMD[@]}"
echo
echo "------------------------------------------------------------"
time "${CMD[@]}"

echo "------------------------------------------------------------"
echo "Marking official metadata and validating outputs"
echo "------------------------------------------------------------"

python - <<PY
from pathlib import Path
import json
import numpy as np
import pandas as pd

out = Path("${PROJECT_ROOT}") / "${OUTPUT_DIR}"
required = [
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

meta_path = out / "run_metadata.json"
with open(meta_path, encoding="utf-8") as f:
    meta = json.load(f)
meta["result_class"] = "official"
meta["formal_benchmark"] = True
meta["dataset"] = "GSE230659"
meta["notes"] = (
    str(meta.get("notes", "")).rstrip()
    + " Formal Scenario C full-data run; no pilot subsampling."
).strip()
meta.pop("result_class_reason", None)
meta.pop("sampling", None)
with open(meta_path, "w", encoding="utf-8") as f:
    json.dump(meta, f, indent=2)

stm = pd.read_csv(out / "state_transition_matrix.csv", index_col=0)
vals = stm.values.astype(float)
print("stm_shape:", stm.shape)
print("stm_nonfinite:", int((~np.isfinite(vals)).sum()))
print("stm_zero_rows:", stm.index[np.isclose(vals.sum(axis=1), 0)].tolist())
if not np.isfinite(vals).all():
    raise SystemExit(1)

for name in ["run_metadata.json", "lineage_metrics.json"]:
    print(f"\\n== {name} ==")
    with open(out / name, encoding="utf-8") as f:
        data = json.load(f)
    print(json.dumps(data, indent=2)[:2400])
PY

echo "============================================================"
echo "Job finished at: $(date)"
echo "============================================================"
