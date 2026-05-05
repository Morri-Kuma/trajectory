#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_build_gse178325_input.$JOB_ID.log
#$ -l s_vmem=128G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"

SCRIPT="scripts/build_gse178325_benchmark_input.py"
OUTPUT_H5AD="benchmark/inputs/gse178325_human_hvg2000/GSE178325_human_HVG2000_benchmark_input.h5ad"
EXPECTED_BATCH="0618"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Job started at : $(date)"
echo "Host           : $(hostname)"
echo "Project root   : ${TRAJ_PROJECT_ROOT}"
echo "Builder script : ${SCRIPT}"
echo "Output h5ad    : ${OUTPUT_H5AD}"
echo "Batch policy   : ${EXPECTED_BATCH} only"
echo "Memory request : s_vmem=128G"
echo "Job ID         : ${JOB_ID:-N/A}"
echo "============================================================"

if [ ! -f "${SCRIPT}" ]; then
  echo "ERROR: builder script not found: ${SCRIPT}" >&2
  exit 2
fi

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi

source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

echo "------------------------------------------------------------"
echo "Environment check"
echo "------------------------------------------------------------"
echo "Python: $(which python)"
python -V
python - <<'PY'
import anndata
import matplotlib
import numpy
import pandas
import scipy
import seaborn
print("imports OK")
print("anndata:", anndata.__version__)
print("numpy:", numpy.__version__)
print("pandas:", pandas.__version__)
print("scipy:", scipy.__version__)
PY

echo "------------------------------------------------------------"
echo "Input 10x manifest check: 0618 batch only"
echo "------------------------------------------------------------"
python - <<'PY'
from pathlib import Path

root = Path.cwd()
data_dir = root / "data" / "gse178325_human" / "rna_seq_10x"
samples = [
    "GSM5683317_StageI_Day0-0618",
    "GSM5683318_StageI_Day0.5-0618",
    "GSM5683319_StageI_Day1-0618",
    "GSM5683320_StageI_Day2-0618",
    "GSM5534144_StageI_Day4-0618",
    "GSM5534145_StageII_Day4-0618",
    "GSM5534146_StageII_Day8-0618",
    "GSM5534147_StageII_Day12-0618",
    "GSM5534148_StageIII_Day4-0618",
    "GSM5534149_StageIII_Day8-0618",
    "GSM5534150_StageIV_Day1-0618",
    "GSM5534151_StageIV_Day2-0618",
    "GSM5534152_StageIV_Day4-0618",
    "GSM5534153_StageIV_Day10-0618",
    "GSM5534154_hCiPSCs-0618",
]
missing = []
for sample_name in samples:
    gsm = sample_name.split("_", 1)[0]
    if "0618" not in sample_name:
        raise SystemExit(f"non-0618 sample in manifest: {sample_name}")
    sample_dir = data_dir / gsm
    for filename in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]:
        path = sample_dir / filename
        if not path.exists():
            missing.append(str(path))
print("samples:", len(samples))
print("missing:", missing)
if missing:
    raise SystemExit(1)
PY

echo "------------------------------------------------------------"
echo "Running GSE178325 benchmark input builder"
echo "------------------------------------------------------------"
time python "${SCRIPT}"

echo "------------------------------------------------------------"
echo "Output validation"
echo "------------------------------------------------------------"
python - <<PY
from pathlib import Path
import anndata as ad
import scipy.sparse as sp

out = Path("${OUTPUT_H5AD}")
if not out.exists():
    raise FileNotFoundError(f"Expected output not found: {out}")

adata = ad.read_h5ad(out, backed="r")
required_obs = [
    "cell_id",
    "dataset_id",
    "batch_id",
    "sample_id",
    "stage",
    "day_within_stage",
    "abs_day",
    "time_label",
    "scTimeBench_timepoint",
    "scTimeBench_cell_type",
]
missing_obs = [c for c in required_obs if c not in adata.obs.columns]
print("shape:", adata.shape)
print("missing_obs:", missing_obs)
print("dataset_id:", adata.uns.get("dataset_id"))
print("time_axis:", adata.uns.get("time_axis"))
print("batch_id values:", sorted(set(adata.obs["batch_id"].astype(str))))
print("timepoints:", sorted(adata.obs["abs_day"].astype(float).unique()))
print("X_pca shape:", adata.obsm["X_pca"].shape if "X_pca" in adata.obsm else None)
print("X sparse:", sp.issparse(adata.X))

if missing_obs:
    raise SystemExit(1)
if adata.uns.get("dataset_id") != "GSE178325":
    raise SystemExit("dataset_id mismatch")
if sorted(set(adata.obs["batch_id"].astype(str))) != ["${EXPECTED_BATCH}"]:
    raise SystemExit("batch_id validation failed")
if adata.n_vars != 2000:
    raise SystemExit(f"expected 2000 HVGs, got {adata.n_vars}")
if "X_pca" not in adata.obsm or adata.obsm["X_pca"].shape[1] != 50:
    raise SystemExit("X_pca validation failed")
PY

echo "============================================================"
echo "Job finished at: $(date)"
echo "Output: ${OUTPUT_H5AD}"
echo "============================================================"
