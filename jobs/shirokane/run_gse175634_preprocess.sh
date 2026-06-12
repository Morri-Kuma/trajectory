#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_gse175634_preprocess.$JOB_ID.log
#$ -l s_vmem=64G
# ── GSE175634 preprocessing job ──────────────────────────────────────────────
#
# DATASET ID NOTE:
#   Canonical GEO accession: GSE175634
#   Source folder:           data/gse175634/
#   The name "GSE174534" that appeared in the original request is a typo
#   (digits transposed). This script uses GSE175634 throughout.
#
# What this script does (TWO sequential steps):
#   Step 1 — Build the HVG2000 benchmark input h5ad from raw MTX files.
#             Skips the silver-standard annotation step (author labels used
#             directly via the 'type' column in GSE175634_cell_metadata.tsv).
#   Step 2 — Build the frozen ground-truth provider
#             (gse175634_cardiac_silver_v1) from the h5ad produced in Step 1.
#
# Memory: up to ~32 GB (raw MTX is ~4 GB; in-memory anndata can be ~20 GB).
# Time  : ~15-30 min depending on IO speed on Shirokane.
#
# Adjust the three variables below to match your Shirokane environment:
#   TRAJ_PROJECT_ROOT — full path to the cloned repository on Shirokane
#   CONDA_SH          — path to conda.sh for your miniconda/anaconda install
#   CONDA_ENV         — name of the conda environment with scanpy/anndata/etc.
# ─────────────────────────────────────────────────────────────────────────────

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "GSE175634 preprocessing"
echo "Started at   : $(date)"
echo "Host         : $(hostname)"
echo "Project root : ${PROJECT_ROOT}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

# ── verify raw data are present ───────────────────────────────────────────────
echo ""
echo "------------------------------------------------------------"
echo "Checking raw data files …"
echo "------------------------------------------------------------"
python - <<'PY'
import os, sys
from pathlib import Path

root = Path(os.environ.get("TRAJ_PROJECT_ROOT", "."))
data_dir = root / "data" / "gse175634"
required = [
    "GSE175634_cell_counts.mtx",
    "gene_indices_counts.tsv",
    "GSE175634_cell_indices.tsv",
    "GSE175634_cell_metadata.tsv",
]
missing = [f for f in required if not (data_dir / f).exists()]
if missing:
    print(f"ERROR: missing raw data files in {data_dir}:", file=sys.stderr)
    for m in missing:
        print(f"  {m}", file=sys.stderr)
    sys.exit(1)
print(f"[OK] all raw data files present in {data_dir}")
PY

# ── Step 1: Build benchmark input h5ad ───────────────────────────────────────
echo ""
echo "------------------------------------------------------------"
echo "Step 1: Building benchmark input h5ad"
echo "  Silver-standard annotation step: SKIPPED"
echo "  (author 'type' labels used directly)"
echo "------------------------------------------------------------"
python scripts/build_gse175634_cardiac_author_input.py

H5AD_PATH="${PROJECT_ROOT}/benchmark/inputs/gse175634_cardiac_author_hvg2000/GSE175634_cardiac_author_HVG2000_benchmark_input.h5ad"
if [ ! -f "${H5AD_PATH}" ]; then
  echo "ERROR: expected h5ad not produced: ${H5AD_PATH}" >&2
  exit 3
fi
echo "[OK] Benchmark h5ad written: ${H5AD_PATH}"

# ── Step 2: Build ground-truth provider ──────────────────────────────────────
echo ""
echo "------------------------------------------------------------"
echo "Step 2: Building ground-truth provider (gse175634_cardiac_silver_v1)"
echo "------------------------------------------------------------"
python scripts/build_gse175634_cardiac_author_provider.py

PROVIDER_DIR="${PROJECT_ROOT}/benchmark/ground_truth/providers/gse175634_cardiac_silver_v1"
if [ ! -f "${PROVIDER_DIR}/reference_graph.json" ]; then
  echo "ERROR: provider reference_graph.json not found: ${PROVIDER_DIR}" >&2
  exit 3
fi
echo "[OK] Provider written: ${PROVIDER_DIR}"

# ── Sanity check ─────────────────────────────────────────────────────────────
echo ""
echo "------------------------------------------------------------"
echo "Running input validation …"
echo "------------------------------------------------------------"
python scripts/validate_gse175634_inputs.py || {
  echo "WARNING: validation reported issues (see above). Review before submitting benchmark." >&2
}

echo "============================================================"
echo "GSE175634 preprocessing DONE at: $(date)"
echo "Next step:"
echo "  bash jobs/shirokane/submit_gse175634_benchmark.sh"
echo "============================================================"
