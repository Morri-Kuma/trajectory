#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -N v21_inject_umap
#$ -o logs/run_v21_inject_umap.$JOB_ID.log
#$ -l s_vmem=64G
# ── Inject a reference UMAP (.obsm['X_umap']) into the two new-dataset benchmark
#    inputs so the official embedding-coherence evaluator has a reference embedding
#    (matching the primary datasets). Runs on a COMPUTE node: the neighbor search
#    spawns threads, which the login node forbids ("can't start new thread").
#    Does NOT re-annotate or re-freeze providers; existing labels/benchmark results
#    are untouched. Idempotent (skips if X_umap already present).
# ─────────────────────────────────────────────────────────────────────────────
set -euo pipefail
PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
cd "${PROJECT_ROOT}"; mkdir -p logs
export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"; export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1
source "${CONDA_SH}"; conda activate "${CONDA_ENV}"

for H in \
  benchmark/inputs/gse218855_marker_fm_transition_silver_hvg2000/GSE218855_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad \
  benchmark/inputs/gse298212_marker_fm_transition_silver_hvg2000/GSE298212_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad
do
  echo "=== inject UMAP -> ${H} at $(date) ==="
  python scripts/inject_reference_umap.py --input-h5ad "${H}"
done

echo "=== verifying obsm ==="
python - <<'PY'
import anndata as ad
for p in ["benchmark/inputs/gse218855_marker_fm_transition_silver_hvg2000/GSE218855_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad",
          "benchmark/inputs/gse298212_marker_fm_transition_silver_hvg2000/GSE298212_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad"]:
    a = ad.read_h5ad(p, backed="r")
    assert "X_umap" in a.obsm, f"X_umap MISSING in {p}"
    print(p.split("/")[-1], "obsm:", list(a.obsm.keys()), "X_umap:", a.obsm["X_umap"].shape)
print("ALL GOOD — submit run_v21_new_datasets_embedding_coherence_array.sh next.")
PY
echo "=== done at $(date) ==="
