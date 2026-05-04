"""
build_gse230659_dataset_pkl.py
Build and serialize the GSE230659 dataset object.

Usage:
    conda activate traj_env
    cd C:\\Users\\37620\\trajectory
    python benchmark/datasets/build_gse230659_dataset_pkl.py

Output:
    benchmark/datasets/gse230659_observed.dataset.pkl

What this script does:
    1. Auto-detects the project root (same logic as 01_load_and_qc.py).
    2. Instantiates GSE230659Dataset pointing at the canonical h5ad.
    3. Smoke-tests load_data() to confirm the split is valid.
    4. Pickles the lightweight dataset object (not the AnnData itself).
    5. Verifies the pkl round-trip by unpickling and calling load_data() again.
"""

from __future__ import annotations

import pickle
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Project-root detection (same logic as 01_load_and_qc.py)
# ---------------------------------------------------------------------------

def _find_project_root() -> Path:
    import os
    # 1. Explicit override: TRAJ_PROJECT_ROOT env var (Shirokane HPC / CI).
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
        raise FileNotFoundError(
            f"TRAJ_PROJECT_ROOT={env!r} does not exist.  "
            "Correct the environment variable and retry."
        )
    # 2. Walk upward from this script until a directory containing 'data/' is found.
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "data").exists():
            return candidate
    # 3. Last-resort fallback.
    return here.parent


PROJECT_ROOT = _find_project_root()
print(f"[build_pkl] Project root: {PROJECT_ROOT}")

# Make sure benchmark/ is importable regardless of how the script is called.
benchmark_root = PROJECT_ROOT / "benchmark"
if str(benchmark_root) not in sys.path:
    sys.path.insert(0, str(benchmark_root))
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from datasets.gse230659_dataset import GSE230659Dataset   # noqa: E402

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

H5AD_PATH  = PROJECT_ROOT / "data" / "processed" / "adata_benchmark.h5ad"
OUTPUT_PKL = PROJECT_ROOT / "benchmark" / "datasets" / "gse230659_observed.dataset.pkl"

# ---------------------------------------------------------------------------
# Step 1: Instantiate
# ---------------------------------------------------------------------------

print(f"\n[build_pkl] Instantiating GSE230659Dataset ...")
dataset = GSE230659Dataset(
    h5ad_path=H5AD_PATH,
    scenario="A",
    project_root=PROJECT_ROOT,
)
print(f"  {dataset}")

# ---------------------------------------------------------------------------
# Step 2: Smoke-test load_data()
# ---------------------------------------------------------------------------

print(f"\n[build_pkl] Smoke-testing load_data() ...")
train_adata, test_adata = dataset.load_data()

assert train_adata.n_obs == 75194, \
    f"Expected 75194 train cells, got {train_adata.n_obs}"
assert test_adata.n_obs > 0, \
    f"test_adata is empty — check terminal_timepoint={dataset.terminal_timepoint}"

# Confirm required obs columns are present in train
REQUIRED_OBS = [
    "cell_id", "dataset_id", "sample_id", "stage",
    "day_within_stage", "abs_day", "time_label",
    "scTimeBench_timepoint", "scTimeBench_cell_type",
]
missing = [c for c in REQUIRED_OBS if c not in train_adata.obs.columns]
if missing:
    raise AssertionError(f"train_adata missing required obs columns: {missing}")

assert "X_pca" in train_adata.obsm, "X_pca missing from train_adata.obsm"

print(f"  train: {train_adata.n_obs} cells × {train_adata.n_vars} genes")
print(f"  test : {test_adata.n_obs}  cells × {test_adata.n_vars}  genes")
print(f"  Smoke test PASSED ✓")

# Free memory before pickling (we pickle only the lightweight class, not the data)
del train_adata, test_adata

# ---------------------------------------------------------------------------
# Step 3: Pickle
# ---------------------------------------------------------------------------

print(f"\n[build_pkl] Writing pkl → {OUTPUT_PKL}")
with open(OUTPUT_PKL, "wb") as f:
    pickle.dump(dataset, f, protocol=pickle.HIGHEST_PROTOCOL)

pkl_size_kb = OUTPUT_PKL.stat().st_size / 1024
print(f"  Written: {pkl_size_kb:.1f} KB")

# ---------------------------------------------------------------------------
# Step 4: Round-trip verification
# ---------------------------------------------------------------------------

print(f"\n[build_pkl] Verifying round-trip ...")
with open(OUTPUT_PKL, "rb") as f:
    dataset_rt = pickle.load(f)

assert dataset_rt.dataset_id  == "GSE230659"
assert dataset_rt.scenario    == "A"
assert dataset_rt.h5ad_path.exists(), \
    f"h5ad path after unpickling does not exist: {dataset_rt.h5ad_path}"

train_rt, test_rt = dataset_rt.load_data()
assert train_rt.n_obs == 75194
assert test_rt.n_obs  > 0
del train_rt, test_rt

print(f"  Round-trip PASSED ✓")
print(f"\n[build_pkl] dataset.pkl ready: {OUTPUT_PKL}")
