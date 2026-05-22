#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_marker_fm_expanded_official_repair_scnode_lineage.$JOB_ID.log
#$ -l s_vmem=64G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
FORCE_REPAIR="${FORCE_REPAIR:-0}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-1}"
export VECLIB_MAXIMUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Repair scNODE expanded-official lineage metrics"
echo "Started at   : $(date)"
echo "Host         : $(hostname)"
echo "Project root : ${PROJECT_ROOT}"
echo "Force repair : ${FORCE_REPAIR}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

python - <<'PY'
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import yaml

from benchmark.evaluation.eval_lineage import run_lineage_evaluation

force = os.environ.get("FORCE_REPAIR", "0") == "1"
configs = []
for p in sorted(Path("benchmark/configs/runtime").glob("scnode_*marker_fm_silver*_formal.yaml")):
    cfg = yaml.safe_load(p.read_text(encoding="utf-8-sig")) or {}
    configs.append((p, cfg))

if not configs:
    raise SystemExit("No scNODE marker_fm_silver configs found.")

completed = []
skipped = []
failed = []

for cfg_path, cfg in configs:
    run_id = cfg.get("run_id", cfg_path.stem)
    dataset_cfg = cfg.get("dataset") or {}
    lineage_cfg = cfg.get("lineage") or {}
    gt = cfg.get("ground_truth") or {}
    sp = cfg.get("scenario_params") or {}
    out = Path((cfg.get("output") or {}).get("base_dir", ""))
    stm = out / "state_transition_matrix.csv"
    edges = out / "lineage_graph_edges.csv"
    metrics = out / "lineage_metrics.json"
    meta_path = out / "run_metadata.json"
    h5ad_path = Path(dataset_cfg.get("h5ad_path", ""))
    time_key = dataset_cfg.get("time_key") or cfg.get("time_key") or "abs_day"
    cell_state_key = (
        gt.get("state_key")
        or lineage_cfg.get("cell_state_key")
        or cfg.get("cell_state_key")
        or "final_milestone_label_expanded"
    )
    ref = Path(lineage_cfg.get("reference_graph_path", ""))
    edge_mode = gt.get("confidence_mode") or lineage_cfg.get("edge_confidence_mode") or "all"
    exclude_uncertain = bool(
        gt.get("exclude_uncertain_states", lineage_cfg.get("exclude_uncertain_states", False))
    )
    train_times = sp.get("train_times")

    print("------------------------------------------------------------")
    print(f"Repairing {run_id}")
    print(f"  output          : {out}")
    print(f"  h5ad            : {h5ad_path}")
    print(f"  cell_state_key  : {cell_state_key}")
    print(f"  reference_graph : {ref}")

    if metrics.exists() and not force:
        print("  lineage_metrics.json already exists; skipping")
        skipped.append(run_id)
        continue
    missing = [str(p) for p in [stm, edges, h5ad_path, ref] if not p.exists()]
    if missing:
        print(f"  ERROR: missing required input(s): {missing}")
        failed.append((run_id, missing))
        continue

    adata = ad.read_h5ad(h5ad_path, backed="r")
    try:
        if train_times:
            train_f = {float(x) for x in train_times}
            t = adata.obs[time_key].astype(float).to_numpy()
            mask = np.array([float(x) in train_f for x in t], dtype=bool)
            print(f"  scenario filter : {int(mask.sum())}/{adata.n_obs} cells")
            eval_adata = adata[mask].to_memory()
            try:
                adata.file.close()
            except Exception:
                pass
        else:
            eval_adata = adata

        out_metrics = run_lineage_evaluation(
            state_transition_matrix_path=str(stm),
            lineage_graph_edges_path=str(edges),
            output_dir=str(out),
            reference_graph_path=str(ref),
            edge_confidence_mode=edge_mode,
            exclude_uncertain_states=exclude_uncertain,
            cell_state_key=cell_state_key,
            time_key=time_key,
            adata=eval_adata,
            ground_truth=gt,
        )
        print(
            "  wrote lineage_metrics.json:",
            out_metrics.get("status"),
            "auroc=",
            out_metrics.get("auroc"),
            "n_reference_edges=",
            out_metrics.get("n_reference_edges"),
        )

        if meta_path.exists():
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        else:
            meta = {}
        original_status = meta.get("status")
        meta["status"] = "completed"
        meta["original_status_before_lineage_repair"] = original_status
        meta["lineage_metrics_repaired_posthoc"] = True
        meta["lineage_metrics_repair_note"] = (
            "Original scNODE run wrote STM/projection outputs but failed while importing "
            "benchmark.evaluation.eval_lineage. This posthoc repair computed official "
            "expanded lineage metrics from the existing scNODE outputs."
        )
        meta["cell_state_key"] = cell_state_key
        meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
        completed.append(run_id)
    finally:
        try:
            if hasattr(adata, "file"):
                adata.file.close()
        except Exception:
            pass

print("============================================================")
print(json.dumps({
    "completed": completed,
    "skipped": skipped,
    "failed": failed,
}, indent=2))
if failed:
    raise SystemExit(1)
PY

echo "============================================================"
echo "scNODE lineage repair finished at: $(date)"
echo "============================================================"
