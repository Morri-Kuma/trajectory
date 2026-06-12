#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_silver_lineage_graphsim.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G
#$ -t 1-40

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"
RUN_TAG="${RUN_TAG:-lineage_graphsim_$(date +%Y%m%d_%H%M%S)}"
RUN_BASELINE="${RUN_BASELINE:-1}"
BACKUP_OLD="${BACKUP_OLD:-1}"

CONFIGS=(
  "benchmark/configs/runtime/scnode_gse175634_cardiac_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse175634_cardiac_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse175634_cardiac_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse175634_cardiac_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse175634_cardiac_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse175634_cardiac_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse175634_cardiac_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse175634_cardiac_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse175634_cardiac_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/cellrank2_gse175634_cardiac_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/wot_gse175634_cardiac_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse178325_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse178325_marker_fm_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse178325_marker_fm_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse178325_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse178325_marker_fm_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse178325_marker_fm_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse178325_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse178325_marker_fm_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse178325_marker_fm_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse230659_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse230659_marker_fm_silver_B_hvg2000_formal.yaml"
  "benchmark/configs/runtime/prescient_gse230659_marker_fm_silver_C_hvg2000_formal.yaml"
  "benchmark/configs/runtime/cellrank2_gse230659_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/runtime/wot_gse230659_marker_fm_silver_A_hvg2000_formal.yaml"
  "benchmark/configs/scnode_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml"
  "benchmark/configs/scnode_gse242424_oskm_ground_truth_B_hvg2000_formal.yaml"
  "benchmark/configs/scnode_gse242424_oskm_ground_truth_C_hvg2000_formal.yaml"
  "benchmark/configs/mioflow_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml"
  "benchmark/configs/mioflow_gse242424_oskm_ground_truth_B_hvg2000_formal.yaml"
  "benchmark/configs/mioflow_gse242424_oskm_ground_truth_C_hvg2000_formal.yaml"
  "benchmark/configs/prescient_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml"
  "benchmark/configs/prescient_gse242424_oskm_ground_truth_B_hvg2000_formal.yaml"
  "benchmark/configs/prescient_gse242424_oskm_ground_truth_C_hvg2000_formal.yaml"
)

if [ "${TASK_ID}" -lt 1 ] || [ "${TASK_ID}" -gt "${#CONFIGS[@]}" ]; then
  echo "ERROR: unsupported task id ${TASK_ID}; expected 1-${#CONFIGS[@]}" >&2
  exit 2
fi

CONFIG="${CONFIGS[$((TASK_ID - 1))]}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export CONFIG
export RUN_TAG
export RUN_BASELINE
export BACKUP_OLD
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Job started at : $(date)"
echo "Task ID        : ${TASK_ID}/${#CONFIGS[@]}"
echo "Config         : ${CONFIG}"
echo "Project root   : ${PROJECT_ROOT}"
echo "Run tag        : ${RUN_TAG}"
echo "Run baseline   : ${RUN_BASELINE}"
echo "Backup old     : ${BACKUP_OLD}"
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
import shutil
from pathlib import Path

import yaml

from benchmark.evaluation.eval_lineage import run_lineage_evaluation


ROOT = Path(os.environ["TRAJ_PROJECT_ROOT"]).resolve()
CONFIG_PATH = Path(os.environ["CONFIG"])
if not CONFIG_PATH.is_absolute():
    CONFIG_PATH = ROOT / CONFIG_PATH

RUN_BASELINE = os.environ.get("RUN_BASELINE", "1") == "1"
BACKUP_OLD = os.environ.get("BACKUP_OLD", "1") == "1"
RUN_TAG = os.environ.get("RUN_TAG", "lineage_graphsim")


def read_yaml(path: Path) -> dict:
    with open(path, encoding="utf-8-sig") as f:
        return yaml.safe_load(f) or {}


def resolve_path(value):
    if value in (None, ""):
        return None
    path = Path(str(value))
    return path if path.is_absolute() else ROOT / path


def first_present(*values, default=None):
    for value in values:
        if value not in (None, ""):
            return value
    return default


def bool_value(value, default=False) -> bool:
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def load_baseline_adata(h5ad_path: Path, time_key: str, train_times):
    import anndata as ad
    import pandas as pd

    print(f"[lineage_graphsim] Loading AnnData for baseline: {h5ad_path}")
    adata = ad.read_h5ad(str(h5ad_path))
    if not train_times:
        return adata
    if time_key not in adata.obs.columns:
        raise KeyError(f"time_key {time_key!r} not found in {h5ad_path}")

    series = adata.obs[time_key]
    mask = series.isin(train_times)
    if int(mask.sum()) == 0:
        mask = series.astype(str).isin({str(t) for t in train_times})
    if int(mask.sum()) == 0:
        numeric = pd.to_numeric(series, errors="coerce")
        train_numeric = [float(t) for t in train_times]
        mask = numeric.isin(train_numeric)
    if int(mask.sum()) == 0:
        raise ValueError(
            f"scenario train_times={train_times!r} matched no cells under {time_key!r}"
        )

    print(
        "[lineage_graphsim] Baseline AnnData filtered to "
        f"{int(mask.sum())}/{adata.n_obs} cells using scenario train_times."
    )
    return adata[mask].copy()


cfg = read_yaml(CONFIG_PATH)
registry_path = ROOT / "benchmark" / "ground_truth" / "registry.yaml"
registry = read_yaml(registry_path).get("providers", {})

dataset = cfg.get("dataset") or {}
lineage = cfg.get("lineage") or {}
cfg_gt = cfg.get("ground_truth") or {}
state_system = cfg.get("state_system") or {}
scenario_params = cfg.get("scenario_params") or {}

provider_id = first_present(
    cfg.get("provider_id"),
    cfg_gt.get("provider_id"),
    lineage.get("state_system_version"),
    state_system.get("version"),
)
registry_gt = dict(registry.get(provider_id, {})) if provider_id else {}

output_dir = resolve_path((cfg.get("output") or {}).get("base_dir"))
if output_dir is None:
    raise SystemExit(f"output.base_dir missing in {CONFIG_PATH}")

state_transition_matrix = output_dir / "state_transition_matrix.csv"
lineage_graph_edges = output_dir / "lineage_graph_edges.csv"
reference_graph_path = resolve_path(
    first_present(
        lineage.get("reference_graph_path"),
        cfg_gt.get("reference_graph_path"),
        registry_gt.get("reference_graph_path"),
    )
)
h5ad_path = resolve_path(
    first_present(
        dataset.get("h5ad_path"),
        state_system.get("source_h5ad"),
        registry_gt.get("source_h5ad"),
    )
)
cell_state_key = first_present(
    lineage.get("cell_state_key"),
    cfg_gt.get("state_key"),
    cfg.get("cell_state_key"),
    registry_gt.get("state_key"),
    default="final_milestone_label_coarse",
)
time_key = first_present(
    dataset.get("time_key"),
    cfg.get("time_key"),
    default="abs_day",
)
edge_confidence_mode = first_present(
    lineage.get("edge_confidence_mode"),
    cfg_gt.get("confidence_mode"),
    registry_gt.get("confidence_mode"),
    default="all",
)
exclude_uncertain_states = bool_value(
    first_present(
        lineage.get("exclude_uncertain_states"),
        cfg_gt.get("exclude_uncertain_states"),
        registry_gt.get("exclude_uncertain_states"),
    ),
    default=False,
)

ground_truth = {**registry_gt, **cfg_gt}
if provider_id:
    ground_truth["provider_id"] = provider_id
ground_truth.setdefault("state_key", cell_state_key)
ground_truth.setdefault("reference_graph_path", str(reference_graph_path.relative_to(ROOT)))
ground_truth.setdefault("confidence_mode", edge_confidence_mode)
ground_truth.setdefault("exclude_uncertain_states", exclude_uncertain_states)

print(json.dumps({
    "config": str(CONFIG_PATH.relative_to(ROOT)),
    "method": cfg.get("method"),
    "scenario": cfg.get("scenario"),
    "dataset": dataset.get("id"),
    "output_dir": str(output_dir.relative_to(ROOT)),
    "provider_id": provider_id,
    "reference_graph_path": str(reference_graph_path.relative_to(ROOT)),
    "cell_state_key": cell_state_key,
    "time_key": time_key,
    "edge_confidence_mode": edge_confidence_mode,
    "exclude_uncertain_states": exclude_uncertain_states,
    "run_baseline": RUN_BASELINE,
}, indent=2))

for required in [CONFIG_PATH, output_dir, state_transition_matrix, lineage_graph_edges, reference_graph_path]:
    if required is None or not Path(required).exists():
        raise SystemExit(f"required path missing: {required}")
if RUN_BASELINE and (h5ad_path is None or not h5ad_path.exists()):
    raise SystemExit(f"baseline requested but h5ad path is missing: {h5ad_path}")

metrics_path = output_dir / "lineage_metrics.json"
if BACKUP_OLD and metrics_path.exists():
    backup_path = output_dir / f"lineage_metrics.pre_graphsim_{RUN_TAG}.json"
    if not backup_path.exists():
        shutil.copy2(metrics_path, backup_path)
        print(f"[lineage_graphsim] Backed up previous metrics to {backup_path}")

adata = None
if RUN_BASELINE:
    adata = load_baseline_adata(
        h5ad_path=h5ad_path,
        time_key=str(time_key),
        train_times=scenario_params.get("train_times"),
    )

metrics = run_lineage_evaluation(
    state_transition_matrix_path=str(state_transition_matrix),
    lineage_graph_edges_path=str(lineage_graph_edges),
    output_dir=str(output_dir),
    reference_graph_path=str(reference_graph_path),
    adata=adata,
    edge_confidence_mode=str(edge_confidence_mode),
    exclude_uncertain_states=exclude_uncertain_states,
    cell_state_key=str(cell_state_key),
    time_key=str(time_key),
    ground_truth=ground_truth,
)

if metrics.get("metric_protocol") != "sctimebench_graph_sim":
    raise SystemExit(f"unexpected metric_protocol: {metrics.get('metric_protocol')!r}")
if metrics.get("status") != "completed":
    raise SystemExit(f"lineage evaluation did not complete cleanly: {metrics.get('status')!r}")
for criterion in ("single_step", "multi_step"):
    values = (metrics.get("graph_metrics") or {}).get(criterion) or {}
    for key in ("auc_roc", "auc_prc", "jaccard_similarity"):
        if values.get(key) is None:
            raise SystemExit(f"missing graph_metrics.{criterion}.{key}")

print(json.dumps({
    "status": metrics.get("status"),
    "metric_protocol": metrics.get("metric_protocol"),
    "single_step": (metrics.get("graph_metrics") or {}).get("single_step"),
    "multi_step": (metrics.get("graph_metrics") or {}).get("multi_step"),
    "baseline_note": (metrics.get("baseline") or {}).get("note"),
}, indent=2))
PY

echo "============================================================"
echo "Lineage graph-sim task finished at: $(date)"
echo "Config: ${CONFIG}"
echo "============================================================"
