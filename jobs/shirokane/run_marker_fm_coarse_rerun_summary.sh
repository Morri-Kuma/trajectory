#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_marker_fm_coarse_rerun_summary.$JOB_ID.log
#$ -l s_vmem=32G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
RUN_OFFICIAL_SUMMARY="${RUN_OFFICIAL_SUMMARY:-1}"
LEIDEN_N_NEIGHBORS="${LEIDEN_N_NEIGHBORS:-15}"
LEIDEN_RESOLUTION="${LEIDEN_RESOLUTION:-0.5}"

cd "${PROJECT_ROOT}"
mkdir -p logs

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"

echo "============================================================"
echo "Marker FM coarse rerun summary"
echo "Started at          : $(date)"
echo "Project root        : ${PROJECT_ROOT}"
echo "Official summary    : ${RUN_OFFICIAL_SUMMARY}"
echo "Leiden n_neighbors  : ${LEIDEN_N_NEIGHBORS}"
echo "Leiden resolution   : ${LEIDEN_RESOLUTION}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

python - <<'PY'
from pathlib import Path
import json

import pandas as pd
import yaml

configs = [
    "benchmark/configs/runtime/scnode_gse178325_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse178325_marker_fm_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse178325_marker_fm_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse178325_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse178325_marker_fm_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse178325_marker_fm_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse178325_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse178325_marker_fm_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse178325_marker_fm_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/scnode_gse230659_marker_fm_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse230659_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse230659_marker_fm_silver_B_hvg2000_formal.yaml",
    "benchmark/configs/runtime/prescient_gse230659_marker_fm_silver_C_hvg2000_formal.yaml",
    "benchmark/configs/runtime/cellrank2_gse230659_marker_fm_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/wot_gse230659_marker_fm_silver_A_hvg2000_formal.yaml",
]

rows = []
errors = []
for cfg_path in configs:
    cfg_path = Path(cfg_path)
    cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8-sig")) or {}
    out_dir = Path((cfg.get("output") or {}).get("base_dir"))
    metrics_path = out_dir / "lineage_metrics.json"
    stm_path = out_dir / "state_transition_matrix.csv"
    edges_path = out_dir / "lineage_graph_edges.csv"
    missing = [str(p) for p in (metrics_path, stm_path, edges_path) if not p.exists()]
    if missing:
        errors.append(f"{cfg_path}: missing {missing}")
        continue

    metrics = json.loads(metrics_path.read_text(encoding="utf-8-sig"))
    stm = pd.read_csv(stm_path, index_col=0)
    labels = [str(x) for x in stm.index] + [str(x) for x in stm.columns]
    stage = sorted({x for x in labels if x.startswith("stage_")})
    if stage:
        errors.append(f"{cfg_path}: STM has stage labels {stage[:20]}")
    if "stage_" in edges_path.read_text(encoding="utf-8", errors="replace"):
        errors.append(f"{cfg_path}: lineage_graph_edges.csv contains stage_ labels")
    if metrics.get("status") != "completed":
        errors.append(f"{cfg_path}: lineage status={metrics.get('status')!r}")
    expected_label_source = "reference_graph_nodes_sctimebench_zero_filled"
    if metrics.get("prediction_label_source") != expected_label_source:
        errors.append(
            f"{cfg_path}: prediction_label_source={metrics.get('prediction_label_source')!r}"
        )
    if metrics.get("cell_state_key") != "final_milestone_label_coarse":
        errors.append(f"{cfg_path}: cell_state_key={metrics.get('cell_state_key')!r}")

    graph = json.loads(Path((cfg.get("lineage") or {}).get("reference_graph_path")).read_text(encoding="utf-8"))
    ref_nodes = [str(n["id"]) for n in graph.get("nodes", [])]
    rows_in_stm = [str(x) for x in stm.index]
    cols_in_stm = [str(x) for x in stm.columns]
    expected_missing = set(x for x in ref_nodes if x not in set(rows_in_stm)) | set(
        x for x in ref_nodes if x not in set(cols_in_stm)
    )
    label_report = metrics.get("prediction_label_report") or {}
    reported_missing = set(label_report.get("missing_reference_nodes") or [])
    if expected_missing != reported_missing:
        errors.append(
            f"{cfg_path}: prediction_label_report missing_reference_nodes mismatch "
            f"expected={sorted(expected_missing)}, reported={sorted(reported_missing)}"
        )

    single = (metrics.get("graph_metrics") or {}).get("single_step") or {}
    multi = (metrics.get("graph_metrics") or {}).get("multi_step") or {}
    rows.append({
        "config": str(cfg_path),
        "method": cfg.get("method"),
        "dataset": (cfg.get("dataset") or {}).get("id"),
        "scenario": cfg.get("scenario"),
        "status": metrics.get("status"),
        "single_step_auc_roc": single.get("auc_roc"),
        "single_step_auc_prc": single.get("auc_prc"),
        "single_step_jaccard": single.get("jaccard_similarity"),
        "multi_step_auc_roc": multi.get("auc_roc"),
        "multi_step_auc_prc": multi.get("auc_prc"),
        "multi_step_jaccard": multi.get("jaccard_similarity"),
        "missing_reference_nodes": ";".join(sorted(reported_missing)),
        "n_missing_reference_nodes": len(reported_missing),
    })

summary_path = Path("benchmark/reports/official_silver/marker_fm_coarse_lineage_rerun_summary.csv")
summary_path.parent.mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(summary_path, index=False)
print(f"Wrote {summary_path} with {len(rows)} rows.")

if errors:
    print("Coarse rerun validation failed:")
    for error in errors:
        print(f"- {error}")
    raise SystemExit(1)

print("All marker FM coarse lineage rerun outputs passed scTimeBench zero-fill validation.")
PY

if [ "${RUN_OFFICIAL_SUMMARY}" = "1" ]; then
  python -m benchmark.evaluation.summarize_official_silver \
    --leiden-n-neighbors "${LEIDEN_N_NEIGHBORS}" \
    --leiden-resolution "${LEIDEN_RESOLUTION}"
fi

echo "============================================================"
echo "Marker FM coarse rerun summary finished at: $(date)"
echo "============================================================"
