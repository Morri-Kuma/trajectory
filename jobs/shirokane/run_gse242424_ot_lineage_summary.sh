#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_gse242424_ot_lineage_summary.$JOB_ID.log
#$ -l s_vmem=16G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"

cd "${PROJECT_ROOT}"
mkdir -p logs benchmark/reports/gse242424_oskm_ground_truth

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"

echo "============================================================"
echo "GSE242424 OT lineage summary started at: $(date)"
echo "Project root: ${PROJECT_ROOT}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

python - <<'PY'
import json
from pathlib import Path

import pandas as pd
import yaml


ROOT = Path.cwd()
CONFIGS = [
    "benchmark/configs/wot_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml",
    "benchmark/configs/cellrank2_gse242424_oskm_ground_truth_A_hvg2000_formal.yaml",
]


def read_yaml(path: Path) -> dict:
    with open(path, encoding="utf-8-sig") as f:
        return yaml.safe_load(f) or {}


def read_json(path: Path) -> dict:
    with open(path, encoding="utf-8-sig") as f:
        return json.load(f) or {}


rows = []
errors = []
for rel_config in CONFIGS:
    config_path = ROOT / rel_config
    cfg = read_yaml(config_path)
    output_dir = ROOT / (cfg.get("output") or {}).get("base_dir")
    metrics_path = output_dir / "lineage_metrics.json"
    if not metrics_path.exists():
        errors.append(f"{rel_config}: missing {metrics_path.relative_to(ROOT)}")
        continue

    metrics = read_json(metrics_path)
    single = (metrics.get("graph_metrics") or {}).get("single_step") or {}
    multi = (metrics.get("graph_metrics") or {}).get("multi_step") or {}
    rows.append(
        {
            "dataset_id": (cfg.get("dataset") or {}).get("id"),
            "scenario": cfg.get("scenario"),
            "method": cfg.get("method"),
            "run_id": cfg.get("run_id"),
            "output_dir": str(output_dir.relative_to(ROOT)),
            "status": metrics.get("status"),
            "metric_protocol": metrics.get("metric_protocol"),
            "single_step_auc_roc": single.get("auc_roc"),
            "single_step_auc_prc": single.get("auc_prc"),
            "single_step_jaccard": single.get("jaccard_similarity"),
            "multi_step_auc_roc": multi.get("auc_roc"),
            "multi_step_auc_prc": multi.get("auc_prc"),
            "multi_step_jaccard": multi.get("jaccard_similarity"),
            "missing_reference_nodes": ";".join(
                map(
                    str,
                    (metrics.get("prediction_label_report") or {}).get(
                        "missing_reference_nodes",
                        [],
                    ),
                )
            ),
        }
    )
    if metrics.get("status") != "completed":
        errors.append(f"{rel_config}: status={metrics.get('status')!r}")
    if metrics.get("metric_protocol") != "sctimebench_graph_sim":
        errors.append(f"{rel_config}: metric_protocol={metrics.get('metric_protocol')!r}")

summary = pd.DataFrame(rows)
out_path = (
    ROOT
    / "benchmark"
    / "reports"
    / "gse242424_oskm_ground_truth"
    / "gse242424_ot_lineage_summary.csv"
)
summary.to_csv(out_path, index=False)
print(f"Wrote {out_path.relative_to(ROOT)} with {len(summary)} rows")

if errors:
    print("Validation errors:")
    for error in errors:
        print(f"- {error}")
    raise SystemExit(1)
PY

echo "============================================================"
echo "GSE242424 OT lineage summary finished at: $(date)"
echo "============================================================"
