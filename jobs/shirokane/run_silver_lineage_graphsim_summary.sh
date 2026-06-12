#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -o logs/run_silver_lineage_graphsim_summary.$JOB_ID.log
#$ -l s_vmem=16G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
RUN_TESTS="${RUN_TESTS:-1}"
RUN_OFFICIAL_SUMMARY="${RUN_OFFICIAL_SUMMARY:-1}"
REQUIRE_ALL_GRAPH_SIM="${REQUIRE_ALL_GRAPH_SIM:-1}"
LEIDEN_N_NEIGHBORS="${LEIDEN_N_NEIGHBORS:-15}"
LEIDEN_RESOLUTION="${LEIDEN_RESOLUTION:-0.5}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"

echo "============================================================"
echo "Summary job started at : $(date)"
echo "Project root           : ${PROJECT_ROOT}"
echo "Run tests              : ${RUN_TESTS}"
echo "Run official summary   : ${RUN_OFFICIAL_SUMMARY}"
echo "Require all graph-sim  : ${REQUIRE_ALL_GRAPH_SIM}"
echo "Leiden n_neighbors     : ${LEIDEN_N_NEIGHBORS}"
echo "Leiden resolution      : ${LEIDEN_RESOLUTION}"
echo "============================================================"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

if [ "${RUN_TESTS}" = "1" ]; then
  echo "------------------------------------------------------------"
  echo "Running scTimeBench graph-sim parity tests"
  echo "------------------------------------------------------------"
  if python - <<'PY'
import importlib.util
raise SystemExit(0 if importlib.util.find_spec("pytest") else 1)
PY
  then
    python -m pytest benchmark/evaluation/tests/test_lineage_graphsim_sctimebench.py -q
  else
    echo "WARNING: pytest is not installed in this environment; skipping parity tests."
    echo "         Install pytest or set RUN_TESTS=0 to silence this warning."
  fi
fi

echo "------------------------------------------------------------"
echo "Writing lineage graph-sim rerun inventory"
echo "------------------------------------------------------------"
python - <<'PY'
import csv
import json
import os
from pathlib import Path

ROOT = Path.cwd()
require_all = os.environ.get("REQUIRE_ALL_GRAPH_SIM", "1") == "1"
report_dir = ROOT / "benchmark" / "reports" / "silver_lineage_graphsim"
report_dir.mkdir(parents=True, exist_ok=True)
rows = []
for path in sorted((ROOT / "benchmark" / "results").glob("*/*_formal/lineage_metrics.json")):
    rel = path.relative_to(ROOT)
    if "_rerun_backup" in rel.parts:
        continue
    with open(path, encoding="utf-8") as f:
        metrics = json.load(f)
    run_dir = path.parent
    gm = metrics.get("graph_metrics") or {}
    single = gm.get("single_step") or {}
    multi = gm.get("multi_step") or {}
    rows.append({
        "result_dir": str(run_dir.relative_to(ROOT)),
        "status": metrics.get("status"),
        "metric_protocol": metrics.get("metric_protocol"),
        "provider_id": metrics.get("provider_id"),
        "cell_state_key": metrics.get("cell_state_key"),
        "edge_confidence_mode": metrics.get("edge_confidence_mode"),
        "n_reference_edges": metrics.get("n_reference_edges"),
        "single_step_auc_roc": single.get("auc_roc"),
        "single_step_auc_prc": single.get("auc_prc"),
        "single_step_jaccard": single.get("jaccard_similarity"),
        "multi_step_auc_roc": multi.get("auc_roc"),
        "multi_step_auc_prc": multi.get("auc_prc"),
        "multi_step_jaccard": multi.get("jaccard_similarity"),
        "baseline_note": (metrics.get("baseline") or {}).get("note"),
    })

out = report_dir / "lineage_graphsim_rerun_inventory.csv"
with open(out, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["result_dir"])
    writer.writeheader()
    writer.writerows(rows)

bad = [
    row for row in rows
    if row["status"] != "completed" or row["metric_protocol"] != "sctimebench_graph_sim"
]
print(f"Wrote {len(rows)} lineage rows to {out.relative_to(ROOT)}")
if bad:
    print("Invalid lineage rows:")
    for row in bad:
        print(f"- {row['result_dir']}: status={row['status']!r}, protocol={row['metric_protocol']!r}")
    if require_all:
        raise SystemExit(1)
PY

if [ "${RUN_OFFICIAL_SUMMARY}" = "1" ]; then
  echo "------------------------------------------------------------"
  echo "Regenerating official-silver marker-FM ranking tables"
  echo "------------------------------------------------------------"
  python -m benchmark.evaluation.summarize_official_silver \
    --leiden-n-neighbors "${LEIDEN_N_NEIGHBORS}" \
    --leiden-resolution "${LEIDEN_RESOLUTION}"
fi

echo "============================================================"
echo "Summary job finished at: $(date)"
echo "Inventory: benchmark/reports/silver_lineage_graphsim/lineage_graphsim_rerun_inventory.csv"
echo "============================================================"
