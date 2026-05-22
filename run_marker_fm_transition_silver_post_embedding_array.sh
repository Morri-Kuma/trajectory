#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -t 1-18
#$ -o logs/run_marker_fm_transition_silver_post_embedding.$JOB_ID.$TASK_ID.log
#$ -l s_vmem=64G

set -euo pipefail

PROJECT_ROOT="${TRAJ_PROJECT_ROOT:-/home/xzy0723/projects/trajectory}"
CONDA_SH="${CONDA_SH:-/home/xzy0723/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-traj_env}"
TASK_ID="${SGE_TASK_ID:-${TASK_ID:-1}}"
REBUILD_PROJECTED_MILESTONE_LABELS="${REBUILD_PROJECTED_MILESTONE_LABELS:-1}"

mkdir -p "${PROJECT_ROOT}/logs"
cd "${PROJECT_ROOT}"

export TRAJ_PROJECT_ROOT="${PROJECT_ROOT}"
export PYTHONPATH="${PROJECT_ROOT}:${PYTHONPATH:-}"
export OMP_NUM_THREADS="${NSLOTS:-1}"
export OPENBLAS_NUM_THREADS="${NSLOTS:-1}"
export MKL_NUM_THREADS="${NSLOTS:-1}"
export NUMEXPR_NUM_THREADS="${NSLOTS:-1}"

if [ ! -f "${CONDA_SH}" ]; then
  echo "ERROR: conda init script not found: ${CONDA_SH}" >&2
  exit 2
fi
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"

mapfile -t CONFIGS < <(python - <<'PY'
import glob
from pathlib import Path
import yaml

for p in sorted(glob.glob("benchmark/configs/runtime/*marker_fm_silver*_formal.yaml")):
    cfg = yaml.safe_load(Path(p).read_text(encoding="utf-8-sig")) or {}
    if str(cfg.get("method", "")).lower() in {"scnode", "mioflow", "prescient"}:
        print(p)
PY
)
N_CONFIGS="${#CONFIGS[@]}"

if [ "${N_CONFIGS}" -eq 0 ]; then
  echo "ERROR: no projection-capable marker_fm_silver configs found." >&2
  exit 3
fi
if [ "${TASK_ID}" -lt 1 ] || [ "${TASK_ID}" -gt "${N_CONFIGS}" ]; then
  echo "Task ${TASK_ID} is outside config range 1-${N_CONFIGS}; nothing to do."
  exit 0
fi

CONFIG="${CONFIGS[$((TASK_ID - 1))]}"
read -r METHOD SCENARIO DATASET H5AD_PATH OUTPUT_DIR RUN_ID PROVIDER_ID < <(python - "${CONFIG}" <<'PY'
import sys
from pathlib import Path
import yaml

cfg = yaml.safe_load(Path(sys.argv[1]).read_text(encoding="utf-8-sig")) or {}
method = str(cfg.get("method", "")).lower()
scenario = str(cfg.get("scenario", ""))
dataset = str((cfg.get("dataset") or {}).get("id", ""))
h5ad = str((cfg.get("dataset") or {}).get("h5ad_path", ""))
out = str((cfg.get("output") or {}).get("base_dir", ""))
run_id = str(cfg.get("run_id", Path(sys.argv[1]).stem))
provider_id = str((cfg.get("ground_truth") or {}).get("provider_id", f"{dataset.lower()}_marker_fm_transition_silver_v1"))
print(method, scenario, dataset, h5ad, out, run_id, provider_id)
PY
)

if [ "${DATASET}" = "GSE178325" ]; then
  EXPECTED_PROVIDER="gse178325_marker_fm_transition_silver_v1"
elif [ "${DATASET}" = "GSE230659" ]; then
  EXPECTED_PROVIDER="gse230659_marker_fm_transition_silver_v1"
else
  echo "ERROR: unsupported dataset for marker_fm post embedding: ${DATASET}" >&2
  exit 4
fi

EMBEDDING_OUT="${OUTPUT_DIR}/embedding_milestone_eval"

echo "============================================================"
echo "Marker-FM post-embedding task started at: $(date)"
echo "Host          : $(hostname)"
echo "Task ID       : ${TASK_ID}/${N_CONFIGS}"
echo "Config        : ${CONFIG}"
echo "Run ID        : ${RUN_ID}"
echo "Dataset       : ${DATASET}"
echo "Method        : ${METHOD}"
echo "Scenario      : ${SCENARIO}"
echo "Provider      : ${PROVIDER_ID}"
echo "Input h5ad    : ${H5AD_PATH}"
echo "Output dir    : ${OUTPUT_DIR}"
echo "Embedding out : ${EMBEDDING_OUT}"
echo "============================================================"

if [ "${PROVIDER_ID}" != "${EXPECTED_PROVIDER}" ]; then
  echo "ERROR: unexpected provider for ${DATASET}: ${PROVIDER_ID}" >&2
  exit 5
fi
if [ ! -f "${H5AD_PATH}" ]; then
  echo "ERROR: reference h5ad missing: ${H5AD_PATH}" >&2
  exit 6
fi
if [ ! -d "${OUTPUT_DIR}" ]; then
  echo "ERROR: method output dir missing: ${OUTPUT_DIR}" >&2
  exit 6
fi
if [ ! -f "${OUTPUT_DIR}/projected_embedding.npy" ]; then
  echo "ERROR: projected_embedding.npy missing in ${OUTPUT_DIR}" >&2
  exit 6
fi
if [ ! -f "${OUTPUT_DIR}/projected_cluster_labels.csv" ]; then
  echo "ERROR: projected_cluster_labels.csv missing in ${OUTPUT_DIR}" >&2
  exit 6
fi

if [ "${REBUILD_PROJECTED_MILESTONE_LABELS}" = "1" ] || [ ! -f "${OUTPUT_DIR}/projected_milestone_labels.csv" ]; then
  echo "------------------------------------------------------------"
  echo "Rebuild projected_milestone_labels.csv from projected_cluster_labels.csv"
  echo "------------------------------------------------------------"
  python - "${OUTPUT_DIR}" "${DATASET}" "${METHOD}" "${SCENARIO}" "${PROVIDER_ID}" <<'PY'
import csv
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

out_dir = Path(sys.argv[1])
dataset_id = sys.argv[2]
method = sys.argv[3]
scenario = sys.argv[4]
provider_id = sys.argv[5]

cluster_path = out_dir / "projected_cluster_labels.csv"
rows = list(csv.DictReader(cluster_path.open(newline="", encoding="utf-8")))
if not rows:
    raise SystemExit(f"empty projected_cluster_labels.csv: {cluster_path}")

label_key = None
for candidate in ("projected_cluster_label", "cluster_label", "label"):
    if candidate in rows[0]:
        label_key = candidate
        break
if label_key is None:
    raise SystemExit(f"no projected label column found in {cluster_path}: {list(rows[0])}")

time_key = None
for candidate in ("timepoint", "projected_timepoint", "target_timepoint"):
    if candidate in rows[0]:
        time_key = candidate
        break
if time_key is None:
    raise SystemExit(f"no timepoint column found in {cluster_path}: {list(rows[0])}")

out_rows = []
for i, row in enumerate(rows):
    label = str(row.get(label_key, "") or "unknown_or_ood")
    tp = str(row.get(time_key, ""))
    projected_cell_id = f"{method}_{dataset_id}_{scenario}_projected_{i:08d}"
    out_rows.append({
        "projected_cell_id": projected_cell_id,
        "cell_idx": row.get("cell_idx", str(i)),
        "timepoint": tp,
        "projected_timepoint": tp,
        "final_milestone_label_coarse": label,
        "final_milestone_label_expanded": label,
        "final_milestone_confidence": "1.0",
        "final_milestone_source": "projected_cluster_label_nearest_centroid_bridge",
        "label_mode": "official_silver",
        "provider_id": provider_id,
        "projected_cluster_label": label,
    })

csv_path = out_dir / "projected_milestone_labels.csv"
with csv_path.open("w", newline="", encoding="utf-8") as f:
    fieldnames = list(out_rows[0].keys())
    writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(out_rows)

counts = Counter(r["final_milestone_label_expanded"] for r in out_rows)
meta = {
    "dataset_id": dataset_id,
    "method": method,
    "scenario": scenario,
    "provider_id": provider_id,
    "label_mode": "official_silver",
    "state_key": "final_milestone_label_expanded",
    "n_projected_cells": len(out_rows),
    "label_counts": dict(counts),
    "annotation_policy": "projected_cluster_label_nearest_centroid_bridge",
    "annotation_policy_note": (
        "This bridge converts the method runner's projected_cluster_label "
        "column into projected_milestone_labels.csv so run-output Leiden "
        "clusters can be compared against projected milestone labels. "
        "Confidence is a numeric placeholder and is not calibrated."
    ),
    "source_projected_cluster_labels": str(cluster_path),
    "output_projected_milestone_labels": str(csv_path),
    "smoke_test": False,
    "formal_benchmark": True,
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
}
(out_dir / "projected_milestone_annotation_metadata.json").write_text(
    json.dumps(meta, indent=2),
    encoding="utf-8",
)
print(json.dumps({
    "projected_milestone_labels": str(csv_path),
    "n_projected_cells": len(out_rows),
    "label_counts": dict(counts),
}, indent=2))
PY
fi

echo "------------------------------------------------------------"
echo "Validate projected milestone labels"
echo "------------------------------------------------------------"
python -m benchmark.annotation.validate_projected_milestone_labels "${OUTPUT_DIR}" --verbose

mkdir -p "${EMBEDDING_OUT}"

echo "------------------------------------------------------------"
echo "Run official_silver milestone embedding coherence"
echo "------------------------------------------------------------"
python -m benchmark.evaluation.eval_embedding_milestone \
  --input-h5ad "${H5AD_PATH}" \
  --run-output-dir "${OUTPUT_DIR}" \
  --output-dir "${EMBEDDING_OUT}" \
  --dataset-id "${DATASET}" \
  --mode run-output \
  --label-mode official_silver \
  --exclude-label ambiguous \
  --exclude-label unknown_or_ood

echo "------------------------------------------------------------"
echo "Validate embedding metric JSON"
echo "------------------------------------------------------------"
python - "${EMBEDDING_OUT}" "${EXPECTED_PROVIDER}" <<'PY'
import json
import sys
from pathlib import Path

out = Path(sys.argv[1])
provider = sys.argv[2]
metrics_path = out / "embedding_metrics_official_silver.json"
if not metrics_path.exists():
    raise SystemExit(f"missing {metrics_path}")
m = json.loads(metrics_path.read_text(encoding="utf-8"))
print(json.dumps({
    "status": m.get("status"),
    "label_mode": m.get("label_mode"),
    "provider_id": m.get("provider_id"),
    "state_key": m.get("state_key"),
    "n_cells_evaluated": m.get("n_cells_evaluated"),
    "adjusted_rand_index": m.get("adjusted_rand_index"),
    "cluster_source": m.get("cluster_source"),
}, indent=2))
if m.get("status") != "completed":
    raise SystemExit(f"embedding metrics not completed: {m.get('status')}")
if m.get("label_mode") != "official_silver":
    raise SystemExit(f"unexpected label_mode: {m.get('label_mode')}")
if m.get("provider_id") != provider:
    raise SystemExit(f"unexpected provider_id: {m.get('provider_id')}")
if m.get("state_key") != "final_milestone_label_expanded":
    raise SystemExit(f"unexpected state_key: {m.get('state_key')}")
if not isinstance(m.get("n_cells_evaluated"), int) or m["n_cells_evaluated"] <= 0:
    raise SystemExit(f"invalid n_cells_evaluated: {m.get('n_cells_evaluated')}")
PY

echo "============================================================"
echo "Marker-FM post-embedding task finished at: $(date)"
echo "Embedding out: ${EMBEDDING_OUT}"
echo "============================================================"
