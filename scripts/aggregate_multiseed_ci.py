#!/usr/bin/env python
"""Aggregate multi-seed runs into per-method mean +/- 95% CI.

Reads the per-seed result directories produced by run_multiseed_benchmark_array.sh
(benchmark/results/<method>/<ds>_marker_fm_silver_<sc>_hvg2000_seed<NN>/) and, for
each (method, dataset, scenario), computes the mean and a 95% confidence interval
(Student-t, df = n-1) across seeds for the three headline metrics:
  lineage  single_step_auc_roc   (lineage_metrics.json)
  forecast wasserstein_distance  (forecast_metrics.json)
  embedding adjusted_rand_index  (embedding_metrics.json)
Writes benchmark/reports/official_silver/multiseed_ci.csv.
"""
from __future__ import annotations

import csv
import json
import math
import os
from collections import defaultdict
from pathlib import Path

ROOT = Path(os.environ.get("TRAJ_PROJECT_ROOT", os.getcwd()))
RES = ROOT / "benchmark" / "results"
OUT = ROOT / "benchmark" / "reports" / "official_silver" / "multiseed_ci.csv"

METHODS = ["scnode", "prescient", "mioflow"]
DATASETS = ["gse178325", "gse230659"]
SCEN = ["A", "B", "C"]
SEEDS = [101, 202, 303, 404, 505]

# Student-t 97.5th percentile by degrees of freedom (n-1), for small n.
T975 = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447,
        7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228}

METRICS = [
    ("lineage_single_step_auroc", "lineage_metrics.json", "single_step_auc_roc"),
    ("forecast_wasserstein", "forecast_metrics.json", "wasserstein_distance"),
    ("embedding_ari", "embedding_metrics.json", "adjusted_rand_index"),
]


def get(path: Path, key: str):
    try:
        d = json.loads(path.read_text())
    except Exception:
        return None
    v = d.get(key)
    # lineage nests single_step under graph_metrics in some schemas; try both
    if v is None and key == "single_step_auc_roc":
        v = d.get("auroc")
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def stats(vals):
    vals = [v for v in vals if v is not None]
    n = len(vals)
    if n == 0:
        return None
    mean = sum(vals) / n
    if n == 1:
        return {"n": 1, "mean": mean, "sd": 0.0, "ci95_half": 0.0,
                "ci_lo": mean, "ci_hi": mean}
    sd = math.sqrt(sum((x - mean) ** 2 for x in vals) / (n - 1))
    sem = sd / math.sqrt(n)
    t = T975.get(n - 1, 1.96)
    half = t * sem
    return {"n": n, "mean": mean, "sd": sd, "ci95_half": half,
            "ci_lo": mean - half, "ci_hi": mean + half}


def main() -> None:
    rows = []
    for method in METHODS:
        for ds in DATASETS:
            for sc in SCEN:
                for mlabel, fname, key in METRICS:
                    vals = []
                    for seed in SEEDS:
                        rundir = RES / method / f"{ds}_marker_fm_silver_{sc}_hvg2000_seed{seed}"
                        vals.append(get(rundir / fname, key))
                    s = stats(vals)
                    if s is None:
                        continue
                    rows.append({
                        "method": method, "dataset": ds.upper(), "scenario": sc,
                        "metric": mlabel, "n_seeds": s["n"],
                        "mean": round(s["mean"], 4), "sd": round(s["sd"], 4),
                        "ci95_half_width": round(s["ci95_half"], 4),
                        "ci95_lo": round(s["ci_lo"], 4), "ci95_hi": round(s["ci_hi"], 4),
                    })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    keys = ["method", "dataset", "scenario", "metric", "n_seeds",
            "mean", "sd", "ci95_half_width", "ci95_lo", "ci95_hi"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"[aggregate_multiseed_ci] wrote {len(rows)} rows -> {OUT}")
    if not rows:
        print("[aggregate_multiseed_ci] no seed runs found yet — run the array first.")


if __name__ == "__main__":
    main()
