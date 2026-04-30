"""
Summarize core benchmark outputs across Forecast Accuracy, Embedding Coherence,
and Lineage Fidelity.

The script is intentionally tolerant of partially implemented dimensions:
WOT/CellRank2 have lineage-only outputs, while early scNODE runs may have
forecast/lineage metrics but deferred embedding metrics.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple


VALID_RESULT_CLASSES = {
    "official",
    "pilot",
    "pilot_backup",
    "reduced_validation",
    "smoke",
    "hpc_validation",
    "diagnostic",
}


def _project_root() -> Path:
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "benchmark").exists() and (
            candidate / "benchmark" / "results"
        ).exists():
            return candidate
    return here.parent.parent


def _read_json(path: Path) -> Dict:
    if not path.exists():
        return {}
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f) or {}
    except Exception:
        return {}


def _discover_runs(root: Path) -> List[Tuple[str, str, Path]]:
    runs: List[Tuple[str, str, Path]] = []
    results_root = root / "benchmark" / "results"
    if not results_root.exists():
        return runs
    for method_dir in sorted(p for p in results_root.iterdir() if p.is_dir()):
        method = method_dir.name
        for run_dir in sorted(p for p in method_dir.iterdir() if p.is_dir()):
            if not run_dir.name.startswith("scenario_"):
                continue
            has_any_metrics = any(
                (run_dir / name).exists()
                for name in (
                    "forecast_metrics.json",
                    "embedding_metrics.json",
                    "lineage_metrics.json",
                    "run_metadata.json",
                )
            )
            if not has_any_metrics:
                continue
            scenario = run_dir.name[len("scenario_"):]
            runs.append((method, scenario, run_dir))
    return runs


def _infer_result_class(meta: Dict) -> str:
    declared = meta.get("result_class")
    if declared in VALID_RESULT_CLASSES:
        return declared
    sampling = meta.get("sampling") or {}
    if sampling.get("subsample_per_timepoint") not in (None, "", False):
        return "pilot"
    return "official"


def _ground_truth(metrics: Dict, meta: Dict) -> Dict:
    gt = metrics.get("ground_truth") or meta.get("ground_truth") or {}
    return gt if isinstance(gt, dict) else {}


def _row(method: str, scenario: str, run_dir: Path) -> Dict:
    forecast = _read_json(run_dir / "forecast_metrics.json")
    embedding = _read_json(run_dir / "embedding_metrics.json")
    lineage = _read_json(run_dir / "lineage_metrics.json")
    meta = _read_json(run_dir / "run_metadata.json")
    gt = _ground_truth(lineage, meta)
    baseline = lineage.get("baseline") or {}
    return {
        "method": method,
        "scenario": scenario,
        "result_class": _infer_result_class(meta),
        "ground_truth_provider": gt.get("provider_id"),
        "ground_truth_status": gt.get("status"),
        "run_status": meta.get("status"),
        "forecast_status": forecast.get("status"),
        "forecast_wasserstein": forecast.get("wasserstein_distance"),
        "forecast_gaussian_mmd": forecast.get("gaussian_mmd"),
        "forecast_energy_distance_mmd": forecast.get("energy_distance_mmd"),
        "forecast_hausdorff_loss": forecast.get("hausdorff_loss"),
        "forecast_n_eval_timepoints": forecast.get("n_eval_timepoints"),
        "forecast_scenario_type": forecast.get("scenario_type"),
        "embedding_status": embedding.get("status"),
        "embedding_ari": embedding.get("adjusted_rand_index"),
        "embedding_entropy": embedding.get("avg_normalized_classifier_entropy"),
        "embedding_n_eval_timepoints": embedding.get("n_eval_timepoints"),
        "lineage_status": lineage.get("status"),
        "lineage_auroc": lineage.get("auroc"),
        "lineage_auprc": lineage.get("auprc"),
        "lineage_jaccard": lineage.get("jaccard_similarity"),
        "lineage_jaccard_topk": lineage.get("jaccard_similarity_topk"),
        "lineage_single_step_recovery": lineage.get("single_step_recovery"),
        "lineage_multi_step_recovery": lineage.get("multi_step_recovery"),
        "baseline_auroc": baseline.get("auroc"),
        "baseline_auprc": baseline.get("auprc"),
        "edge_confidence_mode": lineage.get("edge_confidence_mode"),
        "n_reference_edges": lineage.get("n_reference_edges"),
        "train_times": meta.get("train_times"),
        "heldout_times": meta.get("heldout_times"),
        "n_sim_cells": meta.get("n_sim_cells"),
        "runtime_seconds": meta.get("runtime_seconds"),
        "source": str(run_dir),
    }


def _write_csv(rows: List[Dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "method",
        "scenario",
        "result_class",
        "ground_truth_provider",
        "ground_truth_status",
        "run_status",
        "forecast_status",
        "forecast_wasserstein",
        "forecast_gaussian_mmd",
        "forecast_energy_distance_mmd",
        "forecast_hausdorff_loss",
        "forecast_n_eval_timepoints",
        "forecast_scenario_type",
        "embedding_status",
        "embedding_ari",
        "embedding_entropy",
        "embedding_n_eval_timepoints",
        "lineage_status",
        "lineage_auroc",
        "lineage_auprc",
        "lineage_jaccard",
        "lineage_jaccard_topk",
        "lineage_single_step_recovery",
        "lineage_multi_step_recovery",
        "baseline_auroc",
        "baseline_auprc",
        "edge_confidence_mode",
        "n_reference_edges",
        "train_times",
        "heldout_times",
        "n_sim_cells",
        "runtime_seconds",
        "source",
    ]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k) for k in fieldnames})


def _fmt(value) -> str:
    if value is None:
        return "-"
    if isinstance(value, float):
        if value != value:
            return "-"
        return f"{value:.4f}"
    return str(value)


def _print_table(rows: List[Dict]) -> None:
    cols = [
        "method",
        "scenario",
        "ground_truth_provider",
        "forecast_wasserstein",
        "forecast_gaussian_mmd",
        "forecast_energy_distance_mmd",
        "forecast_hausdorff_loss",
        "embedding_ari",
        "embedding_entropy",
        "lineage_auroc",
        "lineage_auprc",
        "baseline_auroc",
    ]
    if not rows:
        print("[summarize_core] No rows.")
        return
    widths = {c: max(len(c), max(len(_fmt(r.get(c))) for r in rows)) for c in cols}
    print("  ".join(c.ljust(widths[c]) for c in cols))
    print("  ".join("-" * widths[c] for c in cols))
    for row in sorted(rows, key=lambda r: (r["method"], r["scenario"])):
        print("  ".join(_fmt(row.get(c)).ljust(widths[c]) for c in cols))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        default="benchmark/reports/core_summary.csv",
        help="Path for the CSV summary output.",
    )
    parser.add_argument(
        "--official-only",
        action="store_true",
        help="Only include runs classified as official.",
    )
    args = parser.parse_args()

    root = _project_root()
    rows = [_row(method, scenario, run_dir) for method, scenario, run_dir in _discover_runs(root)]
    if args.official_only:
        rows = [row for row in rows if row["result_class"] == "official"]

    out_path = Path(args.output)
    if not out_path.is_absolute():
        out_path = root / out_path
    _write_csv(rows, out_path)
    _print_table(rows)
    print(f"\n[summarize_core] Wrote {len(rows)} rows -> {out_path}")


if __name__ == "__main__":
    main()
