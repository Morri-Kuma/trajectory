"""
make_representation_reports.py
==============================

Assemble the supplementary representation-dynamics report (Work-Plan Step 7)
from per-run evaluator outputs.

Scans a results root for:
    <results>/<method>/<dataset>_<rep>_<scenario>/
        representation_forecast_metrics.json
        representation_temporal_signal.json
        run_metadata.json (optional)
    plus representation/scfm metadata sidecars when present.

Writes the six convention files into the report directory:
    representation_method_summary.md
    representation_forecast_rankings.csv
    representation_temporal_signal.csv
    representation_lineage_metrics.csv
    representation_embedding_coherence.csv
    representation_metadata_manifest.csv

Usage:
    python -m benchmark.reports.representation_dynamics.make_representation_reports \
        --results-root benchmark/results/representation_dynamics \
        --report-dir   benchmark/reports/representation_dynamics
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

RUN_DIR_RE = re.compile(r"(?P<dataset>[^_]+)_(?P<rep>rep_[a-z0-9_]+?)_(?P<scenario>[A-F])$")

FORECAST_METRICS = (
    "rep_wasserstein_distance",
    "rep_gaussian_mmd",
    "rep_energy_distance_mmd",
    "rep_hausdorff_loss",
)
TEMPORAL_METRICS = (
    "temporal_signal_ratio",
    "adjacent_time_centroid_distance",
    "non_adjacent_time_centroid_distance",
    "time_prediction_macro_f1",
    "time_prediction_r2",
)
LINEAGE_METRICS = (
    "forward_mass_fraction",
    "backward_mass_fraction",
    "forbidden_edge_mass",
)
COHERENCE_METRICS = (
    "state_centroid_separation",
    "state_silhouette_score",
)


def _load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f) or {}
    except Exception:
        return {}


def _scan_runs(results_root: Path):
    runs = []
    if not results_root.exists():
        return runs
    for method_dir in sorted(p for p in results_root.iterdir() if p.is_dir()):
        method = method_dir.name
        for run_dir in sorted(p for p in method_dir.iterdir() if p.is_dir()):
            m = RUN_DIR_RE.search(run_dir.name)
            rep = m.group("rep") if m else run_dir.name
            scenario = m.group("scenario") if m else "?"
            dataset = m.group("dataset") if m else "?"
            runs.append({
                "method": method,
                "dataset": dataset,
                "representation": rep,
                "scenario": scenario,
                "dir": run_dir,
                "forecast": _load_json(run_dir / "representation_forecast_metrics.json"),
                "temporal": _load_json(run_dir / "representation_temporal_signal.json"),
                "run_metadata": _load_json(run_dir / "run_metadata.json"),
            })
    return runs


def _write_csv(path: Path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def build_reports(results_root: Path, report_dir: Path) -> dict:
    report_dir.mkdir(parents=True, exist_ok=True)
    runs = _scan_runs(results_root)

    base_cols = ["method", "dataset", "representation", "scenario"]

    # forecast rankings (ranked per scenario by rep_wasserstein_distance asc)
    forecast_rows = []
    for r in runs:
        forecast_rows.append(
            [r[c] for c in base_cols]
            + [r["forecast"].get(k) for k in FORECAST_METRICS]
            + [r["forecast"].get("metric_backend")]
        )
    # rank within scenario by primary metric (wasserstein), Nones last
    def _rank_key(row):
        scenario = row[3]
        val = row[4]
        return (scenario, float("inf") if val is None else val)
    forecast_rows.sort(key=_rank_key)
    # add a rank column per scenario
    ranked = []
    counters: dict = {}
    for row in forecast_rows:
        sc = row[3]
        counters[sc] = counters.get(sc, 0) + 1
        ranked.append(row + [counters[sc]])
    _write_csv(
        report_dir / "representation_forecast_rankings.csv",
        base_cols + list(FORECAST_METRICS) + ["metric_backend", "rank_in_scenario"],
        ranked,
    )

    _write_csv(
        report_dir / "representation_temporal_signal.csv",
        base_cols + list(TEMPORAL_METRICS),
        [[r[c] for c in base_cols] + [r["temporal"].get(k) for k in TEMPORAL_METRICS]
         for r in runs],
    )
    _write_csv(
        report_dir / "representation_lineage_metrics.csv",
        base_cols + list(LINEAGE_METRICS),
        [[r[c] for c in base_cols] + [r["temporal"].get(k) for k in LINEAGE_METRICS]
         for r in runs],
    )
    _write_csv(
        report_dir / "representation_embedding_coherence.csv",
        base_cols + list(COHERENCE_METRICS),
        [[r[c] for c in base_cols] + [r["temporal"].get(k) for k in COHERENCE_METRICS]
         for r in runs],
    )

    # metadata manifest: pull provenance from forecast/temporal embedded metadata
    manifest_cols = base_cols + [
        "reducer_fit_scope", "metric_backend", "model_name", "model_checkpoint",
        "paper_reference", "code_reference", "upstream_commit_or_release",
        "embedding_dim", "preprocessing", "pooling_or_cell_embedding_method",
    ]
    manifest_rows = []
    for r in runs:
        rep_meta = (r["forecast"].get("representation_metadata")
                    or r["temporal"].get("representation_metadata") or {})
        # rep_meta may be a build-inputs block and/or carry an scfm sub-block
        scfm = rep_meta.get("scfm", {}) if isinstance(rep_meta, dict) else {}
        manifest_rows.append([
            r["method"], r["dataset"], r["representation"], r["scenario"],
            rep_meta.get("reducer_fit_scope"),
            r["forecast"].get("metric_backend"),
            scfm.get("model_name"), scfm.get("model_checkpoint"),
            scfm.get("paper_reference"), scfm.get("code_reference"),
            scfm.get("upstream_commit_or_release"),
            scfm.get("embedding_dim"), scfm.get("preprocessing"),
            scfm.get("pooling_or_cell_embedding_method"),
        ])
    _write_csv(report_dir / "representation_metadata_manifest.csv",
               manifest_cols, manifest_rows)

    # human-readable summary
    summary = report_dir / "representation_method_summary.md"
    with open(summary, "w", encoding="utf-8") as f:
        f.write("# Representation-Dynamics Summary (supplementary)\n\n")
        f.write(f"Runs discovered: **{len(runs)}** under `{results_root}`.\n\n")
        f.write("> Supplementary representation-space results. NOT the official "
                "gene-expression silver benchmark. scFM embeddings are frozen and "
                "were not used to define the official silver labels.\n\n")
        f.write("## Forecast rankings (lower is better)\n\n")
        f.write("See `representation_forecast_rankings.csv`. Primary ranking "
                "metric: `rep_wasserstein_distance` (per scenario).\n\n")
        if not runs:
            f.write("_No runs found yet. Run the representation configs and "
                    "evaluators first._\n")
        else:
            f.write("| method | rep | scenario | rep_wasserstein | "
                    "temporal_signal_ratio | state_silhouette |\n")
            f.write("| --- | --- | --- | --- | --- | --- |\n")
            for r in runs:
                f.write(
                    f"| {r['method']} | {r['representation']} | {r['scenario']} | "
                    f"{r['forecast'].get('rep_wasserstein_distance')} | "
                    f"{r['temporal'].get('temporal_signal_ratio')} | "
                    f"{r['temporal'].get('state_silhouette_score')} |\n"
                )
        f.write("\n## Interpretation guide\n\n")
        f.write("- scFM better forecast **and** comparable temporal signal -> "
                "useful state-aware geometry.\n")
        f.write("- scFM better state separation but worse Scenario B forecast / "
                "lower temporal signal -> identity preserved, temporal signal "
                "compressed.\n")
        f.write("- HVG-PCA best across A/B/C -> local expression variation "
                "preserves trajectory signal better than generic pretrained reps.\n")

    return {"n_runs": len(runs), "report_dir": str(report_dir)}


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Assemble representation-dynamics report.")
    p.add_argument("--results-root",
                   default="benchmark/results/representation_dynamics")
    p.add_argument("--report-dir",
                   default="benchmark/reports/representation_dynamics")
    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    info = build_reports(Path(args.results_root), Path(args.report_dir))
    print(f"[make_representation_reports] {info}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
