"""
summarize_lineage.py
--------------------
Minimal Lineage Fidelity summary / ranking script.

Reads one or more ``lineage_metrics.json`` files produced by
``benchmark/evaluation/eval_lineage.py`` and produces a compact
summary table. Rank aggregation is computed ONLY within Lineage Fidelity 鈥?this script never pretends methods have an overall benchmark score
(WOT and CellRank2 are OT-style lineage methods and are not eligible for
Forecast Accuracy or Embedding Coherence; see the framework doc 搂6).

Official vs diagnostic metric classification (v2, 2026-04-20)
-------------------------------------------------------------
Official Lineage Fidelity rank metrics (framework v2 搂9.5):
  - AUROC
  - AUPRC
  - Jaccard Similarity
  - Single-step recovery
  - Multi-step recovery

Diagnostic / supplementary (NOT part of the official rank):
  - jaccard_similarity_topk   鈥?useful diagnostic for dense-STM adapters
                                 (see benchmark/docs/jaccard_singlestep_diagnosis.md)
                                 but not framework-defined. Shown in the table,
                                 never summed into lineage_rank.

Scenario-level baseline
-----------------------
The baseline uses whatever scenario-filtered ``adata`` is handed to the
evaluator and computes a Spearman/maximum vote-based pseudo-STM against the
same reference graph. If WOT and CellRank2 are evaluated on the same adata,
their baselines are identical by
construction. This script emits ONE baseline row per scenario (not one per
method). The baseline is picked from the highest-confidence method run in the
scenario; any disagreement between method-level baseline copies is reported
(``baseline_consistent`` flag) instead of silently de-duplicated.

Result classification
---------------------
Each row carries a ``result_class`` column:
  - ``official`` 鈥?full-data run, countable as a formal benchmark result.
  - ``pilot``    鈥?subsampled / debug / partial run. Reported separately.
  - ``reduced_validation`` 鈥?reduced feature/training setting used to validate
    a new method or dimension before formal runs.
  - ``smoke`` / ``hpc_validation`` / ``pilot_backup`` / ``diagnostic`` 鈥?
    non-formal runs retained for traceability but excluded by
    ``--official-only``.
The pilot/official classification is taken from the run's
``run_metadata.json`` ``result_class`` field if present. If absent, it is
inferred: any run whose ``sampling.subsample_per_timepoint`` is set is
classified as ``pilot``.

Usage
-----
    # Default discovery: scan every benchmark/results/*/scenario_*
    python benchmark/reports/summarize_lineage.py

    # Explicit list
    python benchmark/reports/summarize_lineage.py \
        --inputs \
            wot=A_scgpt_v1=benchmark/results/wot/scenario_A_scgpt_v1 \
            cellrank2=A_scgpt_v1=benchmark/results/cellrank2/scenario_A_scgpt_v1 \
            wot=B_scgpt_v1=benchmark/results/wot/scenario_B_scgpt_v1 \
        --output benchmark/reports/lineage_summary.csv
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Metric classification 鈥?see module docstring.
# ---------------------------------------------------------------------------

# Official Lineage Fidelity rank metric set (framework v2 搂9.5).
# Do NOT silently add/remove from this tuple 鈥?it defines the benchmark.
RANK_METRICS = (
    "auroc",
    "auprc",
    "jaccard_similarity",
    "single_step_recovery",
    "multi_step_recovery",
)

# Diagnostic columns shown in the table but NOT used to compute lineage_rank.
SUPPLEMENTARY_METRICS = (
    "jaccard_similarity_topk",
)


VALID_RESULT_CLASSES = {
    "official",
    "pilot",
    "pilot_backup",
    "reduced_validation",
    "smoke",
    "hpc_validation",
    "diagnostic",
}


# ---------------------------------------------------------------------------
# Path discovery
# ---------------------------------------------------------------------------

def _project_root() -> Path:
    here = Path(__file__).resolve().parent
    for c in [here, *here.parents]:
        if (c / "benchmark").exists() and (c / "benchmark" / "evaluation").exists():
            return c
    return here.parent.parent


def _discover_runs(root: Path) -> List[Tuple[str, str, Path]]:
    """
    Auto-discover (method, scenario, run_dir) triples by scanning
    ``benchmark/results/<method>/scenario_<id>*/lineage_metrics.json``.
    """
    runs: List[Tuple[str, str, Path]] = []
    results_root = root / "benchmark" / "results"
    if not results_root.exists():
        return runs
    for method_dir in sorted(p for p in results_root.iterdir() if p.is_dir()):
        method = method_dir.name
        for run_dir in sorted(p for p in method_dir.iterdir() if p.is_dir()):
            metrics_file = run_dir / "lineage_metrics.json"
            if not metrics_file.exists():
                continue
            name = run_dir.name
            # Preserve the full scenario id including state-system suffix
            # (e.g. "A_scgpt_v1") so that legacy "A" runs don't mix with
            # scGPT-v1 results in the per-scenario ranking bucket.
            scenario = "?"
            if name.startswith("scenario_") and len(name) > len("scenario_"):
                scenario = name[len("scenario_"):]
            runs.append((method, scenario, run_dir))
    return runs


# ---------------------------------------------------------------------------
# File readers
# ---------------------------------------------------------------------------

def _read_metrics(run_dir: Path) -> Optional[Dict]:
    mpath = run_dir / "lineage_metrics.json"
    if not mpath.exists():
        return None
    try:
        with open(mpath, encoding="utf-8") as f:
            return json.load(f)
    except Exception as exc:
        print(f"[summarize] Failed to read {mpath}: {exc}", file=sys.stderr)
        return None


def _read_run_metadata(run_dir: Path) -> Dict:
    """Read run_metadata.json; return {} if absent or unreadable."""
    mpath = run_dir / "run_metadata.json"
    if not mpath.exists():
        return {}
    try:
        with open(mpath, encoding="utf-8") as f:
            return json.load(f) or {}
    except Exception as exc:
        print(f"[summarize] Failed to read {mpath}: {exc}", file=sys.stderr)
        return {}


def _infer_result_class(meta: Dict) -> str:
    """
    Pull result_class from run_metadata if the runner declared it.
    Otherwise infer: anything with a non-null subsample is pilot; else official.
    """
    declared = meta.get("result_class")
    if declared in VALID_RESULT_CLASSES:
        return declared
    sampling = meta.get("sampling") or {}
    if sampling.get("subsample_per_timepoint") not in (None, "", False):
        return "pilot"
    return "official"


def _ground_truth_info(metrics: Dict, meta: Optional[Dict] = None) -> Dict:
    """Extract ground-truth provider fields from metrics or run metadata."""
    gt = metrics.get("ground_truth") or (meta or {}).get("ground_truth") or {}
    if not isinstance(gt, dict):
        gt = {}
    return {
        "ground_truth_provider": gt.get("provider_id"),
        "ground_truth_status": gt.get("status"),
        "ground_truth_state_key": gt.get("state_key") or metrics.get("cell_state_key"),
        "ground_truth_reference_graph": gt.get("reference_graph_path"),
    }


# ---------------------------------------------------------------------------
# Row construction
# ---------------------------------------------------------------------------

def _method_row(method: str, scenario: str, run_dir: Path,
                metrics: Dict, result_class: str,
                sampling: Optional[Dict]) -> Dict:
    """One row per (method, scenario). Baseline is handled separately."""
    row = {
        "kind": "method",
        "result_class": result_class,
        "method": method,
        "scenario": scenario,
        "auroc": metrics.get("auroc"),
        "auprc": metrics.get("auprc"),
        "jaccard_similarity": metrics.get("jaccard_similarity"),
        "single_step_recovery": metrics.get("single_step_recovery"),
        "multi_step_recovery": metrics.get("multi_step_recovery"),
        "jaccard_similarity_topk": metrics.get("jaccard_similarity_topk"),
        "n_reference_edges": metrics.get("n_reference_edges"),
        "edge_confidence_mode": metrics.get("edge_confidence_mode"),
        "status": metrics.get("status"),
        "subsample_per_timepoint": (sampling or {}).get("subsample_per_timepoint"),
        "source": str(run_dir),
    }
    row.update(_ground_truth_info(metrics))
    return row


def _baseline_row_from(scenario: str, result_class: str,
                       source_run_dir: Path,
                       baseline_block: Dict,
                       metrics_envelope: Dict,
                       sampling: Optional[Dict]) -> Dict:
    """Build a scenario-level baseline row from one method's embedded baseline."""
    row = {
        "kind": "baseline",
        "result_class": result_class,
        "method": "correlation_baseline",
        "scenario": scenario,
        "auroc": baseline_block.get("auroc"),
        "auprc": baseline_block.get("auprc"),
        "jaccard_similarity": baseline_block.get("jaccard_similarity"),
        "single_step_recovery": baseline_block.get("single_step_recovery"),
        "multi_step_recovery": baseline_block.get("multi_step_recovery"),
        # The formal baseline uses thresholded graph Jaccard; top-k Jaccard
        # remains a method-output diagnostic column, so leave it blank here.
        "jaccard_similarity_topk": None,
        "n_reference_edges": metrics_envelope.get("n_reference_edges"),
        "edge_confidence_mode": metrics_envelope.get("edge_confidence_mode"),
        "status": "baseline",
        "subsample_per_timepoint": (sampling or {}).get("subsample_per_timepoint"),
        "source": str(source_run_dir),
    }
    row.update(_ground_truth_info(metrics_envelope))
    return row


def _baseline_signature(b: Dict) -> Tuple:
    """
    Stable equality signature for a baseline block. Used to detect
    disagreement between the baseline copies stored in different method runs
    of the same scenario. Rounded to 10 decimal places to absorb harmless
    float jitter.
    """
    def _r(v):
        try:
            return round(float(v), 10) if v is not None else None
        except (TypeError, ValueError):
            return v
    return (
        _r(b.get("auroc")),
        _r(b.get("auprc")),
        _r(b.get("jaccard_similarity")),
        _r(b.get("single_step_recovery")),
        _r(b.get("multi_step_recovery")),
        b.get("cell_state_key"),
        b.get("time_key"),
        b.get("baseline_style"),
        b.get("correlation_method"),
        b.get("averaging_method"),
        (b.get("ground_truth") or {}).get("provider_id"),
        b.get("n_states_used"),
        b.get("n_timepoint_pairs_used"),
        b.get("n_source_cell_votes"),
    )


def _collapse_baselines(per_run: List[Dict]) -> List[Dict]:
    """
    Collapse N method-level baseline entries (one per (method, scenario) run)
    into one scenario-level baseline row per scenario.

    Input entries are the dicts returned by ``_assemble_run`` below, carrying
    the per-run baseline block plus scenario metadata. Entries without a
    baseline are ignored here.
    """
    by_scenario: Dict[Tuple[str, Optional[str]], List[Dict]] = {}
    for e in per_run:
        if not e.get("baseline"):
            continue
        gt_provider = (
            (e.get("metrics") or {}).get("ground_truth") or {}
        ).get("provider_id")
        by_scenario.setdefault((e["scenario"], gt_provider), []).append(e)

    rows: List[Dict] = []
    for (scenario, _gt_provider), entries in by_scenario.items():
        # Partition by result_class so pilot and official baselines don't
        # collapse into one another (a pilot baseline was computed on a
        # subsampled training split and is not comparable to the official
        # baseline).
        by_class: Dict[str, List[Dict]] = {}
        for e in entries:
            by_class.setdefault(e["result_class"], []).append(e)

        for rc, items in by_class.items():
            # Check whether the baselines agree numerically. Keep the first
            # entry's baseline as the canonical one and record consistency.
            sigs = {_baseline_signature(it["baseline"]) for it in items}
            canonical = items[0]
            row = _baseline_row_from(
                scenario=scenario,
                result_class=rc,
                source_run_dir=canonical["run_dir"],
                baseline_block=canonical["baseline"],
                metrics_envelope=canonical["metrics"],
                sampling=canonical["sampling"],
            )
            row["baseline_consistent"] = (len(sigs) == 1)
            row["n_method_runs_contributing"] = len(items)
            if len(sigs) > 1:
                # Surface the disagreement in the status field so readers
                # don't silently trust a collapsed row.
                row["status"] = (
                    f"baseline (INCONSISTENT across {len(items)} runs; "
                    f"{len(sigs)} distinct signatures 鈥?see per-run JSON)"
                )
            rows.append(row)
    return rows


def _assemble_run(method: str, scenario: str, run_dir: Path) -> Optional[Dict]:
    """
    Parse one run directory into a bundle of info used downstream:
      - method row
      - embedded baseline block (if any)
      - run metadata (for result_class / sampling)
    """
    metrics = _read_metrics(run_dir)
    if metrics is None:
        return None
    meta = _read_run_metadata(run_dir)
    sampling = meta.get("sampling") or None
    result_class = _infer_result_class(meta)
    return {
        "method": method,
        "scenario": scenario,
        "run_dir": run_dir,
        "metrics": metrics,
        "baseline": metrics.get("baseline") or None,
        "meta": meta,
        "sampling": sampling,
        "result_class": result_class,
        "method_row": _method_row(method, scenario, run_dir,
                                  metrics, result_class, sampling),
    }


# ---------------------------------------------------------------------------
# Ranking
# ---------------------------------------------------------------------------

def _rank_within_group(rows: List[Dict]) -> None:
    """
    Assign ``lineage_rank`` per row, computed separately within each
    (scenario, result_class, ground_truth_provider) bucket. Official and pilot
    rows are ranked in their own buckets, and annotation providers are kept
    separate so ground-truth sensitivity analyses do not mix score scales.

    Rank is the average of per-metric ranks across RANK_METRICS. Lower = better.
    Rows with all-None metrics get ``lineage_rank = None``.
    Per-metric ranking uses dense ranking with tie handling (average of tied
    ranks, 1-based, higher metric value = better).
    """
    by_group: Dict[Tuple[str, str, Optional[str]], List[Dict]] = {}
    for r in rows:
        key = (
            r["scenario"],
            r.get("result_class", "official"),
            r.get("ground_truth_provider"),
        )
        by_group.setdefault(key, []).append(r)

    for group_rows in by_group.values():
        rank_sums = {id(r): 0.0 for r in group_rows}
        rank_counts = {id(r): 0 for r in group_rows}
        for metric in RANK_METRICS:
            values: List[Tuple[float, int]] = []
            for i, r in enumerate(group_rows):
                v = r.get(metric)
                if v is None:
                    continue
                try:
                    fv = float(v)
                except (TypeError, ValueError):
                    continue
                if fv != fv:  # NaN guard
                    continue
                values.append((fv, i))
            if not values:
                continue
            # Rank 1 = largest metric value.
            values.sort(key=lambda pair: pair[0], reverse=True)
            rank_by_index: Dict[int, float] = {}
            j = 0
            while j < len(values):
                k = j
                while k + 1 < len(values) and values[k + 1][0] == values[j][0]:
                    k += 1
                avg_rank = (j + k + 2) / 2.0
                for q in range(j, k + 1):
                    rank_by_index[values[q][1]] = avg_rank
                j = k + 1
            for i, r in enumerate(group_rows):
                rk = rank_by_index.get(i)
                if rk is None:
                    continue
                rank_sums[id(r)] += rk
                rank_counts[id(r)] += 1

        for r in group_rows:
            n = rank_counts[id(r)]
            r["lineage_rank"] = (rank_sums[id(r)] / n) if n > 0 else None
            r["n_metrics_ranked"] = n


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def _write_csv(rows: List[Dict], out_path: Path) -> None:
    import csv
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "scenario", "result_class", "kind", "method",
        "ground_truth_provider", "ground_truth_status",
        "ground_truth_state_key", "ground_truth_reference_graph",
        "auroc", "auprc", "jaccard_similarity",
        "single_step_recovery", "multi_step_recovery",
        "lineage_rank", "n_metrics_ranked",
        "jaccard_similarity_topk",   # supplementary 鈥?after the rank columns
        "n_reference_edges", "edge_confidence_mode",
        "subsample_per_timepoint",
        "status", "baseline_consistent", "n_method_runs_contributing",
        "source",
    ]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow({k: r.get(k) for k in fieldnames})


def _fmt(v) -> str:
    if v is None:
        return "-"
    if isinstance(v, float):
        if v != v:
            return "-"
        return f"{v:.4f}"
    return str(v)


def _print_table(rows: List[Dict]) -> None:
    if not rows:
        print("[summarize] No rows to print.")
        return
    # Two sub-tables, one per result_class, so the reader cannot confuse
    # pilot rows with official rows.
    by_class: Dict[str, List[Dict]] = {}
    for r in rows:
        by_class.setdefault(r.get("result_class", "official"), []).append(r)

    cols = ["scenario", "ground_truth_provider", "method", "kind",
            "auroc", "auprc", "jaccard_similarity",
            "single_step_recovery", "multi_step_recovery",
            "lineage_rank",
            "jaccard_similarity_topk"]

    # Print formal block first, then diagnostic classes.
    for rc in (
        "official",
        "reduced_validation",
        "pilot",
        "pilot_backup",
        "smoke",
        "hpc_validation",
        "diagnostic",
    ):
        group = by_class.get(rc, [])
        if not group:
            continue
        print(f"\n=== {rc.upper()} ===")
        widths = {c: max(len(c), max(len(_fmt(r.get(c))) for r in group)) for c in cols}
        print("  ".join(c.ljust(widths[c]) for c in cols))
        print("  ".join("-" * widths[c] for c in cols))
        sorted_rows = sorted(
            group,
            key=lambda r: (
                str(r.get("scenario")),
                0 if r.get("kind") == "baseline" else 1,
                float(r.get("lineage_rank"))
                if isinstance(r.get("lineage_rank"), (int, float)) else 99.0,
                str(r.get("method")),
            ),
        )
        for r in sorted_rows:
            print("  ".join(_fmt(r.get(c)).ljust(widths[c]) for c in cols))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--inputs", nargs="*", default=None,
        help=(
            "Optional explicit list of method=scenario=output_dir entries. "
            "If omitted, auto-discovers every lineage_metrics.json under "
            "benchmark/results/*/scenario_*/."
        ),
    )
    parser.add_argument(
        "--output", default="benchmark/reports/lineage_summary.csv",
        help="Path for the CSV summary output.",
    )
    parser.add_argument(
        "--no-baseline", action="store_true",
        help="Skip the scenario-level correlation-baseline rows.",
    )
    parser.add_argument(
        "--official-only", action="store_true",
        help="Only include runs classified as official (exclude pilots).",
    )
    args = parser.parse_args()

    root = _project_root()

    if args.inputs:
        runs: List[Tuple[str, str, Path]] = []
        for entry in args.inputs:
            parts = entry.split("=", 2)
            if len(parts) != 3:
                raise SystemExit(
                    f"Bad --inputs entry {entry!r}; expected method=scenario=output_dir"
                )
            method, scenario, out_dir = parts
            out_p = Path(out_dir)
            if not out_p.is_absolute():
                out_p = root / out_p
            runs.append((method, scenario, out_p))
    else:
        runs = _discover_runs(root)

    if not runs:
        print("[summarize] No runs found. Nothing to do.")
        return

    # Assemble per-run info bundles.
    bundles: List[Dict] = []
    for method, scenario, run_dir in runs:
        b = _assemble_run(method, scenario, run_dir)
        if b is None:
            continue
        if args.official_only and b["result_class"] != "official":
            continue
        bundles.append(b)

    # Method rows.
    method_rows: List[Dict] = [b["method_row"] for b in bundles]

    # Scenario-level baseline rows (one per scenario 脳 result_class).
    baseline_rows: List[Dict] = []
    if not args.no_baseline:
        baseline_rows = _collapse_baselines(bundles)

    all_rows = method_rows + baseline_rows

    # Rank within (scenario, result_class). Baseline rows participate in the
    # rank so readers can see whether a method actually beat the baseline.
    _rank_within_group(all_rows)

    out_path = Path(args.output)
    if not out_path.is_absolute():
        out_path = root / out_path
    _write_csv(all_rows, out_path)

    _print_table(all_rows)
    print(f"\n[summarize] Wrote {len(all_rows)} rows 鈫?{out_path}")


if __name__ == "__main__":
    main()
