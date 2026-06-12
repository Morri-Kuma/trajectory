"""validate_milestone_lineage_outputs.py
==========================================
Validate smoke lineage output directories from run_milestone_lineage_smoke.py.

Checks per output directory
---------------------------
  1.  smoke_metadata.json exists and parses cleanly.
  2.  smoke_test == true.
  3.  formal_benchmark == false.
  4.  lineage_metrics.json exists and parses cleanly.
  5.  status == "completed".
  6.  metric_protocol == "sctimebench_graph_sim".
  7.  graph_metrics.single_step: auc_roc, auc_prc, jaccard_similarity are finite.
  8.  graph_metrics.multi_step:  auc_roc, auc_prc, jaccard_similarity are finite.
  9.  lineage metrics use final_milestone_label_coarse from the provider contract.
  10. state_transition_matrix.csv uses coarse labels; missing coarse states are
      zero-filled by the evaluator and recorded in prediction_label_report.
  11. Graph edge count == 2 (n_reference_edges in metadata).
  12. label_mode and analysis_role are non-empty strings.
  13. lineage_graph_edges.csv exists.
  14. provider_id matches gse230659_milestone_* pattern.

Usage
-----
  python benchmark/evaluation/validate_milestone_lineage_outputs.py \
      benchmark/results/smoke/milestone_lineage [--verbose]

Exit codes: 0 = all pass, 1 = errors found.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple

from benchmark.evaluation.lineage_graphsim_sctimebench import (
    PREDICTION_LABEL_SOURCE_SCTIMEBENCH_ZERO_FILL,
)


EXPECTED_STATES = {"epithelial_like", "intermediate_plastic", "hCiPS"}
EXPECTED_STATE_KEY = "final_milestone_label_coarse"
EXPECTED_N_EDGES = 2
REQUIRED_LINEAGE_PROTOCOL = "sctimebench_graph_sim"
REQUIRED_PREDICTION_LABEL_SOURCE = PREDICTION_LABEL_SOURCE_SCTIMEBENCH_ZERO_FILL

# Fields required to be finite floats under graph_metrics.{criterion}
_REQUIRED_GRAPH_METRIC_FIELDS = ("auc_roc", "auc_prc", "jaccard_similarity")


def _check_dir(out_dir: Path, verbose: bool) -> Tuple[List[str], List[str]]:
    """Return (errors, warnings) for one smoke output directory."""
    errors: List[str] = []
    warnings: List[str] = []

    def note(msg: str) -> None:
        if verbose:
            print(f"    {msg}")

    # 1. smoke_metadata.json
    meta_path = out_dir / "smoke_metadata.json"
    if not meta_path.exists():
        errors.append("smoke_metadata.json missing")
        return errors, warnings
    try:
        with open(meta_path, encoding="utf-8") as f:
            meta: Dict[str, Any] = json.load(f)
        note(f"smoke_metadata.json parsed ({meta_path.stat().st_size} bytes)")
    except Exception as exc:
        errors.append(f"smoke_metadata.json parse error: {exc}")
        return errors, warnings

    # 2. smoke_test == true
    if meta.get("smoke_test") is not True:
        errors.append(f"smoke_test != true (got {meta.get('smoke_test')!r})")
    else:
        note("smoke_test: true OK")

    # 3. formal_benchmark == false
    if meta.get("formal_benchmark") is not False:
        errors.append(f"formal_benchmark != false (got {meta.get('formal_benchmark')!r})")
    else:
        note("formal_benchmark: false OK")

    # 4. lineage_metrics.json
    metrics_path = out_dir / "lineage_metrics.json"
    if not metrics_path.exists():
        errors.append("lineage_metrics.json missing")
        metrics: Dict[str, Any] = {}
    else:
        try:
            with open(metrics_path, encoding="utf-8") as f:
                metrics = json.load(f)
            note(f"lineage_metrics.json parsed ({metrics_path.stat().st_size} bytes)")
        except Exception as exc:
            errors.append(f"lineage_metrics.json parse error: {exc}")
            metrics = {}

    # 5. Status == completed
    status = metrics.get("status", "")
    if not isinstance(status, str) or not status.startswith("completed"):
        errors.append(f"lineage_metrics status is not completed: {status!r}")
    else:
        note(f"status: {status}")

    # 6. metric_protocol must be sctimebench_graph_sim
    protocol = metrics.get("metric_protocol")
    if protocol != REQUIRED_LINEAGE_PROTOCOL:
        errors.append(
            f"metric_protocol={protocol!r}, expected {REQUIRED_LINEAGE_PROTOCOL!r}. "
            "Re-run eval_lineage.py to generate scTimeBench graph-sim metrics."
        )
    else:
        note(f"metric_protocol: {protocol!r} OK")

    # 7 & 8. graph_metrics.{single_step,multi_step} must have finite AUC/Jaccard
    gm = metrics.get("graph_metrics") or {}
    for criterion_key, label in (("single_step", "simple"), ("multi_step", "all_paths")):
        criterion = gm.get(criterion_key) or {}
        for field in _REQUIRED_GRAPH_METRIC_FIELDS:
            val = criterion.get(field)
            if val is None or not isinstance(val, (int, float)) or not math.isfinite(float(val)):
                errors.append(
                    f"graph_metrics.{criterion_key}.{field} is not a finite float: {val!r}"
                )
            else:
                note(f"graph_metrics.{criterion_key}.{field}: {float(val):.4f}")

    # 9. State-key contract: official labels are canonical, not STM-derived.
    metrics_state_key = (
        (metrics.get("ground_truth") or {}).get("state_key")
        or metrics.get("cell_state_key")
    )
    prediction_state_key = metrics.get("prediction_state_key")
    if metrics_state_key != EXPECTED_STATE_KEY:
        errors.append(
            f"lineage state_key={metrics_state_key!r}, expected {EXPECTED_STATE_KEY!r}"
        )
    else:
        note(f"lineage state_key: {metrics_state_key!r} OK")
    if prediction_state_key not in (None, EXPECTED_STATE_KEY):
        errors.append(
            f"prediction_state_key={prediction_state_key!r}, expected "
            f"{EXPECTED_STATE_KEY!r}"
        )
    label_source = metrics.get("prediction_label_source")
    if label_source not in (None, REQUIRED_PREDICTION_LABEL_SOURCE):
        errors.append(
            "prediction_label_source is not scTimeBench zero-fill/coarse: "
            f"{label_source!r}"
        )

    # 10. state_transition_matrix.csv exists and uses coarse states.  Missing
    # reference nodes are allowed only because the evaluator zero-fills them.
    stm_path = out_dir / "state_transition_matrix.csv"
    if not stm_path.exists():
        errors.append("state_transition_matrix.csv missing")
    else:
        try:
            first_line = stm_path.read_text(encoding="utf-8").splitlines()[0]
            cols = {c.strip() for c in first_line.split(",")[1:] if c.strip()}
            expanded_cols = sorted(c for c in cols if c.startswith("stage_"))
            if expanded_cols:
                errors.append(
                    "state_transition_matrix.csv uses expanded/stage labels: "
                    f"{expanded_cols[:20]}"
                )
            missing_states = EXPECTED_STATES - cols
            if missing_states:
                report = metrics.get("prediction_label_report") or {}
                reported_missing = set(report.get("missing_reference_nodes") or [])
                if not missing_states <= reported_missing:
                    errors.append(
                        "state_transition_matrix.csv missing states but "
                        "prediction_label_report does not record matching "
                        f"zero-filled reference nodes: {sorted(missing_states)}"
                    )
                else:
                    note(
                        "Missing STM states zero-filled by evaluator: "
                        f"{sorted(missing_states)}"
                    )
            else:
                note(f"All expected states in STM: {sorted(EXPECTED_STATES)}")
        except Exception as exc:
            errors.append(f"state_transition_matrix.csv read error: {exc}")

    # 11. n_reference_edges == 2
    n_edges = meta.get("n_reference_edges")
    if n_edges != EXPECTED_N_EDGES:
        errors.append(f"n_reference_edges is {n_edges!r}, expected {EXPECTED_N_EDGES}")
    else:
        note(f"n_reference_edges: {n_edges} OK")

    # 12. label_mode and analysis_role non-empty
    for field in ("label_mode", "analysis_role"):
        val = meta.get(field, "")
        if not val:
            errors.append(f"metadata.{field} is empty or missing")
        else:
            note(f"{field}: {val!r}")

    # 13. lineage_graph_edges.csv exists (diagnostic artifact, not metric source)
    edges_path = out_dir / "lineage_graph_edges.csv"
    if not edges_path.exists():
        errors.append("lineage_graph_edges.csv missing")
    else:
        note(f"lineage_graph_edges.csv present ({edges_path.stat().st_size} bytes)")

    # 14. provider_id pattern
    pid = meta.get("provider_id", "")
    if not pid.startswith("gse230659_milestone_"):
        warnings.append(
            f"provider_id {pid!r} does not match gse230659_milestone_* pattern"
        )
    else:
        note(f"provider_id: {pid!r}")

    return errors, warnings


def main(argv: List[str] | None = None) -> int:
    import argparse
    parser = argparse.ArgumentParser(
        description="Validate milestone lineage smoke output directories."
    )
    parser.add_argument(
        "root",
        nargs="?",
        default="benchmark/results/smoke/milestone_lineage",
        help="Root containing per-run smoke output subdirs.",
    )
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args(argv)

    root = Path(args.root)
    if not root.is_absolute():
        here = Path.cwd()
        candidate = here / root
        if not candidate.exists():
            for anc in [here, *here.parents]:
                if (anc / "benchmark").exists():
                    candidate = anc / root
                    break
        root = candidate

    if not root.exists():
        print(f"ERROR: smoke output root does not exist: {root}")
        return 1

    run_dirs = sorted(d for d in root.iterdir() if d.is_dir())
    if not run_dirs:
        print(f"No subdirectories found under {root}")
        return 0

    print(f"Validating {len(run_dirs)} smoke output dir(s) under: {root}")
    print()

    total_errors = 0
    for run_dir in run_dirs:
        errors, warnings = _check_dir(run_dir, args.verbose)
        status_label = "PASS" if not errors else "FAIL"
        print(f"[{status_label}] {run_dir.name}")
        for w in warnings:
            print(f"       WARN: {w}")
        for e in errors:
            print(f"       ERROR: {e}")
        total_errors += len(errors)

    print()
    if total_errors == 0:
        print(f"All {len(run_dirs)} smoke output dir(s) passed validation.")
        return 0
    else:
        print(f"{total_errors} error(s) found across {len(run_dirs)} dir(s).")
        return 1


if __name__ == "__main__":
    sys.exit(main())
