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
  6.  AUROC, AUPRC, jaccard_similarity_topk are finite floats.
  7.  state_transition_matrix.csv exists with expected GSE230659 milestone states.
  8.  Expected states present: epithelial_like, intermediate_plastic, hCiPS.
  9.  Graph edge count == 2 (n_reference_edges in metadata).
  10. label_mode and analysis_role are non-empty strings.
  11. lineage_graph_edges.csv exists.
  12. provider_id matches gse230659_milestone_* pattern.

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


EXPECTED_STATES = {"epithelial_like", "intermediate_plastic", "hCiPS"}
EXPECTED_N_EDGES = 2


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

    # 6. Core metrics are finite floats
    for metric_key in ("auroc", "auprc", "jaccard_similarity_topk"):
        val = metrics.get(metric_key)
        if val is None or not isinstance(val, (int, float)) or not math.isfinite(float(val)):
            errors.append(f"{metric_key} is not a finite float: {val!r}")
        else:
            note(f"{metric_key}: {float(val):.4f}")

    # 7 & 8. state_transition_matrix.csv with expected states
    stm_path = out_dir / "state_transition_matrix.csv"
    if not stm_path.exists():
        errors.append("state_transition_matrix.csv missing")
    else:
        try:
            first_line = stm_path.read_text(encoding="utf-8").splitlines()[0]
            cols = {c.strip() for c in first_line.split(",")[1:] if c.strip()}
            missing_states = EXPECTED_STATES - cols
            if missing_states:
                errors.append(
                    f"state_transition_matrix.csv missing states: {missing_states}"
                )
            else:
                note(f"All expected states in STM: {sorted(EXPECTED_STATES)}")
        except Exception as exc:
            errors.append(f"state_transition_matrix.csv read error: {exc}")

    # 9. n_reference_edges == 2
    n_edges = meta.get("n_reference_edges")
    if n_edges != EXPECTED_N_EDGES:
        errors.append(f"n_reference_edges is {n_edges!r}, expected {EXPECTED_N_EDGES}")
    else:
        note(f"n_reference_edges: {n_edges} OK")

    # 10. label_mode and analysis_role non-empty
    for field in ("label_mode", "analysis_role"):
        val = meta.get(field, "")
        if not val:
            errors.append(f"metadata.{field} is empty or missing")
        else:
            note(f"{field}: {val!r}")

    # 11. lineage_graph_edges.csv exists
    edges_path = out_dir / "lineage_graph_edges.csv"
    if not edges_path.exists():
        errors.append("lineage_graph_edges.csv missing")
    else:
        note(f"lineage_graph_edges.csv present ({edges_path.stat().st_size} bytes)")

    # 12. provider_id pattern
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
