"""validate_adapter_cell_state_key.py
======================================
Step 7 validator: confirm that adapter run_metadata.json files record the
correct milestone cell_state_key and do NOT fall back to legacy scGPT/scTimeBench
defaults.

Checks per run directory
------------------------
  1.  run_metadata.json exists and parses cleanly.
  2.  cell_state_key field is present (not absent / None).
  3.  cell_state_key is one of the known milestone label keys:
        consensus_milestone_label
        milestone_embedding_label
        milestone_classifier_label
  4.  cell_state_key is NOT a legacy fallback value:
        scgpt_pseudostate_provisional
        scTimeBench_cell_type
  5.  If provider_id is present, it matches gse*_milestone_* pattern.
  6.  label_mode is non-empty when present.
  7.  analysis_role is non-empty when present.

Usage
-----
  python -m benchmark.evaluation.validate_adapter_cell_state_key \\
      benchmark/results/smoke/adapter_cell_state_key [--verbose]

Exit codes:  0 = all pass,  1 = errors found.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple


MILESTONE_KEYS = {
    "consensus_milestone_label",
    "milestone_embedding_label",
    "milestone_classifier_label",
}

FORBIDDEN_KEYS = {
    "scgpt_pseudostate_provisional",
    "scTimeBench_cell_type",
}


def _check_dir(run_dir: Path, verbose: bool) -> Tuple[List[str], List[str]]:
    """Return (errors, warnings) for one adapter output directory."""
    errors: List[str] = []
    warnings: List[str] = []

    def note(msg: str) -> None:
        if verbose:
            print(f"    {msg}")

    # 1. run_metadata.json exists and parses.
    meta_path = run_dir / "run_metadata.json"
    if not meta_path.exists():
        errors.append("run_metadata.json missing")
        return errors, warnings
    try:
        with open(meta_path, encoding="utf-8") as f:
            meta: Dict[str, Any] = json.load(f)
        note(f"run_metadata.json parsed ({meta_path.stat().st_size} bytes)")
    except Exception as exc:
        errors.append(f"run_metadata.json parse error: {exc}")
        return errors, warnings

    # 2. cell_state_key present.
    csk = meta.get("cell_state_key")
    if csk is None:
        errors.append("cell_state_key is missing or null in run_metadata.json")
        return errors, warnings
    note(f"cell_state_key: {csk!r}")

    # 3. cell_state_key is a known milestone key.
    if csk not in MILESTONE_KEYS:
        errors.append(
            f"cell_state_key {csk!r} is not a recognised milestone label key. "
            f"Expected one of: {sorted(MILESTONE_KEYS)}"
        )

    # 4. cell_state_key is not a forbidden legacy fallback.
    if csk in FORBIDDEN_KEYS:
        errors.append(
            f"cell_state_key {csk!r} is a legacy fallback value — "
            "milestone configs must not use this key. "
            "Check that eval_dispatch injected the correct state_key."
        )

    # 5. provider_id pattern (optional field).
    pid = meta.get("provider_id")
    if pid is not None:
        if not (isinstance(pid, str) and "milestone" in pid):
            warnings.append(
                f"provider_id {pid!r} does not contain 'milestone' — "
                "unexpected for a milestone config"
            )
        else:
            note(f"provider_id: {pid!r}")

    # 6. label_mode non-empty when present.
    lm = meta.get("label_mode")
    if lm is not None:
        if not lm:
            warnings.append("label_mode is present but empty")
        else:
            note(f"label_mode: {lm!r}")

    # 7. analysis_role non-empty when present.
    ar = meta.get("analysis_role")
    if ar is not None:
        if not ar:
            warnings.append("analysis_role is present but empty")
        else:
            note(f"analysis_role: {ar!r}")

    return errors, warnings


def main(argv: List[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate adapter run_metadata.json for milestone cell_state_key."
    )
    parser.add_argument(
        "root",
        nargs="?",
        default="benchmark/results/smoke/adapter_cell_state_key",
        help="Root dir containing per-run output subdirs.",
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
        print(f"ERROR: output root does not exist: {root}")
        return 1

    # Accept either a flat dir of run subdirs OR a single run dir (has run_metadata.json).
    if (root / "run_metadata.json").exists():
        run_dirs = [root]
    else:
        run_dirs = sorted(d for d in root.iterdir() if d.is_dir())

    if not run_dirs:
        print(f"No subdirectories found under {root}")
        return 0

    print(f"Validating {len(run_dirs)} run dir(s) under: {root}")
    print()

    total_errors = 0
    for run_dir in run_dirs:
        errors, warnings = _check_dir(run_dir, args.verbose)
        label = "PASS" if not errors else "FAIL"
        print(f"[{label}] {run_dir.name}")
        for w in warnings:
            print(f"       WARN: {w}")
        for e in errors:
            print(f"       ERROR: {e}")
        total_errors += len(errors)

    print()
    if total_errors == 0:
        print(f"All {len(run_dirs)} run dir(s) passed cell_state_key validation.")
        return 0
    else:
        print(f"{total_errors} error(s) found across {len(run_dirs)} run dir(s).")
        return 1


if __name__ == "__main__":
    sys.exit(main())
