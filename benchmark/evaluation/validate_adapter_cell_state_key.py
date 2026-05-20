"""Validate adapter run metadata for the official silver cell-state key."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple


EXPECTED_CELL_STATE_KEY = "final_milestone_label_coarse"


def _check_dir(run_dir: Path, verbose: bool) -> Tuple[List[str], List[str]]:
    errors: List[str] = []
    warnings: List[str] = []

    def note(msg: str) -> None:
        if verbose:
            print(f"    {msg}")

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

    csk = meta.get("cell_state_key")
    if csk != EXPECTED_CELL_STATE_KEY:
        errors.append(
            f"cell_state_key {csk!r} != {EXPECTED_CELL_STATE_KEY!r}"
        )
    else:
        note(f"cell_state_key: {csk!r}")

    pid = meta.get("provider_id")
    if pid is not None:
        if not (isinstance(pid, str) and "silver" in pid):
            warnings.append(
                f"provider_id {pid!r} does not look like a silver provider"
            )
        else:
            note(f"provider_id: {pid!r}")

    lm = meta.get("label_mode")
    if lm is not None:
        if lm != "official_silver":
            errors.append(f"label_mode {lm!r} != 'official_silver'")
        else:
            note(f"label_mode: {lm!r}")

    ar = meta.get("analysis_role")
    if ar is not None:
        if not ar:
            warnings.append("analysis_role is present but empty")
        else:
            note(f"analysis_role: {ar!r}")

    return errors, warnings


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Validate official silver cell_state_key in adapter run metadata."
    )
    parser.add_argument(
        "root",
        nargs="?",
        default="benchmark/results",
        help="Run directory or root containing per-run output subdirs.",
    )
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args(argv)

    root = Path(args.root)
    if not root.is_absolute():
        root = Path.cwd() / root

    if not root.exists():
        print(f"ERROR: output root does not exist: {root}")
        return 1

    if (root / "run_metadata.json").exists():
        run_dirs = [root]
    else:
        run_dirs = sorted(p.parent for p in root.rglob("run_metadata.json"))

    if not run_dirs:
        print(f"No run_metadata.json files found under {root}")
        return 0

    print(f"Validating {len(run_dirs)} run dir(s) under: {root}")
    total_errors = 0
    for run_dir in run_dirs:
        errors, warnings = _check_dir(run_dir, args.verbose)
        label = "PASS" if not errors else "FAIL"
        print(f"[{label}] {run_dir}")
        for warning in warnings:
            print(f"       WARN: {warning}")
        for error in errors:
            print(f"       ERROR: {error}")
        total_errors += len(errors)

    print()
    if total_errors == 0:
        print(f"All {len(run_dirs)} run dir(s) passed cell_state_key validation.")
        return 0

    print(f"{total_errors} error(s) found across {len(run_dirs)} run dir(s).")
    return 1


if __name__ == "__main__":
    sys.exit(main())
