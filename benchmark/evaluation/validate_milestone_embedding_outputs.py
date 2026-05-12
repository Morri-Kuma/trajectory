"""validate_milestone_embedding_outputs.py
==========================================
Validate output directories produced by eval_embedding_milestone.py.

Checks per output directory
---------------------------
  1.  All four expected JSON files exist and parse.
  2.  Each mode JSON has label_mode in {consensus, embedding_based, classifier_based}.
  3.  ARI is in [-1, 1] or null with explicit reason recorded.
  4.  mean_normalized_entropy and weighted_mean_normalized_entropy in [0, 1]
      or null with explicit note.
  5.  n_cells_evaluated > 0.
  6.  state_key matches expected label mode.
  7.  provider_id matches expected label mode pattern.
  8.  provider_agreement_metrics.json contains all three pairwise comparisons.
  9.  status == "completed" for each file.

Usage
-----
  python -m benchmark.evaluation.validate_milestone_embedding_outputs \\
      benchmark/results/smoke/embedding_milestone/gse230659_observed [--verbose]

Exit codes: 0 = all pass, 1 = errors found.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple


LABEL_MODE_TO_STATE_KEY = {
    "consensus": "consensus_milestone_label",
    "embedding_based": "milestone_embedding_label",
    "classifier_based": "milestone_classifier_label",
}

EXPECTED_PAIRWISE_KEYS = {
    "consensus_vs_embedding_based",
    "consensus_vs_classifier_based",
    "embedding_based_vs_classifier_based",
}

MODE_FILES = [
    ("consensus", "embedding_metrics_consensus.json"),
    ("embedding_based", "embedding_metrics_embedding_based.json"),
    ("classifier_based", "embedding_metrics_classifier_based.json"),
]


def _is_finite_float(v: Any) -> bool:
    if v is None:
        return False
    try:
        return math.isfinite(float(v))
    except (TypeError, ValueError):
        return False


def _check_dir(out_dir: Path, verbose: bool) -> Tuple[List[str], List[str]]:
    errors: List[str] = []
    warnings: List[str] = []

    def note(msg: str) -> None:
        if verbose:
            print(f"    {msg}")

    # 1. Check all mode-specific files
    for label_mode, fname in MODE_FILES:
        fpath = out_dir / fname
        if not fpath.exists():
            errors.append(f"{fname} missing")
            continue

        try:
            with open(fpath, encoding="utf-8") as f:
                m: Dict[str, Any] = json.load(f)
            note(f"{fname} parsed ({fpath.stat().st_size} bytes)")
        except Exception as exc:
            errors.append(f"{fname} parse error: {exc}")
            continue

        # label_mode field
        got_mode = m.get("label_mode")
        if got_mode != label_mode:
            errors.append(f"{fname}: label_mode={got_mode!r}, expected {label_mode!r}")
        else:
            note(f"label_mode: {got_mode!r} OK")

        # status
        status = m.get("status", "")
        if not isinstance(status, str) or not status.startswith("completed"):
            errors.append(f"{fname}: status is not 'completed': {status!r}")
        else:
            note(f"status: {status!r}")

        # n_cells_evaluated > 0
        n_ev = m.get("n_cells_evaluated")
        if not isinstance(n_ev, int) or n_ev <= 0:
            errors.append(f"{fname}: n_cells_evaluated={n_ev!r}, expected > 0")
        else:
            note(f"n_cells_evaluated: {n_ev}")

        # ARI in [-1, 1] or null with ari_status explanation
        ari = m.get("adjusted_rand_index")
        ari_status = m.get("ari_status", "")
        if ari is None:
            if not ari_status or ari_status == "ok":
                errors.append(f"{fname}: adjusted_rand_index is null but no ari_status reason given")
            else:
                note(f"ARI: null (reason: {ari_status!r})")
        elif not _is_finite_float(ari) or not (-1.0 - 1e-9 <= float(ari) <= 1.0 + 1e-9):
            errors.append(f"{fname}: adjusted_rand_index={ari!r} is not in [-1, 1]")
        else:
            note(f"ARI: {float(ari):.4f}")

        # scTimeBench-style prediction entropy in [0, 1] or null with note.
        for ek in ("mean_prediction_entropy", "weighted_prediction_entropy"):
            ev = m.get(ek)
            if ev is None:
                note_msg = m.get("prediction_entropy_note", "")
                if not note_msg:
                    warnings.append(f"{fname}: {ek} is null (no prediction_entropy_note field)")
                else:
                    note(f"{ek}: null (note: {note_msg!r})")
            elif not _is_finite_float(ev) or not (-1e-9 <= float(ev) <= 1.0 + 1e-9):
                errors.append(f"{fname}: {ek}={ev!r} is not in [0, 1]")
            else:
                note(f"{ek}: {float(ev):.4f}")

        # Diagnostic hard-label cluster entropy in [0, 1] or null with note.
        for ek in ("mean_normalized_entropy", "weighted_mean_normalized_entropy"):
            ev = m.get(ek)
            if ev is None:
                note_msg = m.get("note", "")
                if not note_msg:
                    warnings.append(f"{fname}: {ek} is null (no note field)")
                else:
                    note(f"{ek}: null (note: {note_msg!r})")
            elif not _is_finite_float(ev) or not (-1e-9 <= float(ev) <= 1.0 + 1e-9):
                errors.append(f"{fname}: {ek}={ev!r} is not in [0, 1]")
            else:
                note(f"{ek}: {float(ev):.4f}")

        # state_key matches label_mode
        expected_sk = LABEL_MODE_TO_STATE_KEY[label_mode]
        got_sk = m.get("state_key")
        if got_sk != expected_sk:
            errors.append(
                f"{fname}: state_key={got_sk!r}, expected {expected_sk!r} for label_mode={label_mode!r}"
            )
        else:
            note(f"state_key: {got_sk!r} matches label_mode")

        # provider_id pattern check
        pid = m.get("provider_id", "")
        suffix_map = {
            "consensus": "milestone_consensus",
            "embedding_based": "milestone_embedding",
            "classifier_based": "milestone_classifier",
        }
        if suffix_map[label_mode] not in str(pid):
            warnings.append(
                f"{fname}: provider_id={pid!r} does not contain '{suffix_map[label_mode]}'"
            )
        else:
            note(f"provider_id: {pid!r}")

    # 2. Check provider_agreement_metrics.json
    agree_path = out_dir / "provider_agreement_metrics.json"
    if not agree_path.exists():
        errors.append("provider_agreement_metrics.json missing")
    else:
        try:
            with open(agree_path, encoding="utf-8") as f:
                ag: Dict[str, Any] = json.load(f)
            note(f"provider_agreement_metrics.json parsed ({agree_path.stat().st_size} bytes)")
        except Exception as exc:
            errors.append(f"provider_agreement_metrics.json parse error: {exc}")
            ag = {}

        # status
        ag_status = ag.get("status", "")
        if not isinstance(ag_status, str) or not ag_status.startswith("completed"):
            errors.append(f"provider_agreement_metrics.json: status={ag_status!r}")
        else:
            note(f"agreement status: {ag_status!r}")

        # n_cells_compared > 0
        nc = ag.get("n_cells_compared")
        if not isinstance(nc, int) or nc <= 0:
            errors.append(f"provider_agreement_metrics.json: n_cells_compared={nc!r}")
        else:
            note(f"n_cells_compared: {nc}")

        # All three pairwise ARI keys present
        pairwise = ag.get("pairwise_adjusted_rand_index") or {}
        got_pairs = set(pairwise.keys())
        missing_pairs = EXPECTED_PAIRWISE_KEYS - got_pairs
        if missing_pairs:
            errors.append(
                f"provider_agreement_metrics.json: missing pairwise ARI keys: {missing_pairs}"
            )
        else:
            note(f"pairwise ARI keys: {sorted(got_pairs)}")
            for pk, pv in pairwise.items():
                note(f"  {pk}: {pv}")

        # pairwise_exact_match_fraction keys
        pairwise_em = ag.get("pairwise_exact_match_fraction") or {}
        missing_em = EXPECTED_PAIRWISE_KEYS - set(pairwise_em.keys())
        if missing_em:
            errors.append(
                f"provider_agreement_metrics.json: missing exact_match keys: {missing_em}"
            )
        else:
            note("pairwise_exact_match_fraction: all three pairs present")

        # pairwise_crosstabs keys
        pairwise_ct = ag.get("pairwise_crosstabs") or {}
        missing_ct = EXPECTED_PAIRWISE_KEYS - set(pairwise_ct.keys())
        if missing_ct:
            errors.append(
                f"provider_agreement_metrics.json: missing crosstab keys: {missing_ct}"
            )
        else:
            note("pairwise_crosstabs: all three pairs present")

    return errors, warnings


def main(argv: List[str] | None = None) -> int:
    import argparse
    parser = argparse.ArgumentParser(
        description="Validate milestone embedding coherence output directories."
    )
    parser.add_argument(
        "root",
        nargs="?",
        default="benchmark/results/smoke/embedding_milestone/gse230659_observed",
        help="Output directory to validate (or root containing subdirs).",
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
        print(f"ERROR: output directory does not exist: {root}")
        return 1

    # Check if this directory itself has the output files (single dir mode)
    # or if it contains subdirectories (multi-dir mode).
    has_mode_files = any((root / f).exists() for _, f in MODE_FILES)
    if has_mode_files:
        dirs_to_check = [root]
    else:
        dirs_to_check = sorted(d for d in root.iterdir() if d.is_dir())
        if not dirs_to_check:
            print(f"No output files or subdirectories found under {root}")
            return 0

    print(f"Validating {len(dirs_to_check)} output dir(s):")
    print()

    total_errors = 0
    for d in dirs_to_check:
        errors, warnings = _check_dir(d, args.verbose)
        status_label = "PASS" if not errors else "FAIL"
        print(f"[{status_label}] {d.name}")
        for w in warnings:
            print(f"       WARN: {w}")
        for e in errors:
            print(f"       ERROR: {e}")
        total_errors += len(errors)

    print()
    if total_errors == 0:
        print(f"All {len(dirs_to_check)} output dir(s) passed validation.")
        return 0
    else:
        print(f"{total_errors} error(s) found across {len(dirs_to_check)} dir(s).")
        return 1


if __name__ == "__main__":
    sys.exit(main())
