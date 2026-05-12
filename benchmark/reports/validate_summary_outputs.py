"""validate_summary_outputs.py
==============================
Step 11 validator: check that the four milestone summary CSVs produced by
summarize_milestone_results.py satisfy the structural and classification
constraints required by the benchmark framework.

Exit codes:
  0  All checks pass (or only non-structural warnings were raised)
  1  One or more structural errors found

Usage
-----
  python -m benchmark.reports.validate_summary_outputs
  python -m benchmark.reports.validate_summary_outputs --reports-dir benchmark/reports
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set


# ---------------------------------------------------------------------------
# Required columns per CSV
# ---------------------------------------------------------------------------

CORE_REQUIRED_COLS: Set[str] = {
    "dataset_id", "method", "scenario", "run_id", "result_dir",
    "metric_family", "metric_name", "metric_value",
    "provider_id", "label_mode", "cell_state_key", "reference_graph_path",
    "result_class", "formal_benchmark", "analysis_role",
    "is_primary", "is_sensitivity", "is_legacy", "source_json",
}

LINEAGE_REQUIRED_COLS: Set[str] = {
    "dataset_id", "method", "scenario", "run_id", "result_dir",
    "provider_id", "label_mode", "cell_state_key", "reference_graph_path",
    "n_reference_edges", "auroc", "auprc", "jaccard_topk", "status",
    "result_class", "formal_benchmark", "analysis_role",
    "is_primary", "is_sensitivity", "is_legacy", "source_json",
}

EMBEDDING_REQUIRED_COLS: Set[str] = {
    "dataset_id", "method", "scenario", "run_id", "result_dir",
    "provider_id", "label_mode", "state_key", "embedding_key", "cluster_source",
    "adjusted_rand_index", "mean_prediction_entropy", "weighted_prediction_entropy",
    "entropy_basis", "n_probability_classes", "prediction_entropy_note",
    "mean_normalized_entropy", "weighted_mean_normalized_entropy",
    "hard_label_entropy_note", "n_cells_evaluated", "status",
    "result_class", "formal_benchmark", "analysis_role",
    "is_primary", "is_sensitivity", "is_legacy", "source_json",
}

AGREEMENT_REQUIRED_COLS: Set[str] = {
    "dataset_id", "result_dir", "pair",
    "pairwise_adjusted_rand_index", "pairwise_exact_match_fraction",
    "n_cells_compared", "source_json",
}

PRIMARY_LABEL_MODES: Set[str] = {"consensus", "official_silver"}
SENSITIVITY_LABEL_MODES: Set[str] = {"embedding_based", "classifier_based"}
MILESTONE_LABEL_MODES: Set[str] = PRIMARY_LABEL_MODES | SENSITIVITY_LABEL_MODES


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _project_root() -> Path:
    here = Path(__file__).resolve().parent
    for c in [here, *here.parents]:
        if (c / "benchmark").exists() and (c / "benchmark" / "evaluation").exists():
            return c
    return here.parent.parent


def _truthy(val: str) -> bool:
    """Interpret CSV string as boolean True."""
    return str(val).strip().lower() in ("true", "1", "yes")


def _load_csv(path: Path) -> Optional[List[Dict[str, str]]]:
    """Load CSV into list of dicts. Returns None if file does not exist."""
    if not path.exists():
        return None
    try:
        with open(path, newline="", encoding="utf-8") as f:
            return list(csv.DictReader(f))
    except Exception as exc:
        raise RuntimeError(f"Cannot read {path}: {exc}") from exc


# ---------------------------------------------------------------------------
# Check functions
# ---------------------------------------------------------------------------

class ValidationResult:
    def __init__(self) -> None:
        self.errors: List[str] = []
        self.warnings: List[str] = []
        self.info: List[str] = []

    def error(self, msg: str) -> None:
        self.errors.append(f"  ERROR: {msg}")

    def warn(self, msg: str) -> None:
        self.warnings.append(f"  WARN:  {msg}")

    def note(self, msg: str) -> None:
        self.info.append(f"  ok:    {msg}")

    @property
    def ok(self) -> bool:
        return len(self.errors) == 0


def _check_columns(vr: ValidationResult, rows: List[Dict], required: Set[str],
                   csv_name: str) -> None:
    if not rows:
        return
    actual = set(rows[0].keys())
    missing = required - actual
    if missing:
        vr.error(
            f"{csv_name}: missing required columns: {sorted(missing)}"
        )
    else:
        vr.note(f"{csv_name}: all required columns present")


def _check_core_summary(vr: ValidationResult, rows: List[Dict]) -> None:
    """Validate core_summary.csv rows."""
    _check_columns(vr, rows, CORE_REQUIRED_COLS, "core_summary.csv")
    if not rows:
        return

    milestone_modes = MILESTONE_LABEL_MODES
    for i, row in enumerate(rows, start=2):  # row 1 = header
        lm = row.get("label_mode", "")
        pid = row.get("provider_id", "")
        is_primary = _truthy(row.get("is_primary", ""))
        is_sensitivity = _truthy(row.get("is_sensitivity", ""))
        is_legacy = _truthy(row.get("is_legacy", ""))

        # Milestone label modes must not have empty provider_id
        if lm in milestone_modes and not pid.strip():
            vr.error(
                f"core_summary.csv row {i}: label_mode={lm!r} but provider_id is empty"
            )

        # Primary rows must use a formal primary label mode.
        if is_primary and lm and lm not in PRIMARY_LABEL_MODES:
            vr.error(
                f"core_summary.csv row {i}: is_primary=True but label_mode={lm!r} "
                f"(only {sorted(PRIMARY_LABEL_MODES)!r} may be primary)"
            )

        # Sensitivity rows must only be embedding/classifier
        if is_sensitivity and lm not in ("", *SENSITIVITY_LABEL_MODES):
            if lm not in ("", "unknown"):
                vr.error(
                    f"core_summary.csv row {i}: is_sensitivity=True but "
                    f"label_mode={lm!r}"
                )

        # Legacy rows must not be primary
        if is_legacy and is_primary:
            vr.error(
                f"core_summary.csv row {i}: is_legacy=True AND is_primary=True "
                f"(legacy must not be primary)"
            )

        # scgpt provider_id must never be primary
        if "scgpt" in pid.lower() and is_primary:
            vr.error(
                f"core_summary.csv row {i}: provider_id contains 'scgpt' but "
                f"is_primary=True"
            )

        # unknown label_mode must not be primary
        if lm == "unknown" and is_primary:
            vr.error(
                f"core_summary.csv row {i}: label_mode='unknown' but is_primary=True"
            )


def _check_lineage_summary(vr: ValidationResult, rows: Optional[List[Dict]],
                            path: Path) -> None:
    if rows is None:
        vr.note("lineage_summary.csv does not exist (no lineage results yet)")
        return
    if not rows:
        vr.note("lineage_summary.csv exists but has no data rows "
                "(no milestone-tagged lineage results yet)")
        return

    _check_columns(vr, rows, LINEAGE_REQUIRED_COLS, "lineage_summary.csv")

    milestone_modes = MILESTONE_LABEL_MODES
    for i, row in enumerate(rows, start=2):
        lm = row.get("label_mode", "")
        pid = row.get("provider_id", "")
        is_primary = _truthy(row.get("is_primary", ""))
        is_legacy = _truthy(row.get("is_legacy", ""))

        if lm in milestone_modes and not pid.strip():
            vr.error(
                f"lineage_summary.csv row {i}: label_mode={lm!r} but provider_id is empty"
            )

        if is_legacy and is_primary:
            vr.error(
                f"lineage_summary.csv row {i}: is_legacy=True AND is_primary=True"
            )

        if "scgpt" in pid.lower() and is_primary:
            vr.error(
                f"lineage_summary.csv row {i}: scgpt provider_id but is_primary=True"
            )


def _check_embedding_summary(vr: ValidationResult, rows: Optional[List[Dict]],
                              path: Path) -> None:
    if rows is None:
        vr.note("embedding_summary.csv does not exist (no embedding results yet)")
        return
    if not rows:
        vr.note("embedding_summary.csv exists but has no data rows "
                "(no milestone embedding results yet)")
        return

    _check_columns(vr, rows, EMBEDDING_REQUIRED_COLS, "embedding_summary.csv")

    # Check state_key column is present (spec uses state_key not cell_state_key)
    if rows:
        actual_cols = set(rows[0].keys())
        if "state_key" not in actual_cols:
            vr.error("embedding_summary.csv: missing 'state_key' column "
                     "(must not be named 'cell_state_key')")

    milestone_modes = MILESTONE_LABEL_MODES
    for i, row in enumerate(rows, start=2):
        lm = row.get("label_mode", "")
        pid = row.get("provider_id", "")
        is_primary = _truthy(row.get("is_primary", ""))
        is_legacy = _truthy(row.get("is_legacy", ""))

        if lm in milestone_modes and not pid.strip():
            vr.error(
                f"embedding_summary.csv row {i}: label_mode={lm!r} but provider_id is empty"
            )

        if is_legacy and is_primary:
            vr.error(
                f"embedding_summary.csv row {i}: is_legacy=True AND is_primary=True"
            )

        if "scgpt" in pid.lower() and is_primary:
            vr.error(
                f"embedding_summary.csv row {i}: scgpt provider_id but is_primary=True"
            )


def _check_agreement_summary(vr: ValidationResult, rows: Optional[List[Dict]],
                              path: Path) -> None:
    if rows is None:
        vr.note("provider_agreement_summary.csv does not exist (optional)")
        return
    if not rows:
        vr.note("provider_agreement_summary.csv exists but has no data rows")
        return

    _check_columns(vr, rows, AGREEMENT_REQUIRED_COLS, "provider_agreement_summary.csv")

    # Check pairwise column names
    if rows:
        cols = set(rows[0].keys())
        for expected in ("pairwise_adjusted_rand_index", "pairwise_exact_match_fraction"):
            if expected not in cols:
                vr.error(
                    f"provider_agreement_summary.csv: missing column '{expected}'"
                )


def _check_primary_sensitivity_separation(
    vr: ValidationResult,
    core_rows: Optional[List[Dict]],
) -> None:
    """Verify that primary and sensitivity rows do not mix label_modes."""
    if not core_rows:
        return

    primary_modes: Set[str] = set()
    sensitivity_modes: Set[str] = set()

    for row in core_rows:
        lm = row.get("label_mode", "")
        if _truthy(row.get("is_primary", "")):
            primary_modes.add(lm)
        if _truthy(row.get("is_sensitivity", "")):
            sensitivity_modes.add(lm)

    bad_primary = primary_modes - PRIMARY_LABEL_MODES - {""}
    if bad_primary:
        vr.error(
            f"core_summary.csv: unexpected label_modes marked primary: {bad_primary}"
        )

    bad_sensitivity = sensitivity_modes - SENSITIVITY_LABEL_MODES - {""}
    if bad_sensitivity:
        vr.warn(
            f"core_summary.csv: unexpected label_modes in sensitivity rows: {bad_sensitivity}"
        )

    overlap = primary_modes & sensitivity_modes - {""}
    if overlap:
        vr.error(
            f"core_summary.csv: label_modes appear as BOTH primary and sensitivity: {overlap}"
        )

    if not vr.errors:
        vr.note("Primary / sensitivity label_mode separation is correct")


# ---------------------------------------------------------------------------
# Main validator
# ---------------------------------------------------------------------------

def validate(reports_dir: str) -> ValidationResult:
    vr = ValidationResult()
    root = Path(reports_dir)

    # --- core_summary ---
    core_path = root / "core_summary.csv"
    core_rows = _load_csv(core_path)
    if core_rows is None:
        vr.error("core_summary.csv does not exist — run summarize_milestone_results first")
    else:
        n = len(core_rows)
        vr.note(f"core_summary.csv: {n} data row(s)")
        _check_core_summary(vr, core_rows)
        _check_primary_sensitivity_separation(vr, core_rows)

    # --- lineage_summary ---
    lineage_path = root / "lineage_summary.csv"
    lineage_rows = _load_csv(lineage_path)
    _check_lineage_summary(vr, lineage_rows, lineage_path)
    if lineage_rows is not None:
        vr.note(f"lineage_summary.csv: {len(lineage_rows)} data row(s)")

    # --- embedding_summary ---
    emb_path = root / "embedding_summary.csv"
    emb_rows = _load_csv(emb_path)
    _check_embedding_summary(vr, emb_rows, emb_path)
    if emb_rows is not None:
        vr.note(f"embedding_summary.csv: {len(emb_rows)} data row(s)")

    # --- provider_agreement_summary ---
    agree_path = root / "provider_agreement_summary.csv"
    agree_rows = _load_csv(agree_path)
    _check_agreement_summary(vr, agree_rows, agree_path)
    if agree_rows is not None:
        vr.note(f"provider_agreement_summary.csv: {len(agree_rows)} data row(s)")

    return vr


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--reports-dir",
        default=None,
        help="Directory containing the summary CSVs. "
             "Defaults to <project_root>/benchmark/reports.",
    )
    args = parser.parse_args(argv)

    root = _project_root()
    reports_dir = args.reports_dir or str(root / "benchmark" / "reports")

    print(f"\n[validate_summary] Checking CSVs in: {reports_dir}")
    print(f"{'=' * 60}")

    try:
        vr = validate(reports_dir)
    except RuntimeError as exc:
        print(f"\nFATAL: {exc}", file=sys.stderr)
        return 1

    for msg in vr.info:
        print(msg)
    for msg in vr.warnings:
        print(msg)
    for msg in vr.errors:
        print(msg, file=sys.stderr)

    print(f"{'=' * 60}")
    if vr.ok:
        print(f"[validate_summary] PASSED ({len(vr.info)} checks, "
              f"{len(vr.warnings)} warnings, 0 errors)\n")
        return 0
    else:
        print(
            f"[validate_summary] FAILED: {len(vr.errors)} error(s), "
            f"{len(vr.warnings)} warning(s)\n",
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
