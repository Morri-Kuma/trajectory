"""validate_legacy_scgpt_policy.py
===================================
Step 12 validator: enforce that all scGPT pseudostate providers are correctly
marked as legacy and that no scGPT result is classified as primary in any
summary CSV.

Checks
------
Registry checks (benchmark/ground_truth/registry.yaml):
  - Every provider whose provider_id contains 'scgpt' has:
      deprecated: true
      legacy: true
      label_mode: legacy_scgpt_pseudostate
      analysis_role: legacy_reference_only
      status containing 'legacy' or 'deprecated'
  - No milestone consensus provider is marked deprecated or legacy.

Summary CSV checks (benchmark/reports/):
  - core_summary.csv: no row with provider_id containing 'scgpt' has is_primary=True
  - core_summary.csv: no row with label_mode == legacy_scgpt_pseudostate has is_primary=True
  - core_summary.csv: no row with cell_state_key == scgpt_pseudostate_provisional has is_primary=True
  - lineage_summary.csv: same three checks
  - embedding_summary.csv: same, using state_key instead of cell_state_key
  - Primary rows (is_primary=True) must have:
      label_mode in {consensus, official_silver}
      provider_id containing _milestone_consensus_ or _marker_fm_transition_silver_
      is_legacy == False

Exit codes
----------
  0  All checks pass
  1  One or more structural errors found

Usage
-----
  python -m benchmark.reports.validate_legacy_scgpt_policy
  python -m benchmark.reports.validate_legacy_scgpt_policy --registry benchmark/ground_truth/registry.yaml --reports-dir benchmark/reports
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set

try:
    import yaml
    _YAML_AVAILABLE = True
except ImportError:
    _YAML_AVAILABLE = False


PRIMARY_LABEL_MODES: Set[str] = {"consensus", "official_silver"}
PRIMARY_PROVIDER_MARKERS: Set[str] = {
    "_milestone_consensus_",
    "_marker_fm_transition_silver_",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _project_root() -> Path:
    here = Path(__file__).resolve().parent
    for c in [here, *here.parents]:
        if (c / "benchmark").exists() and (c / "benchmark" / "evaluation").exists():
            return c
    return here.parent.parent


def _truthy(val) -> bool:
    return str(val).strip().lower() in ("true", "1", "yes")


def _load_csv(path: Path) -> Optional[List[Dict[str, str]]]:
    if not path.exists():
        return None
    try:
        with open(path, newline="", encoding="utf-8") as f:
            return list(csv.DictReader(f))
    except Exception as exc:
        raise RuntimeError(f"Cannot read {path}: {exc}") from exc


def _load_registry(path: Path) -> Dict:
    if not _YAML_AVAILABLE:
        raise RuntimeError(
            "PyYAML is not installed. Install it with: pip install pyyaml"
        )
    with open(path, encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return data.get("providers", {})


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


# ---------------------------------------------------------------------------
# Registry checks
# ---------------------------------------------------------------------------

def _check_registry(vr: ValidationResult, providers: Dict) -> None:
    """Verify all scGPT providers are correctly flagged as legacy."""
    scgpt_providers = {
        pid: entry for pid, entry in providers.items()
        if "scgpt" in pid.lower()
    }
    milestone_providers = {
        pid: entry for pid, entry in providers.items()
        if "_milestone_" in pid.lower()
    }

    if not scgpt_providers:
        vr.warn("No scGPT providers found in registry — nothing to check")
    else:
        vr.note(f"Found {len(scgpt_providers)} scGPT provider(s) in registry: "
                f"{sorted(scgpt_providers)}")

    for pid, entry in scgpt_providers.items():
        # deprecated: true
        if not _truthy(entry.get("deprecated", False)):
            vr.error(f"registry provider '{pid}': deprecated is not True "
                     f"(got {entry.get('deprecated')!r})")
        # legacy: true
        if not _truthy(entry.get("legacy", False)):
            vr.error(f"registry provider '{pid}': legacy is not True "
                     f"(got {entry.get('legacy')!r})")
        # label_mode: legacy_scgpt_pseudostate
        lm = entry.get("label_mode", "")
        if lm != "legacy_scgpt_pseudostate":
            vr.error(f"registry provider '{pid}': label_mode={lm!r} "
                     f"(expected 'legacy_scgpt_pseudostate')")
        # analysis_role: legacy_reference_only
        ar = entry.get("analysis_role", "")
        if ar != "legacy_reference_only":
            vr.error(f"registry provider '{pid}': analysis_role={ar!r} "
                     f"(expected 'legacy_reference_only')")
        # status contains 'legacy' or 'deprecated'
        status = str(entry.get("status", "")).lower()
        if "legacy" not in status and "deprecated" not in status:
            vr.error(f"registry provider '{pid}': status={entry.get('status')!r} "
                     f"does not contain 'legacy' or 'deprecated'")

    if scgpt_providers and vr.ok:
        vr.note(f"All {len(scgpt_providers)} scGPT registry provider(s) correctly "
                f"marked as legacy")

    # Milestone providers must NOT be marked legacy/deprecated
    for pid, entry in milestone_providers.items():
        if _truthy(entry.get("deprecated", False)):
            vr.error(f"registry milestone provider '{pid}': deprecated=True "
                     f"(milestone providers must not be deprecated)")
        if _truthy(entry.get("legacy", False)):
            vr.error(f"registry milestone provider '{pid}': legacy=True "
                     f"(milestone providers must not be legacy)")
        # label_mode must NOT be legacy
        lm = entry.get("label_mode", "")
        if "legacy" in lm.lower():
            vr.error(f"registry milestone provider '{pid}': label_mode={lm!r} "
                     f"contains 'legacy' (must not)")

    if milestone_providers:
        vr.note(f"All {len(milestone_providers)} milestone provider(s) correctly "
                f"not marked legacy")


# ---------------------------------------------------------------------------
# CSV row checks
# ---------------------------------------------------------------------------

def _check_csv_no_scgpt_primary(
    vr: ValidationResult,
    rows: List[Dict],
    csv_name: str,
    pid_col: str = "provider_id",
    lm_col: str = "label_mode",
    csk_col: str = "cell_state_key",
    primary_col: str = "is_primary",
    legacy_col: str = "is_legacy",
) -> None:
    """Apply all scGPT-no-primary policy checks to a CSV's rows."""
    if not rows:
        vr.note(f"{csv_name}: no data rows — skip row-level checks")
        return

    cols = set(rows[0].keys())
    has_pid = pid_col in cols
    has_lm = lm_col in cols
    has_csk = csk_col in cols
    has_primary = primary_col in cols
    has_legacy = legacy_col in cols

    n_primary = 0
    n_scgpt_primary = 0

    for i, row in enumerate(rows, start=2):
        is_primary = _truthy(row.get(primary_col, "")) if has_primary else False
        is_legacy = _truthy(row.get(legacy_col, "")) if has_legacy else False
        pid = row.get(pid_col, "") if has_pid else ""
        lm = row.get(lm_col, "") if has_lm else ""
        csk = row.get(csk_col, "") if has_csk else ""

        if is_primary:
            n_primary += 1

        # scGPT provider must never be primary
        if has_pid and "scgpt" in pid.lower() and is_primary:
            vr.error(
                f"{csv_name} row {i}: provider_id={pid!r} contains 'scgpt' "
                f"but is_primary=True"
            )
            n_scgpt_primary += 1

        # legacy_scgpt_pseudostate label_mode must never be primary
        if has_lm and lm == "legacy_scgpt_pseudostate" and is_primary:
            vr.error(
                f"{csv_name} row {i}: label_mode='legacy_scgpt_pseudostate' "
                f"but is_primary=True"
            )

        # scgpt_pseudostate_provisional cell_state_key must never be primary
        if has_csk and csk == "scgpt_pseudostate_provisional" and is_primary:
            vr.error(
                f"{csv_name} row {i}: cell_state_key='scgpt_pseudostate_provisional' "
                f"but is_primary=True"
            )

        # Legacy rows must never be primary
        if is_legacy and is_primary:
            vr.error(
                f"{csv_name} row {i}: is_legacy=True AND is_primary=True "
                f"(pid={pid!r}, lm={lm!r})"
            )

        # Primary rows must satisfy the formal primary policy. Current primary
        # label spaces are the original consensus labels and the official silver
        # labels built by the marker/foundation-model transition workflow.
        if is_primary:
            if has_lm and lm not in PRIMARY_LABEL_MODES:
                vr.error(
                    f"{csv_name} row {i}: is_primary=True but label_mode={lm!r} "
                    f"(only {sorted(PRIMARY_LABEL_MODES)!r} may be primary)"
                )
            if has_pid and pid and not any(marker in pid for marker in PRIMARY_PROVIDER_MARKERS):
                vr.error(
                    f"{csv_name} row {i}: is_primary=True but provider_id={pid!r} "
                    f"does not contain one of {sorted(PRIMARY_PROVIDER_MARKERS)!r}"
                )
            if is_legacy:
                vr.error(
                    f"{csv_name} row {i}: is_primary=True but is_legacy=True"
                )

    vr.note(f"{csv_name}: {len(rows)} row(s), {n_primary} primary row(s), "
            f"{n_scgpt_primary} scGPT-primary violations")


# ---------------------------------------------------------------------------
# Main validator
# ---------------------------------------------------------------------------

def validate(registry_path: str, reports_dir: str) -> ValidationResult:
    vr = ValidationResult()

    # --- Registry ---
    reg_path = Path(registry_path)
    if not reg_path.exists():
        vr.error(f"Registry not found: {reg_path}")
    else:
        try:
            providers = _load_registry(reg_path)
            _check_registry(vr, providers)
        except Exception as exc:
            vr.error(f"Cannot load registry: {exc}")

    # --- CSVs ---
    reports = Path(reports_dir)

    core_rows = _load_csv(reports / "core_summary.csv")
    if core_rows is None:
        vr.warn("core_summary.csv not found — run summarize_milestone_results first")
    else:
        _check_csv_no_scgpt_primary(vr, core_rows, "core_summary.csv",
                                    csk_col="cell_state_key")

    lineage_rows = _load_csv(reports / "lineage_summary.csv")
    if lineage_rows is None:
        vr.note("lineage_summary.csv not found (no milestone lineage results yet)")
    elif lineage_rows:
        _check_csv_no_scgpt_primary(vr, lineage_rows, "lineage_summary.csv",
                                    csk_col="cell_state_key")
    else:
        vr.note("lineage_summary.csv: 0 data rows")

    emb_rows = _load_csv(reports / "embedding_summary.csv")
    if emb_rows is None:
        vr.note("embedding_summary.csv not found (no milestone embedding results yet)")
    elif emb_rows:
        # embedding uses 'state_key' not 'cell_state_key'
        _check_csv_no_scgpt_primary(vr, emb_rows, "embedding_summary.csv",
                                    csk_col="state_key")
    else:
        vr.note("embedding_summary.csv: 0 data rows")

    return vr


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--registry",
        default=None,
        help="Path to registry.yaml. Defaults to <project_root>/benchmark/ground_truth/registry.yaml.",
    )
    parser.add_argument(
        "--reports-dir",
        default=None,
        help="Directory with summary CSVs. Defaults to <project_root>/benchmark/reports.",
    )
    args = parser.parse_args(argv)

    root = _project_root()
    registry_path = args.registry or str(root / "benchmark" / "ground_truth" / "registry.yaml")
    reports_dir = args.reports_dir or str(root / "benchmark" / "reports")

    print(f"\n[validate_legacy_scgpt] Registry : {registry_path}")
    print(f"[validate_legacy_scgpt] Reports  : {reports_dir}")
    print(f"{'=' * 60}")

    try:
        vr = validate(registry_path=registry_path, reports_dir=reports_dir)
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
        print(f"[validate_legacy_scgpt] PASSED ({len(vr.info)} checks, "
              f"{len(vr.warnings)} warnings, 0 errors)\n")
        return 0
    else:
        print(
            f"[validate_legacy_scgpt] FAILED: {len(vr.errors)} error(s), "
            f"{len(vr.warnings)} warning(s)\n",
            file=sys.stderr,
        )
        return 1


if __name__ == "__main__":
    sys.exit(main())
