"""validate_output_naming_compatibility.py
==========================================
Step 10 compatibility validator.

Scans one or more result roots and checks that milestone output directories
meet the Step 10 naming policy:

  Embedding outputs
  -----------------
  - mode-specific embedding_metrics_<label_mode>.json files must exist for the
    label modes declared by provider_agreement_metrics.json, or at least one
    mode-specific embedding_metrics_*.json must exist when no agreement file is
    available.
  - Each mode-specific JSON must contain label_mode, provider_id, state_key.
  - provider_agreement_metrics.json must exist alongside mode-specific JSONs.

  Lineage outputs
  ---------------
  - lineage_metrics.json must contain label_mode, provider_id, cell_state_key
    at the top level (promoted by _enrich_lineage_metrics_compat).
  - Generic compatibility files (lineage_metrics.json, state_transition_matrix.csv,
    lineage_graph_edges.csv) are allowed to exist and carry compatibility_note.

  Projected annotation metadata
  ------------------------------
  - projected_milestone_annotation_metadata.json must contain
    label_modes_available, state_keys_available, provider_ids_available, or the
    older equivalent fields label_modes, state_keys, provider_ids.

  General
  -------
  - No two mode-specific output files inside the same directory should have
    identical label_mode values (guards against overwrite collisions).
  - Legacy generic metrics files in milestone-identified directories must
    contain compatibility_note (or at minimum label_mode/provider_id).

Exit codes
----------
  0  all checks passed
  1  one or more required milestone compatibility checks failed

Usage
-----
  python -m benchmark.evaluation.validate_output_naming_compatibility benchmark/results/smoke
  python -m benchmark.evaluation.validate_output_naming_compatibility \\
      benchmark/results/scnode/scenario_A_scgpt_v1_hvg2000_formal
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

EMBEDDING_MODE_PREFIX = "embedding_metrics_"
EMBEDDING_MODE_SUFFIX = ".json"

EMBEDDING_REQUIRED_FIELDS = {"label_mode", "provider_id", "state_key"}

LINEAGE_METRICS_REQUIRED_FIELDS = {"label_mode", "provider_id", "cell_state_key"}

ANNOTATION_META_REQUIRED_FIELDS = {
    "label_modes_available",
    "state_keys_available",
    "provider_ids_available",
}

ANNOTATION_META_FIELD_ALIASES = {
    "label_modes_available": ("label_modes_available", "label_modes", "label_mode"),
    "state_keys_available": ("state_keys_available", "state_keys", "state_key"),
    "provider_ids_available": ("provider_ids_available", "provider_ids", "provider_id"),
}


def _load_json(path: Path) -> Optional[Dict[str, Any]]:
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except Exception as exc:
        return None


def _embedding_mode_files(directory: Path) -> List[Path]:
    """Return mode-specific embedding metric files in a milestone eval dir."""
    return sorted(
        p for p in directory.glob(f"{EMBEDDING_MODE_PREFIX}*{EMBEDDING_MODE_SUFFIX}")
        if p.name != "embedding_metrics.json"
    )


def _mode_filename(label_mode: str) -> str:
    return f"{EMBEDDING_MODE_PREFIX}{label_mode}{EMBEDDING_MODE_SUFFIX}"


def _is_milestone_dir(directory: Path) -> bool:
    """Heuristic: is this directory a milestone output directory?"""
    # Positive signals
    signals = [
        "embedding_milestone_eval" in str(directory),
        "milestone_lineage" in str(directory),
        "projected_milestone_labels" in str(directory),
        "lineage_consensus" in str(directory),
        # Contains any mode-specific embedding JSON
        bool(_embedding_mode_files(directory)),
        # Contains projected annotation metadata
        (directory / "projected_milestone_annotation_metadata.json").exists(),
        # Contains smoke_metadata.json with milestone-related provider
        _has_milestone_smoke_metadata(directory),
    ]
    return any(signals)


def _has_milestone_smoke_metadata(directory: Path) -> bool:
    meta = _load_json(directory / "smoke_metadata.json")
    if meta is None:
        meta = _load_json(directory / "run_metadata.json")
    if meta is None:
        return False
    provider_id = meta.get("provider_id", "")
    return "milestone" in str(provider_id).lower()


# ---------------------------------------------------------------------------
# Individual checks
# ---------------------------------------------------------------------------

def check_embedding_milestone_dir(
    directory: Path,
    results: List[str],
    failures: List[str],
) -> None:
    """Check a directory that appears to contain embedding milestone outputs."""
    agree_path = directory / "provider_agreement_metrics.json"
    agree = _load_json(agree_path) if agree_path.exists() else None
    expected_files = []
    if agree is not None:
        expected_files = [
            _mode_filename(str(mode))
            for mode in agree.get("label_modes", [])
            if str(mode).strip()
        ]

    mode_files_by_name = {p.name: p for p in _embedding_mode_files(directory)}
    files_to_check: List[Path] = []

    if expected_files:
        for fname in expected_files:
            fpath = directory / fname
            if not fpath.exists():
                failures.append(
                    f"[embedding] MISSING mode-specific file: {fpath.relative_to(fpath.parents[3]) if len(fpath.parts) > 3 else fpath}"
                )
                continue
            files_to_check.append(fpath)
    elif mode_files_by_name:
        files_to_check = list(mode_files_by_name.values())
    else:
        failures.append(f"[embedding] MISSING mode-specific embedding_metrics_*.json in {directory}")

    for fpath in files_to_check:
        fname = fpath.name
        data = _load_json(fpath)
        if data is None:
            failures.append(f"[embedding] INVALID JSON: {fpath.name} in {directory}")
            continue

        missing = EMBEDDING_REQUIRED_FIELDS - set(data.keys())
        if missing:
            failures.append(
                f"[embedding] {fpath.name}: missing fields {sorted(missing)} "
                f"in {directory}"
            )
        else:
            results.append(f"[embedding] OK  {fpath.name} has required fields")

    if not agree_path.exists():
        failures.append(f"[embedding] MISSING provider_agreement_metrics.json in {directory}")
    else:
        if agree is None:
            failures.append(f"[embedding] INVALID JSON: provider_agreement_metrics.json in {directory}")
        else:
            results.append(f"[embedding] OK  provider_agreement_metrics.json exists")

    # Check no label_mode collision across all mode-specific files.
    seen_modes: Dict[str, str] = {}
    for fpath in mode_files_by_name.values():
        fname = fpath.name
        data = _load_json(fpath)
        if data is None:
            continue
        mode = data.get("label_mode", "")
        if mode and mode in seen_modes:
            failures.append(
                f"[embedding] COLLISION: label_mode={mode!r} appears in both "
                f"{seen_modes[mode]} and {fname} in {directory}"
            )
        elif mode:
            seen_modes[mode] = fname


def check_lineage_metrics(
    metrics_path: Path,
    results: List[str],
    failures: List[str],
) -> None:
    """Check a lineage_metrics.json for Step 10 required fields."""
    data = _load_json(metrics_path)
    if data is None:
        failures.append(f"[lineage] INVALID JSON: {metrics_path}")
        return

    # Ground truth sub-object provides fallback for older files
    gt = data.get("ground_truth") or {}

    missing_top = []
    for field in sorted(LINEAGE_METRICS_REQUIRED_FIELDS):
        has_top = field in data
        has_gt = field in gt or (field == "cell_state_key" and "state_key" in gt)
        if not has_top and not has_gt:
            missing_top.append(field)

    if missing_top:
        failures.append(
            f"[lineage] {metrics_path.name}: missing fields {missing_top} "
            f"(top-level or ground_truth) in {metrics_path.parent}"
        )
    else:
        lm = data.get("label_mode") or gt.get("label_mode", "?")
        pid = data.get("provider_id") or gt.get("provider_id", "?")
        csk = data.get("cell_state_key") or gt.get("cell_state_key") or gt.get("state_key", "?")
        results.append(
            f"[lineage] OK  {metrics_path.parent.name}/{metrics_path.name} "
            f"label_mode={lm!r} provider_id={pid!r} cell_state_key={csk!r}"
        )


def check_projected_annotation_metadata(
    meta_path: Path,
    results: List[str],
    failures: List[str],
) -> None:
    """Check projected_milestone_annotation_metadata.json for Step 10 fields."""
    data = _load_json(meta_path)
    if data is None:
        failures.append(f"[annotation] INVALID JSON: {meta_path}")
        return

    missing = [
        canonical
        for canonical, aliases in ANNOTATION_META_FIELD_ALIASES.items()
        if not any(alias in data for alias in aliases)
    ]
    if missing:
        failures.append(
            f"[annotation] {meta_path.parent.name}/projected_milestone_annotation_metadata.json: "
            f"missing fields {sorted(missing)}"
        )
    else:
        lma = (
            data.get("label_modes_available")
            or data.get("label_modes")
            or ([data["label_mode"]] if data.get("label_mode") else [])
        )
        results.append(
            f"[annotation] OK  {meta_path.parent.name}/projected_milestone_annotation_metadata.json "
            f"label_modes_available={lma}"
        )


def check_generic_legacy_in_milestone_dir(
    directory: Path,
    results: List[str],
    failures: List[str],
) -> None:
    """If legacy generic metrics JSONs exist in a milestone dir, they must have compat fields."""
    for fname in ["lineage_metrics.json", "embedding_metrics.json"]:
        fpath = directory / fname
        if not fpath.exists():
            continue
        data = _load_json(fpath)
        if data is None:
            failures.append(f"[legacy] INVALID JSON: {fpath}")
            continue

        gt = data.get("ground_truth") or {}
        has_compat = (
            "compatibility_note" in data
            or "label_mode" in data
            or "label_mode" in gt
        )
        if not has_compat:
            failures.append(
                f"[legacy] {fname} in milestone dir {directory.name} "
                f"lacks compatibility_note/label_mode — cannot identify label mode"
            )
        else:
            results.append(f"[legacy] OK  {directory.name}/{fname} has compat metadata")


# ---------------------------------------------------------------------------
# Directory scanner
# ---------------------------------------------------------------------------

def scan_directory(root: Path, results: List[str], failures: List[str]) -> None:
    """Recursively scan root for milestone output directories and validate them."""
    checked_embedding_dirs: set = set()
    checked_lineage_paths: set = set()
    checked_annotation_paths: set = set()

    for path in sorted(root.rglob("*")):
        if not path.is_dir():
            continue
        if not _is_milestone_dir(path):
            continue

        # ── Embedding milestone checks ──
        if _embedding_mode_files(path) or (path / "provider_agreement_metrics.json").exists():
            if path not in checked_embedding_dirs:
                check_embedding_milestone_dir(path, results, failures)
                checked_embedding_dirs.add(path)

        # ── Lineage metrics checks ──
        for lm_path in path.glob("lineage_metrics*.json"):
            if lm_path not in checked_lineage_paths:
                check_lineage_metrics(lm_path, results, failures)
                checked_lineage_paths.add(lm_path)

        # ── Projected annotation metadata ──
        anno_path = path / "projected_milestone_annotation_metadata.json"
        if anno_path.exists() and anno_path not in checked_annotation_paths:
            check_projected_annotation_metadata(anno_path, results, failures)
            checked_annotation_paths.add(anno_path)

        # ── Legacy generic files in milestone dirs ──
        check_generic_legacy_in_milestone_dir(path, results, failures)

    # Also check root itself if it IS a milestone dir
    if _is_milestone_dir(root):
        if root not in checked_embedding_dirs:
            if _embedding_mode_files(root) or (root / "provider_agreement_metrics.json").exists():
                check_embedding_milestone_dir(root, results, failures)
        for lm_path in root.glob("lineage_metrics*.json"):
            if lm_path not in checked_lineage_paths:
                check_lineage_metrics(lm_path, results, failures)
        anno_path = root / "projected_milestone_annotation_metadata.json"
        if anno_path.exists() and anno_path not in checked_annotation_paths:
            check_projected_annotation_metadata(anno_path, results, failures)
        check_generic_legacy_in_milestone_dir(root, results, failures)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    import argparse
    parser = argparse.ArgumentParser(
        description="Step 10 output naming compatibility validator."
    )
    parser.add_argument(
        "roots",
        nargs="+",
        help="Result root directories to scan.",
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Print passing checks as well as failures.",
    )
    args = parser.parse_args(argv)

    all_results: List[str] = []
    all_failures: List[str] = []

    for root_str in args.roots:
        root = Path(root_str)
        if not root.is_absolute():
            root = Path.cwd() / root
        if not root.exists():
            print(f"[validate_naming] WARNING: root does not exist: {root}")
            continue
        print(f"\n[validate_naming] Scanning: {root}")
        scan_directory(root, all_results, all_failures)

    print()
    if args.verbose:
        for msg in all_results:
            print(f"  PASS  {msg}")

    if all_failures:
        print(f"  FAILURES ({len(all_failures)}):")
        for msg in all_failures:
            print(f"  FAIL  {msg}")
        print(f"\n[validate_naming] Result: {len(all_failures)} failure(s), "
              f"{len(all_results)} pass(es). Exit 1.")
        return 1

    print(f"[validate_naming] All checks passed ({len(all_results)} checks). Exit 0.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
