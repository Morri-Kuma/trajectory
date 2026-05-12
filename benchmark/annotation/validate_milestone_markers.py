"""validate_milestone_markers.py
==================================
Lightweight CLI validator for benchmark/annotation/milestone_markers.yaml.

Checks performed
----------------
  1.  YAML parses without error.
  2.  Required dataset IDs are present: GSE178325, GSE230659.
  3.  Each dataset contains required top-level keys:
        primary_milestones, marker_sets, reference_graph_edges.
  4.  Every primary milestone has a corresponding entry in marker_sets.
  5.  Every node appearing in reference_graph_edges belongs to
        primary_milestones.
  6.  Every marker set entry contains: genes, role, notes.
  7.  Every gene symbol is non-empty and fully uppercase.
  8.  If --gene-list is provided: compute overlap counts and fractions per
        marker set; warn when overlap < --min-overlap.  Low overlap is
        a warning, not a structural error.
  9.  Write a JSON report to --output-json if provided.

Return codes
------------
  0   structural validation passed (gene-overlap warnings do not affect exit code)
  1   YAML or structural validation failed

Usage
-----
    python benchmark/annotation/validate_milestone_markers.py --dry-run
    python benchmark/annotation/validate_milestone_markers.py \\
        --dataset-id GSE178325 --dry-run
    python benchmark/annotation/validate_milestone_markers.py \\
        --gene-list path/to/genes.txt \\
        --output-json path/to/report.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

# Default path relative to repo root; script resolves against CWD.
_DEFAULT_MARKERS_YAML = "benchmark/annotation/milestone_markers.yaml"
_REQUIRED_DATASET_IDS = ("GSE178325", "GSE230659")
_REQUIRED_DATASET_KEYS = (
    "primary_milestones",
    "marker_sets",
    "reference_graph_edges",
)
_REQUIRED_MARKER_KEYS = ("genes", "role", "notes")


# ---------------------------------------------------------------------------
# YAML loader (lazy import so standard-library-only usage stays possible)
# ---------------------------------------------------------------------------

def _load_yaml(path: Path) -> Any:
    try:
        import yaml  # PyYAML
    except ImportError as exc:
        raise SystemExit(
            "ERROR: PyYAML is required for YAML parsing.  "
            "Install it with: pip install pyyaml\n"
            f"Original error: {exc}"
        ) from exc
    with path.open(encoding="ascii") as fh:
        return yaml.safe_load(fh)


# ---------------------------------------------------------------------------
# Structural validation helpers
# ---------------------------------------------------------------------------

def _check_dataset(
    dataset_id: str,
    dataset: Any,
    errors: list,
    warnings: list,
    gene_universe: set,
    min_overlap: int,
) -> dict:
    """Validate a single dataset block; populate errors/warnings in-place.

    Returns a per-dataset marker-set summary dict.
    """
    summary: dict = {}

    # Check top-level required keys.
    if not isinstance(dataset, dict):
        errors.append(f"{dataset_id}: dataset block is not a mapping")
        return summary

    for key in _REQUIRED_DATASET_KEYS:
        if key not in dataset:
            errors.append(f"{dataset_id}: missing required key '{key}'")

    # Bail early if critical keys are absent so we avoid cascading KeyErrors.
    if any(k not in dataset for k in ("primary_milestones", "marker_sets",
                                       "reference_graph_edges")):
        return summary

    primary = list(dataset.get("primary_milestones") or [])
    marker_sets = dataset.get("marker_sets") or {}
    edges = dataset.get("reference_graph_edges") or []

    # Check 4: every primary milestone has a marker set.
    for ms in primary:
        if ms not in marker_sets:
            errors.append(
                f"{dataset_id}: primary milestone '{ms}' has no entry in marker_sets"
            )

    # Collect all graph nodes.
    graph_nodes: set = set()
    for edge in edges:
        if not (isinstance(edge, (list, tuple)) and len(edge) == 2):
            errors.append(
                f"{dataset_id}: reference_graph_edges entry is not a 2-element "
                f"list: {edge!r}"
            )
            continue
        graph_nodes.update(edge)

    # Check 5: every graph node is in primary_milestones.
    primary_set = set(primary)
    for node in sorted(graph_nodes):
        if node not in primary_set:
            errors.append(
                f"{dataset_id}: graph node '{node}' appears in "
                f"reference_graph_edges but is not in primary_milestones"
            )

    # Checks 6, 7, 8 per marker set.
    if not isinstance(marker_sets, dict):
        errors.append(f"{dataset_id}: marker_sets is not a mapping")
        return summary

    for ms_name, ms_body in marker_sets.items():
        ms_errors: list = []
        ms_summary: dict = {"milestone": ms_name}

        if not isinstance(ms_body, dict):
            errors.append(
                f"{dataset_id}.marker_sets.{ms_name}: entry is not a mapping"
            )
            continue

        # Check 7: required marker-set keys.
        for key in _REQUIRED_MARKER_KEYS:
            if key not in ms_body:
                ms_errors.append(f"missing key '{key}'")

        genes = ms_body.get("genes") or []
        ms_summary["gene_count"] = len(genes)
        ms_summary["role"] = ms_body.get("role", "")

        # Check 8: gene symbols.
        bad_genes = []
        for g in genes:
            if not isinstance(g, str) or not g.strip():
                bad_genes.append(repr(g))
            elif g != g.upper():
                bad_genes.append(g)
        if bad_genes:
            ms_errors.append(
                f"gene symbols not uppercase or empty: {bad_genes}"
            )

        if ms_errors:
            for msg in ms_errors:
                errors.append(f"{dataset_id}.marker_sets.{ms_name}: {msg}")

        # Check 9: gene-list overlap (warning only).
        if gene_universe:
            gene_set = {g for g in genes if isinstance(g, str) and g.strip()}
            overlap = gene_set & gene_universe
            overlap_count = len(overlap)
            overlap_frac = overlap_count / len(gene_set) if gene_set else 0.0
            ms_summary["overlap_count"] = overlap_count
            ms_summary["overlap_fraction"] = round(overlap_frac, 4)
            ms_summary["overlap_genes"] = sorted(overlap)
            if overlap_count < min_overlap:
                warnings.append(
                    f"{dataset_id}.marker_sets.{ms_name}: overlap with "
                    f"gene list is {overlap_count} < min_overlap={min_overlap} "
                    f"(genes matched: {sorted(overlap)})"
                )

        summary[ms_name] = ms_summary

    return summary


# ---------------------------------------------------------------------------
# Main validation entry point
# ---------------------------------------------------------------------------

def validate(
    markers_yaml: str | Path,
    dataset_ids: list[str] | None,
    gene_list_path: str | Path | None,
    min_overlap: int,
    output_json: str | Path | None,
    dry_run: bool,
) -> int:
    """Run all validation checks; return 0 on success, 1 on failure."""
    markers_yaml = Path(markers_yaml)
    errors: list = []
    warnings: list = []
    report: dict = {
        "status": "valid",
        "markers_yaml": str(markers_yaml),
        "dataset_ids_validated": [],
        "structural_errors": [],
        "warnings": [],
        "marker_set_summary": {},
        "gene_list_used": None,
        "min_overlap": min_overlap,
    }

    # --- Check 1: YAML parses. ---
    if not markers_yaml.exists():
        errors.append(f"markers YAML not found: {markers_yaml}")
        _finish(report, errors, warnings, output_json, dry_run)
        return 1

    try:
        data = _load_yaml(markers_yaml)
    except Exception as exc:  # noqa: BLE001
        errors.append(f"YAML parse error: {exc}")
        _finish(report, errors, warnings, output_json, dry_run)
        return 1

    if not isinstance(data, dict):
        errors.append("Top-level YAML is not a mapping")
        _finish(report, errors, warnings, output_json, dry_run)
        return 1

    # --- Check 2: required dataset IDs. ---
    target_ids = dataset_ids or list(_REQUIRED_DATASET_IDS)
    for did in _REQUIRED_DATASET_IDS:
        if did not in data:
            errors.append(f"Required dataset ID '{did}' not found in YAML")

    # --- Load optional gene universe. ---
    gene_universe: set = set()
    if gene_list_path is not None:
        gene_list_path = Path(gene_list_path)
        if not gene_list_path.exists():
            errors.append(f"--gene-list file not found: {gene_list_path}")
        else:
            raw = gene_list_path.read_text(encoding="ascii").splitlines()
            gene_universe = {line.strip() for line in raw if line.strip()}
            report["gene_list_used"] = str(gene_list_path)
            print(
                f"[validate_milestone_markers] gene list loaded: "
                f"{len(gene_universe)} genes from {gene_list_path}"
            )

    # --- Per-dataset structural validation. ---
    validated_ids = []
    all_summary: dict = {}
    for did in target_ids:
        if did not in data:
            # Already flagged in check 2; skip to avoid KeyError.
            continue
        validated_ids.append(did)
        ds_summary = _check_dataset(
            did, data[did], errors, warnings,
            gene_universe, min_overlap,
        )
        all_summary[did] = ds_summary

    report["dataset_ids_validated"] = validated_ids
    report["marker_set_summary"] = all_summary

    return _finish(report, errors, warnings, output_json, dry_run)


def _finish(
    report: dict,
    errors: list,
    warnings: list,
    output_json: str | Path | None,
    dry_run: bool,
) -> int:
    """Populate report, print summary, optionally write JSON; return exit code."""
    ok = len(errors) == 0
    report["status"] = "valid" if ok else "invalid"
    report["structural_errors"] = errors
    report["warnings"] = warnings

    # Print summary.
    tag = "[validate_milestone_markers]"
    if dry_run:
        print(f"{tag} dry-run mode: validation checks will run; no files written.")
    if errors:
        print(f"{tag} STRUCTURAL ERRORS ({len(errors)}):")
        for e in errors:
            print(f"  ERROR: {e}")
    else:
        print(f"{tag} structural validation PASSED.")
    if warnings:
        print(f"{tag} warnings ({len(warnings)}):")
        for w in warnings:
            print(f"  WARN:  {w}")
    else:
        print(f"{tag} no warnings.")

    if output_json and not dry_run:
        out = Path(output_json)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", encoding="ascii") as fh:
            json.dump(report, fh, indent=2)
        print(f"{tag} report written -> {out}")
    elif output_json and dry_run:
        print(f"{tag} dry-run: skipping JSON write to {output_json}")

    return 0 if ok else 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Validate benchmark/annotation/milestone_markers.yaml "
            "(scTimeBench framework v2, Step 2)"
        )
    )
    parser.add_argument(
        "--markers-yaml",
        default=_DEFAULT_MARKERS_YAML,
        help=(
            f"Path to milestone markers YAML "
            f"(default: {_DEFAULT_MARKERS_YAML})"
        ),
    )
    parser.add_argument(
        "--dataset-id",
        default=None,
        help=(
            "Validate a single dataset ID (e.g. GSE178325).  "
            "If omitted, all required datasets are validated."
        ),
    )
    parser.add_argument(
        "--gene-list",
        default=None,
        dest="gene_list",
        help=(
            "Optional path to a plain-text file with one gene symbol per line.  "
            "Enables overlap checks against each marker set."
        ),
    )
    parser.add_argument(
        "--min-overlap",
        type=int,
        default=2,
        dest="min_overlap",
        help=(
            "Minimum number of genes from --gene-list that must overlap with "
            "a marker set before issuing a warning (default: 2)."
        ),
    )
    parser.add_argument(
        "--output-json",
        default=None,
        dest="output_json",
        help="Optional path for the JSON validation report.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run all validation checks but skip writing any output files.",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    dataset_ids = [args.dataset_id] if args.dataset_id else None

    rc = validate(
        markers_yaml=args.markers_yaml,
        dataset_ids=dataset_ids,
        gene_list_path=args.gene_list,
        min_overlap=args.min_overlap,
        output_json=args.output_json,
        dry_run=args.dry_run,
    )
    sys.exit(rc)


if __name__ == "__main__":
    main()
