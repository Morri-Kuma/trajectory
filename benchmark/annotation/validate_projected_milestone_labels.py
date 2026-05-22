"""validate_projected_milestone_labels.py
==========================================
Validate a directory produced by annotate_projected_cells.py.

Checks
------
  1.  projected_milestone_labels.csv exists and has required columns.
  2.  n_rows > 0.
  3.  final_milestone_label_expanded is not all ambiguous.
  4.  final_milestone_confidence exists and is numeric.
  5.  timepoint column exists and has at least one unique value.
  6.  If projected_expression.npy exists, row count matches CSV length.
  7.  projected_embedding.npy, if present, is NOT checked against CSV row count for
      two_provider_projection runs (it is method-native latent space, not an annotation input).
      For other modes it is still checked if present.
  7b. For two_provider_projection: projected_hvg_embedding.npy must exist and its
      row count must match CSV length.
  8.  If projected_cluster_labels.csv exists, row count matches CSV length.
  9.  projected_milestone_annotation_metadata.json exists and parses.
  10. Metadata smoke_test and formal_benchmark flags are explicit booleans.
Usage
-----
  python -m benchmark.annotation.validate_projected_milestone_labels \\
      benchmark/results/scnode/gse230659_marker_fm_silver_A_hvg2000_formal [--verbose]

Exit codes: 0 = pass, 1 = fail.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple


REQUIRED_COLUMNS = [
    "projected_cell_id",
    "timepoint",
    "final_milestone_label_coarse",
    "final_milestone_label_expanded",
    "final_milestone_confidence",
    "label_mode",
    "provider_id",
]


def _check_dir(out_dir: Path, verbose: bool) -> Tuple[List[str], List[str]]:
    errors: List[str] = []
    warnings: List[str] = []

    def note(msg: str) -> None:
        if verbose:
            print(f"    {msg}")

    # 1. projected_milestone_labels.csv
    csv_path = out_dir / "projected_milestone_labels.csv"
    if not csv_path.exists():
        errors.append("projected_milestone_labels.csv missing")
        return errors, warnings

    try:
        with open(csv_path, newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            fieldnames = reader.fieldnames or []
        note(f"projected_milestone_labels.csv: {len(rows)} rows, "
             f"{len(fieldnames)} cols ({csv_path.stat().st_size} bytes)")
    except Exception as exc:
        errors.append(f"projected_milestone_labels.csv read error: {exc}")
        return errors, warnings

    # Required columns
    missing_cols = [c for c in REQUIRED_COLUMNS if c not in fieldnames]
    if missing_cols:
        errors.append(f"Missing required columns: {missing_cols}")
    else:
        note(f"All required columns present: {REQUIRED_COLUMNS}")

    # 2. n_rows > 0
    if len(rows) == 0:
        errors.append("projected_milestone_labels.csv is empty (0 rows)")
        return errors, warnings
    note(f"n_projected_cells: {len(rows)}")

    # 3. final_milestone_label_expanded not all ambiguous
    if "final_milestone_label_expanded" in fieldnames:
        labels = [r.get("final_milestone_label_expanded", "") for r in rows]
        n_ambig = sum(1 for l in labels if l == "ambiguous")
        if n_ambig == len(rows):
            errors.append("All final_milestone_label_expanded values are 'ambiguous'")
        else:
            from collections import Counter
            counts = Counter(labels)
            note(f"final_milestone_label_expanded distribution: {dict(counts)}")
            if n_ambig > 0:
                warnings.append(f"{n_ambig}/{len(rows)} cells have label 'ambiguous'")

    # 4. final_milestone_confidence is numeric
    if "final_milestone_confidence" in fieldnames:
        non_numeric = 0
        for r in rows[:100]:  # spot-check
            v = r.get("final_milestone_confidence", "")
            try:
                float(v)
            except (ValueError, TypeError):
                non_numeric += 1
        if non_numeric > 0:
            errors.append(f"final_milestone_confidence has {non_numeric} non-numeric values (spot-checked 100)")
        else:
            note("final_milestone_confidence: numeric OK")

    # 5. timepoint has at least one unique value
    if "timepoint" in fieldnames:
        tps = set(r.get("timepoint", "") for r in rows)
        if len(tps) == 0 or (len(tps) == 1 and "" in tps):
            errors.append("timepoint column is empty or all blank")
        else:
            note(f"timepoints ({len(tps)} unique): {sorted(tps)[:10]}")

    n_csv = len(rows)

    # 6. projected_expression.npy row count check
    expr_path = out_dir / "projected_expression.npy"
    if expr_path.exists():
        try:
            import numpy as np
            expr = np.load(expr_path, allow_pickle=False)
            if expr.ndim == 3:
                n_expr = expr.shape[0] * expr.shape[1]
            else:
                n_expr = expr.shape[0]
            if n_expr != n_csv:
                errors.append(
                    f"projected_expression.npy has {n_expr} cells "
                    f"but projected_milestone_labels.csv has {n_csv} rows"
                )
            else:
                note(f"projected_expression.npy n_cells={n_expr} matches CSV")
        except Exception as exc:
            warnings.append(f"Could not load projected_expression.npy: {exc}")
    else:
        note("projected_expression.npy not present (OK for smoke mode)")

    # 7. projected_embedding.npy - method-native latent space only for two_provider_projection.
    #    Row-count check is skipped for formal two_provider_projection runs; the file may exist
    #    with a different row count (e.g. full projection vs. subset annotated).
    #    For other modes it is still checked if present.
    emb_path = out_dir / "projected_embedding.npy"
    _is_two_provider_for_emb = False
    try:
        import json as _json
        _mp_for_emb = out_dir / "projected_milestone_annotation_metadata.json"
        if _mp_for_emb.exists():
            with open(_mp_for_emb, encoding="utf-8") as _fe:
                _meta_emb = _json.load(_fe)
            _is_two_provider_for_emb = _meta_emb.get("annotation_policy") == "two_provider_projection"
    except Exception:
        pass

    if emb_path.exists():
        if _is_two_provider_for_emb:
            note("projected_embedding.npy present (method-native latent; row count not checked for two_provider_projection)")
        else:
            try:
                import numpy as np
                emb = np.load(emb_path, allow_pickle=False)
                n_emb = emb.shape[0]
                if n_emb != n_csv:
                    errors.append(
                        f"projected_embedding.npy has {n_emb} rows "
                        f"but projected_milestone_labels.csv has {n_csv} rows"
                    )
                else:
                    note(f"projected_embedding.npy n_cells={n_emb} matches CSV")
            except Exception as exc:
                warnings.append(f"Could not load projected_embedding.npy: {exc}")
    else:
        note("projected_embedding.npy not present (OK for smoke mode)")

    # 7b. projected_hvg_embedding.npy - required for two_provider_projection.
    hvg_emb_path = out_dir / "projected_hvg_embedding.npy"
    if _is_two_provider_for_emb:
        if not hvg_emb_path.exists():
            errors.append(
                "projected_hvg_embedding.npy missing for two_provider_projection run.  "
                "Re-run annotate_projected_cells.py to regenerate the HVG2000 PCA embedding sidecar."
            )
        else:
            try:
                import numpy as np
                hvg_emb = np.load(hvg_emb_path, allow_pickle=False)
                n_hvg_emb = hvg_emb.shape[0]
                if n_hvg_emb != n_csv:
                    errors.append(
                        f"projected_hvg_embedding.npy has {n_hvg_emb} rows "
                        f"but projected_milestone_labels.csv has {n_csv} rows"
                    )
                else:
                    note(f"projected_hvg_embedding.npy shape={hvg_emb.shape} matches CSV rows")
            except Exception as exc:
                warnings.append(f"Could not load projected_hvg_embedding.npy: {exc}")
    elif hvg_emb_path.exists():
        note("projected_hvg_embedding.npy present (not checked for non-two_provider_projection runs)")

    # 8. projected_cluster_labels.csv row count check
    clust_path = out_dir / "projected_cluster_labels.csv"
    if clust_path.exists():
        try:
            with open(clust_path, newline="", encoding="utf-8") as f:
                n_clust = sum(1 for _ in csv.DictReader(f))
            if n_clust != n_csv:
                errors.append(
                    f"projected_cluster_labels.csv has {n_clust} rows "
                    f"but projected_milestone_labels.csv has {n_csv} rows"
                )
            else:
                note(f"projected_cluster_labels.csv n_cells={n_clust} matches CSV")
        except Exception as exc:
            warnings.append(f"Could not read projected_cluster_labels.csv: {exc}")
    else:
        note("projected_cluster_labels.csv not present (OK for observed-copy-smoke)")

    # 9. Metadata JSON
    meta_path = out_dir / "projected_milestone_annotation_metadata.json"
    if not meta_path.exists():
        errors.append("projected_milestone_annotation_metadata.json missing")
    else:
        try:
            with open(meta_path, encoding="utf-8") as f:
                meta: Dict[str, Any] = json.load(f)
            note(f"metadata JSON parsed ({meta_path.stat().st_size} bytes)")
        except Exception as exc:
            errors.append(f"metadata JSON parse error: {exc}")
            meta = {}

        # 10. smoke_test and formal_benchmark explicit
        if "smoke_test" not in meta:
            errors.append("metadata missing 'smoke_test' key")
        elif not isinstance(meta["smoke_test"], bool):
            errors.append(f"metadata.smoke_test is not bool: {meta['smoke_test']!r}")
        else:
            note(f"smoke_test: {meta['smoke_test']}")

        if "formal_benchmark" not in meta:
            errors.append("metadata missing 'formal_benchmark' key")
        elif not isinstance(meta["formal_benchmark"], bool):
            errors.append(f"metadata.formal_benchmark is not bool: {meta['formal_benchmark']!r}")
        else:
            note(f"formal_benchmark: {meta['formal_benchmark']}")

        if meta.get("smoke_test") is True and meta.get("formal_benchmark") is True:
            errors.append("metadata cannot have both smoke_test == True and formal_benchmark == True")

        # annotation_policy recorded
        policy = meta.get("annotation_policy", "")
        if not policy:
            warnings.append("metadata.annotation_policy is missing or empty")
        else:
            note(f"annotation_policy: {policy!r}")

        # n_projected_cells matches CSV
        n_meta = meta.get("n_projected_cells")
        if n_meta is not None and n_meta != n_csv:
            errors.append(
                f"metadata.n_projected_cells={n_meta} != CSV rows={n_csv}"
            )
        else:
            note(f"n_projected_cells in metadata: {n_meta}")

    return errors, warnings


def main(argv: List[str] | None = None) -> int:
    import argparse
    parser = argparse.ArgumentParser(
        description="Validate projected_milestone_labels output directory."
    )
    parser.add_argument(
        "root",
        nargs="?",
        default="benchmark/results/scnode/gse230659_marker_fm_silver_A_hvg2000_formal",
        help="Output directory to validate.",
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
        print(f"ERROR: directory does not exist: {root}")
        return 1

    # Single dir or multi-dir mode
    has_csv = (root / "projected_milestone_labels.csv").exists()
    if has_csv:
        dirs_to_check = [root]
    else:
        dirs_to_check = sorted(d for d in root.iterdir() if d.is_dir())
        if not dirs_to_check:
            print(f"No projected_milestone_labels.csv or subdirs found in {root}")
            return 0

    print(f"Validating {len(dirs_to_check)} projected-label dir(s):")
    print()
    total_errors = 0
    for d in dirs_to_check:
        errors, warnings = _check_dir(d, args.verbose)
        label = "PASS" if not errors else "FAIL"
        print(f"[{label}] {d.name}")
        for w in warnings:
            print(f"       WARN: {w}")
        for e in errors:
            print(f"       ERROR: {e}")
        total_errors += len(errors)

    print()
    if total_errors == 0:
        print(f"All {len(dirs_to_check)} dir(s) passed validation.")
        return 0
    else:
        print(f"{total_errors} error(s) found across {len(dirs_to_check)} dir(s).")
        return 1


if __name__ == "__main__":
    sys.exit(main())
