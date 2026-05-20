"""build_projected_milestone_transitions.py
============================================
Step 9 bridge: Build state_transition_matrix.csv and lineage_graph_edges.csv
from projected_milestone_labels.csv using aggregate label frequencies at
consecutive timepoints.

Transition policy: aggregate_frequency_consecutive_timepoints
--------------------------------------------------------------
For each pair of consecutive timepoints (t1, t2):
  - Get label distribution at t1 and t2.
  - For each source state s at t1, transition probabilities to target states
    equal the label frequency distribution at t2.
  - This reflects the "what proportion of cells would reach each state" at the
    next timepoint, given no paired single-cell tracking across timepoints.

This is explicitly a smoke/wiring bridge, not a biological inference method.
Record transition_policy in output metadata.

Usage
-----
  python -m benchmark.evaluation.build_projected_milestone_transitions \\
      --projected-labels <csv> \\
      --output-dir <dir> \\
      --label-mode official_silver \\
      --cell-state-key final_milestone_label_coarse \\
      --exclude-label ambiguous

  python -m benchmark.evaluation.build_projected_milestone_transitions \\
      --projected-labels <csv> \\
      --output-dir <dir> \\
      --label-mode official_silver \\
      --cell-state-key final_milestone_label_coarse \\
      --exclude-label ambiguous \\
      --all-timepoint-pairs
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


def _load_projected_labels(csv_path: Path) -> List[Dict[str, str]]:
    with open(csv_path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _build_transition_matrix(
    rows: List[Dict[str, str]],
    label_col: str,
    timepoint_col: str = "timepoint",
    exclude_labels: Optional[List[str]] = None,
    all_timepoint_pairs: bool = False,
) -> Tuple[Any, Any, Dict[str, Any]]:
    """
    Build state transition matrix from projected milestone labels.

    Returns
    -------
    stm_df : dict-of-dicts {source_state: {target_state: weight}}
    edges : list of (source, target, weight) tuples
    metadata : dict
    """
    import pandas as pd

    exclude_labels = list(exclude_labels) if exclude_labels else []

    # Extract timepoints and labels
    timepoints_raw = [r.get(timepoint_col, "") for r in rows]
    labels_raw = [r.get(label_col, "") for r in rows]

    # Parse timepoints as float if possible
    def _try_float(v: str):
        try:
            return float(v)
        except (ValueError, TypeError):
            return v

    timepoints = [_try_float(t) for t in timepoints_raw]
    unique_tps = sorted(set(timepoints))
    n_tps = len(unique_tps)

    # Filter excluded labels
    valid_mask = [lbl not in exclude_labels for lbl in labels_raw]

    # Group by timepoint
    tp_labels: Dict[Any, List[str]] = {tp: [] for tp in unique_tps}
    for tp, lbl, valid in zip(timepoints, labels_raw, valid_mask):
        if valid:
            tp_labels[tp].append(lbl)

    # All observed states across all timepoints
    all_states = sorted(set(
        lbl for lbl_list in tp_labels.values() for lbl in lbl_list
    ))

    if not all_states:
        raise ValueError(
            f"No valid states found in {label_col} after excluding {exclude_labels}. "
            "Cannot build transition matrix."
        )

    # Build transition pairs
    if all_timepoint_pairs:
        tp_pairs = [(unique_tps[i], unique_tps[j])
                    for i in range(n_tps) for j in range(i+1, n_tps)]
    else:
        tp_pairs = list(zip(unique_tps[:-1], unique_tps[1:]))

    if not tp_pairs:
        raise ValueError(
            f"Only one unique timepoint ({unique_tps}). "
            "Cannot build transition matrix — need at least 2 timepoints."
        )

    # Aggregate transitions: average over all consecutive pairs
    # T[s1, s2] = mean over consecutive pairs of freq_t2[s2]
    state_to_idx = {s: i for i, s in enumerate(all_states)}
    n_states = len(all_states)

    accum = np.zeros((n_states, n_states), dtype=np.float64)
    n_pairs_used = 0

    for t1, t2 in tp_pairs:
        labels_t1 = tp_labels.get(t1, [])
        labels_t2 = tp_labels.get(t2, [])
        if not labels_t1 or not labels_t2:
            continue

        # Frequency distribution at t2
        freq_t2 = Counter(labels_t2)
        total_t2 = sum(freq_t2.values())
        if total_t2 == 0:
            continue

        for s2, cnt in freq_t2.items():
            if s2 not in state_to_idx:
                continue
            j = state_to_idx[s2]
            # All source states at t1 get this probability for this pair
            freq_t1 = Counter(labels_t1)
            total_t1 = sum(freq_t1.values())
            for s1, cnt1 in freq_t1.items():
                if s1 not in state_to_idx:
                    continue
                i = state_to_idx[s1]
                # Weight by source frequency
                accum[i, j] += (cnt1 / total_t1) * (cnt / total_t2)

        n_pairs_used += 1

    if n_pairs_used == 0:
        raise ValueError("No valid consecutive timepoint pairs with data on both sides.")

    # Row-normalize
    row_sums = accum.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # avoid division by zero
    T = accum / row_sums

    # Convert to DataFrames / edge list
    stm_df = {
        all_states[i]: {all_states[j]: round(float(T[i, j]), 6)
                        for j in range(n_states)}
        for i in range(n_states)
    }

    # Edges: keep only pairs with weight > epsilon
    edges = [
        (all_states[i], all_states[j], round(float(T[i, j]), 6))
        for i in range(n_states)
        for j in range(n_states)
        if T[i, j] > 1e-6
    ]

    # Per-timepoint label counts for metadata
    tp_label_counts = {
        str(tp): dict(Counter(lbl_list))
        for tp, lbl_list in tp_labels.items()
    }

    meta = {
        "n_timepoints": n_tps,
        "unique_timepoints": [str(tp) for tp in unique_tps],
        "n_pairs_used": n_pairs_used,
        "n_cells_total": len(rows),
        "n_cells_valid": sum(valid_mask),
        "states": all_states,
        "n_states": n_states,
        "tp_label_counts": tp_label_counts,
    }
    return stm_df, edges, meta


def run_transition_build(
    projected_labels_csv: str,
    output_dir: str,
    label_mode: str = "official_silver",
    cell_state_key: str = "final_milestone_label_coarse",
    exclude_labels: Optional[List[str]] = None,
    dataset_id: Optional[str] = None,
    all_timepoint_pairs: bool = False,
    provider_id: Optional[str] = None,
) -> Dict[str, Any]:
    import pandas as pd

    proj_csv = Path(projected_labels_csv)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    exclude_labels = list(exclude_labels) if exclude_labels else []

    print(f"[build_projected_transitions] Loading {proj_csv}")
    rows = _load_projected_labels(proj_csv)
    print(f"[build_projected_transitions] {len(rows)} projected cells, "
          f"label_col={cell_state_key!r}, label_mode={label_mode!r}")

    if cell_state_key not in (rows[0].keys() if rows else {}):
        available = list(rows[0].keys()) if rows else []
        raise ValueError(
            f"Column {cell_state_key!r} not found in CSV. "
            f"Available: {available}"
        )

    stm_dict, edges, build_meta = _build_transition_matrix(
        rows=rows,
        label_col=cell_state_key,
        timepoint_col="timepoint",
        exclude_labels=exclude_labels,
        all_timepoint_pairs=all_timepoint_pairs,
    )

    # Write state_transition_matrix.csv
    stm_path = out_dir / "state_transition_matrix.csv"
    states = build_meta["states"]
    with open(stm_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["source_state"] + states)
        for src in states:
            row_vals = [stm_dict[src].get(tgt, 0.0) for tgt in states]
            writer.writerow([src] + row_vals)
    print(f"[build_projected_transitions] Wrote STM ({len(states)} states) -> {stm_path}")

    # Write lineage_graph_edges.csv
    edges_path = out_dir / "lineage_graph_edges.csv"
    with open(edges_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["source_state", "target_state", "weight"])
        writer.writeheader()
        for src, tgt, w in edges:
            writer.writerow({"source_state": src, "target_state": tgt, "weight": w})
    print(f"[build_projected_transitions] Wrote {len(edges)} edges -> {edges_path}")

    # Write metadata
    # Step 10 naming policy: mode_specific_outputs records the mode-specific
    # file aliases; generic_compatibility_outputs records stable-name
    # generic files. Both sets exist in the same directory — mode-specific
    # files carry the label_mode suffix; generic files are kept for stable names.
    out_dir_str = str(out_dir)
    metadata = {
        "dataset_id": dataset_id,
        "projected_labels_csv": str(proj_csv),
        "output_dir": out_dir_str,
        "label_mode": label_mode,
        "provider_id": provider_id,
        "cell_state_key": cell_state_key,
        "excluded_labels": exclude_labels,
        "exclusion_reason": (
            "label not present in marker-defined milestone graph"
            if exclude_labels else "none"
        ),
        "transition_policy": "aggregate_frequency_consecutive_timepoints",
        "transition_policy_note": (
            "For each consecutive timepoint pair (t1, t2), the transition "
            "probability from source state s1 to target state s2 equals the "
            "frequency of s2 at t2, weighted by the frequency of s1 at t1. "
            "Row-normalized. This is a smoke bridge without single-cell tracking."
        ),
        "smoke_test": True,
        "formal_benchmark": False,
        "all_timepoint_pairs": all_timepoint_pairs,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "mode_specific_outputs": {
            f"state_transition_matrix_{label_mode}.csv": str(out_dir / f"state_transition_matrix_{label_mode}.csv"),
            f"lineage_graph_edges_{label_mode}.csv": str(out_dir / f"lineage_graph_edges_{label_mode}.csv"),
            f"lineage_metrics_{label_mode}.json": str(out_dir / f"lineage_metrics_{label_mode}.json"),
        },
        "generic_compatibility_outputs": {
            "state_transition_matrix.csv": str(out_dir / "state_transition_matrix.csv"),
            "lineage_graph_edges.csv": str(out_dir / "lineage_graph_edges.csv"),
            "projected_transition_metadata.json": str(out_dir / "projected_transition_metadata.json"),
        },
        "compatibility_note": (
            "generic_compatibility_outputs are kept for stable filenames. "
            "mode_specific_outputs are the preferred files for official silver reports."
        ),
        **build_meta,
    }
    meta_path = out_dir / "projected_transition_metadata.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2, default=str)
    print(f"[build_projected_transitions] Wrote metadata -> {meta_path}")

    return {
        "state_transition_matrix": str(stm_path),
        "lineage_graph_edges": str(edges_path),
        "metadata": metadata,
    }


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build projected milestone transition matrix."
    )
    parser.add_argument("--projected-labels", required=True,
                        help="Path to projected_milestone_labels.csv.")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for STM and edge files.")
    parser.add_argument("--label-mode", default="official_silver",
                        help="Label mode (default: official_silver).")
    parser.add_argument("--cell-state-key", default="final_milestone_label_coarse",
                        help="Column name in the CSV to use as state labels.")
    parser.add_argument("--exclude-label", action="append", dest="exclude_labels",
                        default=[], metavar="LABEL",
                        help="Label to exclude (repeatable).")
    parser.add_argument("--dataset-id", default=None)
    parser.add_argument("--provider-id", default=None,
                        help="Provider ID for metadata compatibility "
                             "(e.g. gse230659_marker_fm_transition_silver_v1).")
    parser.add_argument("--all-timepoint-pairs", action="store_true",
                        help="Use all pairs, not just consecutive pairs.")
    args = parser.parse_args(argv)

    try:
        run_transition_build(
            projected_labels_csv=args.projected_labels,
            output_dir=args.output_dir,
            label_mode=args.label_mode,
            cell_state_key=args.cell_state_key,
            exclude_labels=args.exclude_labels,
            dataset_id=args.dataset_id,
            all_timepoint_pairs=args.all_timepoint_pairs,
            provider_id=args.provider_id,
        )
    except Exception as exc:
        print(f"[build_projected_transitions] ERROR: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
