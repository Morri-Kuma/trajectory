"""build_marker_seed_labels.py
==============================
Annotate observed benchmark cells with milestone marker-gene seed labels.

Stage 1 of a two-stage annotation layer: assigns primary/manual/marker-based
coarse labels and explicitly routes uncertain cells to ``ambiguous`` instead
of forcing all cells into a milestone.  Stage 1 output is designed to support
later Stage 2 trajectory-aware resolution.

Step 3 implementation: reads an input .h5ad, scores each cell against curated
milestone marker-gene signatures (mean expression of matched genes), assigns
milestone labels, and writes an annotated .h5ad plus QC tables.

Annotation columns written to adata.obs
----------------------------------------
Per-milestone scores (one column per scored milestone):
  milestone_score_<milestone_name>      float32

Legacy marker-derived columns (backward-compatible; unchanged):
  milestone_marker_label                str   top-scoring milestone label
  milestone_marker_score                float32   top mean-expression score
  milestone_marker_margin               float32   top score - second-best score

Stage 1 authoritative columns (conservative, threshold-gated):
  stage1_label                          str
      Assigned label: a primary milestone name or "ambiguous".
      hADSCs cells are assigned primarily by
      source sample / time == 0 (not by generic MSC marker score).
  stage1_confidence                     float32
      Score driving the assignment (1.0 for sample/time hADSCs).
  stage1_margin                         float32
      Score margin (top - second).
  stage1_status                         str
      "high_confidence" | "ambiguous"
  stage1_source                         str
      "sample_time" | "marker_primary" | "marker_ambiguous"

Step-3 placeholder columns (compatibility only; will be replaced by
embedding/classifier/consensus logic in later steps):
  milestone_embedding_label             str
  milestone_embedding_confidence        float32
  milestone_classifier_label            str
  milestone_classifier_confidence       float32
  consensus_milestone_label             str
  consensus_confidence                  float32

Metadata written to adata.uns["milestone_annotation_step3"]
    See _build_uns_metadata() for full field list.

Output files
------------
  <output-h5ad>                         annotated AnnData
  <output-dir>/marker_seed_labels_metadata.json
  <output-dir>/marker_seed_label_counts.csv
  <output-dir>/marker_seed_marker_overlap.csv
  <output-dir>/stage1_label_counts.csv
  <output-dir>/stage1_status_counts.csv
  <output-dir>/stage1_source_counts.csv

Usage
-----
  python benchmark/annotation/build_marker_seed_labels.py \\
      --dataset-id GSE230659 \\
      --input-h5ad benchmark/inputs/.../GSE230659_...h5ad \\
      --output-h5ad benchmark/tmp_annotation_step3/GSE230659_step3.h5ad \\
      --output-dir  benchmark/tmp_annotation_step3/gse230659 \\
      --use-var-index-if-needed \\
      --max-cells 2000

hADSC assignment policy (Stage 1)
----------------------------------
hADSCs are assigned PRIMARILY by source sample / time == 0, following the
treatment in Liuyang et al. 2023 (Cell Stem Cell,
DOI: 10.1016/j.stem.2023.02.008).  The hADSCs marker-gene set in
milestone_markers.yaml is an AUXILIARY somatic/fibroblast-like support score
only.  If obs[sample_key] or obs[time_key] columns are absent, the script
falls back to marker-only scoring for all milestones (with a warning).

Status: Step 3 -- marker-score seed labels + Stage 1 conservative routing.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


_VALID_DATASET_IDS = ("GSE178325", "GSE230659")
_DEFAULT_MARKERS_YAML = "benchmark/annotation/milestone_markers.yaml"

# Obs columns that must be present after annotation.
_REQUIRED_OBS_COLS = [
    "milestone_marker_label",
    "milestone_marker_score",
    "milestone_marker_margin",
    "milestone_embedding_label",
    "milestone_embedding_confidence",
    "milestone_classifier_label",
    "milestone_classifier_confidence",
    "consensus_milestone_label",
    "consensus_confidence",
    # Stage 1 authoritative columns.
    "stage1_label",
    "stage1_confidence",
    "stage1_margin",
    "stage1_status",
    "stage1_source",
]

_PLACEHOLDER_NOTE = (
    "milestone_embedding_label, milestone_embedding_confidence, "
    "milestone_classifier_label, milestone_classifier_confidence, "
    "consensus_milestone_label, and consensus_confidence are Step-3 "
    "temporary compatibility placeholders copied from marker-score columns.  "
    "They will be replaced by real embedding, classifier, and consensus "
    "logic in Steps 4-5.  "
    "stage1_label, stage1_confidence, stage1_margin, stage1_status, and "
    "stage1_source are the AUTHORITATIVE Stage 1 outputs; downstream "
    "logic should use these columns rather than the placeholder columns."
)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Annotate observed cells with milestone marker-gene seed labels "
            "(scTimeBench annotation Step 3 / Stage 1)"
        )
    )
    parser.add_argument(
        "--input-h5ad",
        required=True,
        help="Path to input AnnData (.h5ad) with expression matrix.",
    )
    parser.add_argument(
        "--output-h5ad",
        required=True,
        help="Path for the annotated output AnnData (.h5ad).",
    )
    parser.add_argument(
        "--dataset-id",
        required=True,
        choices=list(_VALID_DATASET_IDS),
        help="Dataset identifier used to select the marker config block.",
    )
    parser.add_argument(
        "--markers-yaml",
        default=_DEFAULT_MARKERS_YAML,
        help="Path to milestone_markers.yaml (default: %(default)s).",
    )
    parser.add_argument(
        "--gene-symbol-column",
        default="gene_symbol",
        dest="gene_symbol_column",
        help="Column in adata.var to use as gene symbols (default: gene_symbol).",
    )
    parser.add_argument(
        "--use-var-index-if-needed",
        action="store_true",
        dest="use_var_index",
        help=(
            "Fall back to adata.var_names when --gene-symbol-column is "
            "absent from adata.var."
        ),
    )
    parser.add_argument(
        "--min-marker-overlap",
        type=int,
        default=2,
        dest="min_marker_overlap",
        help=(
            "Minimum matched marker genes required per milestone (default: 2).  "
            "Primary milestones below threshold raise an error."
        ),
    )
    parser.add_argument(
        "--min-top-score",
        type=float,
        default=0.0,
        dest="min_top_score",
        help=(
            "Legacy threshold: cells with top marker score below this are "
            "labelled ambiguous in milestone_marker_label (default: 0.0).  "
            "Does NOT affect stage1_label; use --min-primary-score for Stage 1."
        ),
    )
    parser.add_argument(
        "--min-score-margin",
        type=float,
        default=0.05,
        dest="min_score_margin",
        help=(
            "Legacy threshold: cells where top - second-best score < this are "
            "labelled ambiguous in milestone_marker_label (default: 0.05).  "
            "Does NOT affect stage1_label; use --min-primary-margin for Stage 1."
        ),
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=None,
        dest="max_cells",
        help="Subset to first N cells before scoring (smoke tests only).",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        dest="output_dir",
        help=(
            "Directory for QC output files.  "
            "Defaults to parent directory of --output-h5ad."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Run all computations and print summary but write no output files."
        ),
    )

    # ------------------------------------------------------------------
    # Stage 1: hADSC sample/timepoint assignment
    # ------------------------------------------------------------------
    hadsc_grp = parser.add_argument_group(
        "Stage 1 hADSC assignment",
        description=(
            "hADSCs are assigned primarily by source sample / time == 0.  "
            "If the relevant obs columns are absent the script falls back to "
            "marker-only logic with a warning."
        ),
    )
    hadsc_grp.add_argument(
        "--sample-key",
        default="sample",
        dest="sample_key",
        help="obs column name for sample identity (default: sample).",
    )
    hadsc_grp.add_argument(
        "--time-key",
        default="time",
        dest="time_key",
        help="obs column name for timepoint (default: time).",
    )
    hadsc_grp.add_argument(
        "--hadsc-sample-values",
        default="ADSC,hADSC,hADSCs",
        dest="hadsc_sample_values",
        help=(
            "Comma-separated sample values that identify hADSC cells "
            "(default: ADSC,hADSC,hADSCs)."
        ),
    )
    hadsc_grp.add_argument(
        "--hadsc-time-values",
        default="0,0.0",
        dest="hadsc_time_values",
        help=(
            "Comma-separated timepoint values (numeric or string) that "
            "identify hADSC cells (default: 0,0.0)."
        ),
    )
    hadsc_grp.add_argument(
        "--hadsc-label",
        default="hADSCs",
        dest="hadsc_label",
        help="Stage 1 label string assigned to hADSC cells (default: hADSCs).",
    )

    # ------------------------------------------------------------------
    # Stage 1: primary milestone thresholds and routing labels
    # ------------------------------------------------------------------
    s1_grp = parser.add_argument_group(
        "Stage 1 thresholds and routing",
        description=(
            "Thresholds for routing non-hADSC cells.  All default to 0.0 "
            "(permissive) so that existing runs see minimal change; tighten "
            "for conservative silver-standard label production."
        ),
    )
    s1_grp.add_argument(
        "--min-primary-score",
        type=float,
        default=0.0,
        dest="min_primary_score",
        help=(
            "Minimum top primary marker score for a high-confidence Stage 1 "
            "label (default: 0.0)."
        ),
    )
    s1_grp.add_argument(
        "--min-primary-margin",
        type=float,
        default=0.0,
        dest="min_primary_margin",
        help=(
            "Minimum top - second-best primary score margin for a "
            "high-confidence Stage 1 label (default: 0.0)."
        ),
    )
    s1_grp.add_argument(
        "--ambiguous-label",
        default="ambiguous",
        dest="ambiguous_label",
        help="stage1_label assigned to ambiguous cells (default: ambiguous).",
    )

    return parser


# ---------------------------------------------------------------------------
# Stage 1 helpers
# ---------------------------------------------------------------------------

def _resolve_hadsc_mask(
    adata: object,
    args: argparse.Namespace,
    tag: str,
) -> tuple:
    """Return a boolean mask of hADSC cells identified by sample / timepoint.

    Parameters
    ----------
    adata : AnnData
    args  : parsed CLI namespace
    tag   : log prefix string

    Returns
    -------
    hadsc_mask : np.ndarray[bool]  shape (n_cells,)
    meta       : dict with keys sample_key_found, time_key_found,
                 sample_time_used, n_hadsc_by_sample_time
    """
    import numpy as np

    n_cells = adata.n_obs
    hadsc_mask = np.zeros(n_cells, dtype=bool)
    meta: dict = {
        "sample_key_found": False,
        "time_key_found": False,
        "sample_time_used": False,
        "n_hadsc_by_sample_time": 0,
        "hadsc_sample_values": [],
        "hadsc_time_values": [],
    }

    # Parse sample values.
    hadsc_sample_set = {
        v.strip() for v in args.hadsc_sample_values.split(",") if v.strip()
    }
    meta["hadsc_sample_values"] = sorted(hadsc_sample_set)

    # Parse time values: keep both string and float forms for flexible matching.
    hadsc_time_strings: set = set()
    hadsc_time_floats: set = set()
    for v in args.hadsc_time_values.split(","):
        v = v.strip()
        if v:
            hadsc_time_strings.add(v)
            try:
                hadsc_time_floats.add(float(v))
            except ValueError:
                pass
    meta["hadsc_time_values"] = sorted(hadsc_time_strings)

    obs = adata.obs

    # --- Sample column ---
    if args.sample_key in obs.columns:
        meta["sample_key_found"] = True
        sample_col = obs[args.sample_key].astype(str)
        sample_match = sample_col.isin(hadsc_sample_set).to_numpy()
        hadsc_mask |= sample_match
        n_matched = int(sample_match.sum())
        if n_matched > 0:
            print(
                f"{tag} sample-key '{args.sample_key}': "
                f"{n_matched} cells matched {sorted(hadsc_sample_set)}"
            )
    else:
        print(
            f"{tag} WARN: obs column '{args.sample_key}' not found; "
            "sample-based hADSC assignment skipped."
        )

    # --- Time column ---
    if args.time_key in obs.columns:
        meta["time_key_found"] = True
        time_col = obs[args.time_key]

        # Numeric match.
        try:
            time_float_arr = time_col.astype(float).to_numpy()
            if hadsc_time_floats:
                time_match_num = np.isin(
                    time_float_arr,
                    np.array(sorted(hadsc_time_floats), dtype=float),
                )
            else:
                time_match_num = np.zeros(n_cells, dtype=bool)
        except (ValueError, TypeError):
            time_match_num = np.zeros(n_cells, dtype=bool)

        # String match (catches non-numeric time labels like "D0").
        time_str_col = time_col.astype(str)
        time_match_str = time_str_col.isin(hadsc_time_strings).to_numpy()

        time_match = time_match_num | time_match_str
        hadsc_mask |= time_match
        n_matched = int(time_match.sum())
        if n_matched > 0:
            print(
                f"{tag} time-key '{args.time_key}': "
                f"{n_matched} cells matched {sorted(hadsc_time_strings)}"
            )
    else:
        print(
            f"{tag} WARN: obs column '{args.time_key}' not found; "
            "timepoint-based hADSC assignment skipped."
        )

    meta["sample_time_used"] = (
        meta["sample_key_found"] or meta["time_key_found"]
    )
    meta["n_hadsc_by_sample_time"] = int(hadsc_mask.sum())

    return hadsc_mask, meta


def _assign_stage1_labels(
    n_cells: int,
    hadsc_mask: object,
    score_dict: dict,
    primary_for_stage1: list,
    old_margins: object,
    args: argparse.Namespace,
) -> tuple:
    """Compute Stage 1 label arrays for all cells.

    Assignment priority (highest to lowest):
      1. hADSC by sample/time -> stage1_source = "sample_time"
      2. High-confidence primary -> stage1_source = "marker_primary"
      3. Ambiguous -> stage1_source = "marker_ambiguous"

    Parameters
    ----------
    n_cells             : int
    hadsc_mask          : np.ndarray[bool], shape (n_cells,)
    score_dict          : dict[str, np.ndarray] — all scored milestones
    primary_milestones  : list[str] — all primary milestone names
    primary_for_stage1  : list[str] — primary milestones that compete for
                          non-hADSC cells (excludes hADSCs if sample/time used)
    old_margins         : np.ndarray[float] | None — margins from assign_labels
                          (used as fallback margin for hADSC sample/time cells)
    args                : CLI namespace (for threshold and label-string values)

    Returns
    -------
    stage1_label, stage1_confidence, stage1_margin, stage1_status, stage1_source
        Each is a np.ndarray of shape (n_cells,).
    """
    import numpy as np

    # Initialise to ambiguous defaults for all cells.
    stage1_label = np.full(n_cells, args.ambiguous_label, dtype=object)
    stage1_confidence = np.zeros(n_cells, dtype=np.float64)
    stage1_margin = np.zeros(n_cells, dtype=np.float64)
    stage1_status = np.full(n_cells, "ambiguous", dtype=object)
    stage1_source = np.full(n_cells, "marker_ambiguous", dtype=object)

    # ------------------------------------------------------------------
    # Step A: assign hADSC cells by sample / timepoint
    # ------------------------------------------------------------------
    if hadsc_mask.any():
        stage1_label[hadsc_mask] = args.hadsc_label
        stage1_status[hadsc_mask] = "high_confidence"
        stage1_source[hadsc_mask] = "sample_time"
        stage1_confidence[hadsc_mask] = 1.0
        # Margin: use legacy marker margin where available; fall back to 1.0.
        if old_margins is not None:
            _m = np.asarray(old_margins, dtype=np.float64)
            hadsc_margin = np.where(
                np.isnan(_m[hadsc_mask]) | np.isinf(_m[hadsc_mask]),
                1.0,
                np.abs(_m[hadsc_mask]),
            )
            stage1_margin[hadsc_mask] = hadsc_margin
        else:
            stage1_margin[hadsc_mask] = 1.0

    non_hadsc = ~hadsc_mask
    if not non_hadsc.any():
        return (
            stage1_label.astype(str),
            stage1_confidence.astype(np.float32),
            _clamp_margin(stage1_margin).astype(np.float32),
            stage1_status.astype(str),
            stage1_source.astype(str),
        )

    # ------------------------------------------------------------------
    # Step B: build primary score matrix for non-hADSC competition
    # ------------------------------------------------------------------
    primary_stage1_available = [m for m in primary_for_stage1 if m in score_dict]

    if not primary_stage1_available:
        # No primary markers scored: non-hADSC cells remain ambiguous.
        return (
            stage1_label.astype(str),
            stage1_confidence.astype(np.float32),
            _clamp_margin(stage1_margin).astype(np.float32),
            stage1_status.astype(str),
            stage1_source.astype(str),
        )

    prim_mat = np.column_stack(
        [score_dict[m] for m in primary_stage1_available]
    ).astype(np.float64)  # (n_cells, n_prim)

    top_idx = np.argmax(prim_mat, axis=1)                    # (n_cells,)
    top_prim_scores = prim_mat[np.arange(n_cells), top_idx]  # (n_cells,)

    if prim_mat.shape[1] >= 2:
        tmp = prim_mat.copy()
        tmp[np.arange(n_cells), top_idx] = -np.inf
        second_prim_scores = np.max(tmp, axis=1)
        prim_margins = top_prim_scores - second_prim_scores
    else:
        prim_margins = np.full(n_cells, np.inf)

    prim_name_arr = np.array(primary_stage1_available, dtype=object)
    top_prim_labels = prim_name_arr[top_idx]

    # ------------------------------------------------------------------
    # Step C: vectorised routing for non-hADSC cells
    # ------------------------------------------------------------------
    # High confidence: meets score + margin thresholds.
    high_conf = (
        (top_prim_scores >= args.min_primary_score)
        & (prim_margins >= args.min_primary_margin)
        & non_hadsc
    )
    # Ambiguous: everything else that is non-hADSC
    ambig = non_hadsc & ~high_conf

    # Apply high confidence
    stage1_label[high_conf] = top_prim_labels[high_conf]
    stage1_status[high_conf] = "high_confidence"
    stage1_source[high_conf] = "marker_primary"
    stage1_confidence[high_conf] = top_prim_scores[high_conf]
    stage1_margin[high_conf] = prim_margins[high_conf]

    # Apply ambiguous (write scores even though label is ambiguous)
    stage1_label[ambig] = args.ambiguous_label
    stage1_status[ambig] = "ambiguous"
    stage1_source[ambig] = "marker_ambiguous"
    stage1_confidence[ambig] = top_prim_scores[ambig]
    stage1_margin[ambig] = prim_margins[ambig]

    return (
        stage1_label.astype(str),
        stage1_confidence.astype(np.float32),
        _clamp_margin(stage1_margin).astype(np.float32),
        stage1_status.astype(str),
        stage1_source.astype(str),
    )


def _clamp_margin(arr: object) -> object:
    """Replace inf/nan in margin array with 1.0 / 0.0 for safe h5ad storage."""
    import numpy as np
    arr = np.asarray(arr, dtype=np.float64)
    arr = np.where(np.isinf(arr), 1.0, arr)
    arr = np.where(np.isnan(arr), 0.0, arr)
    return arr


# ---------------------------------------------------------------------------
# Core annotation pipeline
# ---------------------------------------------------------------------------

def _run_annotation(args: argparse.Namespace) -> dict:
    """Execute the full annotation pipeline.  Return metadata dict.

    Imports of anndata / numpy / pandas are deferred to here so the module
    remains importable without those packages installed.
    """
    import numpy as np

    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError(
            "anndata is required for h5ad I/O.  Install: pip install anndata"
        ) from exc

    try:
        import pandas as pd
    except ImportError as exc:
        raise ImportError(
            "pandas is required for CSV output.  Install: pip install pandas"
        ) from exc

    from benchmark.annotation.marker_utils import (
        load_markers_yaml,
        build_gene_index_map,
        resolve_milestone_genes,
        compute_mean_marker_scores,
        assign_labels,
    )

    tag = "[build_marker_seed_labels]"

    # ------------------------------------------------------------------
    # 1. Load AnnData
    # ------------------------------------------------------------------
    input_path = Path(args.input_h5ad)
    print(f"{tag} loading {input_path} ...")
    adata = ad.read_h5ad(str(input_path))
    print(f"{tag} loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # ------------------------------------------------------------------
    # 2. Optionally subset cells (smoke test)
    # ------------------------------------------------------------------
    if args.max_cells is not None and adata.n_obs > args.max_cells:
        print(
            f"{tag} --max-cells {args.max_cells}: subsetting from "
            f"{adata.n_obs} to {args.max_cells} cells."
        )
        adata = adata[:args.max_cells, :].copy()

    # ------------------------------------------------------------------
    # 3. Load marker config
    # ------------------------------------------------------------------
    markers_yaml_path = Path(args.markers_yaml)
    print(f"{tag} loading markers from {markers_yaml_path} ...")
    marker_config = load_markers_yaml(markers_yaml_path)

    if args.dataset_id not in marker_config:
        raise ValueError(
            f"Dataset '{args.dataset_id}' not found in {markers_yaml_path}.  "
            f"Available: {list(marker_config.keys())}"
        )

    ds_cfg = marker_config[args.dataset_id]
    primary_milestones = list(ds_cfg.get("primary_milestones") or [])

    print(f"{tag} milestones: primary={primary_milestones}")

    # ------------------------------------------------------------------
    # 4. Resolve gene symbols to var column indices
    # ------------------------------------------------------------------
    gene_index_map = build_gene_index_map(
        var_names=adata.var_names,
        var_df=adata.var,
        gene_symbol_column=args.gene_symbol_column,
        use_var_index=args.use_var_index,
    )
    print(f"{tag} gene index map: {len(gene_index_map)} unique symbols.")

    # ------------------------------------------------------------------
    # 5. Resolve primary milestone gene indices
    # ------------------------------------------------------------------
    primary_indices, primary_overlap = resolve_milestone_genes(
        marker_config=marker_config,
        dataset_id=args.dataset_id,
        milestones=primary_milestones,
        gene_index_map=gene_index_map,
        min_overlap=args.min_marker_overlap,
        strict=True,
    )

    overlap_summary = dict(primary_overlap)

    for ms, summ in primary_overlap.items():
        print(
            f"{tag}   {ms}: {summ['n_matched']}/{summ['n_required']} genes "
            f"matched ({summ['matched_genes']})"
        )

    # ------------------------------------------------------------------
    # 6. Compute marker scores (primary milestones only)
    # ------------------------------------------------------------------
    all_indices_for_scoring = dict(primary_indices)
    print(
        f"{tag} scoring {len(all_indices_for_scoring)} milestones "
        f"over {adata.n_obs} cells ..."
    )
    score_dict = compute_mean_marker_scores(adata.X, all_indices_for_scoring)

    # ------------------------------------------------------------------
    # 7. Legacy label assignment (milestone_marker_label and friends)
    #    Uses all primary milestones; thresholds from --min-top-score /
    #    --min-score-margin.  Backward-compatible; unchanged.
    # ------------------------------------------------------------------
    labels, top_scores, margins = assign_labels(
        score_dict={m: score_dict[m] for m in primary_indices if m in score_dict},
        milestone_names=list(primary_indices.keys()),
        min_top_score=args.min_top_score,
        min_score_margin=args.min_score_margin,
    )

    # ------------------------------------------------------------------
    # 8. Resolve hADSC mask (sample / timepoint first)
    # ------------------------------------------------------------------
    hadsc_mask, hadsc_meta = _resolve_hadsc_mask(adata, args, tag)
    sample_time_used: bool = hadsc_meta["sample_time_used"]

    if not sample_time_used:
        print(
            f"{tag} WARN: neither '{args.sample_key}' nor '{args.time_key}' "
            "obs columns found.  Falling back to marker-only hADSC assignment; "
            "hADSCs will compete as a primary milestone in Stage 1 scoring."
        )

    # Primary milestones competing for non-hADSC cells in Stage 1:
    # exclude hADSCs from the marker competition when sample/time was used.
    if sample_time_used:
        primary_for_stage1 = [
            m for m in primary_milestones if m != args.hadsc_label
        ]
    else:
        primary_for_stage1 = list(primary_milestones)

    print(
        f"{tag} Stage 1: sample_time_used={sample_time_used}, "
        f"hadsc_by_sample_time={hadsc_meta['n_hadsc_by_sample_time']}, "
        f"primary_for_stage1={primary_for_stage1}"
    )

    # ------------------------------------------------------------------
    # 9. Assign Stage 1 labels
    # ------------------------------------------------------------------
    (
        stage1_label,
        stage1_confidence,
        stage1_margin,
        stage1_status,
        stage1_source,
    ) = _assign_stage1_labels(
        n_cells=adata.n_obs,
        hadsc_mask=hadsc_mask,
        score_dict=score_dict,
        primary_for_stage1=primary_for_stage1,
        old_margins=margins,
        args=args,
    )

    # ------------------------------------------------------------------
    # 10. Write obs columns
    # ------------------------------------------------------------------
    # Per-milestone score columns: primary milestones only.
    output_milestones = set(primary_indices.keys())
    for ms, arr in score_dict.items():
        if ms in output_milestones:
            adata.obs[f"milestone_score_{ms}"] = arr.astype(np.float32)

    # Legacy marker columns (backward-compatible).
    adata.obs["milestone_marker_label"] = labels.astype(str)
    adata.obs["milestone_marker_score"] = top_scores.astype(np.float32)
    adata.obs["milestone_marker_margin"] = margins.astype(np.float32)

    # Step-3 compatibility placeholders.
    adata.obs["milestone_embedding_label"] = labels.astype(str)
    adata.obs["milestone_embedding_confidence"] = top_scores.astype(np.float32)
    adata.obs["milestone_classifier_label"] = labels.astype(str)
    adata.obs["milestone_classifier_confidence"] = top_scores.astype(np.float32)
    adata.obs["consensus_milestone_label"] = labels.astype(str)
    adata.obs["consensus_confidence"] = top_scores.astype(np.float32)

    # Stage 1 authoritative columns.
    adata.obs["stage1_label"] = stage1_label
    adata.obs["stage1_confidence"] = stage1_confidence
    adata.obs["stage1_margin"] = stage1_margin
    adata.obs["stage1_status"] = stage1_status
    adata.obs["stage1_source"] = stage1_source

    # ------------------------------------------------------------------
    # 11. Count summaries
    # ------------------------------------------------------------------
    label_counts = adata.obs["milestone_marker_label"].value_counts().to_dict()
    stage1_label_counts = adata.obs["stage1_label"].value_counts().to_dict()
    stage1_status_counts = adata.obs["stage1_status"].value_counts().to_dict()
    stage1_source_counts = adata.obs["stage1_source"].value_counts().to_dict()

    print(f"{tag} legacy marker label counts: {label_counts}")
    print(f"{tag} stage1_label counts:        {stage1_label_counts}")
    print(f"{tag} stage1_status counts:       {stage1_status_counts}")

    n_total = adata.n_obs
    n_sample_time = int(hadsc_meta["n_hadsc_by_sample_time"])
    n_high_conf = int(stage1_status_counts.get("high_confidence", 0))
    n_ambig = int(stage1_status_counts.get("ambiguous", 0))

    # ------------------------------------------------------------------
    # 12. uns metadata
    # ------------------------------------------------------------------
    uns_meta = _build_uns_metadata(
        args=args,
        primary_milestones=primary_milestones,
        overlap_summary=overlap_summary,
        label_counts=label_counts,
        stage1_label_counts=stage1_label_counts,
        stage1_status_counts=stage1_status_counts,
        stage1_source_counts=stage1_source_counts,
        hadsc_meta=hadsc_meta,
        primary_for_stage1=primary_for_stage1,
        n_cells=adata.n_obs,
        n_genes=adata.n_vars,
        markers_yaml_path=markers_yaml_path,
        n_total=n_total,
        n_sample_time=n_sample_time,
        n_high_conf=n_high_conf,
        n_ambig=n_ambig,
    )
    adata.uns["milestone_annotation_step3"] = uns_meta

    # ------------------------------------------------------------------
    # 13. Write outputs (skip in dry-run)
    # ------------------------------------------------------------------
    out_h5ad = Path(args.output_h5ad)
    out_dir = Path(args.output_dir) if args.output_dir else out_h5ad.parent

    if args.dry_run:
        print(f"{tag} dry-run: skipping all file writes.")
        print(f"{tag} would write annotated h5ad -> {out_h5ad}")
        print(f"{tag} would write QC files in  -> {out_dir}")
    else:
        out_dir.mkdir(parents=True, exist_ok=True)
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)

        print(f"{tag} writing annotated h5ad -> {out_h5ad}")
        adata.write_h5ad(str(out_h5ad))

        meta_path = out_dir / "marker_seed_labels_metadata.json"
        with meta_path.open("w", encoding="ascii") as fh:
            json.dump(uns_meta, fh, indent=2, default=_json_default)
        print(f"{tag} metadata           -> {meta_path}")

        # Legacy label counts CSV.
        counts_path = out_dir / "marker_seed_label_counts.csv"
        counts_df = pd.DataFrame(
            list(label_counts.items()), columns=["milestone_label", "cell_count"]
        ).sort_values("cell_count", ascending=False)
        counts_df.to_csv(str(counts_path), index=False)
        print(f"{tag} label counts       -> {counts_path}")

        # Marker overlap CSV.
        overlap_path = out_dir / "marker_seed_marker_overlap.csv"
        overlap_rows = []
        for ms, summ in overlap_summary.items():
            overlap_rows.append({
                "milestone": ms,
                "n_required": summ["n_required"],
                "n_matched": summ["n_matched"],
                "n_column_indices": summ["n_column_indices"],
                "matched_genes": ";".join(summ["matched_genes"]),
                "missing_genes": ";".join(summ["missing_genes"]),
            })
        pd.DataFrame(overlap_rows).to_csv(str(overlap_path), index=False)
        print(f"{tag} marker overlap     -> {overlap_path}")

        # Stage 1 QC CSVs.
        s1_label_path = out_dir / "stage1_label_counts.csv"
        pd.DataFrame(
            list(stage1_label_counts.items()),
            columns=["stage1_label", "cell_count"],
        ).sort_values("cell_count", ascending=False).to_csv(
            str(s1_label_path), index=False
        )
        print(f"{tag} stage1_label counts -> {s1_label_path}")

        s1_status_path = out_dir / "stage1_status_counts.csv"
        pd.DataFrame(
            list(stage1_status_counts.items()),
            columns=["stage1_status", "cell_count"],
        ).sort_values("cell_count", ascending=False).to_csv(
            str(s1_status_path), index=False
        )
        print(f"{tag} stage1_status counts -> {s1_status_path}")

        s1_source_path = out_dir / "stage1_source_counts.csv"
        pd.DataFrame(
            list(stage1_source_counts.items()),
            columns=["stage1_source", "cell_count"],
        ).sort_values("cell_count", ascending=False).to_csv(
            str(s1_source_path), index=False
        )
        print(f"{tag} stage1_source counts -> {s1_source_path}")

        # Post-write validation.
        _validate_output(out_h5ad, tag)

    return uns_meta


def _build_uns_metadata(
    args: argparse.Namespace,
    primary_milestones: list,
    overlap_summary: dict,
    label_counts: dict,
    stage1_label_counts: dict,
    stage1_status_counts: dict,
    stage1_source_counts: dict,
    hadsc_meta: dict,
    primary_for_stage1: list,
    n_cells: int,
    n_genes: int,
    markers_yaml_path: Path,
    n_total: int,
    n_sample_time: int,
    n_high_conf: int,
    n_ambig: int,
) -> dict:
    def _frac(n: int) -> float:
        return round(n / n_total, 4) if n_total > 0 else 0.0

    return {
        "script": "build_marker_seed_labels",
        "status": "completed",
        "annotation_stage": "marker_seed_labels_step3",
        # Stage 1 policy block.
        "stage1_policy": {
            "name": "marker_sample_time_stage1_v1",
            "hadsc_assignment": (
                "primary assignment by source sample/timepoint (obs columns); "
                "marker gene set is auxiliary somatic/fibroblast-like support only"
            ),
            "hadsc_label": args.hadsc_label,
            "sample_key": args.sample_key,
            "time_key": args.time_key,
            "hadsc_sample_values": hadsc_meta["hadsc_sample_values"],
            "hadsc_time_values": hadsc_meta["hadsc_time_values"],
            "sample_key_found": hadsc_meta["sample_key_found"],
            "time_key_found": hadsc_meta["time_key_found"],
            "sample_time_used": hadsc_meta["sample_time_used"],
            "primary_for_stage1": primary_for_stage1,
            "thresholds": {
                "min_primary_score": args.min_primary_score,
                "min_primary_margin": args.min_primary_margin,
                "ambiguous_label": args.ambiguous_label,
            },
        },
        "dataset_id": args.dataset_id,
        "markers_yaml": str(markers_yaml_path),
        "primary_milestones": primary_milestones,
        "min_marker_overlap": args.min_marker_overlap,
        # Legacy thresholds (milestone_marker_label only).
        "min_top_score": args.min_top_score,
        "min_score_margin": args.min_score_margin,
        "max_cells_used": args.max_cells,
        "placeholder_columns_created": True,
        "placeholder_note": _PLACEHOLDER_NOTE,
        "marker_overlap_summary": {
            ms: {
                "n_required": s["n_required"],
                "n_matched": s["n_matched"],
                "matched_genes": s["matched_genes"],
                "missing_genes": s["missing_genes"],
            }
            for ms, s in overlap_summary.items()
        },
        # Legacy label counts (milestone_marker_label).
        "label_counts": {str(k): int(v) for k, v in label_counts.items()},
        # Stage 1 counts.
        "stage1_label_counts": {
            str(k): int(v) for k, v in stage1_label_counts.items()
        },
        "stage1_status_counts": {
            str(k): int(v) for k, v in stage1_status_counts.items()
        },
        "stage1_source_counts": {
            str(k): int(v) for k, v in stage1_source_counts.items()
        },
        "stage1_summary": {
            "n_cells": n_total,
            "n_assigned_by_sample_time": n_sample_time,
            "frac_assigned_by_sample_time": _frac(n_sample_time),
            "n_high_confidence": n_high_conf,
            "frac_high_confidence": _frac(n_high_conf),
            "n_ambiguous": n_ambig,
            "frac_ambiguous": _frac(n_ambig),
        },
        "n_cells": n_cells,
        "n_genes": n_genes,
    }


def _validate_output(out_h5ad: Path, tag: str) -> None:
    """Reload the written h5ad and assert all required obs columns exist."""
    try:
        import anndata as ad
    except ImportError:
        print(f"{tag} WARN: anndata not available for post-write validation.")
        return
    print(f"{tag} validating output h5ad ...")
    a = ad.read_h5ad(str(out_h5ad), backed="r")
    missing = [c for c in _REQUIRED_OBS_COLS if c not in a.obs.columns]
    if missing:
        raise RuntimeError(
            f"Post-write validation FAILED.  "
            f"Required obs columns missing: {missing}"
        )
    print(
        f"{tag} post-write validation PASSED "
        f"({len(_REQUIRED_OBS_COLS)} required columns present)."
    )


def _json_default(obj: object) -> object:
    """JSON serialiser fallback for numpy scalar types."""
    try:
        import numpy as np
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
    except ImportError:
        pass
    return str(obj)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    tag = "[build_marker_seed_labels]"
    if args.dry_run:
        print(f"{tag} dry-run mode active: no output files will be written.")

    try:
        _run_annotation(args)
    except Exception as exc:  # noqa: BLE001
        print(f"{tag} ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"{tag} done.")
    sys.exit(0)


if __name__ == "__main__":
    main()
