"""build_trajectory_aware_labels.py
==================================
Stage 2 of the trajectory annotation layer: resolves ambiguous Stage 1 cells
using trajectory-aware marker-score and optional embedding evidence, and
produces frozen coarse and expanded milestone labels for scTimeBench metrics.

Input
-----
An .h5ad file annotated by build_marker_seed_labels.py (Stage 1), containing:
  stage1_label, stage1_status, stage1_confidence, stage1_margin
  milestone_score_<milestone>  (one column per scored primary milestone)

Output obs columns
------------------
  trajectory_membership_label     str
      in_trajectory | unknown_or_ood | ambiguous
  trajectory_membership_score     float32
  trajectory_membership_margin    float32
  stage2_label_coarse             str
      Stage 2 decision for previously-ambiguous cells; for Stage 1
      high-confidence cells stage2_source is "not_applicable_stage1_final".
  stage2_label_expanded           str
  stage2_confidence               float32
  stage2_source                   str
  final_milestone_label_coarse    str   FROZEN coarse label (official metrics)
  final_milestone_label_expanded  str   FROZEN expanded label (supplementary)
  final_milestone_confidence      float32
  final_milestone_source          str

Design rules
------------
- Do not use embeddings from the model being evaluated.
- Do not force all cells into primary milestones.
- unknown_or_ood and ambiguous are valid final labels.
- final_milestone_label_coarse drives official scTimeBench-style metrics.
- final_milestone_label_expanded is supplementary only.

Marker-only mode
----------------
If --embedding-key is omitted, Stage 2 resolves ambiguous cells using marker
scores only.  This is the default and fully supported mode.

Provenance
----------
  Liuyang et al. 2023, Cell Stem Cell, DOI: 10.1016/j.stem.2023.02.008
  Guan et al. 2022, Nature, DOI: 10.1038/s41586-022-04593-5
  Companion code: https://github.com/sajuukLyu/CSC_2023

Usage
-----
  python benchmark/annotation/build_trajectory_aware_labels.py \\
      --dataset-id GSE230659 \\
      --input-h5ad  path/to/gse230659_stage1.h5ad \\
      --output-h5ad path/to/gse230659_stage2.h5ad \\
      --output-dir  path/to/stage2_output/

Status: Stage 2 -- trajectory-aware label resolution.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


_VALID_DATASET_IDS = ("GSE178325", "GSE230659")
_DEFAULT_MARKERS_YAML = "benchmark/annotation/milestone_markers.yaml"

_OUTPUT_OBS_COLS = [
    "trajectory_membership_label",
    "trajectory_membership_score",
    "trajectory_membership_margin",
    "stage2_label_coarse",
    "stage2_label_expanded",
    "stage2_confidence",
    "stage2_source",
    "final_milestone_label_coarse",
    "final_milestone_label_expanded",
    "final_milestone_confidence",
    "final_milestone_source",
]

_REQUIRED_STAGE1_COLS = [
    "stage1_label",
    "stage1_status",
    "stage1_confidence",
    "stage1_margin",
]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Stage 2 trajectory-aware milestone label resolution "
            "(scTimeBench annotation framework v2)"
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    req = parser.add_argument_group("required")
    req.add_argument(
        "--input-h5ad", required=True,
        help="Path to Stage 1-annotated input AnnData (.h5ad).",
    )
    req.add_argument(
        "--output-h5ad", required=True,
        help="Path for the Stage 2-annotated output AnnData (.h5ad).",
    )
    req.add_argument(
        "--dataset-id", required=True, choices=list(_VALID_DATASET_IDS),
        help="Dataset identifier (selects config block from markers YAML).",
    )
    req.add_argument(
        "--output-dir", required=True,
        help="Directory for QC CSVs, metadata JSON, and other output files.",
    )

    parser.add_argument(
        "--markers-yaml", default=_DEFAULT_MARKERS_YAML,
        help="Path to milestone_markers.yaml.",
    )
    parser.add_argument(
        "--embedding-key", default=None, dest="embedding_key",
        help=(
            "adata.obsm key for an INDEPENDENT reference embedding "
            "(e.g. X_pca, X_scGPT).  Must NOT be from the model being "
            "evaluated.  If omitted or absent, marker-only mode is used."
        ),
    )

    # Stage 1 column keys
    s1g = parser.add_argument_group("Stage 1 column keys")
    s1g.add_argument(
        "--stage1-label-key", default="stage1_label", dest="stage1_label_key",
    )
    s1g.add_argument(
        "--stage1-status-key", default="stage1_status", dest="stage1_status_key",
    )
    s1g.add_argument(
        "--stage1-confidence-key", default="stage1_confidence",
        dest="stage1_confidence_key",
    )
    s1g.add_argument(
        "--stage1-margin-key", default="stage1_margin", dest="stage1_margin_key",
    )

    # Label strings
    lbl = parser.add_argument_group("label strings")
    lbl.add_argument(
        "--ambiguous-label", default="ambiguous", dest="ambiguous_label",
    )
    lbl.add_argument(
        "--unknown-label", default="unknown_or_ood", dest="unknown_label",
    )

    # Marker thresholds
    mth = parser.add_argument_group("marker thresholds")
    mth.add_argument(
        "--min-single-score", type=float, default=0.0, dest="min_single_score",
        help="Minimum top primary score required for single-milestone assignment.",
    )
    mth.add_argument(
        "--min-transition-score", type=float, default=0.0,
        dest="min_transition_score",
        help="Minimum score required for both top1 and top2 in a transition assignment.",
    )
    mth.add_argument(
        "--transition-margin", type=float, default=0.10, dest="transition_margin",
        help=(
            "Scores within this margin are treated as a potential transition; "
            "scores exceeding it trigger single-milestone assignment."
        ),
    )

    # Embedding thresholds
    eth = parser.add_argument_group("embedding thresholds")
    eth.add_argument(
        "--ood-distance-quantile", type=float, default=0.95,
        dest="ood_distance_quantile",
        help=(
            "Percentile of within-centroid reference distances used as "
            "the per-milestone OOD threshold (0-1)."
        ),
    )
    eth.add_argument(
        "--min-reference-cells", type=int, default=5,
        dest="min_reference_cells",
        help=(
            "Minimum Stage 1 high-confidence cells required per milestone "
            "to build a centroid."
        ),
    )
    eth.add_argument(
        "--transition-position-low", type=float, default=0.25,
        dest="transition_position_low",
        help=(
            "Normalised position threshold below which a cell is assigned "
            "to the nearer milestone (0 = at A, 1 = at B)."
        ),
    )
    eth.add_argument(
        "--transition-position-high", type=float, default=0.75,
        dest="transition_position_high",
        help=(
            "Normalised position threshold above which a cell is assigned "
            "to the farther milestone."
        ),
    )

    parser.add_argument(
        "--dry-run", action="store_true",
        help=(
            "Validate inputs and run all computation but do not write "
            "output_h5ad.  Metadata JSON is written with status 'dry_run'."
        ),
    )
    return parser


# ---------------------------------------------------------------------------
# Config loading and trajectory structure building
# ---------------------------------------------------------------------------

def _load_yaml(path: Path) -> dict:
    try:
        import yaml
    except ImportError as exc:
        raise ImportError(
            "PyYAML is required.  Install: pip install pyyaml"
        ) from exc
    with path.open(encoding="ascii") as fh:
        return yaml.safe_load(fh)


def _load_config(markers_yaml: Path, dataset_id: str) -> dict:
    """Load and return the dataset-specific config block from markers YAML."""
    if not markers_yaml.exists():
        raise FileNotFoundError(f"Markers YAML not found: {markers_yaml}")
    data = _load_yaml(markers_yaml)
    if not isinstance(data, dict) or dataset_id not in data:
        raise ValueError(
            f"Dataset '{dataset_id}' not found in {markers_yaml}.  "
            f"Available: {list(data.keys()) if isinstance(data, dict) else '(invalid YAML)'}"
        )
    return data[dataset_id]


def _build_trajectory_structures(ds_cfg: dict) -> dict:
    """Build adjacency sets and expanded-label maps from dataset config.

    Returns a dict with keys:
      primary_milestones        list[str]
      primary_set               set[str]
      adjacency_directed        set of (from, to) tuples
      adjacency_undirected      set of (a, b) tuples (both directions stored)
      coarse_to_expanded        dict[milestone_name, expanded_stage_name]
      transition_to_expanded    dict[(ms1, ms2), expanded_transition_name]
                                  both (a,b) and (b,a) keys stored
      expanded_trajectory_order list[str]
      excluded_from_official_metrics list[str]
    """
    primary = list(ds_cfg.get("primary_milestones") or [])
    edges = ds_cfg.get("reference_graph_edges") or []
    expanded = list(ds_cfg.get("expanded_trajectory_order") or [])
    excluded = list(ds_cfg.get("excluded_from_official_metrics") or [])

    # Directed adjacency from reference_graph_edges.
    adjacency_directed: set = set()
    for edge in edges:
        if isinstance(edge, (list, tuple)) and len(edge) == 2:
            adjacency_directed.add((str(edge[0]), str(edge[1])))

    # Undirected: store both orderings for fast membership test.
    adjacency_undirected: set = set()
    for a, b in adjacency_directed:
        adjacency_undirected.add((a, b))
        adjacency_undirected.add((b, a))

    # Parse expanded_trajectory_order.
    # Convention: "stage_XX_<semantic>" where <semantic> contains "_to_"
    # for transition stages and is a bare milestone name for coarse stages.
    coarse_to_expanded: dict = {}
    transition_to_expanded: dict = {}

    for stage_name in expanded:
        # Strip leading "stage_XX_" prefix (exactly 2 underscore-delimited tokens).
        parts = stage_name.split("_", maxsplit=2)
        if len(parts) < 3:
            continue
        semantic = parts[2]  # e.g. "hADSCs", "hADSCs_to_epithelial_like"
        if "_to_" in semantic:
            idx = semantic.index("_to_")
            ms1 = semantic[:idx]
            ms2 = semantic[idx + 4:]
            transition_to_expanded[(ms1, ms2)] = stage_name
            transition_to_expanded[(ms2, ms1)] = stage_name  # bidirectional
        else:
            coarse_to_expanded[semantic] = stage_name

    return {
        "primary_milestones": primary,
        "primary_set": set(primary),
        "adjacency_directed": adjacency_directed,
        "adjacency_undirected": adjacency_undirected,
        "coarse_to_expanded": coarse_to_expanded,
        "transition_to_expanded": transition_to_expanded,
        "expanded_trajectory_order": expanded,
        "excluded_from_official_metrics": excluded,
    }


# ---------------------------------------------------------------------------
# Small label helpers
# ---------------------------------------------------------------------------

def _are_adjacent(ms1: str, ms2: str, adjacency_undirected: set) -> bool:
    return (ms1, ms2) in adjacency_undirected


def _to_expanded(coarse: str, coarse_to_expanded: dict) -> str:
    """Map a coarse milestone name to its expanded stage name, or return as-is."""
    return coarse_to_expanded.get(coarse, coarse)


def _membership_label(
    coarse: str,
    primary_set: set,
    unknown_label: str,
    ambiguous_label: str,  # noqa: ARG001 (used as fallback)
) -> str:
    """Map a coarse final label to a trajectory membership category."""
    if coarse in primary_set:
        return "in_trajectory"
    if coarse == unknown_label:
        return "unknown_or_ood"
    return "ambiguous"


# ---------------------------------------------------------------------------
# Resolution helper: apply a single cell's resolved dict into output arrays
# ---------------------------------------------------------------------------

def _apply_resolution(
    global_i: int,
    res: dict,
    final_coarse: object,
    final_expanded: object,
    final_conf: object,
    final_source: object,
    stage2_coarse: object,
    stage2_expanded: object,
    stage2_conf: object,
    stage2_source: object,
    traj_label: object,
    traj_score: object,
    traj_margin: object,
    primary_set: set,
    args: argparse.Namespace,
) -> None:
    coarse = res["coarse"]
    expanded = res["expanded"]
    conf = float(res["confidence"])
    margin = float(res["margin"])
    source = res["source"]

    final_coarse[global_i] = coarse
    final_expanded[global_i] = expanded
    final_conf[global_i] = conf
    final_source[global_i] = source

    stage2_coarse[global_i] = coarse
    stage2_expanded[global_i] = expanded
    stage2_conf[global_i] = conf
    stage2_source[global_i] = source

    traj_label[global_i] = _membership_label(
        coarse, primary_set,
        args.unknown_label, args.ambiguous_label,
    )
    traj_score[global_i] = conf
    traj_margin[global_i] = margin


# ---------------------------------------------------------------------------
# Stage 2 resolution: marker-score evidence
# ---------------------------------------------------------------------------

def _resolve_by_markers(
    ambig_idx: object,
    obs: object,
    traj: dict,
    args: argparse.Namespace,
    warnings: list,
) -> dict:
    """Resolve ambiguous cells using milestone_score_* obs columns.

    Returns a dict mapping global cell index -> resolution dict:
      {coarse, expanded, confidence, margin, source}
    Only resolved cells are included; unresolved cells are absent.
    """
    import numpy as np

    resolutions: dict = {}
    if len(ambig_idx) == 0:
        return resolutions

    primary = traj["primary_milestones"]
    adj = traj["adjacency_undirected"]
    c2e = traj["coarse_to_expanded"]
    t2e = traj["transition_to_expanded"]

    # Collect available primary score arrays (written by Stage 1).
    prim_scores: dict = {}
    missing_prim: list = []
    for ms in primary:
        col = f"milestone_score_{ms}"
        if col in obs.columns:
            prim_scores[ms] = np.nan_to_num(
                obs[col].to_numpy(dtype=np.float64), nan=0.0
            )
        else:
            missing_prim.append(ms)
    if missing_prim:
        warnings.append(
            f"Stage 2 marker: milestone_score columns missing for "
            f"{missing_prim}; those milestones will not compete."
        )

    if not prim_scores:
        warnings.append(
            "Stage 2 marker: no milestone_score columns available; "
            "skipping marker-based resolution."
        )
        return resolutions

    prim_names = list(prim_scores.keys())
    n_ambig = len(ambig_idx)
    prim_mat_sub = np.column_stack(
        [prim_scores[m] for m in prim_names]
    )[ambig_idx]  # (n_ambig, n_prim)

    top1_idx = np.argmax(prim_mat_sub, axis=1)
    top1_scores = prim_mat_sub[np.arange(n_ambig), top1_idx]

    if prim_mat_sub.shape[1] >= 2:
        tmp = prim_mat_sub.copy()
        tmp[np.arange(n_ambig), top1_idx] = -np.inf
        top2_idx = np.argmax(tmp, axis=1)
        top2_scores = prim_mat_sub[np.arange(n_ambig), top2_idx]
    else:
        top2_idx = np.zeros(n_ambig, dtype=int)
        top2_scores = np.zeros(n_ambig, dtype=np.float64)

    prim_name_arr = np.array(prim_names, dtype=object)
    top1_labels = prim_name_arr[top1_idx]
    top2_labels = prim_name_arr[top2_idx]
    prim_margins = top1_scores - top2_scores

    for local_i, global_i in enumerate(ambig_idx):
        t1 = str(top1_labels[local_i])
        t2 = str(top2_labels[local_i])
        s1 = float(top1_scores[local_i])
        s2 = float(top2_scores[local_i])
        margin = float(prim_margins[local_i])

        if s1 >= args.min_single_score and margin > args.transition_margin:
            # Clear single-milestone winner.
            resolutions[int(global_i)] = dict(
                coarse=t1,
                expanded=_to_expanded(t1, c2e),
                confidence=s1,
                margin=margin,
                source="stage2_marker_single",
            )
        elif (
            _are_adjacent(t1, t2, adj)
            and s1 >= args.min_transition_score
            and s2 >= args.min_transition_score
            and abs(s1 - s2) <= args.transition_margin
        ):
            # Adjacent milestones with similar scores: transition zone.
            exp_trans = t2e.get(
                (t1, t2),
                t2e.get((t2, t1), f"{t1}_to_{t2}"),
            )
            resolutions[int(global_i)] = dict(
                coarse=t1,
                expanded=exp_trans,
                confidence=s1,
                margin=abs(s1 - s2),
                source="stage2_marker_transition",
            )
        # else: not resolved by marker evidence; leave for embedding or fallback.

    return resolutions


# ---------------------------------------------------------------------------
# Stage 2 resolution: independent embedding evidence
# ---------------------------------------------------------------------------

def _resolve_by_embedding(
    still_ambig_idx: object,
    adata: object,
    embedding_key: str,
    traj: dict,
    stage1_label_arr: object,
    stage1_status_arr: object,
    args: argparse.Namespace,
    warnings: list,
) -> dict:
    """Resolve still-ambiguous cells using an independent reference embedding.

    Centroid-based assignment.  Reference cells = Stage 1 high-confidence
    primary-milestone cells.  OOD threshold = per-milestone percentile of
    within-centroid reference distances.

    Returns a dict mapping global cell index -> resolution dict (same schema
    as _resolve_by_markers).  Only resolved cells are included.
    """
    import numpy as np

    resolutions: dict = {}
    if len(still_ambig_idx) == 0:
        return resolutions

    if embedding_key not in adata.obsm:
        warnings.append(
            f"--embedding-key '{embedding_key}' not found in adata.obsm "
            f"(available: {list(adata.obsm.keys())}).  "
            "Skipping embedding-based resolution."
        )
        return resolutions

    emb = np.asarray(adata.obsm[embedding_key], dtype=np.float64)  # (n_cells, d)

    primary = traj["primary_milestones"]
    adj = traj["adjacency_undirected"]
    c2e = traj["coarse_to_expanded"]
    t2e = traj["transition_to_expanded"]

    # Build per-milestone centroids from Stage 1 high-confidence reference cells.
    centroids: dict = {}
    ood_thresholds: dict = {}

    for ms in primary:
        ref_mask = (stage1_status_arr == "high_confidence") & (stage1_label_arr == ms)
        ref_idx = (ref_mask).nonzero()[0]
        if len(ref_idx) < args.min_reference_cells:
            warnings.append(
                f"Milestone '{ms}': {len(ref_idx)} Stage 1 high-confidence "
                f"reference cells (need >= {args.min_reference_cells}); "
                "centroid not built."
            )
            continue
        ref_emb = emb[ref_idx]                           # (n_ref, d)
        centroid = np.mean(ref_emb, axis=0)
        centroids[ms] = centroid

        # OOD threshold: percentile of within-centroid reference distances.
        dists_ref = np.sqrt(np.sum((ref_emb - centroid) ** 2, axis=1))
        ood_thresholds[ms] = float(
            np.percentile(dists_ref, args.ood_distance_quantile * 100.0)
        )

    if not centroids:
        warnings.append(
            "No milestone centroids built from Stage 1 reference cells; "
            "skipping embedding-based resolution."
        )
        return resolutions

    centroid_names = list(centroids.keys())
    centroid_mat = np.array(
        [centroids[ms] for ms in centroid_names], dtype=np.float64
    )  # (n_centroids, d)

    # Vectorised pairwise distances: ambiguous cells vs. all centroids.
    ambig_emb = emb[still_ambig_idx]                    # (n_ambig, d)
    diff = ambig_emb[:, None, :] - centroid_mat[None, :, :]  # (n_a, n_c, d)
    dist_mat = np.sqrt(np.sum(diff ** 2, axis=2))       # (n_ambig, n_centroids)

    n_centroids = len(centroid_names)

    for local_i, global_i in enumerate(still_ambig_idx):
        row = dist_mat[local_i]                          # (n_centroids,)
        order = np.argsort(row)

        nearest_ms = centroid_names[int(order[0])]
        nearest_d = float(row[order[0]])

        ood_thr = ood_thresholds.get(nearest_ms)
        if ood_thr is not None and nearest_d > ood_thr:
            resolutions[int(global_i)] = dict(
                coarse=args.unknown_label,
                expanded=args.unknown_label,
                confidence=0.0,
                margin=0.0,
                source="stage2_embedding_ood",
            )
            continue

        def _conf(ms: str, d: float) -> float:
            thr = ood_thresholds.get(ms)
            if thr and thr > 0:
                return max(0.0, 1.0 - d / thr)
            return 0.0

        if n_centroids >= 2:
            second_ms = centroid_names[int(order[1])]
            second_d = float(row[order[1]])
            total_d = nearest_d + second_d
            position = nearest_d / total_d if total_d > 0.0 else 0.5

            if _are_adjacent(nearest_ms, second_ms, adj):
                if position <= args.transition_position_low:
                    # Clearly at nearest milestone.
                    resolutions[int(global_i)] = dict(
                        coarse=nearest_ms,
                        expanded=_to_expanded(nearest_ms, c2e),
                        confidence=_conf(nearest_ms, nearest_d),
                        margin=second_d - nearest_d,
                        source="stage2_embedding_single",
                    )
                elif position >= args.transition_position_high:
                    # Clearly at second milestone.
                    resolutions[int(global_i)] = dict(
                        coarse=second_ms,
                        expanded=_to_expanded(second_ms, c2e),
                        confidence=_conf(second_ms, second_d),
                        margin=nearest_d - second_d,
                        source="stage2_embedding_single",
                    )
                else:
                    # Transition zone between adjacent milestones.
                    exp_trans = t2e.get(
                        (nearest_ms, second_ms),
                        t2e.get(
                            (second_ms, nearest_ms),
                            f"{nearest_ms}_to_{second_ms}",
                        ),
                    )
                    resolutions[int(global_i)] = dict(
                        coarse=nearest_ms,
                        expanded=exp_trans,
                        confidence=_conf(nearest_ms, nearest_d),
                        margin=abs(second_d - nearest_d),
                        source="stage2_embedding_transition",
                    )
            else:
                # Not adjacent: assign to nearest.
                resolutions[int(global_i)] = dict(
                    coarse=nearest_ms,
                    expanded=_to_expanded(nearest_ms, c2e),
                    confidence=_conf(nearest_ms, nearest_d),
                    margin=second_d - nearest_d,
                    source="stage2_embedding_single",
                )
        else:
            # Single centroid available.
            resolutions[int(global_i)] = dict(
                coarse=nearest_ms,
                expanded=_to_expanded(nearest_ms, c2e),
                confidence=_conf(nearest_ms, nearest_d),
                margin=0.0,
                source="stage2_embedding_single",
            )

    return resolutions


# ---------------------------------------------------------------------------
# QC output helpers
# ---------------------------------------------------------------------------

def _write_count_csv_pd(
    path: Path, counts: dict, label_col: str, tag: str, pd: object
) -> None:
    df = pd.DataFrame(
        list(counts.items()), columns=[label_col, "cell_count"]
    ).sort_values("cell_count", ascending=False)
    df.to_csv(str(path), index=False)
    print(f"{tag} {path.name:45s} -> {path}")


def _write_qc_csv(
    path: Path,
    adata: object,
    primary: list,
    args: argparse.Namespace,
    tag: str,
    pd: object,
) -> None:
    """Write per-cell QC CSV including Stage 1, Stage 2, and marker scores."""
    obs = adata.obs
    base_cols = [
        args.stage1_label_key,
        args.stage1_status_key,
        args.stage1_confidence_key,
        "stage2_source",
        "final_milestone_label_coarse",
        "final_milestone_label_expanded",
        "final_milestone_confidence",
        "final_milestone_source",
        "trajectory_membership_label",
        "trajectory_membership_score",
    ]
    score_cols = [
        f"milestone_score_{ms}"
        for ms in primary
        if f"milestone_score_{ms}" in obs.columns
    ]
    cols = [c for c in (base_cols + score_cols) if c in obs.columns]
    df = obs[cols].copy()
    df.index.name = "cell_id"
    df.to_csv(str(path))
    print(f"{tag} {path.name:45s} -> {path}")


# ---------------------------------------------------------------------------
# Core annotation pipeline
# ---------------------------------------------------------------------------

def _run_annotation(args: argparse.Namespace) -> dict:
    """Execute the full Stage 2 annotation pipeline.  Return metadata dict."""
    import numpy as np

    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError(
            "anndata is required.  Install: pip install anndata"
        ) from exc
    try:
        import pandas as pd
    except ImportError as exc:
        raise ImportError(
            "pandas is required.  Install: pip install pandas"
        ) from exc

    tag = "[build_trajectory_aware_labels]"
    warnings_list: list = []

    # ------------------------------------------------------------------
    # 1. Load AnnData
    # ------------------------------------------------------------------
    input_path = Path(args.input_h5ad)
    print(f"{tag} loading {input_path} ...")
    adata = ad.read_h5ad(str(input_path))
    print(f"{tag} loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # ------------------------------------------------------------------
    # 2. Load config and build trajectory structures
    # ------------------------------------------------------------------
    markers_yaml = Path(args.markers_yaml)
    ds_cfg = _load_config(markers_yaml, args.dataset_id)
    traj = _build_trajectory_structures(ds_cfg)

    primary = traj["primary_milestones"]
    primary_set = traj["primary_set"]
    c2e = traj["coarse_to_expanded"]

    print(f"{tag} primary_milestones:  {primary}")
    print(f"{tag} expanded_stages: {len(traj['expanded_trajectory_order'])}")

    # ------------------------------------------------------------------
    # 3. Validate Stage 1 input columns
    # ------------------------------------------------------------------
    s1_col_map = {
        "label":      args.stage1_label_key,
        "status":     args.stage1_status_key,
        "confidence": args.stage1_confidence_key,
        "margin":     args.stage1_margin_key,
    }
    missing_cols = [v for v in s1_col_map.values() if v not in adata.obs.columns]
    if missing_cols:
        raise ValueError(
            f"Required Stage 1 obs columns missing from input h5ad: "
            f"{missing_cols}.  Run build_marker_seed_labels.py first."
        )

    stage1_label_arr = adata.obs[s1_col_map["label"]].to_numpy(dtype=str)
    stage1_status_arr = adata.obs[s1_col_map["status"]].to_numpy(dtype=str)
    stage1_conf_arr = np.nan_to_num(
        adata.obs[s1_col_map["confidence"]].to_numpy(dtype=np.float64), nan=0.0
    )
    stage1_margin_arr = np.nan_to_num(
        adata.obs[s1_col_map["margin"]].to_numpy(dtype=np.float64), nan=0.0
    )

    # Prefer existing stage1_source column; fall back to derived labels.
    if "stage1_source" in adata.obs.columns:
        stage1_src_arr = adata.obs["stage1_source"].to_numpy(dtype=str)
    else:
        warnings_list.append(
            "'stage1_source' column not found; using derived source labels."
        )
        _src_map = {
            "high_confidence": "stage1_high_confidence",
            "ambiguous": "stage1_ambiguous",
        }
        stage1_src_arr = np.array(
            [_src_map.get(s, "stage1_unknown") for s in stage1_status_arr],
            dtype=object,
        )

    n_cells = adata.n_obs

    # ------------------------------------------------------------------
    # 4. Initialise output arrays (default = ambiguous everywhere)
    # ------------------------------------------------------------------
    def _obj(fill: str) -> object:
        import numpy as _np
        return _np.full(n_cells, fill, dtype=object)

    def _f64() -> object:
        import numpy as _np
        return _np.zeros(n_cells, dtype=np.float64)

    final_coarse   = _obj(args.ambiguous_label)
    final_expanded = _obj(args.ambiguous_label)
    final_conf     = _f64()
    final_source   = _obj("unset")

    stage2_coarse   = _obj(args.ambiguous_label)
    stage2_expanded = _obj(args.ambiguous_label)
    stage2_conf     = _f64()
    stage2_source   = _obj("not_applicable_stage1_final")

    traj_label  = _obj("ambiguous")
    traj_score  = _f64()
    traj_margin = _f64()

    # ------------------------------------------------------------------
    # 5. Initialise from Stage 1: high-confidence cells
    # ------------------------------------------------------------------
    hc_mask = stage1_status_arr == "high_confidence"
    ambig_mask = ~hc_mask  # ambiguous + any other unrecognised

    n_hc = int(hc_mask.sum())
    n_ambig_s1 = int(ambig_mask.sum())

    print(
        f"{tag} Stage 1 breakdown: "
        f"high_confidence={n_hc}, "
        f"ambiguous={n_ambig_s1}"
    )

    # High-confidence cells.
    if hc_mask.any():
        hc_lbl = stage1_label_arr[hc_mask]
        hc_exp = np.array(
            [_to_expanded(lbl, c2e) for lbl in hc_lbl], dtype=object
        )
        final_coarse[hc_mask]   = hc_lbl
        final_expanded[hc_mask] = hc_exp
        final_conf[hc_mask]     = stage1_conf_arr[hc_mask]
        final_source[hc_mask]   = stage1_src_arr[hc_mask]

        stage2_coarse[hc_mask]   = hc_lbl
        stage2_expanded[hc_mask] = hc_exp
        stage2_conf[hc_mask]     = stage1_conf_arr[hc_mask]
        stage2_source[hc_mask]   = "not_applicable_stage1_final"

        hc_memb = np.array(
            [
                _membership_label(
                    lbl, primary_set,
                    args.unknown_label, args.ambiguous_label,
                )
                for lbl in hc_lbl
            ],
            dtype=object,
        )
        traj_label[hc_mask]  = hc_memb
        traj_score[hc_mask]  = stage1_conf_arr[hc_mask]
        traj_margin[hc_mask] = stage1_margin_arr[hc_mask]

    # ------------------------------------------------------------------
    # 6. Stage 2 marker-score resolution for ambiguous cells
    # ------------------------------------------------------------------
    ambig_idx = np.where(ambig_mask)[0]
    marker_res: dict = {}

    if n_ambig_s1 > 0:
        print(
            f"{tag} Stage 2 marker resolution for "
            f"{n_ambig_s1} ambiguous cells ..."
        )
        marker_res = _resolve_by_markers(
            ambig_idx=ambig_idx,
            obs=adata.obs,
            traj=traj,
            args=args,
            warnings=warnings_list,
        )
        print(
            f"{tag}   marker resolved: "
            f"{len(marker_res)}/{n_ambig_s1}"
        )

    for gi, res in marker_res.items():
        _apply_resolution(
            gi, res,
            final_coarse, final_expanded, final_conf, final_source,
            stage2_coarse, stage2_expanded, stage2_conf, stage2_source,
            traj_label, traj_score, traj_margin,
            primary_set, args,
        )

    # ------------------------------------------------------------------
    # 7. Optional embedding-based resolution for still-ambiguous cells
    # ------------------------------------------------------------------
    still_ambig_idx = np.array(
        [int(i) for i in ambig_idx if i not in marker_res], dtype=int
    )
    emb_res: dict = {}

    if args.embedding_key and len(still_ambig_idx) > 0:
        print(
            f"{tag} Stage 2 embedding resolution "
            f"(key='{args.embedding_key}') for "
            f"{len(still_ambig_idx)} still-ambiguous cells ..."
        )
        emb_res = _resolve_by_embedding(
            still_ambig_idx=still_ambig_idx,
            adata=adata,
            embedding_key=args.embedding_key,
            traj=traj,
            stage1_label_arr=stage1_label_arr,
            stage1_status_arr=stage1_status_arr,
            args=args,
            warnings=warnings_list,
        )
        print(
            f"{tag}   embedding resolved: "
            f"{len(emb_res)}/{len(still_ambig_idx)}"
        )

        for gi, res in emb_res.items():
            _apply_resolution(
                gi, res,
                final_coarse, final_expanded, final_conf, final_source,
                stage2_coarse, stage2_expanded, stage2_conf, stage2_source,
                traj_label, traj_score, traj_margin,
                primary_set, args,
            )

    # ------------------------------------------------------------------
    # 8. Final fallback: unresolved ambiguous cells
    # ------------------------------------------------------------------
    all_resolved = set(marker_res.keys()) | set(emb_res.keys())
    unresolved_idx = np.array(
        [int(i) for i in ambig_idx if i not in all_resolved], dtype=int
    )
    n_unresolved = len(unresolved_idx)

    if n_unresolved > 0:
        # Use the highest available primary score as best-effort confidence.
        best_scores = np.zeros(n_cells, dtype=np.float64)
        for ms in primary:
            col = f"milestone_score_{ms}"
            if col in adata.obs.columns:
                arr = np.nan_to_num(
                    adata.obs[col].to_numpy(dtype=np.float64), nan=0.0
                )
                best_scores = np.maximum(best_scores, arr)

        final_coarse[unresolved_idx]   = args.ambiguous_label
        final_expanded[unresolved_idx] = args.ambiguous_label
        final_conf[unresolved_idx]     = best_scores[unresolved_idx]
        final_source[unresolved_idx]   = "stage2_unresolved_ambiguous"

        stage2_coarse[unresolved_idx]   = args.ambiguous_label
        stage2_expanded[unresolved_idx] = args.ambiguous_label
        stage2_conf[unresolved_idx]     = best_scores[unresolved_idx]
        stage2_source[unresolved_idx]   = "stage2_unresolved_ambiguous"

        traj_label[unresolved_idx]  = "ambiguous"
        traj_score[unresolved_idx]  = best_scores[unresolved_idx]
        traj_margin[unresolved_idx] = 0.0

    # ------------------------------------------------------------------
    # 9. Write obs columns
    # ------------------------------------------------------------------
    def _safe_f32(arr: object) -> object:
        return np.clip(
            np.nan_to_num(np.asarray(arr, dtype=np.float32), nan=0.0),
            -1e6, 1e6,
        )

    adata.obs["trajectory_membership_label"]    = traj_label.astype(str)
    adata.obs["trajectory_membership_score"]    = _safe_f32(traj_score)
    adata.obs["trajectory_membership_margin"]   = _safe_f32(traj_margin)
    adata.obs["stage2_label_coarse"]            = stage2_coarse.astype(str)
    adata.obs["stage2_label_expanded"]          = stage2_expanded.astype(str)
    adata.obs["stage2_confidence"]              = _safe_f32(stage2_conf)
    adata.obs["stage2_source"]                  = stage2_source.astype(str)
    adata.obs["final_milestone_label_coarse"]   = final_coarse.astype(str)
    adata.obs["final_milestone_label_expanded"] = final_expanded.astype(str)
    adata.obs["final_milestone_confidence"]     = _safe_f32(final_conf)
    adata.obs["final_milestone_source"]         = final_source.astype(str)

    # ------------------------------------------------------------------
    # 10. Count summaries
    # ------------------------------------------------------------------
    def _vc(col: str) -> dict:
        return adata.obs[col].value_counts().to_dict()

    final_coarse_counts   = _vc("final_milestone_label_coarse")
    final_expanded_counts = _vc("final_milestone_label_expanded")
    final_source_counts   = _vc("final_milestone_source")
    traj_counts           = _vc("trajectory_membership_label")
    s2_source_counts      = _vc("stage2_source")

    print(f"{tag} final_milestone_label_coarse: {final_coarse_counts}")
    print(f"{tag} trajectory_membership_label:  {traj_counts}")

    def _frac(n: int) -> float:
        return round(n / n_cells, 4) if n_cells > 0 else 0.0

    n_ms  = int(s2_source_counts.get("stage2_marker_single",       0))
    n_mt  = int(s2_source_counts.get("stage2_marker_transition",   0))
    n_es  = int(s2_source_counts.get("stage2_embedding_single",    0))
    n_et  = int(s2_source_counts.get("stage2_embedding_transition",0))
    n_ur  = int(s2_source_counts.get("stage2_unresolved_ambiguous",0))
    n_ood = int(final_coarse_counts.get(args.unknown_label, 0))
    n_cov = sum(v for k, v in final_coarse_counts.items() if k in primary_set)

    # ------------------------------------------------------------------
    # 11. uns metadata
    # ------------------------------------------------------------------
    uns_meta = {
        "script": "build_trajectory_aware_labels",
        "status": "dry_run" if args.dry_run else "completed",
        "dataset_id": args.dataset_id,
        "input_h5ad": str(args.input_h5ad),
        "output_h5ad": str(args.output_h5ad),
        "markers_yaml": str(markers_yaml),
        "embedding_key": args.embedding_key,
        "primary_milestones": primary,
        "expanded_trajectory_order": traj["expanded_trajectory_order"],
        "excluded_from_official_metrics": traj["excluded_from_official_metrics"],
        "thresholds": {
            "min_single_score":          args.min_single_score,
            "min_transition_score":      args.min_transition_score,
            "transition_margin":         args.transition_margin,
            "ood_distance_quantile":     args.ood_distance_quantile,
            "min_reference_cells":       args.min_reference_cells,
            "transition_position_low":   args.transition_position_low,
            "transition_position_high":  args.transition_position_high,
        },
        "stage1_input_counts": {
            "high_confidence":           n_hc,
            "ambiguous":                 n_ambig_s1,
        },
        "final_coarse_label_counts":   {str(k): int(v) for k, v in final_coarse_counts.items()},
        "final_expanded_label_counts": {str(k): int(v) for k, v in final_expanded_counts.items()},
        "final_source_counts":         {str(k): int(v) for k, v in final_source_counts.items()},
        "trajectory_membership_counts":{str(k): int(v) for k, v in traj_counts.items()},
        "coverage_fraction_primary_milestones": _frac(n_cov),
        "resolution_summary": {
            "n_marker_single":            n_ms,
            "frac_marker_single":         _frac(n_ms),
            "n_marker_transition":        n_mt,
            "frac_marker_transition":     _frac(n_mt),
            "n_embedding_single":         n_es,
            "frac_embedding_single":      _frac(n_es),
            "n_embedding_transition":     n_et,
            "frac_embedding_transition":  _frac(n_et),
            "n_unresolved_ambiguous":     n_ur,
            "frac_unresolved_ambiguous":  _frac(n_ur),
            "n_unknown_or_ood":           n_ood,
            "frac_unknown_or_ood":        _frac(n_ood),
        },
        "warnings": warnings_list,
        "n_cells": n_cells,
    }
    adata.uns["trajectory_aware_annotation_stage2"] = uns_meta

    # ------------------------------------------------------------------
    # 12. Write outputs
    # ------------------------------------------------------------------
    out_h5ad = Path(args.output_h5ad)
    out_dir  = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print(f"{tag} dry-run: skipping h5ad write -> {out_h5ad}")
    else:
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)
        print(f"{tag} writing annotated h5ad -> {out_h5ad}")
        adata.write_h5ad(str(out_h5ad))

    # Metadata JSON always written (even in dry-run; status field reflects it).
    meta_path = out_dir / "trajectory_aware_labels_metadata.json"
    with meta_path.open("w", encoding="ascii") as fh:
        json.dump(uns_meta, fh, indent=2, default=_json_default)
    print(f"{tag} metadata -> {meta_path}")

    # Count CSVs.
    _write_count_csv_pd(out_dir / "final_coarse_label_counts.csv",
                        final_coarse_counts,   "final_milestone_label_coarse",   tag, pd)
    _write_count_csv_pd(out_dir / "final_expanded_label_counts.csv",
                        final_expanded_counts, "final_milestone_label_expanded", tag, pd)
    _write_count_csv_pd(out_dir / "final_source_counts.csv",
                        final_source_counts,   "final_milestone_source",         tag, pd)
    _write_count_csv_pd(out_dir / "trajectory_membership_counts.csv",
                        traj_counts,           "trajectory_membership_label",    tag, pd)

    # Detailed per-cell QC CSV.
    _write_qc_csv(
        path=out_dir / "trajectory_aware_annotation_qc.csv",
        adata=adata,
        primary=primary,
        args=args,
        tag=tag,
        pd=pd,
    )

    if warnings_list:
        print(f"{tag} {len(warnings_list)} warning(s):")
        for w in warnings_list:
            print(f"  WARN: {w}")
    else:
        print(f"{tag} no warnings.")

    return uns_meta


# ---------------------------------------------------------------------------
# JSON serialisation helper
# ---------------------------------------------------------------------------

def _json_default(obj: object) -> object:
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

    tag = "[build_trajectory_aware_labels]"
    if args.dry_run:
        print(f"{tag} dry-run mode active: h5ad will not be written.")

    try:
        _run_annotation(args)
    except Exception as exc:  # noqa: BLE001
        print(f"{tag} ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    print(f"{tag} done.")
    sys.exit(0)


if __name__ == "__main__":
    main()
