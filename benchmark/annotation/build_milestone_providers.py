"""build_milestone_providers.py
================================
Build milestone-based ground-truth provider directories and registry entries.

Updated for Stage 2 two-stage annotation layer.  The official provider is now
built from frozen silver-standard Stage 2 labels:
  final_milestone_label_expanded   (official state key)
  final_milestone_label_coarse     (coarse compatibility key)
  final_milestone_confidence
  final_milestone_source
  trajectory_membership_label

Official providers (recommended, default)
-----------------------------------------
  gse178325_marker_fm_transition_silver_v1/
  gse230659_marker_fm_transition_silver_v1/

Each provider directory contains:
  state_labels.tsv              cell-level labels (schema-only unless --input-h5ad)
  state_metadata.tsv            one row per official graph state
  annotation_votes.tsv          per-cell annotation evidence
  reference_graph.json          milestone reference graph (nodes + edges)
  reference_graph_edges.csv     edge list for quick inspection
  ground_truth_metadata.json    full provider metadata

Label modes (--label-mode)
--------------------------
  official_silver   frozen silver-standard provider from Stage 2 (DEFAULT)

Usage
-----
  python benchmark/annotation/build_milestone_providers.py --dry-run
  python benchmark/annotation/build_milestone_providers.py
  python benchmark/annotation/build_milestone_providers.py \
      --input-h5ad path/to/GSE230659_stage2.h5ad \
      --dataset-id GSE230659 \
      --label-mode official_silver

Provenance
----------
  Liuyang et al. 2023, Cell Stem Cell, DOI: 10.1016/j.stem.2023.02.008
  Guan et al. 2022, Nature, DOI: 10.1038/s41586-022-04593-5
  Companion code: https://github.com/sajuukLyu/CSC_2023

Status: Step 4 -- milestone provider builder, Stage 2 official_silver mode.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))


_VALID_DATASET_IDS = ("GSE178325", "GSE230659")

# Individual mode names
_VALID_LABEL_MODES = (
    "official_silver",
)

_VALID_LABEL_MODE_ARGS = (
    "official_silver",
)

_DEFAULT_MARKERS_YAML = "benchmark/annotation/milestone_markers.yaml"
_DEFAULT_REGISTRY_YAML = "benchmark/ground_truth/registry.yaml"
_DEFAULT_PROVIDER_ROOT = "benchmark/ground_truth/providers"

# Per-mode configuration
_MODE_CONFIG = {
    # ---- Official frozen silver-standard provider (recommended) ----
    "official_silver": {
        "state_key":          "final_milestone_label_expanded",
        "coarse_state_key":   "final_milestone_label_coarse",
        "expanded_state_key": "final_milestone_label_expanded",
        "confidence_key":     "final_milestone_confidence",
        "source_key":         "final_milestone_source",
        "membership_key":     "trajectory_membership_label",
        "analysis_role":      "primary_report",
        "provider_suffix":    "marker_fm_transition_silver_v1",
        "annotation_method":  "marker_fm_transition_silver_provider",
        "status":             "frozen_silver_standard_milestone_provider",
        "label_type":         "frozen_silver_standard",
    },
}

# Default excluded labels (overridden by excluded_from_official_metrics in YAML)
_EXCLUDED_FROM_OFFICIAL_METRICS_DEFAULT = [
    "ambiguous",
    "unknown_or_ood",
]

_ANNOTATION_VOTES_COLS_OFFICIAL_SILVER = [
    "cell_id",
    # Stage 1 columns
    "stage1_label",
    "stage1_status",
    "stage1_confidence",
    "stage1_margin",
    "stage1_source",
    # Stage 2 columns
    "stage2_label_coarse",
    "stage2_label_expanded",
    "stage2_confidence",
    "stage2_source",
    # Final frozen labels
    "final_milestone_label_coarse",
    "final_milestone_label_expanded",
    "final_milestone_confidence",
    "final_milestone_source",
    # Trajectory membership
    "trajectory_membership_label",
    "trajectory_membership_score",
    "trajectory_membership_margin",
    # Marker diagnostics.
    "milestone_marker_label",
    "milestone_marker_score",
    # Provider tracking
    "provider_id",
    "label_mode",
]


# ---------------------------------------------------------------------------
# Mode resolution helper
# ---------------------------------------------------------------------------

def _resolve_modes(label_mode_arg: str | None) -> list:
    """Resolve --label-mode. Only official_silver is supported."""
    if label_mode_arg in (None, "official_silver"):
        return ["official_silver"]
    raise ValueError(f"Unsupported label mode: {label_mode_arg!r}")


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _dataset_prefix(dataset_id: str) -> str:
    return dataset_id.lower()


def _provider_id(dataset_id: str, label_mode: str) -> str:
    prefix = _dataset_prefix(dataset_id)
    suffix = _MODE_CONFIG[label_mode]["provider_suffix"]
    return f"{prefix}_{suffix}"


def _path_text(path: str | Path | None) -> str | None:
    if path is None:
        return None
    return Path(path).as_posix()


def _write_tsv(path: Path, rows: list, fieldnames: list) -> None:
    with path.open("w", newline="", encoding="ascii") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=fieldnames, extrasaction="ignore",
            delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_csv(path: Path, rows: list, fieldnames: list) -> None:
    with path.open("w", newline="", encoding="ascii") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=fieldnames, extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


# ---------------------------------------------------------------------------
# Reference graph builder
# ---------------------------------------------------------------------------

def _stage_payload(stage_id: str) -> str:
    parts = str(stage_id).split("_", 2)
    if len(parts) == 3 and parts[0] == "stage" and parts[1].isdigit():
        return parts[2]
    return str(stage_id)


def _stage_marker_names(stage_id: str) -> list[str]:
    payload = _stage_payload(stage_id)
    if "_to_" in payload:
        return [part for part in payload.split("_to_") if part]
    return [payload]


def _is_transition_stage(stage_id: str) -> bool:
    return "_to_" in _stage_payload(stage_id)


def _build_reference_graph(
    provider_id: str,
    dataset_id: str,
    label_mode: str,
    primary_milestones: list,
    raw_edges: list,
    markers_yaml_path: str,
    ds_cfg: dict | None = None,
) -> dict:
    """Build the reference_graph dict (nodes + edges + _meta).

    For official_silver, _meta includes expanded_trajectory_order,
    excluded_from_official_metrics, and a frozen-silver graph_type.
    """
    mode_cfg = _MODE_CONFIG[label_mode]
    state_key = mode_cfg["state_key"]
    analysis_role = mode_cfg["analysis_role"]

    expanded_order = (ds_cfg or {}).get("expanded_trajectory_order") or []
    use_expanded_graph = bool(expanded_order) and (
        state_key == mode_cfg.get("expanded_state_key")
    )
    graph_states = list(expanded_order) if use_expanded_graph else list(primary_milestones)

    nodes = [
        {
            "id": state_id,
            "label": state_id,
            "status": "confirmed",
            "role": (
                "expanded_transition_stage"
                if use_expanded_graph and _is_transition_stage(state_id)
                else "expanded_primary_stage" if use_expanded_graph
                else "primary_milestone"
            ),
            "order": i,
        }
        for i, state_id in enumerate(graph_states)
    ]

    if label_mode != "official_silver":
        raise ValueError(f"Unsupported label_mode after cleanup: {label_mode}")

    if use_expanded_graph:
        edge_pairs = list(zip(graph_states[:-1], graph_states[1:]))
        edge_type = "marker_defined_expanded_stage_order"
        graph_type = "marker_defined_expanded_transition_graph"
        policy_head = "Official frozen silver-standard expanded transition reference graph."
    else:
        edge_pairs = [
            (str(e[0]), str(e[1]))
            for e in raw_edges
            if isinstance(e, (list, tuple)) and len(e) == 2
        ]
        edge_type = "marker_defined_milestone_order"
        graph_type = "marker_defined_coarse_milestone_graph"
        policy_head = "Official frozen silver-standard coarse milestone reference graph."

    edges = [
        {
            "source": str(src),
            "target": str(tgt),
            "weight": 1.0,
            "confidence": "high",
            "source_status": "confirmed",
            "target_status": "confirmed",
            "edge_type": edge_type,
        }
        for src, tgt in edge_pairs
    ]

    excluded = (
        (ds_cfg or {}).get("excluded_from_official_metrics")
        or list(_EXCLUDED_FROM_OFFICIAL_METRICS_DEFAULT)
    )
    graph_meta = {
        "provider_id":              provider_id,
        "dataset_id":               dataset_id,
        "label_mode":               label_mode,
        "label_type":               mode_cfg["label_type"],
        "state_key":                state_key,
        "expanded_state_key":       mode_cfg.get("expanded_state_key"),
        "graph_type":               graph_type,
        "version":                  "v1",
        "generated_by":             "benchmark/annotation/build_milestone_providers.py",
        "source_markers_yaml":      Path(markers_yaml_path).as_posix(),
        "analysis_role":            analysis_role,
        "deprecated":               False,
        "n_states":                 len(graph_states),
        "n_edges":                  len(edges),
        "official_state_ids":       graph_states,
        "primary_milestones":       list(primary_milestones),
        "expanded_trajectory_order":    expanded_order,
        "excluded_from_official_metrics": excluded,
        "policy": [
            policy_head,
            (
                "Nodes are expanded_trajectory_order stages, including transition "
                "labels, because official marker-FM metrics use "
                "final_milestone_label_expanded."
                if use_expanded_graph else
                "Nodes are primary_milestones only."
            ),
            "All edges have confidence=high; all nodes have status=confirmed.",
            "Excluded from official metrics: ambiguous, unknown_or_ood.",
            "Provenance: Liuyang et al. 2023 DOI 10.1016/j.stem.2023.02.008;"
            " Guan et al. 2022 DOI 10.1038/s41586-022-04593-5;"
            " companion code https://github.com/sajuukLyu/CSC_2023",
        ],
    }

    return {"_meta": graph_meta, "nodes": nodes, "edges": edges}


# ---------------------------------------------------------------------------
# State metadata builder
# ---------------------------------------------------------------------------

def _build_state_metadata_rows(state_ids: list, marker_sets: dict) -> list:
    rows = []
    for i, state_id in enumerate(state_ids):
        marker_names = _stage_marker_names(state_id)
        genes: list[str] = []
        notes_parts: list[str] = []
        for marker_name in marker_names:
            ms_cfg = marker_sets.get(marker_name) or {}
            for gene in ms_cfg.get("genes") or []:
                if gene not in genes:
                    genes.append(gene)
            note = ms_cfg.get("notes") or ""
            if note:
                notes_parts.append(f"{marker_name}: {note}")
        role = (
            "expanded_transition_stage"
            if _is_transition_stage(state_id)
            else "expanded_primary_stage" if str(state_id).startswith("stage_")
            else "primary_milestone"
        )
        notes = " | ".join(notes_parts)
        rows.append({
            "state_id":       state_id,
            "label":          state_id,
            "order":          i,
            "role":           role,
            "status":         "confirmed",
            "marker_genes":   ";".join(str(g) for g in genes),
            "n_marker_genes": len(genes),
            "notes":          str(notes),
        })
    return rows


# ---------------------------------------------------------------------------
# h5ad label reader (lazy anndata import)
# ---------------------------------------------------------------------------

def _read_h5ad_labels(
    input_h5ad: str,
    state_key: str,
    confidence_key: str,
    provider_id: str,
    label_mode: str,
    vote_cols: list,
    excluded_labels: list,
    tag: str,
) -> tuple:
    """Read per-cell labels from a Stage 2 h5ad file.

    Returns
    -------
    state_labels_rows    : list[dict]  (one row per cell)
    annotation_votes_rows: list[dict]  (one row per cell)
    state_labels_status  : str
    label_stats          : dict        (counts, coverage; populated for official_silver)
    """
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError(
            "anndata is required when --input-h5ad is provided.  "
            "Install: pip install anndata"
        ) from exc

    print(f"{tag}   loading h5ad for labels: {input_h5ad}")
    adata = ad.read_h5ad(str(input_h5ad), backed="r")

    if state_key not in adata.obs.columns:
        raise ValueError(
            f"h5ad obs does not contain required column '{state_key}'.  "
            f"Available columns: {list(adata.obs.columns)}"
        )

    cell_ids = (
        list(adata.obs["cell_id"].astype(str))
        if "cell_id" in adata.obs.columns
        else list(adata.obs_names)
    )
    state_ids  = list(adata.obs[state_key].astype(str))
    confidences = (
        list(adata.obs[confidence_key].astype(str))
        if confidence_key in adata.obs.columns
        else [""] * len(cell_ids)
    )

    state_labels_rows = [
        {
            "cell_id":    cid,
            "state_id":   sid,
            "confidence": conf,
            "label_mode": label_mode,
            "provider_id": provider_id,
        }
        for cid, sid, conf in zip(cell_ids, state_ids, confidences)
    ]

    # annotation_votes: use only columns that exist in obs
    avail_vote_cols = [
        c for c in vote_cols
        if c not in ("cell_id", "provider_id", "label_mode")
        and c in adata.obs.columns
    ]
    annotation_votes_rows = []
    for i, cid in enumerate(cell_ids):
        row = {"cell_id": cid, "provider_id": provider_id, "label_mode": label_mode}
        for col in avail_vote_cols:
            row[col] = str(adata.obs[col].iloc[i])
        annotation_votes_rows.append(row)

    print(f"{tag}   read {len(state_labels_rows)} cell labels from h5ad.")

    # Compute label statistics for official_silver
    label_stats: dict = {}
    if label_mode == "official_silver":
        def _vc(col: str) -> dict:
            if col in adata.obs.columns:
                return {
                    str(k): int(v)
                    for k, v in adata.obs[col].value_counts().items()
                }
            return {}

        coarse_counts   = _vc("final_milestone_label_coarse")
        expanded_counts = _vc("final_milestone_label_expanded")
        source_counts   = _vc("final_milestone_source")
        memb_counts     = _vc("trajectory_membership_label")

        n_total  = len(cell_ids)
        state_counts = _vc(state_key)
        n_covered = sum(
            v for k, v in state_counts.items()
            if k not in (excluded_labels or [])
        )
        coverage = round(n_covered / n_total, 4) if n_total > 0 else 0.0

        label_stats = {
            "coarse_label_counts":    coarse_counts,
            "expanded_label_counts":  expanded_counts,
            "state_label_counts":     state_counts,
            "source_counts":          source_counts,
            "membership_counts":      memb_counts,
            "label_coverage_fraction": coverage,
            "n_cell_labels":          n_total,
        }

    return state_labels_rows, annotation_votes_rows, "populated_from_input_h5ad", label_stats


# ---------------------------------------------------------------------------
# Single-provider builder
# ---------------------------------------------------------------------------

def _build_one_provider(
    dataset_id: str,
    label_mode: str,
    ds_cfg: dict,
    provider_root: Path,
    markers_yaml_path: Path,
    input_h5ad: str | None,
    dry_run: bool,
    tag: str,
) -> dict:
    """Build one provider directory and return its metadata dict."""
    mode_cfg     = _MODE_CONFIG[label_mode]
    state_key    = mode_cfg["state_key"]
    confidence_key = mode_cfg["confidence_key"]
    analysis_role  = mode_cfg["analysis_role"]
    pid = _provider_id(dataset_id, label_mode)

    primary_milestones = list(ds_cfg.get("primary_milestones") or [])
    marker_sets  = ds_cfg.get("marker_sets") or {}
    raw_edges    = ds_cfg.get("reference_graph_edges") or []
    expanded_order = list(ds_cfg.get("expanded_trajectory_order") or [])
    excluded       = list(
        ds_cfg.get("excluded_from_official_metrics")
        or _EXCLUDED_FROM_OFFICIAL_METRICS_DEFAULT
    )

    provider_dir = provider_root / pid

    vote_cols = _ANNOTATION_VOTES_COLS_OFFICIAL_SILVER

    # Build graph
    graph = _build_reference_graph(
        provider_id=pid,
        dataset_id=dataset_id,
        label_mode=label_mode,
        primary_milestones=primary_milestones,
        raw_edges=raw_edges,
        markers_yaml_path=str(markers_yaml_path),
        ds_cfg=ds_cfg,
    )
    official_state_ids = [node["id"] for node in graph["nodes"]]
    n_states      = len(official_state_ids)
    n_graph_edges = len(graph["edges"])

    # Build state metadata rows
    state_meta_rows = _build_state_metadata_rows(official_state_ids, marker_sets)

    # Build edge rows for CSV
    edge_rows = [
        {
            "source":     e["source"],
            "target":     e["target"],
            "weight":     e["weight"],
            "confidence": e["confidence"],
            "edge_type":  e["edge_type"],
        }
        for e in graph["edges"]
    ]

    # Labels from h5ad or schema-only
    state_labels_rows: list     = []
    annotation_votes_rows: list = []
    state_labels_status = "schema_only_no_input_h5ad"
    label_stats: dict = {}

    if input_h5ad is not None:
        (
            state_labels_rows,
            annotation_votes_rows,
            state_labels_status,
            label_stats,
        ) = _read_h5ad_labels(
            input_h5ad=input_h5ad,
            state_key=state_key,
            confidence_key=confidence_key,
            provider_id=pid,
            label_mode=label_mode,
            vote_cols=vote_cols,
            excluded_labels=excluded,
            tag=tag,
        )

    def _rel(fname: str) -> str:
        return (provider_dir / fname).as_posix()

    metadata: dict = {
        "provider_id":            pid,
        "dataset_id":             dataset_id,
        "label_mode":             label_mode,
        "label_type":             mode_cfg["label_type"],
        "annotation_method":      mode_cfg["annotation_method"],
        "status":                 mode_cfg["status"],
        "state_key":              state_key,
        "coarse_state_key":       mode_cfg.get("coarse_state_key"),
        "expanded_state_key":     mode_cfg.get("expanded_state_key"),
        "confidence_key":         confidence_key,
        "source_key":             mode_cfg.get("source_key"),
        "membership_key":         mode_cfg.get("membership_key"),
        "analysis_role":          analysis_role,
        "version":                "v1",
        "source_markers_yaml":    markers_yaml_path.as_posix(),
        "excluded_from_official_metrics":  excluded,
        "expanded_trajectory_order":       expanded_order,
        "primary_milestones":              primary_milestones,
        "official_state_ids":              official_state_ids,
        "n_states":               n_states,
        "n_graph_edges":          n_graph_edges,
        "n_cell_labels":          label_stats.get("n_cell_labels", 0),
        "coarse_label_counts":    label_stats.get("coarse_label_counts", {}),
        "expanded_label_counts":  label_stats.get("expanded_label_counts", {}),
        "state_label_counts":     label_stats.get("state_label_counts", {}),
        "source_counts":          label_stats.get("source_counts", {}),
        "membership_counts":      label_stats.get("membership_counts", {}),
        "label_coverage_fraction": label_stats.get("label_coverage_fraction", 0.0),
        "input_h5ad":             _path_text(input_h5ad),
        "state_labels_status":    state_labels_status,
        "reference_graph_path":   _rel("reference_graph.json"),
        "reference_edges_path":   _rel("reference_graph_edges.csv"),
        "label_path":             _rel("state_labels.tsv"),
        "metadata_path":          _rel("state_metadata.tsv"),
        "annotation_votes_path":  _rel("annotation_votes.tsv"),
        "generated_by":           "benchmark/annotation/build_milestone_providers.py",
        "provenance": [
            "Liuyang et al. 2023, Cell Stem Cell,"
            " DOI: 10.1016/j.stem.2023.02.008",
            "Guan et al. 2022, Nature,"
            " DOI: 10.1038/s41586-022-04593-5",
            "Companion code: https://github.com/sajuukLyu/CSC_2023",
        ],
    }

    if dry_run:
        print(f"{tag}   [dry-run] {pid}")
        print(f"{tag}     dir={provider_dir}")
        print(f"{tag}     state_key={state_key}")
        print(f"{tag}     n_states={n_states}, n_edges={n_graph_edges}")
        print(f"{tag}     state_labels_status={state_labels_status}")
        if label_mode == "official_silver":
            print(f"{tag}     expanded_trajectory_order ({len(expanded_order)} stages)")
            print(f"{tag}     excluded_from_official_metrics={excluded}")
        return metadata

    # ---- Write files ----
    provider_dir.mkdir(parents=True, exist_ok=True)

    with (provider_dir / "reference_graph.json").open("w", encoding="ascii") as fh:
        json.dump(graph, fh, indent=2)

    _write_csv(
        provider_dir / "reference_graph_edges.csv",
        edge_rows,
        ["source", "target", "weight", "confidence", "edge_type"],
    )
    _write_tsv(
        provider_dir / "state_metadata.tsv",
        state_meta_rows,
        ["state_id", "label", "order", "role", "status",
         "marker_genes", "n_marker_genes", "notes"],
    )
    _write_tsv(
        provider_dir / "state_labels.tsv",
        state_labels_rows,
        ["cell_id", "state_id", "confidence", "label_mode", "provider_id"],
    )
    _write_tsv(
        provider_dir / "annotation_votes.tsv",
        annotation_votes_rows,
        vote_cols,
    )
    with (provider_dir / "ground_truth_metadata.json").open("w", encoding="ascii") as fh:
        json.dump(metadata, fh, indent=2)

    print(f"{tag}   created: {provider_dir}")
    return metadata


# ---------------------------------------------------------------------------
# Registry updater
# ---------------------------------------------------------------------------

def _update_registry(
    registry_path: Path,
    output_registry_path: Path | None,
    provider_metadatas: list,
    dry_run: bool,
    tag: str,
) -> None:
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("PyYAML required: pip install pyyaml") from exc

    with registry_path.open(encoding="utf-8-sig") as fh:
        registry = yaml.safe_load(fh) or {}
    if "providers" not in registry:
        registry["providers"] = {}

    for meta in provider_metadatas:
        pid = meta["provider_id"]
        lm  = meta["label_mode"]

        if lm != "official_silver":
            raise ValueError(f"Unsupported label_mode in registry update: {lm}")

        excluded = meta.get(
            "excluded_from_official_metrics",
            list(_EXCLUDED_FROM_OFFICIAL_METRICS_DEFAULT),
        )
        entry = {
            "provider_id":       pid,
            "annotation_method": meta["annotation_method"],
            "status":            meta["status"],
            "state_key":         meta["state_key"],
            "label_path":        meta["label_path"],
            "metadata_path":     meta["metadata_path"],
            "reference_graph_path":  meta["reference_graph_path"],
            "reference_edges_path":  meta["reference_edges_path"],
            "confidence_mode":   "all",
            "exclude_uncertain_states": True,
            "excluded_labels":   excluded,
            "n_states":          meta["n_states"],
            "n_graph_edges":     meta["n_graph_edges"],
            "n_cell_labels":     meta.get("n_cell_labels", 0),
            "source_h5ad":       meta.get("input_h5ad") or "",
            "label_mode":        lm,
            "label_type":        meta.get("label_type"),
            "analysis_role":     meta["analysis_role"],
            "coarse_label_key":  meta.get("coarse_state_key")
                                 or "final_milestone_label_coarse",
            "expanded_label_key": meta.get("expanded_state_key"),
            "notes": (
                f"Official frozen silver-standard milestone provider"
                f" for {meta['dataset_id']}."
                f" state_key={meta['state_key']}."
                " Frozen expanded transition labels from Stage 2 trajectory-aware"
                " annotation are used for official marker-FM metrics."
                " Coarse labels are retained as compatibility metadata."
                " Excluded from official metrics:"
                " ambiguous, unknown_or_ood."
            ),
        }

        registry["providers"][pid] = entry
        if not dry_run:
            print(f"{tag}   registry entry updated: {pid}")

    out_path = output_registry_path or registry_path
    if dry_run:
        print(f"{tag} [dry-run] would write registry -> {out_path}")
        return

    with out_path.open("w", encoding="utf-8") as fh:
        yaml.dump(
            registry, fh,
            default_flow_style=False,
            sort_keys=False,
            allow_unicode=True,
        )
    print(f"{tag} registry written -> {out_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Build milestone ground-truth providers and update registry.yaml "
            "(scTimeBench annotation Step 4)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Label mode reference:\n"
            "  official_silver   frozen silver-standard Stage 2 provider (DEFAULT)\n"
        ),
    )
    parser.add_argument(
        "--markers-yaml",
        default=_DEFAULT_MARKERS_YAML,
        help="Path to milestone_markers.yaml (default: %(default)s).",
    )
    parser.add_argument(
        "--registry-yaml",
        default=_DEFAULT_REGISTRY_YAML,
        help="Path to ground_truth/registry.yaml (default: %(default)s).",
    )
    parser.add_argument(
        "--provider-root",
        default=_DEFAULT_PROVIDER_ROOT,
        dest="provider_root",
        help="Root directory for provider subdirectories (default: %(default)s).",
    )
    parser.add_argument(
        "--dataset-id",
        choices=list(_VALID_DATASET_IDS),
        default=None,
        help="Build providers for a single dataset only (default: all datasets).",
    )
    parser.add_argument(
        "--label-mode",
        choices=list(_VALID_LABEL_MODE_ARGS),
        default=None,
        dest="label_mode",
        help=(
            "Which label mode(s) to build.  "
            "Omit for default (official_silver only)."
        ),
    )
    parser.add_argument(
        "--input-h5ad",
        default=None,
        dest="input_h5ad",
        help=(
            "Stage 2 h5ad with obs[state_key] to populate state_labels.tsv.  "
            "For official_silver, requires final_milestone_label_expanded.  "
            "If omitted, writes schema-only TSVs."
        ),
    )
    parser.add_argument(
        "--output-registry",
        default=None,
        dest="output_registry",
        help=(
            "Write updated registry to this path instead of overwriting "
            "--registry-yaml (useful for testing)."
        ),
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "Load config, validate providers, print plan, write nothing.  "
            "If --input-h5ad is given, also validates required state_key exists."
        ),
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    tag = "[build_milestone_providers]"
    if args.dry_run:
        print(f"{tag} dry-run mode: no files will be written.")

    from benchmark.annotation.marker_utils import load_markers_yaml

    markers_yaml_path = Path(args.markers_yaml)
    marker_config = load_markers_yaml(markers_yaml_path)

    datasets = [args.dataset_id] if args.dataset_id else list(_VALID_DATASET_IDS)
    modes    = _resolve_modes(args.label_mode)
    provider_root = Path(args.provider_root)

    print(f"{tag} datasets={datasets}")
    print(f"{tag} modes={modes}")

    all_metadatas = []
    for dataset_id in datasets:
        if dataset_id not in marker_config:
            print(
                f"{tag} WARNING: {dataset_id} not in markers YAML; skipping.",
                file=sys.stderr,
            )
            continue
        ds_cfg = marker_config[dataset_id]
        for mode in modes:
            print(f"{tag} {dataset_id} / {mode} ...")
            meta = _build_one_provider(
                dataset_id=dataset_id,
                label_mode=mode,
                ds_cfg=ds_cfg,
                provider_root=provider_root,
                markers_yaml_path=markers_yaml_path,
                input_h5ad=args.input_h5ad,
                dry_run=args.dry_run,
                tag=tag,
            )
            all_metadatas.append(meta)

    _update_registry(
        registry_path=Path(args.registry_yaml),
        output_registry_path=(
            Path(args.output_registry) if args.output_registry else None
        ),
        provider_metadatas=all_metadatas,
        dry_run=args.dry_run,
        tag=tag,
    )

    print(f"{tag} done. {len(all_metadatas)} provider(s) processed.")
    sys.exit(0)


if __name__ == "__main__":
    main()
