#!/usr/bin/env python3
"""
Build the frozen ground-truth provider for GSE175634 cardiac differentiation
(author-annotation silver provider).

# ─────────────────────────────────────────────────────────────────────────────
# DATASET ID NOTE
# ─────────────────────────────────────────────────────────────────────────────
# The canonical GEO accession is GSE175634.
# The user originally requested "GSE174534" — that is a typo (digits transposed).
# All source files in data/gse175634/ are named GSE175634_*, confirming
# GSE175634 is the correct identifier. This script uses GSE175634 throughout.
# ─────────────────────────────────────────────────────────────────────────────

Provider ID : gse175634_cardiac_silver_v1
Label mode  : official_silver  (author-provided 'type' column, no silver-standard
              annotation generation needed)
Input h5ad  : benchmark/inputs/gse175634_cardiac_author_hvg2000/
                GSE175634_cardiac_author_HVG2000_benchmark_input.h5ad
Output dir  : benchmark/ground_truth/providers/gse175634_cardiac_silver_v1/

Silver-standard annotation step: SKIPPED.
  GSE175634 carries complete per-cell state labels from the original authors
  (the 'type' column). No marker-seed or trajectory-aware annotation pass
  is performed.

Reference graph (from file catalog + publication biology):
    IPSC -> MES -> CMES -> PROG -> CM
                              \\-> CF
  6 named states, 5 directed edges. UNK cells are retained in the h5ad
  but excluded from official metrics via excluded_labels.

Run from the repository root:
    python scripts/build_gse175634_cardiac_author_provider.py [--overwrite]

Prerequisite: the benchmark input h5ad must exist first:
    python scripts/build_gse175634_cardiac_author_input.py
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd

# ─── Constants ────────────────────────────────────────────────────────────────

PROVIDER_ID   = "gse175634_cardiac_silver_v1"
DATASET_ID    = "GSE175634"
LABEL_MODE    = "official_silver"
LABEL_TYPE    = "frozen_silver_standard"
STATE_KEY     = "final_milestone_label_coarse"
GRAPH_TYPE    = "author_annotation_cardiac_differentiation_milestone_graph"
ANALYSIS_ROLE = "primary_report"
GENERATED_BY  = "scripts/build_gse175634_cardiac_author_provider.py"
SOURCE_LABEL  = "author_type_column"

SOURCE_H5AD = (
    "benchmark/inputs/gse175634_cardiac_author_hvg2000/"
    "GSE175634_cardiac_author_HVG2000_benchmark_input.h5ad"
)

# Ordered confirmed states (UNK is excluded from official metrics)
MILESTONE_ORDER: list[str] = [
    "IPSC",   # 0  – pluripotent origin
    "MES",    # 1  – mesoderm intermediate
    "CMES",   # 2  – cardiac mesoderm intermediate
    "PROG",   # 3  – progenitor branch point
    "CM",     # 4  – cardiomyocyte endpoint
    "CF",     # 5  – cardiac fibroblast endpoint
]

# State descriptions (for state_metadata.tsv)
MILESTONE_NOTES: dict[str, str] = {
    "IPSC":  "Induced pluripotent stem cell starting state.",
    "MES":   "Mesoderm-like intermediate; follows IPSC at early timepoints.",
    "CMES":  "Cardiac mesoderm intermediate; follows MES.",
    "PROG":  "Progenitor state at bifurcation point before CM/CF resolution.",
    "CM":    "Cardiomyocyte terminal branch endpoint.",
    "CF":    "Cardiac fibroblast-like terminal branch endpoint.",
}

# Directed edges (source, target)
EDGES: list[tuple[str, str]] = [
    ("IPSC",  "MES"),
    ("MES",   "CMES"),
    ("CMES",  "PROG"),
    ("PROG",  "CM"),
    ("PROG",  "CF"),
]

EDGE_WEIGHT     = 1.0
EDGE_CONFIDENCE = "high"
EDGE_TYPE       = "author_annotation_cardiac_differentiation_milestone_order"

EXCLUDED_LABELS = ["UNK"]


# ─── Project-root detection ───────────────────────────────────────────────────

def find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
        raise FileNotFoundError(f"TRAJ_PROJECT_ROOT={env!r} does not exist.")
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "benchmark").is_dir() and (candidate / "scripts").is_dir():
            return candidate
    return here.parent


# ─── Helpers ─────────────────────────────────────────────────────────────────

def _load_h5ad(h5ad_path: Path):
    try:
        import anndata as ad
    except ImportError:
        raise ImportError("anndata is required.  pip install anndata")
    print(f"  Loading h5ad: {h5ad_path}")
    return ad.read_h5ad(h5ad_path, backed="r")


def _build_reference_graph(n_cells_per_state: dict[str, int]) -> dict:
    """Return the reference_graph.json dict."""
    nodes = []
    for order, state_id in enumerate(MILESTONE_ORDER):
        nodes.append({
            "id":     state_id,
            "label":  state_id,
            "status": "confirmed",
            "role":   "primary_milestone",
            "order":  order,
        })

    edges = []
    for src, tgt in EDGES:
        edges.append({
            "source":        src,
            "target":        tgt,
            "weight":        EDGE_WEIGHT,
            "confidence":    EDGE_CONFIDENCE,
            "source_status": "confirmed",
            "target_status": "confirmed",
            "edge_type":     EDGE_TYPE,
        })

    return {
        "_meta": {
            "provider_id":       PROVIDER_ID,
            "dataset_id":        DATASET_ID,
            "label_mode":        LABEL_MODE,
            "label_type":        LABEL_TYPE,
            "state_key":         STATE_KEY,
            "graph_type":        GRAPH_TYPE,
            "version":           "v1",
            "generated_by":      GENERATED_BY,
            "analysis_role":     ANALYSIS_ROLE,
            "deprecated":        False,
            "n_states":          len(MILESTONE_ORDER),
            "n_edges":           len(EDGES),
            "silver_standard_skipped": True,
            "silver_standard_note": (
                "GSE175634 uses author-provided cell-state annotations (type column). "
                "No marker-seed or trajectory-aware annotation generation was performed."
            ),
            "policy": [
                "Official frozen silver provider for GSE175634 cardiac differentiation.",
                "Labels taken directly from the author 'type' column; no annotation generation step.",
                "Nodes: IPSC, MES, CMES, PROG (bifurcation), CM, CF (two terminal branches).",
                "Graph topology: linear IPSC→MES→CMES→PROG with bifurcation to CM and CF.",
                "UNK cells are retained in the h5ad but excluded from official metrics.",
                f"Source: {SOURCE_H5AD}",
            ],
            "n_cells_per_state": n_cells_per_state,
        },
        "nodes": nodes,
        "edges": edges,
    }


def _build_state_metadata(n_cells_per_state: dict[str, int]) -> pd.DataFrame:
    rows = []
    for order, state_id in enumerate(MILESTONE_ORDER):
        rows.append({
            "state_id": state_id,
            "label":    state_id,
            "status":   "confirmed",
            "role":     "primary_milestone",
            "order":    order,
            "n_cells":  n_cells_per_state.get(state_id, 0),
            "notes":    MILESTONE_NOTES[state_id],
        })
    return pd.DataFrame(rows)


# ─── Main ─────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build GSE175634 cardiac author ground-truth provider."
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Allow overwriting existing provider outputs."
    )
    parser.add_argument(
        "--project-root", type=Path, default=None,
        help="Override project root directory (default: auto-detected)."
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    project_root = args.project_root or find_project_root()
    h5ad_path    = project_root / SOURCE_H5AD
    out_dir      = (
        project_root / "benchmark" / "ground_truth" / "providers" / PROVIDER_ID
    )

    print("=" * 64)
    print(f"Provider build: {PROVIDER_ID}")
    print(f"  Dataset      : {DATASET_ID}")
    print(f"  Label mode   : {LABEL_MODE}  (author annotations; silver step SKIPPED)")
    print(f"  Source h5ad  : {h5ad_path}")
    print(f"  Output dir   : {out_dir}")
    print("=" * 64)

    if not h5ad_path.exists():
        print(
            f"\nERROR: benchmark input h5ad not found:\n  {h5ad_path}\n"
            "Run scripts/build_gse175634_cardiac_author_input.py first.",
            file=sys.stderr,
        )
        sys.exit(1)

    if out_dir.exists() and not args.overwrite:
        existing = list(out_dir.iterdir())
        if existing:
            raise FileExistsError(
                f"{out_dir} already exists and is non-empty. "
                "Pass --overwrite to rebuild."
            )
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── load h5ad ─────────────────────────────────────────────────────────────
    print("\n[step 1] Loading benchmark input h5ad …")
    adata = _load_h5ad(h5ad_path)
    print(f"  n_obs={adata.n_obs}, n_vars={adata.n_vars}")

    if STATE_KEY not in adata.obs.columns:
        raise KeyError(
            f"Required obs column '{STATE_KEY}' not found in {h5ad_path}. "
            "Re-run build_gse175634_cardiac_author_input.py."
        )
    if "cell_id" not in adata.obs.columns:
        raise KeyError("Required obs column 'cell_id' missing.")

    # ── state label counts ────────────────────────────────────────────────────
    label_counts = adata.obs[STATE_KEY].value_counts().to_dict()
    n_cells_per_state = {s: int(label_counts.get(s, 0)) for s in MILESTONE_ORDER}
    n_unk  = int(label_counts.get("UNK", 0))
    n_labeled = sum(n_cells_per_state.values())
    print(f"  State label counts (confirmed states): {n_cells_per_state}")
    print(f"  UNK cells (excluded from metrics)    : {n_unk}")
    print(f"  Labeled cells (official metrics)      : {n_labeled}")

    # ── state_labels.tsv ──────────────────────────────────────────────────────
    print("\n[step 2] Writing state_labels.tsv …")
    obs = adata.obs.copy()
    # Only emit rows for confirmed states (not UNK); UNK cells are not
    # represented in the provider's label file.
    mask_confirmed = obs[STATE_KEY].isin(set(MILESTONE_ORDER))
    labels_df = obs.loc[mask_confirmed, []].copy()
    labels_df["cell_id"]    = obs.loc[mask_confirmed, "cell_id"].values
    labels_df["state_id"]   = obs.loc[mask_confirmed, STATE_KEY].values
    labels_df["confidence"] = 1.0
    labels_df["label_mode"] = LABEL_MODE
    labels_df["provider_id"]= PROVIDER_ID
    labels_df = labels_df[
        ["cell_id", "state_id", "confidence", "label_mode", "provider_id"]
    ].reset_index(drop=True)
    labels_path = out_dir / "state_labels.tsv"
    labels_df.to_csv(labels_path, sep="\t", index=False)
    print(f"  Written: {labels_path}  ({len(labels_df):,} rows)")

    # ── state_metadata.tsv ────────────────────────────────────────────────────
    print("\n[step 3] Writing state_metadata.tsv …")
    meta_df = _build_state_metadata(n_cells_per_state)
    meta_path = out_dir / "state_metadata.tsv"
    meta_df.to_csv(meta_path, sep="\t", index=False)
    print(f"  Written: {meta_path}")

    # ── reference_graph.json ──────────────────────────────────────────────────
    print("\n[step 4] Writing reference_graph.json …")
    graph = _build_reference_graph(n_cells_per_state)
    graph_path = out_dir / "reference_graph.json"
    with graph_path.open("w", encoding="utf-8") as f:
        json.dump(graph, f, indent=2, ensure_ascii=True)
    print(f"  Written: {graph_path}")

    # ── reference_graph_edges.csv ─────────────────────────────────────────────
    print("\n[step 5] Writing reference_graph_edges.csv …")
    edges_df = pd.DataFrame(
        [
            {
                "source":     src,
                "target":     tgt,
                "weight":     EDGE_WEIGHT,
                "confidence": EDGE_CONFIDENCE,
                "edge_type":  EDGE_TYPE,
            }
            for src, tgt in EDGES
        ]
    )
    edges_path = out_dir / "reference_graph_edges.csv"
    edges_df.to_csv(edges_path, index=False)
    print(f"  Written: {edges_path}")

    # ── ground_truth_metadata.json ────────────────────────────────────────────
    print("\n[step 6] Writing ground_truth_metadata.json …")
    gt_meta = {
        "provider_id":             PROVIDER_ID,
        "dataset_id":              DATASET_ID,
        "label_mode":              LABEL_MODE,
        "label_type":              LABEL_TYPE,
        "state_key":               STATE_KEY,
        "annotation_method":       "author_annotation_direct",
        "silver_standard_skipped": True,
        "silver_standard_note": (
            "Labels taken directly from the author 'type' column. "
            "No marker-seed or trajectory-aware annotation generation performed."
        ),
        "analysis_role":    ANALYSIS_ROLE,
        "n_confirmed_states": len(MILESTONE_ORDER),
        "n_graph_edges":    len(EDGES),
        "excluded_labels":  EXCLUDED_LABELS,
        "n_cells_labeled":  n_labeled,
        "n_cells_unk":      n_unk,
        "n_cells_total":    int(adata.n_obs),
        "n_cells_per_state": n_cells_per_state,
        "source_h5ad":      SOURCE_H5AD,
        "generated_by":     GENERATED_BY,
        "output_dir":       str(out_dir),
        "outputs": [
            "state_labels.tsv",
            "state_metadata.tsv",
            "reference_graph.json",
            "reference_graph_edges.csv",
            "ground_truth_metadata.json",
        ],
    }
    gt_meta_path = out_dir / "ground_truth_metadata.json"
    with gt_meta_path.open("w", encoding="utf-8") as f:
        json.dump(gt_meta, f, indent=2, ensure_ascii=True)
    print(f"  Written: {gt_meta_path}")

    # ── summary ───────────────────────────────────────────────────────────────
    print()
    print("=" * 64)
    print("Provider build complete.")
    print(f"  Provider ID    : {PROVIDER_ID}")
    print(f"  n_states       : {len(MILESTONE_ORDER)}")
    print(f"  n_edges        : {len(EDGES)}")
    print(f"  n_labeled cells: {n_labeled:,}  (UNK excluded: {n_unk:,})")
    print(f"  Output dir     : {out_dir}")
    print()
    print("Next step:")
    print("  1. Add provider entry to benchmark/ground_truth/registry.yaml")
    print("     (already done if you used setup scripts).")
    print("  2. Run benchmark methods using configs in")
    print("     benchmark/configs/runtime/*gse175634_cardiac_silver*_formal.yaml")
    print("=" * 64)


if __name__ == "__main__":
    main()
