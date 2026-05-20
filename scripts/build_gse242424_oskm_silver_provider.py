#!/usr/bin/env python3
"""
Build the frozen ground-truth provider for GSE242424 OSKM reprogramming
(author-cluster-matched silver subset).

Provider ID : gse242424_oskm_reprogramming_silver_v1
Input h5ad  : benchmark/inputs/gse242424_author_cluster_matched/
                GSE242424_author_cluster_matched_HVG2000_benchmark_input.h5ad
Output dir  : benchmark/ground_truth/providers/
                gse242424_oskm_reprogramming_silver_v1/

Run from the repository root:
    python scripts/build_gse242424_oskm_silver_provider.py

All outputs are written deterministically; running the script twice
produces byte-identical files (cells are sorted by cell_id before writing).
The original full-dataset h5ad is never read or modified.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PROVIDER_ID = "gse242424_oskm_reprogramming_silver_v1"
DATASET_ID = "GSE242424"
LABEL_MODE = "official_silver"
LABEL_TYPE = "frozen_silver_standard"
STATE_KEY = "final_milestone_label_coarse"
GRAPH_TYPE = "author_cluster_defined_oskm_reprogramming_milestone_graph"
ANALYSIS_ROLE = "primary_report"
GENERATED_BY = "scripts/build_gse242424_oskm_silver_provider.py"
FINAL_MILESTONE_SOURCE = "author_atac_to_rna_cluster_transfer"

SOURCE_H5AD = (
    "benchmark/inputs/gse242424_author_cluster_matched/"
    "GSE242424_author_cluster_matched_HVG2000_benchmark_input.h5ad"
)

# Ordered coarse milestone nodes (index == order)
MILESTONE_ORDER: list[str] = [
    "fibroblast",               # 0
    "fibroblast_like_stalled",  # 1
    "keratinocyte_like",        # 2
    "hOSK",                     # 3
    "partial_intermediate",     # 4
    "partially_reprogrammed",   # 5
    "xOSK",                     # 6
    "primary_intermediate",     # 7
    "pre_iPSC",                 # 8
    "iPSC",                     # 9
]

# Directed edges (source, target)
EDGES: list[tuple[str, str]] = [
    ("fibroblast",           "fibroblast_like_stalled"),
    ("fibroblast",           "keratinocyte_like"),
    ("fibroblast",           "hOSK"),
    ("hOSK",                 "partial_intermediate"),
    ("partial_intermediate", "partially_reprogrammed"),
    ("fibroblast",           "xOSK"),
    ("xOSK",                 "primary_intermediate"),
    ("primary_intermediate", "pre_iPSC"),
    ("pre_iPSC",             "iPSC"),
]

EDGE_WEIGHT = 1.0
EDGE_CONFIDENCE = "high"
EDGE_TYPE = "author_cluster_oskm_reprogramming_milestone_order"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _project_root() -> Path:
    """Walk upward from this script until we find the benchmark/ directory."""
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "benchmark").is_dir():
            return candidate
    raise RuntimeError("Cannot locate project root (no benchmark/ directory found).")


def _read_h5ad(path: Path):
    try:
        import anndata as ad
    except ImportError:
        sys.exit("anndata is required. Install with: pip install anndata")
    return ad.read_h5ad(path)


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------

def build_reference_graph_json() -> dict:
    nodes = [
        {
            "id": label,
            "label": label,
            "status": "confirmed",
            "role": "primary_milestone",
            "order": idx,
        }
        for idx, label in enumerate(MILESTONE_ORDER)
    ]
    edges = [
        {
            "source": src,
            "target": tgt,
            "weight": EDGE_WEIGHT,
            "confidence": EDGE_CONFIDENCE,
            "source_status": "confirmed",
            "target_status": "confirmed",
            "edge_type": EDGE_TYPE,
        }
        for src, tgt in EDGES
    ]
    meta = {
        "provider_id": PROVIDER_ID,
        "dataset_id": DATASET_ID,
        "label_mode": LABEL_MODE,
        "label_type": LABEL_TYPE,
        "state_key": STATE_KEY,
        "graph_type": GRAPH_TYPE,
        "version": "v1",
        "generated_by": GENERATED_BY,
        "analysis_role": ANALYSIS_ROLE,
        "deprecated": False,
        "n_states": len(MILESTONE_ORDER),
        "n_edges": len(EDGES),
        "policy": [
            "Official frozen silver-standard coarse milestone reference graph for GSE242424 OSKM reprogramming.",
            "Nodes are the 10 author-defined coarse ATAC/RNA cluster milestones.",
            "Graph topology reflects the OSKM reprogramming branching order from the original publication.",
            "Two divergent arms: fibroblast->hOSK->partial_intermediate->partially_reprogrammed "
            "and fibroblast->xOSK->primary_intermediate->pre_iPSC->iPSC.",
            "keratinocyte_like and fibroblast_like_stalled are stalled/off-target branches from fibroblast.",
            "All edges have confidence=high; all nodes have status=confirmed.",
            f"Source: {SOURCE_H5AD}",
        ],
    }
    return {"_meta": meta, "nodes": nodes, "edges": edges}


def write_reference_graph_json(out_dir: Path) -> None:
    graph = build_reference_graph_json()
    out_path = out_dir / "reference_graph.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(graph, f, indent=2)
        f.write("\n")
    print(f"  Wrote {out_path}  ({len(graph['nodes'])} nodes, {len(graph['edges'])} edges)")


def write_reference_graph_edges_csv(out_dir: Path) -> None:
    lines = ["source,target,weight,confidence,edge_type"]
    for src, tgt in EDGES:
        lines.append(f"{src},{tgt},{EDGE_WEIGHT},{EDGE_CONFIDENCE},{EDGE_TYPE}")
    out_path = out_dir / "reference_graph_edges.csv"
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"  Wrote {out_path}  ({len(EDGES)} edges)")


def write_state_metadata_tsv(out_dir: Path) -> None:
    """One row per milestone state. n_cells populated from h5ad below."""
    # Roles and notes per milestone
    notes_map: dict[str, str] = {
        "fibroblast":              "Starting fibroblast state; principal trajectory origin.",
        "fibroblast_like_stalled": "Stalled/off-target fibroblast-like branch; does not progress to iPSC.",
        "keratinocyte_like":       "Off-target keratinocyte-like branch diverging from fibroblast.",
        "hOSK":                    "Early hOSK-transduced state; entry point of the main reprogramming arm.",
        "partial_intermediate":    "Partial intermediate on the hOSK->partially_reprogrammed arm.",
        "partially_reprogrammed":  "Partially reprogrammed state; stalled on the hOSK arm.",
        "xOSK":                    "Early xOSK-transduced state; entry point of the productive reprogramming arm.",
        "primary_intermediate":    "Primary intermediate on the xOSK->pre_iPSC->iPSC arm.",
        "pre_iPSC":                "Pre-iPSC state; immediate precursor to the iPSC endpoint.",
        "iPSC":                    "iPSC endpoint; fully reprogrammed pluripotent state.",
    }
    header = "state_id\tlabel\tstatus\trole\torder\tnotes"
    rows = [header]
    for idx, sid in enumerate(MILESTONE_ORDER):
        rows.append(
            f"{sid}\t{sid}\tconfirmed\tprimary_milestone\t{idx}\t{notes_map[sid]}"
        )
    out_path = out_dir / "state_metadata.tsv"
    out_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"  Wrote {out_path}  ({len(MILESTONE_ORDER)} states, n_cells added by main)")


def write_state_metadata_tsv_with_counts(out_dir: Path, counts: dict[str, int]) -> None:
    """Rewrite state_metadata.tsv with n_cells populated from h5ad counts."""
    notes_map: dict[str, str] = {
        "fibroblast":              "Starting fibroblast state; principal trajectory origin.",
        "fibroblast_like_stalled": "Stalled/off-target fibroblast-like branch; does not progress to iPSC.",
        "keratinocyte_like":       "Off-target keratinocyte-like branch diverging from fibroblast.",
        "hOSK":                    "Early hOSK-transduced state; entry point of the main reprogramming arm.",
        "partial_intermediate":    "Partial intermediate on the hOSK->partially_reprogrammed arm.",
        "partially_reprogrammed":  "Partially reprogrammed state; stalled on the hOSK arm.",
        "xOSK":                    "Early xOSK-transduced state; entry point of the productive reprogramming arm.",
        "primary_intermediate":    "Primary intermediate on the xOSK->pre_iPSC->iPSC arm.",
        "pre_iPSC":                "Pre-iPSC state; immediate precursor to the iPSC endpoint.",
        "iPSC":                    "iPSC endpoint; fully reprogrammed pluripotent state.",
    }
    header = "state_id\tlabel\tstatus\trole\torder\tn_cells\tnotes"
    rows = [header]
    for idx, sid in enumerate(MILESTONE_ORDER):
        n = counts.get(sid, 0)
        rows.append(
            f"{sid}\t{sid}\tconfirmed\tprimary_milestone\t{idx}\t{n}\t{notes_map[sid]}"
        )
    out_path = out_dir / "state_metadata.tsv"
    out_path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"  Wrote {out_path}  ({len(MILESTONE_ORDER)} states with n_cells)")


def write_state_labels_tsv(out_dir: Path, obs) -> None:
    """One row per cell: cell_id, state_id, confidence, label_mode, provider_id."""
    import pandas as pd
    df = obs[["cell_id", STATE_KEY]].copy()
    df = df.sort_values("cell_id").reset_index(drop=True)
    df.rename(columns={STATE_KEY: "state_id"}, inplace=True)
    df["confidence"] = 1.0
    df["label_mode"] = LABEL_MODE
    df["provider_id"] = PROVIDER_ID
    out_path = out_dir / "state_labels.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  Wrote {out_path}  ({len(df)} rows)")


def write_annotation_votes_tsv(out_dir: Path, obs) -> None:
    """Per-cell provenance table."""
    import pandas as pd
    cols = [
        "cell_id",
        "author_old_atac_cluster",
        "author_cluster_id",
        "author_cluster_label",
        STATE_KEY,
    ]
    df = obs[cols].copy()
    df = df.sort_values("cell_id").reset_index(drop=True)
    df.rename(columns={STATE_KEY: "final_milestone_label_coarse"}, inplace=True)
    df["final_milestone_confidence"] = 1.0
    df["final_milestone_source"] = FINAL_MILESTONE_SOURCE
    df["provider_id"] = PROVIDER_ID
    df["label_mode"] = LABEL_MODE
    out_path = out_dir / "annotation_votes.tsv"
    df.to_csv(out_path, sep="\t", index=False)
    print(f"  Wrote {out_path}  ({len(df)} rows)")


def write_ground_truth_metadata_json(out_dir: Path, counts: dict[str, int]) -> None:
    """Write provider-level provenance and policy metadata."""
    meta = {
        "provider_id": PROVIDER_ID,
        "dataset_id": DATASET_ID,
        "label_mode": LABEL_MODE,
        "label_type": LABEL_TYPE,
        "annotation_method": "author_cluster_oskm_reprogramming_silver_provider",
        "status": "frozen_silver_standard_milestone_provider",
        "state_key": STATE_KEY,
        "confidence_key": "final_milestone_confidence",
        "source_key": "final_milestone_source",
        "analysis_role": ANALYSIS_ROLE,
        "version": "v1",
        "graph_type": GRAPH_TYPE,
        "n_states": len(MILESTONE_ORDER),
        "n_graph_edges": len(EDGES),
        "n_cell_labels": int(sum(counts.values())),
        "excluded_from_official_metrics": [],
        "primary_milestones": MILESTONE_ORDER,
        "coarse_label_counts": {k: int(counts.get(k, 0)) for k in MILESTONE_ORDER},
        "input_h5ad": SOURCE_H5AD,
        "match_summary_path": (
            "benchmark/inputs/gse242424_author_cluster_matched/"
            "gse242424_author_cluster_match_summary.json"
        ),
        "state_labels_status": "populated_from_author_cluster_matched_h5ad",
        "reference_graph_path": (
            "benchmark/ground_truth/providers/"
            f"{PROVIDER_ID}/reference_graph.json"
        ),
        "reference_edges_path": (
            "benchmark/ground_truth/providers/"
            f"{PROVIDER_ID}/reference_graph_edges.csv"
        ),
        "label_path": (
            "benchmark/ground_truth/providers/"
            f"{PROVIDER_ID}/state_labels.tsv"
        ),
        "metadata_path": (
            "benchmark/ground_truth/providers/"
            f"{PROVIDER_ID}/state_metadata.tsv"
        ),
        "annotation_votes_path": (
            "benchmark/ground_truth/providers/"
            f"{PROVIDER_ID}/annotation_votes.tsv"
        ),
        "generated_by": GENERATED_BY,
        "provenance": [
            "GSE242424/GSE242423 OSKM human fibroblast reprogramming time course",
            "Author ATAC-to-RNA cluster transfer table from kundajelab/scATAC-reprog",
            "Author cluster conversion table from kundajelab/scATAC-reprog",
            "Zenodo 10.5281/zenodo.8313962 clusters.tsv",
        ],
        "notes": [
            "The provider uses the 59,187-cell author-cluster-matched subset, not the full 156,969-cell local GSE242424 input.",
            "The reference graph includes productive, partial, stalled, and off-target OSKM reprogramming branches.",
            "No cells are excluded from official metrics because this provider has no ambiguous or unknown labels.",
        ],
    }
    out_path = out_dir / "ground_truth_metadata.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
        f.write("\n")
    print(f"  Wrote {out_path}")


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def validate(out_dir: Path, obs) -> bool:
    import pandas as pd
    ok = True

    # h5ad shape
    n_cells = len(obs)
    if n_cells != 59187:
        print(f"  [FAIL] h5ad has {n_cells} cells, expected 59187")
        ok = False
    else:
        print(f"  [OK] h5ad shape: {n_cells} cells")

    # null labels
    nulls = obs[STATE_KEY].isnull().sum()
    if nulls > 0:
        print(f"  [FAIL] {nulls} null values in {STATE_KEY}")
        ok = False
    else:
        print(f"  [OK] zero null final_milestone_label_coarse values")

    # state_labels.tsv row count
    sl = pd.read_csv(out_dir / "state_labels.tsv", sep="\t")
    if len(sl) != 59187:
        print(f"  [FAIL] state_labels.tsv has {len(sl)} rows, expected 59187")
        ok = False
    else:
        print(f"  [OK] state_labels.tsv: {len(sl)} rows")

    # annotation_votes.tsv row count
    av = pd.read_csv(out_dir / "annotation_votes.tsv", sep="\t")
    if len(av) != 59187:
        print(f"  [FAIL] annotation_votes.tsv has {len(av)} rows, expected 59187")
        ok = False
    else:
        print(f"  [OK] annotation_votes.tsv: {len(av)} rows")

    # reference_graph.json nodes/edges
    with open(out_dir / "reference_graph.json", encoding="utf-8") as f:
        g = json.load(f)
    n_nodes = len(g["nodes"])
    n_edges = len(g["edges"])
    if n_nodes != 10:
        print(f"  [FAIL] reference_graph.json has {n_nodes} nodes, expected 10")
        ok = False
    else:
        print(f"  [OK] reference_graph.json: {n_nodes} nodes")
    if n_edges != 9:
        print(f"  [FAIL] reference_graph.json has {n_edges} edges, expected 9")
        ok = False
    else:
        print(f"  [OK] reference_graph.json: {n_edges} edges")

    # reference_graph_edges.csv row count
    import csv
    with open(out_dir / "reference_graph_edges.csv", encoding="utf-8") as f:
        data_rows = sum(1 for _ in csv.reader(f)) - 1  # subtract header
    if data_rows != 9:
        print(f"  [FAIL] reference_graph_edges.csv has {data_rows} data rows, expected 9")
        ok = False
    else:
        print(f"  [OK] reference_graph_edges.csv: {data_rows} data rows")

    # All graph node IDs appear in h5ad labels
    h5ad_labels = set(obs[STATE_KEY].unique())
    graph_nodes = {n["id"] for n in g["nodes"]}
    missing = graph_nodes - h5ad_labels
    if missing:
        print(f"  [FAIL] graph node IDs not in h5ad labels: {missing}")
        ok = False
    else:
        print(f"  [OK] all 10 graph node IDs present in h5ad labels")

    # metadata
    metadata_path = out_dir / "ground_truth_metadata.json"
    if not metadata_path.exists():
        print("  [FAIL] ground_truth_metadata.json missing")
        ok = False
    else:
        with open(metadata_path, encoding="utf-8") as f:
            meta = json.load(f)
        if meta.get("provider_id") != PROVIDER_ID:
            print("  [FAIL] ground_truth_metadata.json provider_id mismatch")
            ok = False
        elif meta.get("n_cell_labels") != 59187:
            print("  [FAIL] ground_truth_metadata.json n_cell_labels mismatch")
            ok = False
        else:
            print("  [OK] ground_truth_metadata.json")

    return ok


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    root = _project_root()
    h5ad_path = root / SOURCE_H5AD
    out_dir = (
        root
        / "benchmark"
        / "ground_truth"
        / "providers"
        / PROVIDER_ID
    )

    print(f"Project root : {root}")
    print(f"Input h5ad   : {h5ad_path}")
    print(f"Output dir   : {out_dir}")
    print()

    # Safety check – do not read the original full 156 969-cell h5ad
    if not h5ad_path.exists():
        sys.exit(f"Input h5ad not found: {h5ad_path}")

    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading h5ad …")
    adata = _read_h5ad(h5ad_path)
    obs = adata.obs.reset_index(drop=True)  # index == obs_names; cell_id col exists
    print(f"  Loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")

    # Cell-count per state (for state_metadata)
    counts: dict[str, int] = obs[STATE_KEY].value_counts().to_dict()

    print()
    print("Writing provider files …")
    write_reference_graph_json(out_dir)
    write_reference_graph_edges_csv(out_dir)
    write_state_metadata_tsv_with_counts(out_dir, counts)
    write_state_labels_tsv(out_dir, obs)
    write_annotation_votes_tsv(out_dir, obs)
    write_ground_truth_metadata_json(out_dir, counts)

    print()
    print("Validating …")
    ok = validate(out_dir, obs)

    print()
    if ok:
        print("All validations PASSED.")
        print(f"\nProvider written to:\n  {out_dir}")
        print("\nNext step: register the provider in")
        print("  benchmark/ground_truth/registry.yaml")
    else:
        print("One or more validations FAILED — review output above.")
        sys.exit(1)


if __name__ == "__main__":
    main()
