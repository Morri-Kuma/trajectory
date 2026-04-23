"""
export_v1.py
Promote the provisional scGPT annotation to a working v1 benchmark-ready asset.

Purpose
-------
Read the already-annotated AnnData (adata_scgpt_annotated.h5ad) and the
provisional cluster map (cluster_map_provisional.json), then produce:

  1. Cell-level v1 pseudostate label file  (TSV)
  2. State-level v1 metadata / QC summary  (TSV)
  3. Working v1 reference graph            (JSON + CSV edge table)
  4. Export provenance record              (JSON)

This is a working / silver-standard v1 asset. All file names, JSON metadata
blocks, and column headers make this explicit. Nothing here claims final
biological certainty.

This script does NOT
---------------------
  - Rerun embedding or recompute X_scGPT
  - Modify adata_scgpt_full.h5ad or adata_scgpt_annotated.h5ad
  - Merge or split clusters
  - Wire outputs into WOT / CellRank2 configs (manual step, deferred)

Reference graph method
----------------------
Consecutive-timepoint kNN transition counting:
  For each of the 14 consecutive timepoint pairs (t_i → t_{i+1}):
    - Extract X_scGPT embeddings for source-tp cells and target-tp cells
    - Compute full pairwise similarity matrix via dot product
      (valid because all embeddings are exactly L2-normalised: sim = cos)
    - For each source cell, identify its K=10 nearest target-tp neighbours
    - Tally source-state → target-state transition counts
  Aggregate counts across all 14 pairs, then row-normalise per source state.
  Keep edges where transition probability >= MIN_EDGE_WEIGHT.
  Assign edge confidence based on the status of both endpoint states.

Environment
-----------
  conda activate scgpt_env
  cd C:\\Users\\37620\\trajectory
  python benchmark\\scgpt\\export_v1.py

Expected runtime: ~5–10 min  (dominated by 14 pairwise kNN similarity sweeps)

Outputs
-------
  benchmark/datasets/scgpt_pseudostate_v1.tsv
  benchmark/datasets/scgpt_pseudostate_v1_metadata.tsv
  benchmark/datasets/scgpt_reference_graph_v1.json
  benchmark/datasets/scgpt_reference_graph_v1_edges.csv
  benchmark/results/scgpt/full/v1_export_metadata.json
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import entropy as scipy_entropy


# ---------------------------------------------------------------------------
# 0. Project-root detection
# ---------------------------------------------------------------------------

def _find_project_root() -> Path:
    import os
    # 1. Explicit override: TRAJ_PROJECT_ROOT env var (Shirokane HPC / CI).
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
        raise FileNotFoundError(
            f"TRAJ_PROJECT_ROOT={env!r} does not exist.  "
            "Correct the environment variable and retry."
        )
    # 2. Walk upward from this script until a directory containing 'data/' is found.
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "data").exists():
            return candidate
    # 3. Last-resort fallback.
    return here.parent


PROJECT_ROOT = _find_project_root()
print(f"[export_v1] Project root: {PROJECT_ROOT}")


# ---------------------------------------------------------------------------
# 1. Paths and parameters
# ---------------------------------------------------------------------------

ANNOTATED_H5AD = (
    PROJECT_ROOT / "benchmark" / "results" / "scgpt" / "full"
    / "adata_scgpt_annotated.h5ad"
)
MAP_PATH    = PROJECT_ROOT / "benchmark" / "scgpt" / "cluster_map_provisional.json"
DATASET_DIR = PROJECT_ROOT / "benchmark" / "datasets"
RESULT_DIR  = PROJECT_ROOT / "benchmark" / "results" / "scgpt" / "full"

# Column names (from annotate_provisional.py)
LEIDEN_KEY  = "leiden_scgpt_res0.5"
PS_COL      = "scgpt_pseudostate_provisional"
FAMILY_COL  = "scgpt_state_family"
STATUS_COL  = "scgpt_state_status"
TIME_KEY    = "abs_day"

# Graph construction parameters
KNN_K           = 10     # nearest neighbours per source cell
MIN_EDGE_WEIGHT = 0.02   # prune edges with transition probability < 2 %
CHUNK_SIZE      = 4000   # max source cells per dot-product chunk (memory guard)

# Output file names — version-tagged
LABEL_TSV      = DATASET_DIR / "scgpt_pseudostate_v1.tsv"
META_TSV       = DATASET_DIR / "scgpt_pseudostate_v1_metadata.tsv"
GRAPH_JSON     = DATASET_DIR / "scgpt_reference_graph_v1.json"
GRAPH_EDGE_CSV = DATASET_DIR / "scgpt_reference_graph_v1_edges.csv"
EXPORT_META    = RESULT_DIR  / "v1_export_metadata.json"

VERSION        = "v1"


# ---------------------------------------------------------------------------
# 2. Pre-flight checks
# ---------------------------------------------------------------------------

print("\n[export_v1] === Pre-flight checks ===")

for path, label in [
    (ANNOTATED_H5AD, "annotated h5ad"),
    (MAP_PATH,       "cluster map"),
]:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found:\n  {path}")
    print(f"  {label:<20}: {path.name}  ✓")

DATASET_DIR.mkdir(parents=True, exist_ok=True)
RESULT_DIR.mkdir(parents=True, exist_ok=True)
print(f"  Output dirs ready  ✓")


# ---------------------------------------------------------------------------
# 3. Load inputs
# ---------------------------------------------------------------------------

print(f"\n[export_v1] === Loading inputs ===")
t0 = time.time()

adata = sc.read_h5ad(ANNOTATED_H5AD)
print(f"  h5ad loaded in {time.time()-t0:.1f}s  →  {adata.n_obs:,} cells")

for col in [PS_COL, FAMILY_COL, STATUS_COL, TIME_KEY]:
    if col not in adata.obs.columns:
        raise KeyError(f"Required obs column missing: {col!r}")
if "X_scGPT" not in adata.obsm:
    raise KeyError("X_scGPT not in adata.obsm")

emb = adata.obsm["X_scGPT"].astype(np.float32)   # (75194, 512), L2-normalised

with open(MAP_PATH, "r", encoding="utf-8") as f:
    cluster_map = json.load(f)

map_version = cluster_map["_meta"].get("mapping_version", "unknown")
print(f"  Cluster map version: {map_version}")

all_states = sorted(adata.obs[PS_COL].unique())
n_states   = len(all_states)
state_idx  = {s: i for i, s in enumerate(all_states)}
print(f"  States: {n_states}  →  {all_states}")


# ---------------------------------------------------------------------------
# 4. Helper: status-based edge confidence
# ---------------------------------------------------------------------------

def _edge_confidence(status_src: str, status_tgt: str) -> str:
    """Return 'high', 'medium', or 'low' based on endpoint statuses."""
    if status_src == "uncertain" or status_tgt == "uncertain":
        return "low"
    if status_src == "confirmed" and status_tgt == "confirmed":
        return "high"
    return "medium"   # at least one merge_review endpoint


# Retrieve status for each pseudostate from the cluster map
ps_to_status = {
    entry["pseudostate"]: entry["status"]
    for entry in cluster_map["clusters"].values()
}
ps_to_family = {
    entry["pseudostate"]: entry["family"]
    for entry in cluster_map["clusters"].values()
}


# ---------------------------------------------------------------------------
# 5. Part A — Cell-level v1 label export
# ---------------------------------------------------------------------------

print(f"\n[export_v1] === Part A: cell-level label export ===")

label_df = pd.DataFrame({
    "cell_id"              : adata.obs.index,
    "scgpt_pseudostate_v1" : adata.obs[PS_COL].values,
    "scgpt_state_family"   : adata.obs[FAMILY_COL].values,
    "scgpt_state_status"   : adata.obs[STATUS_COL].values,
})

label_df.to_csv(LABEL_TSV, sep="\t", index=False)
print(f"  Saved {LABEL_TSV.name}  ({label_df.shape[0]:,} rows)")


# ---------------------------------------------------------------------------
# 6. Part B — State-level v1 metadata export
# ---------------------------------------------------------------------------

print(f"\n[export_v1] === Part B: state-level metadata export ===")

n_tp = adata.obs[TIME_KEY].nunique()
meta_records = []

for ps in all_states:
    mask = adata.obs[PS_COL] == ps
    sub  = adata.obs[mask]
    n    = mask.sum()

    tp_vc     = sub[TIME_KEY].value_counts().sort_index()
    dom_tp    = float(tp_vc.idxmax())
    dom_frac  = float(tp_vc.max() / n)
    probs     = (tp_vc / n).values
    norm_ent  = float(scipy_entropy(probs[probs > 0]) / np.log(n_tp))

    # Retrieve notes from cluster map (look up by pseudostate value)
    notes = next(
        (e["notes"] for e in cluster_map["clusters"].values()
         if e["pseudostate"] == ps),
        ""
    )

    meta_records.append({
        "scgpt_pseudostate_v1"  : ps,
        "scgpt_state_family"    : ps_to_family.get(ps, ""),
        "scgpt_state_status"    : ps_to_status.get(ps, ""),
        "n_cells"               : int(n),
        "frac_of_total"         : round(n / adata.n_obs, 4),
        "dominant_timepoint"    : dom_tp,
        "dominant_frac"         : round(dom_frac, 3),
        "norm_entropy"          : round(norm_ent, 3),
        "notes"                 : notes,
    })

meta_df = pd.DataFrame(meta_records)
meta_df.to_csv(META_TSV, sep="\t", index=False)
print(f"  Saved {META_TSV.name}  ({len(meta_df)} states)")


# ---------------------------------------------------------------------------
# 7. Part C — Consecutive-timepoint kNN reference graph construction
# ---------------------------------------------------------------------------

print(f"\n[export_v1] === Part C: reference graph construction ===")
print(f"  Method  : consecutive-timepoint kNN transition counting")
print(f"  K       : {KNN_K} nearest neighbours")
print(f"  Min edge: {MIN_EDGE_WEIGHT} (probability threshold)")

tps = sorted(adata.obs[TIME_KEY].unique())
tp_pairs = list(zip(tps[:-1], tps[1:]))
print(f"  Timepoint pairs: {len(tp_pairs)}")

# Accumulate raw transition counts across all timepoint pairs
# counts[i, j] = total kNN transitions from state i to state j
raw_counts = np.zeros((n_states, n_states), dtype=np.float64)
pair_contrib = np.zeros((n_states, n_states), dtype=np.int32)  # n pairs each edge appeared in

t_graph_start = time.time()

for t_src, t_tgt in tp_pairs:
    src_mask = adata.obs[TIME_KEY] == t_src
    tgt_mask = adata.obs[TIME_KEY] == t_tgt

    src_emb    = emb[src_mask.values]              # (n_src, 512)
    tgt_emb    = emb[tgt_mask.values]              # (n_tgt, 512)
    src_states_arr = adata.obs.loc[src_mask, PS_COL].values
    tgt_states_arr = adata.obs.loc[tgt_mask, PS_COL].values

    n_src = src_emb.shape[0]
    n_tgt = tgt_emb.shape[0]
    k_actual = min(KNN_K, n_tgt)  # guard: can't take more neighbours than targets

    # Map state labels to indices
    src_state_idx = np.array([state_idx[s] for s in src_states_arr], dtype=np.int32)
    tgt_state_idx = np.array([state_idx[s] for s in tgt_states_arr], dtype=np.int32)

    # Compute pairwise similarity in chunks to bound peak memory usage.
    # Since embeddings are L2-normalised, dot product = cosine similarity.
    # knn_state_votes[i, k] = state index of the k-th nearest tgt neighbour of src cell i
    knn_votes = np.empty((n_src, k_actual), dtype=np.int32)

    for chunk_start in range(0, n_src, CHUNK_SIZE):
        chunk_end   = min(chunk_start + CHUNK_SIZE, n_src)
        sim_chunk   = src_emb[chunk_start:chunk_end] @ tgt_emb.T   # (chunk, n_tgt)
        # Top-K: argpartition gives unsorted top-K (fast), sufficient for counting
        top_k_local = np.argpartition(sim_chunk, -k_actual, axis=1)[:, -k_actual:]
        knn_votes[chunk_start:chunk_end] = tgt_state_idx[top_k_local]

    # Tally transitions: for each of the K neighbour votes, add 1 to the
    # corresponding (src_state, tgt_state) bin. All K votes per source cell
    # are counted equally — no distance weighting.
    for ki in range(k_actual):
        np.add.at(raw_counts, (src_state_idx, knn_votes[:, ki]), 1)

    # Track which (state-pair) cells contributed to this timepoint pair
    # (used later to count n_timepoint_pairs per edge)
    for ki in range(k_actual):
        contributing = np.zeros((n_states, n_states), dtype=bool)
        np.logical_or.at(contributing, (src_state_idx, knn_votes[:, ki]), True)
    pair_contrib += contributing.astype(np.int32)

    print(f"  {t_src:6.2f} → {t_tgt:6.2f}:  "
          f"n_src={n_src:6,}  n_tgt={n_tgt:6,}  k={k_actual}")

graph_runtime = time.time() - t_graph_start
print(f"  kNN sweeps complete in {graph_runtime:.1f}s")

# Row-normalise: probability of transitioning from state i to state j
row_sums = raw_counts.sum(axis=1, keepdims=True)
row_sums[row_sums == 0] = 1.0   # avoid divide-by-zero for empty source states
prob_matrix = raw_counts / row_sums


# ---------------------------------------------------------------------------
# 8. Build graph nodes and edges
# ---------------------------------------------------------------------------

print(f"\n[export_v1] === Building graph nodes and edges ===")

# Nodes: one per pseudostate
nodes = []
for ps in all_states:
    mask = adata.obs[PS_COL] == ps
    n    = int(mask.sum())
    nodes.append({
        "id"                  : ps,
        "n_cells"             : n,
        "family"              : ps_to_family.get(ps, "unknown"),
        "status"              : ps_to_status.get(ps, "unknown"),
        "version"             : VERSION,
    })

# Edges: threshold, assign confidence, collect diagnostics
edges      = []
edge_rows  = []   # for the flat CSV

for i, ps_src in enumerate(all_states):
    status_src = ps_to_status.get(ps_src, "unknown")
    for j, ps_tgt in enumerate(all_states):
        if i == j:
            continue   # skip self-loops
        prob = float(prob_matrix[i, j])
        if prob < MIN_EDGE_WEIGHT:
            continue   # below threshold — prune

        status_tgt  = ps_to_status.get(ps_tgt, "unknown")
        confidence  = _edge_confidence(status_src, status_tgt)
        raw_count   = int(raw_counts[i, j])
        n_tp_pairs  = int(pair_contrib[i, j])

        edge = {
            "source"                   : ps_src,
            "target"                   : ps_tgt,
            "weight"                   : round(prob, 4),
            "raw_knn_count"            : raw_count,
            "n_timepoint_pairs"        : n_tp_pairs,
            "confidence"               : confidence,
            "source_status"            : status_src,
            "target_status"            : status_tgt,
        }
        edges.append(edge)
        edge_rows.append(edge)

print(f"  Nodes: {len(nodes)}")
print(f"  Edges (after threshold={MIN_EDGE_WEIGHT}): {len(edges)}")

# Confidence breakdown
for conf in ("high", "medium", "low"):
    n_conf = sum(1 for e in edges if e["confidence"] == conf)
    print(f"    {conf:8s}: {n_conf} edges")


# ---------------------------------------------------------------------------
# 9. Assemble and save reference graph JSON
# ---------------------------------------------------------------------------

print(f"\n[export_v1] === Saving reference graph ===")

graph = {
    "_meta": {
        "description"          : "Working v1 scGPT-derived reference lineage graph for GSE230659.",
        "version"              : VERSION,
        "status"               : "silver_standard_working",
        "graph_method"         : "consecutive_timepoint_kNN_transition_counting",
        "knn_k"                : KNN_K,
        "min_edge_weight"      : MIN_EDGE_WEIGHT,
        "chunk_size"           : CHUNK_SIZE,
        "n_timepoint_pairs"    : len(tp_pairs),
        "embedding_source"     : "X_scGPT (scgpt_whole_human pretrained model)",
        "leiden_resolution"    : 0.5,
        "cluster_map_version"  : map_version,
        "n_cells_total"        : int(adata.n_obs),
        "n_states"             : n_states,
        "n_edges"              : len(edges),
        "created_from"         : str(ANNOTATED_H5AD),
        "policy": [
            "This graph is a silver-standard working asset, not final biology.",
            "Edges with confidence='low' involve at least one 'uncertain' state.",
            "Edges with confidence='medium' involve at least one 'merge_review' state.",
            "All weights are row-normalised transition probabilities per source state.",
            "Self-loops are excluded. Edges below MIN_EDGE_WEIGHT are pruned.",
            "To produce v2: revise cluster_map_provisional.json, rerun annotate_provisional.py then export_v1.py with VERSION='v2'.",
        ],
    },
    "nodes": nodes,
    "edges": edges,
}

with open(GRAPH_JSON, "w", encoding="utf-8") as f:
    json.dump(graph, f, indent=2)
print(f"  Saved {GRAPH_JSON.name}  ({GRAPH_JSON.stat().st_size / 1e3:.0f} KB)")

# Flat edge table CSV
edges_df = pd.DataFrame(edge_rows)
edges_df.to_csv(GRAPH_EDGE_CSV, index=False)
print(f"  Saved {GRAPH_EDGE_CSV.name}  ({len(edges_df)} edges)")


# ---------------------------------------------------------------------------
# 10. Export provenance JSON
# ---------------------------------------------------------------------------

total_time = time.time() - t0

export_meta = {
    "run_timestamp"              : time.strftime("%Y-%m-%dT%H:%M:%S"),
    "script"                     : "benchmark/scgpt/export_v1.py",
    "version"                    : VERSION,
    "status"                     : "silver_standard_working",
    "input_h5ad"                 : str(ANNOTATED_H5AD),
    "cluster_map"                : str(MAP_PATH),
    "cluster_map_version"        : map_version,
    "n_cells"                    : int(adata.n_obs),
    "n_states"                   : n_states,
    "states"                     : all_states,
    "graph_method"               : "consecutive_timepoint_kNN_transition_counting",
    "knn_k"                      : KNN_K,
    "min_edge_weight"            : MIN_EDGE_WEIGHT,
    "n_timepoint_pairs"          : len(tp_pairs),
    "n_graph_nodes"              : len(nodes),
    "n_graph_edges"              : len(edges),
    "n_edges_by_confidence"      : {
        conf: sum(1 for e in edges if e["confidence"] == conf)
        for conf in ("high", "medium", "low")
    },
    "confirmed_states"           : [ps for ps, s in ps_to_status.items() if s == "confirmed"],
    "merge_review_states"        : [ps for ps, s in ps_to_status.items() if s == "merge_review"],
    "uncertain_states"           : [ps for ps, s in ps_to_status.items() if s == "uncertain"],
    "output_files": {
        "label_tsv"    : str(LABEL_TSV),
        "meta_tsv"     : str(META_TSV),
        "graph_json"   : str(GRAPH_JSON),
        "graph_csv"    : str(GRAPH_EDGE_CSV),
    },
    "graph_runtime_sec"          : round(graph_runtime, 1),
    "total_runtime_sec"          : round(total_time, 1),
    "total_runtime_min"          : round(total_time / 60, 2),
}

with open(EXPORT_META, "w", encoding="utf-8") as f:
    json.dump(export_meta, f, indent=2)
print(f"  Saved {EXPORT_META.name}")


# ---------------------------------------------------------------------------
# 11. Final summary
# ---------------------------------------------------------------------------

print(f"\n[export_v1] === Export complete ===")
print(f"  Version          : {VERSION}  (silver-standard working)")
print(f"  Cells labelled   : {adata.n_obs:,}")
print(f"  States           : {n_states}")
print(f"  Graph nodes      : {len(nodes)}")
print(f"  Graph edges      : {len(edges)}")
high = export_meta['n_edges_by_confidence']['high']
med  = export_meta['n_edges_by_confidence']['medium']
low  = export_meta['n_edges_by_confidence']['low']
print(f"    high confidence: {high}")
print(f"    medium conf.   : {med}")
print(f"    low confidence : {low}  (involve uncertain states)")
print(f"  Graph runtime    : {graph_runtime:.1f}s")
print(f"  Total runtime    : {total_time:.1f}s  ({total_time/60:.1f} min)")
print(f"\n  Files written to {DATASET_DIR}:")
for fp in [LABEL_TSV, META_TSV, GRAPH_JSON, GRAPH_EDGE_CSV]:
    print(f"    {fp.name:<50}  {fp.stat().st_size / 1e3:>8.1f} KB")
print(f"\n  Provenance: {EXPORT_META}")

print(f"""
[export_v1] PASSED ✓

  To use v1 assets in WOT configuration (manual step, deferred):
    lineage:
      cell_state_key: scgpt_pseudostate_v1
      reference_graph_path: benchmark/datasets/scgpt_reference_graph_v1.json

  To produce a revised v2 later:
    1. Edit benchmark/scgpt/cluster_map_provisional.json
       (bump mapping_version to 'v2')
    2. python benchmark\\scgpt\\annotate_provisional.py
    3. Set VERSION = 'v2' in this script and rerun
""")
