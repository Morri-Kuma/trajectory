"""
annotate_provisional.py
Add provisional pseudo-state annotation layer on top of the stable full embedding.

Purpose
-------
Read the stable embedding artifact (adata_scgpt_full.h5ad) and the editable
provisional cluster mapping (cluster_map_provisional.json), join the
interpretation onto the AnnData, produce QC summaries and UMAP figures, and
save the annotated AnnData.

This script requires NO GPU and NO recomputation of X_scGPT.
It is designed to be rerun quickly (~1 minute) whenever cluster_map_provisional.json
is edited — for example after a merge decision, a status revision, or adding a new
Leiden resolution column to adata_scgpt_full.h5ad.

What this script does NOT do
-----------------------------
  - Does not rerun embedding (X_scGPT is read from the existing h5ad)
  - Does not rerun Leiden (leiden_scgpt_res0.5 is read from the existing h5ad)
  - Does not modify adata_scgpt_full.h5ad (writes a separate annotated copy)
  - Does not hard-freeze pseudo-states (all assignments remain provisional)
  - Does not touch the active WOT / CellRank2 benchmark configuration

Environment
-----------
  conda activate scgpt_env
  cd C:\\Users\\37620\\trajectory
  python benchmark\\scgpt\\annotate_provisional.py

Inputs
------
  benchmark/results/scgpt/full/adata_scgpt_full.h5ad
  benchmark/scgpt/cluster_map_provisional.json

Outputs  (benchmark/results/scgpt/full/)
-----------------------------------------
  adata_scgpt_annotated.h5ad     full AnnData + 3 new obs columns
  cluster_qc.csv                 per-cluster QC table
  annotation_metadata.json       provenance record for this annotation run
  umap_by_pseudostate.png        UMAP coloured by scgpt_pseudostate_provisional
  umap_by_family.png             UMAP coloured by scgpt_state_family
  umap_by_status.png             UMAP coloured by scgpt_state_status
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
    known = Path(r"C:\Users\37620\trajectory")
    if known.exists() and (known / "data").exists():
        return known
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "data").exists():
            return candidate
    return here.parent


PROJECT_ROOT = _find_project_root()
print(f"[annotate] Project root: {PROJECT_ROOT}")


# ---------------------------------------------------------------------------
# 1. Paths
# ---------------------------------------------------------------------------

FULL_H5AD_PATH   = PROJECT_ROOT / "benchmark" / "results" / "scgpt" / "full" / "adata_scgpt_full.h5ad"
MAP_PATH         = PROJECT_ROOT / "benchmark" / "scgpt" / "cluster_map_provisional.json"
OUTPUT_DIR       = PROJECT_ROOT / "benchmark" / "results" / "scgpt" / "full"

LEIDEN_KEY       = "leiden_scgpt_res0.5"     # source column in obs
TIME_KEY         = "abs_day"                  # timepoint column

# Derived obs columns to add
PS_KEY           = "scgpt_pseudostate_provisional"
FAMILY_KEY       = "scgpt_state_family"
STATUS_KEY       = "scgpt_state_status"


# ---------------------------------------------------------------------------
# 2. Pre-flight checks
# ---------------------------------------------------------------------------

print("\n[annotate] === Pre-flight checks ===")

if not FULL_H5AD_PATH.exists():
    raise FileNotFoundError(
        f"Stable embedding artifact not found:\n  {FULL_H5AD_PATH}\n"
        "Run benchmark/scgpt/full_embed.py first."
    )
print(f"  Embedding artifact : {FULL_H5AD_PATH}  ✓")

if not MAP_PATH.exists():
    raise FileNotFoundError(
        f"Cluster map not found:\n  {MAP_PATH}\n"
        "Ensure benchmark/scgpt/cluster_map_provisional.json exists."
    )
print(f"  Cluster map        : {MAP_PATH}  ✓")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"  Output dir         : {OUTPUT_DIR}  ✓")


# ---------------------------------------------------------------------------
# 3. Load inputs
# ---------------------------------------------------------------------------

print(f"\n[annotate] === Loading inputs ===")

t0 = time.time()
adata = sc.read_h5ad(FULL_H5AD_PATH)
print(f"  Loaded h5ad in {time.time()-t0:.1f}s  →  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

if LEIDEN_KEY not in adata.obs.columns:
    raise KeyError(
        f"Expected Leiden column '{LEIDEN_KEY}' not found in obs.\n"
        f"Available obs columns: {list(adata.obs.columns)}"
    )
if "X_scGPT" not in adata.obsm:
    raise KeyError("X_scGPT not found in adata.obsm — is this the correct h5ad?")
if "X_umap" not in adata.obsm:
    raise KeyError("X_umap not found in adata.obsm — UMAP not yet computed in full_embed.py.")

leiden_clusters_in_data = set(adata.obs[LEIDEN_KEY].astype(str).unique())
print(f"  Leiden clusters in obs: {sorted(leiden_clusters_in_data, key=int)}")

with open(MAP_PATH, "r", encoding="utf-8") as f:
    cluster_map = json.load(f)

clusters_in_map = set(cluster_map["clusters"].keys())
map_version     = cluster_map["_meta"].get("mapping_version", "unknown")
print(f"  Map version: {map_version}  |  Clusters in map: {sorted(clusters_in_map, key=int)}")

# Consistency check: every cluster in the data must exist in the map
missing_from_map = leiden_clusters_in_data - clusters_in_map
extra_in_map     = clusters_in_map - leiden_clusters_in_data
if missing_from_map:
    raise ValueError(
        f"Clusters present in obs but missing from cluster_map_provisional.json:\n"
        f"  {sorted(missing_from_map, key=int)}\n"
        "Add entries for these clusters to the mapping file."
    )
if extra_in_map:
    print(f"  NOTE: map contains {len(extra_in_map)} cluster(s) not present in obs "
          f"(will be ignored): {sorted(extra_in_map, key=int)}")

print("  Consistency check passed  ✓")


# ---------------------------------------------------------------------------
# 4. Build lookup tables from the mapping file
# ---------------------------------------------------------------------------

ps_lookup     = {}  # cluster_id (str) → pseudostate label
family_lookup = {}
status_lookup = {}

for cid, entry in cluster_map["clusters"].items():
    ps_lookup[cid]     = entry["pseudostate"]
    family_lookup[cid] = entry["family"]
    status_lookup[cid] = entry["status"]


# ---------------------------------------------------------------------------
# 5. Add derived obs columns
# ---------------------------------------------------------------------------

print(f"\n[annotate] === Adding annotation columns ===")

leiden_str = adata.obs[LEIDEN_KEY].astype(str)

adata.obs[PS_KEY]     = leiden_str.map(ps_lookup).astype("category")
adata.obs[FAMILY_KEY] = leiden_str.map(family_lookup).astype("category")
adata.obs[STATUS_KEY] = leiden_str.map(status_lookup).astype("category")

# Verify no unmapped cells
n_unmapped = adata.obs[PS_KEY].isna().sum()
if n_unmapped > 0:
    raise ValueError(
        f"{n_unmapped} cells have no pseudostate mapping — "
        "check that all Leiden cluster IDs are present in the map file."
    )

print(f"  {PS_KEY}  : {adata.obs[PS_KEY].nunique()} unique labels")
print(f"  {FAMILY_KEY} : {adata.obs[FAMILY_KEY].nunique()} unique families")
print(f"  {STATUS_KEY} : {adata.obs[STATUS_KEY].nunique()} unique statuses")
print(f"  All {adata.n_obs:,} cells mapped  ✓")


# ---------------------------------------------------------------------------
# 6. Compute per-cluster QC table
# ---------------------------------------------------------------------------

print(f"\n[annotate] === Computing cluster QC ===")

n_timepoints = adata.obs[TIME_KEY].nunique()
records = []

for cid in sorted(leiden_clusters_in_data, key=int):
    mask  = adata.obs[LEIDEN_KEY].astype(str) == cid
    sub   = adata.obs[mask]
    n     = mask.sum()

    tp_counts = sub[TIME_KEY].value_counts().sort_index()
    dom_tp    = float(tp_counts.idxmax())
    dom_frac  = float(tp_counts.max() / n)

    probs = (tp_counts / n).values
    probs_nz = probs[probs > 0]
    norm_entropy = float(scipy_entropy(probs_nz) / np.log(n_timepoints))

    entry = cluster_map["clusters"][cid]

    records.append({
        "leiden_cluster"   : cid,
        "pseudostate"      : entry["pseudostate"],
        "family"           : entry["family"],
        "status"           : entry["status"],
        "n_cells"          : int(n),
        "frac_of_total"    : round(n / adata.n_obs, 4),
        "dominant_timepoint": dom_tp,
        "dominant_frac"    : round(dom_frac, 3),
        "norm_entropy"     : round(norm_entropy, 3),
        "notes"            : entry.get("notes", ""),
    })

qc_df = pd.DataFrame(records).set_index("leiden_cluster")
qc_df = qc_df.sort_values("n_cells", ascending=False)

print(f"\n  {'Cluster':>7} | {'Pseudostate':>7} | {'N':>6} | {'entropy':>7} | {'dom_tp':>6} | {'dom_frac':>8} | status")
print("  " + "-" * 75)
for idx, row in qc_df.iterrows():
    print(f"  {idx:>7} | {row['pseudostate']:>7} | {row['n_cells']:>6,} | "
          f"{row['norm_entropy']:>7.3f} | {row['dominant_timepoint']:>6.1f} | "
          f"{row['dominant_frac']:>8.3f} | {row['status']}")

csv_out = OUTPUT_DIR / "cluster_qc.csv"
qc_df.to_csv(csv_out)
print(f"\n  Saved cluster_qc.csv")


# ---------------------------------------------------------------------------
# 7. Save annotated AnnData
# ---------------------------------------------------------------------------

print(f"\n[annotate] === Saving adata_scgpt_annotated.h5ad ===")
annotated_out = OUTPUT_DIR / "adata_scgpt_annotated.h5ad"
t_save = time.time()
adata.write_h5ad(annotated_out)
save_time = time.time() - t_save
file_mb = annotated_out.stat().st_size / 1e6
print(f"  Saved in {save_time:.1f}s  →  {annotated_out.name}  ({file_mb:.0f} MB)")


# ---------------------------------------------------------------------------
# 8. UMAP figures
# ---------------------------------------------------------------------------

print(f"\n[annotate] === Saving UMAP figures ===")
sc.settings.figdir = str(OUTPUT_DIR)
sc.settings.verbosity = 0

# Status palette: green = confirmed, orange = merge_review, red = uncertain
STATUS_PALETTE = {
    "confirmed"   : "#2ca02c",
    "merge_review": "#ff7f0e",
    "uncertain"   : "#d62728",
}

# Figure 1: by pseudostate
sc.pl.umap(
    adata,
    color  = PS_KEY,
    title  = f"scGPT full — provisional pseudo-states  (n={adata.n_obs:,})",
    show   = False,
    save   = "_by_pseudostate.png",
)
print("  Saved umap_by_pseudostate.png")

# Figure 2: by family
sc.pl.umap(
    adata,
    color  = FAMILY_KEY,
    title  = f"scGPT full — state family",
    show   = False,
    save   = "_by_family.png",
)
print("  Saved umap_by_family.png")

# Figure 3: by status (confirmed / merge_review / uncertain)
sc.pl.umap(
    adata,
    color   = STATUS_KEY,
    palette = STATUS_PALETTE,
    title   = f"scGPT full — cluster status",
    show    = False,
    save    = "_by_status.png",
)
print("  Saved umap_by_status.png")


# ---------------------------------------------------------------------------
# 9. Annotation metadata JSON
# ---------------------------------------------------------------------------

status_counts  = adata.obs[STATUS_KEY].value_counts().to_dict()
family_counts  = adata.obs[FAMILY_KEY].value_counts().to_dict()

annotation_metadata = {
    "run_timestamp"             : time.strftime("%Y-%m-%dT%H:%M:%S"),
    "script"                    : "benchmark/scgpt/annotate_provisional.py",
    "input_h5ad"                : str(FULL_H5AD_PATH),
    "cluster_map"               : str(MAP_PATH),
    "map_version"               : map_version,
    "map_status"                : cluster_map["_meta"]["mapping_status"],
    "leiden_key"                : LEIDEN_KEY,
    "n_cells"                   : int(adata.n_obs),
    "n_leiden_clusters"         : int(len(leiden_clusters_in_data)),
    "output_columns_added"      : [PS_KEY, FAMILY_KEY, STATUS_KEY],
    "status_cell_counts"        : {k: int(v) for k, v in status_counts.items()},
    "family_cell_counts"        : {k: int(v) for k, v in family_counts.items()},
    "confirmed_clusters"        : [cid for cid, e in cluster_map["clusters"].items()
                                   if e["status"] == "confirmed"],
    "merge_review_clusters"     : [cid for cid, e in cluster_map["clusters"].items()
                                   if e["status"] == "merge_review"],
    "uncertain_clusters"        : [cid for cid, e in cluster_map["clusters"].items()
                                   if e["status"] == "uncertain"],
    "output_h5ad"               : str(annotated_out),
    "output_h5ad_size_mb"       : round(file_mb, 1),
}

meta_out = OUTPUT_DIR / "annotation_metadata.json"
with open(meta_out, "w", encoding="utf-8") as f:
    json.dump(annotation_metadata, f, indent=2)
print(f"  Saved annotation_metadata.json")


# ---------------------------------------------------------------------------
# 10. Final summary
# ---------------------------------------------------------------------------

total_time = time.time() - t0
print(f"\n[annotate] === Annotation complete ===")
print(f"  Runtime              : {total_time:.1f}s  ({total_time/60:.1f} min)")
print(f"  Cells annotated      : {adata.n_obs:,}")
print(f"  Leiden clusters      : {len(leiden_clusters_in_data)}")
print(f"  Confirmed clusters   : {annotation_metadata['confirmed_clusters']}")
print(f"  Merge-review clusters: {annotation_metadata['merge_review_clusters']}")
print(f"  Uncertain clusters   : {annotation_metadata['uncertain_clusters']}")
print(f"\n  New obs columns:")
print(f"    {PS_KEY}  →  {adata.obs[PS_KEY].nunique()} unique values")
print(f"    {FAMILY_KEY}       →  {adata.obs[FAMILY_KEY].nunique()} unique values")
print(f"    {STATUS_KEY}       →  {adata.obs[STATUS_KEY].nunique()} unique values")
print(f"\n  Outputs → {OUTPUT_DIR}")
for fp in sorted(OUTPUT_DIR.iterdir()):
    if fp.is_file():
        print(f"    {fp.name:<45}  {fp.stat().st_size / 1e6:>7.1f} MB")

print(f"\n[annotate] PASSED ✓")
print(
    "\n  To revise cluster assignments:\n"
    "    1. Edit benchmark/scgpt/cluster_map_provisional.json\n"
    "    2. Rerun:  python benchmark\\scgpt\\annotate_provisional.py\n"
    "    3. No GPU needed; embedding is not recomputed.\n"
    "\n  To add a new Leiden resolution:\n"
    "    1. Load adata_scgpt_full.h5ad\n"
    "    2. Run sc.tl.leiden(..., key_added='leiden_scgpt_res0.3')\n"
    "    3. Save back to adata_scgpt_full.h5ad\n"
    "    4. Update LEIDEN_KEY in this script and update cluster_map_provisional.json"
)
