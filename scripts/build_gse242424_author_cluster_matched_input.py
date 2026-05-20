#!/usr/bin/env python3
"""
Build the GSE242424 author-cluster matched HVG2000 subset.

This script attaches the author's ATAC-to-RNA transferred cluster labels to
the local GSE242424 HVG2000 h5ad by matching sample/day plus barcode core.
Cells without an exact author-label match are excluded from the output h5ad.

Reads:
  benchmark/inputs/gse242424_hvg2000/GSE242424_HVG2000_benchmark_input.h5ad

Downloads small public author annotation tables:
  - ATAC-to-RNA cluster transfer table
  - old-cluster to paper-cluster conversion table
  - paper cluster description table

Writes:
  benchmark/inputs/gse242424_author_cluster_matched/
    GSE242424_author_cluster_matched_HVG2000_benchmark_input.h5ad
    gse242424_author_cluster_matched_cells.tsv
    gse242424_author_cluster_match_summary.json
    gse242424_author_cluster_match_by_time.csv
    gse242424_author_cluster_counts.csv
    source_tables/*

Revision history:
  v1  2026-05-16  Initial author-cluster matched subset builder.
"""

from __future__ import annotations

import argparse
import io
import json
import os
import re
import time as _time
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd


DATASET_ID = "GSE242424"
TIMESTAMP = datetime.now().strftime("%Y%m%d_%H%M")
T0 = _time.time()

TRANSFER_URL = (
    "https://raw.githubusercontent.com/kundajelab/scATAC-reprog/master/"
    "src/analysis/20200828_RNA_Seurat/sessions/20210725_n59378/"
    "atac.20210717_n62599.cluster.transfer.tsv"
)
CLUSTER_CONVERSION_URL = (
    "https://raw.githubusercontent.com/kundajelab/scATAC-reprog/master/"
    "src/figures_factory/configs/cluster.tsv"
)
CLUSTERS_URL = "https://zenodo.org/record/8313962/files/clusters.tsv?download=1"

COARSE_MILESTONE_BY_CLUSTER_ID = {
    1: "fibroblast",
    2: "fibroblast_like_stalled",
    3: "fibroblast_like_stalled",
    4: "fibroblast_like_stalled",
    5: "fibroblast_like_stalled",
    6: "keratinocyte_like",
    7: "hOSK",
    8: "xOSK",
    9: "partial_intermediate",
    10: "partially_reprogrammed",
    11: "primary_intermediate",
    12: "primary_intermediate",
    13: "pre_iPSC",
    14: "pre_iPSC",
    15: "iPSC",
}


def find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        root = Path(env)
        if root.exists():
            return root
        raise FileNotFoundError(f"TRAJ_PROJECT_ROOT={env!r} does not exist.")
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "benchmark").exists() and (candidate / "scripts").exists():
            return candidate
    return here.parent


PROJECT_ROOT = find_project_root()
DEFAULT_INPUT_H5AD = (
    PROJECT_ROOT
    / "benchmark"
    / "inputs"
    / "gse242424_hvg2000"
    / "GSE242424_HVG2000_benchmark_input.h5ad"
)
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT / "benchmark" / "inputs" / "gse242424_author_cluster_matched"
)


def download_text(url: str, timeout_s: int = 60) -> str:
    with urllib.request.urlopen(url, timeout=timeout_s) as handle:
        return handle.read().decode("utf-8")


def read_public_tsv(url: str, output_path: Path) -> pd.DataFrame:
    text = download_text(url)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(text, encoding="utf-8")
    return pd.read_csv(io.StringIO(text), sep="\t")


def barcode_core(series: pd.Series) -> pd.Series:
    return series.astype(str).str.replace(r"-\d+$", "", regex=True)


def as_jsonable(value: Any) -> Any:
    if hasattr(value, "item"):
        return value.item()
    if isinstance(value, dict):
        return {str(k): as_jsonable(v) for k, v in value.items()}
    if isinstance(value, list):
        return [as_jsonable(v) for v in value]
    return value


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build GSE242424 author-cluster matched HVG2000 subset."
    )
    parser.add_argument(
        "--input-h5ad",
        type=Path,
        default=DEFAULT_INPUT_H5AD,
        help="Local GSE242424 HVG2000 h5ad to filter.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for the matched subset h5ad and reports.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting existing output files.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_h5ad = args.input_h5ad
    output_dir = args.output_dir
    source_dir = output_dir / "source_tables"

    if not input_h5ad.exists():
        raise FileNotFoundError(f"Input h5ad not found: {input_h5ad}")

    output_h5ad = (
        output_dir / "GSE242424_author_cluster_matched_HVG2000_benchmark_input.h5ad"
    )
    matched_tsv = output_dir / "gse242424_author_cluster_matched_cells.tsv"
    summary_json = output_dir / "gse242424_author_cluster_match_summary.json"
    by_time_csv = output_dir / "gse242424_author_cluster_match_by_time.csv"
    cluster_counts_csv = output_dir / "gse242424_author_cluster_counts.csv"

    if output_h5ad.exists() and not args.overwrite:
        raise FileExistsError(
            f"{output_h5ad} already exists. Pass --overwrite to rebuild."
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{TIMESTAMP}] reading author annotation tables")
    transfer = read_public_tsv(
        TRANSFER_URL, source_dir / "author_atac_to_rna_cluster_transfer.tsv"
    )
    conversion = read_public_tsv(
        CLUSTER_CONVERSION_URL, source_dir / "author_cluster_conversion.tsv"
    )
    clusters = read_public_tsv(CLUSTERS_URL, source_dir / "author_clusters.tsv")

    required_transfer = {"barcode", "sample", "atac_cluster"}
    missing = required_transfer - set(transfer.columns)
    if missing:
        raise ValueError(f"Transfer table missing columns: {sorted(missing)}")

    conversion_map = dict(
        zip(conversion["cluster"].astype(str), conversion["new_cluster"].astype(int))
    )
    cluster_name_map = dict(
        zip(clusters["id"].astype(int), clusters["description"].astype(str))
    )

    transfer = transfer.copy()
    transfer["author_old_atac_cluster"] = transfer["atac_cluster"].astype(int)
    transfer["author_cluster_id"] = (
        transfer["atac_cluster"].astype(str).map(conversion_map).astype(int)
    )
    transfer["author_cluster_label"] = transfer["author_cluster_id"].map(
        cluster_name_map
    )
    transfer["final_milestone_label_coarse"] = transfer["author_cluster_id"].map(
        COARSE_MILESTONE_BY_CLUSTER_ID
    )
    transfer["author_barcode"] = transfer["barcode"].astype(str)
    transfer["author_sample"] = transfer["sample"].astype(str)
    transfer["barcode_core"] = barcode_core(transfer["author_barcode"])
    transfer["join_key"] = transfer["author_sample"] + "_" + transfer["barcode_core"]

    if transfer["join_key"].duplicated().any():
        examples = transfer.loc[transfer["join_key"].duplicated(), "join_key"].head()
        raise ValueError(f"Duplicate author join keys, examples: {examples.tolist()}")
    if transfer["author_cluster_label"].isna().any():
        raise ValueError("Some author_cluster_id values lack cluster descriptions.")
    if transfer["final_milestone_label_coarse"].isna().any():
        raise ValueError("Some author_cluster_id values lack coarse milestone mapping.")

    print(f"[{TIMESTAMP}] reading local obs from {input_h5ad}")
    adata_backed = ad.read_h5ad(input_h5ad, backed="r")
    obs = adata_backed.obs[["cell_id", "sample_id", "time_label"]].copy()
    obs["obs_name"] = obs.index.astype(str)
    obs["obs_pos"] = range(obs.shape[0])
    obs["local_barcode"] = obs["cell_id"].astype(str).map(lambda x: x.split("_")[-1])
    obs["barcode_core"] = barcode_core(obs["local_barcode"])
    obs["join_key"] = obs["time_label"].astype(str) + "_" + obs["barcode_core"]

    if obs["join_key"].duplicated().any():
        examples = obs.loc[obs["join_key"].duplicated(), "join_key"].head()
        raise ValueError(f"Duplicate local join keys, examples: {examples.tolist()}")

    label_cols = [
        "join_key",
        "author_sample",
        "author_barcode",
        "author_old_atac_cluster",
        "author_cluster_id",
        "author_cluster_label",
        "final_milestone_label_coarse",
    ]
    merged = obs.merge(transfer[label_cols], on="join_key", how="left")
    matched_mask = merged["author_cluster_id"].notna().to_numpy()
    matched = merged.loc[matched_mask].copy()
    unmatched = merged.loc[~matched_mask].copy()

    print(f"[{TIMESTAMP}] local cells: {len(obs):,}")
    print(f"[{TIMESTAMP}] author transfer rows: {len(transfer):,}")
    print(f"[{TIMESTAMP}] matched local cells: {len(matched):,}")
    print(f"[{TIMESTAMP}] excluded unmatched local cells: {len(unmatched):,}")

    for col in ["author_old_atac_cluster", "author_cluster_id"]:
        matched[col] = matched[col].astype(int)

    matched_report_cols = [
        "cell_id",
        "sample_id",
        "time_label",
        "local_barcode",
        "barcode_core",
        "author_sample",
        "author_barcode",
        "author_old_atac_cluster",
        "author_cluster_id",
        "author_cluster_label",
        "final_milestone_label_coarse",
    ]
    matched[matched_report_cols].to_csv(matched_tsv, sep="\t", index=False)

    by_time = (
        merged.assign(matched=merged["author_cluster_id"].notna())
        .groupby("time_label", observed=False)["matched"]
        .agg(matched_cells="sum", local_cells="count")
        .reset_index()
    )
    by_time["excluded_unmatched_cells"] = (
        by_time["local_cells"] - by_time["matched_cells"]
    )
    by_time["match_rate"] = by_time["matched_cells"] / by_time["local_cells"]
    by_time.to_csv(by_time_csv, index=False)

    cluster_counts = (
        matched.groupby(
            [
                "author_cluster_id",
                "author_cluster_label",
                "final_milestone_label_coarse",
            ],
            observed=False,
        )
        .size()
        .reset_index(name="n_matched_cells")
        .sort_values("author_cluster_id")
    )
    cluster_counts.to_csv(cluster_counts_csv, index=False)

    print(f"[{TIMESTAMP}] writing matched h5ad: {output_h5ad}")
    subset = adata_backed[matched_mask, :].to_memory()
    adata_backed.file.close()

    subset.obs = subset.obs.copy()
    label_values = matched.set_index("obs_name")
    for col in [
        "local_barcode",
        "barcode_core",
        "author_sample",
        "author_barcode",
        "author_old_atac_cluster",
        "author_cluster_id",
        "author_cluster_label",
        "final_milestone_label_coarse",
    ]:
        subset.obs[col] = label_values.loc[subset.obs.index.astype(str), col].to_numpy()
    subset.obs["author_label_match_status"] = "matched_author_transfer"
    subset.obs["author_label_source"] = (
        "kundajelab/scATAC-reprog ATAC-to-RNA transfer table"
    )
    subset.obs["author_label_confidence_note"] = (
        "hard transferred cluster label; per-cell transfer scores not used here"
    )

    subset.uns["gse242424_author_cluster_matching"] = {
        "dataset_id": DATASET_ID,
        "timestamp": TIMESTAMP,
        "source_h5ad": str(input_h5ad),
        "output_h5ad": str(output_h5ad),
        "matching_rule": "time_label/sample plus barcode core after removing trailing -N",
        "local_cells": int(len(obs)),
        "author_transfer_rows": int(len(transfer)),
        "matched_local_cells": int(len(matched)),
        "excluded_unmatched_local_cells": int(len(unmatched)),
        "transfer_rows_not_in_local": int(
            (~transfer["join_key"].isin(set(obs["join_key"]))).sum()
        ),
        "source_urls": {
            "transfer_table": TRANSFER_URL,
            "cluster_conversion": CLUSTER_CONVERSION_URL,
            "clusters": CLUSTERS_URL,
        },
        "coarse_milestone_mapping": {
            str(k): v for k, v in COARSE_MILESTONE_BY_CLUSTER_ID.items()
        },
        "note": (
            "Original full GSE242424 HVG2000 input is unchanged; this h5ad "
            "contains only local cells with exact author transfer-table matches."
        ),
    }
    subset.write_h5ad(output_h5ad)

    summary = {
        "dataset_id": DATASET_ID,
        "timestamp": TIMESTAMP,
        "input_h5ad": str(input_h5ad),
        "output_dir": str(output_dir),
        "output_h5ad": str(output_h5ad),
        "matched_cells_tsv": str(matched_tsv),
        "local_cells": int(len(obs)),
        "author_transfer_rows": int(len(transfer)),
        "matched_local_cells": int(len(matched)),
        "excluded_unmatched_local_cells": int(len(unmatched)),
        "match_rate_local": float(len(matched) / len(obs)),
        "transfer_rows_not_in_local": int(
            (~transfer["join_key"].isin(set(obs["join_key"]))).sum()
        ),
        "barcode_matching_rule": (
            "local_join_key=time_label + '_' + barcode_core; "
            "author_join_key=sample + '_' + barcode_core; "
            "barcode_core removes trailing -N suffix"
        ),
        "source_urls": {
            "transfer_table": TRANSFER_URL,
            "cluster_conversion": CLUSTER_CONVERSION_URL,
            "clusters": CLUSTERS_URL,
        },
        "by_time": by_time.to_dict(orient="records"),
        "cluster_counts": cluster_counts.to_dict(orient="records"),
        "runtime_s": round(_time.time() - T0, 2),
    }
    with summary_json.open("w", encoding="utf-8") as handle:
        json.dump(as_jsonable(summary), handle, indent=2, ensure_ascii=True)

    print(f"[{TIMESTAMP}] wrote {matched_tsv}")
    print(f"[{TIMESTAMP}] wrote {summary_json}")
    print(f"[{TIMESTAMP}] done in {summary['runtime_s']} s")


if __name__ == "__main__":
    main()
