"""Transfer full-gene milestone labels onto an HVG benchmark h5ad.

The intended workflow is:

1. Build a full-gene h5ad for marker scoring.
2. Run build_marker_seed_labels.py on the full-gene h5ad.
3. Transfer only milestone obs columns back to the HVG2000 h5ad used for model
   training/evaluation.

Expression, embeddings, and var in the HVG h5ad are preserved.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path


MILESTONE_OBS_COLUMNS = [
    "milestone_marker_label",
    "milestone_marker_score",
    "milestone_marker_margin",
    "milestone_embedding_label",
    "milestone_embedding_confidence",
    "milestone_classifier_label",
    "milestone_classifier_confidence",
    "consensus_milestone_label",
    "consensus_confidence",
]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Transfer milestone obs columns from full-gene h5ad to HVG h5ad."
    )
    parser.add_argument("--full-gene-h5ad", required=True)
    parser.add_argument("--hvg-h5ad", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument(
        "--alignment-key",
        default="auto",
        choices=["auto", "obs_names", "cell_id"],
        help="How to align cells. auto tries obs_names first, then cell_id.",
    )
    parser.add_argument("--dry-run", action="store_true")
    return parser


def _choose_alignment(full, hvg, mode: str) -> tuple[str, dict | None]:
    if mode in ("auto", "obs_names"):
        if list(full.obs_names) == list(hvg.obs_names):
            return "obs_names_ordered", None
        if set(full.obs_names) == set(hvg.obs_names):
            return "obs_names_reindex", {str(x): i for i, x in enumerate(full.obs_names)}
        if mode == "obs_names":
            raise ValueError("obs_names alignment requested but obs_names do not match.")

    if mode in ("auto", "cell_id"):
        if "cell_id" not in full.obs.columns or "cell_id" not in hvg.obs.columns:
            raise ValueError("cell_id alignment requested but cell_id is missing.")
        full_ids = full.obs["cell_id"].astype(str)
        hvg_ids = hvg.obs["cell_id"].astype(str)
        if full_ids.duplicated().any():
            raise ValueError("full-gene obs['cell_id'] contains duplicates.")
        if hvg_ids.duplicated().any():
            raise ValueError("HVG obs['cell_id'] contains duplicates.")
        if set(full_ids) != set(hvg_ids):
            missing_from_full = sorted(set(hvg_ids) - set(full_ids))[:10]
            missing_from_hvg = sorted(set(full_ids) - set(hvg_ids))[:10]
            raise ValueError(
                "cell_id sets do not match. "
                f"missing_from_full(first10)={missing_from_full}; "
                f"missing_from_hvg(first10)={missing_from_hvg}"
            )
        return "cell_id_reindex", {cid: i for i, cid in enumerate(full_ids)}

    raise ValueError("Could not determine a valid alignment mode.")


def _aligned_series(full_obs, hvg_obs, col: str, alignment: str, index_map: dict | None):
    if alignment == "obs_names_ordered":
        values = full_obs[col].to_numpy()
    elif alignment == "obs_names_reindex":
        positions = [index_map[str(idx)] for idx in hvg_obs.index]
        values = full_obs[col].iloc[positions].to_numpy()
    elif alignment == "cell_id_reindex":
        positions = [index_map[str(cid)] for cid in hvg_obs["cell_id"].astype(str)]
        values = full_obs[col].iloc[positions].to_numpy()
    else:
        raise ValueError(f"Unsupported alignment mode: {alignment}")
    return values


def run_transfer(args: argparse.Namespace) -> dict:
    import anndata as ad
    import pandas as pd

    full_path = Path(args.full_gene_h5ad)
    hvg_path = Path(args.hvg_h5ad)
    out_path = Path(args.output_h5ad)
    out_dir = Path(args.output_dir) if args.output_dir else out_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"[transfer_milestone] loading full-gene labels: {full_path}")
    full = ad.read_h5ad(full_path)
    print(f"[transfer_milestone] loading HVG h5ad: {hvg_path}")
    hvg = ad.read_h5ad(hvg_path)

    missing = [c for c in MILESTONE_OBS_COLUMNS if c not in full.obs.columns]
    if missing:
        raise ValueError(f"Full-gene h5ad missing milestone columns: {missing}")

    if full.n_obs != hvg.n_obs:
        raise ValueError(f"Cell count mismatch: full={full.n_obs}, hvg={hvg.n_obs}")

    alignment, index_map = _choose_alignment(full, hvg, args.alignment_key)
    print(f"[transfer_milestone] alignment: {alignment}")

    for col in MILESTONE_OBS_COLUMNS:
        hvg.obs[col] = _aligned_series(full.obs, hvg.obs, col, alignment, index_map)

    hvg.uns["milestone_label_transfer"] = {
        "script": "benchmark/annotation/transfer_milestone_labels_to_hvg.py",
        "dataset_id": args.dataset_id,
        "status": "completed",
        "source_full_gene_h5ad": str(full_path),
        "target_hvg_h5ad": str(hvg_path),
        "output_h5ad": str(out_path),
        "alignment": alignment,
        "columns_transferred": MILESTONE_OBS_COLUMNS,
        "expression_matrix_policy": "preserved_hvg2000_expression",
        "annotation_policy": "full_gene_marker_labels_transferred_to_hvg",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    counts = (
        hvg.obs["consensus_milestone_label"]
        .astype(str)
        .value_counts()
        .rename_axis("label")
        .reset_index(name="count")
    )

    metadata = dict(hvg.uns["milestone_label_transfer"])
    metadata.update(
        {
            "n_cells": int(hvg.n_obs),
            "n_genes_hvg": int(hvg.n_vars),
            "label_counts": dict(zip(counts["label"].astype(str), counts["count"].astype(int))),
            "full_gene_shape": [int(full.n_obs), int(full.n_vars)],
            "hvg_shape": [int(hvg.n_obs), int(hvg.n_vars)],
        }
    )

    if args.dry_run:
        print("[transfer_milestone] dry-run: no files written")
    else:
        print(f"[transfer_milestone] writing output h5ad: {out_path}")
        hvg.write_h5ad(out_path)
        counts.to_csv(out_dir / "transferred_milestone_label_counts.csv", index=False)
        with open(out_dir / "milestone_label_transfer_metadata.json", "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)

    full.file.close() if getattr(full, "isbacked", False) else None
    hvg.file.close() if getattr(hvg, "isbacked", False) else None
    return metadata


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        meta = run_transfer(args)
    except Exception as exc:
        print(f"[transfer_milestone] ERROR: {exc}")
        return 1
    print(json.dumps(meta, indent=2)[:2000])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
