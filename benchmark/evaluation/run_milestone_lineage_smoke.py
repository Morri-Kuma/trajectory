"""run_milestone_lineage_smoke.py
==================================
Smoke-test helper: run lineage evaluation against a milestone provider config
without executing the full scNODE (or any other) training pipeline.

Smoke strategy
--------------
For a smoke test we need predicted state-transition outputs, but we have no
method's trained model at hand.  We therefore construct an **oracle predicted
STM** directly from the reference graph:

  - reference edges → STM weight 1.0 (row-normalised by out-degree)
  - non-reference pairs → weight 0.0
  - diagonal → weight 0.0 (no self-loops)

This gives AUROC=1.0 / Jaccard=1.0 by construction, which is exactly what we
want for a pipeline smoke test: it confirms that load_reference_graph(), the
metric functions, and the file-I/O all work correctly against a milestone
provider without any biological inference.

The same h5ad is used as the adata source for the correlation baseline.
"ambiguous" labelled cells are excluded before baseline computation; a
clear metadata note records the exclusion.

Output layout (never overlaps formal benchmark directories)
-----------------------------------------------------------
  benchmark/results/smoke/milestone_lineage/<run_id>/
    state_transition_matrix.csv   (oracle STM — NOT a method prediction)
    lineage_graph_edges.csv       (oracle edges)
    lineage_metrics.json          (written by eval_lineage.run_lineage_evaluation)
    baseline_state_transition_matrix.csv   (written by baseline inside eval_lineage)
    baseline_lineage_graph_edges.csv       (written by baseline inside eval_lineage)
    smoke_metadata.json           (smoke-run provenance)

Usage
-----
  python -m benchmark.evaluation.run_milestone_lineage_smoke \\
      benchmark/configs/scnode_gse230659_observed_milestone_consensus_v1_A_hvg2000.yaml

  python benchmark/evaluation/run_milestone_lineage_smoke.py \\
      benchmark/configs/scnode_gse230659_observed_milestone_consensus_v1_A_hvg2000.yaml \\
      [--output-root benchmark/results/smoke/milestone_lineage] \\
      [--dry-run]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _project_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in [here.parent, *here.parents]:
        if (candidate / "benchmark").exists() and (
            candidate / "benchmark" / "evaluation"
        ).exists():
            return candidate
    return here.parents[2]


def _load_yaml(path: Path) -> Dict[str, Any]:
    import yaml
    with open(path, encoding="utf-8-sig") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Expected mapping in {path}")
    return data


def _git_commit(root: Path) -> Optional[str]:
    try:
        result = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return None


def _build_oracle_stm(
    node_ids: List[str],
    reference_edges: set,
) -> "pd.DataFrame":
    """
    Build an oracle state-transition matrix from reference graph edges.

    Each reference edge (src, tgt) is given weight 1.0; all other pairs zero.
    Rows are source states, columns are target states.
    Rows are normalised by out-degree so values lie in [0, 1].

    This is NOT a method prediction — it is used only to exercise the
    eval_lineage pipeline in smoke mode.
    """
    import numpy as np
    import pandas as pd

    n = len(node_ids)
    mat = np.zeros((n, n), dtype=float)
    idx = {s: i for i, s in enumerate(node_ids)}

    for (src, tgt) in reference_edges:
        i, j = idx.get(src), idx.get(tgt)
        if i is not None and j is not None:
            mat[i, j] = 1.0

    # Row-normalise (avoid divide-by-zero for sink nodes).
    row_sums = mat.sum(axis=1, keepdims=True)
    nonzero = row_sums[:, 0] > 0
    mat[nonzero] = mat[nonzero] / row_sums[nonzero]

    return pd.DataFrame(mat, index=node_ids, columns=node_ids)


def _build_oracle_edges(stm: "pd.DataFrame") -> "pd.DataFrame":
    """Convert a weighted STM to a sparse edge DataFrame (all nonzero entries)."""
    import pandas as pd
    rows = []
    for src in stm.index:
        for tgt in stm.columns:
            w = stm.loc[src, tgt]
            if w > 0:
                rows.append({"source_state": src, "target_state": tgt, "weight": float(w)})
    return pd.DataFrame(rows, columns=["source_state", "target_state", "weight"])


def _read_h5ad_safe(h5ad_path: Path):
    """
    Load h5ad via h5py, skipping the uns group entirely.

    The milestone-annotated h5ad stores None in
    uns/milestone_annotation_step3/max_cells_used, which older anndata
    versions cannot decode (IORegistryError: null encoding).  Reading
    via h5py and constructing AnnData manually avoids that bug while
    still giving us a fully functional X / obs / var object.
    """
    import h5py
    import scipy.sparse as sp
    import anndata as ad
    import pandas as pd

    with h5py.File(h5ad_path, "r") as f:
        # -- obs --
        obs_group = f["obs"]
        index_name = obs_group.attrs.get("_index", "_index")
        obs_data = {}
        for col in obs_group.keys():
            item = obs_group[col]
            if isinstance(item, h5py.Dataset):
                arr = item[:]
                if arr.dtype.kind in ("S", "O"):
                    arr = arr.astype(str)
                obs_data[col] = arr
            elif isinstance(item, h5py.Group) and "categories" in item:
                cats = item["categories"][:].astype(str)
                codes = item["codes"][:]
                obs_data[col] = pd.Categorical.from_codes(codes, categories=cats)
        idx = obs_data.pop(index_name, None)
        obs_df = pd.DataFrame(obs_data)
        if idx is not None:
            obs_df.index = idx.astype(str) if hasattr(idx, "astype") else idx

        # -- var --
        var_group = f["var"]
        var_index = var_group["_index"][:].astype(str)
        var_df = pd.DataFrame(index=var_index)

        # -- X (CSR sparse) --
        Xg = f["X"]
        n_obs, n_var = obs_df.shape[0], len(var_index)
        if isinstance(Xg, h5py.Group):
            X = sp.csr_matrix(
                (Xg["data"][:], Xg["indices"][:], Xg["indptr"][:]),
                shape=(n_obs, n_var),
            )
        else:
            X = Xg[:]

    return ad.AnnData(X=X, obs=obs_df, var=var_df)


def _load_adata_filtered(
    h5ad_path: Path,
    state_key: str,
    excluded_labels: List[str],
) -> "tuple[anndata.AnnData, dict]":
    """
    Load h5ad (safely via h5py to avoid null-encoding issues) and drop
    cells whose state_key label is in excluded_labels.

    Returns the filtered AnnData and an exclusion summary dict.
    """
    adata = _read_h5ad_safe(h5ad_path)
    labels = adata.obs[state_key].astype(str)
    n_before = adata.n_obs
    keep_mask = ~labels.isin(excluded_labels)
    n_excluded = int((~keep_mask).sum())
    adata_filtered = adata[keep_mask].copy()
    summary = {
        "n_cells_before_exclusion": n_before,
        "n_cells_excluded": n_excluded,
        "n_cells_after_exclusion": adata_filtered.n_obs,
        "excluded_labels": excluded_labels,
        "exclusion_reason": (
            "label not present in marker-defined milestone graph"
            if excluded_labels else "no exclusions applied"
        ),
    }
    return adata_filtered, summary


# ---------------------------------------------------------------------------
# Main smoke runner
# ---------------------------------------------------------------------------

def run_smoke(
    config_path: Path,
    output_root: Path,
    dry_run: bool = False,
) -> Dict[str, Any]:
    """
    Run a lineage fidelity smoke test for one milestone config.

    Returns the lineage metrics dict (or a dry-run summary dict).
    """
    root = _project_root()

    # -- Load config --
    cfg = _load_yaml(config_path)

    # Guard: refuse to write into formal benchmark directories.
    formal = cfg.get("formal_benchmark", True)  # default-safe: assume formal unless stated
    if formal:
        raise RuntimeError(
            f"Config {config_path.name} has formal_benchmark=true (or unset). "
            "Smoke runs must use configs with formal_benchmark: false."
        )

    run_id = cfg.get("run_id", config_path.stem)
    dataset_id = (cfg.get("dataset") or {}).get("id", "UNKNOWN")
    h5ad_rel = (cfg.get("dataset") or {}).get("h5ad_path", "")
    time_key = (cfg.get("dataset") or {}).get("time_key", "abs_day")

    gt_cfg = cfg.get("ground_truth") or {}
    provider_id = gt_cfg.get("provider_id", "")
    state_key = gt_cfg.get("state_key", "")
    confidence_mode = gt_cfg.get("confidence_mode", "all")
    exclude_uncertain = bool(gt_cfg.get("exclude_uncertain_states", False))

    lin_cfg = cfg.get("lineage") or {}
    graph_rel = lin_cfg.get("reference_graph_path", "")

    out_rel = (cfg.get("output") or {}).get("base_dir", "")

    # Resolve absolute paths.
    h5ad_path = root / h5ad_rel if h5ad_rel else None
    graph_path = root / graph_rel if graph_rel else None

    # Smoke output directory (always under smoke root, never in formal dirs).
    smoke_dir = output_root / run_id
    smoke_dir.mkdir(parents=True, exist_ok=True)

    print(f"[smoke] run_id        : {run_id}")
    print(f"[smoke] provider_id   : {provider_id}")
    print(f"[smoke] state_key     : {state_key}")
    print(f"[smoke] graph_path    : {graph_rel}")
    print(f"[smoke] h5ad_path     : {h5ad_rel}")
    print(f"[smoke] output_dir    : {smoke_dir}")

    if dry_run:
        print("[smoke] DRY RUN — no files written.")
        return {"dry_run": True, "run_id": run_id, "output_dir": str(smoke_dir)}

    # -- Validate inputs --
    if not h5ad_path or not h5ad_path.exists():
        raise FileNotFoundError(f"h5ad not found: {h5ad_path}")
    if not graph_path or not graph_path.exists():
        raise FileNotFoundError(f"Reference graph not found: {graph_path}")

    # -- Load reference graph --
    from benchmark.evaluation.eval_lineage import (
        load_reference_graph,
        run_lineage_evaluation,
    )

    ref_matrix, ref_edges, node_ids = load_reference_graph(
        str(graph_path),
        edge_confidence_mode=confidence_mode,
        exclude_uncertain_states=exclude_uncertain,
    )
    print(f"[smoke] Reference graph: {len(node_ids)} nodes, {len(ref_edges)} edges")

    # -- Labels to exclude: anything not a node in the reference graph --
    graph_node_set = set(node_ids)
    # "ambiguous" is never in the graph; collect any label outside the graph.
    # We will discover which ones after loading adata.

    # -- Build oracle STM from reference graph --
    import pandas as pd
    oracle_stm = _build_oracle_stm(node_ids, ref_edges)
    oracle_edges_df = _build_oracle_edges(oracle_stm)

    # Write STM and edges CSV to smoke dir.
    stm_path = smoke_dir / "state_transition_matrix.csv"
    edges_csv_path = smoke_dir / "lineage_graph_edges.csv"
    oracle_stm.to_csv(stm_path)
    oracle_edges_df.to_csv(edges_csv_path, index=False)
    print(f"[smoke] Oracle STM written  : {stm_path}")
    print(f"[smoke] Oracle edges written: {edges_csv_path}")

    # -- Load h5ad and filter ambiguous / out-of-graph cells --
    # Discover which labels are outside the reference graph.
    import h5py as _h5py
    with _h5py.File(h5ad_path, "r") as _f:
        _item = _f["obs"][state_key]
        if isinstance(_item, _h5py.Dataset):
            _raw = _item[:].astype(str)
        else:
            _cats = _item["categories"][:].astype(str)
            _codes = _item["codes"][:]
            _raw = _cats[_codes]
    all_labels = set(_raw.tolist())

    excluded_labels = sorted(all_labels - graph_node_set)
    print(f"[smoke] Labels in h5ad      : {sorted(all_labels)}")
    print(f"[smoke] Excluding labels    : {excluded_labels}")

    adata_filtered, excl_summary = _load_adata_filtered(h5ad_path, state_key, excluded_labels)
    print(
        f"[smoke] Cells: {excl_summary['n_cells_before_exclusion']} total, "
        f"{excl_summary['n_cells_excluded']} excluded, "
        f"{excl_summary['n_cells_after_exclusion']} retained"
    )

    # -- Load registry entry for metadata --
    try:
        import yaml as _yaml
        _reg_path = root / "benchmark" / "ground_truth" / "registry.yaml"
        _reg_raw = _yaml.safe_load(open(_reg_path, encoding="utf-8")) or {}
        reg_entry = (_reg_raw.get("providers") or {}).get(provider_id) or {}
    except Exception as _exc:
        print(f"[smoke] WARN: could not load registry: {_exc}")
        reg_entry = {}

    label_mode = (
        reg_entry.get("label_mode")
        or cfg.get("state_system", {}).get("label_mode", "")
        or ""
    )
    analysis_role = (
        reg_entry.get("analysis_role", "")
        or cfg.get("state_system", {}).get("analysis_role", "")
        or ""
    )

    ground_truth_meta = {
        "provider_id": provider_id,
        "state_key": state_key,
        "confidence_mode": confidence_mode,
        "exclude_uncertain_states": exclude_uncertain,
        "label_mode": label_mode,
        "analysis_role": analysis_role,
        "n_states": len(node_ids),
        "n_graph_edges": len(ref_edges),
    }

    # -- Run lineage evaluation --
    metrics = run_lineage_evaluation(
        state_transition_matrix_path=str(stm_path),
        lineage_graph_edges_path=str(edges_csv_path),
        output_dir=str(smoke_dir),
        reference_graph_path=str(graph_path),
        adata=adata_filtered,
        edge_confidence_mode=confidence_mode,
        exclude_uncertain_states=exclude_uncertain,
        cell_state_key=state_key,
        time_key=time_key,
        ground_truth=ground_truth_meta,
    )

    # -- Write smoke metadata --
    smoke_meta: Dict[str, Any] = {
        "config_path": str(config_path),
        "run_id": run_id,
        "dataset_id": dataset_id,
        "provider_id": provider_id,
        "label_mode": label_mode,
        "analysis_role": analysis_role,
        "cell_state_key": state_key,
        "reference_graph_path": graph_rel,
        "input_h5ad": h5ad_rel,
        "smoke_test": True,
        "formal_benchmark": False,
        "smoke_stm_type": "oracle_reference_graph",
        "smoke_stm_note": (
            "STM was constructed from the reference graph edges directly "
            "(NOT a method prediction). Used only to exercise the eval pipeline."
        ),
        "excluded_labels": excluded_labels,
        "exclusion_reason": excl_summary["exclusion_reason"],
        "n_cells_before_exclusion": excl_summary["n_cells_before_exclusion"],
        "n_cells_after_exclusion": excl_summary["n_cells_after_exclusion"],
        "n_reference_nodes": len(node_ids),
        "reference_node_ids": node_ids,
        "n_reference_edges": len(ref_edges),
        "output_dir": str(smoke_dir),
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "git_commit": _git_commit(root),
        "metrics_summary": {
            k: metrics.get(k)
            for k in ("auroc", "auprc", "jaccard_similarity", "jaccard_similarity_topk",
                       "single_step_recovery", "multi_step_recovery", "status")
        },
    }

    smoke_meta_path = smoke_dir / "smoke_metadata.json"
    with open(smoke_meta_path, "w", encoding="utf-8") as f:
        json.dump(smoke_meta, f, indent=2)
    print(f"[smoke] Smoke metadata  : {smoke_meta_path}")
    print(f"[smoke] AUROC           : {metrics.get('auroc')}")
    print(f"[smoke] AUPRC           : {metrics.get('auprc')}")
    print(f"[smoke] Jaccard (top-k) : {metrics.get('jaccard_similarity_topk')}")
    print(f"[smoke] Status          : {metrics.get('status')}")
    return metrics


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run milestone lineage fidelity smoke test for a config YAML."
    )
    parser.add_argument("config", help="Path to milestone method config YAML.")
    parser.add_argument(
        "--output-root",
        default="benchmark/results/smoke/milestone_lineage",
        help="Root directory for smoke outputs (default: benchmark/results/smoke/milestone_lineage).",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would be done without writing any files.",
    )
    args = parser.parse_args(argv)

    root = _project_root()
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = root / config_path
    output_root = Path(args.output_root)
    if not output_root.is_absolute():
        output_root = root / output_root

    try:
        run_smoke(config_path, output_root, dry_run=args.dry_run)
    except Exception as exc:
        print(f"[smoke] ERROR: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
