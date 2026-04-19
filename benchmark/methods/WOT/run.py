"""
benchmark/methods/WOT/run.py
WOT smoke-test runner for the scTimeBench-aligned benchmark.
Framework reference: experimental framework v2.md §10.1, §14

Usage:
    conda activate traj_env
    cd C:\\Users\\37620\\trajectory
    python benchmark/methods/WOT/run.py \
        --config benchmark/configs/wot_gse230659_observed.yaml

What this script does:
    1. Reads the YAML config.
    2. Loads the dataset from the pkl (or falls back to the h5ad directly).
    3. Calls dataset.load_data() to obtain (train_adata, test_adata).
    4. Prepares WOT inputs from train_adata:
         - cell_days array  (from obs[dataset.time_key] — "time_label" for the
                             base benchmark adata, "abs_day" for scGPT adatas)
         - growth rates     (uniform prior for smoke test)
    5. Computes consecutive transport maps via wot.ot.OTModel.
       day_field is set to the configured dataset.time_key (not hardcoded).
    6. Aggregates cell-level transport to the cell-state level using
       obs[lineage.cell_state_key] (from config; e.g. "scTimeBench_cell_type"
       for the base config, "scgpt_pseudostate_provisional" for scGPT-v1).
    7. Writes benchmark output files to the configured output directory.
    8. Calls eval_lineage.py to compute Lineage Fidelity metrics
       (active when lineage.reference_graph_path is set in config).

Capability:
    WOT supports_unseen_timepoint_projection = False
    WOT supports_lineage_inference           = True
    → only Lineage Fidelity outputs are produced.
    → Forecast Accuracy and Embedding Coherence are explicitly skipped.
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Project-root / path helpers
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


def _setup_sys_path(project_root: Path):
    for p in [str(project_root), str(project_root / "benchmark")]:
        if p not in sys.path:
            sys.path.insert(0, p)


# ---------------------------------------------------------------------------
# YAML loader (PyYAML if available, else fallback to stdlib)
# ---------------------------------------------------------------------------

def _load_yaml(path: Path) -> dict:
    try:
        import yaml
        with open(path, "r", encoding="utf-8-sig") as f:
            return yaml.safe_load(f)
    except ImportError:
        # Minimal fallback: parse simple key: value lines only.
        # Enough for reading scalar fields needed to locate pkl and output_dir.
        cfg: dict = {}
        with open(path, "r", encoding="utf-8-sig") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#") and ":" in line:
                    k, _, v = line.partition(":")
                    cfg[k.strip()] = v.strip().strip('"').strip("'")
        return cfg


# ---------------------------------------------------------------------------
# WOT availability check
# ---------------------------------------------------------------------------

def _check_wot() -> bool:
    try:
        import wot                   # noqa: F401
        return True
    except ImportError:
        return False


# ---------------------------------------------------------------------------
# Core steps
# ---------------------------------------------------------------------------

def _prepare_wot_inputs(train_adata, time_key: str, growth_rate_source: str):
    """
    Extract WOT-required arrays from train_adata.

    Returns
    -------
    cell_days : np.ndarray, shape (n_cells,)
        Numeric time label for each cell.
    growth_rates : np.ndarray, shape (n_cells,)
        Per-cell growth rate (1.0 = neutral prior for smoke test).
    """
    cell_days = train_adata.obs[time_key].values.astype(np.float64)

    if growth_rate_source == "uniform":
        growth_rates = np.ones(train_adata.n_obs, dtype=np.float64)
        print("  Growth rates: uniform prior (all 1.0)")
    else:
        raise NotImplementedError(
            f"growth_rate_source='{growth_rate_source}' not yet implemented. "
            "Use 'uniform' for the smoke test."
        )

    return cell_days, growth_rates


def _run_wot_transport(train_adata, cell_days, growth_rates, wot_params: dict,
                       output_dir: Path, time_key: str = "time_label"):
    """
    Run WOT transport map computation and save transport maps.

    API used: wot.ot.OTModel (WOT 1.0.x).
      - OTModel takes an AnnData directly; reads time from obs[day_field]
        and growth rates from obs[growth_rate_field].
      - compute_all_transport_maps() writes one h5ad per consecutive pair.
      - TransportMapModel.from_directory() reloads them for downstream use.

    Parameters
    ----------
    time_key : str
        The obs column holding numeric time labels.  Must match the column
        that is actually present in train_adata.obs.  Passed directly as
        day_field to wot.ot.OTModel.  Defaults to "time_label" for backward
        compatibility with the original benchmark adata; use "abs_day" for
        adata_scgpt_annotated.h5ad (scGPT-v1 configs).

    Returns
    -------
    tmap_model : wot.tmap.TransportMapModel
    """
    import wot

    # OTModel reads growth rates from obs[growth_rate_field]; inject before passing.
    train_adata = train_adata.copy()
    train_adata.obs["cell_growth_rate"] = growth_rates

    tmap_dir = output_dir / "tmaps"
    tmap_dir.mkdir(parents=True, exist_ok=True)
    # compute_all_transport_maps writes files as {prefix}_{t0}_{t1}.h5ad
    tmap_prefix = str(tmap_dir / "tmap")

    print(f"  Initializing wot.ot.OTModel (day_field={time_key!r}) ...")
    ot_model = wot.ot.OTModel(
        train_adata,
        day_field=time_key,               # FIX: use configured time_key, not hardcoded "time_label"
        growth_rate_field="cell_growth_rate",
        epsilon=wot_params.get("epsilon", 0.05),
        lambda1=wot_params.get("lambda1", 1.0),
        lambda2=wot_params.get("lambda2", 50.0),
        local_pca=wot_params.get("local_pca", 30),
        growth_iters=wot_params.get("growth_iters", 3),
    )
    print(f"  Timepoints detected: {ot_model.timepoints}")

    print(f"  Computing all transport maps → {tmap_dir} ...")
    ot_model.compute_all_transport_maps(tmap_out=tmap_prefix, overwrite=True)

    print(f"  Loading transport maps from {tmap_dir} ...")
    tmap_model = wot.tmap.TransportMapModel.from_directory(tmap_prefix)

    return tmap_model


def _aggregate_to_state_level(train_adata, tmap_model, cell_state_key: str,
                               time_key: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Aggregate cell-level WOT transport probabilities to the state level.

    Steps (v2 §9.4):
      1. For each consecutive time-pair (t → t+1), retrieve the transport map.
      2. Sum probabilities within each (source_state, target_state) pair.
      3. Row-normalize to get a stochastic state-transition matrix.
      4. Build lineage graph edges from the matrix.

    Returns
    -------
    state_transition_matrix : pd.DataFrame
        Rows = source states, columns = target states, values = transition weights.
    lineage_graph_edges : pd.DataFrame
        Columns: source_state, target_state, weight.
    """
    import scipy.sparse as sp

    obs = train_adata.obs.copy()
    states = sorted(obs[cell_state_key].unique())
    unique_times = sorted(obs[time_key].unique())

    # Accumulate state-level transitions across all consecutive time pairs.
    stm = pd.DataFrame(0.0, index=states, columns=states)

    for t_src, t_tgt in zip(unique_times[:-1], unique_times[1:]):
        src_mask = obs[time_key] == t_src
        tgt_mask = obs[time_key] == t_tgt
        src_ids  = obs.index[src_mask]
        tgt_ids  = obs.index[tgt_mask]

        try:
            tmap = tmap_model.get_coupling(t_src, t_tgt)
        except Exception as e:
            print(f"  Warning: could not get coupling {t_src}→{t_tgt}: {e}")
            continue

        # get_coupling returns AnnData: obs=src cells, var=tgt cells, X=transport matrix
        M = tmap.X
        if sp.issparse(M):
            M = M.toarray()

        # Use cell IDs from the coupling AnnData (authoritative row/col order)
        src_ids = tmap.obs_names
        tgt_ids = tmap.var_names

        # Aggregate: for each (source_state, target_state) pair, sum transport mass
        src_states = obs.loc[src_ids, cell_state_key].values
        tgt_states = obs.loc[tgt_ids, cell_state_key].values

        for si, ss in enumerate(src_states):
            row_sum = M[si, :].sum()
            if row_sum == 0:
                continue
            for tj, ts in enumerate(tgt_states):
                stm.loc[ss, ts] += M[si, tj]

    # Row-normalize
    row_sums = stm.sum(axis=1)
    stm_norm = stm.div(row_sums.replace(0, np.nan), axis=0).fillna(0.0)

    # Build edges
    edges = []
    for src in stm_norm.index:
        for tgt in stm_norm.columns:
            w = stm_norm.loc[src, tgt]
            if w > 0:
                edges.append({"source_state": src, "target_state": tgt, "weight": w})
    edges_df = pd.DataFrame(edges, columns=["source_state", "target_state", "weight"])

    return stm_norm, edges_df


def _write_outputs(stm: pd.DataFrame, edges: pd.DataFrame,
                   output_dir: Path, run_meta: dict):
    """Write the three Lineage Fidelity output files + run_metadata.json."""
    output_dir.mkdir(parents=True, exist_ok=True)

    stm_path   = output_dir / "state_transition_matrix.csv"
    edges_path = output_dir / "lineage_graph_edges.csv"
    meta_path  = output_dir / "run_metadata.json"

    stm.to_csv(stm_path)
    edges.to_csv(edges_path, index=False)

    with open(meta_path, "w") as f:
        json.dump(run_meta, f, indent=2)

    print(f"  state_transition_matrix → {stm_path}")
    print(f"  lineage_graph_edges     → {edges_path}")
    print(f"  run_metadata            → {meta_path}")


def _write_scaffold_outputs(output_dir: Path, run_meta: dict,
                             reason: str = "wot_not_installed"):
    """
    Write placeholder output files when WOT is not installed.
    This keeps the pipeline end-to-end runnable as a scaffold.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    stm = pd.DataFrame(
        index=pd.Index([], name="source_state"),
        columns=pd.Index([], name="target_state"),
    )
    edges = pd.DataFrame(columns=["source_state", "target_state", "weight"])

    run_meta["status"] = f"scaffold_only ({reason})"
    run_meta["dimensions_executed"] = []

    stm.to_csv(output_dir / "state_transition_matrix.csv")
    edges.to_csv(output_dir / "lineage_graph_edges.csv", index=False)
    with open(output_dir / "run_metadata.json", "w") as f:
        json.dump(run_meta, f, indent=2)

    print(f"  Scaffold outputs written to {output_dir}")
    print(f"  Reason: {reason}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="WOT smoke-test runner (scTimeBench v2)"
    )
    parser.add_argument(
        "--config",
        default="benchmark/configs/wot_gse230659_observed.yaml",
        help="Path to YAML config (relative to project root or absolute)",
    )
    parser.add_argument(
        "--subsample",
        type=int,
        default=None,
        help="Subsample to N cells per time point for fast debugging (default: use all)",
    )
    args = parser.parse_args()

    t0             = time.time()
    project_root   = _find_project_root()
    _setup_sys_path(project_root)

    print("=" * 65)
    print("WOT benchmark runner  (scTimeBench v2 — Lineage Fidelity only)")
    print("=" * 65)
    print(f"  Project root : {project_root}")

    # --- Load config ---
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = project_root / config_path
    cfg = _load_yaml(config_path)
    print(f"  Config       : {config_path}")

    # Resolve nested keys robustly
    dataset_cfg  = cfg.get("dataset", cfg)
    output_cfg   = cfg.get("output", cfg)
    wot_cfg      = cfg.get("wot_params", {})
    lineage_cfg  = cfg.get("lineage", {})

    pkl_path_raw  = dataset_cfg.get("pkl_path", "benchmark/datasets/gse230659_observed.dataset.pkl")
    h5ad_fallback = dataset_cfg.get("h5ad_path", "data/processed/adata_benchmark.h5ad")
    output_base   = output_cfg.get("base_dir", "benchmark/results/wot/scenario_A")
    time_key      = dataset_cfg.get("time_key", "time_label")
    cell_state_key          = lineage_cfg.get("cell_state_key", "scTimeBench_cell_type")
    edge_confidence_mode    = lineage_cfg.get("edge_confidence_mode", "all")
    exclude_uncertain_states = lineage_cfg.get("exclude_uncertain_states", False)
    growth_source           = wot_cfg.get("growth_rate_source", "uniform")

    pkl_path    = project_root / pkl_path_raw
    output_dir  = project_root / output_base

    # Build minimal run_metadata now; update at end with final status.
    run_meta = {
        "method":    "wot",
        "dataset":   "GSE230659",
        "scenario":  cfg.get("scenario", "A"),
        "capability_flags": {
            "supports_unseen_timepoint_projection": False,
            "supports_lineage_inference": True,
        },
        "dimensions_executed": ["lineage_fidelity"],
        "config":    str(config_path),
        "runtime_seconds": None,
        "status":    "running",
        "notes": (
            "WOT evaluated on Lineage Fidelity only. "
            "Forecast Accuracy and Embedding Coherence explicitly skipped."
        ),
    }

    # --- Load dataset ---
    print(f"\n[1/5] Loading dataset ...")
    if pkl_path.exists():
        print(f"  Loading from pkl: {pkl_path}")
        with open(pkl_path, "rb") as f:
            dataset = pickle.load(f)
        # Re-resolve h5ad path relative to this machine's project root
        # (the pkl was built on the same machine, but be defensive)
        if not dataset.h5ad_path.exists():
            dataset.h5ad_path = project_root / "data" / "processed" / "adata_benchmark.h5ad"
        train_adata, test_adata = dataset.load_data()
    else:
        # Fallback: load h5ad directly without pkl
        print(f"  pkl not found at {pkl_path}. Falling back to h5ad direct load.")
        import anndata as ad
        adata = ad.read_h5ad(project_root / h5ad_fallback)
        train_adata = adata
        test_adata  = adata[adata.obs[time_key] == 30.0].copy()
        print(f"  Loaded: {train_adata.n_obs} cells × {train_adata.n_vars} genes")

    # Optional subsample for fast local debug
    if args.subsample is not None:
        print(f"\n  [DEBUG] Subsampling to {args.subsample} cells per time point ...")
        import anndata as ad
        groups = []
        for t, grp in train_adata.obs.groupby(time_key, observed=True):
            idx = grp.index[:args.subsample]
            groups.append(train_adata[idx])
        train_adata = ad.concat(groups)
        print(f"  Subsampled: {train_adata.n_obs} cells")

    print(f"\n[2/5] Preparing WOT inputs ...")
    cell_days, growth_rates = _prepare_wot_inputs(
        train_adata, time_key=time_key, growth_rate_source=growth_source
    )
    print(f"  Time points: {sorted(np.unique(cell_days))}")

    # --- WOT execution ---
    wot_available = _check_wot()

    if wot_available:
        print(f"\n[3/5] Running WOT transport maps ...")
        try:
            tmap_model = _run_wot_transport(
                train_adata, cell_days, growth_rates,
                wot_params=wot_cfg, output_dir=output_dir,
                time_key=time_key,            # propagate configured time_key
            )

            print(f"\n[4/5] Aggregating to cell-state level ...")
            stm, edges = _aggregate_to_state_level(
                train_adata, tmap_model,
                cell_state_key=cell_state_key,
                time_key=time_key,
            )
            print(f"  State-transition matrix shape: {stm.shape}")
            print(f"  Lineage graph edges: {len(edges)} edges")

            print(f"\n[5/5] Writing outputs → {output_dir}")
            run_meta["status"]      = "completed"
            run_meta["cell_state_key"] = cell_state_key
            _write_outputs(stm, edges, output_dir, run_meta)

        except Exception as exc:
            print(f"\n  WOT run failed: {exc}")
            print("  Writing scaffold outputs and continuing ...")
            run_meta["error"] = str(exc)
            _write_scaffold_outputs(output_dir, run_meta, reason=f"wot_error: {exc}")

    else:
        print(f"\n[3/5] WOT not installed — writing scaffold outputs ...")
        print("  Install with: pip install wot")
        _write_scaffold_outputs(
            output_dir, run_meta, reason="wot_not_installed"
        )

    # --- Eval ---
    print(f"\n[eval] Running eval_lineage.py ...")
    try:
        # Add project_root to path so evaluation imports work
        eval_path = project_root / "benchmark" / "evaluation" / "eval_lineage.py"
        if eval_path.exists():
            sys.path.insert(0, str(project_root / "benchmark" / "evaluation"))
            from eval_lineage import run_lineage_evaluation   # noqa: E402
            reference_graph_path = lineage_cfg.get("reference_graph_path", None)
            if reference_graph_path == "null" or reference_graph_path == "":
                reference_graph_path = None
            metrics = run_lineage_evaluation(
                state_transition_matrix_path=str(output_dir / "state_transition_matrix.csv"),
                lineage_graph_edges_path=str(output_dir / "lineage_graph_edges.csv"),
                output_dir=str(output_dir),
                reference_graph_path=reference_graph_path,
                adata=None,
                edge_confidence_mode=edge_confidence_mode,
                exclude_uncertain_states=exclude_uncertain_states,
            )
            print(f"  Lineage Fidelity status: {metrics.get('status', 'unknown')}")
        else:
            print(f"  eval_lineage.py not found at {eval_path} — skipping eval.")
    except Exception as exc:
        print(f"  eval_lineage.py failed: {exc} — skipping.")

    # --- Final report ---
    elapsed = time.time() - t0
    run_meta["runtime_seconds"] = round(elapsed, 1)

    # Update run_metadata with final runtime
    meta_path = output_dir / "run_metadata.json"
    if meta_path.exists():
        with open(meta_path) as f:
            meta_existing = json.load(f)
        meta_existing["runtime_seconds"] = run_meta["runtime_seconds"]
        with open(meta_path, "w") as f:
            json.dump(meta_existing, f, indent=2)

    print("\n" + "=" * 65)
    print(f"WOT run complete")
    print(f"  Output dir  : {output_dir}")
    print(f"  Runtime     : {elapsed:.1f}s")
    print(f"  Status      : {run_meta['status']}")
    print("=" * 65)


if __name__ == "__main__":
    main()
