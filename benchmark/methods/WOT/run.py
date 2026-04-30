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
                       output_dir: Path, time_key: str = "time_label",
                       skip_if_exists: bool = False):
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
    skip_if_exists : bool
        If True and the tmap directory already contains at least one tmap_*.h5ad
        file, skip the OT computation entirely and reload the existing files.
        Useful during development to avoid re-running the expensive OT step when
        only the aggregation or evaluation code has changed.
        Default: False (always recompute; safe / reproducible default).

    Returns
    -------
    tmap_model : wot.tmap.TransportMapModel
    """
    import wot

    tmap_dir = output_dir / "tmaps"
    tmap_prefix = str(tmap_dir / "tmap")

    # ------------------------------------------------------------------
    # Caching: reuse existing transport maps if requested and available.
    # ------------------------------------------------------------------
    if skip_if_exists:
        tmap_dir.mkdir(parents=True, exist_ok=True)
        existing = sorted(tmap_dir.glob("tmap_*.h5ad"))
        if existing:
            print(
                f"  [cache] --skip-tmap-if-exists: found {len(existing)} existing tmap(s) "
                f"in {tmap_dir} — skipping OT computation."
            )
            print(f"  Loading cached transport maps from {tmap_dir} ...")
            return wot.tmap.TransportMapModel.from_directory(tmap_prefix)
        else:
            print(
                f"  [cache] --skip-tmap-if-exists set but no existing tmap_*.h5ad files "
                f"found in {tmap_dir} — running OT computation."
            )

    # ------------------------------------------------------------------
    # Inject growth rates into obs in-place (no full copy needed).
    # OTModel reads growth rates from obs[growth_rate_field]; we add the
    # column directly to the caller's AnnData.  The column is only used by
    # OTModel and does not affect downstream aggregation or evaluation.
    # ------------------------------------------------------------------
    train_adata.obs["cell_growth_rate"] = growth_rates

    tmap_dir.mkdir(parents=True, exist_ok=True)
    # compute_all_transport_maps writes files as {prefix}_{t0}_{t1}.h5ad

    print(f"  Initializing wot.ot.OTModel (day_field={time_key!r}) ...")
    ot_model = wot.ot.OTModel(
        train_adata,
        day_field=time_key,               # use configured time_key, not hardcoded "time_label"
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

    Implementation — vectorized sparse aggregation (replaces Python double loop):
    ---
    For each time-pair, the transport matrix M has shape (n_src_cells × n_tgt_cells).
    We want:
        stm[state_s, state_t] = Σ_{i ∈ state_s, j ∈ state_t}  M[i, j]

    This equals  (S_src @ M @ S_tgt.T)[state_s, state_t]  where:
      S_src  (n_states × n_src_cells): indicator matrix, 1 where cell i belongs to state s
      S_tgt  (n_states × n_tgt_cells): indicator matrix, 1 where cell j belongs to state t

    The two sparse matrix multiplications run in milliseconds for 10 K-cell transport
    maps — versus tens of millions of Python iterations in the original loop.
    This is the same mathematical pattern used by CellRank2Adapter._aggregate_cr2_to_state_level().

    Returns
    -------
    state_transition_matrix : pd.DataFrame
        Rows = source states, columns = target states, values = row-normalised weights.
    lineage_graph_edges : pd.DataFrame
        Columns: source_state, target_state, weight.  Only nonzero entries.
    """
    import scipy.sparse as sp

    obs = train_adata.obs.copy()
    states = sorted(obs[cell_state_key].unique())
    n_states = len(states)
    # Integer index for each state label — built once, reused across all time pairs.
    state_to_idx = {s: i for i, s in enumerate(states)}
    unique_times = sorted(obs[time_key].unique())

    # Accumulate raw (un-normalised) state-level mass across all consecutive time pairs.
    # A plain numpy array is faster to write into than a pandas DataFrame with .loc.
    stm_raw = np.zeros((n_states, n_states), dtype=np.float64)

    t_agg_start = time.time()
    pairs_processed = 0

    for t_src, t_tgt in zip(unique_times[:-1], unique_times[1:]):
        try:
            # get_coupling returns AnnData: obs = src cells, var = tgt cells, X = transport matrix.
            tmap = tmap_model.get_coupling(t_src, t_tgt)
        except Exception as e:
            print(f"  Warning: could not get coupling {t_src}→{t_tgt}: {e}")
            continue

        # Keep the transport matrix sparse for efficient multiplication.
        # WOT may return a dense ndarray or a scipy sparse matrix depending on version.
        M = tmap.X
        if not sp.issparse(M):
            M = sp.csr_matrix(np.asarray(M))
        else:
            M = M.tocsr()

        # Row/column ordering is authoritative from the coupling AnnData, not from
        # the time-point mask on obs.  Use tmap.obs_names / tmap.var_names directly.
        src_states_arr = obs.loc[tmap.obs_names, cell_state_key].values
        tgt_states_arr = obs.loc[tmap.var_names, cell_state_key].values
        n_src = len(src_states_arr)
        n_tgt = len(tgt_states_arr)

        # Integer state index for every cell in source / target sets.
        src_idx = np.array([state_to_idx[s] for s in src_states_arr], dtype=np.int32)
        tgt_idx = np.array([state_to_idx[s] for s in tgt_states_arr], dtype=np.int32)

        # Build indicator matrices:
        #   S_src[s, i] = 1  iff source cell i belongs to state s
        #   S_tgt[t, j] = 1  iff target cell j belongs to state t
        S_src = sp.csr_matrix(
            (np.ones(n_src, dtype=np.float64), (src_idx, np.arange(n_src))),
            shape=(n_states, n_src),
        )
        S_tgt = sp.csr_matrix(
            (np.ones(n_tgt, dtype=np.float64), (tgt_idx, np.arange(n_tgt))),
            shape=(n_states, n_tgt),
        )

        # Vectorized aggregation:
        #   block[s, t] = Σ_{i∈s, j∈t} M[i,j]  ==  (S_src @ M @ S_tgt.T)[s, t]
        # .toarray() is safe here: block is (n_states × n_states), always small.
        block = (S_src @ M @ S_tgt.T).toarray()
        stm_raw += block
        pairs_processed += 1

    t_agg_elapsed = time.time() - t_agg_start
    print(f"  [timing] State aggregation — {pairs_processed} pair(s) in {t_agg_elapsed:.3f}s")

    # Wrap accumulated raw counts in a DataFrame, then row-normalise → stochastic matrix.
    stm = pd.DataFrame(stm_raw, index=states, columns=states)
    row_sums = stm.sum(axis=1)
    stm_norm = stm.div(row_sums.replace(0, np.nan), axis=0).fillna(0.0)

    # Build edge list from all nonzero entries of the normalised matrix.
    edges = []
    for src in stm_norm.index:
        for tgt in stm_norm.columns:
            w = float(stm_norm.loc[src, tgt])
            if w > 0.0:
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
    parser.add_argument(
        "--skip-tmap-if-exists",
        action="store_true",
        default=False,
        help=(
            "Reuse existing tmap_*.h5ad files in the output tmaps/ directory instead of "
            "rerunning the WOT OT computation.  Safe to use during development when only "
            "the aggregation or evaluation code has changed.  Default: False (always recompute)."
        ),
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help=(
            "Override the output base directory from the config file. "
            "Accepts an absolute path or a path relative to the project root. "
            "Useful for HPC job arrays that set output paths externally."
        ),
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

    pkl_path = project_root / pkl_path_raw

    # --output-dir CLI flag overrides the config's output.base_dir when provided.
    if args.output_dir:
        cli_out = Path(args.output_dir)
        output_dir = cli_out if cli_out.is_absolute() else project_root / cli_out
    else:
        output_dir = project_root / output_base

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
        # Fallback: load h5ad directly without pkl.
        # Memory-efficient path: when --subsample or scenario_params.train_times
        # is set, use backed='r' so the full gene matrix is NOT loaded into RAM
        # upfront. We materialize only the needed cell subset after filtering.
        # For full-data runs (no subsample, no train_times) backed mode is also
        # fine — the copy() call in scenario-time filtering will materialize.
        print(f"  pkl not found at {pkl_path}. Falling back to h5ad direct load.")
        import anndata as ad
        h5ad_full = project_root / h5ad_fallback

        # Detect whether we'll need only a row-subset.
        _scenario_cfg_pre = cfg.get("scenario_params", cfg.get("scenario_split", {})) or {}
        _train_times_pre  = _scenario_cfg_pre.get("train_times")
        _use_backed = (args.subsample is not None) or bool(_train_times_pre)

        if _use_backed:
            # Backed read: obs is loaded; X stays on disk until explicitly accessed.
            adata = ad.read_h5ad(h5ad_full, backed="r")
            print(f"  Backed read: {adata.n_obs} cells × {adata.n_vars} genes "
                  f"(will materialize filtered subset)")
        else:
            # Full in-memory load for runs that need all cells simultaneously.
            adata = ad.read_h5ad(h5ad_full)
            print(f"  Full load: {adata.n_obs} cells × {adata.n_vars} genes")

        train_adata = adata
        test_adata  = None  # not used by WOT Lineage Fidelity path

    # Optional subsample for fast local debug
    if args.subsample is not None:
        print(f"\n  [DEBUG] Subsampling to {args.subsample} cells per time point ...")
        import anndata as ad
        groups = []
        for t, grp in train_adata.obs.groupby(time_key, observed=True):
            idx = grp.index[:args.subsample]
            groups.append(train_adata[idx].to_memory()
                          if hasattr(train_adata[idx], "to_memory") else train_adata[idx])
        train_adata = ad.concat(groups)
        print(f"  Subsampled: {train_adata.n_obs} cells")

    # --- Scenario-level time-point filter (Scenario B / C) ---------------
    # If the config provides scenario.train_times, keep only cells whose
    # time value is in that list. This is how observed-time extrapolation
    # (Scenario B) is implemented without touching the dataset object —
    # the runner sees a filtered train_adata and WOT's OTModel only computes
    # couplings across the retained time points.
    #
    # Expected YAML shape:
    #   scenario:
    #     train_times: [0.5, 2.0, 4.0, 8.0, 12.0, 16.0]
    #     heldout_times: [...]   # documentation only; not used by WOT
    scenario_cfg = cfg.get("scenario_params", cfg.get("scenario_split", {})) or {}
    train_times = scenario_cfg.get("train_times")
    if train_times:
        import anndata as ad  # noqa: F401
        train_times = [float(t) for t in train_times]
        before = train_adata.n_obs
        mask = train_adata.obs[time_key].astype(float).isin(train_times).values
        # .copy() both applies the boolean mask AND materializes the slice from
        # backed storage (if the h5ad was opened with backed='r').
        train_adata = train_adata[mask].to_memory()
        print(
            f"\n  [scenario] Restricted train_adata to {train_times}: "
            f"{before} → {train_adata.n_obs} cells "
            f"(dropped {before - train_adata.n_obs})"
        )
        remaining_times = sorted(train_adata.obs[time_key].astype(float).unique())
        print(f"  [scenario] Remaining time points: {remaining_times}")

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
            t_ot_start = time.time()
            tmap_model = _run_wot_transport(
                train_adata, cell_days, growth_rates,
                wot_params=wot_cfg, output_dir=output_dir,
                time_key=time_key,            # propagate configured time_key
                skip_if_exists=args.skip_tmap_if_exists,
            )
            print(f"  [timing] Transport map step: {time.time() - t_ot_start:.1f}s")

            print(f"\n[4/5] Aggregating to cell-state level ...")
            t_agg_outer = time.time()
            stm, edges = _aggregate_to_state_level(
                train_adata, tmap_model,
                cell_state_key=cell_state_key,
                time_key=time_key,
            )
            print(f"  [timing] Aggregation step total: {time.time() - t_agg_outer:.3f}s")
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
            # Pass the loaded train_adata plus configured state/time keys so
            # the scTimeBench-style correlation baseline uses the exact same
            # cell universe and adjacent training timepoints as the method.
            metrics = run_lineage_evaluation(
                state_transition_matrix_path=str(output_dir / "state_transition_matrix.csv"),
                lineage_graph_edges_path=str(output_dir / "lineage_graph_edges.csv"),
                output_dir=str(output_dir),
                reference_graph_path=reference_graph_path,
                adata=train_adata,
                edge_confidence_mode=edge_confidence_mode,
                exclude_uncertain_states=exclude_uncertain_states,
                cell_state_key=cell_state_key,
                time_key=time_key,
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
