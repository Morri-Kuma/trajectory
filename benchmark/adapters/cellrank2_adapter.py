"""
cellrank2_adapter.py
CellRank2 method adapter for the scTimeBench-aligned benchmark.
Framework reference: docs/framework/experimental_framework_v2.md §10.2

CellRank2 capability flags (per v2 §6):
  supports_unseen_timepoint_projection = False
  supports_lineage_inference           = True

CellRank2 is evaluated on Lineage Fidelity only.
It is explicitly excluded from Forecast Accuracy and Embedding Coherence.
No OT-projection workaround or synthetic future-cell layer is used.

Kernel choice (per config §cellrank2_params):
  RealTimeKernel.from_wot() — builds a CellRank2 RealTimeKernel from WOT
  optimal transport maps computed internally by this adapter.  This is the
  recommended CellRank2 kernel for observed-time data with measured time points
  and is directly comparable to the WOT adapter's transport-map approach.

Aggregation:
  The cell-level transition matrix T (n_cells × n_cells) produced by the kernel
  is aggregated to the state level via:
      stm_raw = S @ T @ S.T
  where S is the binary cell-state indicator matrix (n_states × n_cells).
  The result is row-normalized to form a stochastic state-transition matrix.
  This aggregation is identical in logic to the WOT adapter's
  _aggregate_to_state_level() but operates on the full T rather than
  per-timepoint couplings.

Configuration (injected into scenario_config by eval_dispatch.py):
  time_key         : obs column for experimental time (default "abs_day")
  cell_state_key   : obs column for cell-state labels
                     (default "final_milestone_label_coarse")
  cellrank2_params : dict — kernel and WOT sub-parameters from the method config
                     (optional; defaults are used if absent)
"""

import traceback

import numpy as np
import pandas as pd
import scipy.sparse as sp
from pathlib import Path

from .base_adapter import BaseAdapter
from benchmark.shared.dataset.preprocessors.scenario_timepoint_split import (
    split_adata_by_timepoints,
)


# ------------------------------------------------------------------
# Module-level aggregation helper
# ------------------------------------------------------------------

def _aggregate_cr2_to_state_level(
    adata,
    transition_matrix,
    cell_state_key: str,
):
    """
    Aggregate a CellRank2 cell-level transition matrix to the state level.

    Parameters
    ----------
    adata : AnnData
        Must have adata.obs[cell_state_key] populated for all cells.
    transition_matrix : scipy sparse or dense array, shape (n_cells, n_cells)
        Cell-level transition probabilities from a CellRank2 kernel.
    cell_state_key : str
        obs column holding the cell-state label for each cell.

    Returns
    -------
    stm_norm : pd.DataFrame
        Row-normalized state-transition matrix.
        Index = source states, columns = target states.
    edges_df : pd.DataFrame
        Columns: source_state, target_state, weight.
        One row per nonzero (source, target) entry in stm_norm.
    """
    cell_states = adata.obs[cell_state_key].values
    states = sorted(set(cell_states))
    n_states = len(states)
    n_cells = len(cell_states)

    state_to_idx = {s: i for i, s in enumerate(states)}
    cell_state_idx = np.array([state_to_idx[s] for s in cell_states], dtype=int)

    # Build cell-state indicator matrix S (n_states × n_cells).
    # S[s, i] = 1  iff  cell i belongs to state s.
    S = sp.csr_matrix(
        (
            np.ones(n_cells, dtype=float),
            (cell_state_idx, np.arange(n_cells)),
        ),
        shape=(n_states, n_cells),
    )

    # Ensure T is sparse for efficient multiplication.
    T = transition_matrix
    if not sp.issparse(T):
        T = sp.csr_matrix(np.asarray(T))

    # Aggregate: stm_raw[s, t] = Σ_{i∈s, j∈t} T[i, j]  =  (S @ T @ S.T)[s, t]
    stm_raw = (S @ T @ S.T).toarray()
    stm_df = pd.DataFrame(stm_raw, index=states, columns=states, dtype=float)

    # Row-normalize to a stochastic matrix.
    row_sums = stm_df.sum(axis=1)
    stm_norm = stm_df.div(row_sums.replace(0, np.nan), axis=0).fillna(0.0)

    # Build edge list from nonzero entries.
    edge_rows = []
    for src in stm_norm.index:
        for tgt in stm_norm.columns:
            w = float(stm_norm.loc[src, tgt])
            if w > 0.0:
                edge_rows.append({"source_state": src, "target_state": tgt, "weight": w})
    edges_df = pd.DataFrame(
        edge_rows if edge_rows else [],
        columns=["source_state", "target_state", "weight"],
    )

    return stm_norm, edges_df


# ------------------------------------------------------------------
# Adapter class
# ------------------------------------------------------------------

class CellRank2Adapter(BaseAdapter):
    """
    Adapter for CellRank2.

    Eligible benchmark dimensions:
      - Lineage Fidelity only

    Not eligible for:
      - Forecast Accuracy
      - Embedding Coherence
    """

    # --- Capability flags ---
    supports_unseen_timepoint_projection: bool = False
    supports_lineage_inference: bool = True

    @property
    def method_id(self) -> str:
        return "cellrank2"

    # ------------------------------------------------------------------
    # Lineage Fidelity implementation
    # ------------------------------------------------------------------

    def _run_lineage_fidelity_impl(self, scenario_id: str) -> dict:
        """
        Run CellRank2 RealTimeKernel inference and construct the predicted lineage.

        Steps (per v2 §9.4):
          0. If scenario_params.train_times is set in the method config, filter
             self.adata to those time points. Mirrors benchmark/methods/WOT/run.py
             so Scenario B (observed-time extrapolation: early-only training) can
             be honored by this adapter. Filtering is done here, once, on the
             entry AnnData — both the WOT sub-fit and the CellRank2 kernel see
             the same filtered cell universe.
          1. Build WOT transport maps using wot.ot.OTModel (same params as WOT adapter).
          2. Initialize cr.kernels.RealTimeKernel.from_wot() from the transport maps.
          3. Compute the cell-level transition matrix.
          4. Aggregate to state level via S @ T @ S.T; row-normalize.
          5. Write state_transition_matrix.csv and lineage_graph_edges.csv.

        Falls back to empty scaffold outputs if CellRank2 or WOT is not installed,
        or if a runtime error occurs. The scaffold path writes empty files so that
        eval_lineage.py reports "completed_with_empty_predictions" in
        lineage_metrics.json — not "completed" — making the failure visible.
        """
        self.start_timer()

        adata = self.adata
        out_dir = self.output_dir

        # --- Read configuration (injected by eval_dispatch.py) ---
        time_key = self.scenario_config.get("time_key", "abs_day")
        cell_state_key = self.scenario_config.get(
            "cell_state_key", "final_milestone_label_coarse"
        )
        cr2_params = self.scenario_config.get("cellrank2_params", {})
        wot_params = cr2_params.get("wot_params", {})

        # -----------------------------------------------------------
        # Step 0: Scenario-level time-point filter (Scenario B / C)
        # -----------------------------------------------------------
        # scenario_params.train_times is injected by eval_dispatch.py when the
        # method config declares it (e.g. Scenario B = [0.5, 2, 4, 8, 12, 16]).
        # If present, restrict the adata used throughout this run to cells at
        # those time points. Equivalent to the block at the top of
        # benchmark/methods/WOT/run.py. Keeps the smallest possible footprint —
        # the adapter still operates on a single AnnData object; only the row
        # set is reduced.
        scenario_cfg = (
            self.scenario_config.get("scenario_params")
            or self.scenario_config.get("scenario_split")
            or {}
        )
        train_times = scenario_cfg.get("train_times") if scenario_cfg else None
        if train_times:
            train_times_f = [float(t) for t in train_times]
            if time_key not in adata.obs.columns:
                raise RuntimeError(
                    f"[CellRank2Adapter] scenario_params.train_times was supplied "
                    f"but adata.obs[{time_key!r}] is missing."
                )
            before = adata.n_obs
            adata, _ = split_adata_by_timepoints(
                adata,
                time_key=time_key,
                train_times=train_times_f,
                heldout_times=None,
                test_includes_start=False,
            )
            self.adata = adata  # keep downstream consistent
            print(
                f"[CellRank2Adapter] scenario_params.train_times applied: "
                f"{before} → {adata.n_obs} cells "
                f"(dropped {before - adata.n_obs}). "
                f"train_times={train_times_f}"
            )

        stm_path = out_dir / "state_transition_matrix.csv"
        edges_path = out_dir / "lineage_graph_edges.csv"

        # --- Check CellRank2 availability ---
        try:
            import cellrank as cr
        except ImportError:
            print(
                "[CellRank2Adapter] cellrank not installed — writing scaffold outputs."
            )
            self._write_scaffold_outputs(
                stm_path, edges_path, scenario_id,
                reason="cellrank2_not_installed",
            )
            return {
                "state_transition_matrix": str(stm_path),
                "lineage_graph_edges": str(edges_path),
            }

        # --- Check WOT availability ---
        try:
            import wot
        except ImportError:
            print(
                "[CellRank2Adapter] wot not installed — writing scaffold outputs."
            )
            self._write_scaffold_outputs(
                stm_path, edges_path, scenario_id,
                reason="wot_not_installed",
            )
            return {
                "state_transition_matrix": str(stm_path),
                "lineage_graph_edges": str(edges_path),
            }

        # --- Validate required obs columns ---
        missing_cols = [
            col for col in [time_key, cell_state_key]
            if col not in adata.obs.columns
        ]
        if missing_cols:
            raise RuntimeError(
                f"[CellRank2Adapter] Required obs columns not found: {missing_cols}. "
                f"Available obs columns: {list(adata.obs.columns)}"
            )

        ncells_subsample = wot_params.get("ncells_subsample")
        if ncells_subsample:
            adata = self._subsample_per_timepoint(
                adata=adata,
                time_key=time_key,
                max_cells=int(ncells_subsample),
                seed=int(wot_params.get("subsample_seed", 42)),
                stratify_key=cell_state_key,
            )
            self.adata = adata

        try:
            # -----------------------------------------------------------
            # Auxiliary integer time key
            # -----------------------------------------------------------
            # Root cause of "more than two extensions" warnings and KeyError
            # inside _restitch_couplings:
            #   WOT derives coupling filenames directly from the day_field
            #   values in obs, producing names like tmap_0.5_2.0.h5ad and
            #   tmap_16.33_16.67.h5ad. Python/OS path libraries treat each
            #   dot-separated suffix as a file extension, so scanpy/AnnData
            #   warns repeatedly. More critically, CellRank's coupling parser
            #   splits on "." to recover time-point values; decimal-valued
            #   names produce ambiguous splits (e.g. "16.33" → ["16", "33"])
            #   causing KeyError when the parser tries to look up a category
            #   by a reconstructed value that does not exist.
            #
            # Fix: map every float time value to a decimal-free integer code
            # (multiply by 100, round to int) BEFORE passing to WOT. This
            # produces safe filenames (tmap_50_200.h5ad, tmap_1600_1633.h5ad).
            # The same integer codes are used as the CellRank time key so that
            # CellRank's filename reconstruction matches the files on disk.
            # The original abs_day column is never modified.
            KERNEL_TIME_KEY = f"{time_key}_cr_code"

            # Sort unique float times to establish a stable, reproducible mapping.
            numeric_times = sorted(adata.obs[time_key].unique())
            # Multiply by 100 and round to int → guaranteed no decimal point.
            time_to_code = {
                t: int(round(t * 100)) for t in numeric_times
            }
            code_to_time = {v: k for k, v in time_to_code.items()}
            print(
                f"[CellRank2Adapter] Integer time-code mapping "
                f"({time_key!r} → {KERNEL_TIME_KEY!r}):"
            )
            for t, c in sorted(time_to_code.items()):
                print(f"  {t} -> {c}")

            # -----------------------------------------------------------
            # Step 1: Build WOT transport maps using integer-coded time key
            # -----------------------------------------------------------
            print(
                f"[CellRank2Adapter] Building WOT transport maps "
                f"(day_field={KERNEL_TIME_KEY!r}, "
                f"epsilon={wot_params.get('epsilon', 0.05)}) ..."
            )

            # Wipe and recreate tmap dir so no stale decimal-named files remain.
            tmap_dir = out_dir / "tmaps_cr2"
            if tmap_dir.exists():
                import shutil
                try:
                    shutil.rmtree(tmap_dir)
                except (PermissionError, OSError) as _rmtree_err:
                    # On some filesystems (e.g. CIFS/SMB mounts) rmtree may be
                    # denied. Fall back to removing individual files so that
                    # stale .h5ad coupling files from previous runs are removed;
                    # the directory itself can remain.
                    print(
                        f"[CellRank2Adapter] Warning: rmtree({tmap_dir}) failed "
                        f"({_rmtree_err!s}). Falling back to per-file removal."
                    )
                    for _f in tmap_dir.iterdir():
                        try:
                            _f.unlink()
                        except OSError:
                            pass
            tmap_dir.mkdir(parents=True, exist_ok=True)
            tmap_prefix = str(tmap_dir / "tmap")

            # WOT copy: inject integer codes + uniform growth rate.
            adata_wot = adata.copy()
            adata_wot.obs[KERNEL_TIME_KEY] = (
                adata_wot.obs[time_key].map(time_to_code).astype(int)
            )
            adata_wot.obs["cell_growth_rate"] = 1.0

            ot_model = wot.ot.OTModel(
                adata_wot,
                day_field=KERNEL_TIME_KEY,          # integer codes → safe filenames
                growth_rate_field="cell_growth_rate",
                epsilon=wot_params.get("epsilon", 0.05),
                lambda1=wot_params.get("lambda1", 1.0),
                lambda2=wot_params.get("lambda2", 50.0),
                local_pca=wot_params.get("local_pca", 30),
                growth_iters=wot_params.get("growth_iters", 3),
            )
            print(f"[CellRank2Adapter] Timepoints (integer codes): {ot_model.timepoints}")

            ot_model.compute_all_transport_maps(
                tmap_out=tmap_prefix, overwrite=True
            )
            n_pairs = max(0, len(ot_model.timepoints) - 1)
            print(
                f"[CellRank2Adapter] Transport maps computed ({n_pairs} consecutive pairs). "
                f"Coupling files in: {tmap_dir}"
            )

            # Release the WOT OT model before building the CellRank2 kernel.
            # wot.ot.OTModel holds an internal reference to adata_wot; deleting
            # it here lets Python reclaim those data structures immediately,
            # avoiding a three-copy peak (original + adata_wot + adata_cr).
            del ot_model

            # -----------------------------------------------------------
            # Step 2: Build CellRank2 RealTimeKernel from WOT coupling dir
            # -----------------------------------------------------------
            # CellRank requires obs[time_key] to be categorical (not float/int).
            # Use the SAME integer codes as the WOT stage so that CellRank's
            # filename reconstruction (from category labels) matches the files
            # already on disk. Use ordered=True so category integer positions
            # follow chronological order — required by _restitch_couplings.
            #
            # Memory optimisation: reuse adata_wot rather than making a second
            # full copy of the 2.4 GB AnnData.  The only structural difference
            # between what WOT needed and what CellRank2 needs is the dtype of
            # KERNEL_TIME_KEY (int vs. ordered categorical) and the presence of
            # cell_growth_rate.  We convert in-place and rebind the variable;
            # no second adata.copy() is required.  Peak RAM goes from ~3× to ~2×
            # the dataset size.
            ordered_codes = sorted(time_to_code.values())
            adata_wot.obs[KERNEL_TIME_KEY] = pd.Categorical(
                adata_wot.obs[KERNEL_TIME_KEY],
                categories=ordered_codes,
                ordered=True,
            )
            adata_cr = adata_wot
            del adata_wot  # rebind only — the object lives on as adata_cr
            print(
                f"[CellRank2Adapter] CellRank time key {KERNEL_TIME_KEY!r} set to "
                f"ordered categorical: {ordered_codes}"
            )
            print(
                f"[CellRank2Adapter] Building RealTimeKernel from coupling dir "
                f"{tmap_dir!r} ..."
            )
            rtk = cr.kernels.RealTimeKernel.from_wot(
                adata_cr,
                path=str(tmap_dir),
                time_key=KERNEL_TIME_KEY,
            )
            rtk.compute_transition_matrix()

            T = rtk.transition_matrix
            print(f"[CellRank2Adapter] Cell-level transition matrix shape: {T.shape}")

            # -----------------------------------------------------------
            # Step 3: Aggregate cell-level transitions to state level
            # -----------------------------------------------------------
            print(
                f"[CellRank2Adapter] Aggregating to state level "
                f"(cell_state_key={cell_state_key!r}) ..."
            )
            stm_norm, edges_df = _aggregate_cr2_to_state_level(
                adata, T, cell_state_key
            )
            print(
                f"[CellRank2Adapter] State-transition matrix: "
                f"{stm_norm.shape[0]} states × {stm_norm.shape[1]} states, "
                f"{len(edges_df)} nonzero edges."
            )

            # -----------------------------------------------------------
            # Step 4: Write outputs
            # -----------------------------------------------------------
            stm_norm.to_csv(stm_path)
            edges_df.to_csv(edges_path, index=False)

            self.write_run_metadata(
                scenario_id=scenario_id,
                status="completed",
                notes=(
                    f"CellRank2 RealTimeKernel.from_wot() with WOT transport maps. "
                    f"time_key={time_key!r}, cell_state_key={cell_state_key!r}, "
                    f"wot_params={wot_params}. "
                    f"State matrix: {stm_norm.shape[0]} × {stm_norm.shape[1]}, "
                    f"{len(edges_df)} nonzero edges."
                ),
            )
            print(f"[CellRank2Adapter] Done. Outputs written to {out_dir}")

        except Exception as exc:
            print(f"[CellRank2Adapter] Runtime error during kernel computation: {exc}")
            traceback.print_exc()
            self._write_scaffold_outputs(
                stm_path, edges_path, scenario_id,
                reason=f"runtime_error: {exc!s}",
            )

        return {
            "state_transition_matrix": str(stm_path),
            "lineage_graph_edges": str(edges_path),
        }

    @staticmethod
    def _subsample_per_timepoint(
        adata,
        time_key: str,
        max_cells: int,
        seed: int,
        stratify_key: str | None = None,
    ):
        """Deterministically cap cells per timepoint before dense WOT maps."""
        if max_cells <= 0:
            return adata

        rng = np.random.default_rng(seed)
        obs = adata.obs
        selected: list = []
        summaries: list[str] = []

        for time_value, time_idx in obs.groupby(time_key, observed=True).groups.items():
            time_idx = pd.Index(time_idx)
            n_cells = len(time_idx)
            if n_cells <= max_cells:
                selected.extend(time_idx.tolist())
                summaries.append(f"{time_value}: kept {n_cells}")
                continue

            if stratify_key and stratify_key in obs.columns:
                chosen = CellRank2Adapter._stratified_sample_index(
                    obs=obs.loc[time_idx],
                    stratify_key=stratify_key,
                    max_cells=max_cells,
                    rng=rng,
                )
            else:
                chosen = rng.choice(time_idx.to_numpy(), size=max_cells, replace=False)

            selected.extend(pd.Index(chosen).tolist())
            summaries.append(f"{time_value}: {n_cells} -> {max_cells}")

        selected_set = set(selected)
        mask = obs.index.isin(selected_set)
        subset = adata[mask]
        filtered = (
            subset.to_memory()
            if getattr(subset, "isbacked", False)
            else subset.copy()
        )
        print(
            "[CellRank2Adapter] ncells_subsample applied per timepoint "
            f"(max={max_cells}, seed={seed}, stratify_key={stratify_key!r}): "
            + "; ".join(summaries)
        )
        print(
            f"[CellRank2Adapter] Subsampled cell count: {adata.n_obs} -> {filtered.n_obs}"
        )
        return filtered

    @staticmethod
    def _stratified_sample_index(
        obs: pd.DataFrame,
        stratify_key: str,
        max_cells: int,
        rng: np.random.Generator,
    ) -> pd.Index:
        groups = [
            pd.Index(group_idx)
            for _, group_idx in obs.groupby(stratify_key, observed=True).groups.items()
        ]
        counts = np.array([len(idx) for idx in groups], dtype=float)

        if len(groups) > max_cells:
            keep_groups = np.argsort(counts)[::-1][:max_cells]
            return pd.Index([rng.choice(groups[i].to_numpy()) for i in keep_groups])

        raw = counts / counts.sum() * max_cells
        alloc = np.maximum(1, np.floor(raw).astype(int))
        alloc = np.minimum(alloc, counts.astype(int))

        while alloc.sum() < max_cells:
            capacity = counts.astype(int) - alloc
            if capacity.max() <= 0:
                break
            residual = raw - np.floor(raw)
            scores = np.where(capacity > 0, residual, -1.0)
            i = int(np.argmax(scores))
            alloc[i] += 1

        while alloc.sum() > max_cells:
            candidates = np.where(alloc > 1)[0]
            if len(candidates) == 0:
                break
            residual = raw[candidates] - np.floor(raw[candidates])
            i = int(candidates[np.argmin(residual)])
            alloc[i] -= 1

        chosen = []
        for idx, n_take in zip(groups, alloc):
            chosen.extend(rng.choice(idx.to_numpy(), size=int(n_take), replace=False))
        return pd.Index(chosen)

    # ------------------------------------------------------------------
    # Scaffold fallback
    # ------------------------------------------------------------------

    def _write_scaffold_outputs(
        self,
        stm_path: Path,
        edges_path: Path,
        scenario_id: str,
        reason: str,
    ) -> None:
        """
        Write empty placeholder outputs and run_metadata with a scaffold status.

        Empty outputs cause eval_lineage.py to report
        "completed_with_empty_predictions" in lineage_metrics.json, which is
        clearly distinguished from a successful "completed" run.
        """
        stm = pd.DataFrame(
            index=pd.Index([], name="source_state"),
            columns=pd.Index([], name="target_state"),
        )
        edges = pd.DataFrame(columns=["source_state", "target_state", "weight"])
        stm.to_csv(stm_path)
        edges.to_csv(edges_path, index=False)
        self.write_run_metadata(
            scenario_id=scenario_id,
            status=f"scaffold_only ({reason})",
            notes=(
                f"CellRank2 adapter fell back to scaffold outputs. "
                f"Reason: {reason}. "
                "Lineage metrics will show 'completed_with_empty_predictions'."
            ),
        )
