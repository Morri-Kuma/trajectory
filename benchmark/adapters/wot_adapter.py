"""
wot_adapter.py
WOT method adapter for the scTimeBench-aligned benchmark.

This adapter is the dispatcher-facing wrapper around the real WOT runner
implementation in ``benchmark/methods/WOT/run.py``.  The standalone runner is
kept for local/HPC jobs; the adapter reuses its core functions so dispatcher
runs produce the same real Lineage Fidelity outputs instead of scaffold files.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from .base_adapter import BaseAdapter
from benchmark.shared.dataset.preprocessors.scenario_timepoint_split import (
    split_adata_by_timepoints,
)
from benchmark.methods.WOT.run import (
    _aggregate_to_state_level,
    _check_wot,
    _prepare_wot_inputs,
    _run_wot_transport,
    _write_outputs,
)


class WOTAdapter(BaseAdapter):
    """
    Adapter for WOT (Waddington Optimal Transport).

    Eligible benchmark dimensions:
      - Lineage Fidelity only

    Not eligible for:
      - Forecast Accuracy
      - Embedding Coherence
    """

    supports_unseen_timepoint_projection: bool = False
    supports_lineage_inference: bool = True

    @property
    def method_id(self) -> str:
        return "wot"

    def _run_lineage_fidelity_impl(self, scenario_id: str) -> dict:
        """
        Run real WOT transport-map inference and state-level aggregation.

        This method mirrors the successful path in
        ``benchmark/methods/WOT/run.py``:
          1. optionally filter to scenario_params.train_times;
          2. prepare WOT time/growth inputs;
          3. compute or reuse transport maps;
          4. aggregate cell-level transport to the state level;
          5. write state_transition_matrix.csv, lineage_graph_edges.csv, and
             run_metadata.json.
        """
        self.start_timer()

        adata = self.adata
        out_dir = self.output_dir

        time_key = self.scenario_config.get("time_key", "time_label")
        cell_state_key = self.scenario_config.get(
            "cell_state_key", "final_milestone_label_coarse"
        )
        wot_params = self.scenario_config.get("wot_params", {}) or {}
        dataset_id = self.scenario_config.get("dataset_id", "unknown")

        scenario_cfg = (
            self.scenario_config.get("scenario_params")
            or self.scenario_config.get("scenario_split")
            or {}
        )
        train_times = scenario_cfg.get("train_times") if scenario_cfg else None

        self._validate_inputs(adata, time_key=time_key, cell_state_key=cell_state_key)

        if train_times:
            adata = self._filter_to_train_times(
                adata=adata,
                time_key=time_key,
                train_times=train_times,
            )
            self.adata = adata

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

        stm_path = out_dir / "state_transition_matrix.csv"
        edges_path = out_dir / "lineage_graph_edges.csv"

        if not _check_wot():
            reason = "wot_not_installed"
            self._write_scaffold_outputs(stm_path, edges_path, scenario_id, reason)
            return {
                "state_transition_matrix": str(stm_path),
                "lineage_graph_edges": str(edges_path),
            }

        growth_source = wot_params.get("growth_rate_source", "uniform")
        skip_if_exists = bool(wot_params.get("skip_tmap_if_exists", False))

        cell_days, growth_rates = _prepare_wot_inputs(
            adata,
            time_key=time_key,
            growth_rate_source=growth_source,
        )
        print(f"[WOTAdapter] Time points: {sorted(np.unique(cell_days))}")

        tmap_model = _run_wot_transport(
            adata,
            cell_days,
            growth_rates,
            wot_params=wot_params,
            output_dir=out_dir,
            time_key=time_key,
            skip_if_exists=skip_if_exists,
        )

        stm, edges = _aggregate_to_state_level(
            adata,
            tmap_model,
            cell_state_key=cell_state_key,
            time_key=time_key,
        )

        elapsed = time.time() - self._start_time if self._start_time else None
        _gt = self.scenario_config.get("ground_truth") or {}
        run_meta = {
            "method": "wot",
            "dataset": dataset_id,
            "scenario": scenario_id,
            "capability_flags": {
                "supports_unseen_timepoint_projection": False,
                "supports_lineage_inference": True,
            },
            "dimensions_executed": ["lineage_fidelity"],
            "runtime_seconds": round(elapsed, 2) if elapsed is not None else None,
            "status": "completed",
            "time_key": time_key,
            "cell_state_key": cell_state_key,
            "provider_id": _gt.get("provider_id") or None,
            "label_mode": _gt.get("label_mode") or None,
            "analysis_role": _gt.get("analysis_role") or None,
            "scenario_params": scenario_cfg,
            "wot_params": wot_params,
            "notes": (
                "WOT adapter reused benchmark/methods/WOT/run.py core functions. "
                "Forecast Accuracy and Embedding Coherence explicitly skipped."
            ),
        }

        _write_outputs(stm, edges, out_dir, run_meta)
        print(
            f"[WOTAdapter] Done. STM shape={stm.shape}, "
            f"edges={len(edges)}, output={out_dir}"
        )

        return {
            "state_transition_matrix": str(stm_path),
            "lineage_graph_edges": str(edges_path),
        }

    @staticmethod
    def _validate_inputs(adata, time_key: str, cell_state_key: str) -> None:
        missing = [
            key for key in (time_key, cell_state_key)
            if key not in adata.obs.columns
        ]
        if missing:
            raise RuntimeError(
                f"[WOTAdapter] Required obs columns not found: {missing}. "
                f"Available obs columns: {list(adata.obs.columns)}"
            )

    @staticmethod
    def _filter_to_train_times(adata, time_key: str, train_times: list):
        train_times_f = [float(t) for t in train_times]
        before = adata.n_obs
        filtered, _ = split_adata_by_timepoints(
            adata,
            time_key=time_key,
            train_times=train_times_f,
            heldout_times=None,
            test_includes_start=False,
        )
        print(
            f"[WOTAdapter] scenario_params.train_times applied: "
            f"{before} -> {filtered.n_obs} cells "
            f"(dropped {before - filtered.n_obs}). train_times={train_times_f}"
        )
        return filtered

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
                chosen = WOTAdapter._stratified_sample_index(
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
            "[WOTAdapter] ncells_subsample applied per timepoint "
            f"(max={max_cells}, seed={seed}, stratify_key={stratify_key!r}): "
            + "; ".join(summaries)
        )
        print(f"[WOTAdapter] Subsampled cell count: {adata.n_obs} -> {filtered.n_obs}")
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

    def _write_scaffold_outputs(
        self,
        stm_path: Path,
        edges_path: Path,
        scenario_id: str,
        reason: str,
    ) -> None:
        """
        Keep the dispatcher failure mode explicit if WOT is unavailable.

        The normal path is no longer scaffolded.  This fallback only exists for
        environments where the WOT package itself cannot be imported.
        """
        stm = pd.DataFrame(
            index=pd.Index([], name="source_state"),
            columns=pd.Index([], name="target_state"),
        )
        edges = pd.DataFrame(columns=["source_state", "target_state", "weight"])
        stm.to_csv(stm_path)
        edges.to_csv(edges_path, index=False)

        _gt = self.scenario_config.get("ground_truth") or {}
        metadata = {
            "method": "wot",
            "dataset": self.scenario_config.get("dataset_id", "unknown"),
            "scenario": scenario_id,
            "capability_flags": {
                "supports_unseen_timepoint_projection": False,
                "supports_lineage_inference": True,
            },
            "dimensions_executed": [],
            "runtime_seconds": (
                round(time.time() - self._start_time, 2)
                if self._start_time else None
            ),
            "status": f"scaffold_only ({reason})",
            "cell_state_key": self.scenario_config.get("cell_state_key"),
            "provider_id": _gt.get("provider_id") or None,
            "label_mode": _gt.get("label_mode") or None,
            "analysis_role": _gt.get("analysis_role") or None,
            "notes": (
                "WOT package is unavailable, so the adapter wrote empty outputs. "
                "Install WOT in the active environment to produce real predictions."
            ),
        }
        with open(self.output_dir / "run_metadata.json", "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2)
