"""
scnode_adapter.py
Dispatcher adapter for scNODE in the scTimeBench-aligned benchmark.

scNODE is the first projection-capable method in this project, so it is
eligible for all three framework dimensions:
  - Forecast Accuracy
  - Embedding Coherence
  - Lineage Fidelity

The adapter reuses the core implementation in ``benchmark/methods/scNODE/run.py``
and avoids Windows-specific paths.  Shirokane runs should set paths through the
method config and/or TRAJ_PROJECT_ROOT, as the standalone runner already does.
"""

from __future__ import annotations

import json
import time
import traceback
from pathlib import Path

import numpy as np

from .base_adapter import BaseAdapter
from benchmark.shared.dataset.preprocessors.scenario_timepoint_split import (
    split_adata_by_timepoints,
)

DEFAULT_N_SIM_CELLS_CAP = 2000


class ScNODEAdapter(BaseAdapter):
    """
    Adapter for scNODE.

    The dispatcher calls each dimension separately.  Training scNODE three
    times would be wasteful, so the adapter has a single _ensure_run() method
    that trains/loads the model once and writes all required output files.
    """

    supports_unseen_timepoint_projection: bool = True
    supports_lineage_inference: bool = True

    @property
    def method_id(self) -> str:
        return "scnode"

    def _run_lineage_fidelity_impl(self, scenario_id: str) -> dict:
        self._ensure_run(scenario_id)
        return {
            "state_transition_matrix": str(self.output_dir / "state_transition_matrix.csv"),
            "lineage_graph_edges": str(self.output_dir / "lineage_graph_edges.csv"),
        }

    def _run_forecast_accuracy_impl(self, scenario_id: str) -> dict:
        self._ensure_run(scenario_id)
        return {
            "projected_expression": str(self.output_dir / "projected_expression.npy"),
            "forecast_metrics": str(self.output_dir / "forecast_metrics.json"),
            "per_timepoint_forecast_metrics": str(
                self.output_dir / "per_timepoint_forecast_metrics.csv"
            ),
        }

    def _run_embedding_coherence_impl(self, scenario_id: str) -> dict:
        self._ensure_run(scenario_id)
        return {
            "projected_embedding": str(self.output_dir / "projected_embedding.npy"),
            "embedding_metrics": str(self.output_dir / "embedding_metrics.json"),
            "projected_cluster_labels": str(
                self.output_dir / "projected_cluster_labels.csv"
            ),
        }

    def _ensure_run(self, scenario_id: str) -> None:
        """Train/load scNODE once and write all benchmark outputs."""
        if getattr(self, "_has_run", False):
            return

        self.start_timer()
        t0 = time.time()
        status = "failed"
        err_notes = ""

        time_key = self.scenario_config.get("time_key", "abs_day")
        dataset_id = self.scenario_config.get("dataset_id", "GSE230659")
        cell_state_key = self.scenario_config.get(
            "cell_state_key", "final_milestone_label_coarse"
        )
        scnode_cfg = self.scenario_config.get("scnode_params", {}) or {}
        scenario_params = self.scenario_config.get("scenario_params", {}) or {}
        train_times_raw = scenario_params.get("train_times")
        heldout_times = [float(t) for t in scenario_params.get("heldout_times", [])]
        train_times = [float(t) for t in train_times_raw] if train_times_raw else None

        self.output_dir.mkdir(parents=True, exist_ok=True)
        model_cache = self.output_dir / "trained_scnode_model.pth"

        try:
            from benchmark.methods.scNODE.run import (
                prepare_data,
                run_embedding_coherence,
                run_forecast_accuracy,
                run_lineage_fidelity,
                train_or_load,
            )

            full_adata = self.adata
            self._validate_inputs(full_adata, time_key, cell_state_key)

            train_adata = self._materialize_training_adata(
                full_adata,
                time_key=time_key,
                train_times=train_times,
            )
            # Downstream lineage baseline must use the same training cell
            # universe as the method.
            self.adata = train_adata

            train_data, train_tps, train_unique_tps, _, _, _ = prepare_data(
                train_adata,
                time_key,
                train_times,
            )
            n_genes = train_data[0].shape[1]

            n_sim_cells_cfg = scnode_cfg.get("n_sim_cells")
            if n_sim_cells_cfg is not None:
                n_sim_cells = int(n_sim_cells_cfg)
            else:
                n_sim_cells = min(
                    int(train_data[0].shape[0]),
                    int(scnode_cfg.get("n_sim_cells_cap", DEFAULT_N_SIM_CELLS_CAP)),
                )
            metric_sample_cells = int(scnode_cfg.get("metric_sample_cells", 1000))
            seed = int(scnode_cfg.get("seed", 42))
            print(
                f"[ScNODEAdapter] n_genes={n_genes}, n_sim_cells={n_sim_cells}, "
                f"metric_sample_cells={metric_sample_cells}, "
                f"train_times={train_times if train_times else 'all'}, "
                f"heldout_times={heldout_times}"
            )

            model = train_or_load(
                train_data,
                train_tps,
                n_genes,
                scnode_cfg,
                model_cache,
            )

            if heldout_times:
                full_eval_adata = self._materialize_full_adata(full_adata)
            else:
                full_eval_adata = train_adata

            run_forecast_accuracy(
                model=model,
                adata_full=full_eval_adata,
                time_key=time_key,
                all_unique_tps=train_unique_tps,
                heldout_tps=heldout_times,
                n_sim_cells=n_sim_cells,
                output_dir=self.output_dir,
                metric_sample_cells=metric_sample_cells,
                seed=seed,
            )
            run_embedding_coherence(
                model=model,
                adata_full=full_eval_adata,
                time_key=time_key,
                all_unique_tps=train_unique_tps,
                heldout_tps=heldout_times,
                n_sim_cells=n_sim_cells,
                cell_state_key=cell_state_key,
                output_dir=self.output_dir,
            )
            run_lineage_fidelity(
                model=model,
                adata_train=train_adata,
                time_key=time_key,
                train_unique_tps=train_unique_tps,
                cell_state_key=cell_state_key,
                output_dir=self.output_dir,
            )

            status = "completed"
            self._has_run = True

        except Exception:
            err_notes = traceback.format_exc()
            print(f"[ScNODEAdapter] ERROR:\n{err_notes}")
            raise

        finally:
            elapsed = time.time() - t0
            _gt = self.scenario_config.get("ground_truth") or {}
            metadata = {
                "method": "scnode",
                "dataset": dataset_id,
                "scenario": scenario_id,
                "capability_flags": {
                    "supports_unseen_timepoint_projection": True,
                    "supports_lineage_inference": True,
                },
                "dimensions_executed": [
                    "forecast_accuracy",
                    "embedding_coherence",
                    "lineage_fidelity",
                ],
                "runtime_seconds": round(elapsed, 2),
                "status": status,
                "time_key": time_key,
                "cell_state_key": cell_state_key,
                "provider_id": _gt.get("provider_id") or None,
                "label_mode": _gt.get("label_mode") or None,
                "analysis_role": _gt.get("analysis_role") or None,
                "train_times": train_times if train_times else "all",
                "heldout_times": heldout_times,
                "n_sim_cells": locals().get("n_sim_cells"),
                "metric_sample_cells": locals().get("metric_sample_cells"),
                "notes": err_notes,
            }
            with open(self.output_dir / "run_metadata.json", "w", encoding="utf-8") as f:
                json.dump(metadata, f, indent=2)

    @staticmethod
    def _validate_inputs(adata, time_key: str, cell_state_key: str) -> None:
        missing = [
            key for key in (time_key, cell_state_key)
            if key not in adata.obs.columns
        ]
        if missing:
            raise RuntimeError(
                f"[ScNODEAdapter] Required obs columns not found: {missing}. "
                f"Available obs columns: {list(adata.obs.columns)}"
            )

    @staticmethod
    def _materialize_training_adata(adata, time_key: str, train_times: list | None):
        if not train_times:
            return adata.to_memory() if getattr(adata, "isbacked", False) else adata

        before = adata.n_obs
        train_adata, _ = split_adata_by_timepoints(
            adata,
            time_key=time_key,
            train_times=train_times,
            heldout_times=None,
            test_includes_start=False,
        )
        print(
            f"[ScNODEAdapter] scenario_params.train_times applied: "
            f"{before} -> {train_adata.n_obs} cells "
            f"(dropped {before - train_adata.n_obs})."
        )
        return train_adata

    @staticmethod
    def _materialize_full_adata(adata):
        return adata.to_memory() if getattr(adata, "isbacked", False) else adata
