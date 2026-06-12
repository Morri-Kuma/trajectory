"""
_generative_adapter.py
Shared base for projection-capable generative method adapters (scNODE-style).

scNODE, scIMF, PI-SDE, and Squidiff all share the same benchmark contract: train a
generative temporal model once, then write Forecast Accuracy, Embedding Coherence, and
Lineage Fidelity outputs. The only differences are the method id, the config params key,
and which ``benchmark/methods/<Method>/run.py`` module is called. This base captures the
common orchestration (scenario split, representation input, single-train, run_metadata)
so each new method needs only a ~10-line subclass.

Each method's ``run.py`` must expose the same callables as
``benchmark/methods/scNODE/run.py``:
    prepare_data(train_adata, time_key, train_times)
        -> (train_data, train_tps, train_unique_tps, *_)
    train_or_load(train_data, train_tps, n_genes, cfg, model_cache) -> model
    run_forecast_accuracy(model, adata_full, time_key, all_unique_tps, heldout_tps,
                          n_sim_cells, output_dir, metric_sample_cells, seed)
    run_embedding_coherence(model, adata_full, time_key, all_unique_tps, heldout_tps,
                            n_sim_cells, cell_state_key, output_dir)
    run_lineage_fidelity(model, adata_train, time_key, train_unique_tps,
                         cell_state_key, output_dir)
"""
from __future__ import annotations

import importlib
import json
import time
import traceback

from .base_adapter import BaseAdapter
from benchmark.shared.dataset.preprocessors.scenario_timepoint_split import (
    split_adata_by_timepoints,
)

DEFAULT_N_SIM_CELLS_CAP = 2000


class GenerativeProjectionAdapter(BaseAdapter):
    """Generic projection-capable generative adapter. Subclasses set the 3 hooks."""

    supports_unseen_timepoint_projection: bool = True
    supports_lineage_inference: bool = True

    # --- subclass hooks ---
    _method_id: str = NotImplemented            # e.g. "scimf"
    _params_key: str = NotImplemented           # e.g. "scimf_params"
    _run_module: str = NotImplemented           # e.g. "benchmark.methods.scIMF.run"

    @property
    def method_id(self) -> str:
        return self._method_id

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
            "per_timepoint_forecast_metrics": str(self.output_dir / "per_timepoint_forecast_metrics.csv"),
        }

    def _run_embedding_coherence_impl(self, scenario_id: str) -> dict:
        self._ensure_run(scenario_id)
        return {
            "embedding": str(self.output_dir / "embedding.npy"),
            "projected_embedding": str(self.output_dir / "projected_embedding.npy"),
            "next_timepoint_embedding": str(self.output_dir / "next_timepoint_embedding.npy"),
            "embedding_metrics": str(self.output_dir / "embedding_metrics.json"),
            "projected_cluster_labels": str(self.output_dir / "projected_cluster_labels.csv"),
        }

    def _ensure_run(self, scenario_id: str) -> None:
        if getattr(self, "_has_run", False):
            return
        self.start_timer()
        t0 = time.time()
        status = "failed"
        err_notes = ""

        time_key = self.scenario_config.get("time_key", "abs_day")
        dataset_id = self.scenario_config.get("dataset_id", "unknown")
        cell_state_key = self.scenario_config.get("cell_state_key", "final_milestone_label_coarse")
        cfg = self.scenario_config.get(self._params_key, {}) or {}
        sp = self.scenario_config.get("scenario_params", {}) or {}
        train_times_raw = sp.get("train_times")
        heldout_times = [float(t) for t in sp.get("heldout_times", [])]
        train_times = [float(t) for t in train_times_raw] if train_times_raw else None

        self.output_dir.mkdir(parents=True, exist_ok=True)
        model_cache = self.output_dir / f"trained_{self._method_id}_model.pth"
        n_sim_cells = metric_sample_cells = None

        try:
            run = importlib.import_module(self._run_module)
            from benchmark.representations.model_input import representation_input_adata

            full_adata = representation_input_adata(self.adata, self.scenario_config)
            self._validate_inputs(full_adata, time_key, cell_state_key)

            train_adata = self._materialize_training_adata(full_adata, time_key, train_times)
            self.adata = train_adata  # lineage baseline uses the same training universe

            train_data, train_tps, train_unique_tps, *_ = run.prepare_data(
                train_adata, time_key, train_times,
            )
            n_genes = train_data[0].shape[1]

            n_sim_cfg = cfg.get("n_sim_cells")
            n_sim_cells = (int(n_sim_cfg) if n_sim_cfg is not None
                           else min(int(train_data[0].shape[0]),
                                    int(cfg.get("n_sim_cells_cap", DEFAULT_N_SIM_CELLS_CAP))))
            metric_sample_cells = int(cfg.get("metric_sample_cells", 1000))
            seed = int(cfg.get("seed", 42))
            print(f"[{self._method_id}] n_genes={n_genes} n_sim_cells={n_sim_cells} "
                  f"train_times={train_times if train_times else 'all'} heldout={heldout_times}")

            model = run.train_or_load(train_data, train_tps, n_genes, cfg, model_cache)
            full_eval_adata = (self._materialize_full_adata(full_adata)
                               if heldout_times else train_adata)

            run.run_forecast_accuracy(model=model, adata_full=full_eval_adata, time_key=time_key,
                                      all_unique_tps=train_unique_tps, heldout_tps=heldout_times,
                                      n_sim_cells=n_sim_cells, output_dir=self.output_dir,
                                      metric_sample_cells=metric_sample_cells, seed=seed)
            run.run_embedding_coherence(model=model, adata_full=full_eval_adata, time_key=time_key,
                                        all_unique_tps=train_unique_tps, heldout_tps=heldout_times,
                                        n_sim_cells=n_sim_cells, cell_state_key=cell_state_key,
                                        output_dir=self.output_dir)
            run.run_lineage_fidelity(model=model, adata_train=train_adata, time_key=time_key,
                                     train_unique_tps=train_unique_tps, cell_state_key=cell_state_key,
                                     output_dir=self.output_dir)
            status = "completed"
            self._has_run = True
        except Exception:
            err_notes = traceback.format_exc()
            print(f"[{self._method_id}] ERROR:\n{err_notes}")
            raise
        finally:
            _gt = self.scenario_config.get("ground_truth") or {}
            metadata = {
                "method": self._method_id, "dataset": dataset_id, "scenario": scenario_id,
                "capability_flags": {"supports_unseen_timepoint_projection": True,
                                     "supports_lineage_inference": True},
                "dimensions_executed": ["forecast_accuracy", "embedding_coherence", "lineage_fidelity"],
                "runtime_seconds": round(time.time() - t0, 2), "status": status,
                "time_key": time_key, "cell_state_key": cell_state_key,
                "provider_id": _gt.get("provider_id") or None, "label_mode": _gt.get("label_mode") or None,
                "analysis_role": _gt.get("analysis_role") or None,
                "train_times": train_times if train_times else "all", "heldout_times": heldout_times,
                "n_sim_cells": n_sim_cells, "metric_sample_cells": metric_sample_cells, "notes": err_notes,
            }
            with open(self.output_dir / "run_metadata.json", "w", encoding="utf-8") as f:
                json.dump(metadata, f, indent=2)

    @staticmethod
    def _validate_inputs(adata, time_key, cell_state_key):
        missing = [k for k in (time_key, cell_state_key) if k not in adata.obs.columns]
        if missing:
            raise RuntimeError(f"Required obs columns not found: {missing}; "
                               f"have {list(adata.obs.columns)}")

    @staticmethod
    def _materialize_training_adata(adata, time_key, train_times):
        if not train_times:
            return adata.to_memory() if getattr(adata, "isbacked", False) else adata
        train_adata, _ = split_adata_by_timepoints(
            adata, time_key=time_key, train_times=train_times,
            heldout_times=None, test_includes_start=False)
        return train_adata

    @staticmethod
    def _materialize_full_adata(adata):
        return adata.to_memory() if getattr(adata, "isbacked", False) else adata
