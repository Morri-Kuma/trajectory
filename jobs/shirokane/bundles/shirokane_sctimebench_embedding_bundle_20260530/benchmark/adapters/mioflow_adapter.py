"""
Dispatcher adapter for MIOFlow.

MIOFlow is treated as a projection-capable generative model. The benchmark
runner trains an ODE in PCA space and inverse-transforms projected cells back
to HVG expression space, matching the scTimeBench treatment of MIOFlow.
"""

from __future__ import annotations

from .base_adapter import BaseAdapter
from benchmark.shared.dataset.preprocessors.scenario_timepoint_split import (
    split_adata_by_timepoints,
)


class MIOFlowAdapter(BaseAdapter):
    supports_unseen_timepoint_projection: bool = True
    supports_lineage_inference: bool = True

    @property
    def method_id(self) -> str:
        return "mioflow"

    def _run_lineage_fidelity_impl(self, scenario_id: str) -> dict:
        result = self._ensure_run(scenario_id)
        return {
            "state_transition_matrix": result["state_transition_matrix"],
            "lineage_graph_edges": result["lineage_graph_edges"],
        }

    def _run_forecast_accuracy_impl(self, scenario_id: str) -> dict:
        result = self._ensure_run(scenario_id)
        return {
            "projected_expression": result["projected_expression"],
            "forecast_metrics": result["forecast_metrics"],
            "per_timepoint_forecast_metrics": result["per_timepoint_forecast_metrics"],
        }

    def _run_embedding_coherence_impl(self, scenario_id: str) -> dict:
        result = self._ensure_run(scenario_id)
        return {
            "embedding": result["embedding"],
            "projected_embedding": result["projected_embedding"],
            "next_timepoint_embedding": result["next_timepoint_embedding"],
            "embedding_metrics": result["embedding_metrics"],
            "projected_cluster_labels": result["projected_cluster_labels"],
        }

    def _ensure_run(self, scenario_id: str) -> dict:
        if getattr(self, "_result", None) is not None:
            return self._result

        from benchmark.methods.MIOFlow.run import run_pipeline

        self._result = run_pipeline(
            adata=self.adata,
            scenario_id=scenario_id,
            scenario_config=self.scenario_config,
            output_dir=self.output_dir,
        )

        scenario_params = self.scenario_config.get("scenario_params", {}) or {}
        train_times = scenario_params.get("train_times")
        time_key = self.scenario_config.get("time_key", "abs_day")
        if train_times:
            self.adata, _ = split_adata_by_timepoints(
                self.adata,
                time_key=time_key,
                train_times=train_times,
                heldout_times=None,
                test_includes_start=False,
            )
        else:
            self.adata = (
                self.adata.to_memory()
                if getattr(self.adata, "isbacked", False)
                else self.adata
            )
        return self._result
