"""
future_model_adapter.py
Template adapter for future generative / forecasting models.
Framework reference: docs/framework/experimental_framework_v2.md §10.3

This template is for models that CAN generate projected cells at unseen future
time points (e.g., neural ODE models, flow-based trajectory models, diffusion
models for single-cell data).

Such models are eligible for all three benchmark dimensions:
  - Forecast Accuracy
  - Embedding Coherence
  - Lineage Fidelity

Rename this file and fill in the method_id and implementations
when a new generative model is added to the benchmark.
"""

import numpy as np
import pandas as pd
import json
from pathlib import Path

from .base_adapter import BaseAdapter


class FutureModelAdapter(BaseAdapter):
    """
    Template adapter for future generative / forecasting models.

    Eligible benchmark dimensions:
      - Forecast Accuracy
      - Embedding Coherence
      - Lineage Fidelity

    To use this template:
      1. Rename the class (e.g., NeuralODEAdapter, TrajectoryNetAdapter).
      2. Set method_id to the model's identifier string.
      3. Implement _run_forecast_accuracy_impl().
      4. Implement _run_embedding_coherence_impl().
      5. Implement _run_lineage_fidelity_impl().
    """

    # --- Capability flags ---
    # These are True because this template is for methods that can
    # generate projected cells at unseen time points.
    supports_unseen_timepoint_projection: bool = True
    supports_lineage_inference: bool = True

    @property
    def method_id(self) -> str:
        # Override with the actual model identifier.
        raise NotImplementedError("Set method_id in the subclass.")

    # ------------------------------------------------------------------
    # Forecast Accuracy implementation
    # ------------------------------------------------------------------

    def _run_forecast_accuracy_impl(self, scenario_id: str) -> dict:
        """
        Project cells from time t to unseen time t+1 and evaluate forecast accuracy.

        Metrics (per v2 §7.3):
          - Wasserstein Distance
          - Gaussian MMD
          - Energy Distance MMD
          - Hausdorff Loss

        Produces:
          - projected_expression.npy
          - forecast_metrics.json
          - per_timepoint_forecast_metrics.csv
        """
        out_dir = self.output_dir

        # TODO: implement projection and metric computation.
        # projected_expr = model.project(adata, from_time=t, to_time=t+1)
        # metrics = compute_forecast_metrics(projected_expr, observed_expr_at_t1)

        proj_path = out_dir / "projected_expression.npy"
        np.save(proj_path, np.array([]))  # placeholder

        metrics = {
            "wasserstein_distance": None,
            "gaussian_mmd": None,
            "energy_distance_mmd": None,
            "hausdorff_loss": None,
        }
        metrics_path = out_dir / "forecast_metrics.json"
        with open(metrics_path, "w") as f:
            json.dump(metrics, f, indent=2)

        per_tp = pd.DataFrame(columns=["timepoint", "wasserstein_distance",
                                        "gaussian_mmd", "energy_distance_mmd",
                                        "hausdorff_loss"])
        per_tp_path = out_dir / "per_timepoint_forecast_metrics.csv"
        per_tp.to_csv(per_tp_path, index=False)

        return {
            "projected_expression": str(proj_path),
            "forecast_metrics": str(metrics_path),
            "per_timepoint_forecast_metrics": str(per_tp_path),
        }

    # ------------------------------------------------------------------
    # Embedding Coherence implementation
    # ------------------------------------------------------------------

    def _run_embedding_coherence_impl(self, scenario_id: str) -> dict:
        """
        Evaluate whether projected cells preserve biologically meaningful structure.

        Metrics (per v2 §8.3):
          - Adjusted Rand Index (ARI)
          - Average normalized classifier entropy

        Produces:
          - projected_embedding.npy
          - embedding_metrics.json
          - projected_cluster_labels.csv
        """
        out_dir = self.output_dir

        # TODO: embed projected cells and compute coherence metrics.

        emb_path = out_dir / "projected_embedding.npy"
        np.save(emb_path, np.array([]))  # placeholder

        metrics = {
            "adjusted_rand_index": None,
            "avg_normalized_classifier_entropy": None,
        }
        metrics_path = out_dir / "embedding_metrics.json"
        with open(metrics_path, "w") as f:
            json.dump(metrics, f, indent=2)

        labels = pd.DataFrame(columns=["cell_id", "projected_cluster_label"])
        labels_path = out_dir / "projected_cluster_labels.csv"
        labels.to_csv(labels_path, index=False)

        return {
            "projected_embedding": str(emb_path),
            "embedding_metrics": str(metrics_path),
            "projected_cluster_labels": str(labels_path),
        }

    # ------------------------------------------------------------------
    # Lineage Fidelity implementation
    # ------------------------------------------------------------------

    def _run_lineage_fidelity_impl(self, scenario_id: str) -> dict:
        """
        Infer lineage structure and produce state-transition outputs.

        Steps (per v2 §9.4):
          1. Infer transitions across benchmark time intervals.
          2. Aggregate to cell-state level.
          3. Construct state-transition matrix.
          4. Convert to lineage graph.

        Produces:
          - state_transition_matrix.csv
          - lineage_graph_edges.csv
          - lineage_metrics.json (populated by eval_lineage.py)
        """
        out_dir = self.output_dir

        # TODO: implement lineage inference for the generative model.

        stm = pd.DataFrame(index=pd.Index([], name="source_state"),
                           columns=pd.Index([], name="target_state"))
        stm_path = out_dir / "state_transition_matrix.csv"
        stm.to_csv(stm_path)

        edges = pd.DataFrame(columns=["source_state", "target_state", "weight"])
        edges_path = out_dir / "lineage_graph_edges.csv"
        edges.to_csv(edges_path, index=False)

        return {
            "state_transition_matrix": str(stm_path),
            "lineage_graph_edges": str(edges_path),
        }
