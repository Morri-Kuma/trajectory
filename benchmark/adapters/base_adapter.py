"""
base_adapter.py
Base class for all method adapters in the scTimeBench-aligned benchmark.
Framework reference: experimental framework v2.md §6

Each method adapter must:
  1. Declare capability flags explicitly.
  2. Implement run_lineage_fidelity() if supports_lineage_inference = True.
  3. Implement run_forecast_accuracy() and run_embedding_coherence() only if
     supports_unseen_timepoint_projection = True.
  4. Never override capability flags to force a method into unsupported dimensions.
"""

from abc import ABC, abstractmethod
import time
import json
from pathlib import Path


class BaseAdapter(ABC):
    """
    Abstract base class for benchmark method adapters.

    Subclasses must declare:
      - supports_unseen_timepoint_projection (bool)
      - supports_lineage_inference (bool)

    and implement the methods appropriate for their capability flags.
    """

    # --- Capability flags (must be overridden in each subclass) ---
    supports_unseen_timepoint_projection: bool = NotImplemented
    supports_lineage_inference: bool = NotImplemented

    def __init__(self, adata, scenario_config: dict, output_dir: str):
        """
        Parameters
        ----------
        adata : AnnData
            Benchmark input object. Required obs fields: cell_id, sample_id, time_label.
            Required var fields: gene_symbol. Required matrix: X.
        scenario_config : dict
            Scenario configuration loaded from scenario_observed.yaml or
            scenario_pseudotime.yaml.
        output_dir : str
            Directory where output files will be written.
        """
        self.adata = adata
        self.scenario_config = scenario_config
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._start_time = None
        self._validate_capability_flags()

    def _validate_capability_flags(self):
        if self.supports_unseen_timepoint_projection is NotImplemented:
            raise NotImplementedError(
                f"{self.__class__.__name__} must declare "
                "'supports_unseen_timepoint_projection'."
            )
        if self.supports_lineage_inference is NotImplemented:
            raise NotImplementedError(
                f"{self.__class__.__name__} must declare "
                "'supports_lineage_inference'."
            )

    @property
    def method_id(self) -> str:
        raise NotImplementedError

    @property
    def eligible_dimensions(self) -> list:
        dims = []
        if self.supports_lineage_inference:
            dims.append("lineage_fidelity")
        if self.supports_unseen_timepoint_projection:
            dims.extend(["forecast_accuracy", "embedding_coherence"])
        return dims

    # ------------------------------------------------------------------
    # Lineage Fidelity interface (required if supports_lineage_inference)
    # ------------------------------------------------------------------

    def run_lineage_fidelity(self, scenario_id: str) -> dict:
        """
        Run lineage fidelity inference and return result paths.

        Must produce:
          - state_transition_matrix.csv
          - lineage_graph_edges.csv
          - lineage_metrics.json

        Returns
        -------
        dict with keys: state_transition_matrix, lineage_graph_edges, lineage_metrics
        """
        if not self.supports_lineage_inference:
            raise RuntimeError(
                f"{self.method_id} does not support lineage inference. "
                "Do not call run_lineage_fidelity() on this adapter."
            )
        return self._run_lineage_fidelity_impl(scenario_id)

    def _run_lineage_fidelity_impl(self, scenario_id: str) -> dict:
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement _run_lineage_fidelity_impl()."
        )

    # ------------------------------------------------------------------
    # Forecast Accuracy interface (only for projection-capable methods)
    # ------------------------------------------------------------------

    def run_forecast_accuracy(self, scenario_id: str) -> dict:
        """
        Run forecast accuracy projection and return result paths.

        Must produce:
          - projected_expression.npy
          - forecast_metrics.json
          - per_timepoint_forecast_metrics.csv

        Only callable if supports_unseen_timepoint_projection = True.
        """
        if not self.supports_unseen_timepoint_projection:
            raise RuntimeError(
                f"{self.method_id} does not support unseen-timepoint projection. "
                "Forecast Accuracy is not applicable to this method. "
                "Do not call run_forecast_accuracy() on this adapter."
            )
        return self._run_forecast_accuracy_impl(scenario_id)

    def _run_forecast_accuracy_impl(self, scenario_id: str) -> dict:
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement _run_forecast_accuracy_impl() "
            "if supports_unseen_timepoint_projection = True."
        )

    # ------------------------------------------------------------------
    # Embedding Coherence interface (only for projection-capable methods)
    # ------------------------------------------------------------------

    def run_embedding_coherence(self, scenario_id: str) -> dict:
        """
        Run embedding coherence evaluation and return result paths.

        Must produce:
          - projected_embedding.npy
          - embedding_metrics.json
          - projected_cluster_labels.csv

        Only callable if supports_unseen_timepoint_projection = True.
        """
        if not self.supports_unseen_timepoint_projection:
            raise RuntimeError(
                f"{self.method_id} does not support unseen-timepoint projection. "
                "Embedding Coherence is not applicable to this method. "
                "Do not call run_embedding_coherence() on this adapter."
            )
        return self._run_embedding_coherence_impl(scenario_id)

    def _run_embedding_coherence_impl(self, scenario_id: str) -> dict:
        raise NotImplementedError(
            f"{self.__class__.__name__} must implement _run_embedding_coherence_impl() "
            "if supports_unseen_timepoint_projection = True."
        )

    # ------------------------------------------------------------------
    # Run metadata
    # ------------------------------------------------------------------

    def write_run_metadata(self, scenario_id: str, status: str, notes: str = ""):
        """
        Write run_metadata.json for this method/scenario run.
        Required fields per framework v2 §12.4.

        Automatically records cell_state_key, provider_id, label_mode, and
        analysis_role from self.scenario_config when available, so that
        milestone-based runs are traceable without each adapter having to
        duplicate this logic.
        """
        elapsed = time.time() - self._start_time if self._start_time else None

        # Extract dynamic cell-state / ground-truth provenance fields.
        # These are injected by eval_dispatch.py from the method config when
        # --method-config is provided. Absent = None so the field is present
        # but clearly null rather than silently missing.
        cell_state_key = self.scenario_config.get("cell_state_key")
        gt = self.scenario_config.get("ground_truth") or {}
        provider_id = gt.get("provider_id") or None
        label_mode = gt.get("label_mode") or None
        analysis_role = gt.get("analysis_role") or None

        metadata = {
            "method": self.method_id,
            "dataset": self.adata.uns.get("dataset_id", "unknown"),
            "scenario": scenario_id,
            "capability_flags": {
                "supports_unseen_timepoint_projection": self.supports_unseen_timepoint_projection,
                "supports_lineage_inference": self.supports_lineage_inference,
            },
            "dimensions_executed": self.eligible_dimensions,
            "runtime_seconds": elapsed,
            "status": status,
            "cell_state_key": cell_state_key,
            "provider_id": provider_id,
            "label_mode": label_mode,
            "analysis_role": analysis_role,
            "notes": notes,
        }
        out_path = self.output_dir / "run_metadata.json"
        with open(out_path, "w") as f:
            json.dump(metadata, f, indent=2)
        return str(out_path)

    def start_timer(self):
        self._start_time = time.time()
