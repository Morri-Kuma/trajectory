"""
eval_dispatch.py
Capability-gated benchmark dispatcher.
Framework reference: experimental framework v2.md §14, Step 1

This dispatcher routes each method to the correct evaluator(s) based on its
declared capability flags. It prevents invalid evaluations by checking flags
before calling any evaluator. No fallback or workaround is applied.

Usage
-----
    python benchmark/evaluation/eval_dispatch.py \
        --method wot \
        --scenario A \
        --adata data/processed/adata_benchmark.h5ad \
        --output-dir benchmark/results/wot/A
"""

import argparse
import importlib
import json
import sys
from pathlib import Path

import anndata


# ------------------------------------------------------------------
# Method registry
# ------------------------------------------------------------------

METHOD_REGISTRY = {
    "wot": "benchmark.adapters.wot_adapter.WOTAdapter",
    "cellrank2": "benchmark.adapters.cellrank2_adapter.CellRank2Adapter",
    # Add future models here:
    # "future_model": "benchmark.adapters.future_model_adapter.FutureModelAdapter",
}


def load_adapter_class(method_id: str):
    """Dynamically import and return the adapter class for a given method."""
    if method_id not in METHOD_REGISTRY:
        raise ValueError(
            f"Unknown method '{method_id}'. "
            f"Registered methods: {list(METHOD_REGISTRY.keys())}"
        )
    module_path, class_name = METHOD_REGISTRY[method_id].rsplit(".", 1)
    module = importlib.import_module(module_path)
    return getattr(module, class_name)


# ------------------------------------------------------------------
# Dispatcher
# ------------------------------------------------------------------

def dispatch(method_id: str, scenario_id: str, adata_path: str,
             scenario_config: dict, output_dir: str,
             reference_graph_path: str = None,
             edge_confidence_mode: str = "all",
             exclude_uncertain_states: bool = False):
    """
    Run the benchmark for one method × scenario.

    Steps:
      1. Load the method adapter.
      2. Check capability flags.
      3. Route to eligible evaluators only.
      4. Write run_metadata.json.

    Parameters
    ----------
    reference_graph_path : str, optional
        Path to the reference lineage graph JSON (scGPT v1 format).
        If provided, Lineage Fidelity metrics are computed against this graph.
        If None (default), metric computation is deferred.
        Typically read from the method config's lineage.reference_graph_path.
    edge_confidence_mode : str
        Edge confidence filter for the reference graph.
        Options: "all" (default), "medium_and_above", "high_only".
        Typically read from the method config's lineage.edge_confidence_mode.
    exclude_uncertain_states : bool
        If True, exclude uncertain-status nodes from the reference graph when
        computing Lineage Fidelity metrics. Default False.
        Typically read from the method config's lineage.exclude_uncertain_states.
    """
    AdapterClass = load_adapter_class(method_id)
    adata = anndata.read_h5ad(adata_path)
    adapter = AdapterClass(adata=adata, scenario_config=scenario_config,
                           output_dir=output_dir)

    print(f"\n{'='*60}")
    print(f"Method       : {method_id}")
    print(f"Scenario     : {scenario_id}")
    print(f"Eligible for : {adapter.eligible_dimensions}")
    print(f"{'='*60}")

    results = {}

    # --- Lineage Fidelity ---
    if adapter.supports_lineage_inference:
        print("[dispatch] Running Lineage Fidelity ...")
        lf_result = adapter.run_lineage_fidelity(scenario_id=scenario_id)
        results["lineage_fidelity"] = lf_result
        print(f"[dispatch] Lineage Fidelity outputs: {lf_result}")

        # Run the lineage evaluator to compute metrics.
        from benchmark.evaluation.eval_lineage import run_lineage_evaluation
        lf_metrics = run_lineage_evaluation(
            state_transition_matrix_path=lf_result["state_transition_matrix"],
            lineage_graph_edges_path=lf_result["lineage_graph_edges"],
            output_dir=output_dir,
            reference_graph_path=reference_graph_path,
            edge_confidence_mode=edge_confidence_mode,
            exclude_uncertain_states=exclude_uncertain_states,
        )
        results["lineage_fidelity"]["metrics"] = lf_metrics
    else:
        print(f"[dispatch] Skipping Lineage Fidelity: {method_id} does not support lineage inference.")

    # --- Forecast Accuracy ---
    if adapter.supports_unseen_timepoint_projection:
        print("[dispatch] Running Forecast Accuracy ...")
        fa_result = adapter.run_forecast_accuracy(scenario_id=scenario_id)
        results["forecast_accuracy"] = fa_result

        from benchmark.evaluation.eval_forecast import run_forecast_evaluation
        fa_metrics = run_forecast_evaluation(
            projected_expression_path=fa_result["projected_expression"],
            adata=adata,
            output_dir=output_dir,
        )
        results["forecast_accuracy"]["metrics"] = fa_metrics
    else:
        print(
            f"[dispatch] Skipping Forecast Accuracy: {method_id} does not support "
            "unseen-timepoint projection (consistent with scTimeBench treatment of OT methods)."
        )

    # --- Embedding Coherence ---
    if adapter.supports_unseen_timepoint_projection:
        print("[dispatch] Running Embedding Coherence ...")
        ec_result = adapter.run_embedding_coherence(scenario_id=scenario_id)
        results["embedding_coherence"] = ec_result

        from benchmark.evaluation.eval_embedding import run_embedding_evaluation
        ec_metrics = run_embedding_evaluation(
            projected_embedding_path=ec_result["projected_embedding"],
            projected_cluster_labels_path=ec_result["projected_cluster_labels"],
            adata=adata,
            output_dir=output_dir,
        )
        results["embedding_coherence"]["metrics"] = ec_metrics
    else:
        print(
            f"[dispatch] Skipping Embedding Coherence: {method_id} does not support "
            "unseen-timepoint projection (consistent with scTimeBench treatment of OT methods)."
        )

    print(f"\n[dispatch] Run complete for {method_id} / {scenario_id}.")
    return results


# ------------------------------------------------------------------
# CLI entry point
# ------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Capability-gated benchmark dispatcher (scTimeBench v2)"
    )
    parser.add_argument("--method", required=True,
                        help="Method ID (e.g. wot, cellrank2)")
    parser.add_argument("--scenario", required=True,
                        help="Scenario ID (A–F)")
    parser.add_argument("--adata", required=True,
                        help="Path to input h5ad (adata_benchmark.h5ad or "
                             "adata_scgpt_annotated.h5ad for scGPT-v1 configs)")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for output files")
    parser.add_argument("--scenario-config", default=None,
                        help="Path to scenario YAML config (optional)")
    parser.add_argument(
        "--method-config", default=None,
        help=(
            "Path to a method-level YAML config "
            "(e.g. cellrank2_gse230659_observed_scgpt_v1.yaml). "
            "When provided, reads lineage.reference_graph_path and "
            "lineage.edge_confidence_mode from the config. "
            "These override the defaults (None and 'all' respectively)."
        ),
    )
    args = parser.parse_args()

    import yaml  # lazy import (PyYAML may not always be installed)

    scenario_config = {}
    if args.scenario_config:
        # Use utf-8-sig to handle optional BOM on Windows; harmless on Linux/Mac.
        with open(args.scenario_config, encoding="utf-8-sig") as f:
            scenario_config = yaml.safe_load(f)

    # Read optional lineage wiring from method-level config.
    reference_graph_path = None
    edge_confidence_mode = "all"
    exclude_uncertain_states = False
    if args.method_config:
        # Use utf-8-sig for the same reason: YAML files may be BOM-encoded on Windows.
        with open(args.method_config, encoding="utf-8-sig") as f:
            method_config = yaml.safe_load(f)
        lineage_cfg = method_config.get("lineage", {})
        rgp = lineage_cfg.get("reference_graph_path")
        # Treat the YAML literal "null" (parsed as None) and empty string as absent.
        if rgp and rgp not in ("null", "~", ""):
            reference_graph_path = rgp
        ecm = lineage_cfg.get("edge_confidence_mode")
        if ecm:
            edge_confidence_mode = ecm
        eus = lineage_cfg.get("exclude_uncertain_states")
        if eus is not None:
            exclude_uncertain_states = bool(eus)

        # Inject method-level dataset/lineage fields into scenario_config so that
        # adapters can read them via self.scenario_config. Adapters receive
        # scenario_config (from --scenario-config) but NOT the method config
        # (from --method-config) directly. Without this injection the adapter
        # cannot know which time/state columns to use.
        dataset_cfg = method_config.get("dataset", {})
        if dataset_cfg.get("time_key") and "time_key" not in scenario_config:
            scenario_config["time_key"] = dataset_cfg["time_key"]
        if lineage_cfg.get("cell_state_key") and "cell_state_key" not in scenario_config:
            scenario_config["cell_state_key"] = lineage_cfg["cell_state_key"]
        if "exclude_uncertain_states" not in scenario_config:
            scenario_config["exclude_uncertain_states"] = exclude_uncertain_states
        # Inject method-specific params (e.g. cellrank2_params) so the adapter
        # can read kernel/WOT configuration from self.scenario_config.
        for top_key in ("cellrank2_params", "wot_params"):
            if method_config.get(top_key) and top_key not in scenario_config:
                scenario_config[top_key] = method_config[top_key]

    dispatch(
        method_id=args.method,
        scenario_id=args.scenario,
        adata_path=args.adata,
        scenario_config=scenario_config,
        output_dir=args.output_dir,
        reference_graph_path=reference_graph_path,
        edge_confidence_mode=edge_confidence_mode,
        exclude_uncertain_states=exclude_uncertain_states,
    )


if __name__ == "__main__":
    main()
