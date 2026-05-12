"""run_adapter_key_smoke.py
============================
Step 7 metadata-only smoke helper.

Purpose
-------
Prove that ``eval_dispatch.py`` injects the correct ``cell_state_key`` from a
milestone method config WITHOUT running the full adapter / training pipeline.

This script replicates only the config-loading and scenario_config-injection
section of ``eval_dispatch.main()`` (lines 263-315 in eval_dispatch.py), then
writes a ``run_metadata.json`` recording the resolved ``cell_state_key``,
``provider_id``, ``label_mode``, and ``analysis_role`` to the smoke output dir.

This proves the wiring is correct and produces a metadata file that
``validate_adapter_cell_state_key.py`` can validate.

Usage
-----
  python -m benchmark.evaluation.run_adapter_key_smoke \\
      --method scnode \\
      --scenario A \\
      --method-config benchmark/configs/scnode_gse230659_observed_milestone_consensus_v1_A_hvg2000.yaml \\
      --output-dir benchmark/results/smoke/adapter_cell_state_key/scnode_consensus_A

  # Or use positional shorthand for a config file:
  python benchmark/evaluation/run_adapter_key_smoke.py \\
      --method wot --scenario A \\
      --method-config benchmark/configs/wot_gse230659_observed_milestone_consensus_v1_A_hvg2000.yaml \\
      --output-dir benchmark/results/smoke/adapter_cell_state_key/wot_consensus_A
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional


def _project_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in [here.parent, *here.parents]:
        if (candidate / "benchmark").exists() and (
            candidate / "benchmark" / "evaluation"
        ).exists():
            return candidate
    return here.parents[2]


def _load_yaml(path: Path) -> Dict[str, Any]:
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("PyYAML required: pip install pyyaml") from exc
    with open(path, encoding="utf-8-sig") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Expected mapping in {path}")
    return data


def resolve_scenario_config(method_config: Dict[str, Any], root: Path) -> Dict[str, Any]:
    """
    Replicate the eval_dispatch.py injection logic for milestone configs.

    Mirrors eval_dispatch.main() lines 263-315: reads ground_truth spec,
    resolves cell_state_key, and injects all relevant fields into
    scenario_config so adapters receive them via self.scenario_config.

    Returns the populated scenario_config dict.
    """
    from benchmark.ground_truth import load_ground_truth

    ground_truth_spec = load_ground_truth(method_config)
    ground_truth_dict = ground_truth_spec.to_dict()

    cell_state_key = ground_truth_spec.state_key
    reference_graph_path = ground_truth_spec.reference_graph_path
    edge_confidence_mode = ground_truth_spec.confidence_mode
    exclude_uncertain_states = ground_truth_spec.exclude_uncertain_states

    if reference_graph_path in ("null", "~", "", None):
        reference_graph_path = None

    scenario_config: Dict[str, Any] = {}

    dataset_cfg = method_config.get("dataset", {})
    if dataset_cfg.get("id"):
        scenario_config["dataset_id"] = dataset_cfg["id"]
    if dataset_cfg.get("time_key"):
        scenario_config["time_key"] = dataset_cfg["time_key"]
    if cell_state_key:
        scenario_config["cell_state_key"] = cell_state_key
    if method_config.get("result_class"):
        scenario_config["result_class"] = method_config["result_class"]
    if "formal_benchmark" in method_config:
        scenario_config["formal_benchmark"] = method_config["formal_benchmark"]
    if ground_truth_dict:
        scenario_config["ground_truth"] = ground_truth_dict
    scenario_config["exclude_uncertain_states"] = exclude_uncertain_states

    for top_key in ("cellrank2_params", "wot_params", "scnode_params",
                    "prescient_params", "mioflow_params"):
        if method_config.get(top_key):
            scenario_config[top_key] = method_config[top_key]

    sp = method_config.get("scenario_params")
    if sp:
        scenario_config["scenario_params"] = sp

    return scenario_config, ground_truth_dict, reference_graph_path, edge_confidence_mode


def run_key_smoke(
    method: str,
    scenario: str,
    method_config_path: Path,
    output_dir: Path,
    root: Path,
) -> Dict[str, Any]:
    """
    Resolve the method config, build scenario_config, and write run_metadata.json
    proving that the correct cell_state_key would be passed to the adapter.
    """
    method_config = _load_yaml(method_config_path)

    # Replicate eval_dispatch injection logic.
    scenario_config, ground_truth_dict, ref_graph_path, edge_mode = (
        resolve_scenario_config(method_config, root)
    )

    cell_state_key = scenario_config.get("cell_state_key")
    gt = scenario_config.get("ground_truth") or {}
    provider_id = gt.get("provider_id")
    label_mode = gt.get("label_mode")
    analysis_role = gt.get("analysis_role")

    output_dir.mkdir(parents=True, exist_ok=True)

    metadata: Dict[str, Any] = {
        "method": method,
        "scenario": scenario,
        "dataset": scenario_config.get("dataset_id", "unknown"),
        "run_id": method_config.get("run_id", method_config_path.stem),
        "status": "key_smoke_only",
        "smoke_note": (
            "Metadata-only smoke run: proves eval_dispatch would inject the correct "
            "cell_state_key into scenario_config without executing adapter training. "
            "Heavy dependencies (scanpy, wot, cellrank, torch) are NOT imported."
        ),
        "cell_state_key": cell_state_key,
        "provider_id": provider_id,
        "label_mode": label_mode,
        "analysis_role": analysis_role,
        "reference_graph_path": ref_graph_path,
        "edge_confidence_mode": edge_mode,
        "formal_benchmark": scenario_config.get("formal_benchmark", False),
        "result_class": scenario_config.get("result_class"),
        "method_config_path": str(method_config_path),
        "scenario_config_keys_injected": sorted(scenario_config.keys()),
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    meta_path = output_dir / "run_metadata.json"
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(metadata, f, indent=2)

    print(f"[key_smoke] method          : {method}")
    print(f"[key_smoke] scenario        : {scenario}")
    print(f"[key_smoke] cell_state_key  : {cell_state_key!r}")
    print(f"[key_smoke] provider_id     : {provider_id!r}")
    print(f"[key_smoke] label_mode      : {label_mode!r}")
    print(f"[key_smoke] analysis_role   : {analysis_role!r}")
    print(f"[key_smoke] run_metadata    : {meta_path}")

    return metadata


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Step 7 metadata-only adapter key smoke runner."
    )
    parser.add_argument("--method", required=True, help="Method ID (scnode/wot/cellrank2)")
    parser.add_argument("--scenario", required=True, help="Scenario ID (A-F)")
    parser.add_argument("--method-config", required=True,
                        help="Path to milestone method config YAML.")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for run_metadata.json output.")
    args = parser.parse_args(argv)

    root = _project_root()
    cfg_path = Path(args.method_config)
    if not cfg_path.is_absolute():
        cfg_path = root / cfg_path
    out_dir = Path(args.output_dir)
    if not out_dir.is_absolute():
        out_dir = root / out_dir

    try:
        meta = run_key_smoke(
            method=args.method,
            scenario=args.scenario,
            method_config_path=cfg_path,
            output_dir=out_dir,
            root=root,
        )
    except Exception as exc:
        print(f"[key_smoke] ERROR: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

    # Fail fast if cell_state_key is wrong.
    csk = meta.get("cell_state_key")
    forbidden = {"scgpt_pseudostate_provisional", "scTimeBench_cell_type"}
    if csk is None:
        print("[key_smoke] FAIL: cell_state_key is None", file=sys.stderr)
        return 1
    if csk in forbidden:
        print(
            f"[key_smoke] FAIL: cell_state_key={csk!r} is a legacy fallback. "
            "Milestone configs must inject a milestone state key.",
            file=sys.stderr,
        )
        return 1

    print(f"[key_smoke] PASS: cell_state_key={csk!r} is not a legacy fallback.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
