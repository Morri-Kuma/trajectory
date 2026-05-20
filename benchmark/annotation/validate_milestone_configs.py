"""validate_milestone_configs.py
================================
Validate milestone-based method config YAML files for Step 5.

Checks performed per config:
  1. YAML parses cleanly.
  2. Required top-level fields present.
  3. ground_truth.provider_id is registered in registry.yaml.
  4. ground_truth.state_key matches the registry entry state_key.
  5. ground_truth.confidence_mode == "all".
  6. dataset.h5ad_path exists on disk.
  7. The milestone obs column (state_key) is present in the h5ad (via h5py).
  8. lineage.reference_graph_path exists on disk.
  9. Reference graph JSON parses; node/edge counts match registry.
  10. lineage.cell_state_key == ground_truth.state_key.
  11. dataset.h5ad_path and state_system.source_h5ad match registry source_h5ad
      when the registry declares one.
  12. Scenario B/C configs have scenario_params with required sub-keys.

Usage
-----
  python benchmark/annotation/validate_milestone_configs.py [config.yaml ...]
  # omit paths to validate all *milestone*.yaml files under benchmark/configs/

Exit codes:  0 = all pass (warnings allowed),  1 = errors found.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Tuple


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


def _load_registry(root: Path) -> Dict[str, Any]:
    reg_path = root / "benchmark" / "ground_truth" / "registry.yaml"
    data = _load_yaml(reg_path)
    return data.get("providers") or {}


def _h5ad_obs_keys(h5ad_path: Path) -> List[str]:
    try:
        import h5py
    except ImportError as exc:
        raise ImportError("h5py required: pip install h5py") from exc
    with h5py.File(h5ad_path, "r") as f:
        if "obs" not in f:
            return []
        return list(f["obs"].keys())


def _load_ref_graph(graph_path: Path) -> Tuple[List[Any], List[Any]]:
    with open(graph_path, encoding="utf-8") as f:
        data = json.load(f)
    return data.get("nodes") or [], data.get("edges") or []


REQUIRED_TOP_KEYS = [
    "run_id", "method", "scenario", "framework_version",
    "ground_truth", "dataset", "output", "lineage", "evaluation",
]
SCENARIO_B_REQUIRED = {"train_times", "heldout_times"}
SCENARIO_C_REQUIRED = {"train_times", "heldout_times", "heldout_interpolation", "heldout_extrapolation"}


def validate_config(
    cfg_path: Path,
    root: Path,
    registry: Dict[str, Any],
) -> Tuple[List[str], List[str]]:
    """Return (errors, warnings) lists for one config file."""
    errors: List[str] = []
    warnings: List[str] = []

    # 1. Parse YAML.
    try:
        cfg = _load_yaml(cfg_path)
    except Exception as exc:
        return [f"YAML parse error: {exc}"], []

    gt = cfg.get("ground_truth") or {}
    dataset = cfg.get("dataset") or {}
    lineage = cfg.get("lineage") or {}
    state_system = cfg.get("state_system") or {}
    scenario = str(cfg.get("scenario", "")).upper()
    provider_id = gt.get("provider_id", "")
    state_key = gt.get("state_key", "")

    # 2. Required top-level keys.
    for key in REQUIRED_TOP_KEYS:
        if key not in cfg:
            errors.append(f"Missing required top-level key: {key!r}")

    # 3. provider_id registered.
    if provider_id not in registry:
        errors.append(f"ground_truth.provider_id {provider_id!r} not in registry.yaml")
        reg_entry: Dict[str, Any] = {}
    else:
        reg_entry = registry[provider_id] or {}

    # 4. state_key matches registry.
    if reg_entry:
        reg_sk = reg_entry.get("state_key", "")
        if state_key != reg_sk:
            errors.append(
                f"ground_truth.state_key {state_key!r} != registry state_key {reg_sk!r}"
            )

    # 5. confidence_mode == all.
    if gt.get("confidence_mode") != "all":
        warnings.append(
            f"ground_truth.confidence_mode is {gt.get('confidence_mode')!r}; milestone providers should use 'all'"
        )

    # 6. h5ad_path exists.
    h5ad_rel = dataset.get("h5ad_path", "")
    h5ad_path = root / h5ad_rel if h5ad_rel else None
    if not h5ad_rel:
        errors.append("dataset.h5ad_path is missing or empty")
    elif h5ad_path and not h5ad_path.exists():
        errors.append(f"dataset.h5ad_path does not exist: {h5ad_rel}")
    elif h5ad_path and state_key:
        # 7. obs column present (via h5py).
        try:
            obs_keys = _h5ad_obs_keys(h5ad_path)
            if state_key not in obs_keys:
                errors.append(
                    f"state_key {state_key!r} not in h5ad obs. Available: {obs_keys}"
                )
        except Exception as exc:
            warnings.append(f"Could not check h5ad obs columns: {exc}")

    # 8 & 9. reference_graph_path.
    graph_rel = lineage.get("reference_graph_path", "")
    graph_path = root / graph_rel if graph_rel else None
    if not graph_rel:
        errors.append("lineage.reference_graph_path is missing or empty")
    elif graph_path and not graph_path.exists():
        errors.append(f"lineage.reference_graph_path does not exist: {graph_rel}")
    else:
        try:
            nodes, edges = _load_ref_graph(graph_path)
            exp_e = reg_entry.get("n_graph_edges")
            exp_n = reg_entry.get("n_states")
            if exp_e is not None and len(edges) != exp_e:
                errors.append(f"Graph has {len(edges)} edges; registry expects {exp_e}")
            if exp_n is not None and len(nodes) != exp_n:
                errors.append(f"Graph has {len(nodes)} nodes; registry expects {exp_n}")
        except Exception as exc:
            errors.append(f"Reference graph JSON error: {exc}")

    # 10. lineage.cell_state_key == ground_truth.state_key.
    lin_sk = lineage.get("cell_state_key", "")
    if lin_sk != state_key:
        errors.append(
            f"lineage.cell_state_key {lin_sk!r} != ground_truth.state_key {state_key!r}"
        )

    # 11. Dataset path consistency with registry source_h5ad.
    if reg_entry:
        reg_source = reg_entry.get("source_h5ad") or ""
        if reg_source:
            if dataset.get("h5ad_path") != reg_source:
                errors.append(
                    f"dataset.h5ad_path {dataset.get('h5ad_path')!r} != "
                    f"registry source_h5ad {reg_source!r}"
                )
            if state_system.get("source_h5ad") != reg_source:
                errors.append(
                    f"state_system.source_h5ad {state_system.get('source_h5ad')!r} != "
                    f"registry source_h5ad {reg_source!r}"
                )

    # 12. Scenario-specific scenario_params.
    sp = cfg.get("scenario_params") or {}
    if scenario == "A" and sp:
        warnings.append("Scenario A config has unexpected scenario_params block")
    elif scenario == "B":
        missing = SCENARIO_B_REQUIRED - set(sp.keys())
        if missing:
            errors.append(f"Scenario B missing scenario_params keys: {sorted(missing)}")
    elif scenario == "C":
        missing = SCENARIO_C_REQUIRED - set(sp.keys())
        if missing:
            errors.append(f"Scenario C missing scenario_params keys: {sorted(missing)}")

    return errors, warnings


def _discover_milestone_configs(root: Path) -> List[Path]:
    return sorted((root / "benchmark" / "configs").glob("*milestone*.yaml"))


def main(argv: List[str] | None = None) -> int:
    import argparse
    parser = argparse.ArgumentParser(
        description="Validate milestone method config YAML files."
    )
    parser.add_argument("configs", nargs="*",
        help="Config YAML paths. If omitted, discovers all *milestone*.yaml under benchmark/configs/.")
    parser.add_argument("--project-root", default=None)
    args = parser.parse_args(argv)

    root = Path(args.project_root) if args.project_root else _project_root()

    if args.configs:
        cfg_paths = [Path(p) for p in args.configs]
    else:
        cfg_paths = _discover_milestone_configs(root)
        if not cfg_paths:
            print("No milestone config files found under benchmark/configs/.")
            return 0

    print(f"Validating {len(cfg_paths)} config file(s) against project root: {root}")
    print()

    try:
        registry = _load_registry(root)
    except Exception as exc:
        print(f"ERROR: Could not load registry.yaml: {exc}")
        return 1

    total_errors = 0
    for cfg_path in cfg_paths:
        errors, warnings = validate_config(cfg_path, root, registry)
        status = "PASS" if not errors else "FAIL"
        print(f"[{status}] {cfg_path.name}")
        for w in warnings:
            print(f"       WARN: {w}")
        for e in errors:
            print(f"       ERROR: {e}")
        total_errors += len(errors)

    print()
    if total_errors == 0:
        print(f"All {len(cfg_paths)} config(s) passed validation.")
        return 0
    else:
        print(f"{total_errors} error(s) found across {len(cfg_paths)} config(s).")
        return 1


if __name__ == "__main__":
    sys.exit(main())
