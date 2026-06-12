#!/usr/bin/env python
"""Generate runtime configs for the new projection-capable methods.

Adds scIMF, PI-SDE, and Squidiff to the two chemical primaries (GSE178325,
GSE230659) across all six scenarios: observed-time A/B/C (HVG2000 input) and
pseudotime D/E/F (DPT-bin input). Mirrors the existing
benchmark/configs/runtime/*.yaml exactly; observed B/C scenario_params are copied
from the matching scNODE configs so the splits stay identical across methods.

    python scripts/generate_method_expansion_configs.py [--out benchmark/configs/runtime]
"""
from __future__ import annotations

import argparse
import os
import yaml

NEW_METHODS = ["scimf", "pisde", "squiddiff"]

# Sensible default param blocks (tune after vendoring each method on Shirokane).
PARAMS = {
    "scimf": {"latent_dim": 50, "hidden_dim": 64, "n_heads": 4, "n_layers": 2,
              "n_epochs": 100, "batch_size": 64, "lr": 0.001, "seed": 42,
              "n_sim_cells": 2000, "metric_sample_cells": 1000, "use_cuda": True},
    "pisde": {"latent_dim": 50, "hidden_dim": 64, "n_epochs": 100, "sde_steps": 100,
              "lr": 0.001, "physics_weight": 1.0, "seed": 42, "n_sim_cells": 2000,
              "metric_sample_cells": 1000, "use_cuda": True},
    "squiddiff": {"pca_dims": 50, "hidden_dim": 128, "n_epochs": 100,
                  "diffusion_steps": 1000, "batch_size": 64, "lr": 0.0001, "seed": 42,
                  "n_sim_cells": 2000, "metric_sample_cells": 1000, "use_cuda": True},
}

DATASETS = {
    "GSE178325": {"provider": "gse178325_marker_fm_transition_silver_v1", "n_states": 5, "n_edges": 4,
        "hvg2000": "benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad",
        "pseudotime": "benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_marker_fm_transition_silver_HVG2000_pseudotime_benchmark_input.h5ad"},
    "GSE230659": {"provider": "gse230659_marker_fm_transition_silver_v1", "n_states": 4, "n_edges": 3,
        "hvg2000": "benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad",
        "pseudotime": "benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/GSE230659_marker_fm_transition_silver_HVG2000_pseudotime_benchmark_input.h5ad"},
}
PSEUDOTIME_SCENARIOS = {
    "D": {"train_times": [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14], "heldout_times": [3, 7, 11]},
    "E": {"train_times": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], "heldout_times": [10, 11, 12, 13, 14]},
    "F": {"train_times": [0, 1, 2, 4, 5, 6, 8, 9], "heldout_times": [3, 7, 10, 11, 12, 13, 14]},
}
STATE_KEY = "final_milestone_label_coarse"
RT = "benchmark/configs/runtime"


def _load(path):
    with open(path, encoding="utf-8-sig") as f:
        return yaml.safe_load(f)


def _observed_scenario_params(dataset_id, scenario):
    """Copy B/C splits from the matching scNODE config so methods stay aligned."""
    if scenario == "A":
        return None
    p = f"{RT}/scnode_{dataset_id.lower()}_marker_fm_silver_{scenario}_hvg2000_formal.yaml"
    if os.path.exists(p):
        return (_load(p) or {}).get("scenario_params")
    return None


def make(method, dataset_id, scenario, out_dir):
    meta = DATASETS[dataset_id]
    pseudo = scenario in ("D", "E", "F")
    tag = "pseudotime" if pseudo else "hvg2000"
    h5ad = meta["pseudotime"] if pseudo else meta["hvg2000"]
    time_key = "dpt_pseudotime_bin_numeric" if pseudo else "abs_day"
    sp = PSEUDOTIME_SCENARIOS[scenario] if pseudo else _observed_scenario_params(dataset_id, scenario)
    prov = meta["provider"]
    ref = f"benchmark/ground_truth/providers/{prov}/reference_graph.json"
    run_id = f"{method}_{dataset_id.lower()}_marker_fm_silver_{scenario}_{tag}_formal"
    cfg = {
        "run_id": run_id, "method": method, "scenario": scenario,
        "scenario_name": f"{dataset_id} {method} scenario {scenario} official_silver marker_fm_transition ({tag})",
        "framework_version": "v2.1", "result_class": "marker_fm_transition_silver_formal",
        "formal_benchmark": True,
        "ground_truth": {"provider_id": prov, "state_key": STATE_KEY,
                         "confidence_mode": "all", "exclude_uncertain_states": False},
        "dataset": {"id": dataset_id, "h5ad_path": h5ad, "time_key": time_key},
        "output": {"base_dir": f"benchmark/results/{method}/{dataset_id.lower()}_marker_fm_silver_{scenario}_{tag}_formal"},
        f"{method}_params": PARAMS[method],
    }
    if sp:
        cfg["scenario_params"] = sp
    cfg["lineage"] = {"cell_state_key": STATE_KEY, "cell_state_annotation_is_proxy": False,
                      "state_system_version": prov, "reference_graph_path": ref,
                      "edge_confidence_mode": "all", "exclude_uncertain_states": False}
    cfg["evaluation"] = {"run_eval_lineage": True,
                         "eval_script": "benchmark/evaluation/eval_lineage.py", "status": "formal"}
    cfg["cell_state_key"] = STATE_KEY
    cfg["time_key"] = time_key
    cfg["state_system"] = {"version": prov, "status": "frozen_silver_standard_milestone_provider",
                           "n_states": meta["n_states"], "n_graph_edges": meta["n_edges"],
                           "source_h5ad": h5ad, "note": "Frozen marker+FM milestone provider."}
    return run_id, cfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=RT)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    n = 0
    for method in NEW_METHODS:
        for ds in DATASETS:
            for scen in ("A", "B", "C", "D", "E", "F"):
                run_id, cfg = make(method, ds, scen, args.out)
                with open(os.path.join(args.out, f"{run_id}.yaml"), "w", encoding="utf-8") as f:
                    yaml.safe_dump(cfg, f, sort_keys=False, default_flow_style=False)
                n += 1
    print(f"Wrote {n} configs to {args.out} (3 methods x 2 datasets x 6 scenarios)")


if __name__ == "__main__":
    main()
