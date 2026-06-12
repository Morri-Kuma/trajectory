#!/usr/bin/env python
"""Generate v2.1 benchmark runtime configs.

Emits scTimeBench-style runtime configs for:
  (1) the two new iPSC datasets GSE298212 (human blood chemical reprogramming) and
      GSE218855 (mouse MEF fast chemical reprogramming): scnode/prescient/mioflow x
      Scenarios A/B/C plus wot/cellrank2 x Scenario A (lineage-only); and
  (2) pseudotime Scenarios D/E/F for the two chemical primaries GSE230659 and
      GSE178325: scnode/prescient/mioflow x D/E/F plus wot/cellrank2 x D.

Configs mirror the existing benchmark/configs/runtime/*.yaml exactly; only the
dataset, scenario, scenario_params, provider, and output fields change. The
method-parameter blocks are copied verbatim from the GSE230659 formal configs.

This generator is deterministic and idempotent: run it locally to populate
benchmark/configs/runtime/, and re-run it on Shirokane after cloning if needed.

    python scripts/generate_v21_benchmark_configs.py [--out benchmark/configs/runtime]
"""
from __future__ import annotations

import argparse
import os
import yaml

# ---------------------------------------------------------------------------
# Method-parameter blocks (verbatim from the GSE230659 formal configs)
# ---------------------------------------------------------------------------
SCNODE_PARAMS = {
    "latent_dim": 50, "drift_latent_size": [50, 50],
    "enc_latent_list": [64, 64], "dec_latent_list": [64, 64],
    "pretrain_iters": 200, "epochs": 10, "iters": 100, "batch_size": 32,
    "lr": 0.001, "latent_coeff": 1.0, "kl_coeff": 0.0, "seed": 42,
    "n_sim_cells": 2000, "n_sim_cells_cap": 2000,
}
PRESCIENT_PARAMS = {
    "num_pcs": 30, "k_dim": 100, "layers": 1, "activation": "softplus",
    "pretrain_epochs": 25, "train_epochs": 100, "train_lr": 0.01, "train_dt": 0.1,
    "train_sd": 0.5, "train_tau": 1.0e-06, "train_batch": 0.03, "train_clip": 0.25,
    "save": 100, "seed": 42, "weight_name": "uniform", "n_sim_cells": 2000,
    "metric_sample_cells": 1000, "assignment_temperature": 1.0, "use_cuda": True,
    "training_protocol": "full_formal_hvg2000",
}
MIOFLOW_PARAMS = {
    "pca_dims": 50, "hidden_dim": 64, "n_epochs": 40, "batch_size": 64,
    "learning_rate": 0.001, "lambda_ot": 1.0, "use_density_loss": False,
    "lambda_density": 5.0, "lambda_energy": 0.01, "energy_time_steps": 10,
    "seed": 42, "n_sim_cells": 2000, "metric_sample_cells": 1000, "use_cuda": True,
}
WOT_PARAMS = {
    "epsilon": 0.05, "lambda1": 1.0, "lambda2": 50.0, "local_pca": 30,
    "growth_iters": 3, "ncells_subsample": None, "growth_rate_source": "uniform",
}
CELLRANK2_PARAMS = {
    "kernel": "RealTimeKernel", "time_key": "abs_day",
    "wot_params": {"epsilon": 0.05, "lambda1": 1.0, "lambda2": 50.0,
                   "local_pca": 30, "growth_iters": 3},
    "gpcca_params": {"n_macrostates": 14, "use_state_labels": True},
}
PARAM_KEY = {
    "scnode": ("scnode_params", SCNODE_PARAMS),
    "prescient": ("prescient_params", PRESCIENT_PARAMS),
    "mioflow": ("mioflow_params", MIOFLOW_PARAMS),
    "wot": ("wot_params", WOT_PARAMS),
    "cellrank2": ("cellrank2_params", CELLRANK2_PARAMS),
}
PROJECTION = ["scnode", "prescient", "mioflow"]
OT = ["wot", "cellrank2"]

# ---------------------------------------------------------------------------
# Dataset metadata
# ---------------------------------------------------------------------------
# Observed-time (A/B/C) new datasets. abs_day timepoints from dataset_inventory.
NEW_DATASETS = {
    "GSE298212": {
        "input": "benchmark/inputs/gse298212_marker_fm_transition_silver_hvg2000/"
                 "GSE298212_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad",
        "provider": "gse298212_marker_fm_transition_silver_v1",
        "n_states": 4, "n_edges": 3,
        "scenarios": {  # A handled by adapter default (no scenario_params)
            "A": None,
            "B": {"train_times": [0.0, 1.0, 3.0], "heldout_times": [6.0, 8.0]},
            "C": {"train_times": [0.0, 3.0, 8.0], "heldout_times": [1.0, 6.0]},
        },
    },
    "GSE218855": {
        "input": "benchmark/inputs/gse218855_marker_fm_transition_silver_hvg2000/"
                 "GSE218855_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad",
        "provider": "gse218855_marker_fm_transition_silver_v1",
        "n_states": 4, "n_edges": 3,
        "scenarios": {
            "A": None,
            "B": {"train_times": [0.0, 4.0, 8.0], "heldout_times": [12.0, 16.0]},
            "C": {"train_times": [0.0, 8.0, 16.0], "heldout_times": [4.0, 12.0]},
        },
    },
}

# Pseudotime (D/E/F) on the two chemical primaries (15 DPT bins, see
# benchmark/configs/scenario_pseudotime.yaml).
PSEUDOTIME_DATASETS = {
    "GSE230659": {
        "input": "benchmark/inputs/gse230659_marker_fm_transition_silver_hvg2000/"
                 "GSE230659_marker_fm_transition_silver_HVG2000_pseudotime_benchmark_input.h5ad",
        "provider": "gse230659_marker_fm_transition_silver_v1", "n_states": 4, "n_edges": 3,
    },
    "GSE178325": {
        "input": "benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/"
                 "GSE178325_marker_fm_transition_silver_HVG2000_pseudotime_benchmark_input.h5ad",
        "provider": "gse178325_marker_fm_transition_silver_v1", "n_states": 5, "n_edges": 4,
    },
}
PSEUDOTIME_TIME_KEY = "dpt_pseudotime_bin_numeric"
PSEUDOTIME_SCENARIOS = {
    "D": {"train_times": [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14], "heldout_times": [3, 7, 11]},
    "E": {"train_times": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9], "heldout_times": [10, 11, 12, 13, 14]},
    "F": {"train_times": [0, 1, 2, 4, 5, 6, 8, 9], "heldout_times": [3, 7, 10, 11, 12, 13, 14]},
}

STATE_KEY = "final_milestone_label_coarse"


def make_config(method, dataset_id, scenario, meta, time_key, scenario_params, tag):
    prov = meta["provider"]
    ref_graph = f"benchmark/ground_truth/providers/{prov}/reference_graph.json"
    run_id = f"{method}_{dataset_id.lower()}_marker_fm_silver_{scenario}_{tag}_formal"
    out_dir = f"benchmark/results/{method}/{dataset_id.lower()}_marker_fm_silver_{scenario}_{tag}_formal"
    cfg = {
        "run_id": run_id,
        "method": method,
        "scenario": scenario,
        "scenario_name": f"{dataset_id} {method} scenario {scenario} official_silver marker_fm_transition ({tag})",
        "framework_version": "v2.1",
        "result_class": "marker_fm_transition_silver_formal",
        "formal_benchmark": True,
        "ground_truth": {
            "provider_id": prov, "state_key": STATE_KEY,
            "confidence_mode": "all", "exclude_uncertain_states": False,
        },
        "dataset": {"id": dataset_id, "h5ad_path": meta["input"], "time_key": time_key},
        "output": {"base_dir": out_dir},
    }
    pkey, pval = PARAM_KEY[method]
    cfg[pkey] = pval
    if method == "cellrank2":
        # keep cellrank2 time_key aligned with the scenario time axis
        cfg[pkey] = dict(pval); cfg[pkey]["time_key"] = time_key
    if scenario_params is not None:
        cfg["scenario_params"] = scenario_params
    cfg["lineage"] = {
        "cell_state_key": STATE_KEY, "cell_state_annotation_is_proxy": False,
        "state_system_version": meta["provider"], "reference_graph_path": ref_graph,
        "edge_confidence_mode": "all", "exclude_uncertain_states": False,
    }
    cfg["evaluation"] = {"run_eval_lineage": True,
                         "eval_script": "benchmark/evaluation/eval_lineage.py", "status": "formal"}
    cfg["cell_state_key"] = STATE_KEY
    cfg["time_key"] = time_key
    cfg["state_system"] = {
        "version": meta["provider"], "status": "frozen_silver_standard_milestone_provider",
        "n_states": meta["n_states"], "n_graph_edges": meta["n_edges"],
        "source_h5ad": meta["input"],
        "note": "Frozen marker+foundation-model milestone provider; official metrics use "
                "final_milestone_label_coarse.",
    }
    return run_id, cfg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="benchmark/configs/runtime")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    written = []

    # (1) New datasets A/B/C (projection) + A (OT)
    for ds, meta in NEW_DATASETS.items():
        for method in PROJECTION:
            for scenario in ("A", "B", "C"):
                run_id, cfg = make_config(method, ds, scenario, meta, "abs_day",
                                          meta["scenarios"][scenario], "hvg2000")
                written.append(_dump(args.out, run_id, cfg))
        for method in OT:
            run_id, cfg = make_config(method, ds, "A", meta, "abs_day", None, "hvg2000")
            written.append(_dump(args.out, run_id, cfg))

    # (2) Pseudotime D/E/F (projection) + D (OT) for the two primaries
    for ds, meta in PSEUDOTIME_DATASETS.items():
        for method in PROJECTION:
            for scenario in ("D", "E", "F"):
                run_id, cfg = make_config(method, ds, scenario, meta, PSEUDOTIME_TIME_KEY,
                                          PSEUDOTIME_SCENARIOS[scenario], "pseudotime")
                written.append(_dump(args.out, run_id, cfg))
        for method in OT:
            run_id, cfg = make_config(method, ds, "D", meta, PSEUDOTIME_TIME_KEY,
                                      PSEUDOTIME_SCENARIOS["D"], "pseudotime")
            written.append(_dump(args.out, run_id, cfg))

    print(f"Wrote {len(written)} configs to {args.out}")
    for w in written:
        print("  ", os.path.basename(w))


def _dump(out_dir, run_id, cfg):
    path = os.path.join(out_dir, f"{run_id}.yaml")
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, sort_keys=False, default_flow_style=False)
    return path


if __name__ == "__main__":
    main()
