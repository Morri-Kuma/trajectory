"""Annotate scNODE reduced-validation run metadata."""

from __future__ import annotations

import json
from pathlib import Path


RUNS = {
    "benchmark/results/scnode/scenario_A_scgpt_v1_full_hvg2000_reduced": {
        "scenario_variant": "A_full_hvg2000_reduced",
        "training_times": "all",
    },
    "benchmark/results/scnode/scenario_B_scgpt_v1_hvg2000_reduced": {
        "scenario_variant": "B_hvg2000_reduced",
    },
    "benchmark/results/scnode/scenario_C_scgpt_v1_hvg2000_reduced": {
        "scenario_variant": "C_hvg2000_reduced",
    },
}


def main() -> None:
    for run_dir_str, extra in RUNS.items():
        run_dir = Path(run_dir_str)
        meta_path = run_dir / "run_metadata.json"
        if not meta_path.exists():
            print(f"missing {meta_path}")
            continue
        with open(meta_path, encoding="utf-8") as f:
            meta = json.load(f)
        meta["result_class"] = "reduced_validation"
        meta["feature_space"] = {
            "type": "HVG",
            "n_genes": 2000,
            "source": "benchmark/results/scnode/hvg_inputs/adata_scnode_full_A_hvg2000.h5ad",
        }
        meta["training_setting"] = {
            "label": "hvg2000_reduced",
            "pretrain_iters": 50,
            "epochs": 3,
            "iters": 20,
            "batch_size": 32,
            "n_sim_cells": meta.get("n_sim_cells"),
        }
        meta["scenario_variant"] = extra.get("scenario_variant")
        with open(meta_path, "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)
        print(f"updated {meta_path}")


if __name__ == "__main__":
    main()
