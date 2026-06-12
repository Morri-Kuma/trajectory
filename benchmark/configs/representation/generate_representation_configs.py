"""
generate_representation_configs.py
==================================

Generate representation-aware method configs (Work-Plan Step 3) for the four
representation arms, derived from the existing *formal* expression-space configs
so that scenario splits (train_times / heldout_times), lineage wiring, and model
hyperparameters stay identical. Only the representation block, dataset path,
result class, and evaluation wiring change.

Run from the repo root:

    python -m benchmark.configs.representation.generate_representation_configs

This writes one YAML per (method, representation, scenario) into this directory.

Each emitted config:
  * sets representation.enabled = true, input_mode = obsm, obsm_key = X_rep,
    final_dim = 50, reducer_fit_scope = train_only;
  * points dataset.h5ad_path at the scenario-specific representation input h5ad
    (built by build_representation_inputs.py with train-only fitting);
  * marks itself supplementary (formal_benchmark: false) and routes evaluation
    to the representation evaluators, NOT the gene-expression eval_dispatch
    forecast (Risk 1: representation forecast != gene-expression forecast).
"""

from __future__ import annotations

import copy
from pathlib import Path

import yaml

METHODS = ("scnode", "mioflow", "prescient")
SCENARIOS = ("A", "B", "C")
REPRESENTATIONS = (
    "rep_hvg_pca50",
    "rep_geneformer_cls_pca50",
    "rep_scgpt_cls_pca50",
    "rep_scfoundation_pca50",
)
DATASET_TAG = "gse230659"
FINAL_DIM = 50


def _repo_root() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "benchmark").is_dir():
            return parent
    raise RuntimeError("Cannot locate repo root.")


def _formal_config_path(root: Path, method: str, scenario: str) -> Path:
    return (
        root / "benchmark" / "configs" / "runtime"
        / f"{method}_{DATASET_TAG}_marker_fm_silver_{scenario}_hvg2000_formal.yaml"
    )


def _rep_input_h5ad(rep: str, scenario: str) -> str:
    return (
        f"benchmark/inputs/representation/{DATASET_TAG}/{rep}/{scenario}/"
        f"{DATASET_TAG.upper()}_{rep}_{scenario}_X_rep.h5ad"
    )


def build_config(base: dict, method: str, rep: str, scenario: str) -> dict:
    cfg = copy.deepcopy(base)
    run_id = f"{method}_{DATASET_TAG}_{rep}_{scenario}"
    cfg["run_id"] = run_id
    cfg["scenario_name"] = (
        f"{DATASET_TAG.upper()} {method} scenario {scenario} "
        f"representation-dynamics {rep} (supplementary)"
    )

    # Point at the scenario-specific representation input (X_rep in obsm + X).
    cfg.setdefault("dataset", {})["h5ad_path"] = _rep_input_h5ad(rep, scenario)
    if "state_system" in cfg:
        cfg["state_system"]["source_h5ad"] = _rep_input_h5ad(rep, scenario)

    cfg["output"] = {
        "base_dir": (
            f"benchmark/results/representation_dynamics/{method}/"
            f"{DATASET_TAG}_{rep}_{scenario}"
        )
    }

    # Representation block (Work-Plan Step 3).
    cfg["representation"] = {
        "enabled": True,
        "input_mode": "obsm",
        "obsm_key": "X_rep",
        "representation_id": rep,
        "input_space": ("expression_pca" if rep == "rep_hvg_pca50"
                        else "scfm_embedding"),
        "output_space": "representation",
        "final_dim": FINAL_DIM,
        "reducer_fit_scope": "train_only",
    }

    # Supplementary classification; keep official labels independent (Risk 4).
    cfg["result_class"] = "representation_dynamics_supplementary"
    cfg["formal_benchmark"] = False

    # Route evaluation to the representation evaluators. The gene-expression
    # eval_dispatch forecast is NOT applicable in representation space.
    cfg["evaluation"] = {
        "mode": "representation",
        "run_eval_lineage": False,
        "forecast_eval_script":
            "benchmark/evaluation/eval_representation_forecast.py",
        "temporal_signal_eval_script":
            "benchmark/evaluation/eval_representation_temporal_signal.py",
        "status": "supplementary",
        "note": (
            "Representation-space forecast (Mode A). Outputs are representation "
            "forecasts, not gene-expression forecasts. Score with the "
            "representation evaluators; do not feed projected outputs to the "
            "official gene-expression eval_forecast."
        ),
    }
    return cfg


def main() -> int:
    root = _repo_root()
    out_dir = root / "benchmark" / "configs" / "representation"
    out_dir.mkdir(parents=True, exist_ok=True)

    written = []
    for method in METHODS:
        for scenario in SCENARIOS:
            base_path = _formal_config_path(root, method, scenario)
            if not base_path.exists():
                print(f"[skip] base config missing: {base_path}")
                continue
            with open(base_path, encoding="utf-8-sig") as f:
                base = yaml.safe_load(f)
            for rep in REPRESENTATIONS:
                cfg = build_config(base, method, rep, scenario)
                out_path = out_dir / f"{method}_{DATASET_TAG}_{rep}_{scenario}.yaml"
                with open(out_path, "w", encoding="utf-8") as f:
                    yaml.safe_dump(cfg, f, sort_keys=False, default_flow_style=False)
                written.append(out_path.name)
    print(f"[generate_representation_configs] wrote {len(written)} configs to {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
