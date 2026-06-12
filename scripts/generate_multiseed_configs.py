#!/usr/bin/env python
"""Generate multi-seed runtime configs for per-method confidence intervals.

For each projection-capable method (scNODE, PRESCIENT, MIOFlow) x primary chemical
dataset (GSE178325, GSE230659) x observed-time scenario (A/B/C), this clones the
existing *_hvg2000_formal.yaml config once per seed, overriding the method's
`<method>_params.seed` and redirecting the output to a per-seed directory
(`..._seed<NN>`) so seed runs never overwrite each other or the canonical run.

Outputs configs to benchmark/configs/runtime/multiseed/ and prints the count.
Run on Shirokane (or locally) before submitting run_multiseed_benchmark_array.sh.
"""
from __future__ import annotations

import copy
import os
from pathlib import Path

import yaml

ROOT = Path(os.environ.get("TRAJ_PROJECT_ROOT", os.getcwd()))
RT = ROOT / "benchmark" / "configs" / "runtime"
OUT = RT / "multiseed"
METHODS = ["scnode", "prescient", "mioflow"]
DATASETS = ["gse178325", "gse230659"]
SCEN = ["A", "B", "C"]
SEEDS = [101, 202, 303, 404, 505]  # 5 seeds -> n=5 per (method,dataset,scenario)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    n = 0
    missing = []
    for method in METHODS:
        pkey = f"{method}_params"
        for ds in DATASETS:
            for sc in SCEN:
                base = RT / f"{method}_{ds}_marker_fm_silver_{sc}_hvg2000_formal.yaml"
                if not base.exists():
                    missing.append(str(base))
                    continue
                cfg0 = yaml.safe_load(base.read_text(encoding="utf-8-sig")) or {}
                for seed in SEEDS:
                    cfg = copy.deepcopy(cfg0)
                    if pkey in cfg and isinstance(cfg[pkey], dict):
                        cfg[pkey]["seed"] = seed
                    else:
                        cfg.setdefault(pkey, {})["seed"] = seed
                    run_id = f"{method}_{ds}_marker_fm_silver_{sc}_hvg2000_seed{seed}"
                    cfg["run_id"] = run_id
                    outdir = f"benchmark/results/{method}/{ds}_marker_fm_silver_{sc}_hvg2000_seed{seed}"
                    cfg.setdefault("output", {})["base_dir"] = outdir
                    dst = OUT / f"{run_id}.yaml"
                    dst.write_text(yaml.safe_dump(cfg, sort_keys=False), encoding="utf-8")
                    n += 1
    print(f"[generate_multiseed_configs] wrote {n} configs to {OUT}")
    print(f"[generate_multiseed_configs] {len(METHODS)*len(DATASETS)*len(SCEN)} cells x {len(SEEDS)} seeds")
    if missing:
        print("[generate_multiseed_configs] WARNING missing base configs:")
        for m in missing:
            print("   ", m)


if __name__ == "__main__":
    main()
