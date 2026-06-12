#!/usr/bin/env python
"""Minimal local smoke test for the v2.1 extensions (code-path only).

Exercises the NEW code paths on tiny synthetic data, with no GPU/torch and no
real benchmark data, to prove the process runs end-to-end before full-data runs
on Shirokane:

  1. registry      : GSE298212Dataset / GSE218855Dataset are registered.
  2. pseudotime    : build_pseudotime_axis() produces DPT bins (Scenarios D-F).
  3. D-split       : the bin axis splits into train/heldout (dispatcher semantics).
  4. factory A/B/C : factory.build_dataset() builds + splits a new dataset (Scenario B).
  5. configs       : generated runtime configs parse and carry the right wiring.

Run from the repo root:
    PYTHONPATH=. python scripts/smoke_test_v21_extensions.py
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import anndata as ad
import yaml

PASS, FAIL = "[PASS]", "[FAIL]"
failures = []


def check(name, cond, detail=""):
    print(f"{PASS if cond else FAIL} {name}" + (f" — {detail}" if detail else ""))
    if not cond:
        failures.append(name)


def synthetic(n_per_t, times, n_genes=40, seed=0):
    rng = np.random.default_rng(seed)
    blocks, days, states = [], [], []
    milestones = ["hADSCs", "epithelial_like", "intermediate_plastic", "hCiPS"]
    for i, t in enumerate(times):
        center = rng.normal(i, 0.2, size=n_genes)
        blocks.append(rng.normal(center, 1.0, size=(n_per_t, n_genes)))
        days += [float(t)] * n_per_t
        states += [milestones[min(i, len(milestones) - 1)]] * n_per_t
    X = np.maximum(np.vstack(blocks), 0).astype("float32")
    obs = {"abs_day": np.array(days, dtype=float),
           "final_milestone_label_coarse": np.array(states, dtype=object)}
    a = ad.AnnData(X=X)
    a.obs["abs_day"] = obs["abs_day"]
    a.obs["final_milestone_label_coarse"] = obs["final_milestone_label_coarse"]
    return a


def main():
    tmp = Path(tempfile.mkdtemp(prefix="v21_smoke_"))

    # 1) registry
    from benchmark.shared.dataset.base import DATASET_REGISTRY
    import benchmark.shared.dataset  # noqa: F401  triggers registry imports
    check("registry: GSE298212Dataset registered", "GSE298212Dataset" in DATASET_REGISTRY)
    check("registry: GSE218855Dataset registered", "GSE218855Dataset" in DATASET_REGISTRY)

    # 2) pseudotime axis builder (Scenarios D-F precompute)
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from build_pseudotime_axis import build as build_pt
    a = synthetic(40, [0, 1, 3, 6, 8])
    in_h5ad = tmp / "obs_input.h5ad"
    a.write_h5ad(in_h5ad)
    pt_h5ad = tmp / "pseudotime_input.h5ad"
    build_pt(str(in_h5ad), str(pt_h5ad), time_key="abs_day", num_bins=15,
             n_top_genes=30, n_neighbors=10)
    apt = ad.read_h5ad(pt_h5ad)
    has_bins = "dpt_pseudotime_bin_numeric" in apt.obs.columns
    n_bins = int(apt.obs["dpt_pseudotime_bin_numeric"].nunique()) if has_bins else 0
    check("pseudotime: bin column present", has_bins)
    check("pseudotime: 15 bins produced", n_bins == 15, f"n_bins={n_bins}")

    # 3) D-scenario split on the bin axis (dispatcher semantics)
    from benchmark.shared.dataset.preprocessors.scenario_timepoint_split import (
        split_adata_by_timepoints,
    )
    train, test = split_adata_by_timepoints(
        apt, time_key="dpt_pseudotime_bin_numeric",
        train_times=[0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14], heldout_times=[3, 7, 11],
    )
    check("D-split: non-empty train", train.n_obs > 0, f"train={train.n_obs}")
    check("D-split: non-empty heldout", test.n_obs > 0, f"test={test.n_obs}")

    # 4) factory build + Scenario B split for a NEW dataset (GSE298212)
    from benchmark.shared.dataset.factory import build_dataset
    b = synthetic(30, [0, 1, 3, 6, 8], seed=1)
    abc_h5ad = tmp / "gse298212_like.h5ad"
    b.write_h5ad(abc_h5ad)
    dataset = build_dataset({
        "id": "GSE298212",
        "h5ad_path": str(abc_h5ad),
        "time_key": "abs_day",
        "cell_state_key": "final_milestone_label_coarse",
        "scenario_params": {"train_times": [0.0, 1.0, 3.0], "heldout_times": [6.0, 8.0]},
    }, output_dir=str(tmp / "out"))
    result = dataset.load_data()
    is_split = isinstance(result, tuple) and len(result) == 2
    check("factory: GSE298212 builds + splits (Scenario B)", is_split)
    if is_split:
        tr, te = result
        check("factory: train holds only train_times", set(np.unique(
            tr.obs["abs_day"].astype(float))).issubset({0.0, 1.0, 3.0}),
            f"train days={sorted(set(tr.obs['abs_day'].astype(float)))}")

    # 5) generated configs parse + correct wiring
    cfg_dir = Path("benchmark/configs/runtime")
    samples = {
        "scnode_gse298212_marker_fm_silver_B_hvg2000_formal.yaml":
            ("GSE298212", "B", "abs_day", "gse298212_marker_fm_transition_silver_v1"),
        "mioflow_gse230659_marker_fm_silver_D_pseudotime_formal.yaml":
            ("GSE230659", "D", "dpt_pseudotime_bin_numeric", "gse230659_marker_fm_transition_silver_v1"),
        "wot_gse218855_marker_fm_silver_A_hvg2000_formal.yaml":
            ("GSE218855", "A", "abs_day", "gse218855_marker_fm_transition_silver_v1"),
    }
    for fname, (ds, scen, tkey, prov) in samples.items():
        p = cfg_dir / fname
        if not p.exists():
            check(f"config: {fname} exists", False)
            continue
        c = yaml.safe_load(p.read_text(encoding="utf-8-sig"))
        ok = (c["dataset"]["id"] == ds and c["scenario"] == scen
              and c["time_key"] == tkey
              and c["ground_truth"]["provider_id"] == prov
              and c["lineage"]["reference_graph_path"].endswith(f"{prov}/reference_graph.json"))
        check(f"config: {fname} wired correctly", ok)
        if scen in ("B", "C", "D", "E", "F"):
            check(f"config: {fname} has scenario_params", "scenario_params" in c)

    print()
    if failures:
        print(f"SMOKE TEST FAILED: {len(failures)} check(s) failed: {failures}")
        sys.exit(1)
    print("SMOKE TEST PASSED: all v2.1 extension code paths run.")


if __name__ == "__main__":
    main()
