"""
recompute_scenario_b_baseline.py
--------------------------------
Scenario-aware recomputation of the correlation baseline for Scenario B
(observed-time extrapolation, scGPT v1 state system).

Why this exists
---------------
The previous `scripts/reeval_lineage_with_baseline.py` script loads the *full*
adata and passes it to `run_lineage_evaluation` for every run directory it is
given. When that was applied to Scenario B, the baseline was silently computed
on all 75k cells instead of only the Scenario B training split
(abs_day in [0.5, 2, 4, 8, 12, 16] with --subsample 500 per time point, i.e.
3000 cells). That made the Scenario B baseline non-scenario-specific and
byte-identical to the Scenario A baseline.

This script reproduces the exact Scenario B training split used by
benchmark/methods/WOT/run.py, then calls compute_correlation_baseline() against
that split. It only overwrites:
  - baseline_state_transition_matrix.csv
  - baseline_lineage_graph_edges.csv
  - the "baseline" block of lineage_metrics.json
Method-level metrics (AUROC/AUPRC/Jaccard/single-step/multi-step) and all other
fields are preserved exactly.

Usage
-----
    python scripts/recompute_scenario_b_baseline.py

By default it targets benchmark/results/wot/scenario_B_scgpt_v1 and uses the
Scenario B training-split rules declared in
benchmark/configs/wot_gse230659_observed_scgpt_v1_scenarioB.yaml.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def _project_root() -> Path:
    here = Path(__file__).resolve().parent
    for c in [here, *here.parents]:
        if (c / "benchmark").exists():
            return c
    return here.parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--adata",
        default="benchmark/results/scgpt/full/adata_scgpt_annotated.h5ad",
    )
    parser.add_argument(
        "--run-dir",
        default="benchmark/results/wot/scenario_B_scgpt_v1",
    )
    parser.add_argument(
        "--reference-graph",
        default="benchmark/datasets/scgpt_reference_graph_v1.json",
    )
    parser.add_argument(
        "--edge-confidence-mode", default="medium_and_above",
    )
    parser.add_argument(
        "--cell-state-key",
        default="scgpt_pseudostate_provisional",
    )
    parser.add_argument(
        "--time-key", default="abs_day",
    )
    parser.add_argument(
        "--train-times",
        nargs="+", type=float,
        default=[0.5, 2.0, 4.0, 8.0, 12.0, 16.0],
        help="Scenario B train time points. Default matches the Scenario B config.",
    )
    parser.add_argument(
        "--subsample", type=int, default=500,
        help=(
            "Head-of-group subsample per time point (default: 500, matches the "
            "Scenario B pilot). Set to 0 to disable."
        ),
    )
    args = parser.parse_args()

    root = _project_root()
    sys.path.insert(0, str(root))

    import anndata as ad  # noqa: E402
    import numpy as np  # noqa: E402
    from benchmark.evaluation.eval_lineage import (  # noqa: E402
        compute_correlation_baseline,
        load_reference_graph,
        compute_jaccard,
        compute_single_step_recovery,
        compute_multi_step_recovery,
    )

    adata_path = (root / args.adata) if not Path(args.adata).is_absolute() else Path(args.adata)
    run_dir = (root / args.run_dir) if not Path(args.run_dir).is_absolute() else Path(args.run_dir)
    ref_path = (root / args.reference_graph) if not Path(args.reference_graph).is_absolute() else Path(args.reference_graph)

    print(f"[recompute_B_baseline] adata       : {adata_path}")
    print(f"[recompute_B_baseline] run_dir     : {run_dir}")
    print(f"[recompute_B_baseline] reference   : {ref_path}")
    print(f"[recompute_B_baseline] train_times : {args.train_times}")
    print(f"[recompute_B_baseline] subsample   : {args.subsample}")

    # ------------------------------------------------------------------
    # Load adata (backed='r' so we don't materialise 75k × 23k in RAM).
    # ------------------------------------------------------------------
    print(f"[recompute_B_baseline] Loading adata (backed='r') ...")
    adata_full = ad.read_h5ad(adata_path, backed="r")
    print(f"[recompute_B_baseline] Full adata: {adata_full.n_obs} × {adata_full.n_vars}")

    # ------------------------------------------------------------------
    # Reproduce the Scenario B training split — identical rule to
    # benchmark/methods/WOT/run.py.
    #   1. Optional: head-of-group subsample per time point.
    #   2. scenario_params.train_times filter on abs_day.
    # Note that benchmark/methods/WOT/run.py applies subsample FIRST and then
    # the train_times filter. We follow that exact order.
    # ------------------------------------------------------------------
    obs = adata_full.obs
    keep_mask = np.zeros(adata_full.n_obs, dtype=bool)

    if args.subsample and args.subsample > 0:
        # Per-group head-of-index selection. Match the .index[:N] rule in
        # benchmark/methods/WOT/run.py lines 406–414.
        for t, grp_idx in obs.groupby(args.time_key, observed=True).groups.items():
            # grp_idx is an Index of obs labels; map back to positional.
            head = grp_idx[: args.subsample]
            positional = obs.index.get_indexer(head)
            keep_mask[positional] = True
    else:
        keep_mask[:] = True

    # Apply keep_mask, then enforce train_times.
    times_arr = obs[args.time_key].astype(float).values
    time_mask = np.isin(times_arr, args.train_times)
    final_mask = keep_mask & time_mask

    # ------------------------------------------------------------------
    # Materialise the Scenario B training split into memory (3000 cells).
    # Backed slicing → ad.to_memory() is cheap at this size.
    # ------------------------------------------------------------------
    train_adata = adata_full[final_mask].to_memory()
    print(
        f"[recompute_B_baseline] Scenario B training split: "
        f"{train_adata.n_obs} cells "
        f"(expected ~3000 for subsample=500 × 6 early time points)"
    )
    print(
        f"[recompute_B_baseline] Time points in split: "
        f"{sorted(train_adata.obs[args.time_key].astype(float).unique())}"
    )
    print(
        f"[recompute_B_baseline] States present: "
        f"{sorted(train_adata.obs[args.cell_state_key].astype(str).unique())}"
    )

    # ------------------------------------------------------------------
    # Load reference graph and recompute baseline.
    # ------------------------------------------------------------------
    reference_matrix, reference_edges, _ = load_reference_graph(
        str(ref_path),
        edge_confidence_mode=args.edge_confidence_mode,
    )
    print(
        f"[recompute_B_baseline] Reference graph: "
        f"{reference_matrix.shape[0]} nodes, {len(reference_edges)} edges "
        f"(mode={args.edge_confidence_mode!r})"
    )

    baseline = compute_correlation_baseline(
        adata=train_adata,
        reference_matrix=reference_matrix,
        reference_edges=reference_edges,
        cell_state_key=args.cell_state_key,
        time_key=args.time_key,
        output_dir=str(run_dir),     # overwrites baseline_*.csv here
    )

    # ------------------------------------------------------------------
    # Rewrite ONLY the baseline block of lineage_metrics.json.
    # Leave the method-level metrics untouched.
    # ------------------------------------------------------------------
    metrics_path = run_dir / "lineage_metrics.json"
    if metrics_path.exists():
        with open(metrics_path, encoding="utf-8") as f:
            lm = json.load(f)
    else:
        print(f"[recompute_B_baseline] WARN: {metrics_path} missing — writing a new one")
        lm = {}

    baseline_with_provenance = dict(baseline)
    baseline_with_provenance["scenario_specific"] = True
    baseline_with_provenance["train_times"] = list(args.train_times)
    baseline_with_provenance["subsample_per_timepoint"] = args.subsample or None
    baseline_with_provenance["n_train_cells"] = int(train_adata.n_obs)
    baseline_with_provenance["recomputed_by"] = "scripts/recompute_scenario_b_baseline.py"

    lm["baseline"] = baseline_with_provenance

    with open(metrics_path, "w", encoding="utf-8") as f:
        json.dump(lm, f, indent=2)
    print(f"[recompute_B_baseline] Updated baseline block → {metrics_path}")

    print("[recompute_B_baseline] Summary of new baseline values:")
    for k in ("auroc", "auprc", "jaccard_similarity",
              "single_step_recovery", "multi_step_recovery",
              "n_states_used", "n_train_cells"):
        print(f"  {k:28s} = {baseline_with_provenance.get(k)}")


if __name__ == "__main__":
    main()
