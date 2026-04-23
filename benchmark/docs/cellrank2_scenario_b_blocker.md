# CellRank2 Scenario B — blocker report

Date: 2026-04-20
Author: cleanup pass (session continuation)

## Status

The adapter-side fix is complete. CellRank2 Scenario B cannot be launched in this sandbox because CellRank2 is not installed here. Every other precondition is met.

## What is now ready

- `benchmark/adapters/cellrank2_adapter.py::CellRank2Adapter._run_lineage_fidelity_impl` has a new Step 0 block that filters `self.adata` to `scenario_params.train_times` when provided. The filter mirrors the one at the top of `benchmark/methods/WOT/run.py` so Scenario B (early-only training) is honored consistently across both methods.
- `benchmark/evaluation/eval_dispatch.py` now injects the method config's `scenario_params` into `scenario_config` so the adapter actually receives `train_times`. Without this, the filter block above would never fire from the CLI path.
- `benchmark/configs/cellrank2_gse230659_observed_scgpt_v1_scenarioB.yaml` already has the correct `scenario_params.train_times` list and is unchanged from the previous session (only the header comment was updated).
- The source AnnData is on disk at `benchmark/results/scgpt/full/adata_scgpt_annotated.h5ad`.
- `wot` is importable in this sandbox.

## What is blocking an actual run

- `cellrank` is not installed in this sandbox. `python -c "import cellrank"` fails with `ModuleNotFoundError`. The previous Scenario A CellRank2 run was done in a different environment.
- Even with CellRank installed, the run is expected to be long (the Scenario A CellRank2 run was ~hours); a subsample pilot (`--subsample 500` equivalent for CR2) would be far quicker but would have to be launched via a wrapper or by extending the adapter to accept a subsample flag.

## How to unblock

Run in the environment that already has CellRank + WOT:

```
pip install cellrank  # only if not already installed
python benchmark/evaluation/eval_dispatch.py \
    --method cellrank2 \
    --scenario B \
    --adata benchmark/results/scgpt/full/adata_scgpt_annotated.h5ad \
    --output-dir benchmark/results/cellrank2/scenario_B_scgpt_v1 \
    --method-config benchmark/configs/cellrank2_gse230659_observed_scgpt_v1_scenarioB.yaml
```

Then, if the run is launched with a subsample for speed, explicitly set `result_class: pilot` in the produced `run_metadata.json` so the summary script keeps it out of the official table (see `benchmark/reports/summarize_lineage.py::_infer_result_class`).

## Correctness verification (without an end-to-end run)

The adapter-side filter is small enough to be reviewable by inspection:

- Reads `self.scenario_config["scenario_params"]["train_times"]` (or falls back to `scenario_split`).
- Raises if `adata.obs[time_key]` is missing.
- Uses `isin(train_times_f)` on a float-cast of the time column, identical to `benchmark/methods/WOT/run.py` lines 427–441.
- Reassigns `self.adata` so the rest of the adapter (WOT sub-fit + `_aggregate_cr2_to_state_level`) sees the filtered cell universe.

No test infrastructure exists in the repo for the adapter path, so a unit test was not added — that would be out of scope for this cleanup pass (the working rules explicitly say "do not redesign"). A regression check should be added in a later session once the first real Scenario B CR2 run is available.
