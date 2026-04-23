# Technical note: is the Scenario B baseline actually scenario-specific?

Date: 2026-04-20
Author: cleanup pass (session continuation)

## Question

The correlation baseline embedded in `benchmark/results/wot/scenario_B_scgpt_v1/lineage_metrics.json` currently reports:

```
auroc              = 0.6742913000977517
auprc              = 0.2339979522788161
jaccard_similarity = 0.26
single_step        = 0.2903225806451613
multi_step         = 0.2884615384615384
n_states_used      = 14
```

These are the **same** numbers as in the Scenario A baseline, byte-for-byte. Is the Scenario B baseline actually scenario-specific, or is it a leak from Scenario A / full-data?

## Answer (short version)

**It is a leak. The Scenario B baseline currently on disk was computed from the full adata (~75k cells across all 15 time points), not from the Scenario B training split (3000 cells restricted to `abs_day ∈ [0.5, 2, 4, 8, 12, 16]`).**

## Evidence

Byte-identical across all three runs:

```
$ diff -q benchmark/results/wot/scenario_A_scgpt_v1/baseline_state_transition_matrix.csv \
          benchmark/results/wot/scenario_B_scgpt_v1/baseline_state_transition_matrix.csv
# silent: files identical
$ diff -q benchmark/results/wot/scenario_A_scgpt_v1/baseline_state_transition_matrix.csv \
          benchmark/results/cellrank2/scenario_A_scgpt_v1/baseline_state_transition_matrix.csv
# silent: files identical
```

The same is true for `baseline_lineage_graph_edges.csv`.

## Code path that produced the leak

The previous session used `scripts/reeval_lineage_with_baseline.py` to populate the baseline block in all existing `lineage_metrics.json` files in place (because the original runs predated the real-baseline implementation). That script:

1. Loads the full adata once (`anndata.read_h5ad(adata_path, backed="r")`).
2. Passes that same full adata into `run_lineage_evaluation(...)` for every run directory listed on the command line — including Scenario B.
3. `compute_correlation_baseline` then iterates over every cell in `adata` to compute per-state mean expression, so for Scenario B it used **all 75k cells**, not the 3000-cell Scenario B training split.

Relevant lines:
- `scripts/reeval_lineage_with_baseline.py` line ~64: `ad.read_h5ad(adata_path, backed="r")` — full adata.
- `benchmark/evaluation/eval_lineage.py::compute_correlation_baseline` — no scenario filter; iterates whatever `adata` is passed in.

By contrast, the fresh-run path in `benchmark/methods/WOT/run.py` at line 502–511 correctly passes `adata=train_adata`, which **has** been filtered by the `scenario_params.train_times` block at lines 427–441. So a brand-new Scenario B run would produce a scenario-specific baseline. The problem is purely that the re-eval script overwrote the correct baseline with a full-adata baseline.

## Why the difference matters

The correlation baseline is meant to answer *"could a naïve non-lineage signal produce these metrics on this run's data?"* If Scenario B is an early-only extrapolation scenario, late-stage states (which only appear in held-out time points) may have fewer or zero cells in the training split, and the per-state mean expressions — and therefore the baseline — should degrade accordingly. Using a full-adata baseline for a Scenario B evaluation systematically **overstates** the baseline for B (it gets all 14 states populated even if most of them are absent from B's training data), which in turn makes WOT Scenario B look worse vs baseline than it really is.

## Fix strategy (what the cleanup session does)

1. Add a scenario-aware re-eval path so the baseline for a given run is computed only on the cells that run actually trained on.
2. For Scenario B (scGPT v1) specifically: reproduce the WOT Scenario B training split (`train_times = [0.5, 2, 4, 8, 12, 16]`, `--subsample 500` head-of-group), run the baseline against that, and overwrite only `baseline_*.csv` and the `baseline` block of `lineage_metrics.json`. The method-level metrics (AUROC/AUPRC/Jaccard/…) stay untouched.
3. After the fix, `baseline_state_transition_matrix.csv` for Scenario B is expected to differ from Scenario A's.

## Why Scenario A's baseline stays the same

Scenario A uses the full adata for training (no `train_times` filter), so its baseline *is* already the full-adata baseline. The Scenario A files on disk are correct by coincidence — they happen to match what a properly scenario-aware computation would produce.

## What a scenario-specific B baseline is expected to look like

Qualitatively:
- Fewer distinct states will have non-zero mean-expression vectors (late-stage states are underrepresented or absent in early-only training data).
- Pearson correlation over the remaining states still gives a dense correlation matrix but with a smaller effective state set; `n_states_used` in the baseline block should drop from 14 to whatever is present in the Scenario B training split.
- AUROC/AUPRC/Jaccard will be different and likely lower than A's (less signal in the states that are present, and fewer positives in the intersection with the reference graph if some reference edges now involve absent states).

Exact numbers will be reported in the Scenario B results directory once the recompute script runs.

## Related files

- `benchmark/methods/WOT/run.py` — the correct path (passes filtered `train_adata` into the evaluator).
- `scripts/reeval_lineage_with_baseline.py` — the leaky path that overwrote Scenario B's baseline with a full-adata baseline.
- `scripts/recompute_scenario_b_baseline.py` — scenario-aware re-computation script created as part of this cleanup pass.
- `benchmark/results/wot/scenario_B_scgpt_v1/baseline_state_transition_matrix.csv` — currently leaking full-adata; will be overwritten by the recompute.
