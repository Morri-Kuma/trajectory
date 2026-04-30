# Diagnosis: identical Jaccard and single-step recovery between WOT and CellRank2 (Scenario A, scGPT v1)

Observed in `benchmark/results/{wot,cellrank2}/scenario_A_scgpt_v1/lineage_metrics.json`:

| Metric               | WOT    | CellRank2 |
|----------------------|--------|-----------|
| AUROC                | 0.8559 | 0.7636    |
| AUPRC                | 0.4564 | 0.3620    |
| `jaccard_similarity` | 0.1582 | 0.1582    |
| `single_step_recovery` | 0.4194 | 0.4194  |
| `multi_step_recovery`  | 0.3654 | 0.3462  |

AUROC, AUPRC, and multi-step differ as expected; Jaccard and single-step are *pixel-identical*. Below is the exact cause of each tie, verified from the STM and edge files on disk.

## 1. Jaccard equality is an ARTIFACT of the evaluator rule

### Where it happens

`benchmark/evaluation/eval_lineage.py::compute_jaccard` operates on the `predicted_edges` set loaded from `lineage_graph_edges.csv` (see `run_lineage_evaluation`, the block that reads `edges_df` and builds `predicted_edges`).

Both adapters emit a **fully dense** edge list:

- `benchmark/methods/WOT/run.py::_aggregate_to_state_level` appends every `(src, tgt)` for which `w > 0`. WOT transport maps produce strictly-positive mass almost everywhere, so every off-diagonal and diagonal cell of the 14×14 state-transition matrix ends up in the edge list.
- `benchmark/adapters/cellrank2_adapter.py::_aggregate_cr2_to_state_level` does the same: `if w > 0.0: edge_rows.append(...)`.

I verified on disk:
- `lineage_graph_edges.csv` has exactly **196 rows** (= 14 × 14) for both methods.
- `wot_edges == cr2_edges` as sets (Jaccard between the two edge sets is 1.0).

Given a fully dense predicted-edge set and a reference set of size `n_ref = 31`:

```
|predicted ∩ reference| = 31         (reference is a subset of the full 14×14)
|predicted ∪ reference| = 196
jaccard                = 31 / 196 = 0.158163…
```

Every method whose STM is numerically dense will hit this exact value. The metric therefore does not discriminate WOT vs CellRank2 — it discriminates *nothing*. This is a real evaluation-rule issue, not a biological property of the two methods.

### Fix applied

I added a second Jaccard field — `jaccard_similarity_topk` — that restricts the predicted edge set to the top-`n_ref` **off-diagonal** entries of the STM (same rule as `single_step_recovery`). This is a method-discriminating similarity and is what the summary/ranking script now uses.

Implementation: `_topk_predicted_edges` in `benchmark/evaluation/eval_lineage.py`, wired into `run_lineage_evaluation` next to `jaccard_similarity`.

After the fix (re-evaluated from the existing CSVs, no method re-run):

| Metric                      | WOT    | CellRank2 |
|-----------------------------|--------|-----------|
| `jaccard_similarity`        | 0.1582 | 0.1582    | ← preserved for traceability; still artifact-bound
| `jaccard_similarity_topk`   | 0.4091 | 0.4091    | ← on top-31 off-diagonal edges (18 of 31 recovered by each)

The top-k Jaccard still ties at 0.4091 here, but for a different and legitimate reason — see §2 below. The *difference* is that top-k Jaccard *can* move with method quality (we verified the top-31 off-diagonal sets have only 25 overlap out of 31 between WOT and CellRank2, so the fact that each recovers 18 is a coincidence, not a forced outcome).

`jaccard_similarity` is kept as-is so historical numbers remain reproducible; downstream ranking uses the top-k version.

## 2. Single-step recovery equality is a COINCIDENCE, not an artifact

### Where the rule lives

`compute_single_step_recovery` in `eval_lineage.py`:

1. Sets `n_ref = number of 1s in the reference matrix` (= 31 under `medium_and_above`).
2. Computes the threshold = the `n_ref`-th largest value in the flattened predicted STM.
3. `pred_edges = {(i,j) : predicted[i,j] >= threshold}`.
4. Recovery = `|ref_edges ∩ pred_edges| / n_ref`.

For both methods the threshold selects **exactly 31 entries** (no tie inflation at the cut point). I verified:

- WOT top-31 threshold = 0.175, 31 entries at or above.
- CR2 top-31 threshold = 0.108, 31 entries at or above.
- The two top-31 sets are NOT identical. There is 25-edge overlap, with 6 method-unique edges on each side.
- Both methods happen to land 13 of their 31 top picks inside the reference edge set.

So the `13 / 31 = 0.4194` equality is a genuine tied count from differently-chosen edge sets. Single-step recovery is behaving as designed; the rule *can* differentiate methods, and the equality here is arithmetic coincidence at this specific reference size.

A note on the diagonal: `compute_single_step_recovery` does NOT mask the STM diagonal, so some of the top-31 picks can be self-loops (which the reference graph explicitly excludes). The top-k Jaccard I added does mask the diagonal and pairs of 18-recovered-edges-out-of-31 fall out naturally; that is the cleaner comparison point for ranking.

### Should single_step_recovery be changed?

Recommendation: **document only, do not change now.**

Reasons:
- The equality here is a coincidence at `n_ref = 31`; it will not generally hold across scenarios.
- Changing the rule would silently re-score every historical result. Nothing forces us to do that today.
- Adding `jaccard_similarity_topk` already gives the summary table a sparse, direction-aware, diagonal-free comparison metric if single-step coincidentally ties.

If later we want single-step to exclude the diagonal too, it should be done as a deliberate, versioned v2.1 metric change with a re-run of all scenarios — not quietly.

## 3. Multi-step recovery is already discriminating

`compute_multi_step_recovery` walks the reference graph with `matrix_power(ref > 0, 2)` and then applies the same top-k rule against that multi-hop reference. WOT (0.365) beats CellRank2 (0.346) here — the two methods differ by one recovered multi-hop path, so this metric already does its job.

## 4. Summary of recommended actions (updated 2026-04-20)

1. **Keep `jaccard_similarity` unchanged** for traceability, and keep it as the
   framework-defined official Jaccard metric. It is artifact-prone for dense
   STMs — the fix belongs at the adapter layer, not in the evaluator (see point 4).
2. **`jaccard_similarity_topk` is DIAGNOSTIC ONLY, not part of the official rank.**
   It is shown in `lineage_summary.csv` as a supplementary column so readers
   can see the top-k comparison, but `benchmark/reports/summarize_lineage.py`
   only ranks on the framework metric set
   `{AUROC, AUPRC, jaccard_similarity, single_step_recovery, multi_step_recovery}`.
   An earlier revision of this note recommended using topk for ranking; that
   recommendation is superseded — silently substituting a metric would have
   changed the benchmark definition without a versioned framework revision.
3. **Do not modify `compute_single_step_recovery` now.** The tie is a
   coincidence; rule-change is a v2.1 task.
4. Longer-term: require adapters to emit a binarised `lineage_graph_edges.csv`
   (top-k or threshold-based) in addition to the dense STM, so `compute_jaccard`
   receives a sparse edge list by construction. This removes the artifact at
   the source. Not done today because it would require changes in both adapters
   and an explicit policy choice for the binarisation threshold. When that
   change lands, `jaccard_similarity` will become method-discriminating on its
   own and `jaccard_similarity_topk` can be removed.

## 5. Exact code references

- `benchmark/evaluation/eval_lineage.py::compute_jaccard` — set-based Jaccard on whatever `lineage_graph_edges.csv` contains.
- `benchmark/evaluation/eval_lineage.py::compute_single_step_recovery` — top-`n_ref` rule on flattened STM (diagonal included).
- `benchmark/evaluation/eval_lineage.py::_topk_predicted_edges` — new helper; top-`n_ref` off-diagonal edges.
- `benchmark/methods/WOT/run.py::_aggregate_to_state_level` — appends every `w > 0` to the edge list (dense emit).
- `benchmark/adapters/cellrank2_adapter.py::_aggregate_cr2_to_state_level` — same (dense emit).
