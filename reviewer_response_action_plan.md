# Reviewer Response — Prioritized Action Plan

**Project:** scTimeBench conditions benchmark · **Date:** 2026-06-12 · **Status:** pre-resubmission triage

Every reviewer claim below was checked against the actual workspace. Each item carries a
**verdict** (Confirmed / Partly / Inaccurate) with the file evidence, then the concrete
tasks needed to close it. Tasks are tagged `[code]` / `[manuscript]` / `[docs]` and sized
**S** (<½ day), **M** (½–2 days), **L** (needs a server/seed run).

The plan is ordered by the cross-review synthesis, not by reviewer number — the synthesis
is what an editor will read against.

---

## Audit verdicts at a glance

| # | Reviewer claim | Verdict | Evidence in repo |
|---|----------------|---------|------------------|
| R1 | OOD gate is post-hoc, not a pipeline gate | **Confirmed** | `src/pipeline.py:180` emits `frac_unknown_ood` from scANVI `scanvi_unknown_flag` (posterior conf.); geometric OOD lives in `scripts/score_ood.py` (k=15, q=0.95) + `src/annotation/ood.py`. Outputs only in `results/annotation_branch/ood_summary_scanvi_latent.json`, `results/test_outputs/ood_recalibration/`. |
| R1 | Reference graphs have only 3–4 edges | **Confirmed** | edges: 178325=4, 230659=3, 218855=3, 298212=3, 175634=5, 242424(OSKM)=9 (`reference_graph_edges.csv`). |
| R1 | Rare states very sparse | **Confirmed (exact)** | 178325 `intermediate_plastic=137`, `xen_like=52`; 230659 `intermediate_plastic=183` (`ground_truth_metadata.json`). |
| R2 | scIMF / PI-SDE / Squidiff are contract stubs, but listed as integrated | **Confirmed** | `benchmark/methods/{scIMF,PISDE,Squidiff}/run.py` all raise `VendoringRequired`; yet `benchmark_manuscript.tex:187–189` lists each "Yes / All three". |
| R2 | New datasets (218855/298212) rely on provisional refs | **Confirmed** | mouse 218855 + human-blood 298212 use `marker_fm_transition_silver_v1` providers, 3-edge graphs; not author-curated. |
| R2 | OSKM vs chemical use different reference protocols | **Confirmed** | 242424 OSKM = 9-edge author-cluster-matched tree; chemical = 3–4-edge marker-silver. Abstract already reports both on one scale (0.71–0.78 vs 0.59–0.66). |
| R3 | `eval_embedding.py` is an inactive stub; official metrics elsewhere | **Confirmed** | `benchmark/evaluation/eval_embedding.py` header "STATUS: INACTIVE"; official numbers come from `eval_embedding_milestone.py` → `embedding_metrics_official_silver.json`. |
| R3 | Multi-seed chain incomplete — seed run dirs missing | **Confirmed** | `benchmark/reports/official_silver/multiseed_ci.csv` exists (54 rows) but no `seed101…505` output dirs under `benchmark/results/` (only config YAMLs in `configs/runtime/multiseed/` + a prescient-internal `seed_42`). |
| R3 | Docs out of sync on test count | **Partly — reviewer mis-stated direction** | Actual suite = **18** `def test_` across `tests/`. **Both** `README.md:258` **and** `docs/reproducibility_statement.md:19` still say **16**. (Reviewer said README shows 18; it doesn't.) |
| R3 | Repro doc says server numbers "in progress" | **Confirmed** | `docs/reproducibility_statement.md:55` says server scANVI "(in progress)", but `models/scanvi_ref_gse178325/model.pt` and `…gse230659/model.pt` exist. Also two repro files diverge: `docs/reproducibility_statement.md` (scoped to the demoted annotation branch) vs `manuscript/REPRODUCIBILITY.md`. |
| R3 | Abstract overloaded | **Confirmed** | Abstract is one ~400-word block built around a single mega-sentence covering all six results. |

Net: the reviewers are right on substance everywhere. The only correction is the test-count
direction (both docs lag at 16; the suite is 18).

---

## Priority 1 — Turn the reference/OOD gate into a real, reproducible main workflow
*Closes: R1 (OOD post-hoc), Risk "gate is calibrated", Synthesis #1*

This is the manuscript's headline contribution ("reference-construction-and-validation
pipeline") and currently its softest claim, because the gate that the abstract calls
"calibrated" is computed by a side script, not the pipeline.

- **P1.1 `[code]` M — Integrate geometric OOD into `src/pipeline.py`.** Fold `score_ood.py`'s
  k-NN-latent flag into the pipeline output next to `frac_unknown_ood`, so a single run emits
  both the classifier-confidence and geometric-OOD gate decisions. Keep `score_ood.py` as a
  thin CLI wrapper around the same function.
- **P1.2 `[code]` M — k / quantile sensitivity.** Sweep k ∈ {5,15,30,50} and quantile ∈
  {0.90,0.95,0.99}; report gate decision (accept/reject) stability per transfer. The collapse
  case (73–76% flagged) and the within-modality agreement (0.76–0.86) should be shown to be
  insensitive to these knobs, or the threshold defended.
- **P1.3 `[code]` M — State-level OOD.** Report `frac_ood` per cell state, not just global, to
  show the gate isn't driven by one abundant population.
- **P1.4 `[code]` S — Positive/negative controls.** Labeled self-transfer (already ~0.94) as
  the negative control; the cross-modality collapse as the positive control. Name them as such.
- **P1.5 `[manuscript]` S — Reframe if not integrated.** If P1.1 doesn't land before
  resubmission, change "reference gate is calibrated" / "accepted-rejected gate" language to
  "post-hoc geometric OOD validation," explicitly.

## Priority 2 — Strengthen statistical robustness (seeds + threshold-free lineage metrics)
*Closes: R1 (lineage metric robustness), R3 (multi-seed provenance), Synthesis #2, Risk "traceable"*

Two separate problems: lineage metrics are fragile on 3–4-edge graphs, and the 5-seed CIs
exist as a CSV but can't be regenerated.

- **P2.1 `[code]` L — Restore seed-run provenance.** Either re-run the `configs/runtime/multiseed/`
  seeds (101–505) to repopulate `benchmark/results/.../seedNNN/`, **or** ship a seed-run manifest
  and label `multiseed_ci.csv` as a synced derived artifact with its source hash. Until then,
  the "all numbers are traceable" claim is unsupported for these rows.
- **P2.2 `[code]` L — Extend CIs to the headline claims.** `multiseed_ci.csv` covers only
  GSE178325 + GSE230659 (scenarios A/B/C). Add 5-seed CIs for the **cross-system** datasets
  (GSE218855, GSE298212) and the **pseudotime** scenarios (D/E/F) that the abstract leans on
  (PRESCIENT 0.64→0.97). The cross-system ranking-flip and the Kendall's W=0.19 claim need CIs.
- **P2.3 `[code]` M — Threshold-free lineage metrics.** On 3–4-edge graphs AUROC/Jaccard are
  threshold-sensitive. Add PR curves + AUPRC, edge-level confusion matrices, and bootstrap CIs
  over edges. State the diagonal-handling and negative-set definition explicitly.
- **P2.4 `[code]` M — Rare-state & label-threshold sensitivity (R1).** Re-score with rare
  states (137/52/183 cells) up-weighted / held-out; sweep the silver label-confidence threshold;
  add edge-drop / edge-confidence ablation. Report how rankings move.
- **P2.5 `[manuscript]` S — Label every leaderboard.** Tag each table "single-seed" or
  "5-seed ± CI" in the caption.

## Priority 3 — Downgrade unfinished content
*Closes: R2 (stub methods, external-dataset framing), Synthesis #3, Risk "three new methods", "six-plus suite"*

- **P3.1 `[manuscript]` S — Move scIMF / PI-SDE / Squidiff to "Planned extensions."** Their
  `run.py` raise `VendoringRequired`; the Methods table (`tex:187–189`) currently marks them
  "Yes / All three." Either implement before resubmission (L, unlikely) or relabel — this is the
  single most concrete factual inconsistency a reviewer can point to.
- **P3.2 `[manuscript]` S — Reframe GSE218855 / GSE298212 as external stress tests**, not
  co-equal validation. Note in-text that the mouse marker set is not author-curated and 298212
  covers Stage 1 only (already true in code comments).
- **P3.3 `[manuscript]` S — Audit "six-plus iPSC suite" (tex:115).** Not all six datasets have
  completed all three tasks. Either state per-dataset task coverage in a matrix or soften to the
  number that genuinely ran the full suite.
- **P3.4 `[manuscript]` S — Mark organoid extension + GSE175634 lineage as preliminary** if not
  fully run (175634 is currently the positive-control / cardiac silver, 5-edge).

## Priority 4 — Unify evaluation protocols & documentation
*Closes: R3 (embedding path, doc sync), R2 (cross-protocol leaderboard), Synthesis #4*

- **P4.1 `[code]` S — Resolve the embedding-evaluator naming collision.** Remove or rename the
  inactive `eval_embedding.py` (e.g. `eval_embedding_FUTURE_generative.py`) and document
  `eval_embedding_milestone.py` → `embedding_metrics_official_silver.json` as the official entry
  point in Methods.
- **P4.2 `[manuscript]` M — Stratify the leaderboard by reference protocol.** Present
  marker-silver (chemical, 3–4 edges) and author-cluster-matched OSKM (9 edges) in separate
  panels, not one unified table. The abstract's "0.71–0.78 vs 0.59–0.66" comparison should be
  explicitly cross-protocol, not cross-system alone.
- **P4.3 `[docs]` S — Fix the test count everywhere.** Suite is **18**; update `README.md:258`
  and `docs/reproducibility_statement.md:19` (both say 16).
- **P4.4 `[docs]` S — Retire the "(in progress)" server status.** `docs/reproducibility_statement.md:55`
  predates the server outputs (`models/scanvi_ref_*/model.pt` exist). Update, and reconcile it
  with `manuscript/REPRODUCIBILITY.md` — the docs/ copy is still scoped to the demoted
  annotation-robustness branch.
- **P4.5 `[docs]` S — Add a result manifest + figure-regeneration map** so README,
  REPRODUCIBILITY, manuscript, and the result files all agree (Synthesis #4).

## Priority 5 — Biological-language hygiene
*Closes: R-synthesis #5, Risk "true lineage recovery"*

- **P5.1 `[manuscript]` S — Replace "lineage recovery" with "recovery against frozen coarse
  reference graphs"** wherever silver-label alignment is described (abstract + §"Pseudotime
  denoises lineage recovery"). Reserve "true lineage" for nothing in this manuscript.
- **P5.2 `[manuscript]` S — Caveat "all numbers are traceable" (tex:80)** until P2.1 lands;
  scope it to "traceable to result files" only for the artifacts that actually regenerate.

## Priority 6 — Rewrite the abstract
*Closes: R3 (overloaded abstract), Synthesis #3*

- **P6.1 `[manuscript]` S — Cut the abstract to 2–3 conclusions:** (1) reference-transfer
  failure detected by the OOD gate, (2) lineage-ranking non-transfer across systems, (3)
  single-seed ranking instability / pseudotime denoising. Drop the per-method numbers from the
  abstract body; move them to Results. Current abstract is one ~400-word block.

---

## Suggested execution order (fastest credibility gain first)

1. **Same-day doc/manuscript fixes (all S):** P3.1 (stub methods → planned), P4.3 (test count),
   P4.4 (server status), P1.5 + P5.1 + P5.2 (claim softening), P6.1 (abstract). These remove
   every *factual inconsistency* a reviewer can verify in minutes, with zero compute.
2. **Framing fixes (S):** P3.2, P3.3, P4.2 (stratified leaderboard), P2.5 (label leaderboards).
3. **Code, no server needed (M):** P1.1–P1.4 (OOD into pipeline + sensitivity), P2.3 (threshold-free
   lineage), P4.1 (embedding stub rename).
4. **Compute-bound (L):** P2.1 (restore seed runs / manifest), P2.2 (extend CIs to cross-system +
   pseudotime), P2.4 (rare-state / label-threshold / edge-drop). These are the long pole — start
   the seed re-runs early since everything in Priority 2 depends on them.

## Open decisions for you
- **scIMF/PI-SDE/Squidiff:** implement the vendored integrations now, or move to planned
  extensions? (Implementing is the only L item in Priority 3 and may not be worth it pre-resubmission.)
- **Seed provenance:** re-run all seeds to restore raw dirs, or ship a manifest + mark the CSV
  as a derived artifact? The honest-and-fast path is the manifest; the bulletproof path is the re-run.
