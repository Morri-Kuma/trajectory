# Response to Reviewers

We thank the three reviewers for an unusually precise and constructive read. The reviews
converged on five issues, and we have addressed each. Below, every point is answered with
the concrete change made and the file(s) touched. We distinguish **[done]** changes (in this
revision) from **[queued]** items that require a server/seed run and are now scaffolded,
labeled, and reproducible from the repository. We have been careful not to add any number we
cannot trace to a result file.

A short summary of the substantive changes:

- The geometric OOD gate is now computed **inside the main pipeline**, not only post-hoc, and
  is reported with a per-state breakdown and a **(k, q) sensitivity sweep** on the real inputs.
- scIMF, PI-SDE, and Squidiff are reclassified as **planned extensions** and removed from every
  results table and figure.
- The OSKM and marker-silver systems are presented **stratified** by reference protocol, and the
  mouse/blood systems are framed as **external stress tests**.
- The abstract is rebuilt around **three conclusions**; biological language is corrected to
  "recovery against frozen coarse reference graphs."
- Documentation, test counts, the embedding entry point, and multi-seed provenance are
  **synced**, with a new seed-run manifest and a derived-artifact label on the CI table.

---

## Reviewer 1 — Technical soundness

**1.1 Weak ground truth; add rare-state, label-threshold, and edge sensitivity.**
We agree the silver references are the central limitation and now say so plainly. The
Limitations paragraph states that lineage fidelity "measures alignment with a coarse reference
graph rather than ground-truth biology," and notes the reference graphs are small (3–4 edges
for the marker-silver systems, 9 for OSKM), so the metrics are read with that sparsity in mind.

- **[done]** Threshold-free / robustness companions added as unit-tested code:
  `benchmark/evaluation/lineage_robustness.py` provides explicit **edge-level confusion** (TP/FP/FN/TN
  over all ordered off-diagonal pairs, i.e. an explicit negative set), **PR-curve points**, and a
  **bootstrap CI over the candidate-edge set** so a 3–4-edge graph reports an interval, not a bare
  point. Tests: `tests/test_lineage_robustness.py` (4 tests). We also note the main evaluator
  already reports the threshold-free **AUROC and AUPRC** (`eval_lineage.py:compute_auroc/compute_auprc`).
- **[queued]** Regenerating per-run bootstrap CIs, rare-state hold-out/up-weighting, and the
  label-confidence-threshold sweep across the benchmark is a downstream run over the saved
  transition matrices; the code path is in place and labeled.

**1.2 The OOD gate is post-hoc, not a pipeline gate.**

- **[done]** The geometric OOD gate is now **integrated into the main pipeline**. `src/pipeline.py`
  calls a shared `src/annotation/ood.py:ood_gate`, emitting `frac_ood_geometric`, the
  accept/reject `gate_decision`, an enrichment-over-reference statistic, and a **per-state OOD
  breakdown** into `summary.json` and a per-query `ood_gate.json`, alongside the existing
  classifier-confidence flag. `scripts/score_ood.py` is now a thin wrapper over the same
  `ood_gate`, so the server (scANVI latent) and pipeline (PCA latent) share one code path.
  Verified end-to-end by the synthetic smoke test and unit tests `tests/test_ood.py` (4 tests).
- **[done]** k / quantile sensitivity, positive/negative controls, state-level OOD: a sweep over
  k∈{5,15,30,50} and q∈{0.90,0.95,0.99} on the real prepped inputs (PCA-latent proxy, the same
  proxy used by the existing recalibration demo) shows the cross-modality **reject decision is not
  an artifact of (k, q)**: the chemical query is flagged **6.5–9.7× more often than the reference's
  own self-flag baseline across every k and across q∈{0.90, 0.95}**, collapsing toward parity only
  in the extreme q=0.99 tail. Script `scripts/ood_sensitivity_sweep.py`; outputs
  `results/test_outputs/ood_recalibration/ood_sensitivity_sweep.{csv,json}`. The Results section
  now names the cross-modality case the **positive control** and the within-modality / labeled
  self-transfer the **negative controls**, and the per-state breakdown is cited
  (`ood_by_group_summary.json`).
- We have softened the claim accordingly: the subsection is retitled "The scANVI reference gate
  flags cross-modality transfers that confidence misses," and we describe it as a *screened,
  controlled* gate, explicitly noting that a continuous calibration curve across intermediate
  transfer distances is future work.

**1.3 Lineage graph metrics need robustness checks (threshold-free metrics, PR curves,
edge-level confusion, bootstrap CIs).**

- **[done]** AUROC and AUPRC (threshold-free) are already computed and reported; the new
  `lineage_robustness.py` adds the explicit edge confusion, PR-curve points, and bootstrap CIs the
  reviewer requested, with unit tests.
- **[queued]** Producing these for every method/dataset/scenario and adding them to the figures is
  a downstream sweep over the saved predictions; the entry points exist and are tested.

## Reviewer 2 — Originality and significance

**2.1 GSE218855/GSE298212 rely on provisional marker references; frame as external stress tests.**

- **[done]** The two systems are now framed as **external stress tests** rather than co-equal
  validation. The Limitations paragraph states the mouse/blood milestone markers are "canonical
  rather than author-curated and so those two systems are framed as external stress tests," and the
  abstract no longer presents them as validation of equal strength. Their role in the cross-system
  *non-transfer* argument (which only needs internal consistency, not author labels) is retained.

**2.2 scIMF, PI-SDE, Squidiff are contract templates that raise `VendoringRequired`, yet listed as
integrated.**

- **[done]** Reclassified as **planned extensions**. The Methods prose now states their adapters
  "raise an explicit `VendoringRequired` contract" and that the methods "are *not* evaluated here and
  are excluded from every results table and figure." Table 2 is split into an **Evaluated roster**
  block and a **Planned extensions (adapter contract only — not evaluated here)** block. The
  Discussion repeats this. They appear in no results table or figure.

**2.3 OSKM and chemical use different reference protocols; present stratified, not unified.**

- **[done]** The OSKM (GSE242424) system is presented **stratified** from the marker-silver systems
  throughout: the cross-system figure marks OSKM as evaluated under the author-cluster ground-truth
  protocol and "not strictly comparable" to the marker-silver bars; the dataset-suite text and
  Limitations state OSKM is "shown stratified rather than pooled with the marker-silver systems." The
  concordance statistic (Kendall's W) is computed over the **four marker-silver systems only**.

## Reviewer 3 — Readability and reproducibility

**3.1 Abstract overloaded.**

- **[done]** The abstract is rebuilt around **three conclusions** — (1) transferred references fail
  detectably across modality (OOD gate), (2) lineage rankings do not transfer across systems, (3)
  single-seed leaderboards are unstable (with the pseudotime caveat) — and the dense per-method
  numbers were moved into Results.

**3.2 Embedding evaluation path is unclear; `eval_embedding.py` says inactive while official reports
read from `embedding_metrics_official_silver.json`.**

- **[done]** The inactive stub was **renamed** `eval_embedding.py` →
  `eval_embedding_generative_inactive.py`, its docstring now states it is *not* the evaluator behind
  the reported numbers, and the dispatcher import and the scNODE comment were updated. The official
  entry point — `eval_embedding_milestone.py` →
  `embedding_milestone_eval/embedding_metrics_official_silver.json` — is now named explicitly in
  Methods and in `REPRODUCIBILITY.md`.

**3.3 Multi-seed reproducibility chain incomplete; seed run dirs not present.**

- **[done]** We added a **seed-run manifest**, `benchmark/reports/official_silver/multiseed_manifest.csv`,
  listing all **90** configured seed-runs (3 methods × 2 datasets × A/B/C × 5 seeds), each with its
  config file, expected output directory, and a `raw_output_present` flag (currently all "no").
  `multiseed_ci.csv` is now **labeled a derived artifact** in `REPRODUCIBILITY.md`, with a
  "single-seed vs multi-seed" note, and re-running the `configs/runtime/multiseed/` sweep regenerates
  the raw outputs. The manuscript's Data-and-code section reflects this.

**3.4 Documentation out of sync (16 vs 18 tests; "server numbers in progress").**

- **[done]** The suite is now **24 tests** (the 18 prior + 6 new — 2 OOD-gate, 4 lineage-robustness;
  precise count verified by `pytest`). `README.md` and `docs/reproducibility_statement.md`
  are updated from "16" to "24". The stale "server scANVI run (in progress)" line is replaced — the
  trained reference models exist (`models/scanvi_ref_gse178325/model.pt`,
  `models/scanvi_ref_gse230659/model.pt`) and the OOD diagnostics are present. The
  `reproducibility_statement.md` scope note now flags that it covers the annotation branch (a side
  sensitivity) and points to `manuscript/REPRODUCIBILITY.md` as the canonical claim map.

## Cross-review synthesis

1. **Strict, reproducible OOD gate** — done: integrated into the pipeline, with k/quantile sweep,
   per-state OOD, and named positive/negative controls (R1.2).
2. **Statistical robustness** — partly done: threshold-free metrics already reported; edge confusion +
   bootstrap CIs added as tested code; the manifest makes the 5-seed CIs reproducible. Extending CIs to
   the cross-system and pseudotime claims is the **main queued run**, and the headline tables are now
   explicitly labeled single-seed.
3. **Downgrade unfinished content** — done: scIMF/PI-SDE/Squidiff → planned extensions; organoid and
   GSE175634 lineage flagged as not-yet-run; mouse/blood framed as stress tests.
4. **Unify evaluation + documentation** — done: official embedding entry point named and stub renamed;
   multi-seed provenance and manifest added; README / REPRODUCIBILITY / manuscript aligned.
5. **Biological caution** — done: "lineage recovery" reframed as "recovery against frozen coarse
   reference graphs"; the gate is described as screened/controlled rather than calibrated.

## On the four "risk / unsupported claims"

- *"All numbers are traceable"* — narrowed: each **single-seed** number is read from a result file
  (`REPRODUCIBILITY.md` maps them); the multi-seed CIs are now explicitly a **derived** artifact with a
  manifest for regeneration.
- *"The scANVI reference gate is calibrated"* — softened to a screened, controlled gate validated by a
  positive and a negative control, with a (k, q) sensitivity sweep; full calibration across intermediate
  distances is named as future work.
- *"Six-plus iPSC suite"* — replaced with a per-dataset run-status statement; no aggregate marketing
  claim remains.
- *"Three newly added methods are not runnable"* — acknowledged: they are now planned extensions,
  excluded from all results.

---

# Response to Reviewers — Round 2

We thank the reviewers for confirming the round-1 improvements. The remaining items were
release-engineering and framing; all are addressed.

**R1 — keep scIMF/PI-SDE/Squidiff out of the formal benchmark.** Confirmed and reinforced.
The three methods appear in **no** title, abstract, results table, or figure (verified
programmatically). They remain only in a clearly separated Table 2 block, "Planned extensions
(adapter contract only — not evaluated here)," and a new Discussion paragraph ("Scope of the
claims") states they "raise an explicit `VendoringRequired` contract, and contribute to no
result in this paper."

**R2 — frame as a silver-reference benchmark, not biological ground truth.** Done. A new
"Scope of the claims" paragraph states the work is "a *silver-reference benchmark and
cautionary analysis*, not a biological-ground-truth study … not a definitive biological
ordering of methods or a mechanism-discovery claim." We also softened the OSKM protocol name
from "author-cluster ground-truth protocol" to "author-cluster reference protocol."

**R3 — reproducibility (the main remaining risk).**

- *README 24-vs-18 inconsistency:* fixed. A second line in `README.md` (and a stale `TODO.md`
  entry) still said 18/16; both now read **24**, matching `pytest`.
- *Raw multi-seed outputs / Fig 8:* strengthened with all three forms of evidence the reviewer
  suggested. (i) **The 90 seed-runs executed** — the run logs are in the repo
  (`logs/run_multiseed_benchmark.122804086.{1..90}.log`), each recording its per-seed output
  path. (ii) **Config checksums:** `multiseed_config_checksums.sha256` pins all 90 seed configs.
  (iii) **Documented HPC rerun recipe** in `REPRODUCIBILITY.md` (generate configs → verify
  checksums → `qsub run_multiseed_benchmark_array.sh` → `aggregate_multiseed_ci.py` →
  `make_multiseed_figure.py`). The manifest now maps each run to its config hash, output dir,
  and log. The raw per-seed metric JSONs remain unsynced for size; the recipe regenerates them.
- *Clean git state + citable release:* the chaotic `git status` traces to a **stale
  `.git/index.lock` left by a crashed git operation on 2026-06-01**, which wedged the index.
  A one-command release script, `scripts/prepare_release.sh`, clears the lock, re-applies the
  hardened `.gitignore` (large regenerable artifacts excluded), refuses to stage any file >5 MB,
  and creates a clean commit + annotated tag. (It must be run on a local machine; the editing
  sandbox cannot delete files inside `.git`.) After tagging, a GitHub release or Zenodo/OSF
  upload gives a citable archive.

**Traceability.** Every figure (Fig 1–8) and headline table maps to a source file and a
regeneration script in `manuscript/REPRODUCIBILITY.md`; single-seed numbers are read directly
from result files, and the multi-seed CIs are labeled a derived artifact with the provenance
above.

---

# Response to Reviewers — Round 3

We thank the reviewers; this round was the most useful yet. The headline addition is a
**negative control** for the lineage metric; the rest sharpen framing and guidance.

## Reviewer 1 — technical evidence

- **Negative control (prioritized #1) — done, and informative.** `scripts/negative_control_lineage.py`
  recomputes the *production* single-step graph-similarity AUROC on every real predicted
  transition matrix (all scenarios A–F, 3 methods × 6 systems, 72 runs) under two null inputs: a
  **label-permutation null** (2000 draws) and a **random row-stochastic matrix** control
  (2000 draws). Result (new Results subsection + Fig. S1): the empirical chance level is
  AUROC ≈ **0.62** (permutation) / ≈ **0.63** (random), *above* the nominal 0.5 because the
  4–5-node reference chains are easy to match partially. Against this empirical null, the
  strong-signal systems (OSKM, cardiac) are clearly above chance (permutation p ≤ 0.005),
  while the human hADSC silver systems frequently are **not** — only **23/72** method×system×scenario
  runs exceed the null's 95th percentile, and scNODE on GSE230659-B (0.44) sits below its own null. Notably PRESCIENT's pseudotime jump on GSE230659 is real (clears the null at p≤0.045, AUROC 0.95–1.00) while most other pseudotime gains do not (4/18). So the metric is
  discriminating (not permissive), and absolute lineage AUROC on the silver systems must be
  read against a permutation floor near 0.6, not 0.5. The tie-correct AUROC helper is unit-
  tested against sklearn (`tests/test_lineage_robustness.py`).
- **Tighten the central message — done.** The abstract now closes on one methodological
  claim: *because rankings move with the task, the reference construction, and the random
  seed, a domain benchmark should report task stratification, an OOD check on reference
  transfer, and seed uncertainty rather than a single leaderboard.*
- **Reduce "ground truth" — done.** The OSKM protocol is now "author-cluster reference
  protocol"; the only remaining uses of the phrase are explicit disclaimers in the
  Limitations/Scope paragraphs. Terminology is consistently silver/proxy/reference.

## Reviewer 2 — novelty and impact

- **Distinguish from prior benchmarks — done.** The Introduction now names the four
  distinctive features: single-domain (pluripotency) focus that makes ranking *transfer*
  testable, task stratification, OOD-gated reference transfer, and quantified single-seed
  instability.
- **"No single best method" as a positive — done.** A new Discussion paragraph ("No single
  best method is the finding, not a gap") frames task-dependent winners as actionable
  guidance for practitioners rather than a shortcoming.
- **Don't let planned methods dilute — confirmed.** scIMF/PI-SDE/Squidiff appear in no
  title/abstract/results table/figure; they remain only in the labeled "Planned extensions"
  block and the Scope paragraph.

## Reviewer 3 — writing and figures

- **Practical guidance — done.** A "Practical guidance" Discussion paragraph now states:
  prioritise OT-loss metrics for forecasting (the most portable task, W = 0.75); do not judge
  lineage by embedding ARI (the tasks rank methods differently); read absolute lineage AUROC
  against the permutation floor; reject cross-modality reference transfer by the enriched OOD
  flag and accept only within-modality gate-passing transfers; treat single-seed leaderboards
  as draws.
- **Abstract less engineering, more conclusions — done.** The abstract is organised around
  three findings plus the methodological central message; component-directory phrasing was
  removed in the earlier round.
- **One-conclusion-per-figure — layout plan provided.** `figure_plan.md` gives a concrete
  consolidation to five main figures (design / reference-OOD / primary tasks / cross-system
  instability / seed robustness + metric validity), with pseudotime and detailed tables moved
  to the Supplement. This is a montage/layout pass over existing PNGs; underlying numbers and
  regeneration scripts are unchanged.

## Status
26 unit tests pass; the manuscript compiles (14 pp, 0 undefined references); every new
number is traceable via `manuscript/REPRODUCIBILITY.md`.

---

# Follow-up: "why don't the tables/figures show all 6 datasets and all 8 methods together?"

Three reasons, two structural and one we have now fixed:

1. **Capability gating.** Of the 8-method roster only scNODE, PRESCIENT and MIOFlow are
   projection-capable *and* completed, so only those three can populate a forecast/embedding
   table. WOT and CellRank2 are lineage-only (optimal transport / fate mapping — no projected
   expression or embedding), and scIMF/PI-SDE/Squidiff are planned (not run). Table 3 therefore
   has 3 methods by construction.
2. **Role / protocol stratification.** Table 3 / Fig. 2 are the *primary* chemical-reprogramming
   benchmark (2 hADSC systems); the other systems appear in the cross-system table/figure, and
   OSKM (GSE242424) is kept on its **author-cluster** reference protocol rather than pooled with
   the marker-silver systems — a stratification an earlier reviewer explicitly requested, and one
   the paper's "rankings don't transfer" thesis requires (a single pooled average would hide the
   main result). Forecast Wasserstein is also on very different scales across systems (≈0.04 vs
   ≈160–260 on GSE230659), so a single forecast axis is not meaningful.
3. **The fixable gap — now addressed.** We added **Table~\ref{tab:master}, a complete results
   matrix**: all five completed methods × all six systems × the three tasks (forecast WD,
   embedding ARI, lineage AUROC), with lineage-only and planned methods marked and OSKM flagged
   as a separate protocol. The main-text tables/figures are task- and role-stratified views of
   this one matrix. Source: `benchmark/reports/official_silver/master_results_matrix.csv`
   (regenerated from the per-run metric JSONs).
