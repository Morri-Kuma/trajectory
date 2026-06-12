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
