# TODO — trajectory project (annotation-robustness MSc slice)

Running task list. Updated 2026-06-09. ✅ done · 🔄 in progress · ⬜ next · 🧊 server-only

## Done this iteration
- ✅ Project audit report — `docs/reports/project_audit_report.md`
- ✅ Literature search + 8 verified new refs — `reference/literature_review_extra.md`
- ✅ Dataset inventory (3 required + 5 additional, verified) — `datasets/dataset_inventory.{md,csv}`
- ✅ Three-month research design — `research_design_three_month_project.md`
- ✅ scANVI annotation-branch design — `scanvi_annotation_branch.md`
- ✅ Modular code under `src/` (data, preprocessing, annotation, trajectory, evaluation, plotting, utils)
- ✅ `config.yaml` (test vs server profiles)
- ✅ Manuscript figure plan — `figure_plan.md`
- ✅ **Small local pipeline test passes** — `scripts/run_smoke_test.py` → `results/test_outputs/annotation_branch/`
- ✅ Deletion proposal (awaiting approval) — `deletion_proposal.md`

## Next (local, no big downloads)
- ✅ Real benchmark inputs wired: `server_dryrun` (HVG2000 subsample) + `local_real` (prepped full-gene inputs) profiles; both run end-to-end locally with the surrogate.
- ✅ `notebooks/01–05` thin drivers over `src/` present and run.
- ✅ Unit-test suite (`tests/`, `pytest.ini`): 24 tests across config, count-reconstruction, comparison metrics, trajectory, surrogate, preprocessing, the geometric OOD gate, and lineage-robustness companions — all pass.
- ✅ Pipeline reads the REAL labels from prepped inputs (`author_cluster_label` ref + `final_milestone_label_coarse` queries) under `local_real`; verified on 3k-cell subsamples.
- ✅ `environment.yml` added (core + scanpy + scvi-tools).

## Server-only (full data, GPU)
- ✅ **Shirokane job chain WRITTEN (2026-06-09):** `jobs/shirokane/submit_scanvi_annotation_branch.sh` runs prep→train→map+compare with `-hold_jid` deps (see `jobs/shirokane/docs/scanvi_annotation_branch_jobs.md`). Remaining server items below are now "submit + verify".
- 🧊 **Prep inputs FIRST (dry-run finding 2026-06-09):** build shared-gene-space scANVI inputs via `src/data/build_scanvi_inputs.py` → job `run_scanvi_prep_inputs.sh` (ready). HVG2000 inputs are gene-disjoint (ref∩query ≈181–238); full-gene shares 17,829 incl. all markers. Feed raw counts to scVI.
- ✅ **scANVI reference trained (2026-06-09): recovery_accuracy = 0.990** (gate 0.90 passed); model saved to `models/scanvi_reference_gse242424/`.
- ⚠️ **scArches map RAN on full data (2026-06-09) but DEGENERATE:** ~76–96% of query cells confidently collapsed to OSKM-only `hOSK` (no chemical analogue); scANVI confidence did NOT flag the OOD. Real, interpretable NEGATIVE result — see `docs/scanvi_run_results_interpretation.md`. NOT a bug (recovery 0.99).
- ✅ **Positive control DONE (GSE175634 self-transfer): test accuracy 0.944**, per-state recall 0.93–0.96 — method works when systems match; the GSE242424→chemical failure is the cross-protocol gap, not a bug.
- ✅ **OOD recalibration DONE**: `src/annotation/ood.py` (unit-tested) flags 38.7%/43.6% of query cells on real data vs scANVI ≤0.3% (`scripts/demo_ood_recalibration.py`). scANVI-latent re-score DONE (2026-06-10): geometric OOD **75.9%/73.0%** vs scANVI ≤0.3% (`results/annotation_branch/ood_summary_scanvi_latent.json`).
- ✅ **Same-protocol DONE: correspondence 0.863 (230→178) / 0.757 (178→230)** vs 0.01–0.04 cross-protocol (20–60×). Within-protocol mapping is trustworthy for abundant states; rare states (hCiPS/intermediate) limited by reference abundance.
- 🧊 Full DPT + PAGA + Slingshot under both annotations; generative-method re-scoring (Fig S2).
- 🧊 Positive control (GSE175634) + external validation (GSE298212) + organoid generalization (E-MTAB-7552).
- 🧊 Differential expression / regulators for annotation-sensitive transitions (Fig 5d).

## Decisions needed from you
- ✅ **Deletion approval executed (2026-06-09):** Tier A+B done — ~34 GB freed (benchmark/archive 32 GB + unreferenced data/gse178325_human/rna_seq/ 2.6 GB + stray logs/.git, _expanded_official_backup, _invalid_lineage_label_space, empty smoke/, caches). Kept canonical rna_seq_10x/. Tier C (stale data/ audit copy; gse230659( rename) NOT done — optional next.
- ✅ **Reference label key (decided 2026-06-09): `author_cluster_label`** (9 author ATAC->RNA clusters). config.yaml updated; correspondence_map already uses these names. `final_milestone_label_coarse` (10) kept as a cheap robustness variant (differs only by splitting the Intermediate cluster).
- ✅ **E-MTAB-7552 (decided 2026-06-09): confirmed real & public** (ArrayExpress + quadbio GitHub); **DEFERRED from the core server run** — organoid biology has no correspondence to the OSKM reference; stays an inventoried dataset for an optional generalization extension (Fig S4) only.

## Notes / caveats
- The local run proves the **pipeline executes**; it does **not** establish any biological result. Biology comes from the server run on full data.
- OSKM (GSE242424 reference) vs chemical (queries) label spaces differ → comparisons use correspondence-aware metrics, not naive accuracy.
- File edits on this working copy are best made via full-file writes (the in-place edit path truncated files mid-write here; resolved).
