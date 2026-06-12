# Release notes — v1.0.0-rc1

**Tag:** `v1.0.0-rc1` · **Date:** 2026-06-12 · **License:** MIT

Code, manuscript, and reproducibility materials for the pluripotency-domain temporal
trajectory-inference benchmark — *"Does it transfer? A pluripotency-domain benchmark of
temporal trajectory-inference methods across iPSC reprogramming systems, with
scANVI-validated cell-state references."* This is the reviewer-revision (round 2) snapshot.

## What this release contains
- The benchmark framework, capability-gated evaluators, and frozen silver-reference
  providers (`benchmark/`, `src/`).
- The geometric out-of-distribution (OOD) reference gate, integrated into the main
  pipeline (`src/pipeline.py` → `src/annotation/ood.py`) and reported with per-state
  breakdown and a k/quantile sensitivity sweep (`scripts/ood_sensitivity_sweep.py`).
- The manuscript source and compiled PDF (`manuscript/`), the point-by-point
  response to reviewers, and the claim-to-source map (`manuscript/REPRODUCIBILITY.md`).
- 24 passing unit tests (`pytest`), including the OOD gate and threshold-free
  lineage-robustness companions.

## Highlights of this revision
- **OOD gate** is now a pipeline gate, not a post-hoc diagnostic; the cross-modality
  reject decision is stable across k∈{5,15,30,50} and q∈{0.90,0.95}.
- **scIMF, PI-SDE, Squidiff** are scoped strictly as *planned extensions* (their adapters
  raise `VendoringRequired`); they appear in no results table or figure.
- Conclusions are framed explicitly as a **silver-reference benchmark and cautionary
  analysis**, not biological ground truth.
- **Multi-seed CIs (Fig. 8)** are labeled a derived artifact, backed by the 90 run logs
  (`logs/run_multiseed_benchmark.122804086.*.log`), SHA-256 config checksums
  (`benchmark/reports/official_silver/multiseed_config_checksums.sha256`), and a documented
  HPC rerun recipe.
- Documentation synced (README/REPRODUCIBILITY/TODO); large regenerable artifacts
  excluded via `.gitignore`.

## Datasets (all public)
NCBI GEO: GSE178325, GSE230659, GSE242424, GSE175634, GSE298212, GSE218855.
ArrayExpress: E-MTAB-7552.

## Reproduce
```bash
pip install -r requirements-annotation-branch.txt pytest
pytest                                   # 24 tests
python scripts/ood_sensitivity_sweep.py  # OOD k/quantile sweep (local CPU)
```
Multi-seed regeneration (HPC) and the full claim-to-source map are documented in
`manuscript/REPRODUCIBILITY.md`.

## Known limitations
Lineage/embedding metrics are scored against marker-defined or scANVI-transferred silver
references (not author-validated biological ground truth). Headline tables are single-seed
point estimates; 5-seed CIs cover the two primary datasets, with cross-system and
pseudotime CIs as outstanding work.

---

## Zenodo / OSF deposit description (copy-paste)

**Title:** Does it transfer? A pluripotency-domain benchmark of temporal trajectory-inference methods across iPSC reprogramming systems (v1.0.0-rc1)

**Authors:** Kuma (Graduate School of Frontier Sciences, The University of Tokyo)

**Description:**
Software and analysis release for a single-domain (induced-pluripotency) benchmark of
temporal single-cell trajectory-inference methods, built on the scTimeBench three-task
design (forecast accuracy, embedding coherence, lineage fidelity). The release contributes
(1) a reference-construction-and-validation pipeline in which silver cell-state references
are built by scANVI/scArches label transfer and screened by a geometric out-of-distribution
gate that quantifies when a transferred reference is trustworthy, and (2) an empirical
analysis showing that method rankings — lineage rankings in particular — do not transfer
across reprogramming systems, and that single-seed leaderboards are unstable. All evaluated
numbers are read from result files; multi-seed confidence intervals are provided as a
derived artifact with run logs, config checksums, and a documented rerun recipe. Results are
framed as a silver-reference benchmark and cautionary analysis rather than biological
ground truth. Public datasets: GEO GSE178325, GSE230659, GSE242424, GSE175634, GSE298212,
GSE218855; ArrayExpress E-MTAB-7552. License: MIT.

**Keywords:** single-cell, trajectory inference, benchmarking, scANVI, scArches,
out-of-distribution, induced pluripotency, reprogramming, reproducibility
