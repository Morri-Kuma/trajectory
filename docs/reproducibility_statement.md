# Reproducibility & Code-Completeness Statement

**Date:** 2026-06-12 · **Scope:** annotation-robustness branch (`src/`, `config.yaml`, jobs, tests)

> The annotation-robustness branch is a **side sensitivity** of the main scTimeBench
> conditions benchmark, not the headline paper. The canonical claim-to-source map for the
> manuscript is `manuscript/REPRODUCIBILITY.md`; this file covers the `src/` annotation
> pipeline and its tests.

This statement documents what a reviewer can reproduce, and the single class of
results that inherently requires GPU/server compute.

## Code architecture (complete)
Modular package `src/` (24 modules, all compile):
`data/` (load, count-reconstruction, scANVI input builder), `preprocessing/`
(normalize→HVG→PCA→kNN, scanpy or sklearn fallback), `annotation/`
(scVI/scANVI/scArches real path + KNN surrogate + dispatcher), `trajectory/`
(DPT/PAGA + dependency-light surrogates), `evaluation/` (the six annotation
comparisons), `plotting/`, `utils/` (config/seed/io), and `pipeline.py`
(config-driven orchestrator). One entry point (`config.yaml`) with profiles
`test`, `server`, `server_dryrun`, `local_real`.

## Reproducible **without** GPU (laptop / CI)
- **Unit tests:** `pytest` → 24 tests across config loading, exact count
  reconstruction (normalize_total→log1p inversion), the comparison metrics
  (composition TV, NMI/ARI/correspondence recall, marker AUROC, pseudotime
  correlation, disagreement), trajectory utilities (graph Jaccard, geodesic
  pseudotime monotonicity, surrogate PAGA), the KNN surrogate, preprocessing,
  the geometric OOD gate (score + transfer-level accept/reject decision +
  per-state breakdown), and the lineage-robustness companions (edge confusion,
  PR curve, bootstrap CI over edges).
- **Synthetic smoke test:** `python scripts/run_smoke_test.py` → full pipeline on
  generated data, 27 output artifacts, asserts presence.
- **Real-data execution (surrogate):** `mode: local_real` runs the entire
  annotation-comparison + trajectory pipeline on the **real** prepped reference
  and query `.h5ad` (subsampled, CPU, KNN surrogate in place of scANVI),
  producing real composition/agreement/pseudotime numbers. This proves the code
  executes correctly on the actual data; the surrogate labels are **not** a
  manuscript result.

## Requires GPU/server (inherent)
- **scANVI/scArches training + mapping** (`mode: server`): scvi-tools + a CUDA
  PyTorch. Produces the reference model, the per-cell `scanvi_label`/confidence
  predictions, and therefore the *real* numbers behind Figures 2–5. Reproduced
  via the Shirokane `qsub` chain (`jobs/shirokane/submit_scanvi_annotation_branch.sh`).
- Full-data DPT/PAGA/Slingshot and any differential-expression panels.

## Determinism & provenance
- Global seed threaded through numpy/torch/scanpy (`utils/seeding.py`).
- Every pipeline run stamps the config hash and annotation method into outputs
  (`uns['scanvi_branch']`, `summary.json`).
- Inputs are never re-downloaded; counts are losslessly reconstructed from the
  stored log-normalized matrices + `obs['total_counts']` (verified corr = 1.0000).

## Environment
- `environment.yml` (conda) and `requirements-annotation-branch.txt` (pip core).
- Core stack (numpy/pandas/scipy/scikit-learn/anndata/matplotlib/pyyaml/h5py/pytest)
  is sufficient for tests + `local_real`; scanpy is auto-used if present;
  scvi-tools + torch only for `server`.

## Honest status line
The **code** is complete, tested, and runs on real data locally. The server scANVI
run is **complete**: the trained reference models are present
(`models/scanvi_ref_gse178325/model.pt`, `models/scanvi_ref_gse230659/model.pt`) and
the geometric OOD diagnostics are in `results/annotation_branch/`. Nothing in the
repository fabricates the paper's numbers — figures are generated from whatever the
run actually produced. The geometric OOD gate is now computed **inside** the main
pipeline (`src/pipeline.py` → `src/annotation/ood.py:ood_gate`), not only post-hoc.
