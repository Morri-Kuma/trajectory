## Per-model wiring status (2026-06-12)

- **Squidiff** — `run.py` now calls the REAL Squidiff API (`Squidiff==1.0.8`, import name
  `Squidiff` capital-S): PCA-space diffusion via `create_model_and_diffusion` /
  `run_training` / `p_sample_loop`, with `obs['Group']` = timepoint index. Needs a GPU
  smoke-test; the one modelling choice to confirm is the perturbation-encoder -> timepoint
  mapping (x_start = t0 cell, Group = target tp). Array tasks 25-36.
- **scIMF** / **PI-SDE** — engines + plumbing done; the `_fit_model/_simulate/_encode`
  hooks still need the authors' source, which is not reachable from the dev sandbox. Get
  the repo URL from each paper's Code Availability:
    * scIMF — Jiang, Li et al. 2026, PLOS Comput Biol, DOI 10.1371/journal.pcbi.1013916
    * PI-SDE — Jiang & Wan 2024, Bioinformatics 40(Suppl.2):ii120-ii127 (ISMB 2024)
  `git clone` it into the `<X>_module/` dir, then wire the three hooks (templates in each
  hook's docstring). Array tasks: scIMF 1-12, PI-SDE 13-24.

---

## UPDATE (run.py glue is now scaffolded)

`benchmark/methods/_generative_common.py` now contains the **model-agnostic** runner (Forecast/Embedding/Lineage artifact construction, ported from scNODE/run.py), and each `benchmark/methods/<X>/run.py` is implemented to delegate to it. The five contract functions are done; what remains per model is wiring **three engine hooks** to the vendored source:

- Squidiff (`SquidiffEngine`): `_fit_model`, `_sample_pcs` (PCA round-trip already done; `pip install squidiff==1.0.8`).
- scIMF (`ScIMFEngine`) / PI-SDE (`PISDEEngine`): `_fit_model`, `_simulate`, `_encode` (git clone the authors' repo into the `<X>_module/` dir).

Until those hooks are wired they raise a precise `NotImplementedError`/`RuntimeError` (no fabricated outputs); the adapter catches it and writes `status=failed`. Each hook carries a `VERIFY:`/template comment with the exact call to confirm against the source. After wiring, smoke-test one task on a GPU node before the full array:

```bash
SGE_TASK_ID=1 bash jobs/shirokane/run_method_expansion_benchmark_array.sh
```

---

# Vendoring the three new generative models (scIMF, PI-SDE, Squidiff)

These three projection-capable models are integrated at the framework level (adapters,
runtime configs, capability gating, dispatch registry, and the benchmark array job all
exist), but their **upstream source is not vendored**, so every run currently fails fast
with a `VendoringRequired` error and writes `status=failed` (no garbage outputs). This
guide is the checklist to make them runnable on Shirokane, where the source and GPU live.

## The contract (identical for all three)

Each `benchmark/methods/<X>/run.py` must implement the same five functions the thin
adapter calls, producing the same standardized artifacts as
`benchmark/methods/scNODE/run.py` — **use that file as the structural template**:

1. `prepare_data(train_adata, time_key, train_times)` → per-timepoint expression arrays
2. `train_or_load(train_data, train_tps, n_genes, cfg, model_cache)` → trained model
3. `run_forecast_accuracy(...)` → writes `projected_expression.npy`, `forecast_metrics.json`, `per_timepoint_forecast_metrics.csv`
4. `run_embedding_coherence(...)` → writes `embedding.npy`, `projected_embedding.npy`, `next_timepoint_embedding.npy`, `embedding_metrics.json`, `projected_cluster_labels.csv`
5. `run_lineage_fidelity(...)` → writes `state_transition_matrix.csv`, `lineage_graph_edges.csv`

The output-writing, forecast-metric, embedding-coherence, and STM→graph helpers in
`scNODE/run.py` are **method-agnostic** once the model can project cells to held-out
time points — copy them verbatim and only swap in the model's train + simulate calls.

## Per-model source

| Model | Vendor under | Source | Notes |
|---|---|---|---|
| scIMF | `benchmark/methods/scIMF/scIMF_module/` | Jiang, Li et al. 2026, *PLOS Comput Biol* 22(1):e1013916 — authors' repo (or scTimeBench `methods/` submodule) | Transformer + interacting mean-field neural SDE; scTimeBench's strongest forecaster |
| PI-SDE | `benchmark/methods/PISDE/PISDE_module/` | Jiang & Wan 2024, *Bioinformatics* 40(Suppl.2):ii120 — authors' repo | Physics-informed neural SDE |
| Squidiff | `benchmark/methods/Squidiff/Squidiff_module/` | He, Zhu et al. 2026, *Nat Methods* 23(1):65 — authors' repo | Conditional diffusion; run in PCA space |

Mirror the layout of `benchmark/methods/scNODE/scNODE_module/` (the working example).

## Steps on Shirokane

```bash
cd /home/xzy0723/projects/trajectory

# 1. Vendor each model's source into its <X>_module/ directory (git clone / copy),
#    then wire the 5 functions in benchmark/methods/<X>/run.py to it.

# 2. Sanity-check the adapters import and the configs build:
python -c "from benchmark.adapters.scimf_adapter import ScIMFAdapter; \
           from benchmark.adapters.pisde_adapter import PISDEAdapter; \
           from benchmark.adapters.squiddiff_adapter import SquidiffAdapter; print('adapters OK')"
python scripts/generate_method_expansion_configs.py     # builds 36 runtime configs

# 3. (Optional) smoke-test one task before the full array:
SGE_TASK_ID=1 bash jobs/shirokane/run_method_expansion_benchmark_array.sh

# 4. Submit the full array (3 methods x 2 primaries x 6 scenarios = 36 tasks):
qsub jobs/shirokane/run_method_expansion_benchmark_array.sh

# 5. After it completes, add the three methods to the manifest + re-summarize
#    (see SUBMIT_publication.md), then regenerate figures/tables.
```

Until step 1 is done, the array job is safe to leave unsubmitted — it will not
produce misleading partial results. The pseudotime tasks (scenarios D/E/F) also
require `run_pseudotime_axis_build.sh` to have produced the `*_pseudotime_benchmark_input.h5ad`.
