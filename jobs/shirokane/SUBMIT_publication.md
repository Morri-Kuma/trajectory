# SUBMIT_publication.md — Shirokane runbook for the remaining publication items

Four compute items remain to take the benchmark from "publication-ready modulo compute"
to fully complete. Each is independent; run any subset. All paths are relative to
`/home/xzy0723/projects/trajectory`. **Sync the repo to Shirokane first** (the new
scripts + jobs below must be present).

New files this round:
- `scripts/generate_multiseed_configs.py`, `jobs/shirokane/run_multiseed_benchmark_array.sh`, `scripts/aggregate_multiseed_ci.py`
- `scripts/compute_spearman_baseline.py`, `jobs/shirokane/run_spearman_baseline.sh`
- `benchmark/methods/VENDORING.md`
- (existing) `jobs/shirokane/run_gse175634_benchmark_array.sh`, `jobs/shirokane/run_method_expansion_benchmark_array.sh`

---

## 1. Multi-seed confidence intervals  (GPU; ~90 short runs)

Gives per-method mean ± 95% CI on the primary datasets — closes the manuscript's main
statistical caveat.

```bash
cd /home/xzy0723/projects/trajectory
python scripts/generate_multiseed_configs.py          # writes 90 configs under configs/runtime/multiseed/
qsub jobs/shirokane/run_multiseed_benchmark_array.sh   # 3 methods x 2 datasets x 3 scen x 5 seeds
# after the array finishes:
python scripts/aggregate_multiseed_ci.py               # -> benchmark/reports/official_silver/multiseed_ci.csv
```

## 2. Spearman correlation baseline  (CPU; one short job)

The scTimeBench null that methods must beat. Machinery already exists; this just feeds
it the AnnData. No retraining.

```bash
qsub jobs/shirokane/run_spearman_baseline.sh           # -> benchmark/reports/official_silver/spearman_baseline.csv
```

## 3. GSE175634 cardiac lineage  (GPU; existing job)

The cardiac provider (`gse175634_cardiac_silver_v1`) and lineage-enabled configs already
exist. Run the benchmark array (if not already complete) to fill in the cardiac lineage
metrics, then it can be folded into the cross-system tables.

```bash
# preprocess only if the input/provider are not yet built:
# qsub jobs/shirokane/run_gse175634_preprocess.sh
qsub jobs/shirokane/run_gse175634_benchmark_array.sh   # scnode/mioflow/prescient A/B/C + wot/cellrank2 A
```

## 4. Three new models — scIMF / PI-SDE / Squidiff  (GPU; needs vendoring)

Follow `benchmark/methods/VENDORING.md` to drop in each model's upstream source and wire
its `run.py`, then:

```bash
python scripts/generate_method_expansion_configs.py    # 36 configs (idempotent)
qsub jobs/shirokane/run_method_expansion_benchmark_array.sh
```

---

## After any run — fold results back in (local or Shirokane)

1. **Sync results back** to your laptop (the `benchmark/results/**` dirs).
2. **Add new datasets/methods to the manifest** `benchmark/results/result_manifest.yaml`
   (multi-seed and baseline do NOT need manifest edits; the 3 new methods and cardiac do —
   add result_dirs entries mirroring the existing blocks).
3. **Re-summarize**: `python benchmark/evaluation/summarize_official_silver.py`
4. **Regenerate figures + stats**:
   ```bash
   python scripts/make_manuscript_figures.py        # fig1-3
   python scripts/make_fig4_cross_system.py         # fig4
   python scripts/rank_concordance_analysis.py      # fig7 + rank_concordance.csv
   ```
5. **Update the manuscript** (Table 4 / §4.5 with CIs and the baseline row) and recompile.
   NOTE: the workspace mount truncates `manuscript/benchmark_manuscript.tex` past ~27 KB;
   compile by copying the full file to a non-mount path first (see `manuscript/REPRODUCIBILITY.md`).

## Suggested order / dependencies

Items 1–4 are independent and can be submitted together. Only the 3-new-model array (4)
has a prerequisite (vendoring). Multi-seed (1) is the highest-value for publication
(per-method CIs); the Spearman baseline (2) is the cheapest (CPU, minutes).
