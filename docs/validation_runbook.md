# Validation runbook — completing the annotation-robustness story

The first server run (`docs/scanvi_run_results_interpretation.md`) showed a
confident cross-protocol collapse. Three follow-ups complete the manuscript; all
reuse the existing `src/` code + Shirokane jobs. Run on Shirokane in the `scvi`
env (`export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH`).

## A. Geometric OOD on the real scANVI latent (cheap, do first)
Re-scores the completed run's model — no retraining:
```bash
python scripts/score_ood.py --config config.yaml --profile server
# -> results/annotation_branch/<query>/ood_scanvi_latent.csv + ood_summary_scanvi_latent.json
```
Expectation: a large OOD fraction (the local PCA-latent proxy gave 39–44%),
confirming scANVI confidence (≤0.3%) is miscalibrated. This is Result 3.3 on the
true latent.

## B. Positive control — GSE175634 (decisive)
GSE175634 (iPSC→cardiomyocyte) HAS author `type` labels (frozen provider
`gse175634_cardiac_silver_v1`). Build a labelled full-gene input, split train/test,
train scANVI on train, map test, and check recovery + agreement:
1. Prep a labelled input for GSE175634 with `scripts/prep_scanvi_inputs_reconstruct.py`
   pattern (full-gene matrix + author labels by barcode; freeze 2000 HVG).
2. Split obs into train/test (e.g. 80/20 stratified by `type`).
3. Train scANVI on train (`scripts/train_scanvi_reference.py` with a GSE175634 profile);
   map test; recovery + per-state recall should be high if the method works when
   systems match. **If yes, it isolates the GSE242424→chemical failure to the
   OSKM↔chemical gap.**

## C. Same-protocol transfer — GSE230659 <-> GSE178325
Use one chemical dataset's `final_milestone_label_coarse` as the reference label
for the other (both chemical, shared state space). Add a config profile pointing
reference=one chemical dataset, query=the other, `reference_label_key:
final_milestone_label_coarse`, and run train→map→compare. A within-protocol
transfer should agree far better (high NMI/correspondence) and gives the
meaningful annotation-strategy comparison for the trajectories (Figs 3–4).

After A–C, re-render the figures and finalise `manuscript/manuscript_draft.md`
(replace the §3.4 "pending" markers with the results).

---

## One-command submission (all .sh files are ready)

On Shirokane, from the repo root (the `scvi` env + `LD_LIBRARY_PATH` fix are
handled inside each job):

```bash
bash jobs/shirokane/submit_scanvi_validation.sh
```

This submits four independent experiments:

| Job script | What it does | Output |
|---|---|---|
| `run_scanvi_score_ood.sh` | geometric OOD on the **completed** scANVI latent (no retrain) | `results/annotation_branch/<query>/ood_scanvi_latent.csv`, `ood_summary_scanvi_latent.json` |
| `run_scanvi_positive_control.sh` | scANVI self-transfer on GSE175634 (labels known) | `results/positive_control_gse175634/positive_control_summary.json`, `confusion.csv` |
| `run_scanvi_train_reference.sh` + `run_scanvi_map.sh` (`CONFIG=config_same_protocol_230to178.yaml`) | same-protocol GSE230659→GSE178325 | `results/same_protocol_230to178/summary.json` + figures |
| `run_scanvi_train_reference.sh` + `run_scanvi_map.sh` (`CONFIG=config_same_protocol_178to230.yaml`) | same-protocol GSE178325→GSE230659 | `results/same_protocol_178to230/summary.json` + figures |

### Or submit individually
```bash
# A — cheapest, do first (reuses the finished model):
qsub jobs/shirokane/run_scanvi_score_ood.sh

# B — positive control (decisive):
qsub jobs/shirokane/run_scanvi_positive_control.sh

# C — same-protocol 230 -> 178 (train, then map holding on it):
T=$(qsub -v CONFIG=config_same_protocol_230to178.yaml,PROFILE=server jobs/shirokane/run_scanvi_train_reference.sh | awk '/Your job/{print $3}')
qsub -hold_jid "$T" -v CONFIG=config_same_protocol_230to178.yaml,PROFILE=server jobs/shirokane/run_scanvi_map.sh

# D — same-protocol 178 -> 230:
T=$(qsub -v CONFIG=config_same_protocol_178to230.yaml,PROFILE=server jobs/shirokane/run_scanvi_train_reference.sh | awk '/Your job/{print $3}')
qsub -hold_jid "$T" -v CONFIG=config_same_protocol_178to230.yaml,PROFILE=server jobs/shirokane/run_scanvi_map.sh
```

### What each result tells the paper
- **A (OOD on real latent):** confirms Result 3.3 — expect a large OOD fraction (local proxy gave 39–44%), versus scANVI's ≤0.3%.
- **B (positive control):** high `test_accuracy` / per-state recall ⇒ the method works when systems match ⇒ the GSE242424→chemical failure is the cross-protocol gap (the manuscript's decisive control).
- **C/D (same-protocol):** high NMI / correspondence ⇒ within-protocol transfer is trustworthy and gives the meaningful annotation-strategy comparison for the trajectories (Figs 3–4 "what works").

After these finish, sync the outputs and I will drop the numbers into
`manuscript/manuscript_draft.md` §3.4 and finalise the figures.
