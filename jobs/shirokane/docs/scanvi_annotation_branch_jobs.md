# Shirokane jobs — scANVI annotation branch

Three dependent stages that run the full server-side scANVI annotation-robustness
analysis (GSE242424 reference → GSE178325 / GSE230659 queries). All paths and
resources are overridable via environment variables; defaults assume the
`/home/xzy0723/projects/trajectory` checkout.

## Stages
| # | Script | Queue | What it does |
|---|---|---|---|
| 1 | `run_scanvi_prep_inputs.sh` | CPU (`s_vmem=96G`) | `src/data/build_scanvi_inputs.py`: join full-gene **counts** + per-cell labels by barcode → 3 shared-gene-space `.h5ad` in `results/annotation_branch/inputs/`. Needed because the HVG2000 inputs are gene-disjoint (ref∩query ≈ 181–238). |
| 2 | `run_scanvi_train_reference.sh` | GPU (`h100=1`) | `scripts/train_scanvi_reference.py`: train scVI→scANVI on the prepped reference (`author_cluster_label`), assert raw counts, enforce the ≥0.90 recovery gate, save model to `models/scanvi_reference_gse242424/`. |
| 3 | `run_scanvi_map_compare.sh` | GPU (`h100=1`) | scArches-map the queries onto the frozen model and run `src.pipeline` in `server` mode → annotation comparison + trajectory + figures in `results/annotation_branch/`. |

## One-command submission (dependency chain)
```bash
bash jobs/shirokane/submit_scanvi_annotation_branch.sh
# prep (CPU) -> train (GPU, holds on prep) -> map+compare (GPU, holds on train)
```

## Before first run — verify / set
- **Conda envs:** `TRAJ_CONDA_ENV` (scanpy/anndata, default `traj_env`) for prep;
  `SCVI_CONDA_ENV` (scvi-tools + a CUDA torch, default `scvi`) for train/map.
- **Raw counts:** every `*_FULLGENE` input must hold raw integer counts — the
  train step hard-fails otherwise. If your full-gene `.h5ad` is normalized, point
  `REF_FULLGENE` / `Q*_FULLGENE` at a counts source or pass `--counts-layer`.
- **Barcode join:** each `*_FULLGENE` and its `*_LABELS` file must share
  `obs_names`. The prep step errors if there is no overlap.
- **Label keys:** reference = `author_cluster_label` (set in `config.yaml`);
  queries = `final_milestone_label_coarse`.
- **GPU flag:** scripts use `-l h100=1`; change to `a100=1` / `v100=1` for your queue.

## Override examples
```bash
# different checkout + A100 queue + custom reference full-gene source
TRAJ_PROJECT_ROOT=/home/me/trajectory \
TRAIN_RESOURCES="-l s_vmem=128G -l a100=1" MAP_RESOURCES="-l s_vmem=128G -l a100=1" \
REF_FULLGENE=data/processed/gse242424_human/GSE242424_raw_full_gene_benchmark_input.h5ad \
bash jobs/shirokane/submit_scanvi_annotation_branch.sh
```

## Outputs
- `results/annotation_branch/inputs/*_labelled.h5ad` (stage 1)
- `models/scanvi_reference_gse242424/` + `reference_genes.txt` (stage 2)
- `results/annotation_branch/{gse178325,gse230659}/` tables + figures, `summary.json` (stage 3)

See `scanvi_annotation_branch.md` (design) and
`results/test_outputs/SERVER_DRYRUN_NOTES.md` (why full-gene inputs).
