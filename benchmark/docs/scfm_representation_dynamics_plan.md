# scFM Representation for Dynamics: Work Plan

## 1. Background

Current trajectory benchmark results mainly evaluate projected cells in the
original expression or PCA-derived space. This is useful, but it does not
answer whether single-cell foundation model representations are helpful for
temporal dynamics reconstruction.

This supplementary module asks:

> When different fixed embedding methods are used as the input matrix for the
> same trajectory inference model, how do the current trajectory benchmark
> metrics change?

Here, scFM embeddings are not used as cell annotation ground truth. They are
treated as alternative input representations for trajectory models. The primary
comparison is a controlled input-representation ablation: HVG-PCA versus
Geneformer, scGPT, and scFoundation embeddings under the same dataset, cells,
scenario split, trajectory model config, random seed, and metric protocol.

## 2. Scientific Objective

The main objective is to measure the effect of the input embedding method on
trajectory inference performance in the current project.

The controlled variable is the input representation:

```text
same dataset
same cells
same time/state labels
same scenario split
same trajectory model config
same random seed
same evaluation protocol
only input embedding changes
```

The primary question is:

> Relative to an HVG-PCA input baseline, do Geneformer, scGPT, or scFoundation
> embeddings improve, worsen, or leave unchanged the current trajectory metrics?

Secondary diagnostics ask whether each representation preserves temporal and
state structure before model training.

The expected contribution is not simply "we used scGPT/Geneformer". The
contribution is a controlled input-representation benchmark for single-cell
trajectory dynamics.

## 3. Core Hypotheses

### H1. scFM embeddings may improve state-level structure

Because scFM models are pretrained on large single-cell corpora, their cell
embeddings may better capture broad cell states, such as somatic, epithelial,
intermediate, XEN-like, and pluripotent states.

Expected signal:

- higher embedding coherence against frozen milestone states;
- better separation of early, intermediate, and terminal populations;
- lower classifier entropy in projected representation space.

### H2. scFM embeddings may suppress temporal dynamics

Pretrained scFM embeddings may emphasize stable cell identity and dataset-level
structure while compressing weak temporal signals. This may hurt extrapolation,
especially in time-series data with strong time-library confounding.

Expected signal:

- lower temporal signal ratio than HVG-PCA;
- weaker adjacent-time ordering;
- worse Scenario B extrapolation despite good state clustering;
- projected cells collapse toward dominant states.

## 4. Experimental Design

### 4.1 Datasets

Initial MVP dataset:

- GSE230659 marker-FM transition silver benchmark.

Recommended expansion:

- GSE178325 marker-FM transition silver benchmark.
- GSE242424 OSKM ground truth benchmark.

GSE230659 is a good first target because it is central to the current project
and has strong time-library confounding, making it useful for testing whether
scFM embeddings preserve or suppress temporal signal.

### 4.2 Benchmark Scenarios

Use the existing observed-time scenarios:

| Scenario | Role | Example interpretation |
| --- | --- | --- |
| A | observed-time interpolation / all-time evaluation | can the model preserve observed trajectory structure? |
| B | future extrapolation | can the model predict later states from early states? |
| C | mixed interpolation and extrapolation | can the model handle both missing internal and future states? |

The representation-dynamics benchmark should be reported separately from the
official expression-space benchmark.

### 4.3 Methods

Run the same trajectory models used in the current benchmark:

| Method | Include? | Notes |
| --- | --- | --- |
| MIOFlow | yes | trains the ODE directly on the active input representation |
| PRESCIENT | yes | trains/simulates in the active input representation |
| scNODE | yes | trains the VAE/latent dynamics from the active input representation |
| WOT | optional lineage-only | can be evaluated on representation-derived transition structure |
| CellRank2 | optional lineage-only | same caution as WOT |

## 5. Representation Arms

All representation arms should be reduced to a common final dimension, such as
50, before model training. This keeps the trajectory model input dimensionality
fixed across embedding methods.

| ID | Representation | Purpose |
| --- | --- | --- |
| rep_hvg_pca50 | HVG2000 -> PCA50 | expression-derived input baseline |
| rep_scgpt_cls_pca50 | scGPT cell/CLS embedding -> PCA50 | test scGPT representation |
| rep_geneformer_cls_pca50 | Geneformer cell embedding -> PCA50 | test rank/attention-based foundation representation |
| rep_scfoundation_pca50 | scFoundation cell embedding -> PCA50 | test large-scale scFM representation |

This benchmark is intentionally limited to these four arms.

## 6. Primary Experimental Mode: Embedding as Model Input

This is the required implementation for this experiment.

Workflow:

```text
raw expression
  -> fixed embedding method
       - HVG-PCA
       - Geneformer
       - scGPT
       - scFoundation
  -> common PCA/standardization to X_rep
  -> trajectory model trains on X_rep
  -> trajectory model predicts future X_rep
  -> evaluate the prediction with the current trajectory metric families
```

Primary comparison unit:

```text
same trajectory model + same scenario + same metric
HVG-PCA input vs Geneformer input vs scGPT input vs scFoundation input
```

Important interpretation:

- The trajectory model output is in the active input representation space.
- Forecast-distance metrics measure distributional prediction quality in that
  active input space.
- Embedding coherence and lineage metrics remain comparable because labels,
  reference graphs, scenario splits, and trajectory model configs are fixed.
- Do not interpret representation-space forecasts as gene-expression forecasts
  unless a separate train-only decoder is added and evaluated.

Out of scope for the primary experiment:

- using scFM only as a frozen evaluator of gene-expression outputs;
- reusing old expression-space formal results as the primary controlled
  baseline;
- mixing representation-space rankings into the official expression-space
  benchmark tables.

## 7. Implementation Plan

### Step 0. Reference-first implementation gate

The embedding extraction code must be grounded in the original paper and the
official implementation reference for each model. Do not invent preprocessing,
tokenization, model-loading, pooling, or embedding-layer logic from memory.

Before implementing or modifying any scFM extractor, record the exact paper,
code repository, checkpoint, tutorial/example file, and commit or release used
in `adata.uns["scfm_embedding_metadata"]` and in the generated
`representation_metadata.json`.

Required references:

| Model | Paper | Implementation reference |
| --- | --- | --- |
| Geneformer | [Transfer learning enables predictions in network biology](https://www.nature.com/articles/s41586-023-06139-9) | The original paper's code availability points to the [ctheodoris/Geneformer Hugging Face repository](https://huggingface.co/ctheodoris/Geneformer), including the pretrained model, transcriptome tokenizer, and code. Treat this as the authoritative implementation source; do not substitute an unofficial GitHub recreation without explicit justification and metadata. |
| scGPT | [scGPT: toward building a foundation model for single-cell multi-omics using generative AI](https://www.nature.com/articles/s41592-024-02201-0) | [bowang-lab/scGPT](https://github.com/bowang-lab/scGPT), the GitHub repository linked by the paper's code availability statement. |
| scFoundation | [Large-scale foundation model on single-cell transcriptomics](https://www.nature.com/articles/s41592-024-02305-7) | [biomap-research/scFoundation](https://github.com/biomap-research/scFoundation), the GitHub repository linked by the paper's code availability statement. |

Implementation rules:

1. Reuse the official model package, tokenizer, vocabulary, and inference
   utilities whenever available.
2. If local wrapper code is needed, keep it thin and cite the exact upstream
   function, notebook, or README section it mirrors.
3. If the reference paper and the upstream repository disagree, stop and record
   the discrepancy before running benchmarks.
4. Every embedding artifact must record `paper_reference`, `code_reference`,
   `upstream_commit_or_release`, `checkpoint`, `tokenizer_or_gene_vocab`,
   `preprocessing`, `embedding_layer`, and `pooling_or_cell_embedding_method`.

### Step 1. Extract frozen scFM embeddings

Create a script such as:

```text
benchmark/representations/extract_scfm_embeddings.py
```

Expected output:

```text
adata.obsm["X_scgpt_cls"]
adata.obsm["X_geneformer_cls"]
adata.obsm["X_scfoundation"]
adata.uns["scfm_embedding_metadata"]
```

Metadata should record:

```text
model_name
model_checkpoint
gene_vocabulary
input_gene_count
matched_gene_count
preprocessing
embedding_layer
embedding_dim
date_created
code_version
paper_reference
code_reference
upstream_commit_or_release
tokenizer_or_gene_vocab
pooling_or_cell_embedding_method
```

### Step 2. Build common 50-dimensional representation inputs

Create:

```text
benchmark/representations/build_representation_inputs.py
```

For the HVG-PCA baseline, follow the reference-paper style preprocessing:

```text
raw counts
  -> library-size normalization
  -> log1p transform
  -> select 2,000 highly variable genes
  -> PCA to 50 dimensions
  -> store as adata.obsm["X_rep"]
```

For each scFM representation:

```text
raw representation
  -> StandardScaler fitted on training cells only
  -> PCA fitted on training cells only
  -> transform all cells
  -> store as adata.obsm["X_rep"]
```

Important anti-leakage rule:

> The scaler and PCA reducer must be fit only on training cells for each
> scenario, then applied to held-out cells.

If scenario-specific reducers are too complex for the first implementation,
start with an explicitly labeled exploratory version:

```text
reducer_fit_scope: all_cells_exploratory
```

Do not use all-cell reducers for formal supplementary rankings.

### Step 3. Add representation-aware configs

Example config block:

```yaml
representation:
  enabled: true
  input_mode: obsm
  obsm_key: X_rep
  representation_id: rep_scgpt_cls_pca50
  input_space: scfm_embedding
  output_space: representation
  final_dim: 50
  reducer_fit_scope: train_only
```

Suggested runtime config names:

```text
scnode_gse230659_rep_scgpt_cls_A.yaml
mioflow_gse230659_rep_scgpt_cls_A.yaml
prescient_gse230659_rep_scgpt_cls_A.yaml
```

### Step 4. Modify method runners minimally

Add a shared helper:

```python
def get_model_input_matrix(adata, cfg):
    rep = cfg.get("representation", {})
    if rep.get("enabled") and rep.get("input_mode") == "obsm":
        return np.asarray(adata.obsm[rep["obsm_key"]], dtype=np.float32)
    return adata.X
```

Method-specific interpretation:

- MIOFlow: train ODE directly in `X_rep` space.
- PRESCIENT: use `X_rep` as the PRESCIENT simulation input space.
- scNODE: use `X_rep` as input to its VAE; output is an active-input-space
  forecast.

### Step 5. Apply current trajectory metrics to the active input space

The experiment should reuse the current trajectory metric families wherever
possible. The difference is that, for representation-input runs, the model
output and the observed target are both in `X_rep` space.

Create:

```text
benchmark/evaluation/eval_representation_forecast.py
```

This script is a representation-aware adapter for the current forecast metric
family. It should mirror the existing forecast-distance definitions while
making the active input space explicit in names and metadata.

Inputs:

```text
projected_representation.npy
observed_representation.npy
cell/time metadata
representation_metadata.json
```

Core metrics:

```text
rep_wasserstein_distance
rep_gaussian_mmd
rep_energy_distance_mmd
rep_hausdorff_loss
```

These correspond to the current forecast-distance metrics, evaluated on the
common 50-dimensional input representation rather than on gene-expression
features. They are intended for controlled comparisons across embedding
methods because all arms use the same cells, splits, model configs, seeds, and
final dimensionality.

### Step 6. Add temporal-signal diagnostics

Create:

```text
benchmark/evaluation/eval_representation_temporal_signal.py
```

Recommended diagnostics:

```text
temporal_signal_ratio
adjacent_time_centroid_distance
non_adjacent_time_centroid_distance
time_prediction_macro_f1
time_prediction_r2
state_centroid_separation
state_silhouette_score
forward_mass_fraction
backward_mass_fraction
forbidden_edge_mass
```

These diagnostics answer whether the representation itself retains dynamic
signal before any trajectory model is trained.

### Step 7. Generate separate reports

Do not mix representation results into the official silver report.

Suggested output directory:

```text
benchmark/reports/representation_dynamics/
```

Suggested report files:

```text
representation_method_summary.md
representation_forecast_rankings.csv
representation_temporal_signal.csv
representation_lineage_metrics.csv
representation_embedding_coherence.csv
representation_metadata_manifest.csv
```

## 8. Metrics

The metric families remain the current trajectory benchmark families:

- forecast accuracy / distributional prediction quality;
- embedding coherence / state separability;
- lineage fidelity / transition-graph consistency;
- temporal-signal diagnostics.

For representation-input runs, forecast-distance metrics are computed in the
active input representation space. This is the correct space for the primary
ablation because the trajectory model was trained and simulated in that space.

### 8.1 Representation-space forecast metrics

Compare predicted representation distributions with observed held-out
representation distributions at each target time point.

| Metric | Direction | Interpretation |
| --- | --- | --- |
| rep_wasserstein_distance | lower better | average marginal distribution shift |
| rep_gaussian_mmd | lower better | kernel distribution mismatch |
| rep_energy_distance_mmd | lower better | energy-distance mismatch |
| rep_hausdorff_loss | lower better | worst-case nearest-neighbor mismatch |

### 8.2 Temporal-signal metrics

These are computed on observed cells before model training.

| Metric | Direction | Interpretation |
| --- | --- | --- |
| temporal_signal_ratio | higher usually better | fraction of representation variance explained by timepoint centroids |
| adjacent_time_centroid_distance | context-dependent | whether neighboring timepoints are distinguishable |
| time_prediction_macro_f1 | higher better | whether timepoint is recoverable from representation |
| state_silhouette_score | higher better | whether state labels are separated |

The key comparison is not simply "higher is always better". A useful
representation should preserve both temporal order and biological state
structure without collapsing all cells into batch or endpoint classes.

### 8.3 Lineage-oriented representation metrics

If projected representations are assigned to state centroids, compute:

```text
forward_mass_fraction
backward_mass_fraction
forbidden_edge_mass
transition_matrix_cosine
transition_matrix_spearman
edge_weight_auprc
```

These metrics test whether the representation-space dynamics move in the
expected biological direction.

## 9. Expected Results

### Expected positive result

scFM embeddings improve representation-space forecast and state coherence:

```text
rep_scgpt_cls_pca50 < rep_hvg_pca50 in forecast loss
rep_geneformer_cls_pca50 or rep_scfoundation_pca50 improves state coherence
temporal_signal_ratio remains comparable to HVG-PCA
```

Interpretation:

> The foundation model embedding provides a useful state-aware geometry for
> trajectory dynamics.

### Expected mixed result

scFM embeddings improve state separation but worsen extrapolation:

```text
state_silhouette_score improves
Scenario B forecast loss worsens
temporal_signal_ratio decreases
```

Interpretation:

> The foundation model preserves broad cellular identity but compresses or
> distorts weak temporal signals needed for dynamic extrapolation.

This is still a valuable result.

### Expected negative result

HVG-PCA remains better than scFM embeddings for temporal dynamics:

```text
rep_hvg_pca50 performs best in A/B/C
scFM arms show lower temporal_signal_ratio
timepoint ordering is weaker in scFM spaces
```

Interpretation:

> Pretrained scFM embeddings are not automatically suitable for temporal
> dynamics reconstruction; local expression variation may preserve trajectory
> signal better than generic pretrained representations.

This is publishable as a benchmark finding if carefully controlled.

## 10. Risk Control

### Risk 1. Confusing active-input-space forecast with gene-expression forecast

Mitigation:

- report representation-input results separately from official expression-space
  rankings;
- use metric names prefixed with `rep_`;
- state clearly that representation-input runs do not predict gene expression
  unless a decoder is added.

### Risk 2. Leakage from fitting reducers on all cells

Mitigation:

- fit scaler/PCA on training cells only;
- record `reducer_fit_scope`;
- keep all-cell reducer runs as exploratory only.

### Risk 3. scFM embedding reflects batch or dataset identity

Mitigation:

- report timepoint and batch predictability;
- compare to the HVG-PCA baseline;
- evaluate across GSE178325 and GSE230659 when possible.

### Risk 4. scFM model is also used in label construction

Mitigation:

- do not use the same scFM embedding as official ground truth;
- keep `official_silver` labels independent;
- label all representation-dynamics results as a separate input-ablation
  analysis.

## 11. Minimal Viable Experiment

MVP scope:

```text
Dataset: GSE230659
Scenarios: A, B, C
Methods: MIOFlow, PRESCIENT, scNODE
Representations:
  - rep_hvg_pca50
  - rep_scgpt_cls_pca50
  - rep_geneformer_cls_pca50
  - rep_scfoundation_pca50
Metrics:
  - forecast accuracy / distributional prediction loss in X_rep space
  - embedding coherence / state separability
  - lineage fidelity / transition consistency
  - temporal_signal_ratio
  - forward/backward/forbidden transition mass
```

Success criteria:

1. All representation arms can run through the same scenario split.
2. HVG-PCA is run through the same representation-input pipeline as the scFM
   arms; old expression-space formal results are not used as the primary
   controlled baseline.
3. Results are reported separately from official expression-space rankings.
4. Each representation has complete metadata.
5. The analysis can answer whether Geneformer, scGPT, or scFoundation improves,
   worsens, or leaves unchanged dynamic forecasting relative to HVG-PCA.

## 12. Full Experiment

Expanded scope:

```text
Datasets:
  - GSE230659
  - GSE178325
  - GSE242424

Methods:
  - MIOFlow
  - PRESCIENT
  - scNODE
  - WOT lineage-only
  - CellRank2 lineage-only

Representations:
  - HVG-PCA50
  - scGPT CLS PCA50
  - Geneformer PCA50
  - scFoundation PCA50
```

Final report questions:

1. Which representation best supports interpolation?
2. Which representation best supports extrapolation?
3. Do scFM embeddings improve state coherence at the cost of temporal signal?
4. Are effects consistent across chemical reprogramming and OSKM reprogramming?
5. As a secondary analysis only, do representation-input results correlate with
   historical gene-expression forecast results?

## 13. Suggested Manuscript Framing

Suggested wording:

> We performed a representation-dynamics ablation to test whether pretrained
> single-cell foundation model embeddings provide better input features for
> trajectory inference than HVG-PCA. The foundation model embeddings were
> frozen, were not used to define the official silver labels, and were evaluated
> with the current trajectory metric families under matched scenario splits,
> model configs, random seeds, and final input dimensionality. This analysis
> distinguishes static cell-state separability from dynamic temporal
> predictability.

Avoid:

> scFM embeddings are biological ground truth.

Avoid:

> scFM embeddings replace expression-space forecast metrics.

## 14. Reference Papers

- [scGPT: toward building a foundation model for single-cell multi-omics using generative AI](https://www.nature.com/articles/s41592-024-02201-0)
- [Geneformer: Transfer learning enables predictions in network biology](https://www.nature.com/articles/s41586-023-06139-9)
- [scFoundation: Large-scale foundation model on single-cell transcriptomics](https://www.nature.com/articles/s41592-024-02305-7)
- [TIGON: Reconstructing growth and dynamic trajectories from single-cell transcriptomics data](https://www.nature.com/articles/s42256-023-00763-w)
- [Benchmarking zero-shot single-cell foundation model embeddings for cellular dynamics reconstruction](https://gist.science/fr/paper/bio/10.64898/2026.03.10.710748)
