# Projected HVG2000 Embedding Annotation Plan

## Purpose

This note defines the implementation plan for formal projected-cell milestone annotation after the metrics-system revision.

The key constraint is that `projected_expression.npy` must remain the method-generated expression output. It should not be overwritten, normalized in place, replaced by embeddings, or reinterpreted as another representation. Any embedding used for annotation must be generated as a separate sidecar artifact with explicit metadata.

## Current State

GSE230659 observed cells already have the real-sensitivity milestone annotation system:

- `milestone_marker_label`
- `milestone_embedding_label`
- `milestone_classifier_label`
- `consensus_milestone_label`

The formal projected-cell annotation policy now allows only two automatic providers:

- `embedding_based`
- `classifier_based`

Marker-score projection is deprecated for formal projected-cell annotation. Marker scores may be retained only as diagnostics.

For scNODE, the current projected outputs are:

- `projected_expression.npy`: generated projected-cell expression, shaped as timepoint by simulated cell by gene.
- `projected_embedding.npy`: scNODE internal latent representation, currently 50-dimensional because `latent_dim: 50`.

The current observed-cell embedding milestone classifier expects 512-dimensional `X_scGPT` features. Therefore, scNODE's 50-dimensional `projected_embedding.npy` is not compatible with the official embedding milestone model and must not be treated as `milestone_embedding_label` input.

## Formal Route Decision

For GSE230659, the formal projected-cell Embedding Coherence annotation route is
now fixed as an HVG2000-consistent route.

This is a project-level rule, not a temporary implementation preference:

- Formal Embedding Coherence uses the HVG2000 gene universe.
- Full-gene scGPT embeddings are allowed only as observed-cell reference,
  provenance evidence, or sensitivity analysis.
- A full-gene-trained scGPT embedding classifier must not be used to annotate
  HVG2000-only projected expression for a formal projected-cell benchmark.
- All formal projected labels must come from providers trained and applied in
  the same gene universe. If a provider uses an internal latent representation,
  that representation must be derived from the same HVG2000 expression contract
  for both observed and projected cells.

The practical consequence is that the current full-gene scGPT embedding model is
not a formal projected-cell provider for scNODE GSE230659 outputs. A new
HVG2000-consistent embedding provider must be registered before complete formal
two-provider projected annotation can be produced.

## Design Principle

The formal route must preserve a strict separation:

```text
projected_expression.npy
  = method-generated expression; primary projected-cell output; immutable input

projected_hvg_embedding.npy
  = sidecar representation derived from projected_expression.npy under the
    same HVG2000 provider contract used for observed cells

projected_milestone_labels.csv
  = annotation table derived from HVG2000-compatible embedding and classifier
    predictions
```

This means the expression matrix semantics do not change. The sidecar embedding
is a downstream representation, not a replacement for the expression output.

## Required Artifacts

For each projected-cell result directory, the target complete formal annotation set is:

- `projected_expression.npy`
- `projected_hvg_embedding.npy`
- `projected_hvg_embedding_metadata.json`
- `projected_milestone_labels.csv`
- `projected_milestone_label_counts.csv`
- `projected_milestone_annotation_metadata.json`

The existing method-native `projected_embedding.npy` may remain, but it must be described as method latent space, not as scGPT or formal milestone-provider input.

## Implementation Steps

### 1. Freeze Input Contract

For each run directory, identify:

- result directory
- `projected_expression.npy`
- `run_metadata.json`
- reference h5ad path from the official GSE230659 real-sensitivity milestone input
- gene order used by `projected_expression.npy`
- target timepoints represented by the first axis of `projected_expression.npy`

The implementation must verify that the gene axis of `projected_expression.npy` matches the benchmark input h5ad var order. If it does not match, the script must stop unless an explicit, metadata-recorded gene-alignment file is provided.

### 2. Build Projected AnnData Without Mutating Expression

Create a temporary or sidecar AnnData for projected cells:

- flatten `projected_expression.npy` from `(n_timepoints, n_cells, n_genes)` to `(n_projected_cells, n_genes)`
- assign stable projected cell ids
- attach projected timepoint metadata
- use the official benchmark input var names and var metadata

The raw array in `projected_expression.npy` must not be rewritten.

### 3. Reuse the Observed HVG2000 Provider Contract

Before generating `projected_hvg_embedding.npy`, inspect and document the
observed-cell route used to train the selected HVG2000 embedding provider.

The projected route must match the observed route:

- same gene vocabulary policy
- same gene-symbol casing policy
- same normalization/log transform/binning policy
- same embedding algorithm and model asset
- same feature ordering and feature count

The current full-gene scGPT `X_scGPT` route does not satisfy this requirement
for HVG2000-only projected expression, so it is not the formal route for
projected-cell Embedding Coherence.

### 4. Generate Sidecar HVG2000 Embeddings

Run the projected AnnData through the selected HVG2000-compatible embedding
provider, producing:

- `projected_hvg_embedding.npy`
- `projected_hvg_embedding_metadata.json`

The metadata must include:

- source expression file path
- source expression shape
- flattened projected-cell count
- gene count
- embedding provider identifier
- preprocessing parameters
- embedding dimension
- confirmation that the route is HVG2000-consistent
- any missing or padded genes
- creation timestamp
- code/script name and relevant git status if available

### 5. Run Embedding-Based Provider

Use the registered HVG2000-compatible embedding milestone model on
`projected_hvg_embedding.npy`.

Output columns:

- `milestone_embedding_label`
- `milestone_embedding_confidence`
- `embedding_prob_<milestone>`

The script must verify that the projected sidecar embedding dimension equals
the model's expected feature count. If not, formal complete annotation must
fail.

### 6. Run Classifier-Based Provider

Use `projected_expression.npy` as the expression input for the classifier-based provider.

Output columns:

- `milestone_classifier_label`
- `milestone_classifier_confidence`
- `classifier_prob_<milestone>`

This route may align genes to the CellTypist/classifier feature order, but that alignment must be performed in memory or into a derived temporary matrix. It must not modify `projected_expression.npy`.

### 7. Build Consensus

Build `consensus_milestone_label` from the two provider outputs:

- if labels agree, use the shared label
- if labels disagree, mark `ambiguous` unless a predeclared confidence rule is used
- write `consensus_status` and `consensus_confidence`

The exact consensus rule must be written into `projected_milestone_annotation_metadata.json`.

### 8. Validation

Update or run projected-label validation to enforce:

- `formal_two_provider_complete: true` only when both providers are available
- no `milestone_embedding_label == unavailable` in formal complete outputs
- `projected_expression.npy` row count after flattening equals projected label count
- `projected_hvg_embedding.npy` row count equals projected label count
- probability columns are present and numerically valid
- marker-score diagnostics, if present, use only `diag_marker_*` columns
- formal label columns are not populated from marker-score projection

Degraded classifier-only outputs may exist for debugging, but they must be marked:

```json
{
  "formal_two_provider_complete": false,
  "annotation_status": "degraded_classifier_only"
}
```

They must not be used for formal Embedding Coherence summaries.

### 9. Re-run Embedding Coherence Only After Complete Annotation

Only after complete projected annotation exists should GSE230659 Embedding Coherence be rerun.

Expected run order:

1. Generate `projected_hvg_embedding.npy`.
2. Generate `projected_milestone_labels.csv`.
3. Validate projected labels.
4. Run milestone Embedding Coherence.
5. Regenerate embedding summary.
6. Update report text.

Forecast Accuracy does not need to be rerun unless `projected_expression.npy` changes, which this plan explicitly forbids.

Lineage Fidelity on observed/reference transition methods is not affected by this projected-cell annotation step.

## Non-Goals

Do not:

- overwrite `projected_expression.npy`
- rename method-native `projected_embedding.npy` to scGPT embedding
- feed 50-dimensional scNODE latent vectors into the 512-dimensional scGPT milestone classifier
- copy marker-score labels into formal projected label columns
- treat classifier-only fallback as a complete formal projected annotation result
- mix full-gene observed scGPT embeddings with HVG2000 projected expression in
  a formal projected-cell provider

## Selected HVG2000 Embedding Provider

The first formal HVG2000-compatible embedding provider is:

```text
observed HVG2000 expression
-> PCA latent embedding
-> scaled logistic-regression milestone classifier
-> milestone_embedding_label
```

Training entry point:

```text
benchmark/annotation/build_hvg_pca_embedding_milestone_model.py
```

Shirokane run script:

```text
run_gse230659_hvg_pca_embedding_provider.sh
```

Expected model directory:

```text
benchmark/annotation_runs/gse230659_real_sensitivity/embedding_model_pca_hvg2000/
```

Expected model asset:

```text
hvg_pca_embedding_milestone_model.joblib
```

The saved joblib asset stores the PCA transformer, post-PCA scaler,
logistic-regression classifier, input HVG gene order, classes, and provider
metadata.  Downstream projected annotation must align `projected_expression.npy`
to that stored gene order before transforming expression to
`projected_hvg_embedding.npy`.

## Open Questions

Before implementation continues, resolve:

- Should the sidecar embedding file be named `projected_hvg_embedding.npy` or
  should provider-specific filenames be required?
- Should legacy `projected_scgpt_embedding.npy` be reserved only for sensitivity
  runs where observed and projected routes both use the same scGPT contract?

## Provenance Audit: 2026-05-10

The first audit question has been resolved for GSE230659.

Observed `X_scGPT` was generated from a full-gene scGPT route, not from the
HVG2000 expression matrix.

Evidence:

- `benchmark/results/scgpt/full/full_embed_metadata.json` records:
  - `script`: `benchmark/scgpt/full_embed.py`
  - `input_h5ad`: `data/processed/20260329_1343_GSE230659_raw.h5ad`
  - `model_dir`: `models/scgpt_whole_human`
  - `gene_col`: `index`
  - `n_genes_raw`: 27,267
  - `n_genes_in_vocab`: 23,112
  - `embedding_shape`: `[75194, 512]`
- `benchmark/results/scgpt/full/adata_scgpt_annotated.h5ad` has:
  - 75,194 cells
  - 23,112 genes
  - `obsm["X_scGPT"]` with shape `(75194, 512)`
- `benchmark/inputs/gse230659_scgpt_hvg2000/GSE230659_scGPT_annotated_HVG2000_benchmark_input.h5ad` has:
  - 75,194 cells
  - 2,000 genes
  - retained `obsm["X_scGPT"]` with shape `(75194, 512)`
  - `uns["benchmark_hvg_source"]` pointing to `benchmark/results/scgpt/full/adata_scgpt_annotated.h5ad`
  - `uns["benchmark_hvg_selection_method"]` equal to top-gene dispersion selection from the scGPT-annotated full-gene AnnData
- The real-sensitivity annotation job uses:
  - `benchmark/inputs/gse230659_milestone_fullgene_labels_hvg2000/GSE230659_milestone_fullgene_labels_HVG2000_benchmark_input.h5ad`
  - `--embedding-key X_scGPT`
  - this HVG2000 h5ad retains full-route `X_scGPT`

Therefore the current official `gse230659_milestone_embedding_v1` model is:

```text
full-gene observed expression
-> scGPT whole-human encoder
-> 512D X_scGPT
-> logistic regression milestone classifier
```

It is not:

```text
observed HVG2000 expression
-> scGPT encoder
-> 512D X_scGPT
```

This matters because current projected cells from scNODE are available as
HVG2000 expression only. Directly embedding those HVG2000 projected cells and
feeding them to the existing full-route embedding classifier would mix input
semantics unless a comparability experiment justifies it.

## Model Asset Audit: 2026-05-10

The current GSE230659 real-sensitivity model assets have different projected
annotation eligibility.

### Classifier Provider

Asset:

```text
benchmark/annotation_runs/gse230659_real_sensitivity/
  classifier_model_celltypist/classifier_milestone_model_celltypist.joblib
```

Metadata:

- `provider_family`: `classifier_based`
- `backend`: `celltypist`
- `input_h5ad`: `benchmark/inputs/gse230659_milestone_real_sensitivity_hvg2000/GSE230659_embedding_based_milestone_HVG2000_benchmark_input.h5ad`
- model features: 1,998 genes
- real-sensitivity reference h5ad: 2,000 genes
- all model features are present in the real-sensitivity HVG2000 h5ad
- two h5ad genes are not used by the CellTypist model:
  - `ERVV-1`
  - `RP1-63G5.7`

Decision:

```text
CellTypist/classifier provider is HVG2000-compatible and can continue to be
used for formal projected-cell annotation, with explicit in-memory gene
alignment to the model feature order.
```

### Embedding Provider

Asset:

```text
benchmark/annotation_runs/gse230659_real_sensitivity/
  embedding_model/embedding_milestone_model.joblib
```

Metadata:

- `provider_family`: `embedding_based`
- `method`: `scGPT_embedding_logistic_regression`
- `input_h5ad`: `benchmark/inputs/gse230659_milestone_fullgene_labels_hvg2000/GSE230659_milestone_fullgene_labels_HVG2000_benchmark_input.h5ad`
- `embedding_key`: `X_scGPT`
- model feature count: 512

The 512D `X_scGPT` stored in the HVG2000 h5ad was generated upstream from
full-gene observed expression:

```text
full-gene observed expression
-> scGPT whole-human encoder
-> 512D X_scGPT
-> logistic regression milestone classifier
```

Decision:

```text
The current embedding provider is a full-gene scGPT route. It must not be used
as the formal embedding-based provider for HVG2000-only projected expression.
It remains valid only for observed-cell reference/provenance and sensitivity
analyses that explicitly preserve the same full-gene route.
```

## Superseded Recommendation

The previous recommendation is now the selected formal policy. If projected
cells are available only as HVG2000 expression, and observed `X_scGPT` was
generated from full-gene expression, the clean formal option is to train a new
HVG2000-consistent embedding milestone provider:

```text
observed HVG2000 expression
-> HVG2000-compatible embedding route
-> new embedding milestone classifier
-> projected HVG2000 expression uses the same route
```

This would preserve comparable input semantics across observed and projected cells, at the cost of registering a new provider version.
