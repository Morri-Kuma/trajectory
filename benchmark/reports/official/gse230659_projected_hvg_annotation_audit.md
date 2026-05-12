# GSE230659 Projected HVG Annotation Audit

Generated: 2026-05-10T16:00:50.929894+00:00

## Scope

This audit covers scNODE A/B/C projected-cell annotation after switching the formal Embedding Coherence route to the HVG2000 PCA embedding provider plus CellTypist classifier provider.

Formal route:

```text
projected_expression.npy
-> hvg2000_pca_logistic_regression embedding provider
-> CellTypist classifier provider
-> consensus_milestone_label
```

## Scenario Summary

| Scenario | Cells | HVG Embedding | Provider Match | Ambiguous | Consensus ARI | Embedding ARI | Classifier ARI |
|---|---:|---|---:|---:|---:|---:|---:|
| A | 30000 | 30000x50 | 0.814 | 5593 (0.186) | 0.033077 | 0.025481 | 0.040033 |
| B | 18000 | 18000x50 | 0.658 | 6149 (0.342) | 0.000935 | 1.5e-05 | 0.008364 |
| C | 10000 | 10000x50 | 0.962 | 380 (0.038) | 0.003935 | 0.002605 | 0.008095 |

## Label Distributions

| Scenario | Mode | Label | Count | Fraction |
|---|---|---|---:|---:|
| A | consensus | epithelial_like | 19227 | 0.641 |
| A | consensus | hCiPS | 2562 | 0.085 |
| A | consensus | intermediate_plastic | 2618 | 0.087 |
| A | consensus | ambiguous | 5593 | 0.186 |
| A | embedding_based | epithelial_like | 23429 | 0.781 |
| A | embedding_based | hCiPS | 3072 | 0.102 |
| A | embedding_based | intermediate_plastic | 3499 | 0.117 |
| A | classifier_based | epithelial_like | 19473 | 0.649 |
| A | classifier_based | hCiPS | 5021 | 0.167 |
| A | classifier_based | intermediate_plastic | 5506 | 0.184 |
| B | consensus | epithelial_like | 11661 | 0.648 |
| B | consensus | hCiPS | 27 | 0.002 |
| B | consensus | intermediate_plastic | 163 | 0.009 |
| B | consensus | ambiguous | 6149 | 0.342 |
| B | embedding_based | epithelial_like | 17746 | 0.986 |
| B | embedding_based | hCiPS | 88 | 0.005 |
| B | embedding_based | intermediate_plastic | 166 | 0.009 |
| B | classifier_based | epithelial_like | 11687 | 0.649 |
| B | classifier_based | hCiPS | 5057 | 0.281 |
| B | classifier_based | intermediate_plastic | 1256 | 0.070 |
| C | consensus | epithelial_like | 9428 | 0.943 |
| C | consensus | hCiPS | 3 | 0.000 |
| C | consensus | intermediate_plastic | 189 | 0.019 |
| C | consensus | ambiguous | 380 | 0.038 |
| C | embedding_based | epithelial_like | 9664 | 0.966 |
| C | embedding_based | hCiPS | 90 | 0.009 |
| C | embedding_based | intermediate_plastic | 246 | 0.025 |
| C | classifier_based | epithelial_like | 9479 | 0.948 |
| C | classifier_based | hCiPS | 70 | 0.007 |
| C | classifier_based | intermediate_plastic | 451 | 0.045 |

## Provider Disagreement

Rows below show nonzero embedding-provider vs classifier-provider label combinations across all cells. Off-diagonal rows become `ambiguous` under the current consensus rule.

| Scenario | Embedding Label | Classifier Label | Count | Fraction | Agreement |
|---|---|---|---:|---:|---|
| A | epithelial_like | epithelial_like | 19227 | 0.641 | yes |
| A | epithelial_like | hCiPS | 1582 | 0.053 | no |
| A | epithelial_like | intermediate_plastic | 2620 | 0.087 | no |
| A | hCiPS | epithelial_like | 242 | 0.008 | no |
| A | hCiPS | hCiPS | 2562 | 0.085 | yes |
| A | hCiPS | intermediate_plastic | 268 | 0.009 | no |
| A | intermediate_plastic | epithelial_like | 4 | 0.000 | no |
| A | intermediate_plastic | hCiPS | 877 | 0.029 | no |
| A | intermediate_plastic | intermediate_plastic | 2618 | 0.087 | yes |
| B | epithelial_like | epithelial_like | 11661 | 0.648 | yes |
| B | epithelial_like | hCiPS | 5028 | 0.279 | no |
| B | epithelial_like | intermediate_plastic | 1057 | 0.059 | no |
| B | hCiPS | epithelial_like | 25 | 0.001 | no |
| B | hCiPS | hCiPS | 27 | 0.002 | yes |
| B | hCiPS | intermediate_plastic | 36 | 0.002 | no |
| B | intermediate_plastic | epithelial_like | 1 | 0.000 | no |
| B | intermediate_plastic | hCiPS | 2 | 0.000 | no |
| B | intermediate_plastic | intermediate_plastic | 163 | 0.009 | yes |
| C | epithelial_like | epithelial_like | 9428 | 0.943 | yes |
| C | epithelial_like | hCiPS | 17 | 0.002 | no |
| C | epithelial_like | intermediate_plastic | 219 | 0.022 | no |
| C | hCiPS | epithelial_like | 44 | 0.004 | no |
| C | hCiPS | hCiPS | 3 | 0.000 | yes |
| C | hCiPS | intermediate_plastic | 43 | 0.004 | no |
| C | intermediate_plastic | epithelial_like | 7 | 0.001 | no |
| C | intermediate_plastic | hCiPS | 50 | 0.005 | no |
| C | intermediate_plastic | intermediate_plastic | 189 | 0.019 | yes |

## Interpretation

- The technical route is coherent: all three scenarios report `formal_two_provider_complete=True`, `gene_universe=HVG2000`, exact feature-order matches, and 50-dimensional projected HVG PCA sidecar embeddings.
- Consensus ARI is low after the route switch, especially for B and C. This appears driven by provider disagreement and the conservative current rule that off-diagonal provider pairs become `ambiguous`.
- Provider agreement should be reviewed before changing consensus policy. A confidence-based tie-break may recover more evaluated cells, but it would change the formal interpretation and should be treated as a versioned policy decision.

