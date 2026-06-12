# scANVI server run — first results & honest interpretation

**Date:** 2026-06-09 (run synced from Shirokane)
**Run:** `mode: server`, reference = GSE242424 `author_cluster_label` (2000 shared HVG, reconstructed counts), queries = GSE178325 (80,475 cells) + GSE230659 (75,194 cells), full data, real scVI→scANVI→scArches.

## What worked
- **Reference recovery = 0.990** (gate 0.90): scANVI learns the OSKM reference's own labels near-perfectly. The training machinery and the reconstructed-count inputs are sound.
- The full pipeline executed end-to-end on the cluster: model saved to `models/scanvi_reference_gse242424/`, both queries mapped, comparison tables + figures written to `results/annotation_branch/`.

## The key finding (a degenerate cross-protocol mapping)
scANVI/scArches **collapsed almost all query cells onto a single out-of-distribution reference label, `hOSK`** (the OSKM-transgene state), with high confidence:

| query | dominant scANVI label | share | OOD-flagged (conf<0.5) |
|---|---|---:|---:|
| GSE178325 | `hOSK` | 77,486 / 80,475 (**96.3%**) | 0.04% |
| GSE230659 | `hOSK` | 57,606 / 75,194 (**76.6%**) | 0.30% |

`hOSK` was only **498 cells (0.8%)** of the reference, and it is exactly one of the two states (`hOSK`/`xOSK`) the pre-registered correspondence map flagged as having **no chemical analogue**. So:

- The low agreement vs marker-silver (NMI 0.027 / 0.101; ARI 0.010 / 0.077; correspondence-overall 0.04 / 0.01; ~99.8% of cells "disagree") is **driven by this collapse**, not by a biologically meaningful state-by-state difference.
- The near-zero OOD fraction shows **scANVI's posterior confidence does NOT detect the cross-protocol mismatch** — it assigns the chemical cells to the nearest reference latent region (`hOSK`) and reports high confidence anyway. This is a known failure mode of classifier-style reference mapping on out-of-distribution populations.
- The near-perfect pseudotime correlation (Spearman 0.9997 / 0.996) is **confounded**: with almost no query cells mapping to `hADSCs`, the annotation-B pseudotime root degenerated, so the two pseudotimes are not a fair "robustness" comparison here.

## Honest conclusion
This run does **not** confirm the original hypotheses; it produces a **clear, interpretable negative result**: *naive cross-protocol scANVI reference mapping (OSKM → chemical reprogramming) is not a valid annotation strategy for these queries — it yields confident-but-spurious labels and its confidence is miscalibrated for OOD.* This is itself publishable and is exactly the risk the design's correspondence-map caveat and risk-management section anticipated. It is **not** a code bug (reference recovery 0.99; valid label categories; structured, non-uniform output).

## Recommended next steps (in priority order)
1. **Positive control — GSE175634** (Fig S3): run the same scANVI machinery where the query *has* author labels (train on part, map the rest). If recovery+mapping there are sensible, it proves the method works when systems match, isolating the failure to the OSKM↔chemical gap. **This is the decisive diagnostic.**
2. **Recalibrate OOD detection:** replace/augment scANVI posterior confidence with a **latent-distance / k-NN-to-reference OOD score** (or scArches reconstruction error). Expectation: the `hOSK`-collapsed cells should score as OOD. This turns the negative result into a method contribution.
3. **Same-protocol transfer (sanity + alternative):** use one chemical dataset's marker-silver labels as the reference for the other (GSE230659 ↔ GSE178325). A within-protocol scANVI transfer should agree far better and gives a meaningful annotation comparison.
4. **Reframe the manuscript:** the story becomes *"when is reference-mapping annotation trustworthy for reprogramming trajectories?"* — with the OSKM→chemical collapse as the cautionary core, the positive control + same-protocol transfer as the "what works," and the OOD recalibration as the fix.

## Artifacts
`results/annotation_branch/{gse178325,gse230659}/` — `agreement_contingency.csv`,
`agreement_metrics.json`, `composition_compare.csv`, `marker_validation_*.csv`,
`pseudotime_corr.json`, `paga_graph_compare.json`, `disagreement_cells.csv`,
figures; `summary.json`; reference model under `models/scanvi_reference_gse242424/`.


---

## Validation results (completed 2026-06-09)
- **Positive control (GSE175634 self-transfer):** test accuracy **0.944**, per-state recall 0.93–0.96 → the machinery is sound; the cross-protocol failure is real, not a bug.
- **Same-protocol (chemical↔chemical):** correspondence **0.863** (230→178) / **0.757** (178→230) vs **0.01–0.04** cross-protocol — a 20–60× improvement. Carried by abundant somatic/epithelial states; rare terminal/intermediate states remain limited by reference abundance.
- **Money figure:** `results/figures/validation_summary.png`. Manuscript updated (`manuscript/manuscript_draft.{md,docx}`, §3.4–3.5). The article now has a complete positive+negative+remedy story with all real numbers.

- **scANVI-latent OOD (re-run 2026-06-10):** geometric OOD on the actual scANVI latent flags **75.9%/73.0%** of query cells vs scANVI confidence ≤0.3% — even stronger than the PCA proxy (38.7%/43.6%). Fig `results/figures/ood_confidence_vs_geometric.png`; `results/annotation_branch/ood_summary_scanvi_latent.json`.
