#!/usr/bin/env python
"""Positive control: scANVI self-transfer on GSE175634 (query HAS author labels).

Trains scANVI on an 80% split of GSE175634 (iPSC->cardiomyocyte; author
`final_milestone_label_coarse`, UNK excluded), maps the held-out 20% via scArches,
and reports test accuracy + per-state recall. If the method works when reference
and query are the same system, this isolates the GSE242424->chemical failure to
the cross-protocol gap. Uses the dataset's own `counts` layer; `sample_id` batch.
"""
from __future__ import annotations
import os, sys, json
import numpy as np, pandas as pd, anndata as ad

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
INPUT = "benchmark/inputs/gse175634_cardiac_author_hvg2000/GSE175634_cardiac_author_HVG2000_benchmark_input.h5ad"
LABEL, BATCH = "final_milestone_label_coarse", "sample_id"
OUT = "results/positive_control_gse175634"; TEST_FRAC = 0.2


def main():
    import scvi
    os.makedirs(OUT, exist_ok=True)
    A = ad.read_h5ad(INPUT)
    A = A[A.obs[LABEL].astype(str) != "UNK"].copy()
    if "counts" not in A.layers:
        raise SystemExit("GSE175634 input has no counts layer")
    A.obs["ref_label"] = A.obs[LABEL].astype(str).values
    rng = np.random.default_rng(0)
    test = np.zeros(A.n_obs, bool)
    lab = A.obs["ref_label"].values
    for s in np.unique(lab):
        idx = np.where(lab == s)[0]
        sel = rng.choice(idx, size=max(1, int(TEST_FRAC * len(idx))), replace=False)
        test[sel] = True
    train, query = A[~test].copy(), A[test].copy()

    scvi.model.SCVI.setup_anndata(train, layer="counts", batch_key=BATCH)
    m = scvi.model.SCVI(train, n_latent=30, n_layers=2); m.train(max_epochs=200)
    sa = scvi.model.SCANVI.from_scvi_model(m, labels_key="ref_label", unlabeled_category="unknown")
    sa.train(max_epochs=40)
    recovery = float(np.mean(np.asarray(sa.predict(train)) == train.obs["ref_label"].values))
    model_dir = os.path.join(OUT, "model"); sa.save(model_dir, overwrite=True)

    scvi.model.SCANVI.prepare_query_anndata(query, model_dir)
    qm = scvi.model.SCANVI.load_query_data(query, model_dir)
    qm.train(max_epochs=100, plan_kwargs={"weight_decay": 0.0})
    pred = np.asarray(qm.predict(query)); true = query.obs["ref_label"].values
    acc = float(np.mean(pred == true))
    per_state = {s: float(np.mean(pred[true == s] == s)) for s in np.unique(true)}
    pd.crosstab(pd.Series(true, name="true"), pd.Series(pred, name="pred")).to_csv(os.path.join(OUT, "confusion.csv"))
    json.dump({"reference_recovery": recovery, "test_accuracy": acc, "per_state_recall": per_state,
               "n_train": int(train.n_obs), "n_test": int(query.n_obs)},
              open(os.path.join(OUT, "positive_control_summary.json"), "w"), indent=2)
    print(f"POSITIVE CONTROL  recovery={recovery:.3f}  test_accuracy={acc:.3f}")
    print("per-state recall:", per_state, flush=True)


if __name__ == "__main__":
    main()
