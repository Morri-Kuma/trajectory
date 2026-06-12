#!/usr/bin/env python
"""Re-score scANVI query labels with a geometric OOD flag, on the scANVI latent.

scANVI posterior confidence misses cross-protocol out-of-distribution cells (the
`hOSK`-collapse). This computes, per query cell, the distance to the nearest
REFERENCE cells in the scANVI latent and flags cells beyond the reference's 95th
percentile. Server step (needs scvi-tools); reuses src/annotation/ood.py.

Usage: python scripts/score_ood.py --config config.yaml --profile server
"""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np, pandas as pd, anndata as ad, yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from src.annotation import ood_score_and_flag, ood_gate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--profile", default="server")
    ap.add_argument("--k", type=int, default=15)
    ap.add_argument("--quantile", type=float, default=0.95)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config)); prof = cfg[a.profile]; ann = cfg["annotation"]
    model_dir = prof["paths"]["models_dir"]; out_dir = prof["paths"]["output_dir"]

    import scvi
    ref = ad.read_h5ad(prof["paths"]["reference_h5ad"])
    ref.obs["ref_label"] = ref.obs[ann["reference_label_key"]].astype(str).values
    ref_model = scvi.model.SCANVI.load(model_dir, adata=ref)
    ref_model.is_trained_ = True  # loaded weights are trained; restore flag for get_latent
    ref_lat = ref_model.get_latent_representation(ref)

    summary = {}
    for q, path in prof["paths"]["query_h5ad"].items():
        A = ad.read_h5ad(path)
        scvi.model.SCANVI.prepare_query_anndata(A, model_dir)
        qm = scvi.model.SCANVI.load_query_data(A, model_dir)
        qm.is_trained_ = True  # project query through the reference-trained encoder (OOD probe; no surgery retrain)
        q_lat = qm.get_latent_representation(A)
        res = ood_score_and_flag(ref_lat, q_lat, k=a.k, quantile=a.quantile)
        # transfer-level decision via the SAME shared gate the pipeline uses
        gate = ood_gate(ref_lat, q_lat, k=a.k, quantile=a.quantile)
        os.makedirs(os.path.join(out_dir, q), exist_ok=True)
        pd.DataFrame({"ood_score": res["ood_score"], "ood_flag": res["ood_flag"]},
                     index=A.obs_names).to_csv(os.path.join(out_dir, q, "ood_scanvi_latent.csv"))
        summary[q] = {"frac_ood": res["frac_ood"],
                      "query_median_knn_dist": res["query_median_knn_dist"],
                      "ref_median_knn_dist": res["ref_median_knn_dist"],
                      "enrichment_vs_ref": gate["enrichment_vs_ref"],
                      "gate_decision": gate["gate_decision"]}
        print(f"[{q}] geometric OOD on scANVI latent: {res['frac_ood']*100:.1f}%"
              f"  -> {gate['gate_decision']} ({gate['enrichment_vs_ref']:.1f}x ref)", flush=True)
    json.dump(summary, open(os.path.join(out_dir, "ood_summary_scanvi_latent.json"), "w"), indent=2)


if __name__ == "__main__":
    main()
