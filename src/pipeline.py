"""End-to-end annotation-robustness pipeline.

Runs, per query dataset:
  load -> joint preprocess with reference (shared embedding) -> annotate
  (scANVI or surrogate) -> 6 comparisons (composition, agreement, markers,
  pseudotime, graph, disagreement) -> tables + figures.

Run:
    python -m src.pipeline --config config.yaml
"""
from __future__ import annotations

import argparse
import os
from typing import Dict

import numpy as np

try:
    import anndata as ad
except Exception:  # pragma: no cover
    ad = None

from .utils import load_config, resolve_path, set_global_seed, ensure_dir, save_json, save_csv
from .data import load_anndata, subsample_anndata, intersect_genes
from .preprocessing import basic_preprocess
from .annotation import annotate_query, ood_gate
from .evaluation import (
    composition_compare, agreement_matrix, marker_validation,
    pseudotime_compare, disagreement_cells,
)
from .trajectory import compute_pseudotime, compute_state_graph, graph_similarity
from .plotting import (
    plot_composition, plot_agreement_heatmap, plot_embedding, plot_pseudotime_scatter,
)


def _time_key(adata):
    for k in ("time_label", "stage_day_label", "sample_id"):
        if k in adata.obs:
            return k
    return adata.obs.columns[0]


def _load_inputs(cfg) -> Dict:
    active = cfg["active"]
    if cfg["mode"] == "test":
        sd = resolve_path(cfg, active["paths"]["synthetic_dir"])
        ref = load_anndata(os.path.join(sd, "reference.h5ad"))
        queries = {q: load_anndata(os.path.join(sd, f"{q}.h5ad"))
                   for q in cfg["annotation"]["query_datasets"]}
    else:
        ref = load_anndata(resolve_path(cfg, active["paths"]["reference_h5ad"]))
        queries = {q: load_anndata(resolve_path(cfg, p))
                   for q, p in active["paths"]["query_h5ad"].items()}
    return {"reference": ref, "queries": queries}


def _harmonize(reference, query):
    genes, harmonized = intersect_genes([reference, query])
    return genes, (harmonized[0], harmonized[1])


def run(config_path: str = "config.yaml") -> Dict:
    cfg = load_config(config_path)
    seed = cfg["project"]["random_seed"]
    set_global_seed(seed)
    active = cfg["active"]
    out_dir = ensure_dir(resolve_path(cfg, active["paths"]["output_dir"]))

    max_cells = active.get("max_cells_per_dataset")
    max_genes = active.get("max_genes")
    n_pcs = active.get("trajectory", {}).get("n_pcs", 30)
    n_neighbors = active.get("trajectory", {}).get("n_neighbors", 15)
    ms_key = cfg["annotation"]["marker_silver_label_key"]
    ref_label_key = cfg["annotation"]["reference_label_key"]
    root_state = cfg["trajectory"]["root_state"]
    corr = cfg["annotation"]["correspondence_map"]
    markers = cfg["markers"]

    inputs = _load_inputs(cfg)
    reference = subsample_anndata(inputs["reference"], max_cells, seed)
    if ref_label_key in reference.obs:
        reference.obs["ref_label"] = reference.obs[ref_label_key].astype(str).values

    summary = {"mode": cfg["mode"], "config_hash": cfg["_config_hash"], "queries": {}}

    for qname, q in inputs["queries"].items():
        q = subsample_anndata(q, max_cells, seed)
        qdir = ensure_dir(os.path.join(out_dir, qname))

        # joint preprocessing for a shared embedding (surrogate needs it)
        genes, (ref_h, q_h) = _harmonize(reference, q)
        joint = ad.concat({"ref": ref_h, "query": q_h}, label="origin", index_unique="-")
        joint = basic_preprocess(joint, n_pcs=n_pcs, n_neighbors=n_neighbors,
                                 max_genes=max_genes, seed=seed)
        ref_sub = joint[joint.obs["origin"] == "ref"].copy()
        q_sub = joint[joint.obs["origin"] == "query"].copy()
        # drop joint-indexed graph artefacts so trajectory code rebuilds clean
        # per-subset kNN graphs from the shared-space X_pca.
        for a in (ref_sub, q_sub):
            a.uns.pop("knn_indices", None)
            a.uns.pop("knn_n_neighbors", None)
            a.uns.pop("neighbors", None)
            try:
                del a.obsp["distances"]
            except Exception:
                pass
        q_sub.obs[ms_key] = q_h.obs[ms_key].astype(str).values
        if "ref_label" in ref_h.obs:
            ref_sub.obs["ref_label"] = ref_h.obs["ref_label"].astype(str).values
        elif ref_label_key in ref_h.obs:
            ref_sub.obs["ref_label"] = ref_h.obs[ref_label_key].astype(str).values

        # annotation (scANVI or surrogate)
        annotate_query(ref_sub, q_sub, cfg)
        q_sub.obs["scanvi_mapped"] = [corr.get(x, "unknown_or_ood")
                                      for x in q_sub.obs["scanvi_label"].astype(str)]

        # Geometric OOD reference gate (integrated into the pipeline; the same
        # routine runs on the server scANVI latent via scripts/score_ood.py).
        # Prefer an scANVI latent if present, else the shared PCA latent.
        ood_cfg = cfg["annotation"].get("ood", {})
        _ood_lat = "X_scANVI" if ("X_scANVI" in ref_sub.obsm and "X_scANVI" in q_sub.obsm) else "X_pca"
        ood_gate_res = ood_gate(
            ref_sub.obsm[_ood_lat], q_sub.obsm[_ood_lat],
            k=int(ood_cfg.get("k", 15)),
            quantile=float(ood_cfg.get("quantile", 0.95)),
            state_labels=q_sub.obs[ms_key].values,
        )
        ood_gate_res["latent"] = _ood_lat
        save_json(ood_gate_res, os.path.join(qdir, "ood_gate.json"))

        tkey = _time_key(q_sub)

        # 1. composition
        comp = composition_compare(q_sub, ms_key, "scanvi_mapped", tkey)
        save_csv(comp["table"], os.path.join(qdir, "composition_compare.csv"), index=False)

        # 2/3. agreement
        agr = agreement_matrix(q_sub, ms_key, "scanvi_label", corr)
        save_csv(agr["contingency"], os.path.join(qdir, "agreement_contingency.csv"))
        save_json(agr["metrics"], os.path.join(qdir, "agreement_metrics.json"))

        # 4. marker validation (both annotations)
        mv_a = marker_validation(q_sub, markers, ms_key)
        mv_b = marker_validation(q_sub, markers, "scanvi_label")
        save_csv(mv_a, os.path.join(qdir, "marker_validation_marker_silver.csv"), index=False)
        save_csv(mv_b, os.path.join(qdir, "marker_validation_scanvi.csv"), index=False)

        # 5. pseudotime under each annotation's root
        pt_a = compute_pseudotime(q_sub, root_state, ms_key, n_dcs=cfg["trajectory"]["n_dcs"])
        pt_b = compute_pseudotime(q_sub, root_state, "scanvi_mapped", n_dcs=cfg["trajectory"]["n_dcs"])
        q_sub.obs["pseudotime_marker_silver"] = pt_a
        q_sub.obs["pseudotime_scanvi"] = pt_b
        pt_cmp = pseudotime_compare(pt_a, pt_b)
        save_json(pt_cmp, os.path.join(qdir, "pseudotime_corr.json"))

        # 5b. state graphs + similarity (same label space via mapping)
        g_a = compute_state_graph(q_sub, ms_key)
        g_b = compute_state_graph(q_sub, "scanvi_mapped")
        gsim = graph_similarity(g_a["edges"], g_b["edges"])
        save_json({"graph_similarity": gsim,
                   "edges_marker_silver": g_a["edges"],
                   "edges_scanvi_mapped": g_b["edges"],
                   "backend": g_a["backend"]},
                  os.path.join(qdir, "paga_graph_compare.json"))

        # 6. disagreement
        dis = disagreement_cells(q_sub, ms_key, "scanvi_label", corr)
        save_csv(dis, os.path.join(qdir, "disagreement_cells.csv"))

        # figures
        figs = {}
        try:
            figs["composition"] = plot_composition(comp["table"], os.path.join(qdir, "fig_composition.png"))
            figs["agreement"] = plot_agreement_heatmap(agr["contingency"], os.path.join(qdir, "fig_agreement.png"))
            figs["embed_silver"] = plot_embedding(q_sub, ms_key, os.path.join(qdir, "fig_embed_marker_silver.png"), title=f"{qname}: marker-silver")
            figs["embed_scanvi"] = plot_embedding(q_sub, "scanvi_label", os.path.join(qdir, "fig_embed_scanvi.png"), title=f"{qname}: scANVI")
            figs["pseudotime"] = plot_pseudotime_scatter(pt_a, pt_b, os.path.join(qdir, "fig_pseudotime.png"), title=f"{qname}: pseudotime A vs B")
        except Exception as exc:
            figs["error"] = str(exc)

        summary["queries"][qname] = {
            "n_cells": int(q_sub.n_obs),
            "n_genes_shared": int(len(genes)),
            "annotation_method": str(q_sub.obs["annotation_method"].iloc[0]),
            "overall_tv_composition": comp["overall_tv"],
            "agreement_nmi": agr["metrics"]["nmi"],
            "agreement_ari": agr["metrics"]["ari"],
            "correspondence_overall_agreement": agr["metrics"].get("correspondence_overall_agreement"),
            "pseudotime_spearman": pt_cmp["spearman"],
            "graph_edge_jaccard": gsim["edge_jaccard"],
            "frac_unknown_ood": float(np.mean(q_sub.obs["scanvi_unknown_flag"].values)),
            "frac_ood_geometric": ood_gate_res["frac_ood_geometric"],
            "ood_gate_decision": ood_gate_res["gate_decision"],
            "ood_enrichment_vs_ref": ood_gate_res["enrichment_vs_ref"],
            "n_disagreement_cells": int(dis.shape[0]),
            "figures": figs,
        }

    save_json(summary, os.path.join(out_dir, "summary.json"))
    return summary


def main():
    ap = argparse.ArgumentParser(description="Annotation-robustness pipeline")
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()
    s = run(args.config)
    print("Pipeline complete. Summary:")
    import json
    print(json.dumps({k: v for k, v in s.items() if k != "queries"}, indent=2))
    for q, info in s["queries"].items():
        print(f"  [{q}] cells={info['n_cells']} method={info['annotation_method']} "
              f"NMI={info['agreement_nmi']:.3f} pt_spearman={info['pseudotime_spearman']:.3f} "
              f"graphJaccard={info['graph_edge_jaccard']}")


if __name__ == "__main__":
    main()
