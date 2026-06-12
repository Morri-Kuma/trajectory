"""Annotation dispatcher: choose scANVI/scArches (server) or KNN surrogate
(test_mode / no scvi), and write results into query .obs.
"""
from __future__ import annotations

from typing import Dict

import numpy as np

from .knn_surrogate import knn_surrogate_annotate


def annotate_query(ref_adata, query_adata, cfg) -> Dict:
    """Annotate ``query_adata`` against ``ref_adata`` per config.

    For the surrogate route both AnnData must already carry a SHARED-space
    ``obsm['X_pca']`` (the pipeline achieves this via joint preprocessing).
    Writes obs columns: scanvi_label, scanvi_confidence, scanvi_unknown_flag,
    annotation_method. Returns the raw result dict.
    """
    active = cfg["active"]
    scvi_cfg = active.get("scvi", {})
    tau = cfg["annotation"]["confidence_threshold"]
    ref_label_key = cfg["annotation"]["reference_label_key"]

    use_scvi = bool(scvi_cfg.get("use_scvi", False))
    result = None

    if use_scvi:
        try:
            from .scanvi_map import map_query_scanvi

            model_dir = active["paths"]["models_dir"]
            result = map_query_scanvi(
                query_adata,
                ref_model_path=model_dir,
                batch_key=scvi_cfg.get("batch_key", "sample_id"),
                query_max_epochs=scvi_cfg.get("scanvi_max_epochs", 100),
                confidence_threshold=tau,
            )
        except ImportError:
            print("[annotate_query] scvi-tools unavailable -> KNN surrogate.")
            use_scvi = False

    if not use_scvi:
        if "X_pca" not in ref_adata.obsm or "X_pca" not in query_adata.obsm:
            raise ValueError(
                "Surrogate annotation needs shared-space obsm['X_pca'] on both "
                "AnnData (run joint preprocessing first)."
            )
        ref_labels = ref_adata.obs[ref_label_key].astype(str).values \
            if ref_label_key in ref_adata.obs else ref_adata.obs["ref_label"].astype(str).values
        result = knn_surrogate_annotate(
            ref_adata.obsm["X_pca"],
            ref_labels,
            query_adata.obsm["X_pca"],
            k=active.get("trajectory", {}).get("n_neighbors", 15),
            confidence_threshold=tau,
        )

    query_adata.obs["scanvi_label"] = np.asarray(result["label"]).astype(str)
    query_adata.obs["scanvi_confidence"] = np.asarray(result["confidence"])
    query_adata.obs["scanvi_unknown_flag"] = np.asarray(result["unknown_flag"])
    query_adata.obs["annotation_method"] = np.asarray(result["method"]).astype(str)
    query_adata.uns["scanvi_branch"] = {
        "reference_label_key": ref_label_key,
        "confidence_threshold": tau,
        "config_hash": cfg.get("_config_hash"),
        "method": str(result["method"][0]) if len(result["method"]) else "none",
    }
    return result
