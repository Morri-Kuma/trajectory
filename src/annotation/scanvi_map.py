"""scArches surgery: map a query onto the frozen reference scANVI model and
predict labels with confidence (production / server path).
"""
from __future__ import annotations

from typing import Dict

import numpy as np


def map_query_scanvi(
    query_adata,
    ref_model_path: str,
    batch_key: str = "sample_id",
    query_max_epochs: int = 100,
    confidence_threshold: float = 0.5,
) -> Dict[str, np.ndarray]:
    """Load reference, run surgery on the query, predict labels + confidence.

    Returns dict with 'label', 'confidence', 'unknown_flag', 'latent', 'method'.
    """
    try:
        import scvi
    except Exception as exc:  # pragma: no cover
        raise ImportError("scvi-tools required for map_query_scanvi") from exc

    scvi.model.SCANVI.prepare_query_anndata(query_adata, ref_model_path)
    q_model = scvi.model.SCANVI.load_query_data(query_adata, ref_model_path)
    q_model.train(max_epochs=query_max_epochs, plan_kwargs={"weight_decay": 0.0})

    pred = np.asarray(q_model.predict(query_adata))
    proba = q_model.predict(query_adata, soft=True)
    conf = np.asarray(proba.max(axis=1)).ravel()
    latent = q_model.get_latent_representation(query_adata)

    unknown = conf < confidence_threshold
    label = pred.astype(object).copy()
    label[unknown] = "unknown_or_ood"
    return {
        "label": label,
        "confidence": conf,
        "unknown_flag": unknown,
        "latent": latent,
        "method": np.array(["scanvi_scarches"] * len(pred)),
    }
