"""Train scVI then scANVI on the reference (production / server path).

Follows the scArches scANVI surgery tutorial:
https://docs.scarches.org/en/latest/scanvi_surgery_pipeline.html

scvi-tools / torch are imported lazily; if absent this raises a clear error so
callers fall back to the KNN surrogate in test_mode.
"""
from __future__ import annotations

from typing import Optional


def _require_scvi():
    try:
        import scvi  # noqa: F401
    except Exception as exc:  # pragma: no cover
        raise ImportError(
            "scvi-tools is required for the real scANVI path. "
            "Install scvi-tools + a torch build, or run in test_mode "
            "(uses the KNN surrogate)."
        ) from exc


def train_reference_scanvi(
    ref_adata,
    labels_key: str = "ref_label",
    batch_key: Optional[str] = "sample_id",
    counts_layer: str = "counts",
    unlabeled_category: str = "unknown",
    n_latent: int = 30,
    n_layers: int = 2,
    scvi_max_epochs: int = 300,
    scanvi_max_epochs: int = 50,
    model_dir: Optional[str] = None,
):
    """Return a trained SCANVI model (and save it if model_dir given)."""
    _require_scvi()
    import scvi

    scvi.model.SCVI.setup_anndata(
        ref_adata, layer=counts_layer, batch_key=batch_key
    )
    scvi_model = scvi.model.SCVI(
        ref_adata, n_latent=n_latent, n_layers=n_layers
    )
    scvi_model.train(max_epochs=scvi_max_epochs)

    scanvi_model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category,
    )
    scanvi_model.train(max_epochs=scanvi_max_epochs)

    if model_dir:
        scanvi_model.save(model_dir, overwrite=True)
    return scanvi_model


def reference_recovery_accuracy(scanvi_model, ref_adata, labels_key="ref_label"):
    """Sanity gate: how well scANVI re-predicts the reference's own labels."""
    import numpy as np

    pred = scanvi_model.predict(ref_adata)
    truth = ref_adata.obs[labels_key].astype(str).values
    return float(np.mean(np.asarray(pred) == truth))
