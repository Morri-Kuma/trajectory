#!/usr/bin/env python
"""Train the reference scANVI model on the prepped GSE242424 reference and save it.

Server / GPU step. Reads a profile from config.yaml (default: `server`), loads the
prepped, labelled reference .h5ad, hard-checks that scVI is being fed raw integer
counts, trains scVI -> scANVI, verifies the reference label-recovery gate, and
saves the model to the configured models_dir.

Usage:
    python scripts/train_scanvi_reference.py --config config.yaml --profile server
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)


def _counts_like(X) -> bool:
    import scipy.sparse as sp
    data = X.data if sp.issparse(X) else np.asarray(X).ravel()
    if data.size == 0:
        return False
    sample = data[:200000]
    return bool(np.all(sample >= 0) and np.allclose(sample, np.round(sample)))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--profile", default="server")
    args = ap.parse_args()

    cfg = yaml.safe_load(open(args.config, encoding="utf-8"))
    prof = cfg[args.profile]
    ann = cfg["annotation"]
    scvi_cfg = prof.get("scvi", {})
    label_key = ann["reference_label_key"]
    ref_path = prof["paths"]["reference_h5ad"]
    model_dir = prof["paths"]["models_dir"]

    import anndata as ad
    ref = ad.read_h5ad(ref_path)
    if label_key not in ref.obs:
        raise KeyError(f"reference_label_key '{label_key}' not in {ref_path} obs")
    ref.obs["ref_label"] = ref.obs[label_key].astype(str).values

    # Hard gate: scVI requires raw integer counts.
    counts = ref.layers["counts"] if "counts" in ref.layers else ref.X
    if not _counts_like(counts):
        raise SystemExit(
            "ERROR: scVI needs raw integer counts, but the reference matrix is not "
            "count-like. Re-build the prepped input from a raw-count source "
            "(build_scanvi_inputs.py --counts-layer <counts>) before training."
        )
    if "counts" not in ref.layers:
        ref.layers["counts"] = ref.X.copy()

    from src.annotation.scanvi_train import (
        train_reference_scanvi, reference_recovery_accuracy,
    )
    print(f"[train] reference={ref_path} cells={ref.n_obs} genes={ref.n_vars} "
          f"label_key={label_key} -> model_dir={model_dir}", flush=True)

    model = train_reference_scanvi(
        ref,
        labels_key="ref_label",
        batch_key=scvi_cfg.get("batch_key", "sample_id"),
        n_latent=scvi_cfg.get("n_latent", 30),
        n_layers=scvi_cfg.get("n_layers", 2),
        scvi_max_epochs=scvi_cfg.get("scvi_max_epochs", 300),
        scanvi_max_epochs=scvi_cfg.get("scanvi_max_epochs", 50),
        model_dir=model_dir,
    )
    acc = reference_recovery_accuracy(model, ref)
    gate = scvi_cfg.get("reference_recovery_min_accuracy", 0.90)
    # also persist the frozen gene list for scArches query reindexing
    os.makedirs(model_dir, exist_ok=True)
    with open(os.path.join(model_dir, "reference_genes.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(map(str, ref.var_names)))
    print(json.dumps({"recovery_accuracy": acc, "gate": gate,
                      "model_dir": model_dir, "n_genes": int(ref.n_vars)}))
    if acc < gate:
        raise SystemExit(f"Recovery accuracy {acc:.3f} < gate {gate}: "
                         f"do NOT proceed to query mapping; inspect the reference.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
