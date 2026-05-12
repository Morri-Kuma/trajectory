"""Train an HVG2000 expression -> PCA embedding milestone provider.

This provider is intended for projected-cell formal annotation when projected
expression is available only in the benchmark HVG gene universe.  It keeps the
embedding route inside the same observed/projected expression space:

    observed HVG expression -> PCA latent embedding -> milestone classifier

The saved joblib asset is a dictionary containing the fitted PCA transformer,
post-PCA scaler, classifier, feature order, and metadata.  Downstream projected
annotation should align projected expression to ``var_names`` before calling the
transformer/classifier.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Optional


EXCLUDE_DEFAULT = ("ambiguous", "unknown", "NA", "nan", "")


def _json_default(obj):
    try:
        import numpy as np

        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.ndarray,)):
            return obj.tolist()
    except Exception:
        pass
    return str(obj)


def _to_dense_float32(X):
    import numpy as np

    try:
        import scipy.sparse as sp

        if sp.issparse(X):
            return X.toarray().astype("float32", copy=False)
    except Exception:
        pass
    return np.asarray(X, dtype="float32")


def _balanced_indices(labels, max_per_class: Optional[int], random_state: int):
    import numpy as np

    rng = np.random.default_rng(random_state)
    keep: List[int] = []
    for label in sorted(set(labels.tolist())):
        idx = np.where(labels == label)[0]
        if max_per_class and len(idx) > max_per_class:
            idx = rng.choice(idx, size=max_per_class, replace=False)
        keep.extend(idx.tolist())
    return np.array(sorted(keep), dtype=int)


def _sample_indices(n: int, max_cells: Optional[int], random_state: int):
    import numpy as np

    if not max_cells or max_cells <= 0 or n <= max_cells:
        return np.arange(n, dtype=int)
    rng = np.random.default_rng(random_state)
    return np.sort(rng.choice(np.arange(n), size=max_cells, replace=False))


def _write_predictions_csv(path: Path, cell_ids, labels, conf, classes, probs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "cell_id",
        "milestone_embedding_label",
        "milestone_embedding_confidence",
    ] + [f"embedding_prob_{c}" for c in classes]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, cid in enumerate(cell_ids):
            row = {
                "cell_id": str(cid),
                "milestone_embedding_label": str(labels[i]),
                "milestone_embedding_confidence": float(conf[i]),
            }
            for j, cls in enumerate(classes):
                row[f"embedding_prob_{cls}"] = float(probs[i, j])
            writer.writerow(row)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Train HVG expression PCA milestone embedding provider."
    )
    p.add_argument("--input-h5ad", required=True)
    p.add_argument("--output-h5ad", default=None)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--dataset-id", required=True)
    p.add_argument("--label-key", default="milestone_marker_label")
    p.add_argument("--exclude-label", action="append", default=list(EXCLUDE_DEFAULT))
    p.add_argument("--n-components", type=int, default=50)
    p.add_argument(
        "--max-pca-fit-cells",
        type=int,
        default=0,
        help="Maximum cells for unsupervised PCA fit; 0 means all cells.",
    )
    p.add_argument("--max-train-cells-per-class", type=int, default=12000)
    p.add_argument("--test-size", type=float, default=0.2)
    p.add_argument("--random-state", type=int, default=42)
    p.add_argument("--C", type=float, default=1.0)
    p.add_argument("--max-iter", type=int, default=1000)
    p.add_argument("--embedding-key", default="X_hvg_pca_milestone")
    p.add_argument("--dry-run", action="store_true")
    return p


def main() -> int:
    args = build_parser().parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        meta = {
            "script": "build_hvg_pca_embedding_milestone_model",
            "status": "dry_run",
            "provider_family": "embedding_based",
            "method": "hvg2000_pca_logistic_regression",
            "dataset_id": args.dataset_id,
            "input_h5ad": args.input_h5ad,
            "label_key": args.label_key,
            "n_components": args.n_components,
        }
        with open(
            out_dir / "hvg_pca_embedding_milestone_model_metadata.json",
            "w",
            encoding="utf-8",
        ) as f:
            json.dump(meta, f, indent=2)
            f.write("\n")
        print("[build_hvg_pca_embedding_milestone_model] dry-run metadata written")
        return 0

    import anndata as ad
    import joblib
    import numpy as np
    import pandas as pd
    from sklearn.decomposition import PCA
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import classification_report, confusion_matrix
    from sklearn.model_selection import train_test_split
    from sklearn.preprocessing import StandardScaler

    print(f"[build_hvg_pca_embedding_milestone_model] loading {args.input_h5ad}")
    adata = ad.read_h5ad(args.input_h5ad)
    if args.label_key not in adata.obs:
        raise KeyError(f"obs label key {args.label_key!r} not found")
    if args.n_components <= 0 or args.n_components > min(adata.n_obs, adata.n_vars):
        raise ValueError(
            f"--n-components must be in [1, {min(adata.n_obs, adata.n_vars)}], "
            f"got {args.n_components}"
        )

    var_names = adata.var_names.astype(str).to_numpy()
    y_all = adata.obs[args.label_key].astype(str).to_numpy()
    excluded = set(args.exclude_label or [])
    train_mask = ~np.isin(y_all, list(excluded))
    if train_mask.sum() == 0:
        raise ValueError("No labelled training cells after excluding labels")

    print("[build_hvg_pca_embedding_milestone_model] materializing HVG expression")
    X_all = _to_dense_float32(adata.X)
    pca_fit_indices = _sample_indices(
        adata.n_obs, args.max_pca_fit_cells, args.random_state
    )

    pca = PCA(
        n_components=args.n_components,
        svd_solver="randomized",
        random_state=args.random_state,
    )
    print(
        "[build_hvg_pca_embedding_milestone_model] fitting PCA on "
        f"{len(pca_fit_indices)} cells x {adata.n_vars} genes"
    )
    pca.fit(X_all[pca_fit_indices])
    X_pca_all = pca.transform(X_all).astype("float32", copy=False)

    train_indices_all = np.where(train_mask)[0]
    sampled_rel = _balanced_indices(
        y_all[train_indices_all],
        max_per_class=args.max_train_cells_per_class,
        random_state=args.random_state,
    )
    train_indices = train_indices_all[sampled_rel]
    X = X_pca_all[train_indices]
    y = y_all[train_indices]
    classes = sorted(set(y.tolist()))
    if len(classes) < 2:
        raise ValueError(f"Need at least two classes to train classifier; got {classes}")

    stratify = y if min(pd.Series(y).value_counts()) >= 2 else None
    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=args.test_size,
        random_state=args.random_state,
        stratify=stratify,
    )

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    classifier = LogisticRegression(
        C=args.C,
        max_iter=args.max_iter,
        class_weight="balanced",
        random_state=args.random_state,
        n_jobs=1,
    )
    print(
        "[build_hvg_pca_embedding_milestone_model] training classifier on "
        f"{len(y_train)} cells; validating on {len(y_test)}"
    )
    classifier.fit(X_train_scaled, y_train)

    y_pred_test = classifier.predict(X_test_scaled)
    report = classification_report(
        y_test, y_pred_test, output_dict=True, zero_division=0
    )
    cm = confusion_matrix(y_test, y_pred_test, labels=classifier.classes_)

    X_pca_all_scaled = scaler.transform(X_pca_all)
    probs_all = classifier.predict_proba(X_pca_all_scaled)
    labels_all = classifier.classes_[np.argmax(probs_all, axis=1)]
    conf_all = np.max(probs_all, axis=1)

    adata.obsm[args.embedding_key] = X_pca_all
    adata.obs["milestone_embedding_label"] = labels_all.astype(str)
    adata.obs["milestone_embedding_confidence"] = conf_all.astype("float32")
    for j, cls in enumerate(classifier.classes_):
        adata.obs[f"embedding_prob_{cls}"] = probs_all[:, j].astype("float32")

    model_path = out_dir / "hvg_pca_embedding_milestone_model.joblib"
    model_asset: Dict = {
        "provider_family": "embedding_based",
        "method": "hvg2000_pca_logistic_regression",
        "dataset_id": args.dataset_id,
        "pca": pca,
        "scaler": scaler,
        "classifier": classifier,
        "var_names": var_names.tolist(),
        "embedding_key": args.embedding_key,
        "label_key": args.label_key,
        "excluded_labels": sorted(excluded),
        "n_features_in_": int(adata.n_vars),
        "n_components": int(args.n_components),
        "classes_": classifier.classes_.astype(str).tolist(),
    }
    joblib.dump(model_asset, model_path)

    pred_path = out_dir / "hvg_pca_embedding_milestone_predictions.csv"
    _write_predictions_csv(
        pred_path, adata.obs_names, labels_all, conf_all, classifier.classes_, probs_all
    )

    cm_path = out_dir / "hvg_pca_embedding_milestone_confusion_matrix.csv"
    pd.DataFrame(cm, index=classifier.classes_, columns=classifier.classes_).to_csv(
        cm_path
    )

    if args.output_h5ad:
        print(f"[build_hvg_pca_embedding_milestone_model] writing {args.output_h5ad}")
        adata.write_h5ad(args.output_h5ad)

    meta: Dict = {
        "script": "build_hvg_pca_embedding_milestone_model",
        "status": "completed",
        "provider_family": "embedding_based",
        "method": "hvg2000_pca_logistic_regression",
        "dataset_id": args.dataset_id,
        "input_h5ad": args.input_h5ad,
        "output_h5ad": args.output_h5ad,
        "gene_universe": "HVG2000",
        "expression_source": "adata.X",
        "embedding_key": args.embedding_key,
        "training_label_key": args.label_key,
        "excluded_labels": sorted(excluded),
        "classes": list(map(str, classifier.classes_)),
        "n_cells_total": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_pca_components": int(args.n_components),
        "n_pca_fit_cells": int(len(pca_fit_indices)),
        "pca_explained_variance_ratio_sum": float(
            np.sum(pca.explained_variance_ratio_)
        ),
        "n_training_pool": int(train_mask.sum()),
        "n_training_sampled": int(len(train_indices)),
        "n_train": int(len(y_train)),
        "n_validation": int(len(y_test)),
        "validation_report": report,
        "model_path": str(model_path),
        "predictions_csv": str(pred_path),
        "confusion_matrix_csv": str(cm_path),
        "feature_order": "input_h5ad.var_names",
        "formal_projected_route": "HVG2000-consistent",
        "notes": (
            "This provider is intended to replace the full-gene scGPT embedding "
            "provider for HVG2000-only projected-cell formal annotation."
        ),
    }
    with open(
        out_dir / "hvg_pca_embedding_milestone_model_metadata.json",
        "w",
        encoding="utf-8",
    ) as f:
        json.dump(meta, f, indent=2, default=_json_default)
        f.write("\n")

    print("[build_hvg_pca_embedding_milestone_model] completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
