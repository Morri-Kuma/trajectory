"""Train a classifier-based milestone annotation provider.

Primary backend implemented here is CellTypist, used as a lightweight
logistic-regression-style cell annotation tool.  A sklearn expression fallback
is also available for environments where CellTypist cannot train on the input.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Optional


EXCLUDE_DEFAULT = ("ambiguous", "unknown", "NA", "nan", "")
SUPPORTED_BACKENDS = ("celltypist", "sklearn_expression")


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


def _write_predictions_csv(path: Path, cell_ids, labels, conf, classes, probs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "cell_id",
        "milestone_classifier_label",
        "milestone_classifier_confidence",
    ] + [f"classifier_prob_{c}" for c in classes]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for i, cid in enumerate(cell_ids):
            row = {
                "cell_id": str(cid),
                "milestone_classifier_label": str(labels[i]),
                "milestone_classifier_confidence": float(conf[i]),
            }
            for j, cls in enumerate(classes):
                row[f"classifier_prob_{cls}"] = float(probs[i, j])
            writer.writerow(row)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Train classifier-based milestone annotation provider.")
    p.add_argument("--input-h5ad", required=True)
    p.add_argument("--query-h5ad", default=None)
    p.add_argument("--output-h5ad", default=None)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--dataset-id", required=True)
    p.add_argument("--label-key", default="milestone_marker_label")
    p.add_argument("--backend", choices=SUPPORTED_BACKENDS, default="celltypist")
    p.add_argument("--exclude-label", action="append", default=list(EXCLUDE_DEFAULT))
    p.add_argument("--max-train-cells-per-class", type=int, default=12000)
    p.add_argument("--random-state", type=int, default=42)
    p.add_argument("--celltypist-check-expression", action="store_true")
    p.add_argument("--max-iter", type=int, default=300)
    p.add_argument("--dry-run", action="store_true")
    return p


def _fit_predict_sklearn(adata_train, labels_train, adata_query, args):
    import joblib
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    X_train = _to_dense_float32(adata_train.X)
    X_query = _to_dense_float32(adata_query.X)
    model = make_pipeline(
        StandardScaler(with_mean=False),
        LogisticRegression(
            max_iter=args.max_iter,
            class_weight="balanced",
            random_state=args.random_state,
            n_jobs=1,
        ),
    )
    model.fit(X_train, labels_train)
    probs = model.predict_proba(X_query)
    labels = model.classes_[np.argmax(probs, axis=1)]
    return model, model.classes_, labels, probs


def _fit_predict_celltypist(adata_train, labels_train, adata_query, args):
    import numpy as np
    import celltypist

    adata_train = _prepare_celltypist_adata(adata_train)
    adata_query = _prepare_celltypist_adata(adata_query)
    model = celltypist.train(
        X=adata_train,
        labels=labels_train,
        check_expression=bool(args.celltypist_check_expression),
        max_iter=args.max_iter,
        n_jobs=1,
    )
    result = celltypist.annotate(adata_query, model=model, mode="best match")
    pred_df = result.predicted_labels
    prob_df = result.probability_matrix

    if "predicted_labels" in pred_df.columns:
        labels = pred_df["predicted_labels"].astype(str).to_numpy()
    else:
        labels = pred_df.iloc[:, 0].astype(str).to_numpy()
    classes = prob_df.columns.astype(str).to_numpy()
    probs = prob_df.to_numpy(dtype=float)
    # Align hard labels with probability argmax if CellTypist returns a mode-specific
    # column name unexpectedly.
    if len(labels) != probs.shape[0]:
        labels = classes[np.argmax(probs, axis=1)]
    return model, classes, labels, probs


def _prepare_celltypist_adata(adata):
    """Return a CellTypist-compatible AnnData copy.

    CellTypist expects nonnegative log1p-normalized values with max <= 9.22.
    The benchmark h5ads are model-facing matrices and may not satisfy that
    strict check, so the wrapper prepares an in-memory normalized copy.
    """
    import numpy as np
    import scanpy as sc

    a = adata.copy()
    try:
        sample = a.X[: min(1000, a.n_obs)]
        max_val = float(sample.max())
        min_val = float(sample.min())
    except Exception:
        max_val = 0.0
        min_val = 0.0

    if min_val < 0 or max_val > 9.22:
        # Treat the matrix as a nonnegative expression-like matrix when possible.
        # Negative values are clipped because CellTypist cannot consume them.
        try:
            import scipy.sparse as sp

            if sp.issparse(a.X):
                a.X = a.X.copy()
                a.X.data = np.clip(a.X.data, 0, None)
            else:
                a.X = np.clip(np.asarray(a.X), 0, None)
        except Exception:
            a.X = np.clip(np.asarray(a.X), 0, None)
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    return a


def main() -> int:
    args = build_parser().parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        meta = {
            "script": "build_classifier_milestone_model",
            "status": "dry_run",
            "dataset_id": args.dataset_id,
            "backend": args.backend,
            "input_h5ad": args.input_h5ad,
            "query_h5ad": args.query_h5ad,
            "label_key": args.label_key,
        }
        with open(out_dir / "classifier_milestone_model_metadata.json", "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=2)
        print("[build_classifier_milestone_model] dry-run metadata written")
        return 0

    import anndata as ad
    import joblib
    import numpy as np
    import pandas as pd
    from sklearn.metrics import classification_report, confusion_matrix
    from sklearn.model_selection import train_test_split

    print(f"[build_classifier_milestone_model] loading train {args.input_h5ad}")
    adata_ref = ad.read_h5ad(args.input_h5ad)
    adata_query = ad.read_h5ad(args.query_h5ad) if args.query_h5ad else adata_ref
    if args.label_key not in adata_ref.obs:
        raise KeyError(f"obs label key {args.label_key!r} not found")

    y_all = adata_ref.obs[args.label_key].astype(str).to_numpy()
    excluded = set(args.exclude_label or [])
    train_mask = ~np.isin(y_all, list(excluded))
    train_indices_all = np.where(train_mask)[0]
    sampled_rel = _balanced_indices(
        y_all[train_indices_all],
        max_per_class=args.max_train_cells_per_class,
        random_state=args.random_state,
    )
    train_indices = train_indices_all[sampled_rel]
    y = y_all[train_indices]
    if len(set(y.tolist())) < 2:
        raise ValueError("Need at least two milestone classes for classifier training")

    # Holdout validation uses the same backend family where practical.
    stratify = y if min(pd.Series(y).value_counts()) >= 2 else None
    tr_rel, va_rel = train_test_split(
        np.arange(len(train_indices)),
        test_size=0.2,
        random_state=args.random_state,
        stratify=stratify,
    )
    adata_train = adata_ref[train_indices[tr_rel]].copy()
    labels_train = y[tr_rel]
    adata_val = adata_ref[train_indices[va_rel]].copy()
    labels_val = y[va_rel]

    if args.backend == "celltypist":
        model, classes, val_labels, val_probs = _fit_predict_celltypist(
            adata_train, labels_train, adata_val, args
        )
        model_full, classes_full, labels_all_query, probs_all_query = _fit_predict_celltypist(
            adata_ref[train_indices].copy(), y, adata_query, args
        )
        model_to_save = model_full
        classes = classes_full
    else:
        model, classes, val_labels, val_probs = _fit_predict_sklearn(
            adata_train, labels_train, adata_val, args
        )
        model_full, classes, labels_all_query, probs_all_query = _fit_predict_sklearn(
            adata_ref[train_indices].copy(), y, adata_query, args
        )
        model_to_save = model_full

    conf_all = np.max(probs_all_query, axis=1)
    adata_query.obs["milestone_classifier_label"] = labels_all_query.astype(str)
    adata_query.obs["milestone_classifier_confidence"] = conf_all.astype("float32")
    for j, cls in enumerate(classes):
        adata_query.obs[f"classifier_prob_{cls}"] = probs_all_query[:, j].astype("float32")

    model_path = out_dir / f"classifier_milestone_model_{args.backend}.joblib"
    joblib.dump(model_to_save, model_path)
    pred_path = out_dir / "classifier_milestone_predictions.csv"
    _write_predictions_csv(
        pred_path,
        adata_query.obs_names,
        labels_all_query,
        conf_all,
        classes,
        probs_all_query,
    )

    report = classification_report(labels_val, val_labels, output_dict=True, zero_division=0)
    cm = confusion_matrix(labels_val, val_labels, labels=classes)
    cm_path = out_dir / "classifier_milestone_confusion_matrix.csv"
    pd.DataFrame(cm, index=classes, columns=classes).to_csv(cm_path)

    if args.output_h5ad:
        print(f"[build_classifier_milestone_model] writing {args.output_h5ad}")
        adata_query.write_h5ad(args.output_h5ad)

    meta: Dict = {
        "script": "build_classifier_milestone_model",
        "status": "completed",
        "provider_family": "classifier_based",
        "backend": args.backend,
        "dataset_id": args.dataset_id,
        "input_h5ad": args.input_h5ad,
        "query_h5ad": args.query_h5ad,
        "output_h5ad": args.output_h5ad,
        "training_label_key": args.label_key,
        "excluded_labels": sorted(excluded),
        "classes": list(map(str, classes)),
        "n_reference_cells": int(adata_ref.n_obs),
        "n_query_cells": int(adata_query.n_obs),
        "n_training_pool": int(train_mask.sum()),
        "n_training_sampled": int(len(train_indices)),
        "n_train": int(len(labels_train)),
        "n_validation": int(len(labels_val)),
        "validation_report": report,
        "model_path": str(model_path),
        "predictions_csv": str(pred_path),
        "confusion_matrix_csv": str(cm_path),
    }
    with open(out_dir / "classifier_milestone_model_metadata.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, default=_json_default)
        f.write("\n")

    print("[build_classifier_milestone_model] completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
