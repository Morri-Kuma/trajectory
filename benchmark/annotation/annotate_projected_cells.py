"""annotate_projected_cells.py
================================
Step 9: Assign milestone labels to projected (simulated) cells.

Formal two-provider projected annotation (policy D)
----------------------------------------------------
two_provider_projection  (policy D -- production mode)
    Apply the two official observed-cell trained providers to projected outputs.
    Both providers consume projected_expression.npy.

        projected_expression.npy ->  HVG2000 PCA embedding provider
                                     (gene align -> pca.transform -> scaler ->
                                      classifier.predict_proba)
                                 ->  milestone_embedding_label
                                    milestone_embedding_confidence
                                    embedding_prob_<milestone> (per class)
                                    sidecar: projected_hvg_embedding.npy (50-dim PCA)

        projected_expression.npy ->  CellTypist classifier milestone model
                                 ->  milestone_classifier_label
                                    milestone_classifier_confidence
                                    classifier_prob_<milestone> (per class)

    Consensus contract (deterministic, written into metadata):
        - Agreement (emb == clf):
              consensus_milestone_label  = shared label
              consensus_confidence       = mean(emb_conf, clf_conf)
              consensus_status           = "agreement"
        - Disagreement:
              consensus_milestone_label  = "ambiguous"
              consensus_confidence       = min(emb_conf, clf_conf)
              consensus_status           = "disagreement_ambiguous"

    A formal complete run requires both providers to succeed
    (formal_two_provider_complete: true).  projected_embedding.npy
    stores the method-native latent space only and is not used here.

    Marker-score diagnostics may be included via --include-marker-diagnostics
    but are written ONLY as diag_marker_* columns and MUST NOT populate
    milestone_embedding_label, milestone_classifier_label, or
    consensus_milestone_label.

    Official models for GSE230659 are encoded in OFFICIAL_ANNOTATION_ASSETS
    below; the --embedding-model-path / --classifier-model-path flags default
    to those canonical paths when --dataset-id is GSE230659.

    annotation_policy: two_provider_projection

Legacy smoke policies (A / B / C)
----------------------------------
marker-smoke  (policy A)
    Marker-gene scoring on projected_expression.npy.  Smoke/diagnostic only.
    annotation_policy: marker_score_projection_smoke

cluster-majority-smoke  (policy B)
    Majority-vote from observed cluster->milestone map.  Smoke only.
    annotation_policy: cluster_majority_projection_smoke

observed-copy-smoke  (policy C)
    Sample observed cells and copy their milestone labels.  Smoke only.
    annotation_policy: observed_cell_copy_smoke

All policies write:
    projected_milestone_labels.csv
    projected_milestone_annotation_metadata.json
    projected_milestone_label_counts.csv

CLI
---
    # Production two-provider annotation (GSE230659 defaults):
    python -m benchmark.annotation.annotate_projected_cells \\
        --reference-h5ad benchmark/inputs/gse230659_milestone_real_sensitivity_hvg2000/GSE230659_real_sensitivity_milestone_HVG2000_benchmark_input.h5ad \\
        --run-output-dir benchmark/results/scnode/gse230659_milestone_consensus_A_hvg2000_formal \\
        --output-csv benchmark/results/scnode/gse230659_milestone_consensus_A_hvg2000_formal/projected_milestone_labels.csv \\
        --dataset-id GSE230659 \\
        --mode two_provider_projection

    # Legacy smoke:
    python -m benchmark.annotation.annotate_projected_cells \\
        --reference-h5ad <path> \\
        --run-output-dir <smoke_dir> \\
        --output-csv <smoke_dir>/projected_milestone_labels.csv \\
        --dataset-id GSE230659 \\
        --mode observed-copy-smoke \\
        --max-cells 2000 \\
        --exclude-label ambiguous
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Official annotation asset registry
# ---------------------------------------------------------------------------
# These paths are relative to the repository root and define the ONLY approved
# annotation inputs for the GSE230659 projected-cell two-provider workflow.
# They must not be changed without a corresponding model retraining step.

OFFICIAL_ANNOTATION_ASSETS: Dict[str, Any] = {
    "GSE230659": {
        "reference_h5ad": (
            "benchmark/inputs/gse230659_milestone_real_sensitivity_hvg2000"
            "/GSE230659_real_sensitivity_milestone_HVG2000_benchmark_input.h5ad"
        ),
        "embedding_model": (
            "benchmark/annotation_runs/gse230659_real_sensitivity"
            "/embedding_model_pca_hvg2000/hvg_pca_embedding_milestone_model.joblib"
        ),
        "classifier_model": (
            "benchmark/annotation_runs/gse230659_real_sensitivity"
            "/classifier_model_celltypist/classifier_milestone_model_celltypist.joblib"
        ),
        "embedding_provider_family": "embedding_based",
        "embedding_provider_method": "hvg2000_pca_logistic_regression",
        "embedding_provider_input": "projected_expression.npy",
        "gene_universe": "HVG2000",
        "classifier_provider_family": "classifier_based",
        "notes": (
            "GSE230659 formal projected-cell annotation uses the HVG2000 PCA "
            "embedding provider (gene align -> pca.transform -> scaler -> "
            "classifier.predict_proba) applied to projected_expression.npy.  "
            "The full-gene scGPT embedding provider is NOT used for "
            "HVG2000-only projected outputs.  Both trained provider models "
            "were built from the real-sensitivity observed h5ad."
        ),
    }
}

# Deterministic consensus rule -- written verbatim into every metadata file so
# that interpretation cannot drift between runs or scripts.
TWO_PROVIDER_CONSENSUS_RULE: Dict[str, Any] = {
    "rule_version": "v2",
    "embedding_provider": "hvg2000_pca_logistic_regression",
    "embedding_provider_input": "projected_expression.npy",
    "inputs": ["milestone_embedding_label", "milestone_classifier_label"],
    "cases": [
        {
            "condition": "milestone_embedding_label == milestone_classifier_label",
            "consensus_milestone_label": "shared_label",
            "consensus_confidence": "mean(milestone_embedding_confidence, milestone_classifier_confidence)",
            "consensus_status": "agreement",
        },
        {
            "condition": "labels disagree",
            "consensus_milestone_label": "ambiguous",
            "consensus_confidence": "min(milestone_embedding_confidence, milestone_classifier_confidence)",
            "consensus_status": "disagreement_ambiguous",
        },
    ],
    "note": (
        "Both providers consume projected_expression.npy.  The HVG2000 PCA "
        "embedding provider applies gene alignment -> pca.transform -> "
        "scaler.transform -> classifier.predict_proba.  A formal complete "
        "run requires both providers to succeed."
    ),
    "marker_score_role": (
        "Marker scores are diagnostic only.  They are written as diag_marker_* "
        "columns and MUST NOT populate milestone_embedding_label, "
        "milestone_classifier_label, or consensus_milestone_label."
    ),
}


# ---------------------------------------------------------------------------
# Label-mode column map
# ---------------------------------------------------------------------------

LABEL_MODES: Dict[str, str] = {
    "consensus": "consensus_milestone_label",
    "embedding_based": "milestone_embedding_label",
    "classifier_based": "milestone_classifier_label",
}

# Default thresholds (match build_marker_seed_labels.py defaults)
DEFAULT_MIN_TOP_SCORE = 0.0
DEFAULT_MIN_SCORE_MARGIN = 0.05

REQUIRED_OUTPUT_COLS = [
    "projected_cell_id",
    "timepoint",
    "milestone_embedding_label",
    "milestone_classifier_label",
    "consensus_milestone_label",
    "consensus_confidence",
]

# Extended column order for the two_provider_projection mode.
# cell_id is listed first (per spec); projected_cell_id follows for compat.
TWO_PROVIDER_OUTPUT_COLS = [
    "cell_id",
    "projected_cell_id",
    "timepoint",
    "milestone_embedding_label",
    "milestone_embedding_confidence",
    "milestone_classifier_label",
    "milestone_classifier_confidence",
    "consensus_milestone_label",
    "consensus_confidence",
    "consensus_status",
]


def _prob_col(milestone: str) -> str:
    return f"prob_{milestone}"


def _softmax_marker_probabilities(
    scores: Dict[str, Any],
    milestone_names: List[str],
    n_cells: int,
) -> np.ndarray:
    """Convert per-milestone marker scores into a probability matrix."""
    if not milestone_names:
        return np.zeros((n_cells, 0), dtype=np.float64)

    cols: List[np.ndarray] = []
    for ms in milestone_names:
        if ms in scores:
            arr = np.asarray(scores[ms], dtype=np.float64)
        else:
            arr = np.full(n_cells, np.nan, dtype=np.float64)
        cols.append(arr)
    score_mat = np.column_stack(cols)

    finite = np.isfinite(score_mat)
    if not finite.any():
        return np.full((n_cells, len(milestone_names)), 1.0 / len(milestone_names))

    row_all_missing = ~finite.any(axis=1)
    min_finite = float(np.nanmin(score_mat[finite]))
    filled = np.where(finite, score_mat, min_finite - 1.0)
    shifted = filled - np.max(filled, axis=1, keepdims=True)
    exp = np.exp(shifted)
    denom = exp.sum(axis=1, keepdims=True)
    probs = exp / np.where(denom == 0, 1.0, denom)
    if row_all_missing.any():
        probs[row_all_missing, :] = 1.0 / len(milestone_names)
    return probs


# ---------------------------------------------------------------------------
# h5py reference loader (avoids anndata null-encoding issue)
# ---------------------------------------------------------------------------

def _load_reference_h5py(h5ad_path: Path) -> Tuple[Dict[str, np.ndarray], np.ndarray, int]:
    """Return (obs_dict, gene_names, n_obs).  Loads obs and var/_index."""
    try:
        import h5py
    except ImportError as exc:
        raise ImportError("h5py required: pip install h5py") from exc

    obs_dict: Dict[str, np.ndarray] = {}
    gene_names: np.ndarray = np.array([], dtype=object)

    with h5py.File(h5ad_path, "r") as f:
        # --- obs ---
        for key in f["obs"].keys():
            obj = f["obs"][key]
            import h5py as _h5
            if isinstance(obj, _h5.Dataset):
                raw = obj[:]
                obs_dict[key] = np.array(
                    [v.decode() if isinstance(v, bytes) else v for v in raw]
                )
            elif isinstance(obj, _h5.Group) and "categories" in obj and "codes" in obj:
                cats = [
                    c.decode() if isinstance(c, bytes) else str(c)
                    for c in obj["categories"][:]
                ]
                codes = obj["codes"][:]
                obs_dict[key] = np.array([
                    cats[int(c)] if (c >= 0 and int(c) < len(cats)) else "NA"
                    for c in codes
                ])
        # --- var gene names ---
        if "var" in f and "_index" in f["var"]:
            raw_genes = f["var"]["_index"][:]
            gene_names = np.array(
                [g.decode() if isinstance(g, bytes) else str(g) for g in raw_genes]
            )

    n_obs = len(next(iter(obs_dict.values()))) if obs_dict else 0
    return obs_dict, gene_names, n_obs


def _load_projected_cluster_labels(csv_path: Path) -> List[Dict[str, str]]:
    with open(csv_path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------------------
# Two-provider helpers (Policy D)
# ---------------------------------------------------------------------------

def _load_pca_embedding_model(model_path: Path) -> Dict[str, Any]:
    """Load the HVG2000 PCA embedding model dict from a joblib file.

    The asset is a dict containing:
        pca           - fitted sklearn PCA
        scaler        - fitted StandardScaler (applied to PCA space)
        classifier    - fitted LogisticRegression
        var_names     - list of gene names (length n_features_in_)
        classes_      - list of class label strings
        n_features_in_ - int (2000 for HVG2000)
        n_components   - int (50)
        method        - str ('hvg2000_pca_logistic_regression')
    """
    try:
        import joblib
    except ImportError as exc:
        raise ImportError("joblib required: pip install joblib") from exc
    try:
        model = joblib.load(model_path)
    except Exception as exc:
        raise RuntimeError(
            f"Failed to load PCA embedding model from {model_path}: {exc}"
        ) from exc
    if not isinstance(model, dict):
        raise RuntimeError(
            f"Expected a dict from {model_path}; got {type(model).__name__}.  "
            "Ensure the correct hvg_pca_embedding_milestone_model.joblib is used."
        )
    for key in ("pca", "scaler", "classifier", "var_names", "classes_"):
        if key not in model:
            raise RuntimeError(
                f"PCA embedding model dict is missing key {key!r}.  "
                "Expected keys: pca, scaler, classifier, var_names, classes_."
            )
    return model


def _load_classifier_model(model_path: Path):
    """Load the CellTypist Model from a joblib file.

    Returns the model object.  The object exposes:
        model.features      - np.ndarray of gene names (length n_model_genes)
        model.scaler        - fitted StandardScaler
        model.classifier    - fitted LogisticRegression (classes_ available)
    """
    try:
        import joblib
    except ImportError as exc:
        raise ImportError("joblib required: pip install joblib") from exc
    try:
        model = joblib.load(model_path)
    except Exception as exc:
        raise RuntimeError(
            f"Failed to load classifier model from {model_path}: {exc}. "
            "Ensure celltypist is installed in the current environment."
        ) from exc
    # Verify expected attributes are present.
    for attr in ("features", "scaler", "classifier"):
        if not hasattr(model, attr):
            raise RuntimeError(
                f"Loaded classifier model is missing attribute {attr!r}. "
                f"Expected a CellTypist Model object; got {type(model).__name__}."
            )
    return model


def _predict_pca_embedding_provider(
    expr_2d: np.ndarray,
    expr_gene_names: np.ndarray,
    model: Dict[str, Any],
    batch_size: int = 4096,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], np.ndarray, Dict[str, Any]]:
    """Apply the HVG2000 PCA embedding provider to projected expression.

    Route: gene-align -> pca.transform -> scaler.transform -> classifier.predict_proba

    Parameters
    ----------
    expr_2d:
        Expression matrix, shape (n_cells, n_expr_genes).
    expr_gene_names:
        Gene names corresponding to columns of expr_2d.
    model:
        Dict loaded from hvg_pca_embedding_milestone_model.joblib.
    batch_size:
        Number of cells per PCA/predict batch.

    Returns
    -------
    labels, confidences, proba_matrix, classes, pca_embedding, gene_alignment_summary
        pca_embedding: shape (n_cells, n_components) -- saved as projected_hvg_embedding.npy
    """
    pca = model["pca"]
    scaler = model["scaler"]
    classifier = model["classifier"]
    model_var_names = np.asarray(model["var_names"])
    classes: List[str] = list(model["classes_"])
    n_components: int = int(model.get("n_components", pca.n_components_))

    # Gene alignment: expression columns -> model var_names order
    expr_col_indices, fill_mask, missing_genes = _build_gene_alignment_index(
        expr_gene_names, model_var_names
    )
    n_missing = int(fill_mask.sum())
    n_found = int(len(model_var_names)) - n_missing
    feature_order_exact = n_missing == 0
    gene_alignment_summary: Dict[str, Any] = {
        "n_model_genes": int(len(model_var_names)),
        "n_expr_genes": int(len(expr_gene_names)),
        "n_genes_matched": n_found,
        "n_genes_missing_in_expr": n_missing,
        "missing_genes_sample": missing_genes[:20],
        "extra_genes_in_expr": int(len(expr_gene_names)) - n_found,
        "feature_order_exact": feature_order_exact,
    }
    print(
        f"[annotate_projected_cells] PCA embedding gene alignment: "
        f"{n_found}/{len(model_var_names)} model genes found in expression "
        f"({n_missing} filled with 0.0)."
    )
    if n_found == 0:
        raise RuntimeError(
            "No PCA embedding model genes found in the expression matrix.  "
            "Check that projected_expression.npy and the reference h5ad use "
            "matching gene names."
        )

    n_cells = expr_2d.shape[0]
    n_classes = len(classes)
    all_proba = np.empty((n_cells, n_classes), dtype=np.float64)
    pca_embedding = np.empty((n_cells, n_components), dtype=np.float64)

    for start in range(0, n_cells, batch_size):
        end = min(start + batch_size, n_cells)
        batch = _align_expression_to_model(
            expr_2d[start:end], expr_col_indices, fill_mask
        )
        batch_pca = pca.transform(batch)              # (batch, n_components)
        pca_embedding[start:end] = batch_pca
        batch_scaled = scaler.transform(batch_pca)    # (batch, n_components)
        all_proba[start:end] = classifier.predict_proba(batch_scaled)

    predicted_idx = np.argmax(all_proba, axis=1)
    labels = np.array([classes[i] for i in predicted_idx], dtype=object)
    confidences = all_proba[np.arange(n_cells), predicted_idx]
    return labels, confidences, all_proba, classes, pca_embedding, gene_alignment_summary


def _build_gene_alignment_index(
    expr_gene_names: np.ndarray,
    model_features: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, List[str]]:
    """Map model feature genes to column indices in the expression matrix.

    Parameters
    ----------
    expr_gene_names:
        Gene names from the benchmark h5ad var axis (length n_expr_genes).
    model_features:
        Gene names the classifier model was trained on (length n_model_genes).

    Returns
    -------
    expr_col_indices:
        Int array of length n_model_genes.  expr_col_indices[j] is the
        column in the expression matrix for model gene j, or -1 if missing.
    fill_mask:
        Bool array of length n_model_genes; True where gene is missing.
    missing_genes:
        List of gene names not found in the expression matrix.
    """
    gene_to_col: Dict[str, int] = {
        str(g): i for i, g in enumerate(expr_gene_names)
    }
    expr_col_indices = np.empty(len(model_features), dtype=np.intp)
    fill_mask = np.zeros(len(model_features), dtype=bool)
    missing_genes: List[str] = []

    for j, gene in enumerate(model_features):
        col = gene_to_col.get(str(gene), -1)
        expr_col_indices[j] = col
        if col == -1:
            fill_mask[j] = True
            missing_genes.append(str(gene))

    return expr_col_indices, fill_mask, missing_genes


def _align_expression_to_model(
    expr_2d: np.ndarray,
    expr_col_indices: np.ndarray,
    fill_mask: np.ndarray,
) -> np.ndarray:
    """Reorder expression columns to match model feature order.

    Missing genes are filled with 0.0 (treated as absent / not detected).
    Returns shape (n_cells, n_model_genes) as float64.
    """
    n_cells = expr_2d.shape[0]
    n_model_genes = len(expr_col_indices)
    aligned = np.zeros((n_cells, n_model_genes), dtype=np.float64)

    present_mask = ~fill_mask
    present_model_cols = np.where(present_mask)[0]
    present_expr_cols = expr_col_indices[present_mask]
    aligned[:, present_model_cols] = expr_2d[:, present_expr_cols].astype(np.float64)
    return aligned


def _predict_classifier_provider(
    expr_2d: np.ndarray,
    expr_gene_names: np.ndarray,
    model,
    batch_size: int = 4096,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, List[str], Dict[str, Any]]:
    """Apply the CellTypist expression classifier to projected expression.

    Parameters
    ----------
    expr_2d:
        Expression matrix, shape (n_cells, n_expr_genes), already 2-D.
    expr_gene_names:
        Gene names corresponding to the columns of expr_2d.
    model:
        CellTypist Model object (with .features, .scaler, .classifier).
    batch_size:
        Number of cells to process per batch to limit peak memory.

    Returns
    -------
    labels, confidences, proba_matrix, classes, gene_alignment_summary
    """
    model_features: np.ndarray = np.asarray(model.features)
    expr_col_indices, fill_mask, missing_genes = _build_gene_alignment_index(
        expr_gene_names, model_features
    )
    n_missing = int(fill_mask.sum())
    n_found = int(len(model_features)) - n_missing
    gene_alignment_summary: Dict[str, Any] = {
        "n_model_genes": int(len(model_features)),
        "n_expr_genes": int(len(expr_gene_names)),
        "n_genes_matched": n_found,
        "n_genes_missing_in_expr": n_missing,
        "missing_genes_sample": missing_genes[:20],
    }
    print(
        f"[annotate_projected_cells] Classifier gene alignment: "
        f"{n_found}/{len(model_features)} model genes found in expression "
        f"({n_missing} filled with 0.0)."
    )
    if n_found == 0:
        raise RuntimeError(
            "No classifier model genes found in the expression matrix. "
            "Check that the expression file and reference h5ad use matching gene names."
        )

    n_cells = expr_2d.shape[0]
    n_classes = len(model.classifier.classes_)
    classes: List[str] = list(model.classifier.classes_)
    all_proba = np.empty((n_cells, n_classes), dtype=np.float64)

    for start in range(0, n_cells, batch_size):
        end = min(start + batch_size, n_cells)
        batch = _align_expression_to_model(
            expr_2d[start:end], expr_col_indices, fill_mask
        )
        batch_scaled = model.scaler.transform(batch)
        all_proba[start:end] = model.classifier.predict_proba(batch_scaled)

    predicted_idx = np.argmax(all_proba, axis=1)
    labels = np.array([classes[i] for i in predicted_idx], dtype=object)
    confidences = all_proba[np.arange(n_cells), predicted_idx]
    return labels, confidences, all_proba, classes, gene_alignment_summary


def _build_consensus_two_provider(
    emb_labels: np.ndarray,
    emb_conf: np.ndarray,
    clf_labels: np.ndarray,
    clf_conf: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Derive consensus from two providers using the documented rule.

    See TWO_PROVIDER_CONSENSUS_RULE for the specification encoded here.
    Both providers are always expected to produce results; neither is Optional.

    Returns
    -------
    consensus_labels, consensus_confidences, consensus_statuses  (all length n_cells)
    """
    n_cells = len(clf_labels)
    consensus_labels = np.empty(n_cells, dtype=object)
    consensus_conf = np.empty(n_cells, dtype=np.float64)
    consensus_status = np.empty(n_cells, dtype=object)

    for i in range(n_cells):
        e_lbl = str(emb_labels[i])
        c_lbl = str(clf_labels[i])
        e_cf = float(emb_conf[i])
        c_cf = float(clf_conf[i])

        if e_lbl == c_lbl:
            consensus_labels[i] = e_lbl
            consensus_conf[i] = (e_cf + c_cf) / 2.0
            consensus_status[i] = "agreement"
        else:
            consensus_labels[i] = "ambiguous"
            consensus_conf[i] = min(e_cf, c_cf)
            consensus_status[i] = "disagreement_ambiguous"

    return consensus_labels, consensus_conf, consensus_status


# ---------------------------------------------------------------------------
# Policy D: two-provider projected annotation (production)
# ---------------------------------------------------------------------------

def _policy_two_provider_projection(
    run_output_dir: Path,
    projected_expression_path: Optional[Path],
    reference_h5ad: Path,
    embedding_model_path: Path,
    classifier_model_path: Path,
    dataset_id: str,
    exclude_labels: List[str],
    max_cells: Optional[int],
    include_marker_diagnostics: bool,
    time_key: str = "abs_day",
    batch_size: int = 4096,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Apply HVG2000 PCA embedding + CellTypist classifier to projected outputs.

    Both providers consume projected_expression.npy.
    milestone_embedding_label and milestone_classifier_label are populated
    exclusively from the two trained provider models.  Marker scores, if
    requested, are written only as diag_marker_* diagnostic columns.

    The HVG2000 PCA embedding sidecar (projected_hvg_embedding.npy) and its
    metadata JSON are saved to run_output_dir alongside the annotation CSV.
    projected_embedding.npy (method-native latent) is not used here.
    """
    # --- Locate expression file ---
    expr_path = projected_expression_path or (run_output_dir / "projected_expression.npy")
    if not expr_path.exists():
        raise FileNotFoundError(
            f"projected_expression.npy not found at {expr_path}."
        )

    # --- Load expression array ---
    print(f"[annotate_projected_cells] Loading {expr_path} ...")
    expr_raw = np.load(expr_path, allow_pickle=False)
    print(f"[annotate_projected_cells] expression shape: {expr_raw.shape}")

    # --- Reshape expression to 2-D and derive timepoint labels ---
    if expr_raw.ndim == 3:
        n_tp, n_cells_per_tp, n_expr_genes = expr_raw.shape
        expr_2d = expr_raw.reshape(-1, n_expr_genes)
    elif expr_raw.ndim == 2:
        n_tp, n_cells_per_tp = 1, expr_raw.shape[0]
        n_expr_genes = expr_raw.shape[1]
        expr_2d = expr_raw
    else:
        raise ValueError(f"Unexpected expression array shape: {expr_raw.shape}")

    # --- Derive timepoint labels from h5ad unique values ---
    print(f"[annotate_projected_cells] Loading reference h5ad for gene names and timepoints ...")
    obs_dict, gene_names, n_obs = _load_reference_h5py(reference_h5ad)
    print(f"[annotate_projected_cells] n_obs={n_obs}, n_genes={len(gene_names)}")

    if time_key in obs_dict and expr_raw.ndim == 3:
        unique_tps = sorted(set(float(v) for v in obs_dict[time_key]))
        if len(unique_tps) == n_tp:
            tp_labels_per_cell = [str(unique_tps[t]) for t in range(n_tp) for _ in range(n_cells_per_tp)]
        else:
            tp_labels_per_cell = [str(t) for t in range(n_tp) for _ in range(n_cells_per_tp)]
    elif expr_raw.ndim == 2:
        tp_labels_per_cell = ["0"] * n_cells_per_tp
    else:
        tp_labels_per_cell = [str(t) for t in range(n_tp) for _ in range(n_cells_per_tp)]

    n_total = expr_2d.shape[0]

    # --- Optional cell subsetting ---
    if max_cells and n_total > max_cells:
        rng = np.random.default_rng(42)
        idx = np.sort(rng.choice(n_total, size=max_cells, replace=False))
        expr_2d = expr_2d[idx]
        tp_labels_per_cell = [tp_labels_per_cell[i] for i in idx]
        n_total = max_cells
    else:
        idx = np.arange(n_total)

    # --- Load models ---
    print(f"[annotate_projected_cells] Loading PCA embedding model: {embedding_model_path}")
    emb_model = _load_pca_embedding_model(embedding_model_path)

    print(f"[annotate_projected_cells] Loading classifier model: {classifier_model_path}")
    clf_model = _load_classifier_model(classifier_model_path)

    # --- Run HVG2000 PCA embedding provider ---
    print(f"[annotate_projected_cells] Running PCA embedding provider on {n_total} cells ...")
    (emb_labels, emb_conf, emb_proba, emb_classes,
     pca_embedding, emb_gene_align_summary) = _predict_pca_embedding_provider(
        expr_2d, gene_names, emb_model, batch_size=batch_size
    )

    # --- Save projected_hvg_embedding.npy sidecar ---
    hvg_emb_path = run_output_dir / "projected_hvg_embedding.npy"
    np.save(hvg_emb_path, pca_embedding)
    print(f"[annotate_projected_cells] Saved projected_hvg_embedding.npy shape={pca_embedding.shape} -> {hvg_emb_path}")

    # --- Save projected_hvg_embedding_metadata.json ---
    hvg_emb_meta: Dict[str, Any] = {
        "source_expression_path": str(expr_path),
        "source_expression_shape": list(expr_raw.shape),
        "flattened_projected_cell_count": n_total,
        "n_genes": n_expr_genes,
        "provider_model_path": str(embedding_model_path),
        "provider_method": str(emb_model.get("method", "hvg2000_pca_logistic_regression")),
        "n_components": int(emb_model.get("n_components", pca_embedding.shape[1])),
        "gene_universe": "HVG2000",
        "feature_order_matched_exactly": emb_gene_align_summary.get("feature_order_exact", False),
        "n_genes_missing_in_expr": emb_gene_align_summary.get("n_genes_missing_in_expr", 0),
        "missing_genes_sample": emb_gene_align_summary.get("missing_genes_sample", []),
        "extra_genes_in_expr": emb_gene_align_summary.get("extra_genes_in_expr", 0),
        "output_shape": list(pca_embedding.shape),
        "creation_timestamp": datetime.now(timezone.utc).isoformat(),
    }
    hvg_emb_meta_path = run_output_dir / "projected_hvg_embedding_metadata.json"
    with open(hvg_emb_meta_path, "w", encoding="utf-8") as f:
        json.dump(hvg_emb_meta, f, indent=2)
    print(f"[annotate_projected_cells] Saved projected_hvg_embedding_metadata.json -> {hvg_emb_meta_path}")

    # --- Run CellTypist classifier provider ---
    print(f"[annotate_projected_cells] Running classifier provider on {n_total} cells ...")
    clf_labels, clf_conf, clf_proba, clf_classes, clf_gene_align_summary = _predict_classifier_provider(
        expr_2d, gene_names, clf_model, batch_size=batch_size
    )

    # --- Build consensus ---
    print(f"[annotate_projected_cells] Building consensus labels ...")
    consensus_labels, consensus_conf, consensus_statuses = _build_consensus_two_provider(
        emb_labels, emb_conf, clf_labels, clf_conf,
    )

    # --- Optional marker-score diagnostics (NOT formal labels) ---
    marker_diag_cols: Dict[str, List[float]] = {}
    marker_overlap_summary: Dict[str, Any] = {}
    if include_marker_diagnostics:
        try:
            from benchmark.annotation.marker_utils import (
                build_gene_index_map,
                compute_mean_marker_scores,
                load_markers_yaml,
                resolve_milestone_genes,
            )
            markers_yaml_path = Path(__file__).parent / "milestone_markers.yaml"
            markers_data = load_markers_yaml(markers_yaml_path)
            ds_markers = markers_data.get(dataset_id, {})
            primary_milestones = ds_markers.get("primary_milestones", [])
            gene_index_map = build_gene_index_map(var_names=gene_names, use_var_index=True)
            milestone_gene_indices, marker_overlap_summary = resolve_milestone_genes(
                marker_config=markers_data,
                dataset_id=dataset_id,
                milestones=primary_milestones,
                gene_index_map=gene_index_map,
                min_overlap=1,
                strict=False,
            )
            scores = compute_mean_marker_scores(expr_2d, milestone_gene_indices)
            for ms in primary_milestones:
                col_key = f"diag_marker_score_{ms}"
                arr = scores.get(ms)
                if arr is not None:
                    marker_diag_cols[col_key] = [round(float(v), 6) for v in arr]
                else:
                    marker_diag_cols[col_key] = [float("nan")] * n_total
            print(
                f"[annotate_projected_cells] Marker diagnostics computed for "
                f"{len(primary_milestones)} milestones (diag_marker_score_* only)."
            )
        except Exception as exc:
            print(
                f"[annotate_projected_cells] WARNING: marker diagnostic failed: {exc}. "
                "Skipping marker diagnostics."
            )

    # --- Assemble rows ---
    rows: List[Dict[str, Any]] = []
    for i in range(n_total):
        row: Dict[str, Any] = {
            "cell_id": f"proj_{i}",
            "projected_cell_id": f"proj_{i}",  # backward-compat alias
            "timepoint": tp_labels_per_cell[i],
            # HVG2000 PCA embedding provider outputs
            "milestone_embedding_label": str(emb_labels[i]),
            "milestone_embedding_confidence": round(float(emb_conf[i]), 6),
            # CellTypist classifier provider outputs
            "milestone_classifier_label": str(clf_labels[i]),
            "milestone_classifier_confidence": round(float(clf_conf[i]), 6),
            # Consensus
            "consensus_milestone_label": str(consensus_labels[i]),
            "consensus_confidence": round(float(consensus_conf[i]), 6),
            "consensus_status": str(consensus_statuses[i]),
            # Bookkeeping
            "annotation_policy": "two_provider_projection",
            "label_is_placeholder": "False",
        }
        # Per-class embedding probabilities
        for j, cls in enumerate(emb_classes):
            row[f"embedding_prob_{cls}"] = round(float(emb_proba[i, j]), 8)

        # Per-class classifier probabilities
        for j, cls in enumerate(clf_classes):
            row[f"classifier_prob_{cls}"] = round(float(clf_proba[i, j]), 8)

        # Marker diagnostics (diagnostic columns only, clearly prefixed)
        for col_key, values in marker_diag_cols.items():
            row[col_key] = values[i]

        rows.append(row)

    # --- Build policy metadata ---
    policy_meta: Dict[str, Any] = {
        "annotation_policy": "two_provider_projection",
        "consensus_rule": TWO_PROVIDER_CONSENSUS_RULE,
        "official_annotation_assets": OFFICIAL_ANNOTATION_ASSETS.get(dataset_id, {}),
        "embedding_model_path": str(embedding_model_path),
        "embedding_provider_method": str(emb_model.get("method", "hvg2000_pca_logistic_regression")),
        "gene_universe": "HVG2000",
        "embedding_model_n_features": int(emb_model.get("n_features_in_", len(emb_model["var_names"]))),
        "embedding_model_n_components": int(emb_model.get("n_components", pca_embedding.shape[1])),
        "embedding_model_classes": list(emb_model["classes_"]),
        "embedding_provider_gene_alignment": emb_gene_align_summary,
        "projected_hvg_embedding_path": str(hvg_emb_path),
        "projected_hvg_embedding_shape": list(pca_embedding.shape),
        "classifier_model_path": str(classifier_model_path),
        "classifier_model_classes": clf_classes,
        "classifier_gene_alignment": clf_gene_align_summary,
        "expression_path": str(expr_path),
        "expression_shape_raw": list(expr_raw.shape),
        "n_cells_after_subset": n_total,
        "marker_diagnostics_included": include_marker_diagnostics,
        "marker_overlap_summary": marker_overlap_summary,
        "marker_score_column_prefix": "diag_marker_score_",
        "marker_score_role": TWO_PROVIDER_CONSENSUS_RULE["marker_score_role"],
        "time_key_used": time_key,
        "probability_columns_embedding": [f"embedding_prob_{c}" for c in emb_classes],
        "probability_columns_classifier": [f"classifier_prob_{c}" for c in clf_classes],
        "probability_basis_embedding": "hvg2000_pca_logistic_regression_predict_proba",
        "probability_basis_classifier": "celltypist_logistic_regression_predict_proba",
        # Both providers succeeded: formal_two_provider_complete is always True here.
        "formal_two_provider_complete": True,
    }
    return rows, policy_meta



# ---------------------------------------------------------------------------
# Policy A: marker-score on projected expression
# ---------------------------------------------------------------------------

def _policy_marker_smoke(
    run_output_dir: Path,
    projected_expression_path: Optional[Path],
    gene_names: np.ndarray,
    dataset_id: str,
    exclude_labels: List[str],
    projected_cluster_csv: Optional[Path],
    max_cells: Optional[int],
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Run marker scoring on projected_expression.npy."""
    from pathlib import Path as _Path
    import sys as _sys
    _sys.path.insert(0, str(_Path(__file__).resolve().parents[2]))
    from benchmark.annotation.marker_utils import (
        build_gene_index_map,
        compute_mean_marker_scores,
        assign_labels,
        load_markers_yaml,
        resolve_milestone_genes,
    )

    expr_path = projected_expression_path or (run_output_dir / "projected_expression.npy")
    if not expr_path.exists():
        raise FileNotFoundError(
            f"projected_expression.npy not found at {expr_path}. "
            "Use --projected-expression to specify the path explicitly, "
            "or switch to --mode cluster-majority-smoke."
        )

    print(f"[annotate_projected_cells] Loading {expr_path} ...")
    expr = np.load(expr_path, allow_pickle=False)
    print(f"[annotate_projected_cells] expression shape: {expr.shape}")

    if expr.ndim == 3:
        n_tp, n_cells_tp, n_genes = expr.shape
        expr_2d = expr.reshape(-1, n_genes)
    elif expr.ndim == 2:
        expr_2d = expr
        n_tp, n_cells_tp = 1, expr.shape[0]
    else:
        raise ValueError(f"Unexpected expression shape: {expr.shape}")

    if max_cells and expr_2d.shape[0] > max_cells:
        rng = np.random.default_rng(42)
        idx = rng.choice(expr_2d.shape[0], size=max_cells, replace=False)
        idx = np.sort(idx)
        expr_2d = expr_2d[idx]
    else:
        idx = np.arange(expr_2d.shape[0])

    n_total = expr_2d.shape[0]

    markers_yaml_path = Path(__file__).parent / "milestone_markers.yaml"
    markers_data = load_markers_yaml(markers_yaml_path)
    ds_markers = markers_data.get(dataset_id, {})
    primary_milestones = ds_markers.get("primary_milestones", [])

    gene_index_map = build_gene_index_map(var_names=gene_names, use_var_index=True)
    milestone_gene_indices, marker_overlap_summary = resolve_milestone_genes(
        marker_config=markers_data,
        dataset_id=dataset_id,
        milestones=primary_milestones,
        gene_index_map=gene_index_map,
        min_overlap=1,
        strict=False,
    )

    scores = compute_mean_marker_scores(expr_2d, milestone_gene_indices)
    milestone_probabilities = _softmax_marker_probabilities(
        scores=scores,
        milestone_names=primary_milestones,
        n_cells=n_total,
    )
    labels, top_scores, margins = assign_labels(
        scores,
        primary_milestones,
        min_top_score=DEFAULT_MIN_TOP_SCORE,
        min_score_margin=DEFAULT_MIN_SCORE_MARGIN,
    )

    if projected_cluster_csv and projected_cluster_csv.exists():
        cluster_rows = _load_projected_cluster_labels(projected_cluster_csv)
        if max_cells and len(cluster_rows) > max_cells:
            cluster_rows = [cluster_rows[i] for i in idx]
        timepoints = [r.get("projected_timepoint", "") for r in cluster_rows]
        cluster_labels = [r.get("projected_cluster_label", "") for r in cluster_rows]
    else:
        if expr.ndim == 3:
            timepoints = [str(t) for t in range(n_tp) for _ in range(n_cells_tp)]
            timepoints = [timepoints[i] for i in idx]
        else:
            timepoints = ["0"] * n_total
        cluster_labels = [""] * n_total

    rows = []
    for i in range(n_total):
        lbl = str(labels[i])
        conf = float(margins[i]) if not np.isnan(margins[i]) else 0.0
        row: Dict[str, Any] = {
            "projected_cell_id": f"proj_marker_{i}",
            "timepoint": timepoints[i],
            "milestone_marker_label": lbl,
            "milestone_embedding_label": lbl,  # placeholder copy
            "milestone_classifier_label": lbl,  # placeholder copy
            "consensus_milestone_label": lbl,
            "consensus_confidence": round(conf, 4),
            "annotation_policy": "marker_score_projection_smoke",
            "label_is_placeholder": "True",
            "label_source": "marker_score",
            "projected_cluster_label": cluster_labels[i] if cluster_labels else "",
        }
        for j, ms in enumerate(primary_milestones):
            row[_prob_col(ms)] = round(float(milestone_probabilities[i, j]), 8)
        rows.append(row)

    meta = {
        "annotation_policy": "marker_score_projection_smoke",
        "expression_path": str(expr_path),
        "n_genes": int(n_genes if expr.ndim == 3 else expr.shape[1]),
        "primary_milestones": primary_milestones,
        "min_top_score": DEFAULT_MIN_TOP_SCORE,
        "min_score_margin": DEFAULT_MIN_SCORE_MARGIN,
        "marker_overlap_summary": marker_overlap_summary,
        "probability_columns": [_prob_col(ms) for ms in primary_milestones],
        "probability_basis": "softmax_marker_scores",
        "n_probability_classes": int(len(primary_milestones)),
        "placeholder_warning": (
            "milestone_embedding_label and milestone_classifier_label are "
            "copies of milestone_marker_label (placeholder policy)."
        ),
    }
    return rows, meta


# ---------------------------------------------------------------------------
# Policy B: cluster-majority vote
# ---------------------------------------------------------------------------

def _policy_cluster_majority(
    run_output_dir: Path,
    projected_cluster_path: Optional[Path],
    obs_dict: Dict[str, np.ndarray],
    exclude_labels: List[str],
    obs_cluster_col: str = "scgpt_pseudostate_provisional",
    obs_milestone_col: str = "consensus_milestone_label",
    max_cells: Optional[int] = None,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Map projected clusters to milestone labels by majority vote on observed cells."""

    csv_path = projected_cluster_path or (run_output_dir / "projected_cluster_labels.csv")
    if not csv_path.exists():
        raise FileNotFoundError(
            f"projected_cluster_labels.csv not found at {csv_path}. "
            "Use --projected-cluster-labels to specify explicitly."
        )

    if obs_cluster_col not in obs_dict:
        raise ValueError(
            f"Reference h5ad does not have obs column {obs_cluster_col!r}. "
            f"Available: {sorted(obs_dict)}"
        )
    if obs_milestone_col not in obs_dict:
        raise ValueError(
            f"Reference h5ad does not have obs column {obs_milestone_col!r}."
        )

    obs_clusters = obs_dict[obs_cluster_col]
    obs_labels = obs_dict[obs_milestone_col]

    cluster_label_map: Dict[str, str] = {}
    cluster_vote_counts: Dict[str, Counter] = {}
    for cl, lbl in zip(obs_clusters.tolist(), obs_labels.tolist()):
        if lbl in exclude_labels:
            continue
        if cl not in cluster_vote_counts:
            cluster_vote_counts[cl] = Counter()
        cluster_vote_counts[cl][lbl] += 1

    for cl, counts in cluster_vote_counts.items():
        if counts:
            cluster_label_map[cl] = counts.most_common(1)[0][0]

    print(f"[annotate_projected_cells] cluster->milestone map ({obs_cluster_col}):")
    for cl in sorted(cluster_label_map):
        cnt = cluster_vote_counts[cl]
        print(f"    {cl} -> {cluster_label_map[cl]}  (votes: {dict(cnt.most_common())})")

    cluster_rows = _load_projected_cluster_labels(csv_path)
    print(f"[annotate_projected_cells] Loaded {len(cluster_rows)} projected cell rows")

    if max_cells and len(cluster_rows) > max_cells:
        rng = np.random.default_rng(42)
        idx = sorted(rng.choice(len(cluster_rows), size=max_cells, replace=False))
        cluster_rows = [cluster_rows[i] for i in idx]

    rows = []
    n_unmapped = 0
    for i, r in enumerate(cluster_rows):
        cl = r.get("projected_cluster_label", "")
        tp = r.get("projected_timepoint", "")
        lbl = cluster_label_map.get(cl, "ambiguous")
        if cl and cl not in cluster_label_map:
            n_unmapped += 1
        row: Dict[str, Any] = {
            "projected_cell_id": f"proj_cluster_{i}",
            "timepoint": tp,
            "milestone_marker_label": lbl,
            "milestone_embedding_label": lbl,   # placeholder copy
            "milestone_classifier_label": lbl,  # placeholder copy
            "consensus_milestone_label": lbl,
            "consensus_confidence": 0.5,
            "annotation_policy": "cluster_majority_projection_smoke",
            "label_is_placeholder": "True",
            "label_source": f"majority_vote_from_{obs_cluster_col}",
            "projected_cluster_label": cl,
        }
        rows.append(row)

    meta = {
        "annotation_policy": "cluster_majority_projection_smoke",
        "projected_cluster_labels_path": str(csv_path),
        "obs_cluster_col": obs_cluster_col,
        "obs_milestone_col": obs_milestone_col,
        "cluster_label_map": cluster_label_map,
        "n_unmapped_clusters": n_unmapped,
        "placeholder_warning": (
            "milestone_embedding_label and milestone_classifier_label are "
            "copies of cluster-majority consensus label (placeholder policy)."
        ),
    }
    return rows, meta


# ---------------------------------------------------------------------------
# Policy C: observed-cell copy smoke
# ---------------------------------------------------------------------------

def _policy_observed_copy(
    obs_dict: Dict[str, np.ndarray],
    exclude_labels: List[str],
    max_cells: int,
    time_key: str = "abs_day",
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Sample observed cells and copy their milestone labels."""

    for col in ["consensus_milestone_label", time_key]:
        if col not in obs_dict:
            raise ValueError(
                f"Reference h5ad missing required obs column {col!r}. "
                f"Available: {sorted(obs_dict)}"
            )

    consensus = obs_dict["consensus_milestone_label"]
    emb = obs_dict.get("milestone_embedding_label", consensus)
    clf = obs_dict.get("milestone_classifier_label", consensus)
    confidence = obs_dict.get("consensus_confidence")
    timepoints = obs_dict[time_key]
    cell_ids = obs_dict.get("_index", np.arange(len(consensus)).astype(str))

    all_labels = np.stack([consensus, emb, clf], axis=1)
    if exclude_labels:
        keep = ~np.any(np.isin(all_labels, exclude_labels), axis=1)
    else:
        keep = np.ones(len(consensus), dtype=bool)

    valid_idx = np.where(keep)[0]
    if len(valid_idx) == 0:
        raise ValueError("All cells excluded by exclude_labels. Cannot create smoke fixture.")

    unique_tps = np.unique(timepoints[valid_idx])
    n_per_tp = max(1, max_cells // len(unique_tps))
    sampled_indices: List[int] = []
    rng = np.random.default_rng(42)
    for tp in unique_tps:
        tp_idx = valid_idx[timepoints[valid_idx] == tp]
        n_take = min(n_per_tp, len(tp_idx))
        chosen = rng.choice(tp_idx, size=n_take, replace=False)
        sampled_indices.extend(chosen.tolist())
    sampled_indices = sorted(sampled_indices[:max_cells])

    rows = []
    for new_i, orig_i in enumerate(sampled_indices):
        lbl = str(consensus[orig_i])
        e_lbl = str(emb[orig_i])
        c_lbl = str(clf[orig_i])
        conf = float(confidence[orig_i]) if confidence is not None else 0.5
        src_id = (
            str(cell_ids[orig_i])
            if hasattr(cell_ids[orig_i], "__str__") else str(orig_i)
        )
        row: Dict[str, Any] = {
            "projected_cell_id": f"obs_copy_{new_i}",
            "timepoint": str(timepoints[orig_i]),
            "milestone_embedding_label": e_lbl,
            "milestone_classifier_label": c_lbl,
            "consensus_milestone_label": lbl,
            "consensus_confidence": round(conf, 4),
            "milestone_marker_label": lbl,
            "annotation_policy": "observed_cell_copy_smoke",
            "label_is_placeholder": "True",
            "label_source": "observed_cell_copy",
            "source_cell_id": src_id,
            "projected_cluster_label": "",
        }
        rows.append(row)

    meta = {
        "annotation_policy": "observed_cell_copy_smoke",
        "time_key": time_key,
        "n_original_valid_cells": int(len(valid_idx)),
        "placeholder_warning": (
            "All label columns are copied from observed reference cells. "
            "This is a smoke fixture only; labels are not derived from "
            "projected model outputs."
        ),
    }
    return rows, meta


# ---------------------------------------------------------------------------
# CSV + metadata writers
# ---------------------------------------------------------------------------

def _write_csv(rows: List[Dict[str, Any]], output_csv: Path) -> None:
    if not rows:
        raise ValueError("No rows to write to projected_milestone_labels.csv")
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    policy = rows[0].get("annotation_policy", "")
    if policy == "two_provider_projection":
        base_cols = TWO_PROVIDER_OUTPUT_COLS
    else:
        base_cols = REQUIRED_OUTPUT_COLS

    extra = [k for k in rows[0] if k not in base_cols]
    fieldnames = base_cols + extra
    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    print(f"[annotate_projected_cells] Wrote {len(rows)} rows -> {output_csv}")


def _write_metadata(
    meta: Dict[str, Any],
    output_csv: Path,
    reference_h5ad: str,
    run_output_dir: Optional[str],
    projected_expression_path: Optional[str],
    projected_embedding_path: Optional[str],
    projected_cluster_labels_path: Optional[str],
    dataset_id: str,
    mode: str,
    rows: List[Dict[str, Any]],
    exclude_labels: List[str],
) -> None:
    n_cells = len(rows)
    timepoints = sorted(set(r["timepoint"] for r in rows))
    label_counts = Counter(r["consensus_milestone_label"] for r in rows)
    is_smoke = mode in (
        "observed-copy-smoke", "observed_cell_copy_smoke",
        "cluster-majority-smoke", "cluster_majority_projection_smoke",
        "marker-smoke", "marker_score_projection_smoke",
    )
    meta_path = output_csv.parent / "projected_milestone_annotation_metadata.json"
    counts_path = output_csv.parent / "projected_milestone_label_counts.csv"

    label_modes_available = list(LABEL_MODES.keys())
    state_keys_available = list(LABEL_MODES.values())
    _ds_prefix = dataset_id.lower() if dataset_id else "unknown"
    provider_ids_available = [
        f"{_ds_prefix}_milestone_consensus_v1",
        f"{_ds_prefix}_milestone_embedding_v1",
        f"{_ds_prefix}_milestone_classifier_v1",
    ]

    is_two_provider = (meta.get("annotation_policy") == "two_provider_projection")
    label_cols_written = TWO_PROVIDER_OUTPUT_COLS if is_two_provider else REQUIRED_OUTPUT_COLS

    prob_cols_written: List[str] = []
    if is_two_provider:
        prob_cols_written = (
            meta.get("probability_columns_embedding", [])
            + meta.get("probability_columns_classifier", [])
        )
    else:
        prob_cols_written = meta.get("probability_columns", [])

    full_meta: Dict[str, Any] = {
        "dataset_id": dataset_id,
        "label_modes_available": label_modes_available,
        "state_keys_available": state_keys_available,
        "provider_ids_available": provider_ids_available,
        "reference_h5ad": reference_h5ad,
        "run_output_dir": run_output_dir,
        "projected_expression_path": projected_expression_path,
        "projected_embedding_path": projected_embedding_path,
        "projected_cluster_labels_path": projected_cluster_labels_path,
        "output_csv": str(output_csv),
        "n_projected_cells": n_cells,
        "timepoints": timepoints,
        "annotation_policy": meta.get("annotation_policy", mode),
        "smoke_test": is_smoke,
        "formal_benchmark": not is_smoke,
        "label_columns_written": label_cols_written,
        "probability_columns_written": prob_cols_written,
        "probability_basis": meta.get("probability_basis", ""),
        "n_probability_classes": meta.get("n_probability_classes"),
        "excluded_labels": exclude_labels,
        "label_counts_consensus": dict(label_counts),
        "placeholder_warning": meta.get("placeholder_warning", ""),
        "timestamp": datetime.now(timezone.utc).isoformat(),
        **{k: v for k, v in meta.items() if k not in (
            "annotation_policy", "placeholder_warning"
        )},
    }

    if is_two_provider:
        full_meta["consensus_rule"] = TWO_PROVIDER_CONSENSUS_RULE
        full_meta["official_annotation_assets"] = OFFICIAL_ANNOTATION_ASSETS.get(dataset_id, {})
        full_meta.pop("marker_copy_warning", None)
    else:
        full_meta["marker_copy_warning"] = (
            "milestone_embedding_label and milestone_classifier_label are "
            "placeholder copies; independent embedding/classifier models "
            "have not been trained."
        )

    meta_json_str = json.dumps(full_meta, indent=2, default=str)
    with open(meta_path, "w", encoding="utf-8") as f:
        f.write(meta_json_str)
    print(f"[annotate_projected_cells] Wrote metadata -> {meta_path}")

    with open(counts_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["label", "count"])
        writer.writeheader()
        for lbl, cnt in sorted(label_counts.items()):
            writer.writerow({"label": lbl, "count": cnt})
    print(f"[annotate_projected_cells] Wrote label counts -> {counts_path}")


# ---------------------------------------------------------------------------
# Path resolver for official model assets
# ---------------------------------------------------------------------------

def _resolve_official_model_path(
    dataset_id: str,
    asset_key: str,
    override_path: Optional[str],
    repo_root: Optional[Path] = None,
) -> Path:
    """Return model path from explicit override or official asset registry."""
    if override_path:
        return Path(override_path)
    assets = OFFICIAL_ANNOTATION_ASSETS.get(dataset_id)
    if assets is None:
        raise ValueError(
            f"No official annotation assets registered for dataset {dataset_id!r}. "
            "Provide --embedding-model-path / --classifier-model-path explicitly."
        )
    rel_path: str = assets[asset_key]
    if repo_root is None:
        repo_root = Path(__file__).resolve().parents[2]
    return repo_root / rel_path


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_annotation(
    reference_h5ad: str,
    run_output_dir: str,
    output_csv: str,
    dataset_id: str,
    mode: str = "cluster-majority-smoke",
    projected_expression: Optional[str] = None,
    projected_embedding: Optional[str] = None,
    projected_cluster_labels: Optional[str] = None,
    max_cells: Optional[int] = None,
    exclude_labels: Optional[List[str]] = None,
    obs_cluster_col: str = "scgpt_pseudostate_provisional",
    embedding_model_path: Optional[str] = None,
    classifier_model_path: Optional[str] = None,
    include_marker_diagnostics: bool = False,
    time_key: str = "abs_day",
    dry_run: bool = False,
) -> Path:
    ref_path = Path(reference_h5ad)
    out_csv = Path(output_csv)
    rout_dir = Path(run_output_dir) if run_output_dir else out_csv.parent
    exclude_labels = list(exclude_labels) if exclude_labels else []

    print(f"[annotate_projected_cells] mode={mode!r}, dataset={dataset_id!r}")
    print(f"[annotate_projected_cells] reference_h5ad={ref_path}")
    print(f"[annotate_projected_cells] output_csv={out_csv}")

    if dry_run:
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        meta_path = out_csv.parent / "projected_milestone_annotation_metadata.json"
        dry_meta: Dict[str, Any] = {"dry_run": True, "mode": mode, "dataset_id": dataset_id}
        if mode == "two_provider_projection":
            dry_meta["consensus_rule"] = TWO_PROVIDER_CONSENSUS_RULE
            dry_meta["official_annotation_assets"] = OFFICIAL_ANNOTATION_ASSETS.get(dataset_id, {})
        with open(meta_path, "w", encoding="utf-8") as f:
            json.dump(dry_meta, f, indent=2)
        print("[annotate_projected_cells] dry-run: metadata written, no annotation performed.")
        return out_csv

    proj_expr_path = Path(projected_expression) if projected_expression else None
    proj_clust_path = Path(projected_cluster_labels) if projected_cluster_labels else None

    if mode == "two_provider_projection":
        emb_model_path = _resolve_official_model_path(
            dataset_id, "embedding_model", embedding_model_path
        )
        clf_model_path = _resolve_official_model_path(
            dataset_id, "classifier_model", classifier_model_path
        )
        rows, policy_meta = _policy_two_provider_projection(
            run_output_dir=rout_dir,
            projected_expression_path=proj_expr_path,
            reference_h5ad=ref_path,
            embedding_model_path=emb_model_path,
            classifier_model_path=clf_model_path,
            dataset_id=dataset_id,
            exclude_labels=exclude_labels,
            max_cells=max_cells,
            include_marker_diagnostics=include_marker_diagnostics,
            time_key=time_key,
        )
    elif mode in ("marker-smoke", "marker_score_projection_smoke", "marker_score_projection"):
        print(f"[annotate_projected_cells] Loading reference h5ad ...")
        obs_dict, gene_names, n_obs = _load_reference_h5py(ref_path)
        print(f"[annotate_projected_cells] n_obs={n_obs}, n_genes={len(gene_names)}")
        rows, policy_meta = _policy_marker_smoke(
            run_output_dir=rout_dir,
            projected_expression_path=proj_expr_path,
            gene_names=gene_names,
            dataset_id=dataset_id,
            exclude_labels=exclude_labels,
            projected_cluster_csv=proj_clust_path or (rout_dir / "projected_cluster_labels.csv"),
            max_cells=max_cells,
        )
        if mode == "marker_score_projection":
            for row in rows:
                row["annotation_policy"] = "marker_score_projection"
            policy_meta["annotation_policy"] = "marker_score_projection"
            policy_meta["formal_annotation_note"] = (
                "Projected cells were annotated from projected_expression.npy "
                "with marker-score milestone assignment."
            )
    elif mode in ("cluster-majority-smoke", "cluster_majority_projection_smoke", "cluster_majority_projection"):
        print(f"[annotate_projected_cells] Loading reference h5ad ...")
        obs_dict, gene_names, n_obs = _load_reference_h5py(ref_path)
        print(f"[annotate_projected_cells] n_obs={n_obs}, n_genes={len(gene_names)}")
        rows, policy_meta = _policy_cluster_majority(
            run_output_dir=rout_dir,
            projected_cluster_path=proj_clust_path,
            obs_dict=obs_dict,
            exclude_labels=exclude_labels,
            obs_cluster_col=obs_cluster_col,
            max_cells=max_cells,
        )
        if mode == "cluster_majority_projection":
            for row in rows:
                row["annotation_policy"] = "cluster_majority_projection"
            policy_meta["annotation_policy"] = "cluster_majority_projection"
            policy_meta["formal_annotation_note"] = (
                "Projected cells annotated by projected-cluster majority vote."
            )
    elif mode in ("observed-copy-smoke", "observed_cell_copy_smoke"):
        print(f"[annotate_projected_cells] Loading reference h5ad ...")
        obs_dict, gene_names, n_obs = _load_reference_h5py(ref_path)
        print(f"[annotate_projected_cells] n_obs={n_obs}, n_genes={len(gene_names)}")
        rows, policy_meta = _policy_observed_copy(
            obs_dict=obs_dict,
            exclude_labels=exclude_labels,
            max_cells=max_cells or 2000,
        )
    else:
        raise ValueError(
            f"Unknown mode {mode!r}. "
            "Choose: two_provider_projection, marker-smoke, "
            "cluster-majority-smoke, observed-copy-smoke"
        )

    if not rows:
        raise ValueError("No projected cells annotated. Check inputs.")

    print(f"[annotate_projected_cells] Annotated {len(rows)} projected cells.")
    lc = Counter(r["consensus_milestone_label"] for r in rows)
    print(f"[annotate_projected_cells] Label counts: {dict(sorted(lc.items()))}")

    _write_csv(rows, out_csv)
    _write_metadata(
        meta=policy_meta,
        output_csv=out_csv,
        reference_h5ad=str(ref_path),
        run_output_dir=str(rout_dir) if run_output_dir else None,
        projected_expression_path=projected_expression,
        projected_embedding_path=projected_embedding,
        projected_cluster_labels_path=projected_cluster_labels,
        dataset_id=dataset_id,
        mode=mode,
        rows=rows,
        exclude_labels=exclude_labels,
    )
    return out_csv


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Step 9: Assign milestone labels to projected cells."
    )
    parser.add_argument("--reference-h5ad", required=True,
                        help="Path to milestone-annotated benchmark h5ad.")
    parser.add_argument("--run-output-dir", default=None,
                        help="Directory containing method projection outputs.")
    parser.add_argument("--output-csv", required=True,
                        help="Output path for projected_milestone_labels.csv.")
    parser.add_argument("--dataset-id", required=True,
                        help="Dataset identifier (e.g. GSE230659).")
    parser.add_argument("--mode",
                        choices=[
                            "two_provider_projection",
                            "marker-smoke", "marker_score_projection_smoke", "marker_score_projection",
                            "cluster-majority-smoke", "cluster_majority_projection_smoke", "cluster_majority_projection",
                            "observed-copy-smoke", "observed_cell_copy_smoke",
                        ],
                        default="cluster-majority-smoke",
                        help=(
                            "'two_provider_projection' is the production mode applying the two "
                            "trained observed-cell provider models to projected outputs. "
                            "Smoke modes are legacy fixtures."
                        ))
    parser.add_argument("--projected-expression", default=None,
                        help="Path to projected_expression.npy (overrides default location).")
    parser.add_argument("--projected-embedding", default=None,
                        help=(
                            "Path to projected_embedding.npy (method-native latent space).  "
                            "Not used by two_provider_projection; kept for legacy smoke modes."
                        ))
    parser.add_argument("--projected-cluster-labels", default=None,
                        help="Path to projected_cluster_labels.csv (overrides default).")
    parser.add_argument("--obs-cluster-col", default="scgpt_pseudostate_provisional",
                        help="Reference obs column for cluster majority vote.")
    parser.add_argument("--max-cells", type=int, default=None,
                        help="Maximum projected cells to annotate.")
    parser.add_argument("--exclude-label", action="append", dest="exclude_labels",
                        default=[], metavar="LABEL",
                        help="Exclude milestone label (repeatable).")
    parser.add_argument("--embedding-model-path", default=None,
                        help=(
                            "Path to hvg_pca_embedding_milestone_model.joblib.  "
                            "Defaults to the official GSE230659 HVG2000 PCA asset "
                            "when --dataset-id is GSE230659."
                        ))
    parser.add_argument("--classifier-model-path", default=None,
                        help=(
                            "Path to classifier_milestone_model_celltypist.joblib. "
                            "Defaults to the official GSE230659 asset when --dataset-id is GSE230659."
                        ))
    parser.add_argument("--include-marker-diagnostics", action="store_true",
                        help=(
                            "Compute marker-gene scores and write them as diag_marker_score_* columns. "
                            "These NEVER populate formal label columns."
                        ))
    # --allow-classifier-only-fallback is removed: both providers now consume
    # projected_expression.npy and are always expected to succeed.
    parser.add_argument("--time-key", default="abs_day",
                        help="obs key for timepoint values in the reference h5ad (default: abs_day).")
    parser.add_argument("--dry-run", action="store_true",
                        help="Parse arguments and validate paths; do not write outputs.")
    args = parser.parse_args(argv)

    try:
        run_annotation(
            reference_h5ad=args.reference_h5ad,
            run_output_dir=args.run_output_dir,
            output_csv=args.output_csv,
            dataset_id=args.dataset_id,
            mode=args.mode,
            projected_expression=args.projected_expression,
            projected_embedding=args.projected_embedding,
            projected_cluster_labels=args.projected_cluster_labels,
            obs_cluster_col=args.obs_cluster_col,
            max_cells=args.max_cells,
            exclude_labels=args.exclude_labels,
            embedding_model_path=args.embedding_model_path,
            classifier_model_path=args.classifier_model_path,
            include_marker_diagnostics=args.include_marker_diagnostics,
            time_key=args.time_key,
            dry_run=args.dry_run,
        )
    except Exception as exc:
        print(f"[annotate_projected_cells] ERROR: {exc}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
