"""
eval_representation_run_from_config.py
======================================

Backfill helper that computes representation-space metrics for a completed
trajectory result directory, driven by the representation config that produced
it (Work-Plan Phase 4).

It is used for both Geneformer and scGPT representation result directories so
that both arms are scored through the exact same evaluator path:

    representation config (.yaml)  +  result_dir/projected_expression.npy
        -> representation_forecast_metrics.json
        -> representation_forecast_per_timepoint.csv
        -> representation_temporal_signal.json

These configs are representation-mode runs (``evaluation.mode: representation``
and ``representation.enabled: true``), so ``projected_expression.npy`` is the
projected *representation* (not gene-expression). This helper never renames or
mutates that file; it only reads it.

The official expression-space metrics (``forecast_metrics.json`` etc.) are left
untouched. Nothing here feeds the official silver report.

Example
-------
    python -m benchmark.evaluation.eval_representation_run_from_config \
        --config benchmark/configs/representation/mioflow_gse230659_rep_geneformer_cls_pca50_B.yaml \
        --result-dir benchmark/results/representation_dynamics/mioflow/gse230659_rep_geneformer_cls_pca50_B
"""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np

from benchmark.evaluation.eval_representation_forecast import (
    run_representation_forecast_evaluation,
)
from benchmark.evaluation.eval_representation_temporal_signal import (
    run_temporal_signal_evaluation,
)

DEFAULT_STATE_KEY = "final_milestone_label_coarse"
DEFAULT_TIME_KEY = "abs_day"
DEFAULT_PROJECTED_FILE = "projected_expression.npy"


# ---------------------------------------------------------------------------
# Config / path helpers
# ---------------------------------------------------------------------------

def load_config(config_path) -> dict:
    """Read a representation config YAML (tolerant of a UTF-8 BOM)."""
    import yaml

    with open(config_path, encoding="utf-8-sig") as f:
        return yaml.safe_load(f) or {}


def resolve_h5ad_path(config: dict, config_path) -> Path:
    """
    Resolve ``dataset.h5ad_path`` to an existing file.

    Tries, in order: the path as given (absolute or relative to CWD), then the
    path joined against each parent directory of the config file (so the repo
    root is found even when the helper is launched from elsewhere).
    """
    dataset = config.get("dataset", {}) or {}
    raw = dataset.get("h5ad_path")
    if not raw:
        raise KeyError("config is missing dataset.h5ad_path")
    candidates = [Path(raw)]
    cfg_path = Path(config_path).resolve()
    for parent in cfg_path.parents:
        candidates.append(parent / raw)
    for cand in candidates:
        if cand.exists():
            return cand
    raise FileNotFoundError(
        f"could not resolve dataset.h5ad_path={raw!r}; tried {len(candidates)} "
        f"locations including CWD and config-relative roots"
    )


def resolve_representation_id(config: dict, adata) -> Optional[str]:
    rep = config.get("representation", {}) or {}
    rid = rep.get("representation_id")
    if rid:
        return rid
    uns_rid = adata.uns.get("active_representation_id")
    return str(uns_rid) if uns_rid is not None else None


def resolve_state_key(config: dict, override: Optional[str]) -> str:
    if override:
        return override
    return (
        config.get("cell_state_key")
        or (config.get("lineage", {}) or {}).get("cell_state_key")
        or (config.get("ground_truth", {}) or {}).get("state_key")
        or DEFAULT_STATE_KEY
    )


def resolve_time_key(config: dict, override: Optional[str]) -> str:
    if override:
        return override
    return (
        (config.get("dataset", {}) or {}).get("time_key")
        or config.get("time_key")
        or DEFAULT_TIME_KEY
    )


def resolve_reference_graph(config: dict, override: Optional[str]) -> Optional[str]:
    if override:
        return override
    return (config.get("lineage", {}) or {}).get("reference_graph_path")


def resolve_eval_timepoints(config: dict, observed_times: np.ndarray) -> list:
    """
    Eval timepoints, in priority order:
      1. scenario_params.heldout_times
      2. observed unique times NOT in scenario_params.train_times
      3. all observed unique timepoints
    """
    sp = config.get("scenario_params", {}) or {}
    observed_unique = sorted({float(t) for t in np.asarray(observed_times, dtype=float)})

    heldout = sp.get("heldout_times")
    if heldout:
        return [float(t) for t in heldout]

    train = sp.get("train_times")
    if train:
        train_set = {float(t) for t in train}
        remaining = [t for t in observed_unique
                     if not any(np.isclose(t, ts) for ts in train_set)]
        if remaining:
            return remaining

    return observed_unique


# ---------------------------------------------------------------------------
# Metadata helpers
# ---------------------------------------------------------------------------

def _to_jsonable(obj):
    """Recursively convert numpy/anndata scalar & array types to plain Python."""
    if isinstance(obj, dict):
        return {str(k): _to_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [_to_jsonable(v) for v in obj.tolist()]
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def resolve_representation_metadata(adata, representation_id: Optional[str]) -> dict:
    """
    Assemble representation provenance for embedding in metric outputs.

    Pulls ``representation_inputs_metadata[representation_id]`` and, when an scFM
    embedding sidecar is present, nests the matching Geneformer/scGPT metadata
    under a flat ``scfm`` key. The report assembler expects
    ``representation_metadata["scfm"]["model_name"]`` and related fields.
    """
    meta: dict = {}
    reps = dict(adata.uns.get("representation_inputs_metadata", {}) or {})
    if representation_id and representation_id in reps:
        meta = dict(reps[representation_id])
    elif reps:
        # single-arm h5ad: take the only entry
        if len(reps) == 1:
            meta = dict(next(iter(reps.values())))
    if representation_id:
        meta.setdefault("representation_id", representation_id)

    scfm_by_model = dict(adata.uns.get("scfm_embedding_metadata", {}) or {})
    if scfm_by_model:
        chosen_key = None
        chosen_meta = None
        if representation_id:
            for model_name in ("geneformer", "scgpt", "scfoundation"):
                if model_name in representation_id and model_name in scfm_by_model:
                    chosen_key = model_name
                    chosen_meta = dict(scfm_by_model[model_name])
                    break
        if chosen_meta is None and len(scfm_by_model) == 1:
            chosen_key, chosen_meta_raw = next(iter(scfm_by_model.items()))
            chosen_meta = dict(chosen_meta_raw)
        if chosen_meta is not None:
            chosen_meta.setdefault("model_key", chosen_key)
            meta["scfm"] = chosen_meta
        else:
            meta["scfm_by_model"] = scfm_by_model

    return _to_jsonable(meta)


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def run_from_config(
    config_path,
    result_dir,
    *,
    backend: str = "auto",
    state_key: Optional[str] = None,
    time_key: Optional[str] = None,
    reference_graph: Optional[str] = None,
    projected_file: str = DEFAULT_PROJECTED_FILE,
    output_dir: Optional[str] = None,
    n_pred_cells: Optional[int] = None,
    seed: int = 0,
) -> dict:
    import anndata as ad

    config = load_config(config_path)
    result_dir = Path(result_dir)
    out_dir = Path(output_dir) if output_dir else result_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    h5ad_path = resolve_h5ad_path(config, config_path)
    adata = ad.read_h5ad(h5ad_path)

    if "X_rep" not in adata.obsm:
        raise KeyError(
            f"adata.obsm['X_rep'] required but not found in {h5ad_path}"
        )
    X_rep = np.asarray(adata.obsm["X_rep"], dtype=np.float64)

    rep_id = resolve_representation_id(config, adata)
    tkey = resolve_time_key(config, time_key)
    skey = resolve_state_key(config, state_key)
    ref_graph = resolve_reference_graph(config, reference_graph)

    if tkey not in adata.obs.columns:
        raise KeyError(f"time_key {tkey!r} not in adata.obs")
    observed_times = np.asarray(adata.obs[tkey].values, dtype=float)

    states = None
    if skey and skey in adata.obs.columns:
        states = np.asarray(adata.obs[skey].astype(str).values)

    eval_timepoints = resolve_eval_timepoints(config, observed_times)
    rep_meta = resolve_representation_metadata(adata, rep_id)

    projected_path = result_dir / projected_file
    if not projected_path.exists():
        raise FileNotFoundError(f"projected file not found: {projected_path}")

    # The forecast evaluator reads observed representation from a .npy path; stage
    # it in a temp dir so the result directory stays clean.
    with tempfile.TemporaryDirectory() as tmp:
        obs_rep_path = Path(tmp) / "observed_representation.npy"
        np.save(obs_rep_path, X_rep)
        forecast_metrics = run_representation_forecast_evaluation(
            projected_representation_path=str(projected_path),
            observed_representation_path=str(obs_rep_path),
            observed_times=observed_times,
            eval_timepoints=eval_timepoints,
            output_dir=str(out_dir),
            backend=backend,
            n_pred_cells=n_pred_cells,
            representation_metadata=rep_meta,
        )

    temporal_metrics = run_temporal_signal_evaluation(
        X_rep,
        observed_times,
        str(out_dir),
        states=states,
        reference_graph_path=ref_graph if (ref_graph and states is not None) else None,
        representation_id=rep_id,
        representation_metadata=rep_meta,
        seed=seed,
    )

    print(
        "[eval_representation_run_from_config] "
        f"rep_id={rep_id} h5ad={h5ad_path.name}\n"
        f"  forecast: status={forecast_metrics.get('status')} "
        f"backend={forecast_metrics.get('metric_backend')} "
        f"n_eval_tp={forecast_metrics.get('n_eval_timepoints')} "
        f"rep_wasserstein={forecast_metrics.get('rep_wasserstein_distance')}\n"
        f"  temporal: TVR={temporal_metrics.get('temporal_signal_ratio')} "
        f"state_sil={temporal_metrics.get('state_silhouette_score')}\n"
        f"  wrote -> {out_dir}/representation_forecast_metrics.json, "
        f"representation_forecast_per_timepoint.csv, "
        f"representation_temporal_signal.json"
    )
    return {"forecast": forecast_metrics, "temporal_signal": temporal_metrics}


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Backfill representation-space metrics from a config + result dir.",
    )
    p.add_argument("--config", required=True, help="Representation config YAML.")
    p.add_argument("--result-dir", required=True,
                   help="Result directory containing projected_expression.npy.")
    p.add_argument("--backend", default="auto", choices=["auto", "numpy", "geomloss"])
    p.add_argument("--state-key", default=None,
                   help="Override obs state-label column (default from config).")
    p.add_argument("--time-key", default=None,
                   help="Override obs time column (default from config).")
    p.add_argument("--reference-graph", default=None,
                   help="Override reference lineage graph JSON (default from config).")
    p.add_argument("--projected-file", default=DEFAULT_PROJECTED_FILE,
                   help="Projected representation .npy in the result dir.")
    p.add_argument("--output-dir", default=None,
                   help="Where to write metrics (default: result dir).")
    p.add_argument("--n-pred-cells", type=int, default=None,
                   help="Cells per timepoint for a 2D projected array.")
    p.add_argument("--seed", type=int, default=0)
    return p


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    run_from_config(
        args.config,
        args.result_dir,
        backend=args.backend,
        state_key=args.state_key,
        time_key=args.time_key,
        reference_graph=args.reference_graph,
        projected_file=args.projected_file,
        output_dir=args.output_dir,
        n_pred_cells=args.n_pred_cells,
        seed=args.seed,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
