"""Backfill scTimeBench-named embedding outputs for existing formal runs.

This is a compatibility helper for formal result directories produced before
the evaluator required method-side ``embedding.npy`` and
``next_timepoint_embedding.npy``. It reuses cached method artifacts instead of
retraining:

- scNODE: infer observed VAE latents from ``trained_scnode_model.pth``.
- MIOFlow: transform observed cells through the cached PCA in
  ``trained_mioflow_model.pth``.
- PRESCIENT: transform observed cells through the cached scaler/PCA in
  ``prescient_internal/data.pt``.

The projected embedding is already present in old runs as
``projected_embedding.npy``; this helper copies it to the scTimeBench-style
``next_timepoint_embedding.npy`` name when needed.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp
import torch
import yaml


def _find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env and Path(env).exists():
        return Path(env).resolve()
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "benchmark").is_dir():
            return parent
    raise RuntimeError("Cannot locate project root. Set TRAJ_PROJECT_ROOT.")


def _resolve(root: Path, path: str | Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else root / p


def _torch_load(path: Path, **kwargs: Any) -> Any:
    try:
        return torch.load(path, weights_only=False, **kwargs)
    except TypeError:
        return torch.load(path, **kwargs)


def _to_dense_float32(x: Any) -> np.ndarray:
    if sp.issparse(x):
        x = x.toarray()
    return np.asarray(x, dtype=np.float32)


def _load_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8-sig") as f:
        return yaml.safe_load(f) or {}


def _load_full_matrix(h5ad_path: Path) -> np.ndarray:
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError("anndata is required for embedding backfill") from exc
    adata = ad.read_h5ad(str(h5ad_path))
    return _to_dense_float32(adata.X)


def _backfill_mioflow(output_dir: Path, h5ad_path: Path) -> np.ndarray:
    payload_path = output_dir / "trained_mioflow_model.pth"
    if not payload_path.exists():
        raise FileNotFoundError(f"missing cached MIOFlow model: {payload_path}")
    payload = _torch_load(payload_path, map_location="cpu")
    if "pca" not in payload:
        raise KeyError(f"{payload_path} does not contain a cached PCA object")
    x_full = _load_full_matrix(h5ad_path)
    return payload["pca"].transform(x_full).astype(np.float32)


def _backfill_prescient(output_dir: Path, h5ad_path: Path) -> np.ndarray:
    data_path = output_dir / "prescient_internal" / "data.pt"
    if not data_path.exists():
        raise FileNotFoundError(f"missing cached PRESCIENT data.pt: {data_path}")
    data_pt = _torch_load(data_path, map_location="cpu")
    for key in ("scaler", "pca"):
        if key not in data_pt:
            raise KeyError(f"{data_path} does not contain {key!r}")
    x_full = _load_full_matrix(h5ad_path)
    x_scaled = data_pt["scaler"].transform(x_full)
    return data_pt["pca"].transform(x_scaled).astype(np.float32)


def _backfill_scnode(output_dir: Path, h5ad_path: Path, scnode_cfg: dict[str, Any]) -> np.ndarray:
    model_path = output_dir / "trained_scnode_model.pth"
    if not model_path.exists():
        raise FileNotFoundError(f"missing cached scNODE model: {model_path}")

    from benchmark.methods.scNODE.run import build_model

    x_full = _load_full_matrix(h5ad_path)
    model = build_model(x_full.shape[1], scnode_cfg)
    state = _torch_load(model_path, map_location="cpu")
    model.load_state_dict(state)
    model.eval()
    with torch.no_grad():
        latent_list, _ = model.vaeReconstruct([x_full])
    return latent_list[0].detach().cpu().numpy().astype(np.float32)


def backfill_from_config(config_path: Path, overwrite: bool = False) -> dict[str, Any]:
    root = _find_project_root()
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))

    cfg = _load_config(config_path)
    method = str(cfg.get("method", "")).lower()
    dataset = cfg.get("dataset") or {}
    output = cfg.get("output") or {}
    h5ad_path = _resolve(root, dataset.get("h5ad_path", ""))
    output_dir = _resolve(root, output.get("base_dir", ""))

    if method not in {"scnode", "mioflow", "prescient"}:
        raise ValueError(f"unsupported projection-capable method: {method!r}")
    if not h5ad_path.exists():
        raise FileNotFoundError(f"missing h5ad: {h5ad_path}")
    if not output_dir.exists():
        raise FileNotFoundError(f"missing output dir: {output_dir}")

    embedding_path = output_dir / "embedding.npy"
    next_path = output_dir / "next_timepoint_embedding.npy"
    projected_path = output_dir / "projected_embedding.npy"

    actions: list[str] = []
    if embedding_path.exists() and not overwrite:
        observed_shape = list(np.load(embedding_path, mmap_mode="r").shape)
        actions.append("embedding_exists")
    else:
        if method == "mioflow":
            observed = _backfill_mioflow(output_dir, h5ad_path)
        elif method == "prescient":
            observed = _backfill_prescient(output_dir, h5ad_path)
        else:
            observed = _backfill_scnode(output_dir, h5ad_path, cfg.get("scnode_params") or {})
        np.save(embedding_path, observed)
        observed_shape = list(observed.shape)
        actions.append("wrote_embedding")

    if next_path.exists() and not overwrite:
        next_shape = list(np.load(next_path, mmap_mode="r").shape)
        actions.append("next_timepoint_embedding_exists")
    else:
        if not projected_path.exists():
            raise FileNotFoundError(
                f"missing projected embedding fallback: {projected_path}"
            )
        shutil.copyfile(projected_path, next_path)
        next_shape = list(np.load(next_path, mmap_mode="r").shape)
        actions.append("copied_projected_to_next_timepoint")

    metadata = {
        "status": "completed",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "config": str(config_path),
        "method": method,
        "h5ad_path": str(h5ad_path),
        "output_dir": str(output_dir),
        "embedding_path": str(embedding_path),
        "embedding_shape": observed_shape,
        "next_timepoint_embedding_path": str(next_path),
        "next_timepoint_embedding_shape": next_shape,
        "projected_embedding_path": str(projected_path),
        "actions": actions,
    }
    (output_dir / "embedding_backfill_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Runtime YAML config")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing embedding.npy and next_timepoint_embedding.npy",
    )
    args = parser.parse_args()

    metadata = backfill_from_config(Path(args.config), overwrite=args.overwrite)
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
