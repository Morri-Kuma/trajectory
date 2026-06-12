"""
model_input.py
==============

Shared helper for representation-aware method runners (Work-Plan Step 4).

``get_model_input_matrix(adata, cfg)`` returns the matrix a trajectory model
should train on:

  * If ``cfg['representation']`` is enabled with ``input_mode == 'obsm'``, the
    pre-built representation in ``adata.obsm[obsm_key]`` is used (e.g. X_rep).
  * Otherwise the existing expression matrix ``adata.X`` is returned unchanged,
    preserving backward compatibility with all existing expression-space runs.

The representation config block is intentionally small and is validated by
:func:`parse_representation_config` so every runner interprets it identically.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

# Accepted values for representation.reducer_fit_scope (mirrors build module).
TRAIN_ONLY = "train_only"
ALL_CELLS_EXPLORATORY = "all_cells_exploratory"

_VALID_INPUT_MODES = ("obsm", "X")
_VALID_FIT_SCOPES = (TRAIN_ONLY, ALL_CELLS_EXPLORATORY)


@dataclass(frozen=True)
class RepresentationConfig:
    enabled: bool = False
    input_mode: str = "X"
    obsm_key: str = "X_rep"
    representation_id: Optional[str] = None
    input_space: Optional[str] = None
    output_space: Optional[str] = None
    final_dim: Optional[int] = None
    reducer_fit_scope: Optional[str] = None

    @property
    def uses_obsm(self) -> bool:
        return self.enabled and self.input_mode == "obsm"


def parse_representation_config(cfg: Optional[dict]) -> RepresentationConfig:
    """
    Parse and validate the ``representation`` block of a method config.

    A missing or empty block yields a disabled config (expression-space, the
    default), so existing configs keep working untouched.
    """
    rep = {} if not cfg else (cfg.get("representation") or {})
    if not rep:
        return RepresentationConfig()

    enabled = bool(rep.get("enabled", False))
    input_mode = rep.get("input_mode", "X")
    if input_mode not in _VALID_INPUT_MODES:
        raise ValueError(
            f"representation.input_mode must be one of {_VALID_INPUT_MODES}, "
            f"got {input_mode!r}."
        )

    fit_scope = rep.get("reducer_fit_scope")
    if fit_scope is not None and fit_scope not in _VALID_FIT_SCOPES:
        raise ValueError(
            f"representation.reducer_fit_scope must be one of {_VALID_FIT_SCOPES}, "
            f"got {fit_scope!r}."
        )

    final_dim = rep.get("final_dim")
    if final_dim is not None:
        final_dim = int(final_dim)

    return RepresentationConfig(
        enabled=enabled,
        input_mode=input_mode,
        obsm_key=rep.get("obsm_key", "X_rep"),
        representation_id=rep.get("representation_id"),
        input_space=rep.get("input_space"),
        output_space=rep.get("output_space"),
        final_dim=final_dim,
        reducer_fit_scope=fit_scope,
    )


def get_model_input_matrix(adata, cfg: Optional[dict]):
    """
    Return the model input matrix for ``adata`` given a method config ``cfg``.

    Representation mode (cfg.representation.enabled and input_mode == 'obsm'):
        returns ``np.asarray(adata.obsm[obsm_key], dtype=float32)``.
    Otherwise:
        returns ``adata.X`` unchanged (backward compatible; may be sparse).

    Raises a clear error if representation mode is requested but the obsm key is
    absent — never silently falls back to expression space, which would
    invalidate the representation comparison.
    """
    rep = parse_representation_config(cfg)
    if not rep.uses_obsm:
        return adata.X

    if rep.obsm_key not in adata.obsm:
        raise KeyError(
            f"representation.enabled but adata.obsm['{rep.obsm_key}'] is missing. "
            f"Build it with build_representation_inputs.py first. "
            f"Available obsm keys: {list(adata.obsm.keys())}."
        )
    mat = np.asarray(adata.obsm[rep.obsm_key], dtype=np.float32)

    if rep.final_dim is not None and mat.shape[1] != rep.final_dim:
        raise ValueError(
            f"representation.final_dim={rep.final_dim} but "
            f"adata.obsm['{rep.obsm_key}'] has {mat.shape[1]} columns. "
            "Rebuild the representation with matching n_components."
        )
    return mat


def representation_is_active(cfg: Optional[dict]) -> bool:
    """Convenience predicate: True iff cfg enables obsm representation input."""
    return parse_representation_config(cfg).uses_obsm


def representation_input_adata(adata, cfg: Optional[dict]):
    """
    Return an AnnData whose ``.X`` is the model input matrix for ``cfg``.

    This is the uniform, low-churn way to make an existing expression-space
    runner train in representation space without editing the runner internals:
    the adapter swaps the AnnData *once* before invoking the pipeline.

      * Representation disabled  -> the original ``adata`` is returned unchanged
        (backward compatible; zero behavioural change).
      * Representation enabled    -> a NEW AnnData is returned with
        ``X = adata.obsm[obsm_key]`` (dense float32), ``var`` reindexed to
        ``<representation_id>_0 .. _{d-1}``, and ``obs`` / ``obsm`` / ``uns``
        preserved so downstream code (time/state keys, metadata) still works.

    The returned object's ``.uns['model_input_space']`` records whether the
    matrix is ``expression`` or ``representation`` for traceability.
    """
    rep = parse_representation_config(cfg)
    if not rep.uses_obsm:
        try:
            adata.uns["model_input_space"] = "expression"
        except Exception:
            pass
        return adata

    mat = get_model_input_matrix(adata, cfg)  # validates obsm presence/dims

    import anndata as ad
    import pandas as pd

    prefix = rep.representation_id or rep.obsm_key
    var = pd.DataFrame(index=[f"{prefix}_{i}" for i in range(mat.shape[1])])
    new = ad.AnnData(
        X=np.asarray(mat, dtype=np.float32),
        obs=adata.obs.copy(),
        var=var,
    )
    # Preserve obsm/uns so timepoint/state metadata and provenance survive.
    for k in adata.obsm.keys():
        try:
            new.obsm[k] = adata.obsm[k]
        except Exception:
            pass
    for k in adata.uns.keys():
        try:
            new.uns[k] = adata.uns[k]
        except Exception:
            pass
    new.uns["model_input_space"] = "representation"
    new.uns["active_representation_id"] = rep.representation_id
    return new
