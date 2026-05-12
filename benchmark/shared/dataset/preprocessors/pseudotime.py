"""scTimeBench-style pseudotime preprocessors.

This module keeps the pseudotime path intentionally close to scTimeBench:
select/prepare a representation, build a Scanpy neighbor graph, run diffusion
maps, choose a root cell from the earliest observed time point, and compute DPT.
"""

from __future__ import annotations

from enum import Enum
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc

from benchmark.shared.constants import ObservationColumns
from benchmark.shared.dataset.base import BaseDatasetPreprocessor


class PreprocessType(Enum):
    NONE = "none"
    PCA = "pca"
    HVG = "hvg"
    ZHENG_HVG = "zheng_hvg"
    OBSM = "obsm"


class DPTPseudotime(BaseDatasetPreprocessor):
    """Compute Scanpy DPT pseudotime and store it in obs.

    The implementation mirrors scTimeBench's ``Pseudotime`` preprocessor, with
    one trajectory-specific extension: ``preprocess_type: obsm`` can use an
    existing embedding such as ``X_scGPT`` before Scanpy computes neighbors.
    """

    def __init__(
        self,
        dataset_dict,
        preprocess_type="zheng_hvg",
        observed_time_key=None,
        pseudotime_key="dpt_pseudotime",
        n_neighbors=15,
        n_top_genes=1000,
        pca_components=50,
        obsm_key="X_scGPT",
        copy_to_timepoint=False,
        **kwargs,
    ):
        super().__init__(dataset_dict)
        self.preprocess_type = PreprocessType(preprocess_type)
        self.observed_time_key = observed_time_key or dataset_dict.get("time_key", "abs_day")
        self.pseudotime_key = pseudotime_key
        self.n_neighbors = int(n_neighbors)
        self.n_top_genes = int(n_top_genes)
        self.pca_components = int(pca_components)
        self.obsm_key = obsm_key
        self.copy_to_timepoint = bool(copy_to_timepoint)

    def requires_caching(self):
        return True

    def _parameters(self):
        return {
            "preprocess_type": self.preprocess_type.value,
            "observed_time_key": self.observed_time_key,
            "pseudotime_key": self.pseudotime_key,
            "n_neighbors": self.n_neighbors,
            "n_top_genes": self.n_top_genes,
            "pca_components": self.pca_components,
            "obsm_key": self.obsm_key,
            "copy_to_timepoint": self.copy_to_timepoint,
        }

    def _prepare_for_dpt(self, ann_data):
        if self.preprocess_type == PreprocessType.OBSM:
            if self.obsm_key not in ann_data.obsm:
                raise RuntimeError(
                    f"DPTPseudotime requested obsm[{self.obsm_key!r}], "
                    f"available keys: {list(ann_data.obsm.keys())}"
                )
            X = np.asarray(ann_data.obsm[self.obsm_key], dtype=np.float32)
            return ad.AnnData(X=X, obs=ann_data.obs.copy())

        if self.preprocess_type == PreprocessType.HVG:
            tmp = ann_data.copy()
            sc.pp.highly_variable_genes(tmp, n_top_genes=self.n_top_genes, inplace=True)
            return tmp[:, tmp.var.highly_variable].copy()

        if self.preprocess_type == PreprocessType.ZHENG_HVG:
            return sc.pp.recipe_zheng17(
                ann_data, n_top_genes=self.n_top_genes, copy=True
            )

        if self.preprocess_type == PreprocessType.PCA:
            tmp = ann_data.copy()
            sc.tl.pca(tmp, n_comps=self.pca_components)
            return ad.AnnData(X=tmp.obsm["X_pca"], obs=ann_data.obs.copy())

        return ann_data.copy()

    def preprocess(self, ann_data, checkpoint_dir=None, **kwargs):
        if self.observed_time_key not in ann_data.obs.columns:
            raise RuntimeError(
                f"DPTPseudotime needs obs[{self.observed_time_key!r}] for root-cell selection."
            )

        cache_path = None
        if checkpoint_dir is not None:
            cache_path = Path(checkpoint_dir) / f"{self.pseudotime_key}.npy"
            if cache_path.exists():
                pseudotime = np.load(cache_path)
                ann_data.obs[self.pseudotime_key] = pseudotime
                if self.copy_to_timepoint:
                    ann_data.obs[ObservationColumns.TIMEPOINT.value] = pseudotime
                return ann_data

        work = self._prepare_for_dpt(ann_data)
        use_rep = "X" if self.preprocess_type in (PreprocessType.OBSM, PreprocessType.PCA) else None
        sc.pp.neighbors(work, n_neighbors=self.n_neighbors, use_rep=use_rep)
        sc.tl.diffmap(work)

        earliest_tp = ann_data.obs[self.observed_time_key].astype(float).min()
        root_idx = ann_data.obs[
            ann_data.obs[self.observed_time_key].astype(float) == earliest_tp
        ].index[0]
        work.uns["iroot"] = work.obs_names.get_loc(root_idx)
        sc.tl.dpt(work)

        pseudotime = np.asarray(work.obs["dpt_pseudotime"], dtype=np.float32)
        ann_data.obs[self.pseudotime_key] = pseudotime
        if self.copy_to_timepoint:
            ann_data.obs[ObservationColumns.TIMEPOINT.value] = pseudotime

        if cache_path is not None:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(cache_path, pseudotime)
        return ann_data


class RoundPseudotimeToBins(BaseDatasetPreprocessor):
    """Discretize a continuous pseudotime column into ordered bins."""

    def __init__(
        self,
        dataset_dict,
        pseudotime_key="dpt_pseudotime",
        bin_key="dpt_pseudotime_bin",
        bin_numeric_key="dpt_pseudotime_bin_numeric",
        num_bins=15,
        copy_to_timepoint=True,
        **kwargs,
    ):
        super().__init__(dataset_dict)
        self.pseudotime_key = pseudotime_key
        self.bin_key = bin_key
        self.bin_numeric_key = bin_numeric_key
        self.num_bins = int(num_bins)
        self.copy_to_timepoint = bool(copy_to_timepoint)

    def _parameters(self):
        return {
            "pseudotime_key": self.pseudotime_key,
            "bin_key": self.bin_key,
            "bin_numeric_key": self.bin_numeric_key,
            "num_bins": self.num_bins,
            "copy_to_timepoint": self.copy_to_timepoint,
        }

    def preprocess(self, ann_data, **kwargs):
        if self.pseudotime_key not in ann_data.obs.columns:
            raise RuntimeError(
                f"RoundPseudotimeToBins needs obs[{self.pseudotime_key!r}]."
            )
        values = ann_data.obs[self.pseudotime_key].to_numpy(dtype=float)
        if not np.all(np.isfinite(values)):
            raise RuntimeError(f"Non-finite values found in {self.pseudotime_key!r}.")

        order = np.argsort(values, kind="mergesort")
        bins = np.empty(values.shape[0], dtype=np.int32)
        for bin_id, idx in enumerate(np.array_split(order, self.num_bins)):
            bins[idx] = bin_id

        ann_data.obs[self.bin_numeric_key] = bins.astype(float)
        ann_data.obs[self.bin_key] = np.array([f"PT_{b:02d}" for b in bins], dtype=object)
        if self.copy_to_timepoint:
            ann_data.obs[ObservationColumns.TIMEPOINT.value] = ann_data.obs[
                self.bin_numeric_key
            ]
        ann_data.uns["dpt_pseudotime_bins"] = {
            "pseudotime_key": self.pseudotime_key,
            "bin_key": self.bin_key,
            "bin_numeric_key": self.bin_numeric_key,
            "num_bins": self.num_bins,
            "binning": "equal_cell_count_by_sorted_dpt_pseudotime",
        }
        return ann_data
