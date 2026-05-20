"""Trajectory Scenario A/B/C timepoint split preprocessor.

This formalizes the split rules that historically lived in method YAML files
and were reimplemented inside individual adapters.
"""

from __future__ import annotations

import numpy as np

from benchmark.shared.constants import ObservationColumns
from benchmark.shared.dataset.base import BaseDatasetPreprocessor


def _to_float_list(values):
    if values is None:
        return None
    return [float(v) for v in values]


def _materialize_subset(adata, mask):
    subset = adata[mask]
    if hasattr(subset, "to_memory") and getattr(subset, "isbacked", False):
        return subset.to_memory()
    return subset.copy()


def split_adata_by_timepoints(
    adata,
    *,
    time_key,
    train_times=None,
    heldout_times=None,
    test_includes_start=True,
):
    """Return train/test AnnData using trajectory's scenario time semantics."""
    if time_key not in adata.obs.columns:
        raise RuntimeError(
            f"ScenarioTimepointSplit requires adata.obs[{time_key!r}], "
            f"available columns: {list(adata.obs.columns)}"
        )

    obs_times = adata.obs[time_key].astype(float)
    all_times = sorted(float(t) for t in obs_times.unique())
    train_times_f = _to_float_list(train_times)
    heldout_times_f = _to_float_list(heldout_times)

    if train_times_f:
        train_set = set(train_times_f)
    elif heldout_times_f:
        heldout_set = set(heldout_times_f)
        train_set = {t for t in all_times if t not in heldout_set}
        train_times_f = sorted(train_set)
    else:
        train_set = set(all_times)
        train_times_f = list(all_times)

    if heldout_times_f:
        test_set = set(heldout_times_f)
    else:
        test_set = {t for t in all_times if t not in train_set}
        heldout_times_f = sorted(test_set)

    if test_includes_start and test_set:
        test_set.add(all_times[0])

    train_mask = obs_times.isin(train_set).to_numpy()
    if test_set:
        test_mask = obs_times.isin(test_set).to_numpy()
    else:
        test_mask = np.ones(adata.n_obs, dtype=bool)

    train_data = _materialize_subset(adata, train_mask)
    test_data = _materialize_subset(adata, test_mask)

    split_meta = {
        "split_preprocessor": "ScenarioTimepointSplit",
        "time_key": time_key,
        "train_times": train_times_f,
        "heldout_times": heldout_times_f,
        "test_includes_start": bool(test_includes_start),
        "n_train_cells": int(train_data.n_obs),
        "n_test_cells": int(test_data.n_obs),
    }
    train_data.uns["split_role"] = "train"
    test_data.uns["split_role"] = "test"
    train_data.uns["scenario_timepoint_split"] = split_meta
    test_data.uns["scenario_timepoint_split"] = split_meta
    return train_data, test_data


class ScenarioTimepointSplit(BaseDatasetPreprocessor):
    def __init__(
        self,
        dataset_dict,
        train_times=None,
        heldout_times=None,
        time_key=None,
        test_includes_start=True,
        **kwargs,
    ):
        super().__init__(dataset_dict)
        self.train_times = train_times
        self.heldout_times = heldout_times
        self.time_key = time_key
        self.test_includes_start = test_includes_start
        self.splits = True

    def _parameters(self):
        return {
            "train_times": self.train_times,
            "heldout_times": self.heldout_times,
            "time_key": self.time_key,
            "test_includes_start": self.test_includes_start,
        }

    def preprocess(self, ann_data, **kwargs):
        time_key = (
            self.time_key
            or self.dataset_dict.get("time_key")
            or ObservationColumns.TIMEPOINT.value
        )
        return split_adata_by_timepoints(
            ann_data,
            time_key=time_key,
            train_times=self.train_times,
            heldout_times=self.heldout_times,
            test_includes_start=self.test_includes_start,
        )
