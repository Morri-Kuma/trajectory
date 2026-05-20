from __future__ import annotations

from pathlib import Path

import anndata as ad

from benchmark.shared.constants import ObservationColumns
from benchmark.shared.dataset.base import BaseDataset


class TrajectoryH5ADDataset(BaseDataset):
    """Base registry class for trajectory h5ad benchmark inputs."""

    dataset_id = "unknown"
    default_time_key = "abs_day"
    default_cell_state_key = "scTimeBench_cell_type"

    def _resolve_h5ad_path(self):
        raw_path = (
            self.dataset_dict.get("h5ad_path")
            or self.dataset_dict.get("data_path")
            or self.dataset_dict.get("path")
        )
        if raw_path is None:
            raise ValueError(
                f"{type(self).__name__} requires h5ad_path or data_path in dataset config."
            )
        path = Path(raw_path)
        if not path.is_absolute():
            project_root = Path(self.dataset_dict.get("project_root", "."))
            path = project_root / path
        if not path.exists():
            raise FileNotFoundError(f"h5ad file not found: {path}")
        return path

    def _load_data(self):
        h5ad_path = self._resolve_h5ad_path()
        backed = self.dataset_dict.get("backed")
        self.data = ad.read_h5ad(h5ad_path, backed=backed)
        self._standardize_obs_columns()
        self.data.uns["dataset_id"] = self.dataset_dict.get("id", self.dataset_id)
        self.data.uns["dataset_registry"] = type(self).__name__
        self.data.uns["source_h5ad_path"] = str(h5ad_path)

    def _standardize_obs_columns(self):
        time_key = self.dataset_dict.get("time_key", self.default_time_key)
        cell_state_key = self.dataset_dict.get(
            "cell_state_key",
            self.dataset_dict.get("cell_type_key", self.default_cell_state_key),
        )

        if ObservationColumns.TIMEPOINT.value not in self.data.obs.columns:
            if time_key not in self.data.obs.columns:
                raise RuntimeError(
                    f"{type(self).__name__} cannot create "
                    f"{ObservationColumns.TIMEPOINT.value!r}: obs[{time_key!r}] is missing."
                )
            self.data.obs[ObservationColumns.TIMEPOINT.value] = self.data.obs[time_key]

        if ObservationColumns.CELL_TYPE.value not in self.data.obs.columns:
            fallback_keys = [
                cell_state_key,
                "scgpt_pseudostate_provisional",
                "cell_type",
                "stage_day_label",
                "stage",
            ]
            for key in fallback_keys:
                if key in self.data.obs.columns:
                    self.data.obs[ObservationColumns.CELL_TYPE.value] = self.data.obs[key].astype(str)
                    break
            else:
                raise RuntimeError(
                    f"{type(self).__name__} cannot create "
                    f"{ObservationColumns.CELL_TYPE.value!r}; tried {fallback_keys}."
                )


class GSE230659Dataset(TrajectoryH5ADDataset):
    dataset_id = "GSE230659"
    default_time_key = "abs_day"
    default_cell_state_key = "scgpt_pseudostate_provisional"


class GSE178325Dataset(TrajectoryH5ADDataset):
    dataset_id = "GSE178325"
    default_time_key = "abs_day"
    default_cell_state_key = "scgpt_pseudostate_provisional"


class GSE242424Dataset(TrajectoryH5ADDataset):
    """GSE242424 iPSC reprogramming time-course (D0-D14 + iPSC).

    cell_state_key is a time/stage proxy (scTimeBench_cell_type), not a
    biologically annotated milestone system.  lineage_fidelity evaluation
    is deferred until a reference graph is available.
    """

    dataset_id = "GSE242424"
    default_time_key = "abs_day"
    default_cell_state_key = "scTimeBench_cell_type"
