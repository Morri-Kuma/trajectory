"""Shared constants for trajectory's scTimeBench-style benchmark layer."""

from enum import Enum


class ObservationColumns(Enum):
    CELL_TYPE = "scTimeBench_cell_type"
    TIMEPOINT = "scTimeBench_timepoint"


DATASET_DIR = "datasets"
PICKLED_DATASET_FILENAME = "dataset.pkl"
METHOD_CONFIG_FILENAME = "method_config.yaml"
