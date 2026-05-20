"""scTimeBench-style dataset and preprocessor base classes.

This module intentionally mirrors scTimeBench's shared/dataset/base.py while
using trajectory package paths. Dataset-specific innovation stays in registry
classes and preprocessors; this file is just the common pipeline.
"""

from __future__ import annotations

import hashlib
import json
import os

from benchmark.shared.constants import DATASET_DIR, ObservationColumns


DATASET_PREPROCESSOR_REGISTRY = {}


def register_dataset_preprocessor(cls):
    DATASET_PREPROCESSOR_REGISTRY[cls.__name__] = cls
    return cls


class BaseDatasetPreprocessor:
    def __init__(self, dataset_dict):
        self.dataset_dict = dataset_dict
        self.splits = False

    def requires_caching(self):
        return False

    def __init_subclass__(cls):
        register_dataset_preprocessor(cls)

    def _parameters(self):
        return {}

    def preprocess(self, ann_data, **kwargs):
        raise NotImplementedError("Subclasses should implement this method.")


DATASET_REGISTRY = {}


def register_dataset(cls):
    DATASET_REGISTRY[cls.__name__] = cls
    return cls


class BaseDataset:
    def __init__(self, dataset_dict, dataset_preprocessors, output_dir):
        self.dataset_dict = dataset_dict
        self.dataset_preprocessors = dataset_preprocessors
        self.output_dir = output_dir
        self.TRAIN_PROCESSED_DATA_FILE = "train_processed_data.h5ad"
        self.TEST_PROCESSED_DATA_FILE = "test_processed_data.h5ad"

    def __init_subclass__(cls):
        register_dataset(cls)

    def _load_data(self):
        raise NotImplementedError("Subclasses should implement this method.")

    def requires_caching(self):
        return any(f.requires_caching() for f in self.dataset_preprocessors)

    def encode_preprocessors(self, i=None):
        preprocessors_to_encode = (
            self.dataset_preprocessors if i is None else self.dataset_preprocessors[:i]
        )
        preprocessor_names = [
            {"name": type(f).__name__, "parameters": f._parameters()}
            for f in preprocessors_to_encode
        ]
        return json.dumps(preprocessor_names, sort_keys=True)

    def encode_dataset_dict(self):
        blocklist = ["data_path", "h5ad_path", "requires_caching", "data_preprocessing_steps", "tag"]
        return json.dumps(
            {k: v for k, v in self.dataset_dict.items() if k not in blocklist},
            sort_keys=True,
        )

    def get_name(self):
        return self.dataset_dict["name"]

    def load_data(self):
        self._load_data()

        assert hasattr(self, "data"), "Dataset must have a 'data' attribute after loading."
        assert (
            ObservationColumns.CELL_TYPE.value in self.data.obs.columns
        ), f"Dataset must have '{ObservationColumns.CELL_TYPE.value}' in observation metadata."
        assert (
            ObservationColumns.TIMEPOINT.value in self.data.obs.columns
        ), f"Dataset must have '{ObservationColumns.TIMEPOINT.value}' in observation metadata."

        encountered_split = False
        for i, dataset_preprocessor in enumerate(self.dataset_preprocessors):
            if dataset_preprocessor.splits and encountered_split:
                raise ValueError(
                    "Multiple dataset preprocessors producing splits are not supported."
                )

            checkpoint_dir = self.get_checkpoint_dir(i)
            if dataset_preprocessor.splits:
                encountered_split = True
                train_data, test_data = dataset_preprocessor.preprocess(
                    self.data, checkpoint_dir=checkpoint_dir
                )
                self.data = (train_data, test_data)
            elif encountered_split:
                train_data = dataset_preprocessor.preprocess(
                    self.data[0], checkpoint_dir=checkpoint_dir
                )
                test_data = dataset_preprocessor.preprocess(
                    self.data[1], checkpoint_dir=checkpoint_dir
                )
                self.data = (train_data, test_data)
            else:
                self.data = dataset_preprocessor.preprocess(
                    self.data, checkpoint_dir=checkpoint_dir
                )

        assert encountered_split, "At least one dataset preprocessor must produce train-test splits."
        return self.data

    def __str__(self):
        return (
            f"Dataset Name: {self.get_name()}\n"
            f"Dataset Config: {self.dataset_dict}\n"
            f"Applied preprocessors: "
            f"{[type(f).__name__ + ', parameters: ' + str(f._parameters()) for f in self.dataset_preprocessors]}"
        )

    def get_dataset_dir(self):
        unique_string = json.dumps(
            {
                "dataset_dict": self.encode_dataset_dict(),
                "preprocessors": self.encode_preprocessors(),
            },
            sort_keys=True,
        )
        return os.path.join(
            self.output_dir,
            DATASET_DIR,
            hashlib.sha256(unique_string.encode()).hexdigest(),
        )

    def get_checkpoint_dir(self, i):
        unique_string = json.dumps(
            {
                "dataset_dict": self.encode_dataset_dict(),
                "preprocessors": self.encode_preprocessors(i),
            },
            sort_keys=True,
        )
        return os.path.join(
            self.output_dir,
            DATASET_DIR,
            "checkpoints",
            hashlib.sha256(unique_string.encode()).hexdigest(),
        )

    def create_dataset_dir(self):
        os.makedirs(self.get_dataset_dir(), exist_ok=True)
