import random

from benchmark.shared.constants import ObservationColumns
from benchmark.shared.dataset.base import BaseDatasetPreprocessor


class PercentileSplitTrainTest(BaseDatasetPreprocessor):
    def __init__(self, dataset_dict, train_pct, time_key=None, seed=None):
        super().__init__(dataset_dict)
        self.train_pct = train_pct
        self.time_key = time_key
        self.seed = seed
        self.splits = True

    def _parameters(self):
        return {
            "train_pct": self.train_pct,
            "time_key": self.time_key,
            "seed": self.seed,
        }

    def preprocess(self, ann_data, **kwargs):
        rng = random.Random(self.seed)
        train_indices = []
        test_indices = []
        tp_column = self.time_key or ObservationColumns.TIMEPOINT.value
        for tp in ann_data.obs[tp_column].unique():
            tp_data = ann_data[ann_data.obs[tp_column] == tp]
            n_train = int(tp_data.n_obs * self.train_pct)
            indices = list(range(tp_data.n_obs))
            rng.shuffle(indices)
            train_indices.extend(tp_data.obs.index[indices[:n_train]])
            test_indices.extend(tp_data.obs.index[indices[n_train:]])
        train_data = ann_data[train_indices].copy()
        test_data = ann_data[test_indices].copy()
        train_data.uns["split_role"] = "train"
        test_data.uns["split_role"] = "test"
        return train_data, test_data
