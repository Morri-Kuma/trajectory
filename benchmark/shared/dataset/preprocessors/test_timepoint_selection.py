from benchmark.shared.constants import ObservationColumns
from benchmark.shared.dataset.base import BaseDatasetPreprocessor


class TestTimepointSelection(BaseDatasetPreprocessor):
    def __init__(self, dataset_dict, test_tps, time_key=None, use_time_indices=None, **kwargs):
        super().__init__(dataset_dict)
        self.test_tps = test_tps
        self.time_key = time_key
        self.use_time_indices = (
            dataset_dict.get("use_time_indices", False)
            if use_time_indices is None
            else use_time_indices
        )
        self.splits = True

    def _parameters(self):
        return {
            "test_tps": self.test_tps,
            "time_key": self.time_key,
            "use_time_indices": self.use_time_indices,
        }

    def preprocess(self, ann_data, **kwargs):
        tp_column = self.time_key or ObservationColumns.TIMEPOINT.value
        timepoints = sorted(ann_data.obs[tp_column].unique())
        start_tp = timepoints[0]
        test_tps = list(self.test_tps)
        if self.use_time_indices:
            test_tps = [timepoints[tp] for tp in self.test_tps]

        train_tps = [tp for tp in timepoints if tp not in test_tps]
        if start_tp not in test_tps:
            test_tps.append(start_tp)

        train_data = ann_data[ann_data.obs[tp_column].isin(train_tps)].copy()
        test_data = ann_data[ann_data.obs[tp_column].isin(test_tps)].copy()
        train_data.uns["split_role"] = "train"
        test_data.uns["split_role"] = "test"
        train_data.uns["train_times"] = list(train_tps)
        test_data.uns["heldout_times"] = list(test_tps)
        return train_data, test_data
