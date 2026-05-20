from benchmark.shared.dataset.base import BaseDatasetPreprocessor


class CopyTrainTest(BaseDatasetPreprocessor):
    def __init__(self, dataset_dict):
        super().__init__(dataset_dict)
        self.splits = True

    def _parameters(self):
        return {}

    def preprocess(self, ann_data, **kwargs):
        train_data = ann_data.copy()
        test_data = ann_data.copy()
        train_data.uns["split_role"] = "train"
        test_data.uns["split_role"] = "test_copy"
        return train_data, test_data
