"""Factory helpers for trajectory's scTimeBench-style dataset layer."""

from __future__ import annotations

from pathlib import Path

import yaml

import benchmark.shared.dataset  # noqa: F401 - triggers registry imports
from benchmark.shared.dataset.base import (
    DATASET_PREPROCESSOR_REGISTRY,
    DATASET_REGISTRY,
)


def _dataset_class_name(dataset_dict):
    name = dataset_dict.get("name")
    if name:
        return name
    dataset_id = dataset_dict.get("id")
    if dataset_id == "GSE230659":
        return "GSE230659Dataset"
    if dataset_id == "GSE178325":
        return "GSE178325Dataset"
    if dataset_id == "GSE242424":
        return "GSE242424Dataset"
    if dataset_id == "GSE298212":
        return "GSE298212Dataset"
    if dataset_id == "GSE218855":
        return "GSE218855Dataset"
    raise ValueError(f"Cannot infer dataset registry class from: {dataset_dict}")


def build_preprocessors(dataset_dict):
    steps = list(dataset_dict.get("data_preprocessing_steps") or [])
    if not steps:
        scenario_params = dataset_dict.get("scenario_params") or {}
        if scenario_params:
            steps.append(
                {
                    "name": "ScenarioTimepointSplit",
                    "time_key": dataset_dict.get("time_key", "abs_day"),
                    "train_times": scenario_params.get("train_times"),
                    "heldout_times": scenario_params.get("heldout_times"),
                }
            )
        else:
            steps.append({"name": "CopyTrainTest"})

    preprocessors = []
    for step in steps:
        step_name = step["name"]
        if step_name not in DATASET_PREPROCESSOR_REGISTRY:
            raise ValueError(
                f"Unknown dataset preprocessor {step_name!r}. "
                f"Available: {sorted(DATASET_PREPROCESSOR_REGISTRY)}"
            )
        cls = DATASET_PREPROCESSOR_REGISTRY[step_name]
        kwargs = {k: v for k, v in step.items() if k != "name"}
        preprocessors.append(cls(dataset_dict, **kwargs))
    return preprocessors


def build_dataset(dataset_dict, output_dir="benchmark/results"):
    dataset_dict = dict(dataset_dict)
    class_name = _dataset_class_name(dataset_dict)
    if class_name not in DATASET_REGISTRY:
        raise ValueError(
            f"Unknown dataset registry class {class_name!r}. "
            f"Available: {sorted(DATASET_REGISTRY)}"
        )
    preprocessors = build_preprocessors(dataset_dict)
    return DATASET_REGISTRY[class_name](dataset_dict, preprocessors, output_dir)


def dataset_config_from_method_config(method_config, project_root=None):
    dataset_cfg = dict(method_config.get("dataset") or {})
    if not dataset_cfg:
        raise ValueError("Method config does not contain a dataset section.")
    dataset_cfg.setdefault("name", f"{dataset_cfg.get('id')}Dataset")
    if project_root is not None:
        dataset_cfg.setdefault("project_root", str(project_root))
    if method_config.get("scenario_params"):
        dataset_cfg.setdefault("scenario_params", method_config["scenario_params"])
    lineage = method_config.get("lineage") or {}
    if lineage.get("cell_state_key"):
        dataset_cfg.setdefault("cell_state_key", lineage["cell_state_key"])
    return dataset_cfg


def build_dataset_from_method_config(method_config, output_dir="benchmark/results", project_root=None):
    return build_dataset(
        dataset_config_from_method_config(method_config, project_root=project_root),
        output_dir=output_dir,
    )


def load_yaml(path):
    with open(path, encoding="utf-8-sig") as f:
        return yaml.safe_load(f)


def build_dataset_from_method_config_path(path, output_dir="benchmark/results", project_root=None):
    path = Path(path)
    config = load_yaml(path)
    if project_root is None:
        project_root = path.resolve().parents[2]
    return build_dataset_from_method_config(
        config,
        output_dir=output_dir,
        project_root=project_root,
    )
