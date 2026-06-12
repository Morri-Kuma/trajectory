"""Repair official-silver lineage outputs after stage-label contamination.

This script fixes the GSE178325/GSE230659 official-silver lineage evaluation
chain without retraining trajectory models:

1. Projection-capable methods:
   - backup old lineage/projection-derived artifacts;
   - rewrite projected_milestone_labels.csv so final_milestone_label_coarse
     contains only official coarse labels;
   - rebuild state_transition_matrix.csv and lineage_graph_edges.csv.

2. WOT / CellRank2:
   - recover their archived stage-label state_transition_matrix.csv files;
   - coarsen the matrix into the official reference label space.

3. Recompute lineage_metrics.json with eval_lineage.run_lineage_evaluation.

Transition labels such as stage_01_hADSCs_to_epithelial_like are mapped to
ambiguous for official metrics because the frozen coarse reference graph has no
transition-state nodes. This is intentionally conservative.
"""

from __future__ import annotations

import json
import re
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from benchmark.evaluation.build_projected_milestone_transitions import (
    run_transition_build,
)
from benchmark.evaluation.eval_lineage import run_lineage_evaluation

MARKERS_YAML = ROOT / "benchmark" / "annotation" / "milestone_markers.yaml"
REGISTRY_YAML = ROOT / "benchmark" / "ground_truth" / "registry.yaml"
INVALID_ARCHIVE = (
    ROOT
    / "benchmark"
    / "results"
    / "_invalid_lineage_label_space"
    / "20260605_stage_label_outputs"
    / "benchmark"
    / "results"
)
RUN_TAG = datetime.now().strftime("coarse_lineage_repair_%Y%m%d_%H%M%S")

PROJECTION_METHODS = ("mioflow", "prescient", "scnode")
LINEAGE_ONLY_METHODS = ("cellrank2", "wot")
DATASETS = ("gse178325", "gse230659")
SCENARIOS = ("A", "B", "C")

ARTIFACTS_TO_BACKUP = (
    "projected_milestone_labels.csv",
    "state_transition_matrix.csv",
    "lineage_graph_edges.csv",
    "lineage_metrics.json",
    "lineage_diagnostics.json",
    "projected_transition_metadata.json",
)


@dataclass(frozen=True)
class RunSpec:
    dataset: str
    method: str
    scenario: str
    result_dir: Path
    provider_id: str
    reference_graph_path: Path


def _load_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8-sig") as f:
        return yaml.safe_load(f) or {}


def _dataset_key(dataset: str) -> str:
    if dataset == "gse178325":
        return "GSE178325"
    if dataset == "gse230659":
        return "GSE230659"
    raise ValueError(f"Unsupported dataset: {dataset}")


def _run_id(dataset: str, method: str, scenario: str) -> str:
    return f"{dataset}_marker_fm_silver_{scenario}_hvg2000_formal"


def _provider_id(dataset: str) -> str:
    return f"{dataset}_marker_fm_transition_silver_v1"


def _reference_graph_path(dataset: str) -> Path:
    return (
        ROOT
        / "benchmark"
        / "ground_truth"
        / "providers"
        / _provider_id(dataset)
        / "reference_graph.json"
    )


def _result_dir(dataset: str, method: str, scenario: str) -> Path:
    return ROOT / "benchmark" / "results" / method / _run_id(dataset, method, scenario)


def _projection_specs() -> Iterable[RunSpec]:
    for dataset in DATASETS:
        for method in PROJECTION_METHODS:
            for scenario in SCENARIOS:
                yield RunSpec(
                    dataset=dataset,
                    method=method,
                    scenario=scenario,
                    result_dir=_result_dir(dataset, method, scenario),
                    provider_id=_provider_id(dataset),
                    reference_graph_path=_reference_graph_path(dataset),
                )


def _lineage_only_specs() -> Iterable[RunSpec]:
    dataset = "gse230659"
    scenario = "A"
    for method in LINEAGE_ONLY_METHODS:
        yield RunSpec(
            dataset=dataset,
            method=method,
            scenario=scenario,
            result_dir=_result_dir(dataset, method, scenario),
            provider_id=_provider_id(dataset),
            reference_graph_path=_reference_graph_path(dataset),
        )


def _backup_artifacts(result_dir: Path) -> Path:
    backup_dir = result_dir / "_repair_backup" / RUN_TAG
    backup_dir.mkdir(parents=True, exist_ok=True)
    for name in ARTIFACTS_TO_BACKUP:
        path = result_dir / name
        if path.exists():
            dest = backup_dir / name
            if not dest.exists():
                shutil.copy2(path, dest)
    return backup_dir


def _stage_to_coarse_map() -> dict[str, dict[str, str]]:
    markers = _load_yaml(MARKERS_YAML)
    mapping: dict[str, dict[str, str]] = {}
    for dataset_key, cfg in markers.items():
        stage_map: dict[str, str] = {}
        for stage in cfg.get("expanded_trajectory_order") or []:
            parts = str(stage).split("_", maxsplit=2)
            if len(parts) < 3:
                continue
            semantic = parts[2]
            stage_map[stage] = "ambiguous" if "_to_" in semantic else semantic
        stage_map["ambiguous"] = "ambiguous"
        stage_map["unknown_or_ood"] = "unknown_or_ood"
        mapping[dataset_key] = stage_map
    return mapping


def _coarsen_label(label: object, stage_map: dict[str, str]) -> str:
    value = str(label)
    if value in stage_map:
        return stage_map[value]
    if value.startswith("stage_"):
        parts = value.split("_", maxsplit=2)
        if len(parts) >= 3:
            semantic = parts[2]
            return "ambiguous" if "_to_" in semantic else semantic
    return value


def _stage_examples(values: Iterable[object]) -> list[str]:
    return sorted({str(v) for v in values if str(v).startswith("stage_")})[:20]


def _repair_projected_labels(spec: RunSpec, stage_map: dict[str, str]) -> None:
    csv_path = spec.result_dir / "projected_milestone_labels.csv"
    if not csv_path.exists():
        raise FileNotFoundError(csv_path)

    df = pd.read_csv(csv_path)
    if "final_milestone_label_coarse" not in df.columns:
        raise ValueError(f"{csv_path}: missing final_milestone_label_coarse")

    old = df["final_milestone_label_coarse"].astype(str)
    if "final_milestone_label_coarse_pre_repair" not in df.columns:
        df["final_milestone_label_coarse_pre_repair"] = old
    df["final_milestone_label_coarse"] = [
        _coarsen_label(label, stage_map) for label in old
    ]
    expanded = _stage_examples(df["final_milestone_label_coarse"])
    if expanded:
        raise ValueError(f"{csv_path}: stage labels remain after repair: {expanded}")
    df.to_csv(csv_path, index=False)


def _coarsen_stage_matrix(
    source_path: Path,
    output_dir: Path,
    stage_map: dict[str, str],
) -> None:
    if not source_path.exists():
        raise FileNotFoundError(source_path)
    raw = pd.read_csv(source_path, index_col=0)
    raw.index = raw.index.astype(str)
    raw.columns = raw.columns.astype(str)

    row_labels = [_coarsen_label(label, stage_map) for label in raw.index]
    col_labels = [_coarsen_label(label, stage_map) for label in raw.columns]
    values = raw.to_numpy(dtype=float)

    row_groups = pd.Index(row_labels, name="source_state")
    col_groups = pd.Index(col_labels, name="target_state")
    coarsened = pd.DataFrame(values, index=row_groups, columns=col_groups)

    drop_labels = {"ambiguous", "unknown_or_ood"}
    coarsened = coarsened.loc[
        ~coarsened.index.isin(drop_labels),
        ~coarsened.columns.isin(drop_labels),
    ]
    coarsened = coarsened.groupby(level=0).mean()
    coarsened = coarsened.T.groupby(level=0).sum().T
    row_sums = coarsened.sum(axis=1)
    nonzero = row_sums > 0
    coarsened.loc[nonzero] = coarsened.loc[nonzero].div(row_sums[nonzero], axis=0)

    expanded = _stage_examples(list(coarsened.index) + list(coarsened.columns))
    if expanded:
        raise ValueError(f"{source_path}: stage labels remain after coarsening: {expanded}")

    output_dir.mkdir(parents=True, exist_ok=True)
    stm_path = output_dir / "state_transition_matrix.csv"
    edges_path = output_dir / "lineage_graph_edges.csv"
    coarsened.to_csv(stm_path)

    rows = []
    for src in coarsened.index:
        for tgt in coarsened.columns:
            weight = float(coarsened.loc[src, tgt])
            if weight > 1e-6:
                rows.append(
                    {
                        "source_state": src,
                        "target_state": tgt,
                        "weight": round(weight, 6),
                    }
                )
    pd.DataFrame(rows, columns=["source_state", "target_state", "weight"]).to_csv(
        edges_path, index=False
    )


def _archived_stm_path(spec: RunSpec) -> Path:
    return INVALID_ARCHIVE / spec.method / spec.result_dir.name / "state_transition_matrix.csv"


def _ground_truth_dict(spec: RunSpec) -> dict:
    registry = (_load_yaml(REGISTRY_YAML).get("providers") or {})
    provider = dict(registry.get(spec.provider_id) or {})
    provider.setdefault("provider_id", spec.provider_id)
    provider.setdefault("state_key", "final_milestone_label_coarse")
    provider.setdefault("label_mode", "official_silver")
    provider.setdefault("reference_graph_path", str(spec.reference_graph_path.relative_to(ROOT)))
    return provider


def _run_lineage(spec: RunSpec) -> dict:
    metrics = run_lineage_evaluation(
        state_transition_matrix_path=str(spec.result_dir / "state_transition_matrix.csv"),
        lineage_graph_edges_path=str(spec.result_dir / "lineage_graph_edges.csv"),
        output_dir=str(spec.result_dir),
        reference_graph_path=str(spec.reference_graph_path),
        adata=None,
        edge_confidence_mode="all",
        exclude_uncertain_states=False,
        cell_state_key="final_milestone_label_coarse",
        time_key="abs_day",
        ground_truth=_ground_truth_dict(spec),
    )
    if metrics.get("status") != "completed":
        raise RuntimeError(f"{spec.result_dir}: lineage status={metrics.get('status')!r}")
    if metrics.get("prediction_state_key") != "final_milestone_label_coarse":
        raise RuntimeError(f"{spec.result_dir}: bad prediction_state_key")
    return metrics


def _clear_invalid_run_metadata_status(result_dir: Path) -> None:
    metadata_path = result_dir / "run_metadata.json"
    if not metadata_path.exists():
        return
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8-sig"))
    except Exception:
        return
    if metadata.get("status") == "failed_lineage_evaluation":
        metadata["status"] = "completed_after_coarse_lineage_repair"
    notes = str(metadata.get("notes") or "")
    notes = re.sub(
        r"\n?state_transition_matrix\.csv uses expanded/stage labels,.*?(?=\n|$)",
        "",
        notes,
    ).strip()
    if notes:
        metadata["notes"] = notes
    elif "notes" in metadata:
        metadata.pop("notes", None)
    metadata["lineage_repair"] = {
        "run_tag": RUN_TAG,
        "state_key": "final_milestone_label_coarse",
        "transition_labels_mapped_to": "ambiguous",
    }
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def _repair_projection_run(spec: RunSpec, stage_map: dict[str, str]) -> dict:
    _backup_artifacts(spec.result_dir)
    _repair_projected_labels(spec, stage_map)
    run_transition_build(
        projected_labels_csv=str(spec.result_dir / "projected_milestone_labels.csv"),
        output_dir=str(spec.result_dir),
        label_mode="official_silver",
        cell_state_key="final_milestone_label_coarse",
        exclude_labels=["ambiguous", "unknown_or_ood"],
        dataset_id=_dataset_key(spec.dataset),
        all_timepoint_pairs=False,
        provider_id=spec.provider_id,
    )
    metrics = _run_lineage(spec)
    _clear_invalid_run_metadata_status(spec.result_dir)
    return metrics


def _repair_lineage_only_run(spec: RunSpec, stage_map: dict[str, str]) -> dict:
    _backup_artifacts(spec.result_dir)
    _coarsen_stage_matrix(_archived_stm_path(spec), spec.result_dir, stage_map)
    metrics = _run_lineage(spec)
    _clear_invalid_run_metadata_status(spec.result_dir)
    return metrics


def main() -> int:
    stage_maps = _stage_to_coarse_map()
    repaired = []

    for spec in _projection_specs():
        dataset_key = _dataset_key(spec.dataset)
        metrics = _repair_projection_run(spec, stage_maps[dataset_key])
        repaired.append((spec, metrics))
        single = (metrics.get("graph_metrics") or {}).get("single_step") or {}
        print(
            f"[OK] {spec.method}/{spec.result_dir.name}: "
            f"AUROC={single.get('auc_roc')}, AUPRC={single.get('auc_prc')}"
        )

    for spec in _lineage_only_specs():
        dataset_key = _dataset_key(spec.dataset)
        metrics = _repair_lineage_only_run(spec, stage_maps[dataset_key])
        repaired.append((spec, metrics))
        single = (metrics.get("graph_metrics") or {}).get("single_step") or {}
        print(
            f"[OK] {spec.method}/{spec.result_dir.name}: "
            f"AUROC={single.get('auc_roc')}, AUPRC={single.get('auc_prc')}"
        )

    summary_path = (
        ROOT
        / "logs"
        / f"official_silver_lineage_repair_summary_{RUN_TAG}.json"
    )
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "run_tag": RUN_TAG,
        "transition_label_policy": "stage_*_A_to_B -> ambiguous",
        "n_repaired_runs": len(repaired),
        "runs": [
            {
                "dataset": spec.dataset,
                "method": spec.method,
                "scenario": spec.scenario,
                "result_dir": str(spec.result_dir.relative_to(ROOT)),
                "status": metrics.get("status"),
                "auroc": metrics.get("auroc"),
                "auprc": metrics.get("auprc"),
                "jaccard_similarity": metrics.get("jaccard_similarity"),
            }
            for spec, metrics in repaired
        ],
    }
    summary_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"[OK] repair summary -> {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
