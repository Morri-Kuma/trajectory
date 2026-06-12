"""Build official-silver embedding/lineage/forecast reports from validated run outputs.

The report is intentionally fail-closed: official embedding metrics may only
be combined with lineage metrics when both consume the state key declared by
the frozen ground-truth provider registry.

Lineage metrics require metric_protocol == "sctimebench_graph_sim".  Runs that
carry the old jaccard_similarity_topk schema are rejected at collection time so
they cannot silently dilute the official ranking.

Forecast metrics require metric_backend == "scTimeBench_exact" and
metric_protocol == "sctimebench_gex_prediction_otloss".  Runs without those
fields are collected as informational only; only exact-mode runs participate in
official ranking.

Ranking convention (lineage):
  single_step_auc_roc, single_step_auc_prc, single_step_jaccard,
  multi_step_auc_roc,  multi_step_auc_prc,  multi_step_jaccard
  — all ranked high-to-low.  Composite = mean within-dataset/scenario rank.

Ranking convention (forecast):
  wasserstein_distance, gaussian_mmd, energy_distance_mmd, hausdorff_loss
  — all ranked low-to-high (lower loss = better).
  Composite = mean within-dataset/scenario rank.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List

import pandas as pd
import yaml

from benchmark.evaluation.lineage_graphsim_sctimebench import (
    PREDICTION_LABEL_SOURCE_SCTIMEBENCH_ZERO_FILL,
)


ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = ROOT / "benchmark" / "results" / "result_manifest.yaml"
REGISTRY_PATH = ROOT / "benchmark" / "ground_truth" / "registry.yaml"
REPORT_DIR = ROOT / "benchmark" / "reports" / "official_silver"
REPORT_NAME = "official_silver_model_rankings.md"
EMBEDDING_METHODS = {"mioflow", "prescient", "scnode"}
FORECAST_METHODS = {"mioflow", "prescient", "scnode"}

# Metric protocol that lineage_metrics.json must declare.
REQUIRED_LINEAGE_PROTOCOL = "sctimebench_graph_sim"
REQUIRED_PREDICTION_LABEL_SOURCE = PREDICTION_LABEL_SOURCE_SCTIMEBENCH_ZERO_FILL

# Metadata that forecast_metrics.json must declare for official ranking.
REQUIRED_FORECAST_BACKEND = "scTimeBench_exact"
REQUIRED_FORECAST_PROTOCOL = "sctimebench_gex_prediction_otloss"

# Forecast loss metrics — all ranked low-to-high (lower = better).
_FORECAST_RANKING_METRICS = [
    "wasserstein_distance",
    "gaussian_mmd",
    "energy_distance_mmd",
    "hausdorff_loss",
]


def _load_yaml(path: Path) -> Dict[str, Any]:
    with open(path, encoding="utf-8-sig") as f:
        return yaml.safe_load(f) or {}


def _load_json(path: Path) -> Dict[str, Any]:
    with open(path, encoding="utf-8-sig") as f:
        return json.load(f)


def _dataset_id(manifest_key: str) -> str:
    return manifest_key.split("_", 1)[0].upper()


def _official_entries(manifest: Dict[str, Any]) -> Iterable[tuple[str, Dict[str, Any]]]:
    for key, entry in (manifest.get("official_silver") or {}).items():
        if Path(str(entry.get("summary", ""))).name == REPORT_NAME:
            yield key, entry


def _extract_lineage_row(
    lineage: Dict[str, Any],
    dataset_id: str,
    method: str,
    scenario: str,
    result_dir: Path,
    provider_id: str,
    expected_state_key: str,
    errors: List[str],
    relative_dir: str,
) -> Dict[str, Any]:
    """
    Pull all required lineage fields out of a lineage_metrics.json dict,
    validate the metric_protocol, and return a flat row dict.

    Returns {} on hard validation failure (caller should skip this run).
    """
    # Protocol guard — must be sctimebench_graph_sim
    protocol = lineage.get("metric_protocol")
    if protocol != REQUIRED_LINEAGE_PROTOCOL:
        errors.append(
            f"{relative_dir}: lineage metric_protocol={protocol!r}, "
            f"expected {REQUIRED_LINEAGE_PROTOCOL!r}. "
            "Re-run eval_lineage.py to generate scTimeBench graph-sim metrics."
        )
        return {}

    lineage_state_key = (
        lineage.get("ground_truth", {}).get("state_key")
        or lineage.get("cell_state_key")
    )
    if lineage.get("label_mode") != "official_silver":
        errors.append(f"{relative_dir}: lineage label_mode is not official_silver")
    if lineage.get("provider_id") != provider_id:
        errors.append(
            f"{relative_dir}: lineage provider_id={lineage.get('provider_id')!r}, "
            f"registry provider_id={provider_id!r}"
        )
    if lineage_state_key != expected_state_key:
        errors.append(
            f"{relative_dir}: lineage state_key={lineage_state_key!r}, "
            f"registry state_key={expected_state_key!r}"
        )
    if lineage.get("status") != "completed":
        errors.append(
            f"{relative_dir}: lineage status={lineage.get('status')!r}; "
            "completed scTimeBench-aligned lineage metrics are required for "
            "official ranking"
        )
        return {}
    if lineage.get("prediction_state_key") != expected_state_key:
        errors.append(
            f"{relative_dir}: prediction_state_key={lineage.get('prediction_state_key')!r}, "
            f"registry state_key={expected_state_key!r}"
        )
        return {}
    if lineage.get("prediction_label_source") != REQUIRED_PREDICTION_LABEL_SOURCE:
        errors.append(
            f"{relative_dir}: prediction_label_source="
            f"{lineage.get('prediction_label_source')!r}; rerun lineage "
            "evaluation with scTimeBench-style reference-node zero-fill alignment"
        )
        return {}

    gm = lineage.get("graph_metrics") or {}
    single = gm.get("single_step") or {}
    multi  = gm.get("multi_step")  or {}
    label_report = lineage.get("prediction_label_report") or {}
    missing_reference_nodes = label_report.get("missing_reference_nodes") or []

    return {
        "dataset_id":   dataset_id,
        "method":       method,
        "scenario":     scenario,
        "run_id":       result_dir.name,
        "state_key":    lineage_state_key,
        "provider_id":  lineage.get("provider_id"),
        "metric_protocol": protocol,
        # Single-step (simple criterion) metrics
        "single_step_auc_roc":  single.get("auc_roc"),
        "single_step_auc_prc":  single.get("auc_prc"),
        "single_step_jaccard":  single.get("jaccard_similarity"),
        "single_step_precision": single.get("precision"),
        "single_step_recall":   single.get("recall"),
        "single_step_f1":       single.get("f1"),
        # Multi-step (all_paths criterion) metrics
        "multi_step_auc_roc":   multi.get("auc_roc"),
        "multi_step_auc_prc":   multi.get("auc_prc"),
        "multi_step_jaccard":   multi.get("jaccard_similarity"),
        "multi_step_precision": multi.get("precision"),
        "multi_step_recall":    multi.get("recall"),
        "multi_step_f1":        multi.get("f1"),
        # Provenance
        "n_reference_edges": lineage.get("n_reference_edges"),
        "edge_confidence_mode": lineage.get("edge_confidence_mode"),
        "status": lineage.get("status"),
        "prediction_label_source": lineage.get("prediction_label_source"),
        "missing_reference_nodes": ";".join(map(str, missing_reference_nodes)),
        "n_missing_reference_nodes": label_report.get("n_missing_reference_nodes", 0),
        # Backward-compat aliases (for combined summary cross-checks)
        "auroc": lineage.get("auroc"),
        "auprc": lineage.get("auprc"),
    }


def _validate_and_collect(
    expected_n_neighbors: int,
    expected_resolution: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, List[str]]:
    manifest = _load_yaml(MANIFEST_PATH)
    registry = (_load_yaml(REGISTRY_PATH).get("providers") or {})
    embedding_rows: List[Dict[str, Any]] = []
    lineage_rows: List[Dict[str, Any]] = []
    forecast_rows: List[Dict[str, Any]] = []
    errors: List[str] = []
    expected_cluster_fragment = (
        f"n_neighbors={expected_n_neighbors}:resolution={float(expected_resolution)}"
    )

    for key, entry in _official_entries(manifest):
        dataset_id = _dataset_id(key)
        label_mode = entry.get("label_mode")
        providers  = entry.get("providers") or []
        if label_mode != "official_silver":
            errors.append(f"{dataset_id}: manifest label_mode={label_mode!r}")
        if len(providers) != 1 or providers[0] not in registry:
            errors.append(f"{dataset_id}: invalid manifest provider list {providers!r}")
            continue
        provider_id = providers[0]
        expected_state_key = registry[provider_id].get("state_key")

        for method, scenarios in (entry.get("result_dirs") or {}).items():
            for scenario, relative_dir in (scenarios or {}).items():
                result_dir   = ROOT / relative_dir
                lineage_path = result_dir / "lineage_metrics.json"
                embedding_path = (
                    result_dir
                    / "embedding_milestone_eval"
                    / "embedding_metrics_official_silver.json"
                )

                if not lineage_path.exists():
                    errors.append(f"{relative_dir}: missing lineage_metrics.json")
                if method in EMBEDDING_METHODS and not embedding_path.exists():
                    errors.append(
                        f"{relative_dir}: missing official embedding metrics for {method}"
                    )

                if lineage_path.exists():
                    lineage = _load_json(lineage_path)
                    row = _extract_lineage_row(
                        lineage=lineage,
                        dataset_id=dataset_id,
                        method=method,
                        scenario=scenario,
                        result_dir=result_dir,
                        provider_id=provider_id,
                        expected_state_key=expected_state_key,
                        errors=errors,
                        relative_dir=str(relative_dir),
                    )
                    if row:
                        lineage_rows.append(row)

                if embedding_path.exists():
                    embedding = _load_json(embedding_path)
                    if embedding.get("label_mode") != "official_silver":
                        errors.append(f"{relative_dir}: embedding label_mode is not official_silver")
                    if embedding.get("provider_id") != provider_id:
                        errors.append(
                            f"{relative_dir}: embedding provider_id={embedding.get('provider_id')!r}, "
                            f"registry provider_id={provider_id!r}"
                        )
                    if embedding.get("state_key") != expected_state_key:
                        errors.append(
                            f"{relative_dir}: embedding state_key={embedding.get('state_key')!r}, "
                            f"registry state_key={expected_state_key!r}"
                        )
                    if expected_cluster_fragment not in str(embedding.get("cluster_source", "")):
                        errors.append(
                            f"{relative_dir}: embedding cluster_source does not contain "
                            f"{expected_cluster_fragment!r}"
                        )
                    if embedding.get("embedding_metric_protocol") != "sctimebench":
                        errors.append(
                            f"{relative_dir}: embedding_metric_protocol is not sctimebench"
                        )
                    if embedding.get("entropy_basis") != "classifier_probability_vector":
                        errors.append(
                            f"{relative_dir}: entropy_basis is not classifier_probability_vector"
                        )
                    if embedding.get("label_source") != "knn_transfer_from_observed_embedding":
                        errors.append(
                            f"{relative_dir}: projected labels are not from scTimeBench kNN transfer"
                        )
                    embedding_rows.append({
                        "dataset_id": dataset_id,
                        "method": method,
                        "scenario": scenario,
                        "run_id": result_dir.name,
                        "embedding_metric_protocol": embedding.get(
                            "embedding_metric_protocol"
                        ),
                        "state_key": embedding.get("state_key"),
                        "provider_id": embedding.get("provider_id"),
                        "adjusted_rand_index": embedding.get("adjusted_rand_index"),
                        "ari_next_timepoint": embedding.get("ari_next_timepoint"),
                        "reference_adjusted_rand_index": embedding.get(
                            "reference_adjusted_rand_index"
                        ),
                        "ari_ground_truth": embedding.get("ari_ground_truth"),
                        "ari_retention_fraction": embedding.get("ari_retention_fraction"),
                        "pred_tp_avg_normalized_entropy": embedding.get(
                            "pred_tp_avg_normalized_entropy"
                        ),
                        "avg_normalized_entropy": embedding.get("avg_normalized_entropy"),
                        "n_cells_evaluated": embedding.get("n_cells_evaluated"),
                        "cluster_source": embedding.get("cluster_source"),
                        "status": embedding.get("status"),
                    })

                if lineage_path.exists() and embedding_path.exists() and lineage_rows:
                    last_lineage_key = lineage_rows[-1].get("state_key")
                    if embedding.get("state_key") != last_lineage_key:
                        errors.append(
                            f"{relative_dir}: embedding state_key={embedding.get('state_key')!r} "
                            f"does not match lineage state_key={last_lineage_key!r}"
                        )

                # Forecast metrics collection (informational for all methods;
                # only exact-mode runs participate in official ranking).
                if method in FORECAST_METHODS:
                    forecast_path = result_dir / "forecast_metrics.json"
                    if forecast_path.exists():
                        fm = _load_json(forecast_path)
                        exact_mode = (
                            fm.get("metric_backend") == REQUIRED_FORECAST_BACKEND
                            and fm.get("metric_protocol") == REQUIRED_FORECAST_PROTOCOL
                            and fm.get("exact_cell_usage") is True
                            and fm.get("lognorm") is False
                            and fm.get("status") == "completed"
                        )
                        forecast_rows.append({
                            "dataset_id": dataset_id,
                            "method": method,
                            "scenario": scenario,
                            "run_id": result_dir.name,
                            "metric_backend": fm.get("metric_backend"),
                            "metric_protocol": fm.get("metric_protocol"),
                            "exact_cell_usage": fm.get("exact_cell_usage"),
                            "lognorm": fm.get("lognorm"),
                            "aggregate": fm.get("aggregate"),
                            "n_eval_timepoints": fm.get("n_eval_timepoints"),
                            "wasserstein_distance": fm.get("wasserstein_distance"),
                            "gaussian_mmd": fm.get("gaussian_mmd"),
                            "energy_distance_mmd": fm.get("energy_distance_mmd"),
                            "hausdorff_loss": fm.get("hausdorff_loss"),
                            "status": fm.get("status"),
                            "official_ranking_eligible": exact_mode,
                        })

    return (
        pd.DataFrame(embedding_rows),
        pd.DataFrame(lineage_rows),
        pd.DataFrame(forecast_rows),
        errors,
    )


def _rank_embedding(df: pd.DataFrame) -> pd.DataFrame:
    ranked = df[df["status"] == "completed"].copy()
    groups = ranked.groupby(["dataset_id", "scenario"])
    ranked["rank_ari"] = groups["adjusted_rand_index"].rank(ascending=False, method="min")
    ranked["rank_classifier_entropy"] = groups["pred_tp_avg_normalized_entropy"].rank(
        ascending=True, method="min"
    )
    ranked["embedding_composite_rank_score"] = ranked[
        ["rank_ari", "rank_classifier_entropy"]
    ].mean(axis=1)
    ranked["embedding_composite_rank"] = ranked.groupby(
        ["dataset_id", "scenario"]
    )["embedding_composite_rank_score"].rank(ascending=True, method="min")
    return ranked.sort_values(
        ["dataset_id", "scenario", "embedding_composite_rank", "method"]
    )


# scTimeBench graph-sim metrics used for official lineage ranking
_LINEAGE_RANKING_METRICS = [
    "single_step_auc_roc",
    "single_step_auc_prc",
    "single_step_jaccard",
    "multi_step_auc_roc",
    "multi_step_auc_prc",
    "multi_step_jaccard",
]


def _rank_lineage(df: pd.DataFrame) -> pd.DataFrame:
    ranked = df[df["status"] == "completed"].copy()
    groups = ranked.groupby(["dataset_id", "scenario"])
    for metric in _LINEAGE_RANKING_METRICS:
        ranked[f"rank_{metric}"] = groups[metric].rank(ascending=False, method="min")
    ranked["lineage_composite_rank_score"] = ranked[
        [f"rank_{m}" for m in _LINEAGE_RANKING_METRICS]
    ].mean(axis=1)
    ranked["lineage_composite_rank"] = ranked.groupby(
        ["dataset_id", "scenario"]
    )["lineage_composite_rank_score"].rank(ascending=True, method="min")
    return ranked.sort_values(
        ["dataset_id", "scenario", "lineage_composite_rank", "method"]
    )


def _rank_forecast(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rank forecast loss metrics within each dataset + scenario.

    Only runs with official_ranking_eligible=True participate.
    All four metrics are ranked low-to-high (lower loss = better).
    Composite rank score = mean of the four individual ranks.
    Composite rank = rank of composite rank score low-to-high.
    """
    if df.empty or "official_ranking_eligible" not in df.columns:
        return pd.DataFrame()
    eligible = df[df["official_ranking_eligible"] == True].copy()  # noqa: E712
    if eligible.empty:
        return eligible
    groups = eligible.groupby(["dataset_id", "scenario"])
    for metric in _FORECAST_RANKING_METRICS:
        eligible[f"rank_{metric}"] = groups[metric].rank(ascending=True, method="min")
    eligible["forecast_composite_rank_score"] = eligible[
        [f"rank_{m}" for m in _FORECAST_RANKING_METRICS]
    ].mean(axis=1)
    eligible["forecast_composite_rank"] = eligible.groupby(
        ["dataset_id", "scenario"]
    )["forecast_composite_rank_score"].rank(ascending=True, method="min")
    return eligible.sort_values(
        ["dataset_id", "scenario", "forecast_composite_rank", "method"]
    )


def _forecast_summary(ranked: pd.DataFrame) -> pd.DataFrame:
    if ranked.empty:
        return ranked
    summary = ranked.groupby("method", as_index=False).agg(
        forecast_runs=("run_id", "count"),
        mean_forecast_composite_rank_score=("forecast_composite_rank_score", "mean"),
        mean_wasserstein_distance=("wasserstein_distance", "mean"),
        mean_gaussian_mmd=("gaussian_mmd", "mean"),
        mean_energy_distance_mmd=("energy_distance_mmd", "mean"),
        mean_hausdorff_loss=("hausdorff_loss", "mean"),
    )
    summary["overall_forecast_rank"] = summary[
        "mean_forecast_composite_rank_score"
    ].rank(ascending=True, method="min")
    return summary.sort_values(["overall_forecast_rank", "method"])


def _embedding_summary(ranked: pd.DataFrame) -> pd.DataFrame:
    summary = ranked.groupby("method", as_index=False).agg(
        embedding_runs=("run_id", "count"),
        mean_embedding_composite_rank_score=("embedding_composite_rank_score", "mean"),
        mean_ari=("adjusted_rand_index", "mean"),
        median_ari=("adjusted_rand_index", "median"),
        mean_classifier_entropy=("pred_tp_avg_normalized_entropy", "mean"),
    )
    summary["overall_embedding_rank"] = summary[
        "mean_embedding_composite_rank_score"
    ].rank(ascending=True, method="min")
    return summary.sort_values(["overall_embedding_rank", "method"])


def _lineage_summary(ranked: pd.DataFrame) -> pd.DataFrame:
    summary = ranked.groupby("method", as_index=False).agg(
        lineage_runs=("run_id", "count"),
        mean_lineage_composite_rank_score=("lineage_composite_rank_score", "mean"),
        mean_single_step_auc_roc=("single_step_auc_roc", "mean"),
        mean_single_step_auc_prc=("single_step_auc_prc", "mean"),
        mean_single_step_jaccard=("single_step_jaccard", "mean"),
        mean_multi_step_auc_roc=("multi_step_auc_roc", "mean"),
        mean_multi_step_auc_prc=("multi_step_auc_prc", "mean"),
        mean_multi_step_jaccard=("multi_step_jaccard", "mean"),
        # Optional precision/recall/f1 summaries
        mean_single_step_precision=("single_step_precision", "mean"),
        mean_single_step_recall=("single_step_recall", "mean"),
        mean_single_step_f1=("single_step_f1", "mean"),
        mean_multi_step_precision=("multi_step_precision", "mean"),
        mean_multi_step_recall=("multi_step_recall", "mean"),
        mean_multi_step_f1=("multi_step_f1", "mean"),
    )
    summary["overall_lineage_rank"] = summary[
        "mean_lineage_composite_rank_score"
    ].rank(ascending=True, method="min")
    return summary.sort_values(["overall_lineage_rank", "method"])


def _combined_summary(embedding: pd.DataFrame, lineage: pd.DataFrame) -> pd.DataFrame:
    combined = embedding.merge(lineage, on="method", how="inner")
    combined = combined[
        (combined["embedding_runs"] == combined["lineage_runs"])
        & (combined["embedding_runs"] > 1)
    ].copy()
    combined["combined_rank_score"] = combined[
        ["mean_embedding_composite_rank_score", "mean_lineage_composite_rank_score"]
    ].mean(axis=1)
    combined["overall_combined_rank"] = combined["combined_rank_score"].rank(
        ascending=True, method="min"
    )
    return combined.sort_values(["overall_combined_rank", "method"])


def _markdown_table(df: pd.DataFrame, columns: List[str], decimals: int = 4) -> str:
    view = df[columns].copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda value: f"{value:.{decimals}f}")
    header  = "| " + " | ".join(columns) + " |"
    divider = "| " + " | ".join(["---"] * len(columns)) + " |"
    rows = [
        "| " + " | ".join(str(value) for value in row) + " |"
        for row in view.itertuples(index=False, name=None)
    ]
    return "\n".join([header, divider, *rows])


def _write_report(
    embedding_rankings: pd.DataFrame,
    lineage_rankings: pd.DataFrame,
    embedding_summary: pd.DataFrame,
    lineage_summary: pd.DataFrame,
    combined_summary: pd.DataFrame,
    forecast_rankings: pd.DataFrame,
    forecast_summary: pd.DataFrame,
    expected_n_neighbors: int,
    expected_resolution: float,
) -> None:
    REPORT_DIR.mkdir(parents=True, exist_ok=True)
    embedding_rankings.to_csv(REPORT_DIR / "official_silver_embedding_rankings.csv", index=False)
    lineage_rankings.to_csv(REPORT_DIR / "official_silver_lineage_rankings.csv", index=False)
    embedding_summary.to_csv(REPORT_DIR / "official_silver_embedding_method_summary.csv", index=False)
    lineage_summary.to_csv(REPORT_DIR / "official_silver_lineage_method_summary.csv", index=False)
    combined_summary.to_csv(
        REPORT_DIR / "official_silver_combined_method_summary_scnode_mioflow_prescient.csv",
        index=False,
    )
    if not forecast_rankings.empty:
        forecast_rankings.to_csv(
            REPORT_DIR / "official_silver_forecast_rankings.csv", index=False
        )
    if not forecast_summary.empty:
        forecast_summary.to_csv(
            REPORT_DIR / "official_silver_forecast_method_summary.csv", index=False
        )

    emb_winners = embedding_rankings[
        embedding_rankings["embedding_composite_rank"] == 1
    ]
    lin_winners = lineage_rankings[
        lineage_rankings["lineage_composite_rank"] == 1
    ]
    state_keys = ", ".join(sorted(set(embedding_rankings["state_key"].dropna())))
    cluster_source = (
        "leiden_from_next_timepoint_embedding:next_timepoint_embedding.npy:"
        f"n_neighbors={expected_n_neighbors}:resolution={float(expected_resolution)}"
    )

    # Build optional forecast section
    if not forecast_rankings.empty and not forecast_summary.empty:
        fc_winners = forecast_rankings[
            forecast_rankings["forecast_composite_rank"] == 1
        ]
        forecast_section = f"""
## Forecast Accuracy Method Summary (scTimeBench OT-Loss, exact mode)

Metrics: wasserstein_distance, gaussian_mmd, energy_distance_mmd, hausdorff_loss —
all ranked low-to-high (lower loss = better).  Only runs with
`metric_backend=scTimeBench_exact`, `metric_protocol=sctimebench_gex_prediction_otloss`,
`exact_cell_usage=True`, and `lognorm=False` are included.

{_markdown_table(forecast_summary, [
    "overall_forecast_rank", "method", "forecast_runs",
    "mean_forecast_composite_rank_score",
    "mean_wasserstein_distance", "mean_gaussian_mmd",
    "mean_energy_distance_mmd", "mean_hausdorff_loss",
])}

## Forecast Winners By Dataset/Scenario

{_markdown_table(fc_winners, [
    "dataset_id", "scenario", "method",
    "wasserstein_distance", "gaussian_mmd",
    "energy_distance_mmd", "hausdorff_loss",
    "forecast_composite_rank_score",
])}
"""
    else:
        forecast_section = "\n## Forecast Accuracy\n\nNo official-mode forecast runs collected.\n"

    markdown = f"""# Official Silver Model Rankings

Generated by `python -m benchmark.evaluation.summarize_official_silver` from
validated `official_silver` run outputs. The report is written only when
embedding, lineage, and the frozen provider registry use the same state key
and lineage runs declare `metric_protocol: {REQUIRED_LINEAGE_PROTOCOL}`.

Official contract: `state_key={state_keys}`; embedding clustering uses
`{cluster_source}`.

Ranking convention: ARI, AUROC, AUPRC, Jaccard metrics are ranked high-to-low;
scTimeBench projected classifier entropy is ranked low-to-high.
Forecast loss metrics are ranked low-to-high (lower = better).
Composite rank scores are averages of within-dataset/scenario ranks, so lower is
better.

Lineage ranking uses six scTimeBench graph-sim metrics:
single_step_auc_roc, single_step_auc_prc, single_step_jaccard,
multi_step_auc_roc, multi_step_auc_prc, multi_step_jaccard.

Forecast ranking uses four scTimeBench OT-loss metrics:
wasserstein_distance, gaussian_mmd, energy_distance_mmd, hausdorff_loss.

## Combined Summary: scNODE / MIOFlow / PRESCIENT

{_markdown_table(combined_summary, [
    "overall_combined_rank", "method", "embedding_runs", "lineage_runs",
    "combined_rank_score", "mean_ari", "mean_classifier_entropy",
    "mean_single_step_auc_prc", "mean_single_step_auc_roc",
    "mean_single_step_jaccard",
])}

## Embedding Coherence Method Summary

{_markdown_table(embedding_summary, [
    "overall_embedding_rank", "method", "embedding_runs",
    "mean_embedding_composite_rank_score", "mean_ari", "median_ari",
    "mean_classifier_entropy",
])}

## Lineage Fidelity Method Summary (scTimeBench graph-sim)

{_markdown_table(lineage_summary, [
    "overall_lineage_rank", "method", "lineage_runs",
    "mean_lineage_composite_rank_score",
    "mean_single_step_auc_roc", "mean_single_step_auc_prc", "mean_single_step_jaccard",
    "mean_multi_step_auc_roc",  "mean_multi_step_auc_prc",  "mean_multi_step_jaccard",
])}

## Embedding Winners By Dataset/Scenario

{_markdown_table(emb_winners, [
    "dataset_id", "scenario", "method", "adjusted_rand_index",
    "pred_tp_avg_normalized_entropy", "embedding_composite_rank_score",
])}

## Lineage Winners By Dataset/Scenario

{_markdown_table(lin_winners, [
    "dataset_id", "scenario", "method",
    "single_step_auc_roc", "single_step_auc_prc", "single_step_jaccard",
    "multi_step_auc_roc",  "multi_step_auc_prc",  "multi_step_jaccard",
    "lineage_composite_rank_score",
])}
{forecast_section}
## Notes

- WOT and CellRank2 are summarized in lineage tables only when they have no official embedding output.
- Lineage runs with `metric_protocol` other than `{REQUIRED_LINEAGE_PROTOCOL}` are excluded from
  this report. Re-run `eval_lineage.py` on those outputs to upgrade them.
- `legacy_jaccard_similarity_topk` (raw edge Jaccard from lineage_graph_edges.csv) is stored
  in lineage_metrics.json for diagnostic traceability only; it is not used for ranking.
- `expanded_sensitivity` embedding results are excluded from this formal report until the frozen
  provider registry and lineage benchmark are explicitly migrated and rerun under that state system.
- Forecast runs without `metric_backend=scTimeBench_exact` are informational only and excluded from
  official forecast ranking. Re-run `eval_forecast.py` (max_cells_per_timepoint=None) to upgrade.
- The complete per-run ranking tables are saved as CSV files in this same directory.
"""
    (REPORT_DIR / REPORT_NAME).write_text(markdown, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--leiden-n-neighbors", type=int, default=15)
    parser.add_argument("--leiden-resolution", type=float, default=0.5)
    args = parser.parse_args()

    embedding, lineage, forecast, errors = _validate_and_collect(
        expected_n_neighbors=args.leiden_n_neighbors,
        expected_resolution=args.leiden_resolution,
    )
    if errors:
        print("Official-silver report contract validation failed:")
        for error in errors:
            print(f"- {error}")
        print("No report files were written.")
        return 1
    if embedding.empty or lineage.empty:
        print("No report files were written: embedding or lineage result set is empty.")
        return 1

    embedding_rankings = _rank_embedding(embedding)
    lineage_rankings   = _rank_lineage(lineage)
    embedding_method_summary = _embedding_summary(embedding_rankings)
    lineage_method_summary   = _lineage_summary(lineage_rankings)
    combined_method_summary  = _combined_summary(
        embedding_method_summary, lineage_method_summary
    )
    forecast_rankings  = _rank_forecast(forecast)
    forecast_method_summary = _forecast_summary(forecast_rankings)

    _write_report(
        embedding_rankings,
        lineage_rankings,
        embedding_method_summary,
        lineage_method_summary,
        combined_method_summary,
        forecast_rankings,
        forecast_method_summary,
        expected_n_neighbors=args.leiden_n_neighbors,
        expected_resolution=args.leiden_resolution,
    )
    fc_msg = (
        f" and {len(forecast_rankings)} forecast runs"
        if not forecast_rankings.empty
        else ""
    )
    print(
        f"Wrote official-silver report tables for {len(embedding_rankings)} embedding "
        f"runs, {len(lineage_rankings)} lineage runs{fc_msg} to {REPORT_DIR}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
