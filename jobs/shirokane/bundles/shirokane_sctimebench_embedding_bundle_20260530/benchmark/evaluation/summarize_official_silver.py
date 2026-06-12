"""Build official-silver embedding/lineage reports from validated run outputs.

The report is intentionally fail-closed: official embedding metrics may only
be combined with lineage metrics when both consume the state key declared by
the frozen ground-truth provider registry.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List

import pandas as pd
import yaml


ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = ROOT / "benchmark" / "results" / "result_manifest.yaml"
REGISTRY_PATH = ROOT / "benchmark" / "ground_truth" / "registry.yaml"
REPORT_DIR = ROOT / "benchmark" / "reports" / "official_silver"
REPORT_NAME = "official_silver_model_rankings.md"
EMBEDDING_METHODS = {"mioflow", "prescient", "scnode"}


def _load_yaml(path: Path) -> Dict[str, Any]:
    with open(path, encoding="utf-8-sig") as f:
        return yaml.safe_load(f) or {}


def _load_json(path: Path) -> Dict[str, Any]:
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def _dataset_id(manifest_key: str) -> str:
    return manifest_key.split("_", 1)[0].upper()


def _official_entries(manifest: Dict[str, Any]) -> Iterable[tuple[str, Dict[str, Any]]]:
    for key, entry in (manifest.get("official_silver") or {}).items():
        if Path(str(entry.get("summary", ""))).name == REPORT_NAME:
            yield key, entry


def _validate_and_collect(
    expected_n_neighbors: int,
    expected_resolution: float,
) -> tuple[pd.DataFrame, pd.DataFrame, List[str]]:
    manifest = _load_yaml(MANIFEST_PATH)
    registry = (_load_yaml(REGISTRY_PATH).get("providers") or {})
    embedding_rows: List[Dict[str, Any]] = []
    lineage_rows: List[Dict[str, Any]] = []
    errors: List[str] = []
    expected_cluster_fragment = (
        f"n_neighbors={expected_n_neighbors}:resolution={float(expected_resolution)}"
    )

    for key, entry in _official_entries(manifest):
        dataset_id = _dataset_id(key)
        label_mode = entry.get("label_mode")
        providers = entry.get("providers") or []
        if label_mode != "official_silver":
            errors.append(f"{dataset_id}: manifest label_mode={label_mode!r}")
        if len(providers) != 1 or providers[0] not in registry:
            errors.append(f"{dataset_id}: invalid manifest provider list {providers!r}")
            continue
        provider_id = providers[0]
        expected_state_key = registry[provider_id].get("state_key")

        for method, scenarios in (entry.get("result_dirs") or {}).items():
            for scenario, relative_dir in (scenarios or {}).items():
                result_dir = ROOT / relative_dir
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
                    lineage_rows.append({
                        "dataset_id": dataset_id,
                        "method": method,
                        "scenario": scenario,
                        "run_id": result_dir.name,
                        "state_key": lineage_state_key,
                        "provider_id": lineage.get("provider_id"),
                        "auroc": lineage.get("auroc"),
                        "auprc": lineage.get("auprc"),
                        "jaccard_topk": lineage.get("jaccard_similarity_topk"),
                        "single_step_recovery": lineage.get("single_step_recovery"),
                        "multi_step_recovery": lineage.get("multi_step_recovery"),
                        "n_reference_edges": lineage.get("n_reference_edges"),
                        "status": lineage.get("status"),
                    })

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

                if lineage_path.exists() and embedding_path.exists():
                    lineage_key = lineage_rows[-1]["state_key"]
                    if embedding.get("state_key") != lineage_key:
                        errors.append(
                            f"{relative_dir}: embedding state_key={embedding.get('state_key')!r} "
                            f"does not match lineage state_key={lineage_key!r}"
                        )

    return pd.DataFrame(embedding_rows), pd.DataFrame(lineage_rows), errors


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


def _rank_lineage(df: pd.DataFrame) -> pd.DataFrame:
    ranked = df[df["status"] == "completed"].copy()
    metrics = [
        "auroc",
        "auprc",
        "jaccard_topk",
        "single_step_recovery",
        "multi_step_recovery",
    ]
    groups = ranked.groupby(["dataset_id", "scenario"])
    for metric in metrics:
        ranked[f"rank_{metric}"] = groups[metric].rank(ascending=False, method="min")
    ranked["lineage_composite_rank_score"] = ranked[
        [f"rank_{metric}" for metric in metrics]
    ].mean(axis=1)
    ranked["lineage_composite_rank"] = ranked.groupby(
        ["dataset_id", "scenario"]
    )["lineage_composite_rank_score"].rank(ascending=True, method="min")
    return ranked.sort_values(
        ["dataset_id", "scenario", "lineage_composite_rank", "method"]
    )


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
        mean_auprc=("auprc", "mean"),
        mean_auroc=("auroc", "mean"),
        mean_jaccard_topk=("jaccard_topk", "mean"),
        mean_single_step_recovery=("single_step_recovery", "mean"),
        mean_multi_step_recovery=("multi_step_recovery", "mean"),
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
    header = "| " + " | ".join(columns) + " |"
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
    markdown = f"""# Official Silver Model Rankings

Generated by `python -m benchmark.evaluation.summarize_official_silver` from
validated `official_silver` run outputs. The report is written only when
embedding, lineage, and the frozen provider registry use the same state key.

Official contract: `state_key={state_keys}`; embedding clustering uses
`{cluster_source}`.

Ranking convention: ARI, AUROC, AUPRC, Jaccard, and recovery metrics are ranked
high-to-low; scTimeBench projected classifier entropy is ranked low-to-high.
Composite rank scores are averages of within-dataset/scenario ranks, so lower is
better.

## Combined Summary: scNODE / MIOFlow / PRESCIENT

{_markdown_table(combined_summary, [
    "overall_combined_rank", "method", "embedding_runs", "lineage_runs",
    "combined_rank_score", "mean_ari", "mean_classifier_entropy", "mean_auprc",
    "mean_auroc", "mean_jaccard_topk",
])}

## Embedding Coherence Method Summary

{_markdown_table(embedding_summary, [
    "overall_embedding_rank", "method", "embedding_runs",
    "mean_embedding_composite_rank_score", "mean_ari", "median_ari",
    "mean_classifier_entropy",
])}

## Lineage Fidelity Method Summary

{_markdown_table(lineage_summary, [
    "overall_lineage_rank", "method", "lineage_runs",
    "mean_lineage_composite_rank_score", "mean_auprc", "mean_auroc",
    "mean_jaccard_topk", "mean_single_step_recovery",
    "mean_multi_step_recovery",
])}

## Embedding Winners By Dataset/Scenario

{_markdown_table(emb_winners, [
    "dataset_id", "scenario", "method", "adjusted_rand_index",
    "pred_tp_avg_normalized_entropy", "embedding_composite_rank_score",
])}

## Lineage Winners By Dataset/Scenario

{_markdown_table(lin_winners, [
    "dataset_id", "scenario", "method", "auprc", "auroc", "jaccard_topk",
    "single_step_recovery", "multi_step_recovery",
    "lineage_composite_rank_score",
])}

## Notes

- WOT and CellRank2 are summarized in lineage tables only when they have no official embedding output.
- `expanded_sensitivity` embedding results are excluded from this formal report until the frozen provider registry and lineage benchmark are explicitly migrated and rerun under that state system.
- The complete per-run ranking tables are saved as CSV files in this same directory.
"""
    (REPORT_DIR / REPORT_NAME).write_text(markdown, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--leiden-n-neighbors", type=int, default=15)
    parser.add_argument("--leiden-resolution", type=float, default=0.5)
    args = parser.parse_args()

    embedding, lineage, errors = _validate_and_collect(
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
    lineage_rankings = _rank_lineage(lineage)
    embedding_method_summary = _embedding_summary(embedding_rankings)
    lineage_method_summary = _lineage_summary(lineage_rankings)
    combined_method_summary = _combined_summary(
        embedding_method_summary, lineage_method_summary
    )
    _write_report(
        embedding_rankings,
        lineage_rankings,
        embedding_method_summary,
        lineage_method_summary,
        combined_method_summary,
        expected_n_neighbors=args.leiden_n_neighbors,
        expected_resolution=args.leiden_resolution,
    )
    print(
        f"Wrote official-silver report tables for {len(embedding_rankings)} embedding "
        f"runs and {len(lineage_rankings)} lineage runs to {REPORT_DIR}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
