"""Build a meeting-ready trajectory result summary figure.

The figure merges the current benchmark reports for GSE178325, GSE230659,
and GSE242424, then summarizes all three benchmark dimensions:

1. Lineage Fidelity
2. Embedding Coherence
3. Forecast Accuracy
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib import pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
REPORT_ROOT = ROOT / "benchmark" / "reports"
PRIMARY_REPORT_DIR = REPORT_ROOT / "official_silver"
GSE242424_REPORT_DIR = REPORT_ROOT / "gse242424_oskm_ground_truth"
OUT_DIR = ROOT / "docs" / "meetings"
PNG_OUT = OUT_DIR / "trajectory_results_visual_summary.png"
PDF_OUT = OUT_DIR / "trajectory_results_visual_summary.pdf"
MD_OUT = OUT_DIR / "trajectory_results_visual_summary.md"
TABLE_OUT = OUT_DIR / "trajectory_results_visual_summary_tables.csv"

DATASET_SCENARIO_ORDER = [
    "GSE178325-A",
    "GSE178325-B",
    "GSE178325-C",
    "GSE230659-A",
    "GSE230659-B",
    "GSE230659-C",
    "GSE242424-A",
    "GSE242424-B",
    "GSE242424-C",
]
PROJECTION_METHOD_ORDER = ["mioflow", "prescient", "scnode"]
LINEAGE_METHOD_ORDER = ["mioflow", "prescient", "scnode", "cellrank2", "wot"]
METHOD_LABELS = {
    "mioflow": "MIOFlow",
    "prescient": "PRESCIENT",
    "scnode": "scNODE",
    "cellrank2": "CellRank2",
    "wot": "WOT",
}


def _set_fonts() -> None:
    available = {font.name for font in font_manager.fontManager.ttflist}
    for name in ("Microsoft YaHei", "SimHei", "Arial Unicode MS", "DejaVu Sans"):
        if name in available:
            plt.rcParams["font.family"] = name
            break
    plt.rcParams["axes.unicode_minus"] = False


def _method_label(method: str) -> str:
    return METHOD_LABELS.get(str(method), str(method))


def _fmt_float(value: float | int | str, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return str(value)


def _fmt_rank(value: float | int | str) -> str:
    if pd.isna(value):
        return "NA"
    try:
        f = float(value)
        return str(int(f)) if f.is_integer() else f"{f:.1f}"
    except (TypeError, ValueError):
        return str(value)


def _load_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "method" in df.columns:
        df = df.dropna(subset=["method"])
    return df.dropna(how="all")


def _numeric(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for column in columns:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    return out


def _rank_group_mean(df: pd.DataFrame, group_cols: list[str], rank_col: str, out_col: str) -> pd.DataFrame:
    out = df.copy()
    out[out_col] = out.groupby(group_cols)[rank_col].rank(ascending=True, method="min")
    return out


def _rerank_lineage(df: pd.DataFrame) -> pd.DataFrame:
    out = _numeric(
        df,
        [
            "single_step_auc_roc",
            "single_step_auc_prc",
            "single_step_jaccard",
        ],
    )
    for metric in [
        "single_step_auc_roc",
        "single_step_auc_prc",
        "single_step_jaccard",
    ]:
        out[f"rank_{metric}"] = out.groupby(["dataset_id", "scenario"])[metric].rank(
            ascending=False,
            method="min",
        )
    out["lineage_composite_rank_score"] = out[
        [
            "rank_single_step_auc_roc",
            "rank_single_step_auc_prc",
            "rank_single_step_jaccard",
        ]
    ].mean(axis=1)
    out["lineage_composite_rank"] = out.groupby(["dataset_id", "scenario"])[
        "lineage_composite_rank_score"
    ].rank(ascending=True, method="min")
    return out


def _rerank_embedding(df: pd.DataFrame) -> pd.DataFrame:
    out = _numeric(
        df,
        [
            "adjusted_rand_index",
            "pred_tp_avg_normalized_entropy",
        ],
    )
    out["rank_ari"] = out.groupby(["dataset_id", "scenario"])[
        "adjusted_rand_index"
    ].rank(ascending=False, method="min")
    out["rank_classifier_entropy"] = out.groupby(["dataset_id", "scenario"])[
        "pred_tp_avg_normalized_entropy"
    ].rank(ascending=True, method="min")
    out["embedding_composite_rank_score"] = out[
        ["rank_ari", "rank_classifier_entropy"]
    ].mean(axis=1)
    out["embedding_composite_rank"] = out.groupby(["dataset_id", "scenario"])[
        "embedding_composite_rank_score"
    ].rank(ascending=True, method="min")
    return out


def _rerank_forecast(df: pd.DataFrame) -> pd.DataFrame:
    out = _numeric(
        df,
        [
            "wasserstein_distance",
            "gaussian_mmd",
            "energy_distance_mmd",
            "hausdorff_loss",
        ],
    )
    rank_cols = []
    for metric in [
        "wasserstein_distance",
        "gaussian_mmd",
        "energy_distance_mmd",
        "hausdorff_loss",
    ]:
        rank_col = f"rank_{metric}"
        out[rank_col] = out.groupby(["dataset_id", "scenario"])[metric].rank(
            ascending=True,
            method="min",
        )
        rank_cols.append(rank_col)
    out["forecast_composite_rank_score"] = out[rank_cols].mean(axis=1)
    out["forecast_composite_rank"] = out.groupby(["dataset_id", "scenario"])[
        "forecast_composite_rank_score"
    ].rank(ascending=True, method="min")
    return out


def _load_primary_rankings() -> dict[str, pd.DataFrame]:
    return {
        "lineage_rankings": _load_csv(PRIMARY_REPORT_DIR / "official_silver_lineage_rankings.csv"),
        "embedding_rankings": _load_csv(PRIMARY_REPORT_DIR / "official_silver_embedding_rankings.csv"),
        "forecast_rankings": _load_csv(PRIMARY_REPORT_DIR / "official_silver_forecast_rankings.csv"),
    }


def _load_gse242424_rankings() -> dict[str, pd.DataFrame]:
    rankings = _load_csv(GSE242424_REPORT_DIR / "gse242424_oskm_ground_truth_rankings.csv")

    lineage = pd.DataFrame(
        {
            "dataset_id": rankings["dataset_id"],
            "method": rankings["method"],
            "scenario": rankings["scenario"],
            "run_id": rankings["run_id"],
            "single_step_auc_roc": rankings["auroc"],
            "single_step_auc_prc": rankings["auprc"],
            "single_step_jaccard": rankings["jaccard_topk"],
            "lineage_composite_rank_score": rankings["lineage_composite_rank_score"],
            "status": rankings["lineage_status"],
        }
    )
    lineage = _numeric(
        lineage,
        [
            "single_step_auc_roc",
            "single_step_auc_prc",
            "single_step_jaccard",
            "lineage_composite_rank_score",
        ],
    )
    lineage = _rank_group_mean(
        lineage,
        ["dataset_id", "scenario"],
        "lineage_composite_rank_score",
        "lineage_composite_rank",
    )
    ot_summary_path = GSE242424_REPORT_DIR / "gse242424_ot_lineage_summary.csv"
    if ot_summary_path.exists():
        ot = _load_csv(ot_summary_path)
        ot = pd.DataFrame(
            {
                "dataset_id": ot["dataset_id"],
                "method": ot["method"],
                "scenario": ot["scenario"],
                "run_id": ot["run_id"],
                "single_step_auc_roc": ot["single_step_auc_roc"],
                "single_step_auc_prc": ot["single_step_auc_prc"],
                "single_step_jaccard": ot["single_step_jaccard"],
                "status": ot["status"],
            }
        )
        ot = _numeric(
            ot,
            [
                "single_step_auc_roc",
                "single_step_auc_prc",
                "single_step_jaccard",
            ],
        )
        for metric in [
            "single_step_auc_roc",
            "single_step_auc_prc",
            "single_step_jaccard",
        ]:
            ot[f"rank_{metric}"] = ot.groupby(["dataset_id", "scenario"])[metric].rank(
                ascending=False,
                method="min",
            )
        ot["lineage_composite_rank_score"] = ot[
            [
                "rank_single_step_auc_roc",
                "rank_single_step_auc_prc",
                "rank_single_step_jaccard",
            ]
        ].mean(axis=1)
        ot = _rank_group_mean(
            ot,
            ["dataset_id", "scenario"],
            "lineage_composite_rank_score",
            "lineage_composite_rank",
        )
        lineage = pd.concat([lineage, ot], ignore_index=True, sort=False)

    embedding = pd.DataFrame(
        {
            "dataset_id": rankings["dataset_id"],
            "method": rankings["method"],
            "scenario": rankings["scenario"],
            "run_id": rankings["run_id"],
            "adjusted_rand_index": rankings["embedding_ari"],
            "pred_tp_avg_normalized_entropy": rankings["embedding_entropy"],
            "embedding_composite_rank_score": rankings["embedding_composite_rank_score"],
            "status": "completed",
        }
    )
    embedding = _numeric(
        embedding,
        [
            "adjusted_rand_index",
            "pred_tp_avg_normalized_entropy",
            "embedding_composite_rank_score",
        ],
    )
    embedding = _rank_group_mean(
        embedding,
        ["dataset_id", "scenario"],
        "embedding_composite_rank_score",
        "embedding_composite_rank",
    )

    forecast_inventory = _load_csv(
        PRIMARY_REPORT_DIR / "official_silver_forecast_exact_rerun_inventory.csv"
    )
    forecast = forecast_inventory[
        (forecast_inventory["dataset_id"] == "GSE242424")
        & (forecast_inventory["official_forecast_exact"].astype(str).str.lower() == "true")
    ].copy()
    forecast = _numeric(
        forecast,
        [
            "wasserstein_distance",
            "gaussian_mmd",
            "energy_distance_mmd",
            "hausdorff_loss",
        ],
    )
    rank_cols = []
    for metric in [
        "wasserstein_distance",
        "gaussian_mmd",
        "energy_distance_mmd",
        "hausdorff_loss",
    ]:
        rank_col = f"rank_{metric}"
        forecast[rank_col] = forecast.groupby(["dataset_id", "scenario"])[metric].rank(
            ascending=True,
            method="min",
        )
        rank_cols.append(rank_col)
    forecast["forecast_composite_rank_score"] = forecast[rank_cols].mean(axis=1)
    forecast = _rank_group_mean(
        forecast,
        ["dataset_id", "scenario"],
        "forecast_composite_rank_score",
        "forecast_composite_rank",
    )

    return {
        "lineage_rankings": lineage,
        "embedding_rankings": embedding,
        "forecast_rankings": forecast,
    }


def _summarize_lineage(df: pd.DataFrame) -> pd.DataFrame:
    work = _numeric(
        df,
        [
            "lineage_composite_rank_score",
            "single_step_auc_roc",
            "single_step_auc_prc",
            "single_step_jaccard",
        ],
    )
    summary = work.groupby("method", as_index=False).agg(
        lineage_runs=("run_id", "count"),
        mean_lineage_composite_rank_score=("lineage_composite_rank_score", "mean"),
        mean_single_step_auc_roc=("single_step_auc_roc", "mean"),
        mean_single_step_auc_prc=("single_step_auc_prc", "mean"),
        mean_single_step_jaccard=("single_step_jaccard", "mean"),
    )
    summary["overall_lineage_rank"] = summary["mean_lineage_composite_rank_score"].rank(
        ascending=True,
        method="min",
    )
    return summary.sort_values(["overall_lineage_rank", "method"])


def _summarize_embedding(df: pd.DataFrame) -> pd.DataFrame:
    work = _numeric(
        df,
        [
            "embedding_composite_rank_score",
            "adjusted_rand_index",
            "pred_tp_avg_normalized_entropy",
        ],
    )
    summary = work.groupby("method", as_index=False).agg(
        embedding_runs=("run_id", "count"),
        mean_embedding_composite_rank_score=("embedding_composite_rank_score", "mean"),
        mean_ari=("adjusted_rand_index", "mean"),
        median_ari=("adjusted_rand_index", "median"),
        mean_classifier_entropy=("pred_tp_avg_normalized_entropy", "mean"),
    )
    summary["overall_embedding_rank"] = summary["mean_embedding_composite_rank_score"].rank(
        ascending=True,
        method="min",
    )
    return summary.sort_values(["overall_embedding_rank", "method"])


def _summarize_forecast(df: pd.DataFrame) -> pd.DataFrame:
    work = _numeric(
        df,
        [
            "forecast_composite_rank_score",
            "wasserstein_distance",
            "gaussian_mmd",
            "energy_distance_mmd",
            "hausdorff_loss",
        ],
    )
    summary = work.groupby("method", as_index=False).agg(
        forecast_runs=("run_id", "count"),
        mean_forecast_composite_rank_score=("forecast_composite_rank_score", "mean"),
        mean_wasserstein_distance=("wasserstein_distance", "mean"),
        mean_gaussian_mmd=("gaussian_mmd", "mean"),
        mean_energy_distance_mmd=("energy_distance_mmd", "mean"),
        mean_hausdorff_loss=("hausdorff_loss", "mean"),
    )
    summary["overall_forecast_rank"] = summary["mean_forecast_composite_rank_score"].rank(
        ascending=True,
        method="min",
    )
    return summary.sort_values(["overall_forecast_rank", "method"])


def _summarize_combined(
    embedding_summary: pd.DataFrame,
    lineage_summary: pd.DataFrame,
) -> pd.DataFrame:
    combined = embedding_summary.merge(lineage_summary, on="method", how="inner")
    combined["combined_rank_score"] = combined[
        ["mean_embedding_composite_rank_score", "mean_lineage_composite_rank_score"]
    ].mean(axis=1)
    combined["overall_combined_rank"] = combined["combined_rank_score"].rank(
        ascending=True,
        method="min",
    )
    return combined.sort_values(["overall_combined_rank", "method"])


def _load_tables() -> dict[str, pd.DataFrame]:
    primary = _load_primary_rankings()
    gse242424 = _load_gse242424_rankings()

    lineage_rankings = pd.concat(
        [primary["lineage_rankings"], gse242424["lineage_rankings"]],
        ignore_index=True,
        sort=False,
    )
    embedding_rankings = pd.concat(
        [primary["embedding_rankings"], gse242424["embedding_rankings"]],
        ignore_index=True,
        sort=False,
    )
    forecast_rankings = pd.concat(
        [primary["forecast_rankings"], gse242424["forecast_rankings"]],
        ignore_index=True,
        sort=False,
    )
    lineage_rankings = _rerank_lineage(lineage_rankings)
    embedding_rankings = _rerank_embedding(embedding_rankings)
    forecast_rankings = _rerank_forecast(forecast_rankings)

    lineage_summary = _summarize_lineage(lineage_rankings)
    embedding_summary = _summarize_embedding(embedding_rankings)
    forecast_summary = _summarize_forecast(forecast_rankings)

    return {
        "lineage_rankings": lineage_rankings,
        "lineage_summary": lineage_summary,
        "embedding_rankings": embedding_rankings,
        "embedding_summary": embedding_summary,
        "forecast_rankings": forecast_rankings,
        "forecast_summary": forecast_summary,
        "combined_summary": _summarize_combined(embedding_summary, lineage_summary),
    }


def _make_table(
    ax,
    df: pd.DataFrame,
    title: str,
    *,
    col_widths: list[float] | None = None,
    fontsize: float = 8.5,
) -> None:
    ax.axis("off")
    ax.text(0.0, 1.08, title, transform=ax.transAxes, fontsize=13, weight="bold", color="#222222")
    table = ax.table(
        cellText=df.values,
        colLabels=df.columns,
        loc="upper left",
        cellLoc="center",
        colLoc="center",
        colWidths=col_widths,
        bbox=[0.0, 0.0, 1.0, 0.94],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(fontsize)
    for (row, _col), cell in table.get_celld().items():
        cell.set_edgecolor("#d5dbe5")
        cell.set_linewidth(0.7)
        if row == 0:
            cell.set_facecolor("#243447")
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        else:
            cell.set_facecolor("#f7f9fc" if row % 2 == 0 else "white")


def _lineage_table(lineage_summary: pd.DataFrame) -> pd.DataFrame:
    df = lineage_summary.sort_values(["overall_lineage_rank", "method"])
    return pd.DataFrame(
        {
            "Method": df["method"].map(_method_label),
            "Rank": df["overall_lineage_rank"].map(_fmt_rank),
            "AUROC": df["mean_single_step_auc_roc"].map(_fmt_float),
            "AUPRC": df["mean_single_step_auc_prc"].map(_fmt_float),
            "Jaccard": df["mean_single_step_jaccard"].map(_fmt_float),
        }
    )


def _embedding_table(embedding_summary: pd.DataFrame) -> pd.DataFrame:
    df = embedding_summary.sort_values(["overall_embedding_rank", "method"])
    return pd.DataFrame(
        {
            "Method": df["method"].map(_method_label),
            "Rank": df["overall_embedding_rank"].map(_fmt_rank),
            "ARI": df["mean_ari"].map(_fmt_float),
            "Median ARI": df["median_ari"].map(_fmt_float),
            "Entropy": df["mean_classifier_entropy"].map(_fmt_float),
        }
    )


def _forecast_table(forecast_summary: pd.DataFrame) -> pd.DataFrame:
    df = forecast_summary.sort_values(["overall_forecast_rank", "method"])
    return pd.DataFrame(
        {
            "Method": df["method"].map(_method_label),
            "Rank": df["overall_forecast_rank"].map(_fmt_rank),
            "WD": df["mean_wasserstein_distance"].map(lambda x: _fmt_float(x, 1)),
            "Energy MMD": df["mean_energy_distance_mmd"].map(_fmt_float),
            "Hausdorff": df["mean_hausdorff_loss"].map(lambda x: _fmt_float(x, 1)),
        }
    )


def _projection_summary_table(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    combined = tables["combined_summary"].copy()
    forecast = tables["forecast_summary"][
        ["method", "overall_forecast_rank", "mean_wasserstein_distance"]
    ]
    merged = combined.merge(forecast, on="method", how="left").sort_values(
        "overall_combined_rank"
    )
    return pd.DataFrame(
        {
            "Method": merged["method"].map(_method_label),
            "Emb+Lin rank": merged["overall_combined_rank"].map(_fmt_rank),
            "Embedding rank": merged["overall_embedding_rank"].map(_fmt_rank),
            "Lineage rank": merged["overall_lineage_rank"].map(_fmt_rank),
            "Forecast rank": merged["overall_forecast_rank"].map(_fmt_rank),
            "Embedding ARI": merged["mean_ari"].map(_fmt_float),
            "Lineage AUROC": merged["mean_single_step_auc_roc"].map(_fmt_float),
            "Forecast WD": merged["mean_wasserstein_distance"].map(lambda x: _fmt_float(x, 1)),
        }
    )


def _metric_heatmap(
    df: pd.DataFrame,
    *,
    value_col: str,
    method_order: list[str],
) -> pd.DataFrame:
    work = df.copy()
    work[value_col] = pd.to_numeric(work[value_col], errors="coerce")
    work["Method"] = work["method"].map(_method_label)
    work["Dataset/Scenario"] = work["dataset_id"] + "-" + work["scenario"]
    pivot = work.pivot_table(
        index="Method",
        columns="Dataset/Scenario",
        values=value_col,
        aggfunc="mean",
    )
    return pivot.reindex(
        index=[_method_label(m) for m in method_order],
        columns=DATASET_SCENARIO_ORDER,
    )


def _draw_heatmap(
    ax,
    heatmap: pd.DataFrame,
    *,
    title: str,
    cmap_name: str,
    vmin: float,
    vmax: float,
    text_digits: int = 3,
    good_high: bool = True,
) -> None:
    values = heatmap.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(values)
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(color="#f1f3f6")
    ax.imshow(masked, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_title(title, fontsize=13, weight="bold", loc="left", pad=12)
    ax.set_xticks(np.arange(len(heatmap.columns)))
    ax.set_xticklabels(heatmap.columns, rotation=35, ha="right", fontsize=7.5)
    ax.set_yticks(np.arange(len(heatmap.index)))
    ax.set_yticklabels(heatmap.index, fontsize=9.5)
    ax.tick_params(length=0)

    threshold = (vmin + vmax) / 2
    for i in range(heatmap.shape[0]):
        for j in range(heatmap.shape[1]):
            value = heatmap.iloc[i, j]
            if pd.isna(value):
                text = "NA"
                color = "#1f2933"
            else:
                text = f"{float(value):.{text_digits}f}"
                if good_high:
                    color = "white" if value >= threshold else "#1f2933"
                else:
                    color = "white" if value < threshold else "#1f2933"
            ax.text(j, i, text, ha="center", va="center", fontsize=7.2, color=color, weight="bold")

    for spine in ax.spines.values():
        spine.set_visible(False)


def _write_markdown(
    *,
    projection_table: pd.DataFrame,
    lineage_table: pd.DataFrame,
    embedding_table: pd.DataFrame,
    forecast_table: pd.DataFrame,
    tables: dict[str, pd.DataFrame],
) -> None:
    n_lineage = len(tables["lineage_rankings"])
    n_embedding = len(tables["embedding_rankings"])
    n_forecast = len(tables["forecast_rankings"])

    content = f"""# Trajectory 结果可视化汇总

生成自当前本地 benchmark report CSVs。

![Trajectory results visual summary](trajectory_results_visual_summary.png)

## 当前汇总规模

- Lineage Fidelity: {n_lineage} runs
- Embedding Coherence: {n_embedding} runs
- Forecast Accuracy: {n_forecast} exact-mode runs

## 主要结论

1. 当前图已合并 GSE178325、GSE230659 和 GSE242424 三组数据。
2. WOT 和 CellRank2 仍为 lineage-only 方法，因此只进入 Lineage Fidelity，不进入 Embedding、Forecast 或 Emb+Lineage combined rank。
3. Forecast heatmap 使用每个 dataset/scenario 内四个 forecast loss 排名的平均值，即 scTimeBench-style averaged rank。
4. 顶部 method summary 是合并后重新计算的跨场景均值和总体排名。

## Projection-capable Combined Summary

{projection_table.to_markdown(index=False)}

## Lineage Fidelity Summary

{lineage_table.to_markdown(index=False)}

## Embedding Coherence Summary

{embedding_table.to_markdown(index=False)}

## Forecast Accuracy Summary

{forecast_table.to_markdown(index=False)}
"""
    MD_OUT.write_text(content, encoding="utf-8")


def main() -> None:
    _set_fonts()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tables = _load_tables()

    projection_table = _projection_summary_table(tables)
    lineage_table = _lineage_table(tables["lineage_summary"])
    embedding_table = _embedding_table(tables["embedding_summary"])
    forecast_table = _forecast_table(tables["forecast_summary"])

    lineage_hm = _metric_heatmap(
        tables["lineage_rankings"],
        value_col="single_step_auc_roc",
        method_order=LINEAGE_METHOD_ORDER,
    )
    embedding_hm = _metric_heatmap(
        tables["embedding_rankings"],
        value_col="adjusted_rand_index",
        method_order=PROJECTION_METHOD_ORDER,
    )
    forecast_hm = _metric_heatmap(
        tables["forecast_rankings"],
        value_col="forecast_composite_rank_score",
        method_order=PROJECTION_METHOD_ORDER,
    )

    export_table = pd.concat(
        [
            projection_table.assign(Table="Projection-capable combined"),
            lineage_table.assign(Table="Lineage Fidelity"),
            embedding_table.assign(Table="Embedding Coherence"),
            forecast_table.assign(Table="Forecast Accuracy"),
        ],
        ignore_index=True,
        sort=False,
    )
    export_table.to_csv(TABLE_OUT, index=False, encoding="utf-8-sig")

    n_lineage = len(tables["lineage_rankings"])
    n_embedding = len(tables["embedding_rankings"])
    n_forecast = len(tables["forecast_rankings"])

    fig = plt.figure(figsize=(21, 12.5), dpi=180)
    gs = fig.add_gridspec(
        nrows=3,
        ncols=3,
        height_ratios=[0.38, 1.85, 3.25],
        hspace=0.45,
        wspace=0.24,
    )

    ax_title = fig.add_subplot(gs[0, :])
    ax_title.axis("off")
    ax_title.text(
        0.0,
        0.8,
        "Trajectory benchmark summary",
        fontsize=23,
        weight="bold",
        color="#15202b",
        transform=ax_title.transAxes,
    )

    _make_table(
        fig.add_subplot(gs[1, 0]),
        lineage_table,
        "Lineage Fidelity",
        col_widths=[0.24, 0.11, 0.18, 0.18, 0.18],
        fontsize=8.1,
    )
    _make_table(
        fig.add_subplot(gs[1, 1]),
        embedding_table,
        "Embedding Coherence",
        col_widths=[0.24, 0.11, 0.16, 0.21, 0.18],
        fontsize=8.1,
    )
    _make_table(
        fig.add_subplot(gs[1, 2]),
        forecast_table,
        "Forecast Accuracy",
        col_widths=[0.24, 0.11, 0.13, 0.22, 0.19],
        fontsize=8.1,
    )

    _draw_heatmap(
        fig.add_subplot(gs[2, 0]),
        lineage_hm,
        title="Lineage AUROC",
        cmap_name="YlGnBu",
        vmin=0.4,
        vmax=max(0.9, float(np.nanmax(lineage_hm.to_numpy(dtype=float)))),
        text_digits=3,
        good_high=True,
    )
    _draw_heatmap(
        fig.add_subplot(gs[2, 1]),
        embedding_hm,
        title="Embedding ARI",
        cmap_name="YlGnBu",
        vmin=0.0,
        vmax=max(0.35, float(np.nanmax(embedding_hm.to_numpy(dtype=float)))),
        text_digits=3,
        good_high=True,
    )
    _draw_heatmap(
        fig.add_subplot(gs[2, 2]),
        forecast_hm,
        title="Forecast average rank",
        cmap_name="RdYlGn_r",
        vmin=1.0,
        vmax=max(3.0, float(np.nanmax(forecast_hm.to_numpy(dtype=float)))),
        text_digits=2,
        good_high=False,
    )

    fig.text(
        0.01,
        0.01,
        "Source: benchmark report CSVs. Rank: lower is better. ARI/AUROC: higher is better. Forecast losses: lower is better.",
        fontsize=8.5,
        color="#52616b",
    )
    fig.savefig(PNG_OUT, bbox_inches="tight", facecolor="white")
    fig.savefig(PDF_OUT, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    _write_markdown(
        projection_table=projection_table,
        lineage_table=lineage_table,
        embedding_table=embedding_table,
        forecast_table=forecast_table,
        tables=tables,
    )
    print(f"Wrote {PNG_OUT}")
    print(f"Wrote {PDF_OUT}")
    print(f"Wrote {MD_OUT}")
    print(f"Wrote {TABLE_OUT}")


if __name__ == "__main__":
    main()
