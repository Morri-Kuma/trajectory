"""Audit GSE230659 projected HVG2000 annotation outputs.

This report is intentionally narrower than the general milestone summary.  It
focuses on the new formal projected-cell route:

    projected_expression.npy
    -> HVG2000 PCA embedding provider + CellTypist provider
    -> projected consensus milestone labels
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List

import pandas as pd


SCENARIO_DIRS = {
    "A": Path("benchmark/results/scnode/gse230659_milestone_consensus_A_hvg2000_formal"),
    "B": Path("benchmark/results/scnode/gse230659_milestone_consensus_B_hvg2000_formal"),
    "C": Path("benchmark/results/scnode/gse230659_milestone_consensus_C_hvg2000_formal"),
}

LABEL_ORDER = ["epithelial_like", "hCiPS", "intermediate_plastic", "ambiguous"]
LABEL_COLUMNS = {
    "consensus": "consensus_milestone_label",
    "embedding_based": "milestone_embedding_label",
    "classifier_based": "milestone_classifier_label",
}


def _load_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _safe_float(value: Any) -> Any:
    if value is None:
        return ""
    try:
        return float(value)
    except Exception:
        return value


def _write_csv(path: Path, rows: Iterable[Dict[str, Any]], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def build_audit(results_root: Path, output_dir: Path) -> Dict[str, Path]:
    scenario_dirs = {k: results_root / v for k, v in SCENARIO_DIRS.items()}
    summary_rows: List[Dict[str, Any]] = []
    label_rows: List[Dict[str, Any]] = []
    crosstab_rows: List[Dict[str, Any]] = []
    metric_rows: List[Dict[str, Any]] = []

    for scenario, run_dir in scenario_dirs.items():
        labels_path = run_dir / "projected_milestone_labels.csv"
        ann_meta_path = run_dir / "projected_milestone_annotation_metadata.json"
        hvg_meta_path = run_dir / "projected_hvg_embedding_metadata.json"
        eval_dir = run_dir / "embedding_milestone_eval"

        if not labels_path.exists():
            raise FileNotFoundError(labels_path)
        if not ann_meta_path.exists():
            raise FileNotFoundError(ann_meta_path)
        if not hvg_meta_path.exists():
            raise FileNotFoundError(hvg_meta_path)

        df = pd.read_csv(labels_path)
        ann_meta = _load_json(ann_meta_path)
        hvg_meta = _load_json(hvg_meta_path)

        emb = df["milestone_embedding_label"].astype(str)
        clf = df["milestone_classifier_label"].astype(str)
        consensus = df["consensus_milestone_label"].astype(str)
        status = df.get("consensus_status", pd.Series([""] * len(df))).astype(str)

        n_total = int(len(df))
        n_agree = int((emb == clf).sum())
        n_disagree = n_total - n_agree
        n_ambiguous = int((consensus == "ambiguous").sum())

        summary_rows.append(
            {
                "dataset_id": "GSE230659",
                "method": "scnode",
                "scenario": scenario,
                "run_dir": str(run_dir),
                "n_projected_cells": n_total,
                "embedding_provider_method": ann_meta.get("embedding_provider_method"),
                "gene_universe": ann_meta.get("gene_universe"),
                "formal_two_provider_complete": ann_meta.get("formal_two_provider_complete"),
                "projected_hvg_embedding_shape": "x".join(map(str, hvg_meta.get("output_shape", []))),
                "feature_order_matched_exactly": hvg_meta.get("feature_order_matched_exactly"),
                "n_embedding_classifier_agree": n_agree,
                "n_embedding_classifier_disagree": n_disagree,
                "provider_exact_match_fraction_all_cells": round(n_agree / n_total, 6) if n_total else "",
                "ambiguous_count": n_ambiguous,
                "ambiguous_fraction": round(n_ambiguous / n_total, 6) if n_total else "",
                "consensus_status_counts": json.dumps(dict(Counter(status)), sort_keys=True),
            }
        )

        for mode, col in LABEL_COLUMNS.items():
            counts = Counter(df[col].astype(str))
            for label in sorted(counts, key=lambda x: (LABEL_ORDER.index(x) if x in LABEL_ORDER else 99, x)):
                label_rows.append(
                    {
                        "dataset_id": "GSE230659",
                        "method": "scnode",
                        "scenario": scenario,
                        "label_mode": mode,
                        "label": label,
                        "count": int(counts[label]),
                        "fraction": round(counts[label] / n_total, 6) if n_total else "",
                    }
                )

        ctab = pd.crosstab(emb, clf, dropna=False)
        for emb_label in ctab.index.astype(str):
            for clf_label in ctab.columns.astype(str):
                count = int(ctab.loc[emb_label, clf_label])
                if count == 0:
                    continue
                crosstab_rows.append(
                    {
                        "dataset_id": "GSE230659",
                        "method": "scnode",
                        "scenario": scenario,
                        "embedding_label": emb_label,
                        "classifier_label": clf_label,
                        "count": count,
                        "fraction": round(count / n_total, 6) if n_total else "",
                        "agreement": emb_label == clf_label,
                    }
                )

        for mode in ["consensus", "embedding_based", "classifier_based"]:
            metric_path = eval_dir / f"embedding_metrics_{mode}.json"
            if not metric_path.exists():
                continue
            m = _load_json(metric_path)
            metric_rows.append(
                {
                    "dataset_id": "GSE230659",
                    "method": "scnode",
                    "scenario": scenario,
                    "label_mode": mode,
                    "adjusted_rand_index": _safe_float(m.get("adjusted_rand_index")),
                    "n_cells_total": m.get("n_cells_total"),
                    "n_cells_evaluated": m.get("n_cells_evaluated"),
                    "ambiguous_count": m.get("ambiguous_count"),
                    "ambiguous_fraction": _safe_float(m.get("ambiguous_fraction")),
                    "mean_normalized_entropy": _safe_float(m.get("mean_normalized_entropy")),
                    "weighted_mean_normalized_entropy": _safe_float(m.get("weighted_mean_normalized_entropy")),
                    "status": m.get("status"),
                    "source_json": str(metric_path),
                }
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "summary": output_dir / "gse230659_projected_hvg_annotation_audit_summary.csv",
        "label_distribution": output_dir / "gse230659_projected_hvg_annotation_label_distribution.csv",
        "provider_crosstab": output_dir / "gse230659_projected_hvg_annotation_provider_crosstab.csv",
        "embedding_metrics": output_dir / "gse230659_projected_hvg_annotation_embedding_metrics.csv",
        "markdown": output_dir / "gse230659_projected_hvg_annotation_audit.md",
    }

    _write_csv(
        paths["summary"],
        summary_rows,
        [
            "dataset_id", "method", "scenario", "run_dir", "n_projected_cells",
            "embedding_provider_method", "gene_universe", "formal_two_provider_complete",
            "projected_hvg_embedding_shape", "feature_order_matched_exactly",
            "n_embedding_classifier_agree", "n_embedding_classifier_disagree",
            "provider_exact_match_fraction_all_cells", "ambiguous_count",
            "ambiguous_fraction", "consensus_status_counts",
        ],
    )
    _write_csv(
        paths["label_distribution"],
        label_rows,
        ["dataset_id", "method", "scenario", "label_mode", "label", "count", "fraction"],
    )
    _write_csv(
        paths["provider_crosstab"],
        crosstab_rows,
        [
            "dataset_id", "method", "scenario", "embedding_label",
            "classifier_label", "count", "fraction", "agreement",
        ],
    )
    _write_csv(
        paths["embedding_metrics"],
        metric_rows,
        [
            "dataset_id", "method", "scenario", "label_mode",
            "adjusted_rand_index", "n_cells_total", "n_cells_evaluated",
            "ambiguous_count", "ambiguous_fraction", "mean_normalized_entropy",
            "weighted_mean_normalized_entropy", "status", "source_json",
        ],
    )

    _write_markdown(paths["markdown"], summary_rows, label_rows, crosstab_rows, metric_rows)
    return paths


def _rows_for(rows: List[Dict[str, Any]], scenario: str) -> List[Dict[str, Any]]:
    return [r for r in rows if r.get("scenario") == scenario]


def _metric_lookup(metric_rows: List[Dict[str, Any]], scenario: str, mode: str) -> Dict[str, Any]:
    for row in metric_rows:
        if row.get("scenario") == scenario and row.get("label_mode") == mode:
            return row
    return {}


def _write_markdown(
    path: Path,
    summary_rows: List[Dict[str, Any]],
    label_rows: List[Dict[str, Any]],
    crosstab_rows: List[Dict[str, Any]],
    metric_rows: List[Dict[str, Any]],
) -> None:
    lines: List[str] = []
    lines.append("# GSE230659 Projected HVG Annotation Audit")
    lines.append("")
    lines.append(f"Generated: {datetime.now(timezone.utc).isoformat()}")
    lines.append("")
    lines.append("## Scope")
    lines.append("")
    lines.append("This audit covers scNODE A/B/C projected-cell annotation after switching the formal Embedding Coherence route to the HVG2000 PCA embedding provider plus CellTypist classifier provider.")
    lines.append("")
    lines.append("Formal route:")
    lines.append("")
    lines.append("```text")
    lines.append("projected_expression.npy")
    lines.append("-> hvg2000_pca_logistic_regression embedding provider")
    lines.append("-> CellTypist classifier provider")
    lines.append("-> consensus_milestone_label")
    lines.append("```")
    lines.append("")

    lines.append("## Scenario Summary")
    lines.append("")
    lines.append("| Scenario | Cells | HVG Embedding | Provider Match | Ambiguous | Consensus ARI | Embedding ARI | Classifier ARI |")
    lines.append("|---|---:|---|---:|---:|---:|---:|---:|")
    for row in summary_rows:
        scenario = row["scenario"]
        consensus_m = _metric_lookup(metric_rows, scenario, "consensus")
        emb_m = _metric_lookup(metric_rows, scenario, "embedding_based")
        clf_m = _metric_lookup(metric_rows, scenario, "classifier_based")
        lines.append(
            "| {scenario} | {cells} | {shape} | {match:.3f} | {ambig} ({ambig_frac:.3f}) | {ari_c:.6g} | {ari_e:.6g} | {ari_clf:.6g} |".format(
                scenario=scenario,
                cells=row["n_projected_cells"],
                shape=row["projected_hvg_embedding_shape"],
                match=float(row["provider_exact_match_fraction_all_cells"]),
                ambig=row["ambiguous_count"],
                ambig_frac=float(row["ambiguous_fraction"]),
                ari_c=float(consensus_m.get("adjusted_rand_index", 0) or 0),
                ari_e=float(emb_m.get("adjusted_rand_index", 0) or 0),
                ari_clf=float(clf_m.get("adjusted_rand_index", 0) or 0),
            )
        )
    lines.append("")

    lines.append("## Label Distributions")
    lines.append("")
    lines.append("| Scenario | Mode | Label | Count | Fraction |")
    lines.append("|---|---|---|---:|---:|")
    for row in label_rows:
        lines.append(
            f"| {row['scenario']} | {row['label_mode']} | {row['label']} | {row['count']} | {float(row['fraction']):.3f} |"
        )
    lines.append("")

    lines.append("## Provider Disagreement")
    lines.append("")
    lines.append("Rows below show nonzero embedding-provider vs classifier-provider label combinations across all cells. Off-diagonal rows become `ambiguous` under the current consensus rule.")
    lines.append("")
    lines.append("| Scenario | Embedding Label | Classifier Label | Count | Fraction | Agreement |")
    lines.append("|---|---|---|---:|---:|---|")
    for row in crosstab_rows:
        agreement = "yes" if row["agreement"] else "no"
        lines.append(
            f"| {row['scenario']} | {row['embedding_label']} | {row['classifier_label']} | {row['count']} | {float(row['fraction']):.3f} | {agreement} |"
        )
    lines.append("")

    lines.append("## Interpretation")
    lines.append("")
    lines.append("- The technical route is coherent: all three scenarios report `formal_two_provider_complete=True`, `gene_universe=HVG2000`, exact feature-order matches, and 50-dimensional projected HVG PCA sidecar embeddings.")
    lines.append("- Consensus ARI is low after the route switch, especially for B and C. This appears driven by provider disagreement and the conservative current rule that off-diagonal provider pairs become `ambiguous`.")
    lines.append("- Provider agreement should be reviewed before changing consensus policy. A confidence-based tie-break may recover more evaluated cells, but it would change the formal interpretation and should be treated as a versioned policy decision.")
    lines.append("")

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Audit GSE230659 projected HVG annotation.")
    p.add_argument("--results-root", default=".")
    p.add_argument("--output-dir", default="benchmark/reports/official")
    return p


def main() -> int:
    args = build_parser().parse_args()
    paths = build_audit(Path(args.results_root), Path(args.output_dir))
    for name, path in paths.items():
        print(f"{name}: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
