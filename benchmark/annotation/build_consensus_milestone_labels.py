"""Merge embedding/classifier milestone labels into consensus labels."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
from typing import Any, Dict, List


TIE_BREAK_CHOICES = ("majority", "embedding", "classifier", "marker", "abstain")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Build consensus milestone labels.")
    p.add_argument("--input-h5ad", required=True)
    p.add_argument("--output-h5ad", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--dataset-id", required=True)
    p.add_argument("--embedding-label-key", default="milestone_embedding_label")
    p.add_argument("--classifier-label-key", default="milestone_classifier_label")
    p.add_argument("--marker-label-key", default="milestone_marker_label")
    p.add_argument("--embedding-confidence-key", default="milestone_embedding_confidence")
    p.add_argument("--classifier-confidence-key", default="milestone_classifier_confidence")
    p.add_argument("--tie-break", choices=TIE_BREAK_CHOICES, default="majority")
    p.add_argument("--ambiguous-label", default="ambiguous")
    p.add_argument("--dry-run", action="store_true")
    return p


def _as_float_array(series, default: float = 0.0):
    import numpy as np

    try:
        return series.astype(float).to_numpy()
    except Exception:
        return np.full(len(series), default, dtype=float)


def _write_counts(path: Path, labels) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["label", "count"])
        writer.writeheader()
        for label, count in sorted(Counter(labels).items()):
            writer.writerow({"label": label, "count": count})


def _write_qc(path: Path, rows: List[Dict[str, Any]]) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = build_parser().parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        with open(out_dir / "consensus_milestone_labels_metadata.json", "w", encoding="utf-8") as f:
            json.dump({"status": "dry_run", "dataset_id": args.dataset_id}, f, indent=2)
        print("[build_consensus_milestone_labels] dry-run metadata written")
        return 0

    import anndata as ad
    import numpy as np

    print(f"[build_consensus_milestone_labels] loading {args.input_h5ad}")
    adata = ad.read_h5ad(args.input_h5ad)
    for key in (args.embedding_label_key, args.classifier_label_key):
        if key not in adata.obs:
            raise KeyError(f"required obs column missing: {key}")

    emb = adata.obs[args.embedding_label_key].astype(str).to_numpy()
    clf = adata.obs[args.classifier_label_key].astype(str).to_numpy()
    marker = (
        adata.obs[args.marker_label_key].astype(str).to_numpy()
        if args.marker_label_key in adata.obs
        else np.array([args.ambiguous_label] * adata.n_obs, dtype=object)
    )
    emb_conf = (
        _as_float_array(adata.obs[args.embedding_confidence_key])
        if args.embedding_confidence_key in adata.obs
        else np.zeros(adata.n_obs)
    )
    clf_conf = (
        _as_float_array(adata.obs[args.classifier_confidence_key])
        if args.classifier_confidence_key in adata.obs
        else np.zeros(adata.n_obs)
    )

    consensus: List[str] = []
    confidence: List[float] = []
    source: List[str] = []
    for e, c, m, ec, cc in zip(emb, clf, marker, emb_conf, clf_conf):
        votes = Counter([e, c, m])
        top_label, top_count = votes.most_common(1)[0]
        if args.tie_break == "majority" and top_count >= 2:
            consensus.append(top_label)
            conf_parts = []
            if e == top_label:
                conf_parts.append(ec)
            if c == top_label:
                conf_parts.append(cc)
            confidence.append(float(sum(conf_parts) / len(conf_parts)) if conf_parts else float(max(ec, cc)))
            source.append("three_provider_majority")
        elif e == c:
            consensus.append(e)
            confidence.append(float((ec + cc) / 2.0))
            source.append("embedding_classifier_agree")
        elif args.tie_break == "embedding":
            consensus.append(e)
            confidence.append(float(ec))
            source.append("tie_break_embedding")
        elif args.tie_break == "classifier":
            consensus.append(c)
            confidence.append(float(cc))
            source.append("tie_break_classifier")
        elif args.tie_break == "marker":
            consensus.append(m)
            confidence.append(float(max(ec, cc)))
            source.append("tie_break_marker")
        else:
            consensus.append(args.ambiguous_label)
            confidence.append(float(max(ec, cc)))
            source.append("discordant_abstain")

    adata.obs["consensus_milestone_label"] = consensus
    adata.obs["consensus_confidence"] = np.asarray(confidence, dtype="float32")
    adata.obs["consensus_label_source"] = source

    agreement = emb == clf
    qc_rows = [
        {
            "cell_id": str(cid),
            "milestone_embedding_label": str(e),
            "milestone_classifier_label": str(c),
            "milestone_marker_label": str(m),
            "consensus_milestone_label": str(co),
            "consensus_confidence": float(cf),
            "consensus_label_source": str(src),
            "embedding_classifier_agree": bool(ag),
        }
        for cid, e, c, m, co, cf, src, ag in zip(
            adata.obs_names, emb, clf, marker, consensus, confidence, source, agreement
        )
    ]

    print(f"[build_consensus_milestone_labels] writing {args.output_h5ad}")
    adata.write_h5ad(args.output_h5ad)
    _write_qc(out_dir / "consensus_label_agreement_qc.csv", qc_rows)
    _write_counts(out_dir / "consensus_label_counts.csv", consensus)

    meta: Dict[str, Any] = {
        "script": "build_consensus_milestone_labels",
        "status": "completed",
        "dataset_id": args.dataset_id,
        "input_h5ad": args.input_h5ad,
        "output_h5ad": args.output_h5ad,
        "embedding_label_key": args.embedding_label_key,
        "classifier_label_key": args.classifier_label_key,
        "marker_label_key": args.marker_label_key,
        "tie_break": args.tie_break,
        "n_cells": int(adata.n_obs),
        "n_embedding_classifier_agree": int(np.sum(agreement)),
        "agreement_fraction": float(np.mean(agreement)),
        "consensus_label_counts": dict(Counter(consensus)),
        "source_counts": dict(Counter(source)),
    }
    with open(out_dir / "consensus_milestone_labels_metadata.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
        f.write("\n")

    print("[build_consensus_milestone_labels] completed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
