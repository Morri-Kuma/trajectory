"""
Backfill ground-truth provider metadata into existing result JSON files.

This is for runs produced before the ground-truth provider layer existed. It
does not change metric values; it only adds the provider block used by newer
dispatcher/evaluator runs.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path


def _project_root() -> Path:
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "benchmark").exists() and (
            candidate / "benchmark" / "results"
        ).exists():
            return candidate
    return here.parent


ROOT = _project_root()
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from benchmark.ground_truth import load_ground_truth


def _read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with open(path, encoding="utf-8") as f:
        return json.load(f) or {}


def _write_json(path: Path, data: dict) -> None:
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)


def _looks_like_scgpt_v1(run_dir: Path, metrics: dict, meta: dict) -> bool:
    haystack = " ".join(
        str(x)
        for x in (
            run_dir.name,
            metrics.get("cell_state_key"),
            metrics.get("edge_confidence_mode"),
            meta.get("cell_state_key"),
        )
    )
    return (
        "scgpt_v1" in haystack
        or "scgpt_pseudostate_provisional" in haystack
    )


def main() -> None:
    root = ROOT
    gt = load_ground_truth(provider_id="scgpt_v1").to_dict()
    results_root = root / "benchmark" / "results"
    updated = 0

    for lineage_path in results_root.glob("*/*/lineage_metrics.json"):
        run_dir = lineage_path.parent
        metrics = _read_json(lineage_path)
        meta_path = run_dir / "run_metadata.json"
        meta = _read_json(meta_path)
        if metrics.get("ground_truth"):
            continue
        if not _looks_like_scgpt_v1(run_dir, metrics, meta):
            continue

        metrics["ground_truth"] = gt
        _write_json(lineage_path, metrics)
        if meta_path.exists():
            meta["ground_truth"] = gt
            _write_json(meta_path, meta)
        updated += 1
        print(f"updated {lineage_path}")

    print(f"Backfilled ground_truth metadata for {updated} runs.")


if __name__ == "__main__":
    main()
