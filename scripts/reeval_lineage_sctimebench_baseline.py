"""
Re-evaluate existing Lineage Fidelity outputs with the formal scTimeBench-style
correlation baseline.

This script does not rerun WOT, CellRank2, or scNODE. It reads existing
state_transition_matrix.csv and lineage_graph_edges.csv files, rebuilds the
scenario-filtered AnnData used by the run, and calls run_lineage_evaluation().

Example
-------
python scripts/reeval_lineage_sctimebench_baseline.py \
  --run benchmark/configs/wot_gse230659_observed_scgpt_v1.yaml=benchmark/results/wot/scenario_A_scgpt_v1 \
  --run benchmark/configs/wot_gse230659_observed_scgpt_v1_scenarioB.yaml=benchmark/results/wot/scenario_B_scgpt_v1 \
  --run benchmark/configs/wot_gse230659_observed_scgpt_v1_scenarioC.yaml=benchmark/results/wot/scenario_C_scgpt_v1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable


def _project_root() -> Path:
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "benchmark").exists() and (candidate / "scripts").exists():
            return candidate
    return here.parent


def _resolve(root: Path, path: str | Path) -> Path:
    p = Path(path)
    return p if p.is_absolute() else root / p


def _load_ground_truth_from_config(method_config: dict):
    from benchmark.ground_truth import load_ground_truth

    spec = load_ground_truth(method_config)
    ref_path = spec.reference_graph_path
    if ref_path in ("null", "~", ""):
        ref_path = None
    return spec, ref_path


def _scenario_filtered_adata(adata, method_config: dict, time_key: str):
    """Return the AnnData cell universe used for this scenario's lineage eval."""
    scenario_params = method_config.get("scenario_params") or {}
    train_times = scenario_params.get("train_times") or []
    if not train_times:
        return adata

    train_times_float = {float(t) for t in train_times}
    obs_times = adata.obs[time_key].astype(float)
    mask = obs_times.isin(train_times_float).to_numpy()
    return adata[mask].to_memory()


def _iter_runs(entries: Iterable[str]):
    for entry in entries:
        if "=" not in entry:
            raise SystemExit(f"Bad --run entry {entry!r}; expected CONFIG=RUN_DIR")
        config, run_dir = entry.split("=", 1)
        yield config, run_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--adata",
        default="benchmark/results/scgpt/full/adata_scgpt_annotated.h5ad",
        help="Annotated full AnnData used by WOT/CellRank2 official runs.",
    )
    parser.add_argument(
        "--run",
        action="append",
        required=True,
        help="CONFIG=RUN_DIR pair. May be repeated.",
    )
    parser.add_argument(
        "--source-chunk-size",
        type=int,
        default=512,
        help="Source-cell chunk size for Spearman baseline computation.",
    )
    args = parser.parse_args()

    root = _project_root()
    sys.path.insert(0, str(root))

    import anndata as ad  # noqa: E402
    import yaml  # noqa: E402
    from benchmark.evaluation import eval_lineage  # noqa: E402

    adata_path = _resolve(root, args.adata)
    print(f"[reeval-sctimebench] Loading adata backed='r': {adata_path}")
    adata_full = ad.read_h5ad(adata_path, backed="r")
    print(f"[reeval-sctimebench] full adata: {adata_full.n_obs} x {adata_full.n_vars}")

    # Patch chunk size for this process without changing the evaluator API.
    original_compute = eval_lineage.compute_correlation_baseline

    def compute_with_chunk(*a, **kw):
        kw.setdefault("source_chunk_size", args.source_chunk_size)
        return original_compute(*a, **kw)

    eval_lineage.compute_correlation_baseline = compute_with_chunk

    for config_entry, run_entry in _iter_runs(args.run):
        config_path = _resolve(root, config_entry)
        run_dir = _resolve(root, run_entry)
        print("\n" + "=" * 72)
        print(f"[reeval-sctimebench] config : {config_path}")
        print(f"[reeval-sctimebench] run dir: {run_dir}")

        with open(config_path, encoding="utf-8-sig") as f:
            method_config = yaml.safe_load(f)

        dataset_cfg = method_config.get("dataset") or {}
        time_key = dataset_cfg.get("time_key", "abs_day")
        ground_truth_spec, reference_graph_path = _load_ground_truth_from_config(method_config)
        reference_graph = _resolve(root, reference_graph_path)
        eval_adata = _scenario_filtered_adata(adata_full, method_config, time_key)
        print(f"[reeval-sctimebench] eval adata: {eval_adata.n_obs} x {eval_adata.n_vars}")
        print(f"[reeval-sctimebench] time_key={time_key!r}, state_key={ground_truth_spec.state_key!r}")

        stm = run_dir / "state_transition_matrix.csv"
        edges = run_dir / "lineage_graph_edges.csv"
        if not stm.exists() or not edges.exists():
            raise FileNotFoundError(f"Missing state_transition_matrix.csv or lineage_graph_edges.csv in {run_dir}")

        metrics = eval_lineage.run_lineage_evaluation(
            state_transition_matrix_path=str(stm),
            lineage_graph_edges_path=str(edges),
            output_dir=str(run_dir),
            reference_graph_path=str(reference_graph),
            adata=eval_adata,
            edge_confidence_mode=ground_truth_spec.confidence_mode,
            exclude_uncertain_states=ground_truth_spec.exclude_uncertain_states,
            cell_state_key=ground_truth_spec.state_key,
            time_key=time_key,
            ground_truth=ground_truth_spec.to_dict(),
        )
        baseline = metrics.get("baseline") or {}
        print(
            "[reeval-sctimebench] baseline: "
            f"style={baseline.get('baseline_style')}, "
            f"auroc={baseline.get('auroc')}, "
            f"auprc={baseline.get('auprc')}, "
            f"pairs={baseline.get('n_timepoint_pairs_used')}, "
            f"votes={baseline.get('n_source_cell_votes')}"
        )


if __name__ == "__main__":
    main()
