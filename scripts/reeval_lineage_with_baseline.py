"""
reeval_lineage_with_baseline.py
--------------------------------
Re-run eval_lineage.run_lineage_evaluation() on existing method outputs so
that lineage_metrics.json picks up the new correlation baseline and the new
jaccard_similarity_topk field.

This script never calls the method adapters; it only re-evaluates existing
state_transition_matrix.csv / lineage_graph_edges.csv files against the
reference graph. Safe to re-run.

Usage
-----
    python scripts/reeval_lineage_with_baseline.py \
        --adata benchmark/results/scgpt/full/adata_scgpt_annotated.h5ad \
        --runs \
            wot=benchmark/results/wot/scenario_A_scgpt_v1 \
            cellrank2=benchmark/results/cellrank2/scenario_A_scgpt_v1 \
        --reference-graph benchmark/datasets/scgpt_reference_graph_v1.json \
        --edge-confidence-mode medium_and_above \
        --cell-state-key scgpt_pseudostate_provisional
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


def _project_root() -> Path:
    here = Path(__file__).resolve().parent
    for c in [here, *here.parents]:
        if (c / "benchmark").exists():
            return c
    return here.parent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--adata", required=True)
    parser.add_argument("--runs", nargs="+", required=True,
                        help="method=output_dir entries (space-separated)")
    parser.add_argument("--reference-graph", required=True)
    parser.add_argument("--edge-confidence-mode", default="medium_and_above")
    parser.add_argument("--cell-state-key",
                        default="scgpt_pseudostate_provisional")
    parser.add_argument("--time-key", default=None,
                        help="Time column for scTimeBench-style adjacent-pair baseline.")
    parser.add_argument("--exclude-uncertain-states", action="store_true",
                        default=False)
    args = parser.parse_args()

    root = _project_root()
    sys.path.insert(0, str(root))

    import anndata as ad  # noqa: E402
    from benchmark.evaluation.eval_lineage import run_lineage_evaluation  # noqa: E402

    adata_path = (root / args.adata) if not Path(args.adata).is_absolute() else Path(args.adata)
    ref_path = (root / args.reference_graph) if not Path(args.reference_graph).is_absolute() else Path(args.reference_graph)

    # Use backed='r' so we never materialise the full n_cells × n_genes matrix.
    # _state_mean_expression chunks over cells, so a backed CSR handle is
    # sufficient and keeps memory flat at the 75k×23k GSE230659 scale.
    print(f"[reeval] Loading adata (backed='r'): {adata_path}")
    adata = ad.read_h5ad(adata_path, backed="r")
    print(f"[reeval] adata: {adata.n_obs} cells × {adata.n_vars} genes")

    for entry in args.runs:
        if "=" not in entry:
            raise SystemExit(f"Bad --runs entry {entry!r}; expected method=output_dir")
        method, out_dir = entry.split("=", 1)
        out_dir_p = (root / out_dir) if not Path(out_dir).is_absolute() else Path(out_dir)
        stm = out_dir_p / "state_transition_matrix.csv"
        edges = out_dir_p / "lineage_graph_edges.csv"
        if not stm.exists() or not edges.exists():
            print(f"[reeval] SKIP {method}: missing STM or edges in {out_dir_p}")
            continue
        print(f"\n[reeval] {method}  →  {out_dir_p}")
        run_lineage_evaluation(
            state_transition_matrix_path=str(stm),
            lineage_graph_edges_path=str(edges),
            output_dir=str(out_dir_p),
            reference_graph_path=str(ref_path),
            adata=adata,
            edge_confidence_mode=args.edge_confidence_mode,
            exclude_uncertain_states=args.exclude_uncertain_states,
            cell_state_key=args.cell_state_key,
            time_key=args.time_key,
        )


if __name__ == "__main__":
    main()
