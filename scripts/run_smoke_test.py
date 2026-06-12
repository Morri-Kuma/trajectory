"""Local smoke test: generate synthetic data, run the annotation-robustness
pipeline in test_mode, and assert that all expected outputs were produced.

Run from the repo root:
    python scripts/run_smoke_test.py

Exits non-zero if any expected artifact is missing.
"""
from __future__ import annotations

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)


def main() -> int:
    from src.utils import load_config, resolve_path
    import scripts.make_synthetic_smoke_data as gen
    from src import pipeline

    cfg = load_config(os.path.join(ROOT, "config.yaml"))
    assert cfg["mode"] == "test", "run_smoke_test expects mode: test in config.yaml"

    syn_dir = resolve_path(cfg, cfg["active"]["paths"]["synthetic_dir"])
    os.makedirs(syn_dir, exist_ok=True)

    # 1) generate tiny inputs if absent
    if not os.path.exists(os.path.join(syn_dir, "reference.h5ad")):
        sys.argv = ["make_synthetic_smoke_data.py", "--out-dir", syn_dir,
                    "--seed", str(cfg["project"]["random_seed"])]
        gen.main()

    # 2) run the pipeline
    summary = pipeline.run(os.path.join(ROOT, "config.yaml"))

    # 3) assert expected outputs
    out_dir = resolve_path(cfg, cfg["active"]["paths"]["output_dir"])
    expected_per_query = [
        "composition_compare.csv", "agreement_contingency.csv",
        "agreement_metrics.json", "marker_validation_marker_silver.csv",
        "marker_validation_scanvi.csv", "pseudotime_corr.json",
        "paga_graph_compare.json", "disagreement_cells.csv",
        "fig_composition.png", "fig_agreement.png", "fig_pseudotime.png",
    ]
    missing = []
    for q in cfg["annotation"]["query_datasets"]:
        for f in expected_per_query:
            p = os.path.join(out_dir, q, f)
            if not os.path.exists(p):
                missing.append(p)
    if not os.path.exists(os.path.join(out_dir, "summary.json")):
        missing.append(os.path.join(out_dir, "summary.json"))

    print("\n================ SMOKE TEST SUMMARY ================")
    print(json.dumps(summary["queries"], indent=2, default=str)[:2000])
    if missing:
        print("\nFAIL — missing outputs:")
        for m in missing:
            print("  -", m)
        return 1
    print("\nPASS — all expected outputs present.")
    print(f"Outputs in: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
