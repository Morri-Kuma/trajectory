#!/usr/bin/env python3
"""
Lightweight sanity-check / validation script for GSE175634 benchmark inputs.

Checks that expression data, metadata, cell annotations, and the ground-truth
provider are present and internally consistent before submitting jobs to
Shirokane.

# ─────────────────────────────────────────────────────────────────────────────
# DATASET ID NOTE
# ─────────────────────────────────────────────────────────────────────────────
# The canonical GEO accession is GSE175634 (not GSE174534 — digits transposed
# in the original request). All source files use the prefix GSE175634_.
# ─────────────────────────────────────────────────────────────────────────────

Run from the repository root:
    python scripts/validate_gse175634_inputs.py

Returns exit code 0 if all checks pass, 1 if any check fails.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

# ─── Config ───────────────────────────────────────────────────────────────────

DATASET_ID    = "GSE175634"
PROVIDER_ID   = "gse175634_cardiac_silver_v1"
STATE_KEY     = "final_milestone_label_coarse"
TIME_KEY      = "abs_day"
EXPECTED_STATES = {"IPSC", "MES", "CMES", "PROG", "CM", "CF"}

RAW_DATA_FILES = [
    "GSE175634_cell_counts.mtx",
    "gene_indices_counts.tsv",
    "GSE175634_cell_indices.tsv",
    "GSE175634_cell_metadata.tsv",
]

BENCHMARK_H5AD = (
    "benchmark/inputs/gse175634_cardiac_author_hvg2000/"
    "GSE175634_cardiac_author_HVG2000_benchmark_input.h5ad"
)

PROVIDER_FILES = [
    f"benchmark/ground_truth/providers/{PROVIDER_ID}/state_labels.tsv",
    f"benchmark/ground_truth/providers/{PROVIDER_ID}/state_metadata.tsv",
    f"benchmark/ground_truth/providers/{PROVIDER_ID}/reference_graph.json",
    f"benchmark/ground_truth/providers/{PROVIDER_ID}/reference_graph_edges.csv",
    f"benchmark/ground_truth/providers/{PROVIDER_ID}/ground_truth_metadata.json",
]

RUNTIME_CONFIGS = [
    f"benchmark/configs/runtime/{method}_gse175634_cardiac_silver_{sc}_hvg2000_formal.yaml"
    for method in ("scnode", "prescient", "mioflow")
    for sc in ("A", "B", "C")
] + [
    "benchmark/configs/runtime/wot_gse175634_cardiac_silver_A_hvg2000_formal.yaml",
    "benchmark/configs/runtime/cellrank2_gse175634_cardiac_silver_A_hvg2000_formal.yaml",
]


# ─── Helpers ─────────────────────────────────────────────────────────────────

PASS = "\033[32m[PASS]\033[0m"
FAIL = "\033[31m[FAIL]\033[0m"
WARN = "\033[33m[WARN]\033[0m"
INFO = "\033[34m[INFO]\033[0m"

_failures: list[str] = []


def ok(msg: str) -> None:
    print(f"  {PASS} {msg}")


def fail(msg: str) -> None:
    print(f"  {FAIL} {msg}")
    _failures.append(msg)


def warn(msg: str) -> None:
    print(f"  {WARN} {msg}")


def info(msg: str) -> None:
    print(f"  {INFO} {msg}")


def find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
        raise FileNotFoundError(f"TRAJ_PROJECT_ROOT={env!r} does not exist.")
    here = Path(__file__).resolve().parent
    for candidate in [here, *here.parents]:
        if (candidate / "benchmark").is_dir() and (candidate / "scripts").is_dir():
            return candidate
    return here.parent


# ─── Checks ──────────────────────────────────────────────────────────────────

def check_raw_data(root: Path) -> None:
    print("\n[1] Raw data files (data/gse175634/)")
    data_dir = root / "data" / "gse175634"
    for fname in RAW_DATA_FILES:
        p = data_dir / fname
        if p.exists():
            size_mb = p.stat().st_size / 1e6
            ok(f"{fname}  ({size_mb:.0f} MB)")
        else:
            fail(f"Missing: {p}")


def check_benchmark_h5ad(root: Path) -> None:
    print(f"\n[2] Benchmark h5ad  ({BENCHMARK_H5AD})")
    h5ad = root / BENCHMARK_H5AD
    if not h5ad.exists():
        fail(f"Benchmark h5ad not found: {h5ad}")
        warn("Run: python scripts/build_gse175634_cardiac_author_input.py")
        return
    size_mb = h5ad.stat().st_size / 1e6
    ok(f"h5ad exists  ({size_mb:.0f} MB)")

    # Deep check with anndata
    try:
        import anndata as ad
    except ImportError:
        warn("anndata not installed — skipping h5ad content checks.")
        return

    try:
        adata = ad.read_h5ad(h5ad, backed="r")
    except Exception as exc:
        fail(f"Could not read h5ad: {exc}")
        return

    ok(f"n_obs={adata.n_obs:,}, n_vars={adata.n_vars:,}")

    # Required obs columns
    required_obs = ["cell_id", "sample_id", "time_label", "abs_day",
                    STATE_KEY, "diffday", "individual", "type"]
    missing_obs = [c for c in required_obs if c not in adata.obs.columns]
    if missing_obs:
        fail(f"Missing obs columns: {missing_obs}")
    else:
        ok(f"All required obs columns present: {required_obs}")

    # layers['counts']
    if "counts" not in adata.layers:
        fail("layers['counts'] missing")
    else:
        ok("layers['counts'] present")

    # State labels
    if STATE_KEY in adata.obs.columns:
        label_counts = adata.obs[STATE_KEY].value_counts().to_dict()
        info(f"State label counts: {label_counts}")
        found_states = set(label_counts.keys()) - {"UNK"}
        missing_states = EXPECTED_STATES - found_states
        if missing_states:
            warn(f"Expected states not found after QC: {missing_states}")
        else:
            ok(f"All expected states present: {sorted(EXPECTED_STATES)}")
        n_unk = label_counts.get("UNK", 0)
        info(f"UNK cells: {n_unk:,} (excluded from official metrics)")

    # Time axis
    if TIME_KEY in adata.obs.columns:
        days = sorted(adata.obs[TIME_KEY].dropna().unique().tolist())
        ok(f"abs_day values: {days}")
        if len(days) < 5:
            warn(f"Fewer time points than expected ({len(days)}); check QC filtering.")

    # dataset_id in uns
    ds = adata.uns.get("dataset_id", "")
    if ds != DATASET_ID:
        warn(f"uns['dataset_id']={ds!r}, expected {DATASET_ID!r}")
    else:
        ok(f"uns['dataset_id'] = {DATASET_ID!r}")

    # silver_standard_skipped flag
    skipped = adata.uns.get("silver_standard_skipped", None)
    if not skipped:
        warn("uns['silver_standard_skipped'] not set or False")
    else:
        ok("uns['silver_standard_skipped'] is set (silver step correctly skipped)")

    try:
        adata.file.close()
    except Exception:
        pass


def check_provider(root: Path) -> None:
    print(f"\n[3] Ground-truth provider  ({PROVIDER_ID})")
    for rel in PROVIDER_FILES:
        p = root / rel
        if p.exists():
            ok(str(rel))
        else:
            fail(f"Missing: {rel}")
            if "state_labels" in rel or "reference_graph" in rel:
                warn(
                    "Run: python scripts/build_gse175634_cardiac_author_provider.py"
                )

    # Reference graph structural check
    graph_path = root / f"benchmark/ground_truth/providers/{PROVIDER_ID}/reference_graph.json"
    if graph_path.exists():
        try:
            g = json.loads(graph_path.read_text(encoding="utf-8"))
            n_nodes = len(g.get("nodes", []))
            n_edges = len(g.get("edges", []))
            ok(f"reference_graph.json: {n_nodes} nodes, {n_edges} edges")
            if n_nodes != 6:
                warn(f"Expected 6 nodes (IPSC MES CMES PROG CM CF), got {n_nodes}")
            if n_edges != 5:
                warn(f"Expected 5 edges, got {n_edges}")
            node_ids = {nd["id"] for nd in g.get("nodes", [])}
            missing = EXPECTED_STATES - node_ids
            if missing:
                fail(f"Missing expected nodes in graph: {missing}")
        except Exception as exc:
            fail(f"Could not parse reference_graph.json: {exc}")

    # Registry entry
    registry_path = root / "benchmark/ground_truth/registry.yaml"
    if registry_path.exists():
        text = registry_path.read_text(encoding="utf-8")
        if PROVIDER_ID in text:
            ok(f"Provider '{PROVIDER_ID}' found in registry.yaml")
        else:
            fail(f"Provider '{PROVIDER_ID}' NOT found in registry.yaml — add it.")


def check_runtime_configs(root: Path) -> None:
    print(f"\n[4] Runtime benchmark configs  (benchmark/configs/runtime/)")
    try:
        import yaml
    except ImportError:
        yaml = None
        warn("PyYAML not installed; checking config presence only.")

    for rel in RUNTIME_CONFIGS:
        p = root / rel
        if p.exists():
            ok(Path(rel).name)
            if yaml is None:
                continue
            try:
                cfg = yaml.safe_load(p.read_text(encoding="utf-8-sig")) or {}
            except Exception as exc:
                fail(f"Could not parse YAML config {rel}: {exc}")
                continue

            dataset_id = (cfg.get("dataset") or {}).get("id")
            if dataset_id != DATASET_ID:
                fail(f"{rel}: dataset.id={dataset_id!r}, expected {DATASET_ID!r}")

            provider_id = (
                cfg.get("provider_id")
                or (cfg.get("ground_truth") or {}).get("provider_id")
            )
            if provider_id != PROVIDER_ID:
                fail(f"{rel}: provider_id={provider_id!r}, expected {PROVIDER_ID!r}")

            scenario = cfg.get("scenario")
            if scenario in {"B", "C"}:
                scenario_params = cfg.get("scenario_params") or {}
                if not scenario_params.get("train_times"):
                    fail(f"{rel}: scenario {scenario} missing train_times")
                if not scenario_params.get("heldout_times"):
                    fail(f"{rel}: scenario {scenario} missing heldout_times")
        else:
            fail(f"Missing: {rel}")


def check_shirokane_scripts(root: Path) -> None:
    print(f"\n[5] Shirokane submission scripts")
    scripts = [
        "jobs/shirokane/run_gse175634_preprocess.sh",
        "jobs/shirokane/run_gse175634_benchmark_array.sh",
        "jobs/shirokane/submit_gse175634_benchmark.sh",
    ]
    for s in scripts:
        p = root / s
        if p.exists():
            ok(s)
        else:
            fail(f"Missing: {s}")


# ─── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    root = find_project_root()
    print("=" * 64)
    print(f"GSE175634 input validation")
    print(f"  Dataset ID   : {DATASET_ID}")
    print(f"  Provider ID  : {PROVIDER_ID}")
    print(f"  Project root : {root}")
    print("=" * 64)

    check_raw_data(root)
    check_benchmark_h5ad(root)
    check_provider(root)
    check_runtime_configs(root)
    check_shirokane_scripts(root)

    print()
    print("=" * 64)
    if _failures:
        print(f"RESULT: {len(_failures)} check(s) FAILED:")
        for f in _failures:
            print(f"  - {f}")
        print("=" * 64)
        sys.exit(1)
    else:
        print("RESULT: All checks PASSED.")
        print("=" * 64)
        sys.exit(0)


if __name__ == "__main__":
    main()
