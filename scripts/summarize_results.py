import csv
import json
import os
import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import yaml


def find_project_root() -> Path:
    env_root = os.environ.get("TRAJ_PROJECT_ROOT")
    if env_root:
        root = Path(env_root).expanduser().resolve()
        if not root.exists():
            raise FileNotFoundError(f"TRAJ_PROJECT_ROOT does not exist: {root}")
        return root

    here = Path(__file__).resolve()
    for p in [here.parent] + list(here.parents):
        if (p / "benchmark").exists() and (p / "scripts").exists():
            return p

    raise RuntimeError(
        "Could not determine project root. "
        "Set TRAJ_PROJECT_ROOT or place this script inside the trajectory repo."
    )


def load_json(path: Path) -> Dict:
    raw = path.read_text(encoding="utf-8").rstrip("\x00\ufeff\n\r\t ")
    return json.loads(raw)


def guess_platform(result_dir: Path, run_meta: Dict) -> str:
    result_str = str(result_dir)
    config = str(run_meta.get("config", ""))

    if "shirokane_test" in result_str or "/home/xzy0723/" in config:
        return "shirokane"
    if ":\\" in config or config.startswith("C:\\"):
        return "windows"
    return "unknown"


def infer_result_class(result_dir: Path, run_meta: Dict) -> str:
    if run_meta.get("result_class"):
        return str(run_meta["result_class"])

    s = str(result_dir).lower()
    if "pilot" in s:
        return "pilot"
    if "smoke" in s or "reduced" in s or "cpu" in s:
        return "pilot"
    if "shirokane_test" in s:
        return "archive"
    if "backup" in s or "split" in s or "debug" in s:
        return "archive"
    return "official"


def infer_is_full_data(result_dir: Path, run_meta: Dict) -> str:
    s = str(result_dir).lower()
    if "pilot" in s:
        return "no"
    return "yes"


def get_baseline(metrics: Dict) -> Dict:
    return metrics.get("baseline", {}) if isinstance(metrics.get("baseline"), dict) else {}


def safe_float(x):
    return x if isinstance(x, (int, float)) else ""


def safe_str(x):
    return "" if x is None else str(x)


def collect_result_rows(method_root: Path) -> List[Dict]:
    rows: List[Dict] = []

    if not method_root.exists():
        return rows

    for result_dir in sorted(method_root.iterdir()):
        if not result_dir.is_dir():
            continue

        run_file = result_dir / "run_metadata.json"
        metric_file = result_dir / "lineage_metrics.json"

        if not run_file.exists() or not metric_file.exists():
            continue

        try:
            run_meta = load_json(run_file)
            metrics = load_json(metric_file)
        except Exception as e:
            print(f"skip unreadable: {result_dir} ({e})")
            continue

        baseline = get_baseline(metrics)

        method = safe_str(run_meta.get("method")) or method_root.name
        scenario = safe_str(run_meta.get("scenario"))

        if not scenario:
            # fallback from directory name
            name = result_dir.name.lower()
            if "scenario_a" in name:
                scenario = "A"
            elif "scenario_b" in name:
                scenario = "B"
            elif "scenario_c" in name:
                scenario = "C"

        platform_guess = guess_platform(result_dir, run_meta)
        result_class = infer_result_class(result_dir, run_meta)
        is_full_data = infer_is_full_data(result_dir, run_meta)

        row = {
            "method": method,
            "scenario": scenario,
            "result_dir_name": result_dir.name,
            "result_dir": str(result_dir),
            "platform_guess": platform_guess,
            "result_class": result_class,
            "is_full_data": is_full_data,

            "runtime_seconds": safe_float(run_meta.get("runtime_seconds")),
            "status": safe_str(run_meta.get("status")),
            "config": safe_str(run_meta.get("config")),
            "cell_state_key": safe_str(
                run_meta.get("cell_state_key") or metrics.get("cell_state_key")
            ),

            "auroc": safe_float(metrics.get("auroc")),
            "auprc": safe_float(metrics.get("auprc")),
            "jaccard_similarity": safe_float(metrics.get("jaccard_similarity")),
            "jaccard_similarity_topk": safe_float(metrics.get("jaccard_similarity_topk")),
            "single_step_recovery": safe_float(metrics.get("single_step_recovery")),
            "multi_step_recovery": safe_float(metrics.get("multi_step_recovery")),

            "baseline_auroc": safe_float(baseline.get("auroc")),
            "baseline_auprc": safe_float(baseline.get("auprc")),
            "baseline_jaccard_similarity": safe_float(baseline.get("jaccard_similarity")),
            "baseline_single_step_recovery": safe_float(baseline.get("single_step_recovery")),
            "baseline_multi_step_recovery": safe_float(baseline.get("multi_step_recovery")),

            "edge_confidence_mode": safe_str(metrics.get("edge_confidence_mode")),
            "n_reference_edges": safe_float(metrics.get("n_reference_edges")),
            "baseline_method": safe_str(baseline.get("method")),
            "baseline_style": safe_str(baseline.get("baseline_style")),
            "baseline_correlation_method": safe_str(baseline.get("correlation_method")),
            "baseline_averaging_method": safe_str(baseline.get("averaging_method")),
            "baseline_time_key": safe_str(baseline.get("time_key")),
            "baseline_n_states_used": safe_float(baseline.get("n_states_used")),
            "baseline_n_timepoint_pairs_used": safe_float(baseline.get("n_timepoint_pairs_used")),
            "baseline_n_source_cell_votes": safe_float(baseline.get("n_source_cell_votes")),
            "baseline_n_predicted_edges_thresholded": safe_float(
                baseline.get("n_predicted_edges_thresholded")
            ),
            "baseline_n_predicted_edges_topk": safe_float(baseline.get("n_predicted_edges_topk")),
        }
        rows.append(row)

    return rows


def collect_single_result_dir(result_dir: Path, method_hint: str = "", scenario_hint: str = "") -> List[Dict]:
    run_file = result_dir / "run_metadata.json"
    metric_file = result_dir / "lineage_metrics.json"

    if not run_file.exists() or not metric_file.exists():
        print(f"skip incomplete lineage result: {result_dir}")
        return []

    try:
        run_meta = load_json(run_file)
        metrics = load_json(metric_file)
    except Exception as e:
        print(f"skip unreadable: {result_dir} ({e})")
        return []

    baseline = get_baseline(metrics)
    method = safe_str(run_meta.get("method")) or method_hint or result_dir.parent.name
    scenario = safe_str(run_meta.get("scenario")) or scenario_hint
    platform_guess = guess_platform(result_dir, run_meta)
    result_class = infer_result_class(result_dir, run_meta)
    is_full_data = infer_is_full_data(result_dir, run_meta)

    return [{
        "method": method,
        "scenario": scenario,
        "result_dir_name": result_dir.name,
        "result_dir": str(result_dir),
        "platform_guess": platform_guess,
        "result_class": result_class,
        "dataset_role": safe_str(run_meta.get("dataset_role")),
        "benchmark_dataset": safe_str(run_meta.get("benchmark_dataset")),
        "is_full_data": is_full_data,

        "runtime_seconds": safe_float(run_meta.get("runtime_seconds")),
        "status": safe_str(run_meta.get("status")),
        "config": safe_str(run_meta.get("config")),
        "cell_state_key": safe_str(
            run_meta.get("cell_state_key") or metrics.get("cell_state_key")
        ),

        "auroc": safe_float(metrics.get("auroc")),
        "auprc": safe_float(metrics.get("auprc")),
        "jaccard_similarity": safe_float(metrics.get("jaccard_similarity")),
        "jaccard_similarity_topk": safe_float(metrics.get("jaccard_similarity_topk")),
        "single_step_recovery": safe_float(metrics.get("single_step_recovery")),
        "multi_step_recovery": safe_float(metrics.get("multi_step_recovery")),

        "baseline_auroc": safe_float(baseline.get("auroc")),
        "baseline_auprc": safe_float(baseline.get("auprc")),
        "baseline_jaccard_similarity": safe_float(baseline.get("jaccard_similarity")),
        "baseline_single_step_recovery": safe_float(baseline.get("single_step_recovery")),
        "baseline_multi_step_recovery": safe_float(baseline.get("multi_step_recovery")),

        "edge_confidence_mode": safe_str(metrics.get("edge_confidence_mode")),
        "n_reference_edges": safe_float(metrics.get("n_reference_edges")),
        "baseline_method": safe_str(baseline.get("method")),
        "baseline_style": safe_str(baseline.get("baseline_style")),
        "baseline_correlation_method": safe_str(baseline.get("correlation_method")),
        "baseline_averaging_method": safe_str(baseline.get("averaging_method")),
        "baseline_time_key": safe_str(baseline.get("time_key")),
        "baseline_n_states_used": safe_float(baseline.get("n_states_used")),
        "baseline_n_timepoint_pairs_used": safe_float(baseline.get("n_timepoint_pairs_used")),
        "baseline_n_source_cell_votes": safe_float(baseline.get("n_source_cell_votes")),
        "baseline_n_predicted_edges_thresholded": safe_float(
            baseline.get("n_predicted_edges_thresholded")
        ),
        "baseline_n_predicted_edges_topk": safe_float(baseline.get("n_predicted_edges_topk")),
    }]


def iter_manifest_result_dirs(manifest: Dict, groups: Iterable[str], root: Path) -> Iterable[Tuple[Path, str, str]]:
    for group in groups:
        group_items = manifest.get(group, {})
        if not isinstance(group_items, dict):
            continue
        for result_set in group_items.values():
            result_dirs = result_set.get("result_dirs", {})
            if not isinstance(result_dirs, dict):
                continue
            for method, scenarios in result_dirs.items():
                if isinstance(scenarios, dict):
                    for scenario, rel_path in scenarios.items():
                        yield (root / rel_path, str(method), str(scenario))
                elif isinstance(scenarios, str):
                    yield (root / scenarios, str(method), "")


def collect_manifest_rows(manifest_path: Path, groups: List[str], root: Path) -> List[Dict]:
    with open(manifest_path, "r", encoding="utf-8") as f:
        manifest = yaml.safe_load(f) or {}

    rows: List[Dict] = []
    for result_dir, method, scenario in iter_manifest_result_dirs(manifest, groups, root):
        rows.extend(collect_single_result_dir(result_dir, method, scenario))
    return rows


def preferred_rank(row: Dict) -> Tuple[int, int, float]:
    """
    Lower is better.
    Rank priority:
      1) official full-data shirokane
      2) shirokane_test full-data
      3) official full-data windows/local
      4) anything full-data unknown
      5) pilot
    Then prefer completed status, then shorter runtime as tie-breaker.
    """
    result_class = row.get("result_class", "")
    platform = row.get("platform_guess", "")
    is_full_data = row.get("is_full_data", "")
    status = row.get("status", "")

    if result_class == "pilot":
        base = 50
    elif is_full_data == "yes" and platform == "shirokane" and result_class == "official":
        base = 10
    elif is_full_data == "yes" and platform == "shirokane":
        base = 11
    elif is_full_data == "yes" and platform == "windows":
        base = 20
    elif is_full_data == "yes":
        base = 30
    else:
        base = 40

    status_penalty = 0 if status == "completed" else 1
    runtime = row.get("runtime_seconds")
    runtime_val = float(runtime) if runtime not in ("", None) else 1e18

    return (base, status_penalty, runtime_val)


def write_csv(path: Path, rows: List[Dict]) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")

    fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def select_preferred_rows(rows: List[Dict]) -> List[Dict]:
    grouped: Dict[Tuple[str, str, str], List[Dict]] = {}
    for row in rows:
        key = (
            row.get("benchmark_dataset", ""),
            row["method"],
            row["scenario"],
        )
        grouped.setdefault(key, []).append(row)

    preferred: List[Dict] = []
    for key, candidates in sorted(grouped.items()):
        chosen = sorted(candidates, key=preferred_rank)[0]
        preferred.append(chosen)

    return preferred


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize lineage benchmark results.")
    parser.add_argument(
        "--manifest",
        type=Path,
        help="Optional result manifest. When provided, only listed result_dirs from the selected groups are summarized.",
    )
    parser.add_argument(
        "--groups",
        nargs="+",
        default=["official"],
        help="Manifest groups to summarize, for example: official external_validation.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for summary CSVs. Defaults to benchmark/results.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    root = find_project_root()
    results_root = root / "benchmark" / "results"
    output_dir = (root / args.output_dir) if args.output_dir and not args.output_dir.is_absolute() else args.output_dir
    if output_dir is None:
        output_dir = results_root
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.manifest:
        manifest_path = (root / args.manifest) if not args.manifest.is_absolute() else args.manifest
        all_rows = collect_manifest_rows(manifest_path, args.groups, root)
    else:
        all_rows: List[Dict] = []
        all_rows.extend(collect_result_rows(results_root / "wot"))
        all_rows.extend(collect_result_rows(results_root / "cellrank2"))

    if not all_rows:
        raise RuntimeError("No result directories with both run_metadata.json and lineage_metrics.json were found.")

    # Stable sort for readability
    all_rows = sorted(
        all_rows,
        key=lambda r: (
            r["method"],
            r["scenario"],
            r["platform_guess"],
            r["result_class"],
            r["result_dir_name"],
        ),
    )

    preferred_rows = select_preferred_rows(all_rows)

    out_all = output_dir / "summary_lineage_metrics_all.csv"
    out_preferred = output_dir / "summary_lineage_metrics_preferred.csv"

    write_csv(out_all, all_rows)
    write_csv(out_preferred, preferred_rows)

    print(f"Project root: {root}")
    if args.manifest:
        print(f"Manifest    : {manifest_path}")
        print(f"Groups      : {', '.join(args.groups)}")
    print(f"Wrote all runs     : {out_all}")
    print(f"Wrote preferred set: {out_preferred}")
    print(f"Total runs found   : {len(all_rows)}")
    print(f"Preferred rows     : {len(preferred_rows)}")


if __name__ == "__main__":
    main()
