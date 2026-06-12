"""Add GSE178325 WOT/CellRank2 official-silver result dirs to the manifest."""

from __future__ import annotations

import argparse
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "benchmark" / "results" / "result_manifest.yaml"
ENTRY_KEY = "gse178325_marker_fm_transition_silver_hvg2000"
RESULT_DIRS = {
    "cellrank2": {
        "A": "benchmark/results/cellrank2/gse178325_marker_fm_silver_A_hvg2000_formal",
    },
    "wot": {
        "A": "benchmark/results/wot/gse178325_marker_fm_silver_A_hvg2000_formal",
    },
}


def update_manifest(path: Path, dry_run: bool = False) -> bool:
    data = yaml.safe_load(path.read_text(encoding="utf-8-sig")) or {}
    official = data.setdefault("official_silver", {})
    if ENTRY_KEY not in official:
        raise KeyError(f"Missing official_silver entry: {ENTRY_KEY}")

    result_dirs = official[ENTRY_KEY].setdefault("result_dirs", {})
    changed = False
    for method, scenarios in RESULT_DIRS.items():
        method_dirs = result_dirs.setdefault(method, {})
        for scenario, result_dir in scenarios.items():
            if method_dirs.get(scenario) != result_dir:
                method_dirs[scenario] = result_dir
                changed = True

    if changed and not dry_run:
        text = yaml.safe_dump(
            data,
            allow_unicode=True,
            sort_keys=False,
            width=120,
        )
        path.write_text(text, encoding="utf-8")
    return changed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--check", action="store_true", help="Only check whether changes would be needed.")
    args = parser.parse_args()

    changed = update_manifest(args.manifest, dry_run=args.check)
    if args.check:
        print("manifest_needs_update=" + ("1" if changed else "0"))
    elif changed:
        print(f"Updated {args.manifest}")
    else:
        print(f"No manifest changes needed: {args.manifest}")


if __name__ == "__main__":
    main()
