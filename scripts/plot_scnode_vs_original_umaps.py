#!/usr/bin/env python3
"""
Plot original-vs-scNODE cell-label UMAPs for every dataset with a scNODE
Scenario-A runtime config.

The UMAP for each dataset is fit on a combined sample of:
  - original benchmark-input HVG expression
  - scNODE projected_expression

The original and scNODE points are then saved as two separate figures that share
the same UMAP coordinate system and color map. This makes distribution shifts
between observed and inferred cells visually comparable.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from sklearn.decomposition import PCA


@dataclass(frozen=True)
class DatasetSpec:
    dataset_id: str
    run_id: str
    config_path: Path
    h5ad_path: Path
    scnode_dir: Path
    label_key: str


def find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        return Path(env).resolve()
    here = Path(__file__).resolve()
    for candidate in [here.parent, *here.parents]:
        if (candidate / "benchmark").is_dir() and (candidate / "scripts").is_dir():
            return candidate
    return here.parents[1]


def discover_specs(root: Path) -> list[DatasetSpec]:
    specs: list[DatasetSpec] = []
    cfg_paths = []
    for cfg_dir in [
        root / "benchmark" / "configs" / "runtime",
        root / "benchmark" / "configs",
    ]:
        cfg_paths.extend(sorted(cfg_dir.glob("scnode_*_A_hvg2000_formal.yaml")))

    seen_paths: set[Path] = set()
    for cfg_path in cfg_paths:
        if cfg_path in seen_paths:
            continue
        seen_paths.add(cfg_path)
        cfg = yaml.safe_load(cfg_path.read_text(encoding="utf-8-sig")) or {}
        if cfg.get("method") != "scnode" or str(cfg.get("scenario")) != "A":
            continue

        dataset = cfg.get("dataset") or {}
        output = cfg.get("output") or {}
        h5ad_rel = dataset.get("h5ad_path")
        out_rel = output.get("base_dir")
        dataset_id = str(dataset.get("id") or cfg_path.stem)
        if not h5ad_rel or not out_rel:
            continue

        scnode_dir = root / out_rel
        if not (scnode_dir / "projected_expression.npy").exists():
            continue

        projected_milestone = scnode_dir / "projected_milestone_labels.csv"
        label_key = (
            "final_milestone_label_expanded"
            if projected_milestone.exists()
            else "final_milestone_label_coarse"
        )

        specs.append(
            DatasetSpec(
                dataset_id=dataset_id,
                run_id=str(cfg.get("run_id") or cfg_path.stem),
                config_path=cfg_path,
                h5ad_path=root / h5ad_rel,
                scnode_dir=scnode_dir,
                label_key=label_key,
            )
        )

    # Keep only the four formal dataset runs, not smoke/debug configs.
    specs = [
        spec for spec in specs
        if re.match(r"^GSE\d+$", spec.dataset_id)
    ]
    return specs


def stratified_positions(labels: np.ndarray, max_total: int, seed: int) -> np.ndarray:
    labels = np.asarray(labels).astype(str)
    n = labels.shape[0]
    if max_total <= 0 or n <= max_total:
        return np.arange(n, dtype=int)

    rng = np.random.default_rng(seed)
    series = pd.Series(labels)
    groups = {label: idx.to_numpy(dtype=int) for label, idx in series.groupby(series).groups.items()}
    keys = sorted(groups)
    counts = np.array([len(groups[k]) for k in keys], dtype=float)

    if len(keys) > max_total:
        chosen_labels = sorted(keys, key=lambda k: len(groups[k]), reverse=True)[:max_total]
        return np.array([rng.choice(groups[k]) for k in chosen_labels], dtype=int)

    raw = counts / counts.sum() * max_total
    alloc = np.maximum(1, np.floor(raw).astype(int))
    alloc = np.minimum(alloc, counts.astype(int))

    while alloc.sum() < max_total:
        capacity = counts.astype(int) - alloc
        if capacity.max() <= 0:
            break
        residual = raw - np.floor(raw)
        scores = np.where(capacity > 0, residual, -1.0)
        alloc[int(np.argmax(scores))] += 1

    while alloc.sum() > max_total:
        candidates = np.where(alloc > 1)[0]
        if len(candidates) == 0:
            break
        residual = raw[candidates] - np.floor(raw[candidates])
        alloc[int(candidates[np.argmin(residual)])] -= 1

    selected: list[np.ndarray] = []
    for key, n_take in zip(keys, alloc):
        selected.append(rng.choice(groups[key], size=int(n_take), replace=False))
    return np.sort(np.concatenate(selected).astype(int))


def to_dense_float32(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        matrix = matrix.toarray()
    arr = np.asarray(matrix, dtype=np.float32)
    np.nan_to_num(arr, copy=False, nan=0.0, posinf=0.0, neginf=0.0)
    return arr


def load_original_sample(spec: DatasetSpec, max_cells: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    adata = ad.read_h5ad(spec.h5ad_path, backed="r")
    try:
        label_key = spec.label_key
        if label_key not in adata.obs.columns:
            label_key = "final_milestone_label_coarse"
        labels_all = adata.obs[label_key].astype(str).to_numpy()
        pos = stratified_positions(labels_all, max_total=max_cells, seed=seed)
        labels = labels_all[pos]
        X = to_dense_float32(adata[pos, :].X)
        return X, labels
    finally:
        try:
            adata.file.close()
        except Exception:
            pass


def load_scnode_sample(spec: DatasetSpec, max_cells: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    expr = np.load(spec.scnode_dir / "projected_expression.npy", mmap_mode="r")
    X_all = expr.reshape((-1, expr.shape[-1]))

    milestone_path = spec.scnode_dir / "projected_milestone_labels.csv"
    cluster_path = spec.scnode_dir / "projected_cluster_labels.csv"
    if milestone_path.exists():
        labels_df = pd.read_csv(milestone_path)
        if spec.label_key in labels_df.columns:
            labels_all = labels_df[spec.label_key].astype(str).to_numpy()
        elif "final_milestone_label_coarse" in labels_df.columns:
            labels_all = labels_df["final_milestone_label_coarse"].astype(str).to_numpy()
        else:
            labels_all = labels_df["projected_cluster_label"].astype(str).to_numpy()
    else:
        labels_df = pd.read_csv(cluster_path)
        labels_all = labels_df["projected_cluster_label"].astype(str).to_numpy()

    if len(labels_all) != X_all.shape[0]:
        raise ValueError(
            f"{spec.run_id}: label count {len(labels_all)} does not match "
            f"projected_expression rows {X_all.shape[0]}"
        )

    pos = stratified_positions(labels_all, max_total=max_cells, seed=seed + 17)
    labels = labels_all[pos]
    X = to_dense_float32(X_all[pos, :])
    return X, labels


def color_map(labels: np.ndarray) -> dict[str, tuple]:
    unique = sorted(pd.Series(labels.astype(str)).unique())
    palettes = [
        plt.get_cmap("tab20"),
        plt.get_cmap("tab20b"),
        plt.get_cmap("tab20c"),
        plt.get_cmap("Set3"),
    ]
    colors: list[tuple] = []
    for cmap in palettes:
        colors.extend([cmap(i) for i in range(cmap.N)])
    return {label: colors[i % len(colors)] for i, label in enumerate(unique)}


def plot_panel(
    coords: np.ndarray,
    labels: np.ndarray,
    title: str,
    out_path: Path,
    label_to_color: dict[str, tuple],
    point_size: float,
) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 7.2), dpi=180)
    labels = labels.astype(str)

    for label in sorted(label_to_color):
        mask = labels == label
        if not np.any(mask):
            continue
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=point_size,
            c=[label_to_color[label]],
            label=label,
            linewidths=0,
            alpha=0.72,
        )

    for label in sorted(label_to_color):
        mask = labels == label
        if not np.any(mask):
            continue
        center = np.median(coords[mask], axis=0)
        ax.text(
            center[0],
            center[1],
            label,
            fontsize=8,
            ha="center",
            va="center",
            color="black",
            bbox={
                "boxstyle": "round,pad=0.18",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.72,
            },
        )

    ax.set_title(title, fontsize=12)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        markerscale=3.0,
        fontsize=7,
        title="cell label",
        title_fontsize=8,
    )
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def fit_umap(X_orig: np.ndarray, X_scn: np.ndarray, seed: int) -> tuple[np.ndarray, np.ndarray]:
    import umap

    X = np.vstack([X_orig, X_scn]).astype(np.float32, copy=False)
    n_pcs = min(50, X.shape[1], X.shape[0] - 1)
    pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=seed)
    X_pca = pca.fit_transform(X)
    reducer = umap.UMAP(
        n_neighbors=30,
        min_dist=0.25,
        metric="euclidean",
        random_state=seed,
        n_components=2,
    )
    coords = reducer.fit_transform(X_pca)
    return coords[: len(X_orig)], coords[len(X_orig):]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-original-cells", type=int, default=20000)
    parser.add_argument("--max-scnode-cells", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("benchmark/reports/visualizations/scnode_vs_original"),
    )
    args = parser.parse_args()

    root = find_project_root()
    out_dir = args.output_dir if args.output_dir.is_absolute() else root / args.output_dir
    specs = discover_specs(root)
    if not specs:
        raise RuntimeError("No scNODE Scenario-A dataset specs found.")

    manifest = []
    for i, spec in enumerate(specs):
        print("=" * 80)
        print(f"[{i + 1}/{len(specs)}] {spec.dataset_id} ({spec.run_id})")
        print(f"  h5ad      : {spec.h5ad_path}")
        print(f"  scNODE dir: {spec.scnode_dir}")
        print(f"  label key : {spec.label_key}")

        X_orig, y_orig = load_original_sample(
            spec,
            max_cells=args.max_original_cells,
            seed=args.seed + i * 100,
        )
        X_scn, y_scn = load_scnode_sample(
            spec,
            max_cells=args.max_scnode_cells,
            seed=args.seed + i * 100,
        )
        print(f"  original sample: {X_orig.shape}")
        print(f"  scNODE sample  : {X_scn.shape}")

        emb_orig, emb_scn = fit_umap(X_orig, X_scn, seed=args.seed + i)
        union_labels = np.concatenate([y_orig.astype(str), y_scn.astype(str)])
        colors = color_map(union_labels)

        prefix = out_dir / spec.run_id
        original_png = prefix.with_name(prefix.name + "_original_cell_labels_umap.png")
        scnode_png = prefix.with_name(prefix.name + "_scnode_cell_labels_umap.png")

        plot_panel(
            emb_orig,
            y_orig,
            f"{spec.dataset_id} original benchmark input",
            original_png,
            colors,
            point_size=3.5,
        )
        plot_panel(
            emb_scn,
            y_scn,
            f"{spec.dataset_id} scNODE inferred cells",
            scnode_png,
            colors,
            point_size=4.0,
        )

        manifest.append(
            {
                "dataset_id": spec.dataset_id,
                "run_id": spec.run_id,
                "config_path": str(spec.config_path.relative_to(root)),
                "h5ad_path": str(spec.h5ad_path.relative_to(root)),
                "scnode_dir": str(spec.scnode_dir.relative_to(root)),
                "label_key": spec.label_key,
                "n_original_plotted": int(len(y_orig)),
                "n_scnode_plotted": int(len(y_scn)),
                "original_png": str(original_png.relative_to(root)),
                "scnode_png": str(scnode_png.relative_to(root)),
                "labels": sorted(pd.Series(union_labels).unique().tolist()),
            }
        )
        print(f"  wrote: {original_png}")
        print(f"  wrote: {scnode_png}")

    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("=" * 80)
    print(f"Wrote manifest: {manifest_path}")


if __name__ == "__main__":
    main()
