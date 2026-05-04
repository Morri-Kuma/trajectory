"""
Build a scGPT-derived ground-truth provider for GSE178325 0618-only data.

This script follows the project framework used for the GSE230659 scGPT-v1
provider, but keeps GSE178325 self-contained:

1. Load the GSE178325 0618 post-QC raw/full-gene benchmark input.
2. Compute scGPT embeddings into obsm["X_scGPT"] unless already present.
3. Run neighbors, UMAP, and Leiden clustering on X_scGPT.
4. Treat Leiden clusters as scGPT pseudostates.
5. Build a frozen reference graph by consecutive-timepoint kNN transition
   counting in scGPT embedding space.
6. Export provider files under a dataset-specific provider directory, while
   keeping the same artifact names and state-key convention as the GSE230659
   scGPT-v1 workflow.

The generated provider is a silver-standard scGPT-derived benchmark asset, not
final curated biological ground truth.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import scipy.sparse as sp
from scipy.stats import entropy as scipy_entropy
from torch.utils.data import DataLoader, SequentialSampler
from tqdm import tqdm


PROVIDER_ID = "scgpt_v1_gse178325_0618"
VERSION = "v1"


def _find_project_root() -> Path:
    env = os.environ.get("TRAJ_PROJECT_ROOT")
    if env:
        p = Path(env)
        if p.exists():
            return p
        raise FileNotFoundError(f"TRAJ_PROJECT_ROOT={env!r} does not exist")
    here = Path(__file__).resolve()
    for candidate in [here.parent, *here.parents]:
        if (candidate / "benchmark").exists() and (candidate / "data").exists():
            return candidate
    return here.parents[2]


def _parse_args() -> argparse.Namespace:
    root = _find_project_root()
    parser = argparse.ArgumentParser(
        description="Build GSE178325 0618 scGPT-derived ground-truth provider"
    )
    parser.add_argument(
        "--input-h5ad",
        default=str(
            root
            / "data/processed/gse178325_human/"
            / "GSE178325_0618_raw_full_gene_benchmark_input.h5ad"
        ),
    )
    parser.add_argument(
        "--model-dir",
        default=str(root / "models/scgpt_whole_human"),
    )
    parser.add_argument(
        "--output-dir",
        default=str(root / "benchmark/results/scgpt/gse178325_0618/full"),
    )
    parser.add_argument(
        "--provider-dir",
        default=str(root / "benchmark/ground_truth/providers" / PROVIDER_ID),
    )
    parser.add_argument("--time-key", default="abs_day")
    parser.add_argument("--stage-key", default="stage")
    parser.add_argument("--stage-proxy-key", default="stage")
    parser.add_argument("--gene-col", default="index")
    parser.add_argument("--leiden-resolution", type=float, default=0.5)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument("--max-length", type=int, default=1200)
    parser.add_argument("--knn-k", type=int, default=10)
    parser.add_argument("--min-edge-weight", type=float, default=0.02)
    parser.add_argument("--chunk-size", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--use-fast-transformer", action="store_true")
    parser.add_argument(
        "--require-cuda",
        action="store_true",
        help="Fail immediately if CUDA is unavailable.",
    )
    parser.add_argument(
        "--skip-embedding-if-exists",
        action="store_true",
        help="Reuse output adata if it already contains X_scGPT.",
    )
    parser.add_argument(
        "--skip-umap",
        action="store_true",
        help="Skip UMAP plotting if only provider files are needed.",
    )
    return parser.parse_args()


def _add_scgpt_to_path() -> None:
    scgpt_env = os.environ.get("SCGPT_REPO")
    if not scgpt_env:
        return
    scgpt_repo = Path(scgpt_env)
    if scgpt_repo.exists() and str(scgpt_repo) not in sys.path:
        sys.path.insert(0, str(scgpt_repo))


def _load_or_embed(args: argparse.Namespace, out_dir: Path):
    output_h5ad = out_dir / "adata_scgpt_full.h5ad"
    leiden_key = f"leiden_scgpt_res{args.leiden_resolution}"

    if args.skip_embedding_if_exists and output_h5ad.exists():
        adata = sc.read_h5ad(output_h5ad)
        if "X_scGPT" in adata.obsm and leiden_key in adata.obs:
            print(f"[provider] Reusing existing artifact: {output_h5ad}")
            return adata, output_h5ad, {"reused_existing_artifact": True}
        print("[provider] Existing artifact lacks X_scGPT or Leiden; recomputing.")

    input_h5ad = Path(args.input_h5ad)
    model_dir = Path(args.model_dir)
    if not input_h5ad.exists():
        raise FileNotFoundError(input_h5ad)
    for fname in ["vocab.json", "args.json", "best_model.pt"]:
        if not (model_dir / fname).exists():
            raise FileNotFoundError(model_dir / fname)

    print(f"[provider] Loading input: {input_h5ad}")
    adata = sc.read_h5ad(input_h5ad)
    print(f"[provider] Loaded shape: {adata.n_obs:,} x {adata.n_vars:,}")
    if not adata.var_names.is_unique:
        raise ValueError(
            "Input var_names are not unique. Rebuild the raw/full-gene input "
            "with scripts/build_gse178325_raw_full_gene_input.py so the "
            "GSE178325 scGPT path matches the GSE230659 full-gene workflow."
        )

    if args.time_key not in adata.obs:
        raise KeyError(f"time key not found in obs: {args.time_key}")
    if args.stage_proxy_key not in adata.obs:
        raise KeyError(f"stage proxy key not found in obs: {args.stage_proxy_key}")

    with open(model_dir / "vocab.json", encoding="utf-8") as f:
        vocab = json.load(f)
    vocab_genes = set(vocab) - {"<pad>", "<cls>", "<eoc>"}
    gene_symbols = set(map(str, adata.var.index))
    overlap = gene_symbols & vocab_genes
    overlap_frac = len(overlap) / max(len(gene_symbols), 1)
    print(
        "[provider] Vocab overlap: "
        f"{len(overlap):,}/{len(gene_symbols):,} ({100 * overlap_frac:.1f}%)"
    )
    if overlap_frac < 0.3:
        raise RuntimeError("scGPT vocabulary overlap is too low; check var.index.")

    if torch.cuda.is_available():
        device = "cuda"
        gpu_name = torch.cuda.get_device_name(0)
    else:
        if args.require_cuda:
            raise RuntimeError(
                "CUDA is unavailable. This job must run on a GPU node; "
                "resubmit to the H100 queue instead of CPU queues."
            )
        device = "cpu"
        gpu_name = "CPU"
        print("[provider] WARNING: CUDA unavailable; scGPT embedding may be slow.")

    t0 = time.time()
    adata = _embed_data_low_memory(
        adata,
        model_dir,
        gene_col=args.gene_col,
        max_length=args.max_length,
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        device=device,
        use_fast_transformer=args.use_fast_transformer,
    )
    embed_runtime = time.time() - t0
    if "X_scGPT" not in adata.obsm:
        raise RuntimeError("scGPT embedding did not populate obsm['X_scGPT']")

    emb = np.asarray(adata.obsm["X_scGPT"], dtype=np.float32)
    norms = np.linalg.norm(emb, axis=1)
    print(
        "[provider] X_scGPT:",
        emb.shape,
        "norm mean/std:",
        f"{norms.mean():.6f}/{norms.std():.6f}",
    )

    sc.settings.verbosity = 1
    sc.pp.neighbors(
        adata,
        use_rep="X_scGPT",
        n_neighbors=args.n_neighbors,
        random_state=args.seed,
    )
    if not args.skip_umap:
        sc.tl.umap(adata, random_state=args.seed)
    sc.tl.leiden(
        adata,
        resolution=args.leiden_resolution,
        random_state=args.seed,
        key_added=leiden_key,
    )
    print(
        f"[provider] Leiden {leiden_key}: "
        f"{adata.obs[leiden_key].nunique()} clusters"
    )

    adata.write_h5ad(output_h5ad)
    meta = {
        "reused_existing_artifact": False,
        "device": device,
        "gpu_name": gpu_name,
        "embed_runtime_sec": round(embed_runtime, 1),
        "vocab_overlap_frac": round(overlap_frac, 4),
        "n_vocab_overlap_genes": int(len(overlap)),
    }
    return adata, output_h5ad, meta


def _embed_data_low_memory(
    adata,
    model_dir: Path,
    gene_col: str,
    max_length: int,
    batch_size: int,
    num_workers: int,
    device: str,
    use_fast_transformer: bool,
):
    """
    Low-memory scGPT embedding path.

    scGPT's public ``embed_data`` helper densifies ``adata.X`` and may spawn many
    DataLoader workers. For GSE178325 this exceeded 128G virtual memory before
    the first embedding batch. This implementation keeps ``adata.X`` sparse,
    materializes only one row per Dataset item, and uses ``num_workers=0`` by
    default.
    """
    _add_scgpt_to_path()
    from scgpt.data_collator import DataCollator  # noqa: WPS433
    from scgpt.model import TransformerModel  # noqa: WPS433
    from scgpt.tokenizer import GeneVocab  # noqa: WPS433
    from scgpt.utils import load_pretrained  # noqa: WPS433

    model_dir = Path(model_dir)
    vocab_file = model_dir / "vocab.json"
    model_config_file = model_dir / "args.json"
    model_file = model_dir / "best_model.pt"
    pad_token = "<pad>"
    special_tokens = [pad_token, "<cls>", "<eoc>"]

    if gene_col == "index":
        adata.var["index"] = adata.var.index.astype(str)
    elif gene_col not in adata.var:
        raise KeyError(f"gene_col not found in adata.var: {gene_col}")

    vocab = GeneVocab.from_file(vocab_file)
    for token in special_tokens:
        if token not in vocab:
            vocab.append_token(token)
    vocab.set_default_index(vocab[pad_token])

    id_in_vocab = np.array(
        [vocab[gene] if gene in vocab else -1 for gene in adata.var[gene_col]],
        dtype=np.int64,
    )
    keep = id_in_vocab >= 0
    print(
        f"[provider] Low-memory embed keeps {int(keep.sum())}/{adata.n_vars} "
        "genes in scGPT vocab"
    )
    if not np.any(keep):
        raise RuntimeError("No input genes are present in the scGPT vocabulary.")

    # Subset once, preserving sparse X. Copying this 80k x 1860 sparse matrix is
    # far cheaper than the dense conversion done by scGPT's helper.
    adata = adata[:, keep].copy()
    adata.var["id_in_vocab"] = id_in_vocab[keep]
    if sp.issparse(adata.X):
        adata.X = adata.X.tocsr()
    else:
        adata.X = np.asarray(adata.X, dtype=np.float32)

    with open(model_config_file, encoding="utf-8") as f:
        model_configs = json.load(f)

    genes = adata.var[gene_col].astype(str).tolist()
    gene_ids = np.array(vocab(genes), dtype=np.int64)

    model = TransformerModel(
        ntoken=len(vocab),
        d_model=model_configs["embsize"],
        nhead=model_configs["nheads"],
        d_hid=model_configs["d_hid"],
        nlayers=model_configs["nlayers"],
        nlayers_cls=model_configs["n_layers_cls"],
        n_cls=1,
        vocab=vocab,
        dropout=model_configs["dropout"],
        pad_token=model_configs["pad_token"],
        pad_value=model_configs["pad_value"],
        do_mvc=True,
        do_dab=False,
        use_batch_labels=False,
        domain_spec_batchnorm=False,
        explicit_zero_prob=False,
        use_fast_transformer=use_fast_transformer,
        fast_transformer_backend="flash",
        pre_norm=False,
    )
    device_obj = torch.device(device)
    state = torch.load(model_file, map_location=device_obj)
    load_pretrained(model, state, verbose=False)
    del state
    model.to(device_obj)
    model.eval()

    class SparseRowDataset(torch.utils.data.Dataset):
        def __init__(self, X, gene_ids):
            self.X = X
            self.gene_ids = gene_ids

        def __len__(self):
            return self.X.shape[0]

        def __getitem__(self, idx):
            row = self.X[idx]
            if sp.issparse(row):
                row = row.toarray().ravel()
            else:
                row = np.asarray(row).ravel()
            row = row.astype(np.float32, copy=False)
            nonzero_idx = np.nonzero(row)[0]
            values = row[nonzero_idx]
            genes = self.gene_ids[nonzero_idx]
            genes = np.insert(genes, 0, vocab["<cls>"])
            values = np.insert(values, 0, model_configs["pad_value"])
            return {
                "id": idx,
                "genes": torch.from_numpy(genes).long(),
                "expressions": torch.from_numpy(values).float(),
            }

    dataset = SparseRowDataset(adata.X, gene_ids)
    collator = DataCollator(
        do_padding=True,
        pad_token_id=vocab[model_configs["pad_token"]],
        pad_value=model_configs["pad_value"],
        do_mlm=False,
        do_binning=True,
        max_length=max_length,
        sampling=True,
        keep_first_n_tokens=1,
    )
    loader = DataLoader(
        dataset,
        batch_size=batch_size,
        sampler=SequentialSampler(dataset),
        collate_fn=collator,
        drop_last=False,
        num_workers=num_workers,
        pin_memory=(device == "cuda"),
    )

    cell_embeddings = np.zeros(
        (len(dataset), model_configs["embsize"]),
        dtype=np.float32,
    )
    count = 0
    autocast_enabled = device == "cuda"
    with torch.no_grad(), torch.cuda.amp.autocast(enabled=autocast_enabled):
        for data_dict in tqdm(loader, desc="Embedding cells"):
            input_gene_ids = data_dict["gene"].to(device_obj, non_blocking=True)
            src_key_padding_mask = input_gene_ids.eq(vocab[model_configs["pad_token"]])
            embeddings = model._encode(
                input_gene_ids,
                data_dict["expr"].to(device_obj, non_blocking=True),
                src_key_padding_mask=src_key_padding_mask,
                batch_labels=None,
            )
            embeddings = embeddings[:, 0, :].detach().cpu().numpy()
            cell_embeddings[count : count + len(embeddings)] = embeddings
            count += len(embeddings)

    norms = np.linalg.norm(cell_embeddings, axis=1, keepdims=True)
    cell_embeddings = cell_embeddings / np.maximum(norms, 1e-12)
    adata.obsm["X_scGPT"] = cell_embeddings.astype(np.float32, copy=False)
    return adata


def _normalized_entropy(values: pd.Series) -> float:
    counts = values.value_counts()
    if len(counts) <= 1:
        return 0.0
    return float(scipy_entropy(counts.values) / np.log(len(counts)))


def _assign_pseudostates(adata, args: argparse.Namespace) -> tuple[pd.DataFrame, str]:
    leiden_key = f"leiden_scgpt_res{args.leiden_resolution}"
    ps_key = "scgpt_pseudostate_provisional"
    family_key = "scgpt_state_family"
    status_key = "scgpt_state_status"

    clusters = sorted(adata.obs[leiden_key].astype(str).unique(), key=lambda x: int(x))
    cluster_order = []
    for cluster in clusters:
        mask = adata.obs[leiden_key].astype(str) == cluster
        obs = adata.obs.loc[mask]
        dominant_proxy = obs[args.stage_proxy_key].astype(str).value_counts().idxmax()
        median_time = float(pd.to_numeric(obs[args.time_key]).median())
        cluster_order.append((cluster, median_time, dominant_proxy))
    cluster_order.sort(key=lambda x: (x[1], x[0]))

    cluster_to_ps = {
        cluster: f"G178325_PS_{i:02d}"
        for i, (cluster, _, _) in enumerate(cluster_order)
    }

    rows = []
    leiden_str = adata.obs[leiden_key].astype(str)
    adata.obs[ps_key] = leiden_str.map(cluster_to_ps).astype("category")

    family_lookup = {}
    status_lookup = {}
    for cluster, _, _ in cluster_order:
        mask = leiden_str == cluster
        obs = adata.obs.loc[mask]
        proxy_counts = obs[args.stage_proxy_key].astype(str).value_counts()
        dominant_proxy = str(proxy_counts.idxmax())
        purity = float(proxy_counts.max() / proxy_counts.sum())
        time_entropy = _normalized_entropy(pd.to_numeric(obs[args.time_key]).astype(str))
        n_cells = int(mask.sum())
        status = "confirmed"
        if purity < 0.45 or n_cells < max(50, int(0.002 * adata.n_obs)):
            status = "uncertain"
        elif purity < 0.70 or time_entropy > 0.65:
            status = "merge_review"
        ps = cluster_to_ps[cluster]
        family_lookup[ps] = dominant_proxy
        status_lookup[ps] = status
        rows.append(
            {
                "state_id": ps,
                "leiden_cluster": cluster,
                "family": dominant_proxy,
                "status": status,
                "n_cells": n_cells,
                "dominant_stage_proxy": dominant_proxy,
                "stage_proxy_purity": purity,
                "time_entropy_normalized": time_entropy,
                "time_min": float(pd.to_numeric(obs[args.time_key]).min()),
                "time_median": float(pd.to_numeric(obs[args.time_key]).median()),
                "time_max": float(pd.to_numeric(obs[args.time_key]).max()),
                "version": VERSION,
            }
        )

    adata.obs[family_key] = adata.obs[ps_key].astype(str).map(family_lookup).astype("category")
    adata.obs[status_key] = adata.obs[ps_key].astype(str).map(status_lookup).astype("category")

    state_meta = pd.DataFrame(rows).sort_values(["time_median", "state_id"])
    return state_meta, ps_key


def _edge_confidence(src_status: str, tgt_status: str) -> str:
    if src_status == "uncertain" or tgt_status == "uncertain":
        return "low"
    if src_status == "confirmed" and tgt_status == "confirmed":
        return "high"
    return "medium"


def _build_reference_graph(
    adata,
    state_meta: pd.DataFrame,
    ps_key: str,
    args: argparse.Namespace,
) -> tuple[list[dict], pd.DataFrame]:
    emb = np.asarray(adata.obsm["X_scGPT"], dtype=np.float32)
    norms = np.linalg.norm(emb, axis=1, keepdims=True)
    emb = emb / np.maximum(norms, 1e-12)

    states = list(state_meta["state_id"])
    state_idx = {s: i for i, s in enumerate(states)}
    counts = np.zeros((len(states), len(states)), dtype=np.int64)
    pair_counts = np.zeros((len(states), len(states)), dtype=np.int64)

    time_values = pd.to_numeric(adata.obs[args.time_key]).to_numpy(dtype=float)
    state_values = adata.obs[ps_key].astype(str).to_numpy()
    timepoints = sorted(np.unique(time_values))

    for t0, t1 in zip(timepoints[:-1], timepoints[1:]):
        src_idx = np.where(np.isclose(time_values, t0))[0]
        tgt_idx = np.where(np.isclose(time_values, t1))[0]
        if src_idx.size == 0 or tgt_idx.size == 0:
            continue
        tgt_emb = emb[tgt_idx]
        tgt_states = state_values[tgt_idx]
        k = min(args.knn_k, tgt_idx.size)
        pair_seen = set()
        for start in range(0, src_idx.size, args.chunk_size):
            chunk_idx = src_idx[start : start + args.chunk_size]
            sim = emb[chunk_idx] @ tgt_emb.T
            nn_local = np.argpartition(-sim, kth=k - 1, axis=1)[:, :k]
            for row_i, neighbors in enumerate(nn_local):
                src_state = state_values[chunk_idx[row_i]]
                src_pos = state_idx[src_state]
                for local_j in neighbors:
                    tgt_state = tgt_states[local_j]
                    tgt_pos = state_idx[tgt_state]
                    counts[src_pos, tgt_pos] += 1
                    pair_seen.add((src_pos, tgt_pos))
        for src_pos, tgt_pos in pair_seen:
            pair_counts[src_pos, tgt_pos] += 1
        print(
            f"[provider] kNN transitions {t0:g}->{t1:g}: "
            f"{src_idx.size} source, {tgt_idx.size} target"
        )

    rows = []
    status = dict(zip(state_meta["state_id"], state_meta["status"]))
    with np.errstate(divide="ignore", invalid="ignore"):
        weights = counts / counts.sum(axis=1, keepdims=True)
    weights = np.nan_to_num(weights, nan=0.0)

    for i, src in enumerate(states):
        for j, tgt in enumerate(states):
            if src == tgt:
                continue
            weight = float(weights[i, j])
            if weight < args.min_edge_weight:
                continue
            rows.append(
                {
                    "source": src,
                    "target": tgt,
                    "weight": round(weight, 6),
                    "raw_knn_count": int(counts[i, j]),
                    "n_timepoint_pairs": int(pair_counts[i, j]),
                    "confidence": _edge_confidence(status[src], status[tgt]),
                    "source_status": status[src],
                    "target_status": status[tgt],
                    "edge_type": "scgpt_consecutive_timepoint_knn",
                }
            )

    edges_df = pd.DataFrame(rows).sort_values(["source", "target"])
    edges = edges_df.to_dict(orient="records")
    return edges, edges_df


def _write_outputs(
    adata,
    output_h5ad: Path,
    state_meta: pd.DataFrame,
    ps_key: str,
    edges: list[dict],
    edges_df: pd.DataFrame,
    args: argparse.Namespace,
    embed_meta: dict,
    out_dir: Path,
    provider_dir: Path,
) -> None:
    provider_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    annotated_h5ad = out_dir / "adata_scgpt_annotated.h5ad"
    adata.write_h5ad(annotated_h5ad)

    label_path = provider_dir / "state_labels.tsv"
    label_df = pd.DataFrame(
        {
            "cell_id": adata.obs.index.astype(str),
            "scgpt_pseudostate_v1": adata.obs[ps_key].astype(str).values,
            "scgpt_pseudostate_provisional": adata.obs[ps_key].astype(str).values,
        }
    )
    label_df.to_csv(label_path, sep="\t", index=False)

    metadata_path = provider_dir / "state_metadata.tsv"
    state_meta.to_csv(metadata_path, sep="\t", index=False)

    edges_path = provider_dir / "reference_graph_edges.csv"
    edges_df.to_csv(edges_path, index=False)

    nodes = []
    for row in state_meta.to_dict(orient="records"):
        nodes.append(
            {
                "id": row["state_id"],
                "n_cells": int(row["n_cells"]),
                "family": row["family"],
                "status": row["status"],
                "version": VERSION,
                "leiden_cluster": str(row["leiden_cluster"]),
                "time_range": [float(row["time_min"]), float(row["time_max"])],
            }
        )

    graph = {
        "_meta": {
            "description": "Working v1 scGPT-derived reference lineage graph for GSE178325 0618.",
            "version": VERSION,
            "status": "silver_standard_working",
            "dataset": "GSE178325",
            "batch_policy": "0618 only",
            "graph_method": "consecutive_timepoint_kNN_transition_counting",
            "knn_k": args.knn_k,
            "min_edge_weight": args.min_edge_weight,
            "state_key": ps_key,
            "time_key": args.time_key,
            "embedding_source": "X_scGPT (scgpt_whole_human pretrained model)",
            "leiden_resolution": args.leiden_resolution,
            "n_cells_total": int(adata.n_obs),
            "n_states": int(len(nodes)),
            "n_edges": int(len(edges)),
            "created_from": str(args.input_h5ad),
            "policy": [
                "This graph is a silver-standard working asset, not final biology.",
                "States are Leiden clusters in scGPT embedding space ordered by median abs_day.",
                "Edges are consecutive-timepoint kNN transition probabilities in X_scGPT.",
                "Self-loops are excluded and edges below min_edge_weight are pruned.",
            ],
        },
        "nodes": nodes,
        "edges": edges,
    }
    graph_path = provider_dir / "reference_graph.json"
    graph_path.write_text(json.dumps(graph, indent=2), encoding="utf-8")

    cluster_time = pd.crosstab(
        adata.obs[ps_key].astype(str),
        pd.to_numeric(adata.obs[args.time_key]),
    )
    cluster_time.to_csv(out_dir / "pseudostate_timepoint_matrix.csv")

    metadata = {
        "run_timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "script": "benchmark/scgpt/build_gse178325_0618_provider.py",
        "provider_id": PROVIDER_ID,
        "version": VERSION,
        "status": "silver_standard_working",
        "input_h5ad": str(args.input_h5ad),
        "output_h5ad": str(output_h5ad),
        "annotated_h5ad": str(annotated_h5ad),
        "provider_dir": str(provider_dir),
        "state_key": ps_key,
        "time_key": args.time_key,
        "n_cells": int(adata.n_obs),
        "n_states": int(len(nodes)),
        "n_graph_edges": int(len(edges)),
        "n_edges_by_confidence": edges_df["confidence"].value_counts().to_dict()
        if not edges_df.empty
        else {},
        "embedding": embed_meta,
        "outputs": {
            "label_tsv": str(label_path),
            "metadata_tsv": str(metadata_path),
            "graph_json": str(graph_path),
            "graph_csv": str(edges_path),
            "pseudostate_timepoint_matrix": str(out_dir / "pseudostate_timepoint_matrix.csv"),
        },
    }
    (provider_dir / "ground_truth_metadata.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )

    if "X_umap" in adata.obsm:
        sc.settings.figdir = str(out_dir)
        adata.obs["abs_day_str"] = pd.to_numeric(adata.obs[args.time_key]).astype(str)
        sc.pl.umap(adata, color="abs_day_str", show=False, save="_gse178325_by_timepoint.png")
        sc.pl.umap(adata, color=ps_key, show=False, save="_gse178325_by_pseudostate.png")

    print("[provider] Wrote provider files:")
    print(f"  {graph_path}")
    print(f"  {edges_path}")
    print(f"  {metadata_path}")
    print(f"  {label_path}")


def main() -> None:
    args = _parse_args()
    out_dir = Path(args.output_dir)
    provider_dir = Path(args.provider_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    adata, output_h5ad, embed_meta = _load_or_embed(args, out_dir)
    state_meta, ps_key = _assign_pseudostates(adata, args)
    print(f"[provider] Pseudostates: {state_meta.shape[0]}")
    print(state_meta[["state_id", "family", "status", "n_cells", "time_median"]])

    edges, edges_df = _build_reference_graph(adata, state_meta, ps_key, args)
    print(f"[provider] Reference edges after pruning: {len(edges)}")
    if not edges_df.empty:
        print(edges_df.head(20).to_string(index=False))

    _write_outputs(
        adata=adata,
        output_h5ad=output_h5ad,
        state_meta=state_meta,
        ps_key=ps_key,
        edges=edges,
        edges_df=edges_df,
        args=args,
        embed_meta=embed_meta,
        out_dir=out_dir,
        provider_dir=provider_dir,
    )


if __name__ == "__main__":
    main()
