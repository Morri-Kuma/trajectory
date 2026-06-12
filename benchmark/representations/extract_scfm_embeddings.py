"""
extract_scfm_embeddings.py
==========================

CLI to extract frozen single-cell foundation model (scFM) cell embeddings and
write them into an AnnData object, with full provenance metadata.

This is the implementation of Work-Plan Step 1. It is *reference-first*: the
actual embedding computation is delegated to the official package for each
model (see ``scfm_extractors.py``). This script is responsible only for I/O,
argument handling, and metadata assembly.

Outputs (into the output h5ad):
    adata.obsm["X_geneformer_cls"]      (model=geneformer)
    adata.obsm["X_scgpt_cls"]           (model=scgpt)
    adata.obsm["X_scfoundation"]        (model=scfoundation)
    adata.uns["scfm_embedding_metadata"][<model>] = <metadata block>

A standalone ``representation_metadata.json`` is also written next to the
output so provenance survives independently of the h5ad.

Examples (run on the Shirokane GPU host with the official packages installed)
----------------------------------------------------------------------------
    python -m benchmark.representations.extract_scfm_embeddings \
        --model scgpt \
        --input-h5ad  inputs/GSE230659_full_gene.h5ad \
        --output-h5ad inputs/GSE230659_scgpt.h5ad \
        --model-path  checkpoints/scGPT_human \
        --gene-col feature_name --batch-size 64

    python -m benchmark.representations.extract_scfm_embeddings \
        --model geneformer \
        --input-h5ad  inputs/GSE230659_full_gene.h5ad \
        --output-h5ad inputs/GSE230659_geneformer.h5ad \
        --model-path  checkpoints/Geneformer-V2-104M \
        --ensembl-col ensembl_id --counts-col n_counts

    python -m benchmark.representations.extract_scfm_embeddings \
        --model scfoundation \
        --input-h5ad  inputs/GSE230659_full_gene.h5ad \
        --output-h5ad inputs/GSE230659_scfoundation.h5ad \
        --model-path  checkpoints/scFoundation/models.ckpt \
        --scfoundation-repo /path/to/scFoundation/model \
        --gene-index-tsv /path/to/OS_scRNA_gene_index.19264.tsv
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .metadata import build_embedding_metadata
from .references import SUPPORTED_MODELS, get_reference
from .scfm_extractors import (
    EXTRACTORS,
    ScFMDependencyError,
    ScFMInputError,
)


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Extract frozen scFM cell embeddings (reference-first).",
    )
    p.add_argument("--model", required=True, choices=list(SUPPORTED_MODELS))
    p.add_argument("--input-h5ad", required=True, help="Input AnnData (.h5ad).")
    out = p.add_mutually_exclusive_group(required=True)
    out.add_argument("--output-h5ad", help="Output AnnData path.")
    out.add_argument(
        "--output-dir",
        help="Directory for outputs; writes <model>_embeddings.h5ad inside.",
    )

    # checkpoint / model paths
    p.add_argument(
        "--model-path",
        help="Model directory (scgpt/geneformer) or checkpoint file (scfoundation).",
    )
    p.add_argument(
        "--upstream-commit-or-release",
        default="unknown",
        help="Commit hash or release tag of the official code/checkpoint used.",
    )

    # vocab / tokenizer options
    p.add_argument("--gene-col", default="feature_name",
                   help="[scgpt] var column with gene names (or 'index').")
    p.add_argument("--ensembl-col", default="ensembl_id",
                   help="[geneformer] var column with Ensembl IDs.")
    p.add_argument("--counts-col", default="n_counts",
                   help="[geneformer] obs column with per-cell total counts.")
    p.add_argument("--token-dictionary-file", default=None,
                   help="[geneformer] custom token dictionary pickle.")
    p.add_argument("--emb-mode", default="cls", choices=["cls", "cell"],
                   help="[geneformer] CLS-token or mean-pool cell embedding.")
    p.add_argument("--emb-layer", type=int, default=-1, choices=[-1, 0],
                   help="[geneformer] -1: 2nd-to-last layer, 0: last layer.")
    p.add_argument("--model-version", default="V2",
                   help="[geneformer] V1 (~30M) or V2 (~104M).")
    p.add_argument("--scfoundation-repo", default=None,
                   help="[scfoundation] path to official 'model' folder "
                        "(containing get_embedding.py, OS_scRNA_gene_index.19264.tsv, "
                        "models/models.ckpt).")
    p.add_argument("--gene-index-tsv", default=None,
                   help="[scfoundation] OS_scRNA_gene_index.19264.tsv "
                        "(default: <repo>/OS_scRNA_gene_index.19264.tsv).")
    p.add_argument("--gene-name-col", default=None,
                   help="[scfoundation] var column with gene names (default var_names).")
    p.add_argument("--pool-type", default="all", choices=["all", "max"],
                   help="[scfoundation] cell-embedding pooling.")
    p.add_argument("--tgthighres", default="t4",
                   help="[scfoundation] target high-resolution token.")
    p.add_argument("--scfoundation-version", default="ce", choices=["ce", "rde"],
                   help="[scfoundation] official cell-embedding version key.")
    p.add_argument("--scfoundation-python", default=None,
                   help="[scfoundation] python interpreter of the scFoundation "
                        "environment used to launch get_embedding.py "
                        "(default: this interpreter).")
    p.add_argument("--scfoundation-task-name", default="rep",
                   help="[scfoundation] task_name label for the official output file.")
    p.add_argument("--scfoundation-ckpt-name", default="models",
                   help="[scfoundation] ckpt_name label for the official output file.")

    # runtime
    p.add_argument("--batch-size", type=int, default=64)
    p.add_argument("--device", default="cuda")
    p.add_argument("--use-fast-transformer", action="store_true",
                   help="[scgpt] use flash-attn backend.")
    p.add_argument("--nproc", type=int, default=4, help="[geneformer] CPU procs.")
    return p


def _resolve_output(args) -> Path:
    if args.output_h5ad:
        return Path(args.output_h5ad)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    return out_dir / f"{args.model}_embeddings.h5ad"


def _call_extractor(model: str, adata, args):
    if model == "geneformer":
        return EXTRACTORS[model](
            adata,
            model_dir=args.model_path,
            token_dictionary_file=args.token_dictionary_file,
            emb_mode=args.emb_mode,
            emb_layer=args.emb_layer,
            model_version=args.model_version,
            ensembl_col=args.ensembl_col,
            counts_col=args.counts_col,
            forward_batch_size=args.batch_size,
            nproc=args.nproc,
            upstream_commit_or_release=args.upstream_commit_or_release,
        )
    if model == "scgpt":
        return EXTRACTORS[model](
            adata,
            model_dir=args.model_path,
            gene_col=args.gene_col,
            batch_size=args.batch_size,
            device=args.device,
            use_fast_transformer=args.use_fast_transformer,
            upstream_commit_or_release=args.upstream_commit_or_release,
        )
    if model == "scfoundation":
        if not args.scfoundation_repo:
            raise ScFMInputError(
                "[scfoundation] --scfoundation-repo (official model folder) is "
                "required. The extractor invokes the official get_embedding.py."
            )
        return EXTRACTORS[model](
            adata,
            scfoundation_repo=args.scfoundation_repo,
            model_path=args.model_path,
            gene_index_tsv=args.gene_index_tsv,
            gene_name_col=args.gene_name_col,
            pool_type=args.pool_type,
            tgthighres=args.tgthighres,
            version=args.scfoundation_version,
            task_name=args.scfoundation_task_name,
            ckpt_name=args.scfoundation_ckpt_name,
            python_executable=args.scfoundation_python,
            device=args.device,
            upstream_commit_or_release=args.upstream_commit_or_release,
        )
    raise ValueError(f"Unsupported model {model!r}.")


def run(args) -> Path:
    try:
        import anndata
    except Exception as exc:  # pragma: no cover
        raise ScFMDependencyError(
            f"anndata is required to read/write h5ad: {exc}"
        ) from exc

    in_path = Path(args.input_h5ad)
    if not in_path.exists():
        raise ScFMInputError(f"--input-h5ad not found: {in_path}")
    adata = anndata.read_h5ad(in_path)

    ref = get_reference(args.model)
    embeddings, run_meta = _call_extractor(args.model, adata, args)

    if embeddings.shape[0] != adata.n_obs:
        raise ScFMInputError(
            f"[{args.model}] extractor returned {embeddings.shape[0]} embeddings "
            f"for {adata.n_obs} cells; cell-order alignment failed."
        )

    meta = build_embedding_metadata(
        reference=ref,
        model_checkpoint=run_meta["model_checkpoint"],
        upstream_commit_or_release=run_meta["upstream_commit_or_release"],
        input_gene_count=run_meta["input_gene_count"],
        matched_gene_count=run_meta["matched_gene_count"],
        embedding_dim=run_meta["embedding_dim"],
        extra=run_meta.get("extra"),
    )

    obsm_key = run_meta["obsm_key"]
    adata.obsm[obsm_key] = embeddings
    existing = dict(adata.uns.get("scfm_embedding_metadata", {}))
    existing[args.model] = meta
    adata.uns["scfm_embedding_metadata"] = existing

    out_path = _resolve_output(args)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out_path)

    sidecar = out_path.with_name(
        out_path.stem + ".representation_metadata.json"
    )
    with open(sidecar, "w", encoding="utf-8") as f:
        json.dump({args.model: meta}, f, indent=2)

    print(
        f"[extract_scfm_embeddings] {args.model}: wrote obsm['{obsm_key}'] "
        f"({embeddings.shape[0]}x{embeddings.shape[1]}) -> {out_path}\n"
        f"[extract_scfm_embeddings] metadata -> {sidecar}"
    )
    return out_path


def main(argv=None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        run(args)
    except (ScFMDependencyError, ScFMInputError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
