"""
scfm_extractors.py
==================

Thin, reference-grounded adapters that extract *frozen* single-cell foundation
model (scFM) cell embeddings. Each adapter calls the official package / code for
its model and does NOT re-implement tokenization, vocabulary matching, model
loading, pooling, or embedding-layer selection.

If the official dependency or checkpoint is unavailable, the adapter raises a
clear, actionable error (``ScFMDependencyError`` / ``ScFMInputError``) rather
than substituting an invented approximation. This is by design: extraction is
expected to run on the Shirokane GPU server where the official packages and
checkpoints are installed.

Official sources inspected (also recorded in ``references.MODEL_REFERENCES``):
  - Geneformer    : https://www.nature.com/articles/s41586-023-06139-9
                    https://huggingface.co/ctheodoris/Geneformer
                    geneformer.TranscriptomeTokenizer + geneformer.EmbExtractor
                    (https://geneformer.readthedocs.io/en/latest/geneformer.tokenizer.html,
                     https://geneformer.readthedocs.io/en/latest/geneformer.emb_extractor.html)
  - scGPT         : https://www.nature.com/articles/s41592-024-02201-0
                    https://github.com/bowang-lab/scGPT
                    scgpt.tasks.cell_emb.embed_data
                    (https://github.com/bowang-lab/scGPT/blob/main/scgpt/tasks/cell_emb.py)
  - scFoundation  : https://www.nature.com/articles/s41592-024-02305-7
                    https://github.com/biomap-research/scFoundation
                    model/get_embedding.py (invoked as a subprocess; we do NOT
                    re-implement its tokenization / encoder / pooling)
                    (https://github.com/biomap-research/scFoundation/blob/main/model/get_embedding.py,
                     https://github.com/biomap-research/scFoundation/blob/main/model/load.py)
"""

from __future__ import annotations

import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional

import numpy as np

from .references import get_reference


class ScFMDependencyError(RuntimeError):
    """Official scFM package or checkpoint is not available."""


class ScFMInputError(ValueError):
    """Input AnnData is missing fields the official extractor requires."""


def _import_or_raise(module_name: str, model: str, install_hint: str):
    try:
        return __import__(module_name)
    except Exception as exc:  # pragma: no cover - exercised on Shirokane only
        raise ScFMDependencyError(
            f"[{model}] Could not import the official package "
            f"'{module_name}'. This extractor is a thin wrapper over the "
            f"official implementation and will not substitute an approximation. "
            f"Install/enable it on the run host: {install_hint}. "
            f"Underlying import error: {exc}"
        ) from exc


def _require_path(path: Optional[str], model: str, what: str) -> Path:
    if not path:
        raise ScFMInputError(
            f"[{model}] A {what} path is required but none was provided."
        )
    p = Path(path)
    if not p.exists():
        raise ScFMDependencyError(
            f"[{model}] {what} not found at '{p}'. Provide the official "
            f"checkpoint/model files on the run host."
        )
    return p


# ===========================================================================
# Geneformer
# ===========================================================================
#
# Official tokenizer (geneformer.TranscriptomeTokenizer) hard-requires:
#   - row/gene attribute  adata.var["ensembl_id"]  (fixed name)
#   - col/cell attribute  adata.obs["n_counts"]    (fixed name)
# See https://geneformer.readthedocs.io/en/latest/geneformer.tokenizer.html
# If the caller stored these under different column names, we copy them into the
# official fixed names on a temp AnnData *before* tokenization (we never rename
# the caller's object). CLS embeddings (emb_mode="cls") require a V2 (CLS-token)
# checkpoint; V1 has no CLS token.
# ---------------------------------------------------------------------------

def _geneformer_prepare_official_adata(adata, ensembl_col, counts_col, index_col):
    """
    Return a copy of ``adata`` with the official fixed columns populated:
      var["ensembl_id"] <- var[ensembl_col]   (if a custom column was supplied)
      obs["n_counts"]   <- obs[counts_col]     (if a custom column was supplied)
    and a stable ``index_col`` for restoring cell order after extraction.

    This is pure pandas/AnnData column copying (no model dependency) so it is
    unit-testable without the geneformer package.
    """
    a = adata.copy()
    if ensembl_col != "ensembl_id":
        a.var["ensembl_id"] = list(a.var[ensembl_col])
    if counts_col != "n_counts":
        a.obs["n_counts"] = list(a.obs[counts_col])
    n = a.obs.shape[0]
    a.obs[index_col] = np.arange(n)
    return a


def _geneformer_matched_gene_count(ensembl_ids, gene_token_dict) -> int:
    """
    Number of input genes whose Ensembl ID is present in the official token
    dictionary (Ensembl ID -> token). This is the tokenizer-consistent count of
    genes that can actually be tokenized, NOT the raw input gene count.
    """
    keys = set(gene_token_dict.keys())
    return int(sum(1 for g in ensembl_ids if g in keys))


def extract_geneformer(
    adata,
    *,
    model_dir: str,
    token_dictionary_file: Optional[str] = None,
    emb_mode: str = "cls",
    emb_layer: int = -1,
    model_version: str = "V2",
    ensembl_col: str = "ensembl_id",
    counts_col: str = "n_counts",
    forward_batch_size: int = 100,
    nproc: int = 4,
    upstream_commit_or_release: str = "unknown",
):
    """
    Extract Geneformer cell embeddings using the official ``geneformer`` package
    (``TranscriptomeTokenizer`` + ``EmbExtractor``).

    Mirrors the official API exactly:
      TranscriptomeTokenizer(custom_attr_name_dict=..., model_version=...,
                             [token_dictionary_file=...]).tokenize_data(...)
      EmbExtractor(model_type='Pretrained', emb_mode='cls'|'cell',
                   cell_emb_style='mean_pool', emb_layer=-1,
                   model_version=...).extract_embs(...)

    Returns
    -------
    (embeddings, run_meta) : (np.ndarray [n_cells, emb_dim], dict)
    """
    ref = get_reference("geneformer")

    if emb_mode == "cls" and model_version != "V2":
        raise ScFMInputError(
            "[geneformer] emb_mode='cls' requires a V2 (CLS-token) checkpoint "
            f"(model_version='V2'), but model_version={model_version!r}. The V1 "
            "model series has no CLS token; use emb_mode='cell' (mean-pool) for "
            "V1, or select a V2 checkpoint."
        )

    if ensembl_col not in adata.var.columns:
        raise ScFMInputError(
            f"[geneformer] adata.var must contain Ensembl IDs in column "
            f"'{ensembl_col}' (the official TranscriptomeTokenizer keys genes by "
            f"Ensembl ID via adata.var['ensembl_id']). "
            f"Available var columns: {list(adata.var.columns)}."
        )
    if counts_col not in adata.obs.columns:
        raise ScFMInputError(
            f"[geneformer] adata.obs must contain per-cell total counts in column "
            f"'{counts_col}' (the official TranscriptomeTokenizer requires "
            f"adata.obs['n_counts']). Available obs columns: {list(adata.obs.columns)}."
        )

    _import_or_raise(
        "geneformer", "geneformer",
        "pip install git+https://huggingface.co/ctheodoris/Geneformer",
    )
    model_path = _require_path(model_dir, "geneformer", "model directory")

    from geneformer import EmbExtractor, TranscriptomeTokenizer  # type: ignore

    input_gene_count = int(adata.n_vars)

    # Copy custom columns into the official fixed names on a temp copy and add a
    # stable index so we can restore the original cell order from the output.
    index_col = "__rep_cell_index"
    prepared = _geneformer_prepare_official_adata(
        adata, ensembl_col, counts_col, index_col
    )
    ensembl_ids = list(prepared.var["ensembl_id"])

    tk_kwargs = {"model_version": model_version}
    if token_dictionary_file:
        tk_kwargs["token_dictionary_file"] = token_dictionary_file

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        data_dir = tmp / "h5ad"
        data_dir.mkdir()
        token_dir = tmp / "tokenized"
        token_dir.mkdir()
        out_dir = tmp / "embs"
        out_dir.mkdir()

        prepared.write_h5ad(data_dir / "input.h5ad")

        tk = TranscriptomeTokenizer(
            custom_attr_name_dict={index_col: index_col},
            nproc=nproc,
            **tk_kwargs,
        )
        # Tokenizer-consistent matched-gene count (Ensembl IDs present in the
        # official token dictionary). NOT the raw input gene count.
        matched_gene_count = _geneformer_matched_gene_count(
            ensembl_ids, tk.gene_token_dict
        )

        tk.tokenize_data(str(data_dir), str(token_dir), "rep", file_format="h5ad")
        dataset_path = token_dir / "rep.dataset"

        embex = EmbExtractor(
            model_type="Pretrained",
            num_classes=0,
            emb_mode=emb_mode,
            cell_emb_style="mean_pool",
            emb_layer=emb_layer,
            max_ncells=None,  # all cells
            emb_label=[index_col],
            forward_batch_size=forward_batch_size,
            nproc=nproc,
            model_version=model_version,
            **(
                {"token_dictionary_file": token_dictionary_file}
                if token_dictionary_file
                else {}
            ),
        )
        embs_df = embex.extract_embs(
            str(model_path),
            str(dataset_path),
            str(out_dir),
            "rep_embs",
            output_torch_embs=False,
        )

    embs_df = embs_df.sort_values(index_col)
    emb_cols = [c for c in embs_df.columns if c != index_col]
    embeddings = embs_df[emb_cols].to_numpy(dtype=np.float32)

    run_meta = {
        "model_checkpoint": str(model_path),
        "upstream_commit_or_release": upstream_commit_or_release,
        "input_gene_count": input_gene_count,
        "matched_gene_count": matched_gene_count,
        "embedding_dim": int(embeddings.shape[1]),
        "obsm_key": ref["obsm_key"],
        "extra": {
            "emb_mode": emb_mode,
            "emb_layer": emb_layer,
            "model_version": model_version,
            "ensembl_col_source": ensembl_col,
            "counts_col_source": counts_col,
            "matched_gene_count_basis": "ensembl_ids_in_official_token_dictionary",
        },
    }
    return embeddings, run_meta


# ===========================================================================
# scGPT
# ===========================================================================
#
# Official cell-embedding entry point is scgpt.tasks.cell_emb.embed_data, which
# internally performs vocab matching, value binning, the model forward,
# CLS-token (first position) extraction and per-cell L2 normalization. We
# delegate entirely to it and never re-implement that path.
#
# Official embed_data signature (github.com/bowang-lab/scGPT, cell_emb.py):
#   embed_data(adata_or_file, model_dir, gene_col="feature_name", max_length=1200,
#              batch_size=64, obs_to_save=None, device="cuda",
#              use_fast_transformer=True, return_new_adata=False)
# The official default use_fast_transformer=True requires flash-attn + a
# compatible GPU. We default to False for portability on run hosts without
# flash-attn and RECORD this override explicitly in metadata (set
# use_fast_transformer=True via the CLI to match the official default).
# ---------------------------------------------------------------------------

SCGPT_OFFICIAL_DEFAULT_USE_FAST_TRANSFORMER = True


def _scgpt_run_extra(use_fast_transformer: bool, max_length: int, gene_col: str) -> dict:
    """Build the scGPT provenance ``extra`` dict, recording any override of the
    official embed_data default for use_fast_transformer."""
    extra = {
        "max_length": max_length,
        "gene_col": gene_col,
        "l2_normalized": True,  # official get_batch_cell_embeddings L2-normalizes
        "cell_embedding_mode": "cls",
        "use_fast_transformer": bool(use_fast_transformer),
        "official_default_use_fast_transformer":
            SCGPT_OFFICIAL_DEFAULT_USE_FAST_TRANSFORMER,
    }
    if bool(use_fast_transformer) != SCGPT_OFFICIAL_DEFAULT_USE_FAST_TRANSFORMER:
        extra["use_fast_transformer_override"] = (
            "Overrode official embed_data default "
            f"(use_fast_transformer={SCGPT_OFFICIAL_DEFAULT_USE_FAST_TRANSFORMER}) "
            f"-> {bool(use_fast_transformer)} for portability on hosts without "
            "flash-attn. Computation is otherwise the official embed_data path."
        )
    return extra


def extract_scgpt(
    adata,
    *,
    model_dir: str,
    gene_col: str = "feature_name",
    max_length: int = 1200,
    batch_size: int = 64,
    device: str = "cuda",
    use_fast_transformer: bool = False,
    upstream_commit_or_release: str = "unknown",
):
    """
    Extract scGPT CLS cell embeddings via the official
    ``scgpt.tasks.cell_emb.embed_data`` (cell_embedding_mode='cls').

    The official function computes the CLS-position embedding and L2-normalizes
    it. With ``return_new_adata=True`` it returns an AnnData whose ``X`` is the
    cell-embedding matrix; we return that array unchanged and re-key it to
    ``X_scgpt_cls`` at the caller.
    """
    ref = get_reference("scgpt")
    if gene_col != "index" and gene_col not in adata.var.columns:
        raise ScFMInputError(
            f"[scgpt] adata.var must contain gene names in column '{gene_col}' "
            f"(matched against the checkpoint's vocab.json), or pass "
            f"--gene-col index. Available var columns: {list(adata.var.columns)}."
        )
    _import_or_raise(
        "scgpt", "scgpt", "pip install scgpt  (see github.com/bowang-lab/scGPT)"
    )
    model_path = _require_path(model_dir, "scgpt", "model directory")
    # The official loader expects vocab.json / args.json / best_model.pt inside.
    for required in ("vocab.json", "args.json", "best_model.pt"):
        if not (model_path / required).exists():
            raise ScFMDependencyError(
                f"[scgpt] Expected '{required}' inside the scGPT model "
                f"directory '{model_path}' (official embed_data layout)."
            )

    from scgpt.tasks.cell_emb import embed_data  # type: ignore

    input_gene_count = int(adata.n_vars)
    out = embed_data(
        adata,
        str(model_path),
        gene_col=gene_col,
        max_length=max_length,
        batch_size=batch_size,
        device=device,
        use_fast_transformer=use_fast_transformer,
        return_new_adata=True,  # returns AnnData(X=cell_embeddings)
    )
    embeddings = np.asarray(out.X, dtype=np.float32)

    # embed_data drops genes with id_in_vocab < 0 and annotates the input adata
    # in place; recover the matched count from that official annotation.
    matched = input_gene_count
    if "id_in_vocab" in adata.var.columns:
        matched = int(np.sum(np.asarray(adata.var["id_in_vocab"]) >= 0))

    run_meta = {
        "model_checkpoint": str(model_path),
        "upstream_commit_or_release": upstream_commit_or_release,
        "input_gene_count": input_gene_count,
        "matched_gene_count": matched,
        "embedding_dim": int(embeddings.shape[1]),
        "obsm_key": ref["obsm_key"],
        "extra": _scgpt_run_extra(use_fast_transformer, max_length, gene_col),
    }
    return embeddings, run_meta


# ===========================================================================
# scFoundation
# ===========================================================================
#
# scFoundation ships as scripts (no installable package). The official
# cell-embedding entry point is model/get_embedding.py, which performs gene-panel
# alignment (its own local main_gene_selection), per-cell normalization, the two
# read-depth tokens, the xTrimoGene encoder forward, and pooling. We invoke that
# script AS-IS via subprocess and read its .npy output. We deliberately do NOT
# re-implement any of its tokenization / encoder / pooling logic here.
#
# Official invocation (input_type=singlecell, output_type=cell, version=ce):
#   python get_embedding.py --task_name <t> --input_type singlecell \
#       --output_type cell --pool_type all --tgthighres t4 --pre_normalized F \
#       --version ce --data_path <input> --save_path <out> --ckpt_name <name>
# It reads ./OS_scRNA_gene_index.19264.tsv and ./models/models.ckpt relative to
# its own directory (cwd), so we run with cwd set to the official model folder.
# Output file:
#   <save_path>/<task_name>_<ckpt_name>_singlecell_cell_embedding_<tgthighres>_resolution.npy
# Refs: https://github.com/biomap-research/scFoundation/blob/main/model/get_embedding.py
#       https://github.com/biomap-research/scFoundation/blob/main/model/load.py
# ---------------------------------------------------------------------------

def _scfoundation_output_filename(task_name, ckpt_name, tgthighres) -> str:
    """Official get_embedding.py output naming (cell, singlecell)."""
    return (
        f"{task_name}_{ckpt_name}_singlecell_cell_embedding_"
        f"{tgthighres}_resolution.npy"
    )


def _scfoundation_matched_gene_count(gene_names, panel_gene_names) -> int:
    """Genes present in the official 19,264-gene panel (set intersection)."""
    return int(len(set(gene_names) & set(panel_gene_names)))


def extract_scfoundation(
    adata,
    *,
    scfoundation_repo: str,
    model_path: Optional[str] = None,
    gene_index_tsv: Optional[str] = None,
    gene_name_col: Optional[str] = None,
    pool_type: str = "all",
    tgthighres: str = "t4",
    pre_normalized: str = "F",
    version: str = "ce",
    task_name: str = "rep",
    ckpt_name: str = "models",
    python_executable: Optional[str] = None,
    device: str = "cuda",
    upstream_commit_or_release: str = "unknown",
):
    """
    Extract scFoundation cell embeddings by invoking the official
    ``model/get_embedding.py`` as a subprocess (thin wrapper; no re-implemented
    tokenization / encoder / pooling).

    ``scfoundation_repo`` is the official ``model`` folder containing
    ``get_embedding.py``, ``load.py``, ``pretrainmodels.py``,
    ``OS_scRNA_gene_index.19264.tsv``, and ``models/models.ckpt`` (the checkpoint
    layout the official ``version='ce'`` cell path loads from ``./models/models.ckpt``).

    Returns
    -------
    (embeddings, run_meta)
    """
    ref = get_reference("scfoundation")

    repo = _require_path(scfoundation_repo, "scfoundation", "official model folder")
    get_emb = repo / "get_embedding.py"
    if not get_emb.exists():
        raise ScFMDependencyError(
            f"[scfoundation] official get_embedding.py not found in "
            f"'{repo}'. Pass --scfoundation-repo pointing at the official "
            f"scFoundation 'model' folder."
        )

    # Gene panel (read directly; needed for matched-gene count and required by
    # the official script as ./OS_scRNA_gene_index.19264.tsv).
    panel_path = (
        Path(gene_index_tsv) if gene_index_tsv
        else repo / "OS_scRNA_gene_index.19264.tsv"
    )
    if not panel_path.exists():
        raise ScFMInputError(
            "[scfoundation] OS_scRNA_gene_index.19264.tsv not found "
            f"(looked at '{panel_path}'). It must be present in the official "
            "model folder (or pass --gene-index-tsv)."
        )

    # version='ce' (cell embedding) loads ./models/models.ckpt relative to the
    # script's cwd. Validate that official layout.
    expected_ckpt = repo / "models" / "models.ckpt"
    ckpt_for_meta = Path(model_path) if model_path else expected_ckpt
    if version == "ce" and not expected_ckpt.exists():
        raise ScFMDependencyError(
            "[scfoundation] official version='ce' cell embedding loads "
            f"'./models/models.ckpt' relative to the model folder, but "
            f"'{expected_ckpt}' is missing. Place (or symlink) the scFoundation "
            "checkpoint there on the run host."
        )

    _import_or_raise("anndata", "scfoundation", "pip install anndata")
    import anndata  # type: ignore
    import pandas as pd
    import scipy.sparse as sp

    input_gene_count = int(adata.n_vars)
    panel = list(pd.read_csv(panel_path, header=0, delimiter="\t")["gene_name"])
    if gene_name_col and gene_name_col in adata.var.columns:
        gene_names = list(adata.var[gene_name_col])
    else:
        gene_names = list(adata.var_names)
    matched_gene_count = _scfoundation_matched_gene_count(gene_names, panel)

    py = python_executable or sys.executable

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        in_h5ad = tmp / "scf_input.h5ad"
        out_dir = tmp / "out"
        out_dir.mkdir()

        # Hand the raw counts to the official script with gene names it reads via
        # var.gene_name (falling back to var_names). The script itself applies
        # main_gene_selection / normalization / read-depth tokens.
        X = adata.X
        var = pd.DataFrame(index=[str(g) for g in gene_names])
        var["gene_name"] = [str(g) for g in gene_names]
        tmp_adata = anndata.AnnData(
            X=(X.copy() if sp.issparse(X) else np.asarray(X)),
            obs=pd.DataFrame(index=[str(o) for o in adata.obs_names]),
            var=var,
        )
        tmp_adata.write_h5ad(in_h5ad)

        cmd = [
            py, "get_embedding.py",
            "--task_name", task_name,
            "--input_type", "singlecell",
            "--output_type", "cell",
            "--pool_type", pool_type,
            "--tgthighres", tgthighres,
            "--pre_normalized", pre_normalized,
            "--version", version,
            "--data_path", str(in_h5ad),
            "--save_path", str(out_dir),
            "--ckpt_name", ckpt_name,
        ]
        try:
            subprocess.run(
                cmd, cwd=str(repo), check=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
            )
        except FileNotFoundError as exc:  # pragma: no cover - Shirokane only
            raise ScFMDependencyError(
                f"[scfoundation] could not launch the official get_embedding.py "
                f"via '{py}'. Run this extractor with the python interpreter of "
                f"the scFoundation environment. Underlying error: {exc}"
            ) from exc
        except subprocess.CalledProcessError as exc:  # pragma: no cover
            tail = (exc.stdout or "")[-2000:]
            raise ScFMDependencyError(
                "[scfoundation] official get_embedding.py failed "
                f"(exit {exc.returncode}). Last output:\n{tail}"
            ) from exc

        out_npy = out_dir / _scfoundation_output_filename(
            task_name, ckpt_name, tgthighres
        )
        if not out_npy.exists():
            raise ScFMDependencyError(
                f"[scfoundation] expected official output '{out_npy.name}' was "
                f"not produced in '{out_dir}'. Contents: "
                f"{[p.name for p in out_dir.glob('*')]}"
            )
        embeddings = np.load(out_npy)

    embeddings = np.asarray(embeddings, dtype=np.float32)
    if embeddings.ndim == 1:
        embeddings = embeddings.reshape(1, -1)
    if embeddings.shape[0] != int(adata.n_obs):
        raise ScFMDependencyError(
            f"[scfoundation] official output had {embeddings.shape[0]} rows for "
            f"{adata.n_obs} cells; cell-order/count alignment failed."
        )

    run_meta = {
        "model_checkpoint": str(ckpt_for_meta),
        "upstream_commit_or_release": upstream_commit_or_release,
        "input_gene_count": input_gene_count,
        "matched_gene_count": matched_gene_count,
        "embedding_dim": int(embeddings.shape[1]),
        "obsm_key": ref["obsm_key"],
        "extra": {
            "invocation": "official get_embedding.py (subprocess)",
            "pool_type": pool_type,
            "tgthighres": tgthighres,
            "pre_normalized": pre_normalized,
            "version": version,
            "task_name": task_name,
            "ckpt_name": ckpt_name,
            "gene_panel_size": 19264,
            "matched_gene_count_basis": "input_genes_in_19264_panel",
        },
    }
    return embeddings, run_meta


EXTRACTORS = {
    "geneformer": extract_geneformer,
    "scgpt": extract_scgpt,
    "scfoundation": extract_scfoundation,
}
