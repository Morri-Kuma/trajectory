"""
references.py
=============

Canonical paper / code / checkpoint references for each single-cell foundation
model used in the representation-dynamics benchmark.

This module is the single source of truth for the *reference-first* metadata
required by the work plan (Step 0). The values here were transcribed directly
from the official paper code-availability statements and the official
implementation entry points (see ``source_inspected`` for the exact file/class
that each adapter mirrors). They are embedded in every emitted embedding
artifact via ``adata.uns["scfm_embedding_metadata"]`` so that any downstream
result is fully traceable.

If the reference paper and the upstream repository disagree, STOP and record the
discrepancy before running benchmarks (work plan Step 0, rule 3).
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Reference registry
# ---------------------------------------------------------------------------
# Each entry documents:
#   paper_reference            : DOI/URL of the peer-reviewed paper
#   code_reference             : official implementation repository
#   default_checkpoint         : the checkpoint the adapter expects by default
#   tokenizer_or_gene_vocab    : the tokenizer/vocabulary the official code uses
#   preprocessing              : preprocessing performed by the official code
#   embedding_layer            : which layer the cell embedding is read from
#   pooling_or_cell_embedding_method : how the per-cell vector is pooled
#   source_inspected           : the exact official file/function this adapter
#                                mirrors (cited so a reviewer can diff it)
#   obsm_key                   : adata.obsm key the extractor writes
#   notes                      : anything a reviewer must know

MODEL_REFERENCES = {
    "geneformer": {
        "model_name": "Geneformer",
        "paper_reference": "https://www.nature.com/articles/s41586-023-06139-9",
        "code_reference": "https://huggingface.co/ctheodoris/Geneformer",
        "default_checkpoint": "ctheodoris/Geneformer (Geneformer-V2-104M)",
        "tokenizer_or_gene_vocab": (
            "geneformer.TranscriptomeTokenizer rank-value encoding; "
            "token_dictionary.pkl (Ensembl ID -> token), "
            "gene_median_dictionary.pkl"
        ),
        "preprocessing": (
            "Rank-value encoding per cell: counts normalized by each gene's "
            "non-zero median expression across Genecorpus, ranked, truncated to "
            "model input size. Requires raw counts in adata.X, Ensembl IDs in "
            "adata.var['ensembl_id'], and per-cell total counts in "
            "adata.obs['n_counts']. Performed by the official TranscriptomeTokenizer."
        ),
        "embedding_layer": "-1 (2nd-to-last hidden layer; EmbExtractor default)",
        "pooling_or_cell_embedding_method": (
            "CLS-token embedding (emb_mode='cls'): cell embedding read from the "
            "<cls> token prepended to the rank-value encoding"
        ),
        "source_inspected": (
            "geneformer.EmbExtractor (emb_mode='cls', cell_emb_style='mean_pool', "
            "emb_layer=-1, model_version='V2').extract_embs(...) and "
            "geneformer.TranscriptomeTokenizer; "
            "https://geneformer.readthedocs.io/en/latest/geneformer.emb_extractor.html"
        ),
        "obsm_key": "X_geneformer_cls",
        "notes": (
            "emb_mode='cls' requires a V2 (CLS-token) checkpoint. For a V1 "
            "checkpoint without a CLS token, the official fallback is "
            "emb_mode='cell' (mean-pool of gene embeddings); this must be set "
            "explicitly and is recorded in the metadata."
        ),
    },
    "scgpt": {
        "model_name": "scGPT",
        "paper_reference": "https://www.nature.com/articles/s41592-024-02201-0",
        "code_reference": "https://github.com/bowang-lab/scGPT",
        "default_checkpoint": "scGPT whole-human (best_model.pt + vocab.json + args.json)",
        "tokenizer_or_gene_vocab": (
            "scgpt.tokenizer.GeneVocab loaded from the checkpoint's vocab.json; "
            "genes matched to var[gene_col]; unmatched genes (id_in_vocab < 0) dropped"
        ),
        "preprocessing": (
            "Value binning inside scgpt.DataCollator(do_binning=True) applied to "
            "non-zero gene expression after a <cls> token is prepended "
            "(value=pad_value). max_length default 1200, sampling=True, "
            "keep_first_n_tokens=1. Mirrors official embed_data()."
        ),
        "embedding_layer": "final transformer encoder output (model._encode)",
        "pooling_or_cell_embedding_method": (
            "CLS-token embedding: embeddings[:, 0, :] (the <cls> position), then "
            "L2-normalized per cell (official get_batch_cell_embeddings)"
        ),
        "source_inspected": (
            "scgpt.tasks.cell_emb.embed_data() and get_batch_cell_embeddings("
            "cell_embedding_mode='cls'); "
            "https://github.com/bowang-lab/scGPT/blob/main/scgpt/tasks/cell_emb.py"
        ),
        "obsm_key": "X_scgpt_cls",
        "notes": (
            "Official code writes adata.obsm['X_scGPT']; we re-key to "
            "'X_scgpt_cls' for consistency with the plan but preserve the exact "
            "official computation (CLS position + L2 normalization)."
        ),
    },
    "scfoundation": {
        "model_name": "scFoundation",
        "paper_reference": "https://www.nature.com/articles/s41592-024-02305-7",
        "code_reference": "https://github.com/biomap-research/scFoundation",
        "default_checkpoint": "scFoundation models.ckpt (key='cell', version='ce')",
        "tokenizer_or_gene_vocab": (
            "Fixed 19,264-gene panel from OS_scRNA_gene_index.19264.tsv; "
            "main_gene_selection() aligns/zero-pads input genes to this panel"
        ),
        "preprocessing": (
            "Single-cell: per-cell log1p(x / total * 1e4) then append two "
            "read-depth tokens [tgthighres, log10(total_count)] -> length 19266; "
            "data_gene_ids = arange(19266). Mirrors official get_embedding.py "
            "(input_type='singlecell', pre_normalized='F')."
        ),
        "embedding_layer": "xTrimoGene encoder output (pretrainmodel.encoder)",
        "pooling_or_cell_embedding_method": (
            "pool_type='all' (default): concat of [last token, 2nd-to-last token, "
            "max-pool over gene tokens, mean-pool over gene tokens] -> 4*embdim. "
            "pool_type='max': max over all tokens -> embdim."
        ),
        "source_inspected": (
            "scFoundation model/get_embedding.py main() invoked AS-IS via "
            "subprocess (output_type='cell', input_type='singlecell', "
            "pool_type='all', tgthighres='t4', pre_normalized='F', version='ce'); "
            "the script itself performs main_gene_selection (defined inside "
            "get_embedding.py), normalization, read-depth tokens, the xTrimoGene "
            "encoder forward and pooling. We do NOT re-implement that logic. "
            "https://github.com/biomap-research/scFoundation/blob/main/model/get_embedding.py "
            "(load.py: https://github.com/biomap-research/scFoundation/blob/main/model/load.py)"
        ),
        "obsm_key": "X_scfoundation",
        "notes": (
            "Thin subprocess wrapper around the official get_embedding.py. "
            "Requires the official scFoundation 'model' folder containing "
            "get_embedding.py, load.py, pretrainmodels.py, "
            "OS_scRNA_gene_index.19264.tsv, and models/models.ckpt (the official "
            "version='ce' cell path loads ./models/models.ckpt relative to the "
            "folder). Run the extractor with the scFoundation environment's "
            "python so the subprocess inherits torch + deps. tgthighres default "
            "'t4' per the official cell-embedding example. We do NOT import "
            "main_gene_selection from load.py (get_embedding.py defines its own)."
        ),
    },
}

#: scFM model keys supported by this module (and ONLY these).
SUPPORTED_MODELS = tuple(MODEL_REFERENCES.keys())

#: The four representation arms compared in this benchmark (and ONLY these).
REPRESENTATION_ARMS = (
    "rep_hvg_pca50",
    "rep_geneformer_cls_pca50",
    "rep_scgpt_cls_pca50",
    "rep_scfoundation_pca50",
)

#: Mapping representation_id -> the obsm key holding the *raw* scFM embedding.
#: rep_hvg_pca50 has no scFM source (it is built from raw counts).
REP_TO_RAW_OBSM = {
    "rep_hvg_pca50": None,
    "rep_geneformer_cls_pca50": "X_geneformer_cls",
    "rep_scgpt_cls_pca50": "X_scgpt_cls",
    "rep_scfoundation_pca50": "X_scfoundation",
}


def get_reference(model: str) -> dict:
    """Return a *copy* of the canonical reference block for ``model``."""
    key = model.lower()
    if key not in MODEL_REFERENCES:
        raise ValueError(
            f"Unknown scFM model {model!r}. Supported: {SUPPORTED_MODELS}."
        )
    return dict(MODEL_REFERENCES[key])
