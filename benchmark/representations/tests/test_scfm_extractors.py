"""
Tests for the reference-first scFM extractor adapters.

These tests must NOT require real checkpoints, GPUs, or the official scFM
packages. They exercise the pure-Python helpers and assert structural
reference-first guarantees by inspecting the adapter source.
"""

import inspect
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import benchmark.representations.scfm_extractors as ex
from benchmark.representations.metadata import build_embedding_metadata, validate_embedding_metadata
from benchmark.representations.references import get_reference

SOURCE = Path(ex.__file__).read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Minimal AnnData-like stub (no anndata dependency)
# ---------------------------------------------------------------------------

class _FakeAdata:
    def __init__(self, var, obs):
        self.var = var
        self.obs = obs

    @property
    def n_obs(self):
        return self.obs.shape[0]

    @property
    def n_vars(self):
        return self.var.shape[0]

    def copy(self):
        return _FakeAdata(self.var.copy(), self.obs.copy())


# ===========================================================================
# Geneformer: custom columns copied into official fixed names
# ===========================================================================

def test_geneformer_copies_custom_cols_into_official_names():
    var = pd.DataFrame({"my_ensembl": ["ENSG1", "ENSG2", "ENSG3"]},
                       index=["g1", "g2", "g3"])
    obs = pd.DataFrame({"my_counts": [10.0, 20.0]}, index=["c1", "c2"])
    adata = _FakeAdata(var, obs)

    prepared = ex._geneformer_prepare_official_adata(
        adata, ensembl_col="my_ensembl", counts_col="my_counts",
        index_col="__rep_cell_index",
    )
    # Official fixed names now exist and carry the custom column values.
    assert list(prepared.var["ensembl_id"]) == ["ENSG1", "ENSG2", "ENSG3"]
    assert list(prepared.obs["n_counts"]) == [10.0, 20.0]
    # Stable index for cell-order restoration.
    assert list(prepared.obs["__rep_cell_index"]) == [0, 1]
    # Original object untouched (we copied).
    assert "ensembl_id" not in adata.var.columns
    assert "n_counts" not in adata.obs.columns


def test_geneformer_official_names_passthrough_when_already_present():
    var = pd.DataFrame({"ensembl_id": ["ENSG1", "ENSG2"]}, index=["g1", "g2"])
    obs = pd.DataFrame({"n_counts": [5.0, 6.0]}, index=["c1", "c2"])
    adata = _FakeAdata(var, obs)
    prepared = ex._geneformer_prepare_official_adata(
        adata, "ensembl_id", "n_counts", "__rep_cell_index"
    )
    assert list(prepared.var["ensembl_id"]) == ["ENSG1", "ENSG2"]
    assert list(prepared.obs["n_counts"]) == [5.0, 6.0]


# ===========================================================================
# Geneformer: matched_gene_count is tokenizer-consistent, not == input
# ===========================================================================

def test_geneformer_matched_gene_count_uses_token_dictionary():
    token_dict = {"ENSG1": 5, "ENSG2": 6, "<pad>": 0, "<cls>": 1}
    ensembl_ids = ["ENSG1", "ENSG2", "ENSG_MISSING", "ENSG_ALSO_MISSING"]
    matched = ex._geneformer_matched_gene_count(ensembl_ids, token_dict)
    assert matched == 2
    # Crucially NOT blindly equal to the input gene count.
    assert matched != len(ensembl_ids)


def test_geneformer_source_does_not_hardcode_matched_equals_input():
    # The old bug set matched_gene_count = input_gene_count. Ensure it's gone.
    assert "matched_gene_count\": input_gene_count" not in SOURCE
    assert "_geneformer_matched_gene_count(" in SOURCE


def test_geneformer_cls_requires_v2():
    var = pd.DataFrame({"ensembl_id": ["ENSG1"]}, index=["g1"])
    obs = pd.DataFrame({"n_counts": [1.0]}, index=["c1"])
    adata = _FakeAdata(var, obs)
    with pytest.raises(ex.ScFMInputError, match="emb_mode='cls' requires a V2"):
        ex.extract_geneformer(adata, model_dir="/nonexistent",
                              emb_mode="cls", model_version="V1")


def test_geneformer_passes_model_version_to_tokenizer():
    # Tokenizer must receive model_version consistently with EmbExtractor.
    assert 'tk_kwargs = {"model_version": model_version}' in SOURCE
    assert "TranscriptomeTokenizer(" in SOURCE


# ===========================================================================
# scFoundation: no load.py main_gene_selection import, no manual forward
# ===========================================================================

def test_scfoundation_does_not_import_main_gene_selection_from_load():
    assert "main_gene_selection" not in SOURCE or (
        # allowed only as a comment reference, never as an import
        "from load import" not in SOURCE
    )
    assert "from load import" not in SOURCE


def test_scfoundation_has_no_manual_encoder_or_pooling():
    for forbidden in (
        "pretrainmodel.encoder",
        "pretrainmodel.token_emb",
        "load_model_frommmf(",
        "gatherData(",
        "torch.concat([geneemb",
    ):
        assert forbidden not in SOURCE, (
            f"scFoundation adapter must not hand-roll official logic: {forbidden!r}"
        )


def test_scfoundation_invokes_official_get_embedding_subprocess():
    assert "get_embedding.py" in SOURCE
    assert "subprocess.run(" in SOURCE


def test_scfoundation_output_filename_matches_official_pattern():
    name = ex._scfoundation_output_filename("rep", "models", "t4")
    assert name == "rep_models_singlecell_cell_embedding_t4_resolution.npy"


def test_scfoundation_matched_gene_count_panel_intersection():
    panel = ["A", "B", "C", "D"]
    genes = ["A", "B", "Z"]
    assert ex._scfoundation_matched_gene_count(genes, panel) == 2


def test_scfoundation_requires_repo():
    adata = _FakeAdata(pd.DataFrame(index=["g1"]), pd.DataFrame(index=["c1"]))
    with pytest.raises((ex.ScFMInputError, ex.ScFMDependencyError)):
        ex.extract_scfoundation(adata, scfoundation_repo="/definitely/not/here")


# ===========================================================================
# scGPT: still delegates to official embed_data; override recorded
# ===========================================================================

def test_scgpt_delegates_to_official_embed_data():
    assert "from scgpt.tasks.cell_emb import embed_data" in SOURCE
    assert "embed_data(" in SOURCE
    # Must not reimplement binning / forward / cls extraction locally.
    for forbidden in ("do_binning", "DataCollator", "_encode(", "[:, 0, :]"):
        assert forbidden not in SOURCE


def test_scgpt_records_use_fast_transformer_override():
    extra = ex._scgpt_run_extra(use_fast_transformer=False, max_length=1200,
                                gene_col="feature_name")
    assert extra["use_fast_transformer"] is False
    assert extra["official_default_use_fast_transformer"] is True
    assert "use_fast_transformer_override" in extra  # override recorded


def test_scgpt_no_override_recorded_when_matching_official_default():
    extra = ex._scgpt_run_extra(use_fast_transformer=True, max_length=1200,
                                gene_col="index")
    assert extra["use_fast_transformer"] is True
    assert "use_fast_transformer_override" not in extra


# ===========================================================================
# Metadata stays strict for every model's reference + a synthetic run_meta
# ===========================================================================

@pytest.mark.parametrize("model,emb_dim", [
    ("geneformer", 512), ("scgpt", 512), ("scfoundation", 3072),
])
def test_metadata_builds_and_validates_for_each_model(model, emb_dim):
    meta = build_embedding_metadata(
        reference=get_reference(model),
        model_checkpoint="/ckpts/model",
        upstream_commit_or_release="abc1234",
        input_gene_count=2000,
        matched_gene_count=1500,
        embedding_dim=emb_dim,
    )
    validate_embedding_metadata(meta)
    # Reference-first provenance must be present.
    for field in ("paper_reference", "code_reference", "source_inspected",
                  "tokenizer_or_gene_vocab", "preprocessing", "embedding_layer",
                  "pooling_or_cell_embedding_method"):
        assert meta.get(field)
