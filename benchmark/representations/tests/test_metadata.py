"""Tests for scFM embedding metadata schema/validation and the reference registry."""

import pytest

from benchmark.representations.metadata import (
    MetadataValidationError,
    REQUIRED_EMBEDDING_METADATA_FIELDS,
    build_embedding_metadata,
    validate_embedding_metadata,
)
from benchmark.representations.references import (
    MODEL_REFERENCES,
    REP_TO_RAW_OBSM,
    SUPPORTED_MODELS,
    get_reference,
)


def _valid_meta():
    return build_embedding_metadata(
        reference=get_reference("scgpt"),
        model_checkpoint="/ckpts/scGPT_human",
        upstream_commit_or_release="f609711",
        input_gene_count=2000,
        matched_gene_count=1800,
        embedding_dim=512,
    )


def test_build_embedding_metadata_has_all_required_fields():
    meta = _valid_meta()
    for field in REQUIRED_EMBEDDING_METADATA_FIELDS:
        assert field in meta, f"missing {field}"
    assert validate_embedding_metadata(meta) is meta


def test_validate_rejects_missing_field():
    meta = _valid_meta()
    del meta["paper_reference"]
    with pytest.raises(MetadataValidationError, match="missing required fields"):
        validate_embedding_metadata(meta)


def test_validate_rejects_blank_field():
    meta = _valid_meta()
    meta["preprocessing"] = "   "
    with pytest.raises(MetadataValidationError, match="must not be blank"):
        validate_embedding_metadata(meta)


def test_validate_rejects_matched_exceeding_input():
    meta = _valid_meta()
    meta["matched_gene_count"] = meta["input_gene_count"] + 1
    with pytest.raises(MetadataValidationError, match="cannot exceed"):
        validate_embedding_metadata(meta)


def test_validate_rejects_negative_counts():
    meta = _valid_meta()
    meta["embedding_dim"] = -3
    with pytest.raises(MetadataValidationError, match="non-negative integer"):
        validate_embedding_metadata(meta)


def test_references_present_for_all_supported_models():
    assert set(SUPPORTED_MODELS) == {"geneformer", "scgpt", "scfoundation"}
    for model in SUPPORTED_MODELS:
        ref = get_reference(model)
        # reference-first: paper + code + source inspected must be recorded
        assert ref["paper_reference"].startswith("https://www.nature.com/")
        assert ref["code_reference"].startswith("http")
        assert ref["source_inspected"]
        assert ref["obsm_key"] == REP_TO_RAW_OBSM[
            {"geneformer": "rep_geneformer_cls_pca50",
             "scgpt": "rep_scgpt_cls_pca50",
             "scfoundation": "rep_scfoundation_pca50"}[model]
        ]


def test_get_reference_unknown_model_raises():
    with pytest.raises(ValueError, match="Unknown scFM model"):
        get_reference("uce")


def test_reference_block_is_copy_not_alias():
    a = get_reference("geneformer")
    a["paper_reference"] = "mutated"
    assert MODEL_REFERENCES["geneformer"]["paper_reference"] != "mutated"
