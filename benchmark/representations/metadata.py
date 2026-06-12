"""
metadata.py
===========

Schema + validation for ``adata.uns["scfm_embedding_metadata"]`` and the
``representation_metadata.json`` artifacts described in the work plan
(Step 0 / Step 1).

Keeping the schema in one place lets the extractor, the representation builder,
the evaluators, and the tests all agree on which provenance fields are required.
"""

from __future__ import annotations

from datetime import datetime, timezone

# Code version for this representation-dynamics module. Bump on behavioural
# changes to extraction / reduction so artifacts remain traceable.
CODE_VERSION = "representation_dynamics_v1"

#: Fields that MUST be present (and non-empty) in scfm_embedding_metadata.
#: This list matches the work-plan requirements (Step 1 + the implementation
#: goals in the task brief).
REQUIRED_EMBEDDING_METADATA_FIELDS = (
    "model_name",
    "model_checkpoint",
    "paper_reference",
    "code_reference",
    "upstream_commit_or_release",
    "tokenizer_or_gene_vocab",
    "input_gene_count",
    "matched_gene_count",
    "preprocessing",
    "embedding_layer",
    "embedding_dim",
    "pooling_or_cell_embedding_method",
    "date_created",
    "code_version",
)

#: Fields whose values must not be None/empty (subset of the above that can
#: never be legitimately blank). ``upstream_commit_or_release`` is allowed to be
#: the string "unknown" but must still be present.
_NONEMPTY_FIELDS = (
    "model_name",
    "model_checkpoint",
    "paper_reference",
    "code_reference",
    "tokenizer_or_gene_vocab",
    "preprocessing",
    "embedding_layer",
    "pooling_or_cell_embedding_method",
    "date_created",
    "code_version",
)


def utc_now_iso() -> str:
    """Return current UTC time as an ISO-8601 string (date_created)."""
    return datetime.now(timezone.utc).isoformat()


class MetadataValidationError(ValueError):
    """Raised when an embedding metadata block is missing required provenance."""


def validate_embedding_metadata(meta: dict) -> dict:
    """
    Validate a single scFM embedding metadata block.

    Parameters
    ----------
    meta : dict
        Metadata block (one entry from ``adata.uns['scfm_embedding_metadata']``
        or a standalone ``representation_metadata.json`` model record).

    Returns
    -------
    dict
        The same ``meta`` (unchanged) if valid.

    Raises
    ------
    MetadataValidationError
        If any required field is missing, or any non-empty field is blank, or
        integer count fields are negative / non-integer.
    """
    if not isinstance(meta, dict):
        raise MetadataValidationError(
            f"Embedding metadata must be a dict, got {type(meta).__name__}."
        )

    missing = [f for f in REQUIRED_EMBEDDING_METADATA_FIELDS if f not in meta]
    if missing:
        raise MetadataValidationError(
            f"Embedding metadata is missing required fields: {missing}. "
            f"Present: {sorted(meta.keys())}."
        )

    blank = [
        f for f in _NONEMPTY_FIELDS
        if meta.get(f) is None or str(meta.get(f)).strip() == ""
    ]
    if blank:
        raise MetadataValidationError(
            f"Embedding metadata fields must not be blank: {blank}."
        )

    for count_field in ("input_gene_count", "matched_gene_count", "embedding_dim"):
        val = meta.get(count_field)
        if not isinstance(val, (int,)) or isinstance(val, bool) or val < 0:
            raise MetadataValidationError(
                f"{count_field} must be a non-negative integer, got {val!r}."
            )

    if meta["matched_gene_count"] > meta["input_gene_count"]:
        raise MetadataValidationError(
            "matched_gene_count "
            f"({meta['matched_gene_count']}) cannot exceed input_gene_count "
            f"({meta['input_gene_count']})."
        )

    return meta


def build_embedding_metadata(
    *,
    reference: dict,
    model_checkpoint: str,
    upstream_commit_or_release: str,
    input_gene_count: int,
    matched_gene_count: int,
    embedding_dim: int,
    extra: dict | None = None,
) -> dict:
    """
    Assemble a fully-populated, schema-valid embedding metadata block.

    ``reference`` is one of ``benchmark.representations.references.MODEL_REFERENCES``
    values (or a copy from :func:`references.get_reference`). Run-specific fields
    (checkpoint, matched gene counts, embedding dim, upstream commit) are supplied
    by the caller because they are only known at extraction time.
    """
    meta = {
        "model_name": reference["model_name"],
        "model_checkpoint": model_checkpoint,
        "paper_reference": reference["paper_reference"],
        "code_reference": reference["code_reference"],
        "upstream_commit_or_release": upstream_commit_or_release or "unknown",
        "tokenizer_or_gene_vocab": reference["tokenizer_or_gene_vocab"],
        # alias kept for the work-plan's "gene_vocabulary" wording
        "gene_vocabulary": reference["tokenizer_or_gene_vocab"],
        "input_gene_count": int(input_gene_count),
        "matched_gene_count": int(matched_gene_count),
        "preprocessing": reference["preprocessing"],
        "embedding_layer": reference["embedding_layer"],
        "embedding_dim": int(embedding_dim),
        "pooling_or_cell_embedding_method": reference[
            "pooling_or_cell_embedding_method"
        ],
        "source_inspected": reference.get("source_inspected", ""),
        "date_created": utc_now_iso(),
        "code_version": CODE_VERSION,
    }
    if extra:
        meta.update(extra)
    return validate_embedding_metadata(meta)
