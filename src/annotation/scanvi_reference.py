"""Build the scANVI reference AnnData from GSE242424.

Production path (server). Harmonizes the reference gene space to the query
intersection and writes a frozen gene list so scArches surgery can reuse it.
"""
from __future__ import annotations

from typing import List, Sequence

from ..data.load import load_anndata, intersect_genes


def build_reference(
    reference_path: str,
    query_paths: Sequence[str],
    label_key: str,
    out_gene_list: str | None = None,
):
    """Load reference + queries, restrict to common genes (reference order),
    and stamp ``ref_label``. Returns (reference_adata, query_adatas, genes)."""
    ref = load_anndata(reference_path)
    queries = [load_anndata(p) for p in query_paths]
    if label_key not in ref.obs:
        raise KeyError(f"reference label_key '{label_key}' not in obs")
    genes, harmonized = intersect_genes([ref, *queries])
    ref_h, query_hs = harmonized[0], harmonized[1:]
    ref_h.obs["ref_label"] = ref.obs[label_key].astype(str).values

    if out_gene_list:
        from ..utils.io import ensure_dir
        import os

        ensure_dir(os.path.dirname(out_gene_list) or ".")
        with open(out_gene_list, "w", encoding="utf-8") as fh:
            fh.write("\n".join(genes))
    return ref_h, query_hs, genes
