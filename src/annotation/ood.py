"""Geometric out-of-distribution (OOD) scoring for reference mapping.

The first scANVI run showed that the classifier posterior is *miscalibrated* for
cross-protocol query cells: it confidently assigns out-of-distribution cells to a
"sink" reference label (here `hOSK`). A geometric OOD score fixes this without
trusting the classifier — each query cell's distance to its k nearest REFERENCE
cells in the shared latent space, calibrated against the reference's own
kNN-distance distribution. Cells far from the entire reference manifold are
flagged OOD regardless of how confident the classifier is.

Dependency-light (numpy + scikit-learn) so it is unit-testable and can post-hoc
re-score any saved latent (scANVI latent on the server; PCA latent locally).
"""
from __future__ import annotations

import numpy as np


def knn_reference_distance(ref_latent, query_latent, k: int = 15):
    """Mean distance to the k nearest reference cells, for query and (self) reference."""
    from sklearn.neighbors import NearestNeighbors

    k = int(min(k, ref_latent.shape[0] - 1))
    nn = NearestNeighbors(n_neighbors=k).fit(ref_latent)
    d_q, _ = nn.kneighbors(query_latent)
    d_r, _ = NearestNeighbors(n_neighbors=k + 1).fit(ref_latent).kneighbors(ref_latent)
    return d_q.mean(1), d_r[:, 1:].mean(1)  # drop self-neighbour for reference


def ood_score_and_flag(ref_latent, query_latent, k: int = 15, quantile: float = 0.95):
    """Return per-query OOD score (ratio to reference median) + flag (above the
    reference's `quantile` kNN-distance) + summary."""
    dq, dr = knn_reference_distance(ref_latent, query_latent, k)
    thr = float(np.quantile(dr, quantile))
    med = float(np.median(dr))
    return {
        "ood_score": dq / (med + 1e-12),
        "ood_flag": dq > thr,
        "threshold": thr,
        "ref_median_knn_dist": med,
        "frac_ood": float((dq > thr).mean()),
        "query_median_knn_dist": float(np.median(dq)),
    }


def ood_gate(ref_latent, query_latent, k: int = 15, quantile: float = 0.95,
             state_labels=None, reject_enrichment: float = 3.0):
    """Pipeline-facing wrapper around :func:`ood_score_and_flag`.

    Turns the per-cell geometric OOD score into a transfer-level *gate decision*
    so the same routine can be called from the main pipeline (PCA latent) and
    from ``scripts/score_ood.py`` (server scANVI latent).

    The reference self-flags exactly ``(1 - quantile)`` of its own cells by
    construction, so the decision-relevant statistic is the ENRICHMENT
    ``frac_ood / (1 - quantile)`` --- how many times more often a query cell is
    flagged than a reference cell. The transfer is rejected when the query is
    flagged more than ``reject_enrichment`` times the reference baseline. This is
    insensitive to ``k`` and to ``quantile`` over a sensible range (see
    ``scripts/ood_sensitivity_sweep.py``).

    If ``state_labels`` is given, a per-state OOD fraction is also returned.
    """
    res = ood_score_and_flag(ref_latent, query_latent, k=k, quantile=quantile)
    enrichment = float(res["frac_ood"] / (1.0 - quantile))
    out = {
        "frac_ood_geometric": float(res["frac_ood"]),
        "ood_threshold": float(res["threshold"]),
        "ref_median_knn_dist": float(res["ref_median_knn_dist"]),
        "query_median_knn_dist": float(res["query_median_knn_dist"]),
        "k": int(k),
        "quantile": float(quantile),
        "enrichment_vs_ref": enrichment,
        "gate_decision": "reject" if enrichment > reject_enrichment else "accept",
    }
    if state_labels is not None:
        sl = np.asarray(state_labels).astype(str)
        flag = np.asarray(res["ood_flag"])
        out["ood_by_state"] = {
            s: float(flag[sl == s].mean()) for s in sorted(set(sl.tolist()))
        }
    return out
