"""Generate tiny synthetic AnnData inputs for the local smoke test.

Creates a labelled OSKM-like reference and two label-free chemical-like queries
that share a gene panel (incl. canonical markers), with a planted trajectory so
the annotation-comparison pipeline produces meaningful, non-trivial numbers
WITHOUT downloading any real data. This is test plumbing only — never a
manuscript input.

Usage:
    python scripts/make_synthetic_smoke_data.py --out-dir results/test_outputs/synthetic
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
import anndata as ad

MARKERS = ["POU5F1", "SOX2", "NANOG", "LIN28A",       # pluripotency
           "COL1A2", "DCN", "LUM",                    # stromal/fibroblast
           "CDH1", "EPCAM",                           # epithelial/MET
           "GATA6", "SOX17"]                          # XEN
FILLER = [f"G{i:03d}" for i in range(49)]
GENES = MARKERS + FILLER
G = len(GENES)
IDX = {g: i for i, g in enumerate(GENES)}

# filler structure: some track pseudotime, some mark off-target program
T_FILLER = FILLER[:12]
OFF_FILLER = FILLER[12:18]


def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def _program_means(t, is_xen, is_off):
    """Return per-gene Poisson means (length G) for a cell at pseudotime t."""
    lam = np.full(G, 0.4)
    pluri = _sigmoid((t - 0.7) / 0.1) * 6.0
    strom = _sigmoid((0.3 - t) / 0.1) * 6.0
    epi = np.exp(-((t - 0.4) / 0.12) ** 2) * 4.0
    for g in ["POU5F1", "SOX2", "NANOG", "LIN28A"]:
        lam[IDX[g]] += pluri
    for g in ["COL1A2", "DCN", "LUM"]:
        lam[IDX[g]] += strom
    for g in ["CDH1", "EPCAM"]:
        lam[IDX[g]] += epi
    xen = 5.0 if is_xen else 0.2
    for g in ["GATA6", "SOX17"]:
        lam[IDX[g]] += xen
    for g in T_FILLER:
        lam[IDX[g]] += t * 3.0
    if is_off:
        for g in OFF_FILLER:
            lam[IDX[g]] += 5.0
    return lam


def _simulate(n, label_fn, time_fn, seed, xen_frac=0.0, off_frac=0.0, sample_prefix="S"):
    rng = np.random.default_rng(seed)
    t = np.clip(rng.beta(1.4, 1.4, size=n), 0, 1)
    t.sort()
    is_xen = np.zeros(n, dtype=bool)
    is_off = np.zeros(n, dtype=bool)
    if xen_frac > 0:  # a chemical-specific XEN side branch in the mid-trajectory
        mid = np.where((t > 0.4) & (t < 0.7))[0]
        pick = rng.choice(mid, size=int(len(mid) * xen_frac), replace=False) if len(mid) else []
        is_xen[pick] = True
    if off_frac > 0:  # OSKM off-target (hOSK/xOSK) scattered cells
        pick = rng.choice(n, size=int(n * off_frac), replace=False)
        is_off[pick] = True

    X = np.zeros((n, G), dtype=np.int32)
    size = rng.gamma(shape=8.0, scale=1.0 / 8.0, size=n)  # library variation ~1
    for i in range(n):
        lam = _program_means(t[i], is_xen[i], is_off[i]) * size[i]
        X[i] = rng.poisson(lam)

    labels = np.array([label_fn(ti, xi, oi) for ti, xi, oi in zip(t, is_xen, is_off)])
    times = np.array([time_fn(ti) for ti in t])
    obs = pd.DataFrame({
        "latent_t": t,
        "sample_id": [f"{sample_prefix}_{tm}" for tm in times],
        "time_label": times,
        "pct_counts_mt": rng.uniform(0.5, 8.0, size=n),
        "n_genes_by_counts": (X > 0).sum(1),
        "total_counts": X.sum(1),
    }, index=[f"{sample_prefix}cell{i:05d}" for i in range(n)])
    A = ad.AnnData(X=X.astype(np.float32), obs=obs,
                   var=pd.DataFrame(index=GENES))
    return A, labels


def _time_bin(t):
    days = ["D0", "D2", "D4", "D6", "D8", "D10", "D12", "D14", "iPSC"]
    return days[min(int(t * 9), 8)]


def _ref_label(t, is_xen, is_off):
    if is_off:
        return "hOSK" if t < 0.6 else "xOSK"
    if t < 0.2:
        return "Fibroblast"
    if t < 0.4:
        return "Keratinocyte-like"
    if t < 0.55:
        return "Intermediate"
    if t < 0.7:
        return "Partially-reprogrammed"
    if t < 0.85:
        return "Pre-iPSC"
    return "iPSC"


def _query_label(t, is_xen, is_off):
    if is_xen:
        return "xen_like"
    if t < 0.2:
        return "hADSCs"
    if t < 0.45:
        return "epithelial_like"
    if t < 0.75:
        return "intermediate_plastic"
    if t < 0.82:
        return "ambiguous"
    return "hCiPS"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="results/test_outputs/synthetic")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    ref, ref_labels = _simulate(1500, _ref_label, _time_bin, args.seed,
                                off_frac=0.06, sample_prefix="REF")
    ref.obs["final_milestone_label_coarse"] = ref_labels
    ref.obs["author_cluster_label"] = ref_labels
    ref.write_h5ad(os.path.join(args.out_dir, "reference.h5ad"))

    # GSE178325: includes a chemical-specific xen_like branch
    q1, q1_labels = _simulate(1200, _query_label, _time_bin, args.seed + 1,
                              xen_frac=0.25, sample_prefix="G178")
    q1.obs["final_milestone_label_coarse"] = q1_labels
    q1.write_h5ad(os.path.join(args.out_dir, "gse178325.h5ad"))

    # GSE230659: cleaner near-linear trajectory, no xen branch
    q2, q2_labels = _simulate(1100, _query_label, _time_bin, args.seed + 2,
                              xen_frac=0.0, sample_prefix="G230")
    q2.obs["final_milestone_label_coarse"] = q2_labels
    q2.write_h5ad(os.path.join(args.out_dir, "gse230659.h5ad"))

    print(f"Wrote synthetic smoke data to {args.out_dir}")
    for name, A, L in [("reference", ref, ref_labels), ("gse178325", q1, q1_labels),
                       ("gse230659", q2, q2_labels)]:
        uniq, cnt = np.unique(L, return_counts=True)
        print(f"  {name}: {A.n_obs} cells x {A.n_vars} genes | "
              f"labels: {dict(zip(uniq, cnt))}")


if __name__ == "__main__":
    main()
