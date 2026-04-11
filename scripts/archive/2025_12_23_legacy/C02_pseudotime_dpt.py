# C02_pseudotime_dpt.py
# For GSE178325 (human) chemical reprogramming
#
# Requested behavior:
# - ONLY keep donor=0618
# - ONLY keep stages: Somatic, StageI, StageII, StageIII, StageIV, hCiPSC
# - Exclude H1
# - Compute DPT (root=Somatic), but if stage medians are NOT monotonic in expected order,
#   automatically fall back to a stage-aware pseudotime ("stage_progress") to enforce:
#     Somatic -> StageI -> StageII -> StageIII -> StageIV -> hCiPSC
# - Keep DPT result in obs['pseudotime_dpt'] for reference.
# - Optional: isotonic calibration of DPT to progress (if sklearn available).
# - Timestamped outputs + sanity check prints.

import re
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime

from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree

# ---- optional isotonic regression ----
try:
    from sklearn.isotonic import IsotonicRegression
    HAS_SKLEARN = True
except Exception:
    HAS_SKLEARN = False

# =========================
# Config
# =========================
ROOT = Path(r"E:\scgpt")
BASE_RESULTS = ROOT / r"results\C_traj"

DATASET_TAG = "gse178325"

EXCLUDE_STAGES = {"H1"}                 # Plan A
DONOR_KEEP = "0618"                     # ONLY donor=0618
STAGES_KEEP = ["Somatic", "StageI", "StageII", "StageIII", "StageIV", "hCiPSC"]

MIN_CELLS_PER_CLUSTER = 150
N_PCS_FOR_MST = 20
SHOW_CLUSTER_LABELS = True

N_DCS = 10

STAGE_CMAP = {
    "Somatic": "Greys",
    "StageI": "Blues",
    "StageII": "Oranges",
    "StageIII": "Greens",
    "StageIV": "Purples",
    "hCiPSC": "Reds",
    "H1": "Greys",
}
CMAP_START = 0.25
CMAP_END = 0.95


# =========================
# Helpers
# =========================
def get_run_dir() -> Path:
    cand1 = BASE_RESULTS / DATASET_TAG / "_latest_run.txt"
    cand2 = BASE_RESULTS / "_latest_run.txt"

    if cand1.exists():
        latest = cand1
    elif cand2.exists():
        latest = cand2
    else:
        raise FileNotFoundError(
            f"Missing _latest_run.txt under {cand1.parent} or {cand2.parent}. Run C01 first."
        )

    run_dir = Path(latest.read_text(encoding="utf-8").strip())
    if not run_dir.exists():
        raise FileNotFoundError(f"Run dir not found: {run_dir}")
    return run_dir


def standardize_stage_name(x: str) -> str:
    s = str(x).strip()
    low = s.lower()

    if "hcipsc" in low:
        return "hCiPSC"
    if low in {"h1", "esc", "hesc"} or low.startswith("h1"):
        return "H1"
    if "somatic" in low:
        return "Somatic"

    # StageIV before StageI
    if "stageiv" in low or re.search(r"\bstage\s*4\b", low):
        return "StageIV"
    # StageIII before StageII
    if "stageiii" in low or re.search(r"\bstage\s*3\b", low):
        return "StageIII"
    if "stageii" in low or re.search(r"\bstage\s*2\b", low):
        return "StageII"
    if re.search(r"stagei(?!i)(?!v)", low) or re.search(r"\bstage\s*1\b", low):
        return "StageI"

    return s


def normalize_01(x: np.ndarray) -> np.ndarray:
    x = x.astype(float, copy=True)
    x[~np.isfinite(x)] = np.nan
    mn = np.nanmin(x)
    mx = np.nanmax(x)
    if not (np.isfinite(mn) and np.isfinite(mx) and mx > mn):
        return x
    return (x - mn) / (mx - mn)


def build_progress(adata: sc.AnnData) -> np.ndarray:
    rank = {
        "Somatic": 0,
        "StageI": 1,
        "StageII": 2,
        "StageIII": 3,
        "StageIV": 4,
        "hCiPSC": 5,
        "H1": 6,
    }
    st = adata.obs["stage"].astype(str).map(rank)
    day = pd.to_numeric(adata.obs.get("day", np.nan), errors="coerce")
    prog = st.to_numpy(dtype=float) * 1000.0 + day.to_numpy(dtype=float)
    prog[~np.isfinite(prog)] = np.nan
    return prog


def ensure_embeddings(adata: sc.AnnData):
    if "X_pca" not in adata.obsm:
        sc.tl.pca(adata, n_comps=50, svd_solver="randomized", zero_center=False)

    # always rebuild neighbors for filtered subset stability
    sc.pp.neighbors(adata, n_neighbors=25, n_pcs=50, use_rep="X_pca")

    conn = adata.obsp.get("connectivities", None)
    if conn is not None:
        rs = np.ravel(conn.sum(axis=1))
        n_iso = int((rs == 0).sum())
        print(f"[INFO] isolated cells after rebuild neighbors: {n_iso}")
        if n_iso > 0:
            sc.pp.neighbors(adata, n_neighbors=40, n_pcs=50, use_rep="X_pca")
            conn2 = adata.obsp["connectivities"]
            rs2 = np.ravel(conn2.sum(axis=1))
            n_iso2 = int((rs2 == 0).sum())
            print(f"[INFO] isolated cells after rebuild neighbors (n_neighbors=40): {n_iso2}")

    if "X_diffmap" not in adata.obsm:
        sc.tl.diffmap(adata)

# definite start
def pick_somatic_root_cell(adata: sc.AnnData):
    st = adata.obs["stage"].astype(str)
    m = (st == "Somatic").to_numpy()
    idxs = np.where(m)[0]
    if len(idxs) == 0:
        # fallback: smallest day
        day = pd.to_numeric(adata.obs.get("day", np.nan), errors="coerce").to_numpy()
        iroot = int(np.nanargmin(day)) if np.isfinite(day).any() else 0
        root_gsm = str(adata.obs["gsm"].iloc[iroot]) if "gsm" in adata.obs.columns else "NA"
        return iroot, root_gsm

    if "X_pca" in adata.obsm:
        X = adata.obsm["X_pca"][idxs, : min(50, adata.obsm["X_pca"].shape[1])]
        ctr = X.mean(axis=0, keepdims=True)
        d2 = ((X - ctr) ** 2).sum(axis=1)
        iroot = int(idxs[int(np.argmin(d2))])
    else:
        iroot = int(idxs[0])

    root_gsm = str(adata.obs["gsm"].iloc[iroot]) if "gsm" in adata.obs.columns else "NA"
    return iroot, root_gsm

# run DPT from Scanpy
def run_dpt(adata: sc.AnnData, iroot: int, n_dcs: int = 10) -> np.ndarray:
    adata.uns["iroot"] = int(iroot)
    sc.tl.dpt(adata, n_dcs=n_dcs)
    if "dpt_pseudotime" not in adata.obs:
        raise ValueError("scanpy did not produce obs['dpt_pseudotime'].")
    return pd.to_numeric(adata.obs["dpt_pseudotime"], errors="coerce").to_numpy(dtype=float)

# calculate the medians of Pseudotime by Stage
def stage_medians(pt: np.ndarray, stages: np.ndarray, stage_order: list[str]) -> dict:
    out = {}
    for st in stage_order:
        m = (stages == st)
        if m.any():
            out[st] = float(np.nanmedian(pt[m]))
        else:
            out[st] = np.nan
    return out

# Endpoint Sanity Check
def endpoint_ok(med: dict) -> bool:
    a = med.get("Somatic", np.nan)
    b = med.get("hCiPSC", np.nan)
    return bool(np.isfinite(a) and np.isfinite(b) and (b > a))

# Monotonicity Check
def monotonic_ok(med: dict, stage_order: list[str]) -> bool:
    prev = None
    for st in stage_order:
        v = med.get(st, np.nan)
        if not np.isfinite(v):
            continue
        if prev is not None and v <= prev:
            return False
        prev = v
    return True

# compare Pseudotime with experiment progress
def corr_progress(pt: np.ndarray, prog: np.ndarray) -> float:
    m = np.isfinite(pt) & np.isfinite(prog)
    if m.sum() < 20:
        return np.nan
    return float(np.corrcoef(prog[m], pt[m])[0, 1])

# Sort these GSMs in the order of the biological progress sequence.
def auto_build_order_gsms(adata: sc.AnnData):
    """
    Build ORDER_GSMS from obs by sorting:
      stage_rank (Somatic < StageI < StageII < StageIII < StageIV < hCiPSC < others)
      then (baseline-without-Day first)
      then day
      then gsm
    """
    if "gsm" not in adata.obs.columns:
        raise ValueError("obs['gsm'] missing. Please ensure B01 wrote obs['gsm'].")
    if "stage" not in adata.obs.columns:
        raise ValueError("obs['stage'] missing.")
    if "day" not in adata.obs.columns:
        raise ValueError("obs['day'] missing.")

    df = adata.obs[["gsm", "stage", "day"]].copy().dropna(subset=["gsm"])
    df["stage"] = df["stage"].astype(str).map(standardize_stage_name)
    df["day"] = pd.to_numeric(df["day"], errors="coerce")

    # --- NEW: baseline sample (no "Day" in sample_name) should come FIRST within stage ---
    # is_day_sample: 1 if contains "Day", else 0
    if "sample_name" in adata.obs.columns:
        sn = adata.obs.loc[df.index, "sample_name"].astype(str)
        df["is_day_sample"] = sn.str.contains(r"\bday\b", case=False, regex=True).astype(int)
    else:
        # if no sample_name info, fallback: treat all as day-samples (keep old behavior)
        df["is_day_sample"] = 1

    rank_map = {
        "Somatic": 0,
        "StageI": 1,
        "StageII": 2,
        "StageIII": 3,
        "StageIV": 4,
        "hCiPSC": 5,
        "H1": 6,
    }
    df["stage_rank"] = df["stage"].map(rank_map).fillna(99).astype(int)

    # one row per GSM
    df_u = (
        df.groupby("gsm", observed=True)
          .agg(
              stage=("stage", lambda x: str(x.value_counts().index[0])),
              day=("day", "median"),
              stage_rank=("stage_rank", "min"),
              # --- NEW: per-GSM is_day_sample (any cell has Day -> 1) ---
              is_day_sample=("is_day_sample", "max"),
          )
          .reset_index()
    )

    # --- NEW SORT: stage_rank -> is_day_sample(0 first) -> day -> gsm ---
    df_u = df_u.sort_values(["stage_rank", "is_day_sample", "day", "gsm"], ascending=True)

    order_gsms = df_u["gsm"].astype(str).tolist()

    print("[INFO] auto ORDER_GSMS built from obs (top 20 shown):")
    print(df_u.head(20).to_string(index=False))
    return order_gsms, df_u

# Generate a gradient color sequence at the same stage
def sample_cmap(cmap_name: str, n: int, start: float = CMAP_START, end: float = CMAP_END):
    cmap = plt.colormaps.get_cmap(cmap_name)
    if n <= 1:
        return [cmap((start + end) / 2)]
    xs = np.linspace(start, end, n)
    return [cmap(float(x)) for x in xs]

# Match the sample sequence generated by "auto_build_order_gsms" at the front with the gradient color generated by "sample_cmap" precisely.
def build_timepoint_labels_and_palette(adata: sc.AnnData, order_gsms: list[str]):
    present = set(adata.obs["gsm"].astype(str).unique().tolist())
    order_gsms = [g for g in order_gsms if g in present]
    if len(order_gsms) == 0:
        raise ValueError("No GSMs left after intersecting order_gsms with data.")

    gsm2label = {}
    if "sample_name" in adata.obs.columns:
        tmp = adata.obs[["gsm", "sample_name"]].dropna()
        for g, sub in tmp.groupby("gsm", observed=True):
            gsm2label[str(g)] = str(sub["sample_name"].iloc[0])
    else:
        for g in order_gsms:
            gsm2label[g] = g

    gsm2stage = {}
    tmp = adata.obs[["gsm", "stage"]].dropna()
    for g, sub in tmp.groupby("gsm", observed=True):
        st = str(sub["stage"].astype(str).value_counts().index[0])
        gsm2stage[str(g)] = standardize_stage_name(st)

    ordered_labels = [gsm2label.get(g, g) for g in order_gsms]

    stage2labels = {}
    for g, lab in zip(order_gsms, ordered_labels):
        st = gsm2stage.get(g, "Unknown")
        stage2labels.setdefault(st, []).append(lab)

    label2color = {}
    for st, labs in stage2labels.items():
        cmap_name = STAGE_CMAP.get(st, "Greys")
        cols = sample_cmap(cmap_name, len(labs))
        for lab, col in zip(labs, cols):
            label2color[lab] = col

    palette_list = [label2color.get(lab, (0.7, 0.7, 0.7, 1.0)) for lab in ordered_labels]

    mapper = {g: gsm2label.get(g, g) for g in order_gsms}
    tp = adata.obs["gsm"].astype(str).map(mapper)
    adata.obs["timepoint_gsm_ordered"] = pd.Categorical(tp, categories=ordered_labels, ordered=True)

    return ordered_labels, palette_list, order_gsms

#  calculate the geometric center (Centroids) of each sampling time point in the low-dimensional space of UMAP.
def compute_timepoint_centroids(adata: sc.AnnData, ordered_labels: list[str]):
    umap = adata.obsm["X_umap"]
    lab = adata.obs["timepoint_gsm_ordered"].astype(str)

    rows = []
    for i, l in enumerate(ordered_labels):
        m = (lab == str(l)).to_numpy()
        if m.sum() == 0:
            rows.append({"order": i, "label": l, "n": 0, "x": np.nan, "y": np.nan})
        else:
            xy = umap[m]
            rows.append({"order": i, "label": l, "n": int(m.sum()), "x": float(xy[:, 0].mean()), "y": float(xy[:, 1].mean())})
    return pd.DataFrame(rows)

# Cluster and calculate the average position and average development progress of each cluster.
def compute_cluster_stats(adata: sc.AnnData, cluster_key="leiden", pt_key="pseudotime"):
    umap = adata.obsm["X_umap"]
    Xpca = adata.obsm["X_pca"][:, :N_PCS_FOR_MST]
    cl = adata.obs[cluster_key].astype(str).to_numpy()
    pt = pd.to_numeric(adata.obs[pt_key], errors="coerce").to_numpy()

    ok = np.isfinite(pt)
    df = pd.DataFrame({"cluster": cl[ok], "pt": pt[ok], "x": umap[ok, 0], "y": umap[ok, 1]})

    cent = (
        df.groupby("cluster", observed=True)
        .agg(n=("pt", "size"), pt_median=("pt", "median"), x_umap=("x", "mean"), y_umap=("y", "mean"))
        .reset_index()
    )

    cent = cent[cent["n"] >= MIN_CELLS_PER_CLUSTER].copy()
    if cent.shape[0] < 2:
        raise ValueError(f"Too few clusters after filtering (MIN_CELLS_PER_CLUSTER={MIN_CELLS_PER_CLUSTER}). Lower it.")

    clusters_kept = cent["cluster"].astype(str).tolist()

    pca_centroids = []
    for c in clusters_kept:
        m = (cl[ok] == c)
        pca_centroids.append(Xpca[ok][m].mean(axis=0))
    pca_centroids = np.vstack(pca_centroids)

    sort_idx = np.argsort(cent["pt_median"].to_numpy(dtype=float))
    cent = cent.iloc[sort_idx].reset_index(drop=True)
    pca_centroids = pca_centroids[sort_idx]
    clusters_kept = [clusters_kept[i] for i in sort_idx]

    cent["leiden_time"] = [str(i) for i in range(cent.shape[0])]
    return cent, clusters_kept, pca_centroids

# Re-map the sorted "time labels" back onto each individual cell.
def assign_leiden_time(adata: sc.AnnData, cent: pd.DataFrame, cluster_key="leiden"):
    mapper = dict(zip(cent["cluster"].astype(str), cent["leiden_time"].astype(str)))
    new = adata.obs[cluster_key].astype(str).map(mapper)
    cats = list(cent["leiden_time"].astype(str))
    adata.obs["leiden_time"] = pd.Categorical(new, categories=cats, ordered=True)

# Connect the centroids of the clustering clusters with the shortest total distance, thereby deriving the evolutionary path of the cell population.
def mst_edges_from_pca(pca_centroids: np.ndarray):
    D = squareform(pdist(pca_centroids, metric="euclidean"))
    mst = minimum_spanning_tree(D).toarray()
    edges = set()
    n = mst.shape[0]
    for i in range(n):
        for j in range(n):
            if mst[i, j] > 0:
                a, b = (i, j) if i < j else (j, i)
                edges.add((a, b))
    return sorted(edges)

# Draw development trajectory arrows with directions on the UMAP plot
def draw_mst_arrows_on_umap(ax: plt.Axes, cent: pd.DataFrame, edges):
    coords = cent[["x_umap", "y_umap"]].to_numpy(dtype=float)
    ptm = cent["pt_median"].to_numpy(dtype=float)

    for a, b in edges:
        i, j = (a, b) if ptm[a] <= ptm[b] else (b, a)
        x1, y1 = coords[i]
        x2, y2 = coords[j]
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="->", lw=2.2, color="black", alpha=0.85))

    if SHOW_CLUSTER_LABELS:
        for _, r in cent.iterrows():
            ax.text(float(r["x_umap"]), float(r["y_umap"]), str(r["leiden_time"]),
                    fontsize=10, ha="center", va="center", color="black",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.65))


# =========================
# Sanity checks
# =========================
# Evaluates and prints a detailed quality report for a pseudotime candidate.
def print_candidate_report(name: str, pt: np.ndarray, stages: np.ndarray, prog: np.ndarray, stage_order: list[str]):
    med = stage_medians(pt, stages, stage_order)
    ok_end = endpoint_ok(med)
    ok_mono = monotonic_ok(med, stage_order)
    corr = corr_progress(pt, prog)

    print(f"[CANDIDATE] {name}")
    for st in stage_order:
        print(f"  {st:8s} median={med[st]:.4f}")
    print(f"  endpoint_ok={ok_end}, monotonic_ok={ok_mono}, corr(progress,pt)={corr:.4f}")
    return med, ok_end, ok_mono, corr

# Sanity check for the 'Origin' of the trajectory (the earliest 1% of cells).
def early_1pct_report(adata: sc.AnnData, pt_key="pseudotime"):
    o = adata.obs
    pt = pd.to_numeric(o[pt_key], errors="coerce").to_numpy(dtype=float)
    pt_f = pt[np.isfinite(pt)]
    if pt_f.size < 100:
        print("[SANITY] early 1% skipped (too few finite pseudotime)")
        return
    q = float(np.quantile(pt_f, 0.01))
    m = pt <= q
    print(f"[SANITY] {pt_key} 1% quantile = {q:.6f}")
    print("[SANITY] early 1% stage composition:")
    print(o.loc[m, "stage"].astype(str).value_counts().to_string())
    print("[SANITY] early 1% top GSM:")
    print(o.loc[m, "gsm"].astype(str).value_counts().head(10).to_string())


# =========================
# Main
# =========================
def main():
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = get_run_dir()
    adir = run_dir / "adata"
    figdir = run_dir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    in_h5ad = adir / f"C01_umap_cluster_{DATASET_TAG}.h5ad"
    if not in_h5ad.exists():
        in_h5ad2 = adir / "C1_umap_cluster.h5ad"
        if in_h5ad2.exists():
            in_h5ad = in_h5ad2
        else:
            raise FileNotFoundError(f"Missing input: {in_h5ad} (or {in_h5ad2}). Run C01 first.")

    print("[INFO] reading:", in_h5ad)
    adata = sc.read_h5ad(in_h5ad)

    # Standardize stage names
    if "stage" not in adata.obs.columns:
        raise ValueError("obs['stage'] missing.")
    adata.obs["stage"] = adata.obs["stage"].astype(str).map(standardize_stage_name)

    # Exclude H1
    before_n = adata.n_obs
    m_keep = ~adata.obs["stage"].astype(str).isin(EXCLUDE_STAGES)
    adata = adata[m_keep].copy()
    print(f"[INFO] exclude {','.join(sorted(EXCLUDE_STAGES))}: {before_n} -> {adata.n_obs}")

    # Keep donor=0618 only
    if "donor" not in adata.obs.columns:
        raise ValueError("obs['donor'] missing.")
    before_n = adata.n_obs
    adata = adata[adata.obs["donor"].astype(str) == str(DONOR_KEEP)].copy()
    print(f"[INFO] keep donor={DONOR_KEEP}: {before_n} -> {adata.n_obs}")

    # Keep only specified stages
    before_n = adata.n_obs
    adata = adata[adata.obs["stage"].astype(str).isin(STAGES_KEEP)].copy()
    print(f"[INFO] keep stages={STAGES_KEEP}: {before_n} -> {adata.n_obs}")
    print("[INFO] stage counts after filtering:")
    print(adata.obs["stage"].astype(str).value_counts().to_string())

    # Ensure embeddings
    ensure_embeddings(adata)

    # Root = Somatic
    iroot, root_gsm = pick_somatic_root_cell(adata)
    print(f"[INFO] donor={DONOR_KEEP} root iroot={iroot}, root_gsm={root_gsm}")

    # DPT
    dpt_raw = run_dpt(adata, iroot=iroot, n_dcs=N_DCS)

    # Endpoint-first orientation for DPT
    def orient_dpt(raw: np.ndarray) -> np.ndarray:
        pt = normalize_01(raw)
        st = adata.obs["stage"].astype(str).to_numpy()
        ms = (st == "Somatic")
        mh = (st == "hCiPSC")
        med_s = float(np.nanmedian(pt[ms])) if ms.any() else np.nan
        med_h = float(np.nanmedian(pt[mh])) if mh.any() else np.nan
        if np.isfinite(med_s) and np.isfinite(med_h) and (med_h <= med_s):
            pt = 1.0 - pt
        return normalize_01(pt)

    pt_dpt = orient_dpt(dpt_raw)
    adata.obs["pseudotime_dpt"] = pt_dpt

    # Stage-progress pseudotime (guarantees stage order)
    prog = build_progress(adata)
    pt_prog = normalize_01(prog)
    adata.obs["pseudotime_progress"] = pt_prog

    # Optional isotonic calibration of DPT to progress
    pt_iso = None
    if HAS_SKLEARN:
        x = pt_dpt.copy()
        y = prog.copy()
        m = np.isfinite(x) & np.isfinite(y)
        if m.sum() > 200:
            ir = IsotonicRegression(increasing=True, out_of_bounds="clip")
            yhat = ir.fit_transform(x[m], y[m])
            # predict for all points via piecewise mapping:
            # simplest: refit on sorted x and interpolate
            xs = x[m]
            ys = yhat
            order = np.argsort(xs)
            xs2 = xs[order]
            ys2 = ys[order]
            pt_iso = np.interp(x, xs2, ys2, left=ys2[0], right=ys2[-1])
            pt_iso = normalize_01(pt_iso)
            adata.obs["pseudotime_dpt_isotonic"] = pt_iso
        else:
            print("[INFO] sklearn available but too few finite points for isotonic; skip.")
    else:
        print("[INFO] sklearn not available; skip isotonic calibration.")

    # ===== Candidate selection =====
    stages = adata.obs["stage"].astype(str).to_numpy()

    print("\n========== PSEUDOTIME CANDIDATES ==========")
    med_dpt, ok_end_dpt, ok_mono_dpt, corr_dpt = print_candidate_report("dpt", pt_dpt, stages, prog, STAGES_KEEP)

    cands = [("dpt", pt_dpt, ok_end_dpt, ok_mono_dpt, corr_dpt)]

    if pt_iso is not None:
        med_iso, ok_end_iso, ok_mono_iso, corr_iso = print_candidate_report("dpt_isotonic", pt_iso, stages, prog, STAGES_KEEP)
        cands.append(("dpt_isotonic", pt_iso, ok_end_iso, ok_mono_iso, corr_iso))

    med_prog, ok_end_prog, ok_mono_prog, corr_prog = print_candidate_report("stage_progress", pt_prog, stages, prog, STAGES_KEEP)
    cands.append(("stage_progress", pt_prog, ok_end_prog, ok_mono_prog, corr_prog))

    # pick best by (endpoint_ok, monotonic_ok, corr)
    def score(t):
        name, pt, ok_end, ok_mono, corr = t
        c = -np.inf if not np.isfinite(corr) else corr
        return (1 if ok_end else 0, 1 if ok_mono else 0, c)

    best = max(cands, key=score)
    best_name, best_pt, best_end, best_mono, best_corr = best

    adata.obs["pseudotime"] = best_pt
    print("==========================================")
    print(f"[INFO] pseudotime method selected: {best_name} | endpoint_ok={best_end}, monotonic_ok={best_mono}, corr={best_corr:.4f}\n")

    # Quick endpoint report
    ms = (stages == "Somatic")
    mh = (stages == "hCiPSC")
    print(f"[INFO] median(pseudotime) in Somatic = {float(np.nanmedian(best_pt[ms])):.4f}")
    print(f"[INFO] median(pseudotime) in hCiPSC = {float(np.nanmedian(best_pt[mh])):.4f}")

    # Store meta
    adata.uns["iroot"] = int(iroot)
    adata.uns["iroot_somatic"] = int(iroot)
    adata.uns["root_gsm"] = str(root_gsm)
    adata.uns["donor_kept"] = str(DONOR_KEEP)
    adata.uns["stages_kept"] = list(STAGES_KEEP)
    adata.uns["pseudotime_method"] = str(best_name)

    # Cluster order + MST arrows (using FINAL pseudotime)
    cent, _, pca_centroids = compute_cluster_stats(adata, cluster_key="leiden", pt_key="pseudotime")
    assign_leiden_time(adata, cent, cluster_key="leiden")
    edges = mst_edges_from_pca(pca_centroids)

    # GSM order (for timepoint plot)
    order_gsms, _ = auto_build_order_gsms(adata)

    suffix = f"{DATASET_TAG}_donor{DONOR_KEEP}_only_{ts}"

    # ===== Plots =====
    fig1, ax1 = plt.subplots(figsize=(7, 6))
    sc.pl.umap(
        adata,
        color="pseudotime",
        ax=ax1,
        show=False,
        title=f"pseudotime ({best_name}) | donor={DONOR_KEEP}",
    )
    out1 = figdir / f"umap_pseudotime_somaticroot_noH1_{suffix}.png"
    fig1.savefig(out1, dpi=220, bbox_inches="tight")
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(7, 6))
    sc.pl.umap(
        adata,
        color="pseudotime",
        ax=ax2,
        show=False,
        title=f"pseudotime ({best_name}) + MST arrows",
    )
    draw_mst_arrows_on_umap(ax2, cent, edges)
    out2 = figdir / f"umap_pseudotime_arrows_mst_somaticroot_noH1_{suffix}.png"
    fig2.savefig(out2, dpi=220, bbox_inches="tight")
    plt.close(fig2)

    fig3, ax3 = plt.subplots(figsize=(7, 6))
    sc.pl.umap(
        adata,
        color="leiden_time",
        ax=ax3,
        show=False,
        legend_loc="on data",
        title="leiden_time (clusters ordered by pseudotime)",
    )
    out3 = figdir / f"umap_leiden_time_somaticroot_noH1_{suffix}.png"
    fig3.savefig(out3, dpi=220, bbox_inches="tight")
    plt.close(fig3)

    cent_out = figdir / f"cluster_pseudotime_order_somaticroot_noH1_{suffix}.tsv"
    cent.to_csv(cent_out, sep="\t", index=False)

    # Timepoint plot
    ordered_labels, palette_list, _ = build_timepoint_labels_and_palette(adata, order_gsms)
    cent_tp = compute_timepoint_centroids(adata, ordered_labels)
    cent_tp_out = figdir / f"timepoint_umap_centroids_somaticroot_noH1_{suffix}.tsv"
    cent_tp.to_csv(cent_tp_out, sep="\t", index=False)

    fig4, ax4 = plt.subplots(figsize=(14, 6))
    sc.pl.umap(
        adata,
        color="timepoint_gsm_ordered",
        ax=ax4,
        show=False,
        palette=palette_list,
        size=6,
        alpha=0.85,
        legend_loc="right margin",
        title=f"UMAP colored by timepoint (auto order; stage families) | donor={DONOR_KEEP}",
    )
    out4 = figdir / f"umap_timepoint_ordered_somaticroot_noH1_{suffix}.png"
    fig4.savefig(out4, dpi=240, bbox_inches="tight")
    plt.close(fig4)
    print(f"[OK] wrote donor timepoint plot: {out4}")
    print(f"[OK] wrote donor centroids:    {cent_tp_out}")

    # Save h5ad
    out_h5ad = adir / f"C02_pseudotime_somaticroot_noH1_{suffix}.h5ad"
    adata.write_h5ad(out_h5ad)

    # ===== SANITY CHECKS (FINAL) =====
    print("\n========== SANITY CHECKS (FINAL pseudotime) ==========")
    pt_final = pd.to_numeric(adata.obs["pseudotime"], errors="coerce").to_numpy(dtype=float)
    print(f"[SANITY] finite: n_obs={adata.n_obs}, n_nan={int(np.isnan(pt_final).sum())}, n_inf={int(np.isinf(pt_final).sum())}, "
          f"min={np.nanmin(pt_final):.4f}, max={np.nanmax(pt_final):.4f}")

    med_final = stage_medians(pt_final, stages, STAGES_KEEP)
    print("[SANITY] stage medians (expected increasing):")
    for st in STAGES_KEEP:
        n = int((stages == st).sum())
        print(f"  {st:8s} n={n:6d} median={med_final[st]:.4f}")

    if monotonic_ok(med_final, STAGES_KEEP):
        print("[SANITY] stage median monotonicity: OK")
    else:
        print("[WARN] stage median monotonicity: FAILED (check graph / labels / consider strict progress only)")

    corr_f = corr_progress(pt_final, prog)
    print(f"[SANITY] corr(progress, pseudotime) = {corr_f:.4f}")

    early_1pct_report(adata, pt_key="pseudotime")
    print("===============================================\n")

    print(f"[OK] wrote: {out_h5ad}")
    print(f"[OK] wrote figure: {out1}")
    print(f"[OK] wrote figure: {out2}")
    print(f"[OK] wrote figure: {out3}")
    print(f"[OK] wrote cluster order: {cent_out}")
    print(f"[INFO] MIN_CELLS_PER_CLUSTER={MIN_CELLS_PER_CLUSTER}, N_PCS_FOR_MST={N_PCS_FOR_MST}")
    print(f"[INFO] run dir: {run_dir}")


if __name__ == "__main__":
    main()
