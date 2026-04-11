import re
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import matplotlib as mpl
from pathlib import Path

BASE_RESULTS = Path(r"E:\scgpt\results\C_traj")
RUN = Path((BASE_RESULTS / "_latest_run.txt").read_text(encoding="utf-8").strip())

ADIR = RUN / "adata"
FIGDIR = RUN / "figures"
WOTDIR = RUN / "wot"
WOTIN = WOTDIR / "input"

sc.settings.figdir = str(FIGDIR)
sc.settings.autoshow = False

# ---------- tunable end-point parameters----------
TOP_FRAC_TARGET = 0.10  #  the top 10% in Stage III Day 12 is regarded as the "successful end-point subset".
RIGHT_TRAJ_Q = 0.80     # early cells top20% as right_traj (quantile 0.8)

PLURI_GENES = [
    "POU5F1", "NANOG", "SOX2", "DPPA4", "LIN28A", "PRDM14", "TDGF1",
    "DNMT3B", "ZFP42", "UTF1", "NODAL", "TERT", "SALL4", "KLF4", "MYC"
]

CMAP_VMIN = 0.25
CMAP_VMAX = 0.90
# ---------------------------------


# ------------------------------ helpers ------------------------------
def find_tmap_files() -> list[Path]:
    # your tmaps are directly under RUN/wot/
    files = []
    files += list(WOTDIR.glob("tmaps_*.h5ad"))
    files += list((WOTDIR / "tmaps").glob("tmaps_*.h5ad"))
    files += list((WOTDIR / "tmaps").glob("*.h5ad"))
    return sorted(set(files))

# Extract the starting time point ($t_0$) and the ending time point ($t_1$) 
# from the specific file name output by WOT (Waddington-OT).
def parse_pair(fp: Path):
    # tmaps_<t0>_<t1>.h5ad ; t0/t1 can be floats
    m = re.search(r"tmaps_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)\.h5ad$", fp.name)
    if not m:
        return None
    return float(m.group(1)), float(m.group(2))

# Identify the storage location of the cell ID in the WOT output file
def decide_id_scheme(sample_tmap: Path) -> bool:
    tm = sc.read_h5ad(sample_tmap)
    v0 = str(tm.var_names[0]) if tm.n_vars > 0 else ""
    o0 = str(tm.obs_names[0]) if tm.n_obs > 0 else ""
    return v0.startswith("cell_") or o0.startswith("cell_")

# By prefixing "cell_" to the unified ID, all IDs can be forcibly converted to string type, 
# completely eliminating the ambiguity caused by numeric IDs.
def orig_to_wot_id(orig: str, has_prefix: bool) -> str:
    if has_prefix and not orig.startswith("cell_"):
        return "cell_" + orig
    return orig

# unifie the name of cell stage
def normalize_stage_name(s: str) -> str:
    x = str(s).strip()
    low = x.lower()
    if "hcipsc" in low:
        return "hCiPSC"
    if "stageiii" in low:
        return "StageIII"
    if "stageii" in low:
        return "StageII"
    if re.search(r"stagei(?!i)", low):
        return "StageI"
    return x

# Locate the position of the search sample or time point label
def pick_timepoint_label(adata: sc.AnnData) -> str:
    for k in ["sample_name", "gsm", "GSM"]:
        if k in adata.obs:
            return k
    return ""

# Convert the floating-point number representing the "number of days" for cell culture into a string.
def fmt_day(d: float) -> str:
    if pd.isna(d):
        return "NA"
    if abs(d - round(d)) < 1e-9:
        return str(int(round(d)))
    return str(d)

# Uniformly sample $n$ positions within the range [vmin, vmax] from the specified Matplotlib color mapping table, 
# and return the corresponding list of hexadecimal (HEX) colors.
def sample_cmap_hex(cmap_name: str, n: int, vmin: float, vmax: float) -> list[str]:
    cmap = mpl.cm.get_cmap(cmap_name)
    if n <= 1:
        vals = np.array([(vmin + vmax) / 2.0])
    else:
        vals = np.linspace(vmin, vmax, n)
    return [mpl.colors.to_hex(cmap(float(v))) for v in vals]

# construct a gradient color scheme from light to dark for single-cell time-series data
def build_success_timepoint_palette(adata: sc.AnnData):
    stage_order = ["StageI", "StageII", "StageIII"]
    stage_rank = {s: i for i, s in enumerate(stage_order)}
    stage_to_cmap = {"StageI": "Blues", "StageII": "Oranges", "StageIII": "Greens"}

    tp_key = pick_timepoint_label(adata)
    if tp_key:
        tp = adata.obs[tp_key].astype(str)
    else:
        tp = adata.obs["stage_std"].astype(str) + "_Day" + adata.obs["day_in_stage"].map(fmt_day).astype(str)

    adata.obs["timepoint_label"] = tp

    meta = adata.obs[["timepoint_label", "stage_std", "day_in_stage"]].copy()
    meta["day_in_stage"] = pd.to_numeric(meta["day_in_stage"], errors="coerce")

    tab = (
        meta.groupby("timepoint_label", observed=True)
        .agg(
            stage=("stage_std", lambda s: s.mode().iloc[0] if len(s.mode()) else s.iloc[0]),
            day=("day_in_stage", "median"),
            n=("stage_std", "size"),
        )
        .reset_index()
    )

    tab = tab[tab["stage"].isin(stage_order)].copy()
    tab["stage_rank"] = tab["stage"].map(stage_rank).astype(int)
    tab = tab.sort_values(["stage_rank", "day", "timepoint_label"]).reset_index(drop=True)

    color_map = {}
    ordered_labels = []
    for st in stage_order:
        sub = tab[tab["stage"] == st].copy()
        labels = sub["timepoint_label"].tolist()
        ordered_labels.extend(labels)

        cols = sample_cmap_hex(stage_to_cmap[st], len(labels), CMAP_VMIN, CMAP_VMAX)
        for lab, c in zip(labels, cols):
            color_map[lab] = c

    return ordered_labels, color_map

# By rounding off to eliminate minor numerical errors, 
# ensure that the time points can precisely match when used as indices or keys.
def _tkey(x: float, nd: int = 6) -> float:
    return float(np.round(float(x), nd))

# Sync the "continuous time axis" information (day_wot) that was previously calculated 
# from the external mapping table to the current single-cell data object (adata)
def load_time_map_and_attach(adata: sc.AnnData) -> pd.DataFrame:
    """
    Read timepoint_map.tsv, attach:
      stage_std, day_in_stage, day_wot
    """
    time_map = pd.read_csv(WOTIN / "timepoint_map.tsv", sep="\t")
    if not {"stage", "day", "day_wot"}.issubset(set(time_map.columns)):
        raise ValueError(f"timepoint_map.tsv must have columns stage/day/day_wot. got: {time_map.columns.tolist()}")

    # normalize
    time_map["stage_std"] = time_map["stage"].astype(str).map(normalize_stage_name)
    time_map["day_in_stage"] = pd.to_numeric(time_map["day"], errors="coerce").astype(float)
    time_map["day_wot"] = pd.to_numeric(time_map["day_wot"], errors="coerce").astype(float)

    adata.obs["stage_std"] = adata.obs["stage"].astype(str).map(normalize_stage_name)
    if "day" in adata.obs:
        adata.obs["day_in_stage"] = pd.to_numeric(adata.obs["day"], errors="coerce").astype(float)
    elif "day_in_stage" in adata.obs:
        adata.obs["day_in_stage"] = pd.to_numeric(adata.obs["day_in_stage"], errors="coerce").astype(float)
    else:
        raise ValueError("adata.obs must have 'day' or 'day_in_stage'")

    adata.obs = adata.obs.merge(
        time_map[["stage_std", "day_in_stage", "day_wot"]],
        on=["stage_std", "day_in_stage"],
        how="left",
        sort=False,
    )

    if adata.obs["day_wot"].isna().any():
        bad = adata.obs.loc[adata.obs["day_wot"].isna(), ["stage_std", "day_in_stage"]].drop_duplicates()
        raise ValueError(f"Some cells have missing day_wot mapping:\n{bad}")

    return time_map

# Build unique timepoints sorted by day_wot, assign tp_idx (0..K-1).
# This MUST match the idx naming you saw in tmaps (0..14).
def build_tp_table_with_idx(time_map: pd.DataFrame) -> pd.DataFrame:
    tp = time_map[["stage_std", "day_in_stage", "day_wot"]].drop_duplicates().copy()
    tp = tp.sort_values(["day_wot", "stage_std", "day_in_stage"]).reset_index(drop=True)
    tp["tp_idx"] = np.arange(tp.shape[0], dtype=int)
    return tp

# Determine whether the file names generated by WOT use the "original day mode" or the "index number mode".
def detect_idx_mode_from_tmaps(times: list[float], K: int) -> bool:
    """
    idx mode if all times are integer-like and max is close to K-1
    """
    if not times:
        return False
    all_intlike = all(abs(t - round(t)) < 1e-6 for t in times)
    tmax = max(times)
    return all_intlike and abs(tmax - float(K - 1)) < 1e-6


# ------------------------------ main ------------------------------
def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)

    # 1) load main adata
    adata = sc.read_h5ad(ADIR / "C2_pseudotime.h5ad")

    # 2) attach day_wot from timepoint_map.tsv
    time_map = load_time_map_and_attach(adata)
    tp_tbl = build_tp_table_with_idx(time_map)
    K = int(tp_tbl.shape[0])

    # 3) locate target timepoint: StageIII Day12 -> day_wot & tp_idx
    row = tp_tbl[
        (tp_tbl["stage_std"] == "StageIII") &
        (np.isclose(tp_tbl["day_in_stage"].astype(float), 12.0))
    ]
    if row.shape[0] != 1:
        raise ValueError(f"Cannot uniquely locate StageIII Day12 in timepoint_map. Found:\n{row}")

    TARGET_DAYWOT = float(row["day_wot"].iloc[0])  # 28.0
    TARGET_IDX = int(row["tp_idx"].iloc[0])        # should be 13 with your table

    # 4) pluripotency score
    genes_upper = {g.upper(): g for g in adata.var_names}
    used = []
    for g in PLURI_GENES:
        gu = g.upper()
        if gu in genes_upper:
            used.append(genes_upper[gu])
        elif g in adata.var_names:
            used.append(g)
    used = sorted(set(used))
    if len(used) < 4:
        raise ValueError(f"Too few pluripotency genes found (found {len(used)}): {used}")

    sc.tl.score_genes(adata, gene_list=used, score_name="pluri_score", use_raw=False)

    # 5) target subset within StageIII Day12 (use day_in_stage + day_wot)
    mask_T = (
        (adata.obs["stage_std"].astype(str) == "StageIII") &
        (np.isclose(adata.obs["day_in_stage"].astype(float), 12.0)) &
        (np.isclose(adata.obs["day_wot"].astype(float), TARGET_DAYWOT))
    )
    if mask_T.sum() == 0:
        raise ValueError("No cells found for StageIII Day12 in adata.obs")

    scores_T = adata.obs.loc[mask_T, "pluri_score"].astype(float).values
    thr = float(np.quantile(scores_T, 1.0 - TOP_FRAC_TARGET))
    mask_target = mask_T & (adata.obs["pluri_score"].astype(float) >= thr)
    adata.obs["is_target_stageIII12_top"] = mask_target.astype(int)

    print("[INFO] Target timepoint day_wot =", TARGET_DAYWOT, " (tp_idx =", TARGET_IDX, ")")
    print("[INFO] Using pluripotency genes n=", len(used), ":", used)
    print("[INFO] StageIII Day12 cells =", int(mask_T.sum()))
    print("[INFO] Target subset (top {:.0f}%) size = {}".format(TOP_FRAC_TARGET * 100, int(mask_target.sum())))
    print("[INFO] pluri_score threshold =", thr)

    # 6) read tmaps
    tmap_files = find_tmap_files()
    if not tmap_files:
        raise FileNotFoundError(f"No tmaps files found under {WOTDIR}")

    pairs_raw = []
    for fp in tmap_files:
        pr = parse_pair(fp)
        if pr is None:
            continue
        t0, t1 = pr
        pairs_raw.append((_tkey(t0), _tkey(t1), fp))

    if not pairs_raw:
        raise ValueError("No tmaps files matched pattern tmaps_<t0>_<t1>.h5ad")

    times = sorted(set([t for (a, b, _) in pairs_raw for t in (a, b)]))
    idx_mode = detect_idx_mode_from_tmaps(times, K)

    print("[INFO] Detected tmaps mode =", ("idx(0..K-1)" if idx_mode else "day_wot(named)"))
    if idx_mode:
        # Use idx target; in your case should stop at 13, not 14
        TARGET_KEY = float(TARGET_IDX)

        # build pair_map and keep consecutive 0->1->...->TARGET_IDX
        pair_map = {(t0, t1): fp for (t0, t1, fp) in pairs_raw}
        pairs = []
        missing = []
        for i in range(0, TARGET_IDX):
            a = float(i)
            b = float(i + 1)
            fp = pair_map.get((_tkey(a), _tkey(b)))
            if fp is None:
                missing.append((a, b))
            else:
                pairs.append((a, b, fp))

        if missing:
            raise ValueError(
                "Missing tmaps for adjacent index pairs up to target:\n" +
                "\n".join([f"  {int(a)}->{int(b)}" for a, b in missing[:30]])
            )
        if not pairs or pairs[-1][1] != float(TARGET_IDX):
            raise ValueError("Consecutive idx tmaps do not reach the target index.")

    else:
        # day_wot named tmaps (not your current case)
        TARGET_KEY = _tkey(TARGET_DAYWOT)

        pairs_all = [(t0, t1, fp) for (t0, t1, fp) in pairs_raw if t1 <= TARGET_KEY + 1e-9]
        pair_map = {(t0, t1): fp for (t0, t1, fp) in pairs_all}
        times2 = sorted(set([t for (a, b, _) in pairs_all for t in (a, b)]))
        if TARGET_KEY not in times2:
            raise ValueError(f"TARGET day_wot={TARGET_DAYWOT} not found in tmaps. Max t1={max(times2)}")

        pairs = []
        missing = []
        for a, b in zip(times2[:-1], times2[1:]):
            if b > TARGET_KEY + 1e-9:
                break
            fp = pair_map.get((a, b))
            if fp is None:
                missing.append((a, b))
            else:
                pairs.append((a, b, fp))

        if missing:
            raise ValueError(
                "Missing tmaps for some adjacent day_wot pairs:\n" +
                "\n".join([f"  {a}->{b}" for a, b in missing[:30]])
            )
        if not pairs or pairs[-1][1] != TARGET_KEY:
            raise ValueError("Consecutive day_wot tmaps do not reach the target day_wot.")

    print("[INFO] Using tmaps pairs:", len(pairs), "last =", pairs[-1][2].name)

    # 7) id prefix scheme
    has_prefix = decide_id_scheme(pairs[-1][2])
    print("[INFO] tmaps id prefix cell_ =", has_prefix)

    # 8) indicator vector on TARGET cells (var axis of last pair)
    last_fp = pairs[-1][2]
    tm_last = sc.read_h5ad(last_fp)
    var_ids = tm_last.var_names.astype(str).tolist()

    target_orig = set(adata.obs_names[mask_target].astype(str))
    target_wot = set(orig_to_wot_id(x, has_prefix) for x in target_orig)

    v = np.array([1.0 if cid in target_wot else 0.0 for cid in var_ids], dtype=np.float64)
    w = np.ones_like(v, dtype=np.float64)

    next_cells = var_ids
    p_by_wot = {cid: float(val) for cid, val in zip(var_ids, v)}  # store target time too

    # 9) backward propagate
    for (t0, t1, fp) in pairs[::-1]:
        tm = sc.read_h5ad(fp)
        M = tm.X
        if not sp.issparse(M):
            M = sp.csr_matrix(M)
        else:
            M = M.tocsr()

        obs_ids = tm.obs_names.astype(str).tolist()
        var_ids = tm.var_names.astype(str).tolist()

        if var_ids == next_cells:
            v_aligned = v
            w_aligned = w
        else:
            pos = {c: i for i, c in enumerate(next_cells)}
            v_aligned = np.zeros(len(var_ids), dtype=np.float64)
            w_aligned = np.zeros(len(var_ids), dtype=np.float64)
            for j, cid in enumerate(var_ids):
                i = pos.get(cid)
                if i is not None:
                    v_aligned[j] = v[i]
                    w_aligned[j] = w[i]

        v0 = M @ v_aligned
        w0 = M @ w_aligned
        p0 = np.divide(v0, w0, out=np.zeros_like(v0), where=(w0 > 0))
        p0 = np.clip(p0, 0.0, 1.0)

        for cid, pv in zip(obs_ids, p0):
            p_by_wot[cid] = float(pv)

        v, w, next_cells = v0, w0, obs_ids
        print(f"[OK] backprop {t0}->{t1}  n_obs={len(obs_ids)} n_var={len(var_ids)}")

    # 10) map back to original ids
    p_vec = np.full(adata.n_obs, np.nan, dtype=np.float64)
    for i, orig in enumerate(adata.obs_names.astype(str)):
        wot_id = orig_to_wot_id(orig, has_prefix)
        p_vec[i] = p_by_wot.get(wot_id, np.nan)

    # set hCiPSC to 1 for visualization
    is_h = (adata.obs["stage_std"].astype(str) == "hCiPSC").values
    p_vec[is_h] = 1.0
    adata.obs["p_iPSC"] = p_vec

    # 11) right_traj threshold on early cells (by continuous day_wot)
    mask_early = (
        (adata.obs["day_wot"].astype(float) < TARGET_DAYWOT) &
        (adata.obs["stage_std"].astype(str) != "hCiPSC") &
        np.isfinite(adata.obs["p_iPSC"].values)
    )
    qthr = float(np.nanquantile(adata.obs.loc[mask_early, "p_iPSC"].values, RIGHT_TRAJ_Q))
    adata.obs["right_traj"] = (adata.obs["p_iPSC"].values >= qthr).astype(int)
    print("[DONE] right_traj threshold (early, top20%) =", qthr)

    # 12) coloring
    terminal = (adata.obs["stage_std"].astype(str) == "hCiPSC").to_numpy()
    fail = (adata.obs["right_traj"].astype(int).to_numpy() == 0) & (~terminal)
    success = (adata.obs["right_traj"].astype(int).to_numpy() == 1) & (~terminal)

    ordered_tp_labels, tp_color_map = build_success_timepoint_palette(adata)

    final_lab = np.full(adata.n_obs, "Fail", dtype=object)
    final_lab[success] = adata.obs["timepoint_label"].astype(str).to_numpy()[success]
    final_lab[terminal] = "Terminal"

    categories = ["Fail"] + ordered_tp_labels + ["Terminal"]
    adata.obs["traj_timepoint_color"] = pd.Categorical(final_lab, categories=categories, ordered=True)

    palette = ["#BDBDBD"]  # Fail grey
    for tp in ordered_tp_labels:
        palette.append(tp_color_map.get(tp, "#BDBDBD"))
    palette.append("#d62728")  # Terminal red

    # 13) plots
    sc.pl.umap(adata, color=["pluri_score"], save="_pluri_score.png")
    sc.pl.umap(adata, color=["is_target_stageIII12_top"], save="_target_stageIII12_top.png")
    sc.pl.umap(adata, color=["p_iPSC"], save="_p_iPSC_scoreTarget.png")
    sc.pl.umap(adata, color=["right_traj"], save="_right_traj_scoreTarget.png")
    sc.pl.umap(
        adata,
        color=["traj_timepoint_color"],
        palette=palette,
        legend_loc="right margin",
        size=8,
        alpha=0.85,
        title="Fail=grey; Success colored by timepoint shades; Terminal(hCiPSC)=red",
        save="_traj_timepoint_failgrey.png",
    )

    # 14) save
    out = ADIR / "C4_wot_labels_scoreTarget.h5ad"
    adata.write_h5ad(out)
    print("[DONE] wrote:", out)
    print("[DONE] figures in:", FIGDIR)


if __name__ == "__main__":
    main()