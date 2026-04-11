import re
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from pathlib import Path

BASE_RESULTS = Path(r"E:\scgpt\results\C_traj")
RUN = Path((BASE_RESULTS / "_latest_run.txt").read_text(encoding="utf-8").strip())

ADIR = RUN / "adata"
FIGDIR = RUN / "figures"
WOTDIR = RUN / "wot"
WOTIN = WOTDIR / "input"

sc.settings.figdir = str(FIGDIR)
sc.settings.autoshow = False


def find_tmap_files() -> list[Path]:
    """WOT outputs may be in RUN/wot/tmaps_*.h5ad or in RUN/wot/tmaps/*.h5ad (if --out used)."""
    files = []
    files += list(WOTDIR.glob("tmaps_*.h5ad"))          # your current case
    files += list((WOTDIR / "tmaps").glob("*.h5ad"))    # if user used --out
    # de-dup
    files = sorted(set(files))
    return files


def parse_pair(fp: Path):
    m = re.search(r"tmaps_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)\.h5ad$", fp.name)
    if not m:
        return None
    return float(m.group(1)), float(m.group(2))


def load_wot_id_mapping():
    """Return (orig->wot, wot->orig) mappings if wot_id_map.tsv exists; else None."""
    mp = WOTIN / "wot_id_map.tsv"
    if not mp.exists():
        return None, None
    df = pd.read_csv(mp, sep="\t")
    # expected columns: wot_id, orig_id
    if not {"wot_id", "orig_id"}.issubset(df.columns):
        return None, None
    o2w = dict(zip(df["orig_id"].astype(str), df["wot_id"].astype(str)))
    w2o = dict(zip(df["wot_id"].astype(str), df["orig_id"].astype(str)))
    return o2w, w2o


def main():
    FIGDIR.mkdir(parents=True, exist_ok=True)

    # Load original adata (keep UMAP, clusters, etc.)
    adata = sc.read_h5ad(ADIR / "C2_pseudotime.h5ad")

    # Make sure stage/day are present
    if "stage" not in adata.obs or "day" not in adata.obs:
        raise ValueError("adata.obs must contain 'stage' and 'day'")

    # Add t_wot to adata (same mapping as C03)
    time_map = pd.read_csv(WOTIN / "timepoint_map.tsv", sep="\t")
    adata.obs["day"] = adata.obs["day"].astype(float)
    adata.obs = adata.obs.merge(time_map[["stage", "day", "day_wot"]], on=["stage", "day"], how="left", sort=False)
    adata.obs["day_wot"] = adata.obs["day_wot"].astype(float)  
    if adata.obs["t_wot"].isna().any():
        bad = adata.obs[adata.obs["t_wot"].isna()][["stage", "day"]].drop_duplicates()
        raise ValueError(f"Some cells have missing t_wot mapping:\n{bad}")
    adata.obs["t_wot"] = adata.obs["t_wot"].astype(int)

    # Find transport maps
    tmap_files = find_tmap_files()
    if not tmap_files:
        raise FileNotFoundError(f"No transport maps found under {WOTDIR}")

    pairs = []
    for fp in tmap_files:
        pr = parse_pair(fp)
        if pr is not None:
            t0, t1 = pr
            pairs.append((t0, t1, fp))
    pairs = sorted(pairs, key=lambda x: (x[0], x[1]))
    if not pairs:
        raise ValueError("Found tmaps files but could not parse (t0,t1) from filenames.")

    print("[INFO] Found tmaps:", len(pairs))
    print("[INFO] First/last:", pairs[0][2].name, " ... ", pairs[-1][2].name)

    # Determine whether WOT ids are prefixed, and load optional mapping table
    o2w, w2o = load_wot_id_mapping()

    # We will build probabilities in WOT-id space (the ids inside tmaps matrices)
    # target: all hCiPSC cells at final t_wot
    T = float(adata.obs["day_wot"].max())
    target_orig = set(adata.obs_names[(adata.obs["stage"] == "hCiPSC") & (adata.obs["t_wot"] == T)].astype(str))

    # helper: orig_id -> wot_id (best-effort)
    def orig_to_wot(orig: str, sample_prefix: bool) -> str:
        if o2w is not None and orig in o2w:
            return o2w[orig]
        if sample_prefix and not orig.startswith("cell_"):
            return "cell_" + orig
        return orig

    # Start from last map (t=T-1 -> T)
    last_fp = pairs[-1][2]
    tm_last = sc.read_h5ad(last_fp)
    var_ids = tm_last.var_names.astype(str).tolist()

    sample_prefix = (len(var_ids) > 0 and var_ids[0].startswith("cell_"))

    target_wot = set(orig_to_wot(x, sample_prefix) for x in target_orig)

    # v_T: indicator on var_ids (cells at time T)
    v = np.array([1.0 if cid in target_wot else 0.0 for cid in var_ids], dtype=np.float64)
    w = np.ones_like(v, dtype=np.float64)
    next_cells = var_ids  # names aligned to v,w

    # store p for final time T cells too
    p_by_wot = {cid: float(val) for cid, val in zip(var_ids, v)}

    # Walk backward through all pairs
    for (t0, t1, fp) in pairs[::-1]:
        tm = sc.read_h5ad(fp)
        M = tm.X
        if not sp.issparse(M):
            M = sp.csr_matrix(M)
        else:
            M = M.tocsr()

        obs_ids = tm.obs_names.astype(str).tolist()  # time t0
        var_ids = tm.var_names.astype(str).tolist()  # time t1

        # align current v,w (for next_cells) onto this tm.var_names
        if var_ids == next_cells:
            v_aligned = v
            w_aligned = w
        else:
            pos = {c: i for i, c in enumerate(next_cells)}
            v_aligned = np.zeros(len(var_ids), dtype=np.float64)
            w_aligned = np.zeros(len(var_ids), dtype=np.float64)
            for j, cid in enumerate(var_ids):
                i = pos.get(cid, None)
                if i is not None:
                    v_aligned[j] = v[i]
                    w_aligned[j] = w[i]
                else:
                    # if missing, leave 0 (no mass / no info)
                    pass

        v0 = M @ v_aligned
        w0 = M @ w_aligned
        p0 = np.divide(v0, w0, out=np.zeros_like(v0), where=(w0 > 0))
        p0 = np.clip(p0, 0.0, 1.0)

        # record probabilities for obs_ids (time t0 cells)
        for cid, pv in zip(obs_ids, p0):
            p_by_wot[cid] = float(pv)

        # update for next iteration (moving backward)
        v, w, next_cells = v0, w0, obs_ids

        print(f"[OK] backprop {t0}->{t1}  n_obs={len(obs_ids)} n_var={len(var_ids)}")

    # Map p back to original adata.obs_names
    def wot_to_orig(wot_id: str) -> str:
        if w2o is not None and wot_id in w2o:
            return w2o[wot_id]
        if wot_id.startswith("cell_"):
            return wot_id.replace("cell_", "", 1)
        return wot_id

    # build p_iPSC for each original cell
    p_vec = np.full(adata.n_obs, np.nan, dtype=np.float64)
    for i, orig in enumerate(adata.obs_names.astype(str)):
        wot_id = orig_to_wot(orig, sample_prefix)
        p_vec[i] = p_by_wot.get(wot_id, np.nan)

    adata.obs["p_iPSC"] = p_vec

    # right_traj: top 20% among non-hCiPSC (recommended)
    mask_train = (adata.obs["stage"] != "hCiPSC") & np.isfinite(adata.obs["p_iPSC"].values)
    thr = float(np.nanquantile(adata.obs.loc[mask_train, "p_iPSC"].values, 0.80))
    adata.obs["right_traj"] = (adata.obs["p_iPSC"].values >= thr).astype(int)

    # Plot
    sc.pl.umap(adata, color=["p_iPSC"], save="_p_iPSC.png")
    sc.pl.umap(adata, color=["right_traj"], save="_right_traj.png")

    # Save
    out = ADIR / "C4_wot_labels.h5ad"
    adata.write_h5ad(out)
    print("[DONE] wrote:", out)
    print("[DONE] threshold thr(top20% non-hCiPSC) =", thr)
    print("[DONE] figures in:", FIGDIR)


if __name__ == "__main__":
    main()
