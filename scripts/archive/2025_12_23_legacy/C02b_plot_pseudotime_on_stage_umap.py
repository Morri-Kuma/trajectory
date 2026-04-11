# C02b_plot_pseudotime_on_stage_umap.py
# Plot UMAP with:
# 1) stage (fixed colors)
# 2) timepoint ordered (auto ORDER_GSMS from obs: stage_rank + day)
# 3) timepoint-ordered smooth curve + arrows (centroid polyline/spline)
#
# Input : results/C_traj/<dataset>/<run>/adata/C02_pseudotime*.h5ad  (auto pick latest)
# Output: results/C_traj/<dataset>/<run>/figures/*.png + centroids.tsv (ALL with timestamp, no overwrite)

import re
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
from pathlib import Path
from datetime import datetime

# ---- optional: spline smoothing ----
try:
    from scipy.interpolate import splprep, splev
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False

# =========================
# Config
# =========================
BASE_RESULTS = Path(r"E:\scgpt\results\C_traj")
DATASET_TAG = "gse178325"   # 你的当前数据集（用于优先读取 dataset 的 _latest_run.txt）

# ===== Stage-specific colormap families (for timepoint gradient) =====
STAGE_CMAP = {
    "Somatic": "Greys",
    "StageI": "Blues",
    "StageII": "Oranges",
    "StageIII": "Greens",
    "StageIV": "Purples",
    "hCiPSC": "Reds",
    "H1": "Reds",   # 如果不想让 H1 像终点，可改成 "Greys"
}
CMAP_START = 0.25
CMAP_END = 0.95

# ===== Fixed stage colors for stage plot =====
STAGE_FIXED_PALETTE = {
    "Somatic": "#7f7f7f",  # grey
    "StageI":  "#1f77b4",  # blue
    "StageII": "#ff7f0e",  # orange
    "StageIII":"#2ca02c",  # green
    "StageIV": "#9467bd",  # purple
    "hCiPSC":  "#d62728",  # red
    "H1":      "#000000",  # black
}

# ===== Plot parameters =====
POINT_SIZE = 6
POINT_ALPHA = 0.85

# Curve / arrows
DRAW_SMOOTH_CURVE = True
SPLINE_N_POINTS = 350
SPLINE_SMOOTH_S = 0.0   # 0.0 = pass through all centroids; >0 = smoother
ARROW_STYLE = dict(arrowstyle="->", lw=3.0, color="black", alpha=0.9)
ARROW_U_STEP = 0.03
SHOW_TIMEPOINT_LABELS = True

# =========================


def get_run_dir() -> Path:
    """
    优先读取:
      E:\scgpt\results\C_traj\gse178325\_latest_run.txt
    否则 fallback:
      E:\scgpt\results\C_traj\_latest_run.txt
    """
    candidates = [
        BASE_RESULTS / DATASET_TAG / "_latest_run.txt",
        BASE_RESULTS / "_latest_run.txt",
    ]
    for latest in candidates:
        if latest.exists():
            run_dir = Path(latest.read_text(encoding="utf-8").strip())
            if run_dir.exists():
                return run_dir
    raise FileNotFoundError(
        "Cannot find a valid _latest_run.txt. "
        f"Tried: {[str(x) for x in candidates]}"
    )


def pick_latest_c02_file(adir: Path) -> Path:
    """
    在 run_dir/adata 下选择最新的 C02 输出:
      C02_pseudotime*.h5ad
    """
    cands = sorted(adir.glob("C02_pseudotime*.h5ad"),
                   key=lambda p: p.stat().st_mtime,
                   reverse=True)
    if not cands:
        # 兼容旧命名
        p = adir / "C2_pseudotime.h5ad"
        if p.exists():
            return p
        raise FileNotFoundError(f"No C02_pseudotime*.h5ad found under: {adir}")
    return cands[0]


def standardize_stage_name(x: str) -> str:
    """
    Robust stage normalization.
    CRITICAL: check StageIII BEFORE StageII, because 'stageiii' contains substring 'stageii'.
    """
    s = str(x).strip()
    low = s.lower()

    if "hcipsc" in low:
        return "hCiPSC"
    if low == "h1" or " h1" in low:
        return "H1"
    if "somatic" in low:
        return "Somatic"
    if "stageiv" in low:
        return "StageIV"
    if "stageiii" in low:
        return "StageIII"
    if "stageii" in low:
        return "StageII"
    if re.search(r"stagei(?!i)", low):
        return "StageI"
    return s


def detect_gsm_column(adata: sc.AnnData) -> str:
    for c in ["gsm", "sample_name", "sample", "batch"]:
        if c in adata.obs.columns:
            return c
    raise ValueError("Cannot find GSM column in adata.obs. Need one of: gsm / sample_name / sample / batch")


def extract_gsm(series: pd.Series) -> pd.Series:
    s = series.astype(str)

    def _get(x: str):
        m = re.search(r"(GSM\d+)", x)
        return m.group(1) if m else x

    return s.map(_get)


def get_cmap_safe(name: str):
    # 兼容新版 matplotlib（避免 get_cmap deprecated warning）
    if hasattr(mpl, "colormaps"):
        return mpl.colormaps.get_cmap(name)
    return cm.get_cmap(name)


def sample_cmap(cmap_name: str, n: int, start: float = CMAP_START, end: float = CMAP_END):
    cmap = get_cmap_safe(cmap_name)
    if n <= 1:
        return [cmap((start + end) / 2)]
    xs = np.linspace(start, end, n)
    return [cmap(float(x)) for x in xs]


def auto_build_order_gsms(ad: sc.AnnData):
    """
    从 obs 自动构建 ORDER_GSMS，排序规则：
      stage_rank (Somatic=0, StageI=1, StageII=2, StageIII=3, StageIV=4, hCiPSC=5, H1=6, else=99)
      然后 day（数值）
      然后 gsm
    """
    if "gsm" not in ad.obs.columns:
        raise ValueError("Need adata.obs['gsm'] to auto-build ORDER_GSMS.")

    tmp = ad.obs[["gsm", "stage", "day"]].copy()
    tmp["stage"] = tmp["stage"].astype(str).map(standardize_stage_name)
    tmp["day"] = pd.to_numeric(tmp["day"], errors="coerce")

    # 每个 GSM 取 stage 众数、day 中位数
    def mode_str(s):
        vc = s.astype(str).value_counts()
        return str(vc.index[0]) if len(vc) else "Unknown"

    meta = (
        tmp.groupby("gsm", observed=True)
        .agg(stage=("stage", mode_str), day=("day", "median"))
        .reset_index()
    )

    stage_rank_map = {
        "Somatic": 0,
        "StageI": 1,
        "StageII": 2,
        "StageIII": 3,
        "StageIV": 4,
        "hCiPSC": 5,
        "H1": 6,
    }
    meta["stage_rank"] = meta["stage"].map(stage_rank_map).fillna(99).astype(int)
    meta = meta.sort_values(["stage_rank", "day", "gsm"], kind="mergesort").reset_index(drop=True)

    order_gsms = meta["gsm"].astype(str).tolist()

    print("[INFO] auto ORDER_GSMS built from obs (top 20 shown):")
    print(meta.head(20).to_string(index=False))
    return order_gsms, meta


def build_timepoint_labels_and_palette(ad: sc.AnnData, ORDER_GSMS):
    """
    Build:
    - ad.obs['timepoint_gsm_ordered'] : ordered categorical by ORDER_GSMS (each GSM has its own color)
    - ordered_labels, palette_list    : stage families + gradient within stage
    """
    gsm_col = detect_gsm_column(ad)
    ad.obs["gsm_id"] = extract_gsm(ad.obs[gsm_col]).astype(str)

    # label: prefer sample_name
    gsm2label = {}
    if "sample_name" in ad.obs.columns:
        tmp = ad.obs[["gsm_id", "sample_name"]].dropna()
        for g, sub in tmp.groupby("gsm_id"):
            gsm2label[str(g)] = str(sub["sample_name"].iloc[0])
    else:
        for g in ORDER_GSMS:
            gsm2label[g] = g

    # stage per GSM: prefer ad.obs['stage']
    gsm2stage = {}
    if "stage" in ad.obs.columns:
        tmp = ad.obs[["gsm_id", "stage"]].dropna()
        for g, sub in tmp.groupby("gsm_id"):
            st = str(sub["stage"].astype(str).value_counts().index[0])
            gsm2stage[str(g)] = standardize_stage_name(st)
    else:
        for g in ORDER_GSMS:
            gsm2stage[g] = "Unknown"

    ordered_labels = [gsm2label.get(g, g) for g in ORDER_GSMS]

    # stage -> labels（保持全局时间顺序）
    stage2labels = {}
    for g, lab in zip(ORDER_GSMS, ordered_labels):
        st = gsm2stage.get(g, "Unknown")
        stage2labels.setdefault(st, []).append(lab)

    # 为每个 stage 抽色
    label2color = {}
    for st, labs in stage2labels.items():
        cmap_name = STAGE_CMAP.get(st, "Greys")
        cols = sample_cmap(cmap_name, len(labs))
        for lab, col in zip(labs, cols):
            label2color[lab] = col

    palette_list = [label2color.get(lab, (0.7, 0.7, 0.7, 1.0)) for lab in ordered_labels]

    # 写入有序分类列
    mapper = {g: gsm2label.get(g, g) for g in ORDER_GSMS}
    tp = ad.obs["gsm_id"].map(mapper)
    ad.obs["timepoint_gsm_ordered"] = pd.Categorical(tp, categories=ordered_labels, ordered=True)

    # debug prints
    print("[INFO] gsm_col used:", gsm_col)
    if "stage" in ad.obs.columns:
        print("[INFO] stage counts (after standardize):")
        print(ad.obs["stage"].value_counts().to_string())
    print("[INFO] timepoint order preview (first 15):")
    for g, lab in list(zip(ORDER_GSMS, ordered_labels))[:15]:
        st = gsm2stage.get(g, "Unknown")
        cmap_name = STAGE_CMAP.get(st, "Greys")
        print(f"  {g} -> {lab} | stage={st} | cmap={cmap_name}")

    return ordered_labels, palette_list


def compute_timepoint_centroids(ad: sc.AnnData, ordered_labels):
    if "X_umap" not in ad.obsm:
        raise ValueError("Missing adata.obsm['X_umap'].")
    umap = ad.obsm["X_umap"]
    lab = ad.obs["timepoint_gsm_ordered"]

    rows = []
    for i, l in enumerate(ordered_labels):
        m = (lab.astype(str) == str(l)).to_numpy()
        if m.sum() == 0:
            rows.append({"order": i, "label": l, "n": 0, "x": np.nan, "y": np.nan})
        else:
            xy = umap[m]
            rows.append(
                {"order": i, "label": l, "n": int(m.sum()),
                 "x": float(xy[:, 0].mean()), "y": float(xy[:, 1].mean())}
            )
    return pd.DataFrame(rows)


def draw_polyline_arrows(ax, cent_df: pd.DataFrame):
    df = cent_df.copy()
    df = df[np.isfinite(df["x"]) & np.isfinite(df["y"]) & (df["n"] > 0)].sort_values("order")
    if df.shape[0] < 2:
        print("[WARN] Too few timepoints with cells to draw arrows.")
        return

    x = df["x"].to_numpy(float)
    y = df["y"].to_numpy(float)
    orders = df["order"].to_numpy(int)

    for k in range(len(x) - 1):
        ax.annotate("", xy=(x[k + 1], y[k + 1]), xytext=(x[k], y[k]), arrowprops=ARROW_STYLE)

    if SHOW_TIMEPOINT_LABELS:
        for xi, yi, o in zip(x, y, orders):
            ax.text(
                xi, yi, str(o + 1),
                fontsize=11, ha="center", va="center", color="black",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7),
            )


def draw_smooth_curve_arrows(ax, cent_df: pd.DataFrame, n_points: int, smooth_s: float):
    df = cent_df.copy()
    df = df[np.isfinite(df["x"]) & np.isfinite(df["y"]) & (df["n"] > 0)].sort_values("order")
    if df.shape[0] < 2:
        print("[WARN] Too few timepoints with cells to draw curve.")
        return

    x = df["x"].to_numpy(float)
    y = df["y"].to_numpy(float)
    orders = df["order"].to_numpy(int)

    if (not HAS_SCIPY) or (len(x) < 4):
        if not HAS_SCIPY:
            print("[WARN] scipy not available; fallback to polyline arrows.")
        else:
            print("[WARN] too few points for spline; fallback to polyline arrows.")
        draw_polyline_arrows(ax, cent_df)
        return

    d = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)
    u = np.r_[0.0, np.cumsum(d)]
    if u[-1] == 0:
        u = np.linspace(0, 1, len(x))
    else:
        u = u / u[-1]

    k = min(3, len(x) - 1)
    tck, _ = splprep([x, y], u=u, s=smooth_s, k=k)
    uu = np.linspace(0, 1, n_points)
    xs, ys = splev(uu, tck)

    ax.plot(xs, ys, color="black", lw=3.0, alpha=0.85)

    for i in range(len(u) - 1):
        umid = (u[i] + u[i + 1]) / 2
        u2 = min(1.0, umid + ARROW_U_STEP)
        x1, y1 = splev(umid, tck)
        x2, y2 = splev(u2, tck)
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=ARROW_STYLE)

    if SHOW_TIMEPOINT_LABELS:
        for xi, yi, o in zip(x, y, orders):
            ax.text(
                xi, yi, str(o + 1),
                fontsize=11, ha="center", va="center", color="black",
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.7),
            )


def main():
    run_dir = get_run_dir()
    figdir = run_dir / "figures"
    adir = run_dir / "adata"
    figdir.mkdir(parents=True, exist_ok=True)

    in_c02 = pick_latest_c02_file(adir)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    stem = in_c02.stem  # 用输入文件名做区分
    print("[INFO] using C02 input:", in_c02)

    # 读入（backed）并构建轻量对象（只保留 obs + X_umap）
    ad0 = sc.read_h5ad(in_c02, backed="r")
    if "X_umap" not in ad0.obsm:
        raise ValueError("Input file has no X_umap. Run C01 before C02b.")
    obs = ad0.obs.copy()
    umap = np.array(ad0.obsm["X_umap"])
    ad = sc.AnnData(obs=obs, obsm={"X_umap": umap})
    # 关闭文件句柄
    ad0.file.close()

    # 标准化 stage
    if "stage" in ad.obs.columns:
        ad.obs["stage"] = ad.obs["stage"].astype(str).map(standardize_stage_name)

    # 自动 ORDER_GSMS（适配 gse178325）
    ORDER_GSMS, _meta = auto_build_order_gsms(ad)

    # -------- Plot stage with fixed palette (exact colors) --------
    if "stage" in ad.obs.columns:
        # 只保留当前数据里出现的类别，并按推荐顺序排序
        prefer_order = ["Somatic", "StageI", "StageII", "StageIII", "StageIV", "hCiPSC", "H1"]
        present = [x for x in prefer_order if x in set(ad.obs["stage"].astype(str).unique())]
        others = sorted(set(ad.obs["stage"].astype(str).unique()) - set(present))
        stage_order = present + others

        ad.obs["stage"] = pd.Categorical(ad.obs["stage"].astype(str), categories=stage_order, ordered=True)
        palette_list = [
            STAGE_FIXED_PALETTE.get(c, "#808080") for c in stage_order
        ]

        fig, ax = plt.subplots(figsize=(10, 7))
        sc.pl.umap(
            ad,
            color="stage",
            palette=palette_list,
            size=POINT_SIZE,
            alpha=POINT_ALPHA,
            legend_loc="right margin",
            ax=ax,
            show=False,
            title=f"stage (fixed colors) | source={stem}",
        )
        out_stage = figdir / f"umap_stage_from_{stem}_{ts}.png"
        fig.savefig(out_stage, dpi=240, bbox_inches="tight")
        plt.close(fig)
        print("[OK] wrote:", out_stage)
    else:
        print("[WARN] adata.obs['stage'] not found, skip stage plot.")

    # -------- Timepoint ordered palette (C02-consistent) --------
    ordered_labels, palette_tp = build_timepoint_labels_and_palette(ad, ORDER_GSMS)

    fig, ax = plt.subplots(figsize=(14, 6))
    sc.pl.umap(
        ad,
        color="timepoint_gsm_ordered",
        palette=palette_tp,
        size=POINT_SIZE,
        alpha=POINT_ALPHA,
        legend_loc="right margin",
        ax=ax,
        show=False,
        title=f"UMAP colored by timepoint (auto order) | source={stem}",
    )
    out_tp = figdir / f"umap_timepoint_ordered_from_C02b_{stem}_{ts}.png"
    fig.savefig(out_tp, dpi=240, bbox_inches="tight")
    plt.close(fig)
    print("[OK] wrote:", out_tp)

    # -------- Compute centroids & save --------
    cent_df = compute_timepoint_centroids(ad, ordered_labels)
    cent_out = figdir / f"timepoint_umap_centroids_{stem}_{ts}.tsv"
    cent_df.to_csv(cent_out, sep="\t", index=False)
    print("[OK] wrote:", cent_out)

    # -------- Plot timepoint + smooth curve arrows --------
    fig, ax = plt.subplots(figsize=(14, 6))
    sc.pl.umap(
        ad,
        color="timepoint_gsm_ordered",
        palette=palette_tp,
        size=POINT_SIZE,
        alpha=POINT_ALPHA,
        legend_loc="right margin",
        ax=ax,
        show=False,
        title=f"Timepoint-ordered curve + arrows (auto order) | source={stem}",
    )

    if DRAW_SMOOTH_CURVE:
        draw_smooth_curve_arrows(ax, cent_df, n_points=SPLINE_N_POINTS, smooth_s=SPLINE_SMOOTH_S)
    else:
        draw_polyline_arrows(ax, cent_df)

    out_curve = figdir / f"umap_timepoint_ordered_smoothcurve_arrows_{stem}_{ts}.png"
    fig.savefig(out_curve, dpi=260, bbox_inches="tight")
    plt.close(fig)
    print("[OK] wrote:", out_curve)

    if DRAW_SMOOTH_CURVE and (not HAS_SCIPY):
        print("[NOTE] smooth curve requested but scipy not found; used polyline arrows instead.")
    if DRAW_SMOOTH_CURVE and HAS_SCIPY:
        print(f"[INFO] used spline smoothing: s={SPLINE_SMOOTH_S}, n_points={SPLINE_N_POINTS}")

    print("[OK] figures in:", figdir)


if __name__ == "__main__":
    main()
