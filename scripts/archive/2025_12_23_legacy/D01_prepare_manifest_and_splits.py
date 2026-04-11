# E:\scgpt\scripts\2025_12_28\D01_prepare_manifest_and_splits.py
import os
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


# =========================
# 0) Paths (EDIT IF NEEDED)
# =========================
MODEL_DIR = Path(r"E:\scgpt\models\scgpt_whole_human")
ADATA_IN  = Path(r"E:\scgpt\results\C_traj\20251224_153430\adata\C4_wot_labels_scoreTarget.h5ad")
OUTDIR    = Path(r"E:\scgpt\results\D_finetune\run_20251228_1655")

# =========================
# 1) Settings
# =========================
SEED = 0
MAX_CELLS_PER_GSM = 5000

# Regression label: MVP 默认用 p_iPSC（更贴近“诱导成功概率”）
REG_KEY_PRIMARY = "p_iPSC"
REG_KEY_FALLBACK = "pseudotime"

# Time split rule:
#   train: day_numeric <= 2
#   val:   day_numeric in {4, 8}
#   test:  day_numeric >= 12   (包含 12 以及你设置的终点占位如 999)
TIME_VAL_DAYS = {4.0, 8.0}
TIME_TEST_MIN = 12.0


def assert_exists(p: Path, name: str):
    if not p.exists():
        raise FileNotFoundError(f"[ERROR] Missing {name}: {p}")


def downsample_by_gsm(obs: pd.DataFrame, cell_ids, max_cells_per_gsm: int, rng: np.random.Generator):
    """Downsample each GSM to at most max_cells_per_gsm. Keep all if <= max."""
    sub = obs.loc[cell_ids, ["gsm"]].copy()
    out = []
    for gsm, idx in sub.groupby("gsm").groups.items():
        idx = np.array(list(idx))
        if len(idx) > max_cells_per_gsm:
            idx = rng.choice(idx, size=max_cells_per_gsm, replace=False)
        out.append(idx)
    if len(out) == 0:
        return []
    return np.concatenate(out).tolist()


def main():
    # ---------- check inputs ----------
    assert_exists(MODEL_DIR / "vocab.json", "MODEL vocab.json")
    assert_exists(MODEL_DIR / "args.json",  "MODEL args.json")
    assert_exists(MODEL_DIR / "best_model.pt", "MODEL best_model.pt")
    assert_exists(ADATA_IN, "Input h5ad")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    (OUTDIR / "splits").mkdir(exist_ok=True)
    (OUTDIR / "data").mkdir(exist_ok=True)
    (OUTDIR / "logs").mkdir(exist_ok=True)

    rng = np.random.default_rng(SEED)

    print("[INFO] Reading:", ADATA_IN)
    adata = sc.read_h5ad(str(ADATA_IN))
    print(f"[INFO] Loaded adata: n_cells={adata.n_obs}, n_genes={adata.n_vars}")

    # ---------- required columns ----------
    required = ["stage", "gsm", "day_numeric"]
    for k in required:
        if k not in adata.obs.columns:
            raise ValueError(f"[ERROR] Missing adata.obs['{k}']")

    # ---------- choose regression key ----------
    reg_key = None
    if REG_KEY_PRIMARY in adata.obs.columns:
        reg_key = REG_KEY_PRIMARY
    elif REG_KEY_FALLBACK in adata.obs.columns:
        reg_key = REG_KEY_FALLBACK
    else:
        raise ValueError(f"[ERROR] Need regression label: '{REG_KEY_PRIMARY}' or '{REG_KEY_FALLBACK}' in adata.obs")

    print(f"[INFO] Regression label key = {reg_key}")

    # ---------- intersect genes with pretrained vocab ----------
    with open(MODEL_DIR / "vocab.json", "r", encoding="utf-8") as f:
        vocab = json.load(f)
    vocab_genes = set(vocab.keys())

    keep_genes = [g for g in adata.var_names if g in vocab_genes]
    print(f"[INFO] Gene alignment: {len(keep_genes)} / {adata.n_vars} genes kept (intersect with vocab)")

    if len(keep_genes) == 0:
        raise RuntimeError("[ERROR] No genes overlap between adata.var_names and pretrained vocab.json keys")

    adata = adata[:, keep_genes].copy()

    # ---------- write fixed inputs ----------
    used_genes_path = OUTDIR / "data" / "used_genes.tsv"
    pd.Series(adata.var_names, name="gene").to_csv(used_genes_path, sep="\t", index=False, header=False)

    shutil.copy2(MODEL_DIR / "vocab.json", OUTDIR / "data" / "vocab.json")
    shutil.copy2(MODEL_DIR / "args.json",  OUTDIR / "data" / "pretrain_args.json")

    # record label keys for reproducibility
    with open(OUTDIR / "data" / "label_keys.json", "w", encoding="utf-8") as f:
        json.dump(
            {"cls_key": "stage", "reg_key": reg_key, "seed": SEED, "max_cells_per_gsm": MAX_CELLS_PER_GSM},
            f, ensure_ascii=False, indent=2
        )

    # ---------- manifest (cells per gsm/stage/day) ----------
    df = adata.obs[["gsm", "stage", "day_numeric"]].copy()
    df["day_numeric"] = pd.to_numeric(df["day_numeric"], errors="coerce")
    df["n_cells"] = 1
    manifest = df.groupby(["gsm", "stage", "day_numeric"], as_index=False)["n_cells"].sum()
    manifest_path = OUTDIR / "data" / "manifest_cells_per_gsm.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    print("[OK] Wrote manifest:", manifest_path)

    # ---------- split 1: GSM holdout (MAIN) ----------
    gsms = sorted(df["gsm"].astype(str).unique().tolist())
    rng.shuffle(gsms)

    n_test = max(1, int(round(0.2 * len(gsms))))
    test_gsms = set(gsms[:n_test])
    train_pool = gsms[n_test:]
    n_val = max(1, int(round(0.1 * len(train_pool))))
    val_gsms = set(train_pool[:n_val])
    train_gsms = set(train_pool[n_val:])

    train_ids = adata.obs_names[adata.obs["gsm"].astype(str).isin(train_gsms)].tolist()
    val_ids   = adata.obs_names[adata.obs["gsm"].astype(str).isin(val_gsms)].tolist()
    test_ids  = adata.obs_names[adata.obs["gsm"].astype(str).isin(test_gsms)].tolist()

    train_ids = downsample_by_gsm(adata.obs, train_ids, MAX_CELLS_PER_GSM, rng)
    val_ids   = downsample_by_gsm(adata.obs, val_ids,   MAX_CELLS_PER_GSM, rng)
    test_ids  = downsample_by_gsm(adata.obs, test_ids,  MAX_CELLS_PER_GSM, rng)

    gsm_split = {
        "train": train_ids,
        "val": val_ids,
        "test": test_ids,
        "train_gsms": sorted(list(train_gsms)),
        "val_gsms": sorted(list(val_gsms)),
        "test_gsms": sorted(list(test_gsms)),
        "seed": SEED,
        "max_cells_per_gsm": MAX_CELLS_PER_GSM,
        "note": "MAIN split: hold out GSMs to reduce sample overfitting. Each GSM downsampled to <= MAX_CELLS_PER_GSM."
    }

    gsm_split_path = OUTDIR / "splits" / "gsm_split.json"
    with open(gsm_split_path, "w", encoding="utf-8") as f:
        json.dump(gsm_split, f, ensure_ascii=False, indent=2)
    print("[OK] Wrote split:", gsm_split_path)
    print(f"[INFO] gsm_split sizes: train={len(train_ids)} val={len(val_ids)} test={len(test_ids)} (GSM test={len(test_gsms)})")

    # ---------- split 2: time split (SUPPLEMENT) ----------
    day = pd.to_numeric(adata.obs["day_numeric"], errors="coerce").astype(float)

    train_ids_t = adata.obs_names[day <= 2.0].tolist()
    val_ids_t   = adata.obs_names[day.isin(list(TIME_VAL_DAYS))].tolist()
    test_ids_t  = adata.obs_names[day >= TIME_TEST_MIN].tolist()

    train_ids_t = downsample_by_gsm(adata.obs, train_ids_t, MAX_CELLS_PER_GSM, rng)
    val_ids_t   = downsample_by_gsm(adata.obs, val_ids_t,   MAX_CELLS_PER_GSM, rng)
    test_ids_t  = downsample_by_gsm(adata.obs, test_ids_t,  MAX_CELLS_PER_GSM, rng)

    time_split = {
        "train": train_ids_t,
        "val": val_ids_t,
        "test": test_ids_t,
        "seed": SEED,
        "max_cells_per_gsm": MAX_CELLS_PER_GSM,
        "rule": f"train day<=2; val day in {sorted(list(TIME_VAL_DAYS))}; test day>={TIME_TEST_MIN} (includes terminal placeholders like 999)"
    }

    time_split_path = OUTDIR / "splits" / "time_split.json"
    with open(time_split_path, "w", encoding="utf-8") as f:
        json.dump(time_split, f, ensure_ascii=False, indent=2)
    print("[OK] Wrote split:", time_split_path)
    print(f"[INFO] time_split sizes: train={len(train_ids_t)} val={len(val_ids_t)} test={len(test_ids_t)}")

    # ---------- write a small summary table ----------
    summary = []
    summary.append(["gsm_split", "train", len(train_ids), len(train_gsms)])
    summary.append(["gsm_split", "val",   len(val_ids),   len(val_gsms)])
    summary.append(["gsm_split", "test",  len(test_ids),  len(test_gsms)])
    summary.append(["time_split","train", len(train_ids_t), np.nan])
    summary.append(["time_split","val",   len(val_ids_t),   np.nan])
    summary.append(["time_split","test",  len(test_ids_t),  np.nan])

    summary_df = pd.DataFrame(summary, columns=["split", "subset", "n_cells", "n_gsms"])
    summary_path = OUTDIR / "data" / "split_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    print("[OK] Wrote summary:", summary_path)

    print("\n[DONE] OUTDIR =", OUTDIR)


if __name__ == "__main__":
    main()
