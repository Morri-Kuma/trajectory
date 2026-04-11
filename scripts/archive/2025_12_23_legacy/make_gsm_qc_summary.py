import re
import csv
import gzip
from pathlib import Path
import pandas as pd


def count_lines_tsv(path: Path) -> int:
    """Count lines in (possibly gzipped) tsv quickly."""
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as f:
            return sum(1 for _ in f)
    else:
        with path.open("r", encoding="utf-8", errors="ignore") as f:
            return sum(1 for _ in f)


def parse_stage_day(sample_name: str):
    """
    Parse Stage/Day from folder name like:
    GSM7229998_StageI_Day0.5-0618
    GSM5534145_StageII_Day4-0618
    """
    stage = None
    day = None
    m_stage = re.search(r"(Stage[IVX]+)", sample_name, flags=re.IGNORECASE)
    if m_stage:
        stage = m_stage.group(1)
    m_day = re.search(r"Day([0-9]+(?:\.[0-9]+)?)", sample_name, flags=re.IGNORECASE)
    if m_day:
        day = m_day.group(1)
    return stage, day


def read_metrics(metrics_csv: Path) -> dict:
    """Read cellranger metrics_summary.csv if present."""
    key_map = {
        "estimated number of cells": "estimated_cells",
        "number of reads": "number_of_reads",
        "mean reads per cell": "mean_reads_per_cell",
    }
    out = {"estimated_cells": None, "number_of_reads": None, "mean_reads_per_cell": None}

    if not metrics_csv.exists():
        return out

    with metrics_csv.open("r", encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.reader(f)
        rows = list(reader)

    if len(rows) > 0 and len(rows[0]) == 2:
        for k, v in rows:
            lk = k.strip().lower()
            for patt, stdkey in key_map.items():
                if patt in lk:
                    out[stdkey] = v.strip()
    else:
        df = pd.read_csv(metrics_csv)
        cols_lower = {c.lower(): c for c in df.columns}
        if "metric name" in cols_lower and "value" in cols_lower:
            kcol = cols_lower["metric name"]
            vcol = cols_lower["value"]
            for _, r in df.iterrows():
                lk = str(r[kcol]).strip().lower()
                for patt, stdkey in key_map.items():
                    if patt in lk:
                        out[stdkey] = str(r[vcol]).strip()

    return out


def get_cells_filtered_from_gsm_dir(gsm_dir: Path):
    """
    Layout (you made rna_seq_10x):
      GSMxxxx/
        barcodes.tsv(.gz)
        features.tsv(.gz) or genes.tsv(.gz)
        matrix.mtx(.gz)
    """
    for fname in ["barcodes.tsv.gz", "barcodes.tsv"]:
        p = gsm_dir / fname
        if p.exists():
            return count_lines_tsv(p), str(p)

    # Recursive fallback
    for p in gsm_dir.rglob("barcodes.tsv.gz"):
        return count_lines_tsv(p), str(p)
    for p in gsm_dir.rglob("barcodes.tsv"):
        return count_lines_tsv(p), str(p)

    return None, None


def find_matrix_features(gsm_dir: Path):
    """Locate matrix/features paths (non-recursive first, then recursive fallback)."""
    matrix_path = None
    features_path = None

    for fname in ["matrix.mtx.gz", "matrix.mtx"]:
        p = gsm_dir / fname
        if p.exists():
            matrix_path = str(p)
            break

    # ✅ features 优先，其次 genes（兼容老 10x）
    for fname in ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"]:
        p = gsm_dir / fname
        if p.exists():
            features_path = str(p)
            break

    # Recursive fallback
    if matrix_path is None:
        for p in gsm_dir.rglob("matrix.mtx.gz"):
            matrix_path = str(p)
            break
        if matrix_path is None:
            for p in gsm_dir.rglob("matrix.mtx"):
                matrix_path = str(p)
                break

    if features_path is None:
        for patt in ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"]:
            for p in gsm_dir.rglob(patt):
                features_path = str(p)
                break
            if features_path is not None:
                break

    return matrix_path, features_path


def build_gsm_label_map_from_tars(tar_dir: Path) -> dict:
    """
    Map:
      GSM5387815_StageII-JNKIN8-0330.tar.gz  ->  GSM5387815 : StageII-JNKIN8-0330
      GSM5534158_H1.tar.gz -> GSM5534158 : H1
    """
    m = {}
    for p in tar_dir.glob("GSM*_*.tar.gz"):
        name = p.name
        gsm = name.split("_", 1)[0]              # GSM5387815
        label = name.split("_", 1)[1].removesuffix(".tar.gz")
        m[gsm] = label
    return m


def main():
    # =========================
    # ✅ 你只需要改这里三行
    # =========================
    DATASET_TAG = "gse178325"  # 用来区分输出
    tar_dir = Path(r"E:\scgpt\data\gse178325_human").resolve()  # tar.gz 所在目录
    raw_root = (tar_dir / "rna_seq_10x").resolve()             # GSM 文件夹所在目录（你扁平化后的）
    # 如果你没做 rna_seq_10x，而是直接用 rna_seq 也行：
    # raw_root = (tar_dir / "rna_seq").resolve()

    # 输出目录：放在 raw_root 下，且带上 gse178325 标记（不影响 B 部分找文件名）
    qc_dir = raw_root / f"qc_{DATASET_TAG}"
    qc_dir.mkdir(parents=True, exist_ok=True)

    # 用 tar 文件名补回 sample_label，保证 stage/day 能解析出来
    gsm2label = build_gsm_label_map_from_tars(tar_dir)

    gsm_dirs = sorted([p for p in raw_root.glob("GSM*") if p.is_dir()])
    if len(gsm_dirs) == 0:
        print("[ERROR] No GSM directories found under:", raw_root)
        print("        Check raw_root path.")
        return

    records = []
    for d in gsm_dirs:
        gsm = d.name  # 这里就是 GSM5387815（因为目录名就是 GSM）
        label = gsm2label.get(gsm)  # StageII-JNKIN8-0330 / HEFs-0330 / H1 ...

        # ✅ 关键：拼回旧格式，复用你原来的解析逻辑
        sample_name = f"{gsm}_{label}" if label else gsm
        stage, day = parse_stage_day(sample_name)

        cells_filtered, barcode_path = get_cells_filtered_from_gsm_dir(d)
        matrix_path, features_path = find_matrix_features(d)

        metrics_csv = d / "metrics_summary.csv"
        metrics = read_metrics(metrics_csv)

        records.append(
            {
                "GSM": gsm,
                "sample_name": sample_name,
                "stage": stage,
                "day": day,
                "cells_filtered": cells_filtered,
                "estimated_cells": metrics["estimated_cells"],
                "number_of_reads": metrics["number_of_reads"],
                "mean_reads_per_cell": metrics["mean_reads_per_cell"],
                "barcode_path": barcode_path,
                "features_path": features_path,
                "matrix_path": matrix_path,
                "metrics_path": str(metrics_csv) if metrics_csv.exists() else None,
            }
        )

    df = pd.DataFrame(records)

    # --- write full debug table ---
    out_tsv = qc_dir / "gsm_qc_summary.tsv"  # 文件名保持不变，B 部分通常会依赖这个名字
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[OK] Wrote: {out_tsv}")

    # ---- Route A: only cells_filtered is required ----
    df["day_numeric"] = pd.to_numeric(df["day"], errors="coerce")

    # Label iPSC sample(s)
    mask_ipsc = df["sample_name"].str.contains("hCiPSC", case=False, na=False)
    df.loc[mask_ipsc, "stage"] = "hCiPSCs"
    df.loc[mask_ipsc, "day"] = "iPSC"
    df.loc[mask_ipsc, "day_numeric"] = 999.0

    # Label H1 (ESC) if present
    mask_h1 = df["sample_name"].str.contains(r"(_H1$|_H1-)", case=False, na=False)
    df.loc[mask_h1, "stage"] = "H1"
    df.loc[mask_h1, "day"] = "ESC"
    df.loc[mask_h1, "day_numeric"] = 1000.0

    # File integrity flags
    df["has_barcodes"] = df["barcode_path"].notna()
    df["has_features"] = df["features_path"].notna()
    df["has_matrix"] = df["matrix_path"].notna()

    # Stage ordering (✅ 加上 StageIV，避免排序奇怪)
    stage_order = {"StageI": 1, "StageII": 2, "StageIII": 3, "StageIV": 4, "hCiPSCs": 5, "H1": 6}
    df["stage_order"] = df["stage"].map(stage_order).fillna(99).astype(int)
    df = df.sort_values(["stage_order", "day_numeric"], na_position="last")

    # Write the clean Methods/QC table
    methods_tsv = qc_dir / "gsm_qc_methods_table.tsv"
    df[
        [
            "GSM",
            "sample_name",
            "stage",
            "day",
            "cells_filtered",
            "has_barcodes",
            "has_features",
            "has_matrix",
        ]
    ].to_csv(methods_tsv, sep="\t", index=False)
    print(f"[OK] Wrote methods/QC table: {methods_tsv}")

    # Debug prints
    print("[DEBUG] df shape:", df.shape)
    print("[DEBUG] df columns:", list(df.columns))
    cols_preview = [
        "GSM",
        "sample_name",
        "stage",
        "day",
        "cells_filtered",
        "has_barcodes",
        "has_features",
        "has_matrix",
    ]
    print(df.reindex(columns=cols_preview).head(30))


if __name__ == "__main__":
    main()
