# B01_build_per_gsm_h5ad.py
import re
import pandas as pd
import scanpy as sc
from pathlib import Path

# =========================
# Config (edit if needed)
# =========================
ROOT = Path(r"E:\scgpt")  # 工程根目录

DATASET_TAG = "gse178325"

# ✅ 你的 QC 表路径（A 部分生成的）
META = Path(r"E:\scgpt\data\gse178325_human\rna_seq_10x\qc_gse178325\gsm_qc_summary.tsv")

# ✅ 输出单样本 h5ad：放到带 dataset 标记的目录，避免和旧文件混
OUTDIR = ROOT / rf"data\processed_{DATASET_TAG}\per_gsm_h5ad"
OUTDIR.mkdir(parents=True, exist_ok=True)

# 基因命名体系：二选一
# "symbol": 使用 gene symbols（HGNC）
# "ensembl": 使用 Ensembl gene id
GENE_NAMING = "symbol"  # "symbol" or "ensembl"

# donor/batch 默认值（更建议从 sample_name 里解析 0330/0618/0605...）
DEFAULT_DONOR = "unknown"


def read_10x_mtx(tenx_dir: Path):
    """用 scanpy 读 10x mtx（支持 .gz）"""
    if GENE_NAMING == "symbol":
        return sc.read_10x_mtx(tenx_dir, var_names="gene_symbols", make_unique=True)
    else:
        return sc.read_10x_mtx(tenx_dir, var_names="gene_ids", make_unique=True)


def is_missing(x) -> bool:
    """把 NaN 和字符串 'nan' 都视为缺失"""
    if x is None:
        return True
    if isinstance(x, float) and pd.isna(x):
        return True
    s = str(x).strip()
    return (s == "") or (s.lower() == "nan")


def is_hcipsc_sample(sample_name: str) -> bool:
    s = str(sample_name)
    return ("hCiPSC" in s) or ("hCiPSCs" in s)


def is_h1_sample(sample_name: str) -> bool:
    s = str(sample_name)
    return ("_H1" in s) or (s.endswith("H1"))


def parse_batch_from_sample_name(sample_name: str) -> str:
    """
    从 sample_name 末尾解析批次：-0330 / -0618 / -0605 / -0809 / -1117 / -1230
    """
    s = str(sample_name).strip()
    m = re.search(r"-(\d{4})$", s)
    if m:
        return m.group(1)
    return DEFAULT_DONOR


def tenx_dir_from_row(row) -> Path:
    """
    ✅ 最稳：直接用 A 部分写进 tsv 的 matrix_path 的父目录
    这样不依赖目录命名（GSMxxxx 或 GSMxxxx_Stage... 都没关系）
    """
    mp = row.get("matrix_path", None)
    if not is_missing(mp):
        p = Path(str(mp)).resolve()
        if p.exists():
            return p.parent
    raise FileNotFoundError(f"matrix_path missing or not exists for GSM={row.get('GSM')} sample={row.get('sample_name')}")


def parse_stage_day(row, hcipsc_day: float):
    """
    从 gsm_qc_summary.tsv 的 stage/day 字段解析并兜底修复：
    - hCiPSC/hCiPSCs：stage='hCiPSC', day=hcipsc_day（= max_non_hcipsc_day + 9）
    - H1：stage='H1', day=hcipsc_day + 1（放在最后）
    - somatic（HEFs/hASFs/hADSCs）：day = -1
    - stage-only（StageI-0618 这种 day 缺失）：用 stage anchor day 占位，保证后续不崩
    """
    sample_name = str(row["sample_name"]).strip()
    stage_raw = row.get("stage", None)
    day_raw = row.get("day", None)

    # hCiPSC
    if is_hcipsc_sample(sample_name):
        return "hCiPSC", float(hcipsc_day)

    # H1
    if is_h1_sample(sample_name):
        return "H1", float(hcipsc_day + 1.0)

    # stage
    stage = str(stage_raw).strip() if not is_missing(stage_raw) else "Unknown"

    # day
    if not is_missing(day_raw):
        # 正常可转 float
        try:
            return stage, float(day_raw)
        except Exception:
            pass

    # day 缺失：兜底策略
    # 1) somatic
    if re.search(r"(HEFs|hASFs|hADSCs)", sample_name, flags=re.IGNORECASE):
        return stage if stage != "Unknown" else "Somatic", -1.0

    # 2) StageI/II/III/IV 但没有 day：用 anchor day
    stage_anchor = {"StageI": 0.0, "StageII": 4.0, "StageIII": 8.0, "StageIV": 12.0}
    if stage in stage_anchor:
        return stage, stage_anchor[stage]

    # 3) 实在不行：Unknown + NaN（不建议，但至少不崩）
    return stage, float("nan")


def main():
    df = pd.read_csv(META, sep="\t")

    required_cols = {"GSM", "sample_name", "day", "matrix_path"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"gsm_qc_summary.tsv missing columns: {missing}")

    # 计算：hCiPSC 的占位 day = (非 hCiPSC 的最大 day) + 9
    non_hcipsc_mask = ~df["sample_name"].astype(str).apply(is_hcipsc_sample)
    max_day = pd.to_numeric(df.loc[non_hcipsc_mask, "day"], errors="coerce").max()
    if pd.isna(max_day):
        # 如果全是 NaN，就给一个保守默认
        max_day = 12.0
    hcipsc_day = float(max_day) + 9.0
    print(f"[INFO] computed hCiPSC day = {hcipsc_day} (max non-hCiPSC day={max_day})")

    for _, row in df.iterrows():
        gsm = str(row["GSM"]).strip()
        sample_name = str(row["sample_name"]).strip()

        # donor/batch：优先从 sample_name 末尾解析 4 位数字
        donor = parse_batch_from_sample_name(sample_name)
        if "donor" in df.columns and not is_missing(row.get("donor", None)):
            donor = str(row["donor"]).strip()

        stage, day = parse_stage_day(row, hcipsc_day)

        tenx_dir = tenx_dir_from_row(row)
        print(f"[INFO] {gsm}  tenx_dir={tenx_dir}  stage={stage}  day={day}  donor={donor}")

        adata = read_10x_mtx(tenx_dir)

        # 固定 raw counts
        adata.layers["counts"] = adata.X.copy()

        # 写 obs 字段
        adata.obs["gsm"] = gsm
        adata.obs["sample_name"] = sample_name
        adata.obs["stage"] = str(stage)
        adata.obs["day"] = float(day)
        adata.obs["donor"] = str(donor)
        adata.obs["dataset"] = DATASET_TAG

        # cell id 全局唯一
        adata.obs_names = [f"{gsm}_{bc}" for bc in adata.obs_names]

        out = OUTDIR / f"{gsm}.h5ad"
        adata.write_h5ad(out)
        print(f"[OK] wrote {out} shape={adata.shape}")

    print(f"\n[DONE] All GSM processed. Output dir: {OUTDIR}")


if __name__ == "__main__":
    main()
