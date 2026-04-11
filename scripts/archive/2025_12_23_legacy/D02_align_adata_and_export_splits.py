# e:/scgpt/scripts/2025_12_23/D02_align_adata_and_export_splits.py
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Tuple, Set

import numpy as np
import pandas as pd
import scanpy as sc

#
@dataclass
class Cfg:
    # D01 输出目录（你当前这次 run）
    run_dir: Path = Path(r"E:\scgpt\results\D_finetune\run_20251228_1655")

    # C 阶段的带标签 adata
    adata_path: Path = Path(
        r"E:\scgpt\results\C_traj\20251224_153430\adata\C4_wot_labels_scoreTarget.h5ad"
    )

    # 预训练模型 vocab（whole-human）
    vocab_json: Path = Path(r"E:\scgpt\models\scgpt_whole_human\vocab.json")

    # 输出文件名
    out_adata_name: str = "adata_aligned_vocab.h5ad"


def read_vocab_genes(vocab_json: Path) -> Set[str]:
    """
    兼容三种常见 vocab.json:
      1) {gene: id}
      2) {"itos":[...]}
      3) {"stoi":{...}}
    """
    obj = json.loads(Path(vocab_json).read_text(encoding="utf-8"))
    if not isinstance(obj, dict):
        raise ValueError(f"vocab.json must be a dict: {vocab_json}")

    if len(obj) > 0 and all(isinstance(v, int) for v in obj.values()):
        genes = set(obj.keys())
    elif "itos" in obj and isinstance(obj["itos"], list):
        genes = set(obj["itos"])
    elif "stoi" in obj and isinstance(obj["stoi"], dict):
        genes = set(obj["stoi"].keys())
    else:
        raise ValueError(f"Unrecognized vocab structure: {vocab_json}")

    # 去掉 special tokens（<pad> 等）
    genes = {g for g in genes if not (isinstance(g, str) and g.startswith("<") and g.endswith(">"))}
    return genes


def load_label_keys(path: Path) -> Dict[str, str]:
    """
    label_keys.json 允许你以后改键名而不改代码。
    若不存在则回退默认键。
    """
    default = {
        "stage_key": "stage",
        "gsm_key": "gsm",
        "day_key": "day_numeric",
        "pseudotime_key": "pseudotime",
        "reg_key": "p_iPSC",
        "target_flag_key": "is_target_stageIII12_top",
    }
    if path.exists():
        obj = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(obj, dict):
            default.update(obj)
    return default


def load_split_indices(adata: sc.AnnData, split_json: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    兼容多种 split.json 格式。
    你当前 D01 的格式是：
      {"train":[cell_id...], "val":[...], "test":[...], ...}
    也兼容：
      {"train_cells":[...]} / {"splits":{"train":[...]}} / 以及直接给 indices 的情况
    """
    sp = json.loads(split_json.read_text(encoding="utf-8"))

    def _pick(name: str):
        # 兼容常见命名
        candidates = [
            f"{name}_cells", f"{name}Cells",
            f"{name}_cell_ids", f"{name}CellIds",
            f"{name}_obs_names", f"{name}ObsNames",
            name,  # 你当前就是 "train"/"val"/"test"
            f"{name}_idx", f"{name}Idx",
            f"{name}_indices", f"{name}Indices",
        ]
        for k in candidates:
            if k in sp:
                return sp[k]

        # 兼容嵌套结构
        for container_key in ["splits", "split", "data"]:
            if container_key in sp and isinstance(sp[container_key], dict) and name in sp[container_key]:
                return sp[container_key][name]

        raise KeyError(
            f"{split_json.name}: cannot find split list for '{name}'. "
            f"Top-level keys={list(sp.keys())[:30]}"
        )

    def _normalize(x):
        # 可能是 dict 包着 list：{"cells":[...]} / {"obs_names":[...]} 之类
        if isinstance(x, dict):
            for kk in ["cells", "cell_ids", "obs_names", "ids", "idx", "indices"]:
                if kk in x:
                    x = x[kk]
                    break

        # 可能 list 元素是 dict：{"cell_id": "..."} 之类
        if isinstance(x, list) and len(x) > 0 and isinstance(x[0], dict):
            for kk in ["cell_id", "obs_name", "id"]:
                if kk in x[0]:
                    x = [d[kk] for d in x]
                    break
        return x

    train_raw = _normalize(_pick("train"))
    val_raw = _normalize(_pick("val"))
    test_raw = _normalize(_pick("test"))

    def _to_indices(raw, split_name: str) -> np.ndarray:
        # 如果已经是整数索引
        if isinstance(raw, (list, tuple)) and len(raw) > 0 and isinstance(raw[0], (int, np.integer)):
            idx = np.asarray(raw, dtype=np.int32)
        else:
            # 当作 cell_id / obs_name 字符串列表
            raw = [str(x) for x in raw]
            idx = adata.obs_names.get_indexer(raw).astype(np.int32)

        miss = int(np.sum(idx < 0))
        if miss > 0:
            raise ValueError(f"[{split_json.name}] {split_name}: {miss} items not found in adata.obs_names")
        return idx

    idx_train = _to_indices(train_raw, "train")
    idx_val = _to_indices(val_raw, "val")
    idx_test = _to_indices(test_raw, "test")
    return idx_train, idx_val, idx_test


def main():
    cfg = Cfg()

    data_dir = cfg.run_dir / "data"
    splits_dir = cfg.run_dir / "splits"
    data_dir.mkdir(parents=True, exist_ok=True)

    # label key mapping
    keys = load_label_keys(data_dir / "label_keys.json")
    stage_key = keys["stage_key"]
    gsm_key = keys["gsm_key"]
    day_key = keys["day_key"]
    pseudotime_key = keys["pseudotime_key"]
    reg_key = keys["reg_key"]
    target_flag_key = keys["target_flag_key"]

    print(f"[INFO] Reading adata: {cfg.adata_path}")
    adata = sc.read_h5ad(cfg.adata_path)
    print(f"[INFO] Loaded: n_cells={adata.n_obs}, n_genes={adata.n_vars}")

    # check required columns
    need_cols = [stage_key, gsm_key, day_key, pseudotime_key, reg_key]
    for c in need_cols:
        if c not in adata.obs.columns:
            raise KeyError(f"adata.obs missing required column: {c}")

    # vocab align
    print(f"[INFO] Reading vocab: {cfg.vocab_json}")
    vocab_genes = read_vocab_genes(cfg.vocab_json)

    keep_genes = [g for g in adata.var_names.astype(str).tolist() if g in vocab_genes]
    adata2 = adata[:, keep_genes].copy()
    print(f"[INFO] Gene alignment: kept {len(keep_genes)} / {adata.n_vars} genes")

    # export gene list and vocab copy
    (data_dir / "genes_kept_ordered.txt").write_text("\n".join(keep_genes) + "\n", encoding="utf-8")
    (data_dir / "vocab.json").write_text(cfg.vocab_json.read_text(encoding="utf-8"), encoding="utf-8")

    # encode stage -> int
    stages = adata2.obs[stage_key].astype(str)
    stage_categories = sorted(stages.unique().tolist())
    stage2id = {s: i for i, s in enumerate(stage_categories)}
    id2stage = {i: s for s, i in stage2id.items()}
    y_stage = stages.map(stage2id).astype(np.int16).values

    (data_dir / "stage_label_map.json").write_text(
        json.dumps({"stage2id": stage2id, "id2stage": id2stage}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    # export label arrays
    y_reg = adata2.obs[reg_key].astype(float).values.astype(np.float32)
    y_time = adata2.obs[pseudotime_key].astype(float).values.astype(np.float32)

    np.savez_compressed(
        data_dir / "labels_arrays.npz",
        y_stage=y_stage,
        y_reg=y_reg,
        y_pseudotime=y_time,
    )

    # load split json and convert to index arrays
    gsm_json = splits_dir / "gsm_split.json"
    time_json = splits_dir / "time_split.json"
    if not gsm_json.exists() or not time_json.exists():
        raise FileNotFoundError("Missing gsm_split.json or time_split.json in splits/")

    gsm_train, gsm_val, gsm_test = load_split_indices(adata2, gsm_json)
    time_train, time_val, time_test = load_split_indices(adata2, time_json)

    np.savez_compressed(splits_dir / "gsm_split_indices.npz", train=gsm_train, val=gsm_val, test=gsm_test)
    np.savez_compressed(splits_dir / "time_split_indices.npz", train=time_train, val=time_val, test=time_test)
    print(f"[OK] Wrote: {splits_dir / 'gsm_split_indices.npz'}")
    print(f"[OK] Wrote: {splits_dir / 'time_split_indices.npz'}")

    # meta table for reproducibility/paper
    def mark_split(n: int, idx_train: np.ndarray, idx_val: np.ndarray, idx_test: np.ndarray) -> np.ndarray:
        arr = np.array(["none"] * n, dtype=object)
        arr[idx_train] = "train"
        arr[idx_val] = "val"
        arr[idx_test] = "test"
        return arr

    meta = pd.DataFrame(
        {
            "cell_id": adata2.obs_names.astype(str),
            "gsm": adata2.obs[gsm_key].astype(str).values,
            "stage": adata2.obs[stage_key].astype(str).values,
            "day_numeric": adata2.obs[day_key].astype(float).values,
            "pseudotime": y_time,
            "p_iPSC": y_reg,
            "split_gsm": mark_split(adata2.n_obs, gsm_train, gsm_val, gsm_test),
            "split_time": mark_split(adata2.n_obs, time_train, time_val, time_test),
        }
    )
    if target_flag_key in adata2.obs.columns:
        meta[target_flag_key] = adata2.obs[target_flag_key].astype(bool).values

    meta_path = data_dir / "meta_cells.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)
    print(f"[OK] Wrote: {meta_path}")

    # save aligned adata
    out_adata = data_dir / cfg.out_adata_name
    adata2.write_h5ad(out_adata)
    print(f"[OK] Wrote aligned adata: {out_adata}")

    print("\n[DONE] D02 finished.")


if __name__ == "__main__":
    main()
