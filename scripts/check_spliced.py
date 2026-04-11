from pathlib import Path
import scanpy as sc

root = Path(r"C:\Users\37620\trajectory\data\processed")

# 递归查找所有 h5ad 文件
h5ad_files = list(root.rglob("*.h5ad"))

print("找到的 .h5ad 文件：")
for f in h5ad_files:
    print(" -", f)

if not h5ad_files:
    print("\n这个目录下没有找到 .h5ad 文件。")
    print("你可以先看看里面是不是只有 matrix.mtx.gz / barcodes.tsv.gz / features.tsv.gz 这类原始文件。")
else:
    for f in h5ad_files:
        print("\n" + "=" * 80)
        print("文件:", f)

        adata = sc.read_h5ad(f)

        layer_keys = list(adata.layers.keys())
        obs_cols = list(adata.obs.columns)

        print("adata shape:", adata.shape)
        print("layers:", layer_keys)
        print("obs columns (前50个):", obs_cols[:50])

        print("是否有 spliced:", "spliced" in layer_keys)
        print("是否有 unspliced:", "unspliced" in layer_keys)

        # 顺手检查一下有没有可用的聚类列
        candidate_cluster_cols = [c for c in obs_cols if "cluster" in c.lower() or "leiden" in c.lower() or "louvain" in c.lower()]
        print("可能的聚类列:", candidate_cluster_cols)