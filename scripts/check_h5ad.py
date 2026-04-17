import scanpy as sc
import pandas as pd

# 读取 h5ad 文件
file_path = r"C:\Users\37620\Documents\GitHub\scTimeBench\data\human_germ.h5ad"
adata = sc.read_h5ad(file_path)

# 基本信息
print("adata shape:", adata.shape)
print("\nobs 前10行：")
print(adata.obs.head(10))

print("\nvar 前10行：")
print(adata.var.head(10))

# 如果想看表达矩阵前10行前10列
print("\nX 前10行前10列：")
print(pd.DataFrame(adata.X[:10, :10].toarray() if hasattr(adata.X, "toarray") else adata.X[:10, :10]))