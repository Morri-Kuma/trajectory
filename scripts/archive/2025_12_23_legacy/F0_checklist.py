import os
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

h5ad_path = r"E:\scgpt\results\C_traj\20251224_153430\adata\C4_wot_labels_scoreTarget.h5ad"
adata = sc.read_h5ad(h5ad_path)

adata.obs_names_make_unique()

# 这个文件里有 p_iPSC
y = adata.obs["p_iPSC"].to_numpy()

X = adata.X
if not isinstance(X, np.ndarray):
    X = X.toarray()  # 数据很大时会吃内存；先跑通为主

genes = np.array(adata.var_names)

rho = []
for i in range(X.shape[1]):
    r, _ = spearmanr(X[:, i], y)
    rho.append(r)
rho = np.array(rho)

pos = genes[np.argsort(-rho)[:200]]
neg = genes[np.argsort(rho)[:200]]

os.makedirs(r"results\tables", exist_ok=True)
pd.Series(pos).to_csv(r"results\tables\cand_pos_top200.txt", index=False, header=False)
pd.Series(neg).to_csv(r"results\tables\cand_neg_top200.txt", index=False, header=False)
print("[OK] wrote results/tables/cand_pos_top200.txt and cand_neg_top200.txt")
