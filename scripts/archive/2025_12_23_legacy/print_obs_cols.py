import scanpy as sc

p = r"E:\scgpt\results\C_traj\20251224_153430\adata\C4_wot_labels_scoreTarget.h5ad"
adata = sc.read_h5ad(p)
print("[obs columns]", list(adata.obs.columns))
