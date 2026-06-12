"""02 - preprocessing: normalize -> PCA -> kNN."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.utils import load_config
from src import pipeline
from src.preprocessing import basic_preprocess
cfg = load_config("config.yaml")
q = list(pipeline._load_inputs(cfg)["queries"].values())[0]
q = basic_preprocess(q, n_pcs=30, n_neighbors=15, max_genes=2000, seed=0)
print("X_pca:", q.obsm["X_pca"].shape, "| backend:", q.uns["preprocess"]["backend"])
