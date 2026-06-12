"""01 - data overview: load inputs and summarize."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.utils import load_config, resolve_path
from src import pipeline
cfg = load_config("config.yaml")
inp = pipeline._load_inputs(cfg)
for name, a in {"reference": inp["reference"], **inp["queries"]}.items():
    print(name, a.shape, "| obs:", [c for c in a.obs.columns][:8])
