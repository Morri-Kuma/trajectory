"""05 - figures: run the full pipeline and list generated figures."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src import pipeline
s = pipeline.run("config.yaml")
for q, info in s["queries"].items():
    print(q, "figures:", list(info["figures"].values()))
