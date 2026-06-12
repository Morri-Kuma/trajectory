"""04 - trajectory: pseudotime correlation under the two annotations."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src import pipeline
s = pipeline.run("config.yaml")
for q, info in s["queries"].items():
    print(q, "pseudotime spearman=%.3f" % info["pseudotime_spearman"],
          "| graph edge Jaccard=%.3f" % info["graph_edge_jaccard"])
