"""03 - annotation: surrogate scANVI + agreement on one query (via pipeline)."""
import os, sys, json; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src import pipeline
s = pipeline.run("config.yaml")
for q, info in s["queries"].items():
    print(q, "NMI=%.3f" % info["agreement_nmi"], "corr_agree=%.3f" % (info["correspondence_overall_agreement"] or float("nan")))
