"""
src — annotation-robustness pipeline for the trajectory project.

Modular, reproducible code for the three-month MSc slice:
    data -> preprocessing -> annotation (marker-silver vs scANVI/surrogate)
         -> trajectory -> evaluation (annotation comparison) -> plotting

Heavy dependencies (scanpy, scvi-tools, torch) are imported lazily so the
package is importable and the smoke test runs on a CPU-only laptop. See
config.yaml and research_design_three_month_project.md.
"""

__version__ = "0.1.0"
