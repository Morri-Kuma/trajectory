# Notebooks

Thin, runnable drivers over `src/` mirroring the suggested notebook sequence.
They are plain `.py` so they run anywhere and stay diff-friendly; open them as
notebooks via jupytext if preferred. Run from the repo root, e.g.
`python notebooks/03_scanvi_annotation_test.py`. All default to `mode: test`
(synthetic data); switch `config.yaml` to `mode: server` for real inputs.

- `01_data_overview.py` — load inputs, print shapes/labels.
- `02_preprocessing_test.py` — preprocess + PCA/kNN sanity check.
- `03_scanvi_annotation_test.py` — annotate a query (surrogate) + agreement.
- `04_trajectory_test.py` — pseudotime under both annotations + correlation.
- `05_figure_generation_test.py` — full pipeline + figures (== smoke test).
