# MIOFlow

This directory contains the local vendored MIOFlow implementation used by the
trajectory benchmark, plus the project-specific runner at `run.py`.

Upstream MIOFlow models single-cell time-series data with neural ODEs and
optimal transport. In this benchmark it is treated as a generative/projected-cell
method and is eligible for:

- Forecast Accuracy
- Embedding Coherence
- Lineage Fidelity

Capability flags are declared in
`benchmark/configs/method_capabilities.yaml`.

## Local Installation

From the repository root:

```bash
pip install -e benchmark/methods/MIOFlow
```

Or install the upstream package from PyPI:

```bash
pip install MIOFlow
```

The local `pyproject.toml` declares the package dependencies used by this
vendored copy. Benchmark jobs are normally launched through the repository-level
run scripts rather than upstream tutorials.

## Benchmark Usage

Run a single config from the repository root:

```bash
python benchmark/methods/MIOFlow/run.py \
    --config benchmark/configs/mioflow_gse242424_oskm_silver_A_hvg2000_formal.yaml
```

Marker-FM transition silver runtime configs live under
`benchmark/configs/runtime/`, for example:

```bash
python benchmark/methods/MIOFlow/run.py \
    --config benchmark/configs/runtime/mioflow_gse230659_marker_fm_silver_A_hvg2000_formal.yaml
```

Array-style benchmark runs are usually started from the project root with helper
scripts such as:

```bash
bash run_marker_fm_transition_silver_primary_array.sh
bash run_gse242424_oskm_silver_formal_array.sh
```

## Output Contract

The benchmark runner writes method outputs under
`benchmark/results/mioflow/<run_id>/`. Depending on the config, expected files
include:

- `run_metadata.json`
- `forecast_metrics.json`
- `embedding_metrics.json`
- `lineage_metrics.json`
- `projected_expression.npy`
- `projected_embedding.npy`
- `state_transition_matrix.csv`
- `lineage_graph_edges.csv`

Not every historical run contains every file. Interpret missing files through
the method capability flags and the run metadata.

## Upstream Project

- Upstream repository: https://github.com/KrishnaswamyLab/MIOFlow
- Paper: https://arxiv.org/abs/2206.14928
- PyPI: https://pypi.org/project/mioflow/

## Citation

If you use MIOFlow in research, cite:

```bibtex
@misc{https://doi.org/10.48550/arxiv.2206.14928,
  doi = {10.48550/ARXIV.2206.14928},
  url = {https://arxiv.org/abs/2206.14928},
  author = {Huguet, Guillaume and Magruder, D. S. and Tong, Alexander and Fasina, Oluwadamilola and Kuchroo, Manik and Wolf, Guy and Krishnaswamy, Smita},
  title = {Manifold Interpolating Optimal-Transport Flows for Trajectory Inference},
  publisher = {arXiv},
  year = {2022}
}
```

## License

This vendored copy follows the license included in `LICENSE.md`.
