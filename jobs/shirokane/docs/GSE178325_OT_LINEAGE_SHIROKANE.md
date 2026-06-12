# GSE178325 WOT/CellRank2 Lineage Run on Shirokane

This project-level job adds the missing GSE178325 scenario-A official-silver
lineage runs for WOT and CellRank2. It is meant to be run directly from the
current `trajectory` project checkout on Shirokane. No separate zip bundle is
required.

## What It Runs

- Task 1: `wot_gse178325_marker_fm_silver_A_hvg2000_formal`
- Task 2: `cellrank2_gse178325_marker_fm_silver_A_hvg2000_formal`

Both tasks use:

- input: `benchmark/inputs/gse178325_marker_fm_transition_silver_hvg2000/GSE178325_marker_fm_transition_silver_HVG2000_benchmark_input.h5ad`
- state key: `final_milestone_label_coarse`
- reference graph: `benchmark/ground_truth/providers/gse178325_marker_fm_transition_silver_v1/reference_graph.json`
- output dirs:
  - `benchmark/results/wot/gse178325_marker_fm_silver_A_hvg2000_formal`
  - `benchmark/results/cellrank2/gse178325_marker_fm_silver_A_hvg2000_formal`

## Required Project Files

Make sure the Shirokane checkout contains these files:

- `benchmark/configs/runtime/wot_gse178325_marker_fm_silver_A_hvg2000_formal.yaml`
- `benchmark/configs/runtime/cellrank2_gse178325_marker_fm_silver_A_hvg2000_formal.yaml`
- `jobs/shirokane/run_gse178325_ot_lineage_array.sh`
- `jobs/shirokane/run_gse178325_ot_lineage_summary.sh`
- `jobs/shirokane/submit_gse178325_ot_lineage.sh`
- `scripts/update_gse178325_ot_manifest.py`

## Submit Directly

```bash
cd /home/xzy0723/projects/trajectory
bash jobs/shirokane/submit_gse178325_ot_lineage.sh
```

The submit script first runs:

```bash
python scripts/update_gse178325_ot_manifest.py
```

so `benchmark/results/result_manifest.yaml` will collect the new WOT/CellRank2
GSE178325-A result directories when the official summary is regenerated.

For more memory:

```bash
S_VMEM=256G bash jobs/shirokane/submit_gse178325_ot_lineage.sh
```

Run only WOT first:

```bash
TASK_RANGE=1 RUN_SUMMARY=0 bash jobs/shirokane/submit_gse178325_ot_lineage.sh
```

Run only CellRank2:

```bash
TASK_RANGE=2 RUN_SUMMARY=0 bash jobs/shirokane/submit_gse178325_ot_lineage.sh
```

## Monitor

```bash
qstat -u "$USER"
ls logs/run_gse178325_ot_lineage.*.log
tail -n 80 logs/run_gse178325_ot_lineage.<JOB_ID>.<TASK_ID>.log
```

The submit script also queues a held summary job by default. It validates the
two new outputs and then runs:

```bash
python -m benchmark.evaluation.summarize_official_silver \
  --leiden-n-neighbors 15 \
  --leiden-resolution 0.5
```

## Expected Output Files

Each method should write:

- `state_transition_matrix.csv`
- `lineage_graph_edges.csv`
- `lineage_metrics.json`
- `run_metadata.json`

The summary job writes:

- `benchmark/reports/official_silver/gse178325_ot_lineage_summary.csv`
- updated `benchmark/reports/official_silver/official_silver_lineage_rankings.csv`
- updated `benchmark/reports/official_silver/official_silver_lineage_method_summary.csv`
- updated `benchmark/reports/official_silver/official_silver_model_rankings.md`

## Failure Conditions Checked

The run script fails if:

- the input h5ad is missing;
- `final_milestone_label_coarse` is missing;
- h5ad/reference/method outputs contain `stage_*` labels;
- the adapter falls back to scaffold outputs;
- `lineage_metrics.json` status is not `completed`;
- `prediction_label_source` is not `reference_graph_nodes_sctimebench_zero_filled`.
