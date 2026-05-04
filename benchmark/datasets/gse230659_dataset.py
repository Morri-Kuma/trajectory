"""
gse230659_dataset.py
Dataset class for GSE230659 (human chemical iPSC reprogramming, Liuyang 2023).
Framework reference: experimental framework v2.md §5

This module provides a pickle-able dataset object that loads
``data/processed/adata_benchmark.h5ad`` and returns a (train, test) split
suitable for the scTimeBench benchmark pipeline.

Scenario A — Observed-time smoke test (current implementation)
--------------------------------------------------------------
For the first runnable WOT smoke test we use the simplest valid split:

  train  : all cells across all 15 observed time points
            (WOT needs every consecutive time-pair to infer transport maps)
  test   : last observed time point only (abs_day == 30.0, i.e. hCiPSC day 30)
            This is the terminal state that Lineage Fidelity evaluation targets.

NOTE: returning ``train == full_adata`` and ``test == day-30 subset`` is a
scaffold-level split.  Once the reference lineage graph is frozen (v2 §9.3),
the test split will be replaced with the evaluation harness that compares
predicted state transitions against the reference graph.  Until then, the
``test`` object is provided purely to allow the pipeline to be end-to-end
runnable without raising errors at the WOT run step.

Pickle compatibility
--------------------
The class stores only paths and configuration at init time.  AnnData objects
are loaded on-demand inside ``load_data()`` and are NOT stored on ``self``.
This keeps the pickled object small (a few hundred bytes) and avoids h5py
file-handle issues when the pickle is re-opened in a new process.
"""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Tuple, Optional


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATASET_ID       = "GSE230659"
TIME_KEY         = "time_label"        # obs column holding numeric abs_day
TERMINAL_TIMEPOINT = 30.0              # hCiPSC day — used as the test split
SCENARIO_ID      = "A"                 # Observed-time interpolation (scaffold)


# ---------------------------------------------------------------------------
# Dataset class
# ---------------------------------------------------------------------------

class GSE230659Dataset:
    """
    Pickle-able dataset for GSE230659.

    Parameters
    ----------
    h5ad_path : str or Path
        Absolute path to ``adata_benchmark.h5ad``.
        Defaults to the canonical location under the project root.
    scenario : str
        Benchmark scenario ID (currently only ``"A"`` is implemented).
    project_root : str or Path, optional
        Project root directory.  Auto-detected from ``h5ad_path`` if not given.

    Attributes
    ----------
    dataset_id : str
        Always ``"GSE230659"``.
    scenario : str
        Scenario ID this instance was built for.
    h5ad_path : Path
        Absolute path to the source h5ad file.
    time_key : str
        obs column used as the benchmark time axis.
    terminal_timepoint : float
        The terminal time point used as the ``test`` split.
    """

    dataset_id: str = DATASET_ID

    def __init__(
        self,
        h5ad_path: Optional[str | Path] = None,
        scenario: str = SCENARIO_ID,
        project_root: Optional[str | Path] = None,
    ):
        # Resolve project root
        if project_root is not None:
            self._project_root = Path(project_root)
        else:
            # Auto-detect project root.
            # 1. Explicit override: TRAJ_PROJECT_ROOT env var (Shirokane HPC / CI).
            import os
            env = os.environ.get("TRAJ_PROJECT_ROOT")
            if env:
                p = Path(env)
                if not p.exists():
                    raise FileNotFoundError(
                        f"TRAJ_PROJECT_ROOT={env!r} does not exist.  "
                        "Correct the environment variable and retry."
                    )
                self._project_root = p
            else:
                # 2. Walk upward from this file until a 'data/' directory is found.
                here = Path(__file__).resolve().parent
                self._project_root = here.parent  # fallback
                for candidate in [here, *here.parents]:
                    if (candidate / "data").exists():
                        self._project_root = candidate
                        break

        # Resolve h5ad path
        if h5ad_path is not None:
            self.h5ad_path = Path(h5ad_path)
        else:
            self.h5ad_path = (
                self._project_root / "data" / "processed" / "adata_benchmark.h5ad"
            )

        self.scenario           = scenario
        self.time_key           = TIME_KEY
        self.terminal_timepoint = TERMINAL_TIMEPOINT

        # Validate the h5ad exists at construction time (fast check, no read).
        if not self.h5ad_path.exists():
            raise FileNotFoundError(
                f"adata_benchmark.h5ad not found at {self.h5ad_path}.\n"
                "Run scripts/01_load_and_qc.py first."
            )

    # ------------------------------------------------------------------
    # Pickle protocol — only serialize lightweight attributes
    # ------------------------------------------------------------------

    def __getstate__(self):
        return {
            "_project_root":    str(self._project_root),
            "h5ad_path":        str(self.h5ad_path),
            "scenario":         self.scenario,
            "time_key":         self.time_key,
            "terminal_timepoint": self.terminal_timepoint,
        }

    def __setstate__(self, state):
        self._project_root    = Path(state["_project_root"])
        self.h5ad_path        = Path(state["h5ad_path"])
        self.scenario         = state["scenario"]
        self.time_key         = state["time_key"]
        self.terminal_timepoint = state["terminal_timepoint"]

    # ------------------------------------------------------------------
    # Main interface
    # ------------------------------------------------------------------

    def load_data(self):
        """
        Load the benchmark h5ad and return a (train, test) AnnData pair.

        Returns
        -------
        train_adata : AnnData
            All 75,194 cells across all 15 observed time points.
            This is the input to WOT transport map inference.

        test_adata : AnnData
            Cells at the terminal time point only (abs_day == 30.0, n ≈ 9,142).

        Notes
        -----
        Scaffold split (scenario A smoke test):
            train = full dataset    — WOT needs all consecutive time pairs
            test  = day-30 subset   — terminal state; will be used for
                                      Lineage Fidelity evaluation once a
                                      reference lineage graph is frozen.

        Once the reference graph is available, ``test_adata`` should carry
        the reference lineage edges rather than (or in addition to) raw cells.
        """
        import anndata as ad

        print(f"[GSE230659Dataset] Loading {self.h5ad_path} ...")
        adata = ad.read_h5ad(self.h5ad_path)
        print(f"  Loaded: {adata.n_obs} cells × {adata.n_vars} genes")
        print(f"  Time points: {sorted(adata.obs[self.time_key].unique())}")

        # train: full dataset (all time points — required by WOT)
        train_adata = adata.copy()

        # test: terminal time point only
        # SCAFFOLD NOTE: this is not a held-out "unseen" split in the generative
        # sense. For Lineage Fidelity (the only active WOT dimension), the
        # evaluation compares predicted state transitions against a reference
        # lineage, not against these held-out cells.  The test_adata object is
        # provided here so the pipeline has a concrete second object to pass
        # around; it will be replaced once the reference graph is frozen.
        terminal_mask = adata.obs[self.time_key] == self.terminal_timepoint
        test_adata    = adata[terminal_mask].copy()
        print(f"  train: {train_adata.n_obs} cells | test (day {self.terminal_timepoint:.0f}): {test_adata.n_obs} cells")

        # Annotate uns with split provenance
        train_adata.uns["split_role"]    = "train"
        train_adata.uns["split_scenario"] = self.scenario
        test_adata.uns["split_role"]     = "test_scaffold"
        test_adata.uns["split_scenario"] = self.scenario
        test_adata.uns["split_note"]     = (
            "Scaffold split: day-30 terminal cells. Replace with reference "
            "lineage graph edges when the benchmark reference is frozen."
        )

        return train_adata, test_adata

    # ------------------------------------------------------------------
    # Convenience helpers
    # ------------------------------------------------------------------

    def time_points(self) -> list:
        """Return the sorted list of observed time points without loading X."""
        import anndata as ad
        # Read obs only (fast path — avoids loading the full expression matrix)
        adata = ad.read_h5ad(self.h5ad_path, backed="r")
        times = sorted(adata.obs[self.time_key].unique().tolist())
        adata.file.close()
        return times

    def cell_counts_per_timepoint(self) -> "pd.Series":
        """Return a Series of cell counts indexed by time_label."""
        import anndata as ad
        adata = ad.read_h5ad(self.h5ad_path, backed="r")
        counts = adata.obs[self.time_key].value_counts().sort_index()
        adata.file.close()
        return counts

    def __repr__(self):
        return (
            f"GSE230659Dataset("
            f"scenario={self.scenario!r}, "
            f"h5ad={self.h5ad_path.name!r}, "
            f"terminal={self.terminal_timepoint})"
        )
