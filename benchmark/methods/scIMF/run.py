"""benchmark/methods/scIMF/run.py — scIMF benchmark runner.

scIMF (Jiang, Li et al. 2026, PLOS Comput Biol 22(1):e1013916) — Transformer + interacting mean-field neural SDE.

The model-AGNOSTIC Forecast/Embedding/Lineage artifact construction lives in
benchmark/methods/_generative_common.py. This file implements only the scIMF engine:
data plumbing + caching are real; the neural-SDE train / simulate / encode calls are
confined to _fit_model() / _simulate() / _encode() and MUST be confirmed against the
vendored source under benchmark/methods/scIMF/scIMF_module/ (git clone the authors' repo). They
fail loudly with precise guidance rather than emit anything fabricated.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_MOD = Path(__file__).parent / "scIMF_module"
if _MOD.exists() and str(_MOD) not in sys.path:
    sys.path.insert(0, str(_MOD))

from benchmark.methods import _generative_common as _common  # noqa: E402

_DEFAULTS = dict(latent_dim=50, epochs=10, iters=100, pretrain_iters=200, batch_size=32,
                 lr=1e-3, seed=42, n_sim_cells=None, n_sim_cells_cap=2000,
                 metric_sample_cells=1000, use_cuda=True)


def _g(cfg, k):
    return cfg.get(k, _DEFAULTS[k])


class ScIMFEngine:
    """Neural-SDE engine implementing the _generative_common Engine protocol.
    The model trains on per-timepoint expression and can (a) encode cells to a
    latent and (b) simulate the SDE forward from t0 to arbitrary target tps."""

    def __init__(self, n_genes, cfg, cache_path):
        self.n_genes = n_genes
        self.cfg = cfg
        self.cache_path = Path(cache_path)
        self.model = None
        self.train_tps = None

    def fit(self, train_data, train_tps):
        self.train_tps = [float(t) for t in train_tps]
        self.model = self._fit_model([np.asarray(x, np.float32) for x in train_data],
                                     self.train_tps)
        return self

    # ---- Engine protocol -----------------------------------------------------
    def embed(self, X):
        return np.asarray(self._encode(np.asarray(X, np.float32)))

    def project(self, X0, target_tps, n_sim_cells):
        target_tps = [float(t) for t in target_tps]
        n = int(n_sim_cells) if n_sim_cells else int(np.asarray(X0).shape[0])
        # _simulate returns gene-space recon (n, T, g) and latent (n, T, d)
        recon, latent = self._simulate(np.asarray(X0, np.float32), target_tps, n)
        return np.asarray(recon), np.asarray(latent)

    # ---- model-specific hooks: CONFIRM against the vendored scimf source ------
    def _require_source(self):
        try:
            import scimf  # noqa: F401
            return
        except Exception:
            pass
        if not _MOD.exists():
            raise RuntimeError(
                "scIMF source not vendored. git clone the authors' repository into "
                "benchmark/methods/scIMF/scIMF_module/ (see benchmark/methods/VENDORING.md), "
                "then wire ScIMFEngine._fit_model/_simulate/_encode to it.")

    def _fit_model(self, data_by_tp, sorted_tps):
        """Train the scIMF SDE on per-timepoint expression. VERIFY against vendored source.
        Template:
            from scimf import Model, train         # confirm names against the repo
            model = Model(n_genes=self.n_genes, latent_dim=_g(self.cfg,'latent_dim'),
                          use_cuda=_g(self.cfg,'use_cuda'))
            train(model, data_by_tp, sorted_tps, epochs=_g(self.cfg,'epochs'),
                  iters=_g(self.cfg,'iters'), lr=_g(self.cfg,'lr'),
                  pretrain_iters=_g(self.cfg,'pretrain_iters'), seed=_g(self.cfg,'seed'))
            return model
        """
        self._require_source()
        raise NotImplementedError(
            "VERIFY: wire ScIMFEngine._fit_model() to the vendored scIMF training API "
            "(template in this method's docstring). data_by_tp + sorted_tps are ready.")

    def _simulate(self, X0, target_tps, n):
        """Simulate the SDE forward from t0 cells X0 to each target tp.
        Must return (recon[n,T,g] in gene space, latent[n,T,d]). VERIFY:
            states = self.model.simulate(x0=X0, t_grid=[t0]+target_tps, n=n)
            recon  = states.gene_space            # decode if model is latent-space
            latent = states.latent
        """
        self._require_source()
        raise NotImplementedError(
            "VERIFY: wire ScIMFEngine._simulate() to the vendored scIMF forward-simulation API "
            f"(target_tps ready; return gene-space recon + latent).")

    def _encode(self, X):
        """Encode cells to the model latent (for centroids/ARI). VERIFY:
            return self.model.encode(X)
        If scIMF has no encoder, fall back to PCA(50) on X (document the choice)."""
        self._require_source()
        raise NotImplementedError(
            "VERIFY: wire ScIMFEngine._encode() to the vendored scIMF encoder, or use a PCA(50) "
            "fallback if the model is gene-space only.")


# --------------------------------------------------------------------------- #
# 5-function contract (delegates to _generative_common)
# --------------------------------------------------------------------------- #
def prepare_data(adata, time_key, train_times=None):
    return _common.prepare_data(adata, time_key, train_times)


def train_or_load(train_data, train_tps, n_genes, cfg, model_cache):
    return ScIMFEngine(n_genes, cfg or {}, model_cache).fit(train_data, train_tps)


def run_forecast_accuracy(model, **kw):
    return _common.run_forecast_accuracy(model, **kw)


def run_embedding_coherence(model, **kw):
    return _common.run_embedding_coherence(model, **kw)


def run_lineage_fidelity(model, **kw):
    return _common.run_lineage_fidelity(model, **kw)
