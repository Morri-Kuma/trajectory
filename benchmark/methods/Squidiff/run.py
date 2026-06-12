"""benchmark/methods/Squidiff/run.py — Squidiff benchmark runner.

Squidiff (He, Zhu et al. 2026, Nat Methods 23(1):65) is an encoder-conditioned guided
diffusion model (OpenAI improved-diffusion lineage). Install: `pip install Squidiff==1.0.8`
(import name is `Squidiff`, capital S). It trains from an on-disk h5ad whose obs['Group']
holds the conditioning class and X holds the features, and samples with
diffusion.p_sample_loop. It has no native gene decoder, so we diffuse in PCA space and
inverse-transform projected PCs back to expression (the MIOFlow/PRESCIENT convention).

Model-agnostic Forecast/Embedding/Lineage artifact construction lives in
benchmark/methods/_generative_common.py. This file implements the Squidiff engine. The
PCA round-trip is real; the train/sample calls target Squidiff's real API
(create_model_and_diffusion / run_training / p_sample_loop).

MODELLING NOTE (VERIFY on GPU): Squidiff was designed for control->perturbed prediction
with an encoder over a start cell + a Group label. We map that onto timepoint progression
by setting Group = timepoint index and, at sampling, conditioning on x_start = a t0 cell
and Group = the target timepoint. Confirm this mapping (and the x_start injection) gives
sensible forecasts on a smoke-test before the full array; it is the one genuine modelling
choice here.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_MOD = Path(__file__).parent / "Squidiff_module"
if _MOD.exists() and str(_MOD) not in sys.path:
    sys.path.insert(0, str(_MOD))

from benchmark.methods import _generative_common as _common  # noqa: E402

# keys match benchmark/configs/runtime/squiddiff_*.yaml -> squiddiff_params
_DEFAULTS = dict(pca_dims=50, hidden_dim=128, n_epochs=100, diffusion_steps=1000,
                 num_layers=3, batch_size=64, lr=1e-4, noise_schedule="linear", seed=42,
                 n_sim_cells=None, n_sim_cells_cap=2000, metric_sample_cells=1000)


def _g(cfg, k):
    return cfg.get(k, _DEFAULTS[k])


class SquidiffEngine:
    """PCA-space Squidiff (guided-diffusion) engine: embed + project."""

    def __init__(self, n_genes, cfg, cache_path):
        self.n_genes = n_genes
        self.cfg = cfg
        self.cache_path = Path(cache_path)
        self.workdir = self.cache_path.parent / "squidiff_work"
        self.pca = None
        self.scaler = None
        self.tp_index = {}
        self.ckpt = None
        self._model = None
        self._diffusion = None

    # ---- real PCA round-trip -------------------------------------------------
    def _to_pcs(self, X):
        return self.pca.transform(self.scaler.transform(np.asarray(X, np.float32))).astype(np.float32)

    def _from_pcs(self, Z):
        return self.scaler.inverse_transform(self.pca.inverse_transform(np.asarray(Z, np.float32)))

    def fit(self, train_data, train_tps):
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        Xall = np.concatenate([np.asarray(x, np.float32) for x in train_data], axis=0)
        self.scaler = StandardScaler().fit(Xall)
        self.pca = PCA(n_components=int(_g(self.cfg, "pca_dims")),
                       random_state=int(_g(self.cfg, "seed"))).fit(self.scaler.transform(Xall))
        sorted_tps = sorted(float(t) for t in train_tps)
        self.tp_index = {t: i for i, t in enumerate(sorted_tps)}
        self._fit_model(train_data, sorted_tps)
        return self

    # ---- Engine protocol -----------------------------------------------------
    def embed(self, X):
        return self._to_pcs(X)

    def project(self, X0, target_tps, n_sim_cells):
        z0 = self._to_pcs(X0)
        n = int(n_sim_cells) if n_sim_cells else z0.shape[0]
        recon, latent = [], []
        for t in (float(t) for t in target_tps):
            z_t = self._sample_pcs(z0, t, n)
            latent.append(z_t)
            recon.append(self._from_pcs(z_t))
        return np.stack(recon, axis=1), np.stack(latent, axis=1)

    # ---- Squidiff real-API hooks --------------------------------------------
    def _squidiff(self):
        try:
            import Squidiff  # capital S (pip install Squidiff==1.0.8)
            return Squidiff
        except Exception as exc:
            raise RuntimeError(
                "Squidiff not importable. `pip install Squidiff==1.0.8` in traj_env "
                "(import name is `Squidiff`, capital S). " + repr(exc))

    def _fit_model(self, train_data, sorted_tps):
        """Train Squidiff in PCA space via its on-disk run_training entrypoint.
        Builds a temp h5ad: X = PCA(cells), obs['Group'] = timepoint index."""
        import anndata as ad
        import pandas as pd
        Sq = self._squidiff()
        from Squidiff import dist_util
        from Squidiff.train_squidiff import run_training

        self.workdir.mkdir(parents=True, exist_ok=True)
        Z = np.concatenate([self._to_pcs(x) for x in train_data], axis=0).astype(np.float32)
        groups = np.concatenate([[self.tp_index[t]] * len(train_data[i])
                                 for i, t in enumerate(sorted_tps)]).astype(int)
        adata = ad.AnnData(X=Z, obs=pd.DataFrame({"Group": groups}))
        h5 = self.workdir / "train_pca.h5ad"
        adata.write_h5ad(str(h5))

        n_pcs = int(_g(self.cfg, "pca_dims"))
        bs = int(_g(self.cfg, "batch_size"))
        steps_per_epoch = max(1, -(-Z.shape[0] // bs))   # ceil(n_cells / batch)
        train_steps = max(1, int(_g(self.cfg, "n_epochs")) * steps_per_epoch)
        args = dict(Sq.model_and_diffusion_defaults())
        args.update(dict(
            data_path=str(h5), logger_path=str(self.workdir),
            gene_size=n_pcs, output_dim=n_pcs, num_layers=int(_g(self.cfg, "num_layers")),
            num_channels=int(_g(self.cfg, "hidden_dim")),
            class_cond=False, use_encoder=True, use_drug_structure=False, comb_num=1,
            diffusion_steps=int(_g(self.cfg, "diffusion_steps")),
            noise_schedule=_g(self.cfg, "noise_schedule"),
            batch_size=bs, lr=float(_g(self.cfg, "lr")),
            lr_anneal_steps=train_steps,
            microbatch=-1, ema_rate="0.9999", schedule_sampler="uniform",
            log_interval=max(1, train_steps // 5), save_interval=train_steps,
            resume_checkpoint="", use_fp16=False, fp16_scale_growth=1e-3, weight_decay=0.0,
        ))
        dist_util.setup_dist()
        run_training(args)
        ckpts = sorted(self.workdir.glob("model*.pt"))
        if not ckpts:
            raise RuntimeError(f"Squidiff training wrote no model*.pt under {self.workdir}")
        self.ckpt = ckpts[-1]

    def _load(self):
        if self._model is not None:
            return
        import torch as th
        Sq = self._squidiff()
        from Squidiff import dist_util
        n_pcs = int(_g(self.cfg, "pca_dims"))
        args = dict(Sq.model_and_diffusion_defaults())
        args.update(dict(gene_size=n_pcs, output_dim=n_pcs,
                         num_layers=int(_g(self.cfg, "num_layers")),
                         num_channels=int(_g(self.cfg, "hidden_dim")),
                         class_cond=False, use_encoder=True,
                         diffusion_steps=int(_g(self.cfg, "diffusion_steps")),
                         noise_schedule=_g(self.cfg, "noise_schedule")))
        model, diffusion = Sq.create_model_and_diffusion(
            **{k: args[k] for k in Sq.model_and_diffusion_defaults().keys()})
        model.load_state_dict(th.load(str(self.ckpt), map_location=dist_util.dev()))
        model.to(dist_util.dev()); model.eval()
        self._model, self._diffusion = model, diffusion

    def _sample_pcs(self, z0, target_tp, n):
        """Conditional reverse-diffusion sample of n PCA cells at target_tp.
        VERIFY: conditions on x_start = resampled t0 cells + Group = tp index."""
        import torch as th
        from Squidiff import dist_util
        self._load()
        n_pcs = int(_g(self.cfg, "pca_dims"))
        idx = np.random.default_rng(int(_g(self.cfg, "seed"))).integers(0, z0.shape[0], n)
        x_start = th.tensor(z0[idx], dtype=th.float32, device=dist_util.dev())
        grp = int(self.tp_index.get(target_tp, min(self.tp_index.values(),
                  key=lambda i: abs(i - target_tp))))
        model_kwargs = {"x_start": x_start,
                        "group": th.tensor([grp] * n, device=dist_util.dev())}
        with th.no_grad():
            sample = self._diffusion.p_sample_loop(
                self._model, (n, n_pcs), model_kwargs=model_kwargs, clip_denoised=False)
        return sample.detach().cpu().numpy()


# --------------------------------------------------------------------------- #
# 5-function contract
# --------------------------------------------------------------------------- #
def prepare_data(adata, time_key, train_times=None):
    return _common.prepare_data(adata, time_key, train_times)


def train_or_load(train_data, train_tps, n_genes, cfg, model_cache):
    return SquidiffEngine(n_genes, cfg or {}, model_cache).fit(train_data, train_tps)


def run_forecast_accuracy(model, **kw):
    return _common.run_forecast_accuracy(model, **kw)


def run_embedding_coherence(model, **kw):
    return _common.run_embedding_coherence(model, **kw)


def run_lineage_fidelity(model, **kw):
    return _common.run_lineage_fidelity(model, **kw)
