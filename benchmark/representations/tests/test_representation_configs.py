"""
Tests for the generated representation configs and scope constraints.

Guarantees:
  - every generated config carries a valid representation block parseable by
    parse_representation_config and is supplementary (formal_benchmark: false);
  - only the four allowed representation arms appear (no GeneCompass / UCE /
    hybrid / random-projection controls);
  - configs route to the representation evaluators, not the gene-expression one.
"""

import re
from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

from benchmark.representations.model_input import parse_representation_config
from benchmark.representations.references import REPRESENTATION_ARMS

CONFIG_DIR = (
    Path(__file__).resolve().parents[2] / "configs" / "representation"
)
FORBIDDEN_TOKENS = ("genecompass", "uce", "random_projection", "hybrid")


def _config_files():
    return sorted(CONFIG_DIR.glob("*.yaml"))


def test_configs_exist():
    files = _config_files()
    assert len(files) == 36, f"expected 36 configs, found {len(files)}"


@pytest.mark.parametrize("path", _config_files(), ids=lambda p: p.name)
def test_each_config_has_valid_representation_block(path):
    with open(path, encoding="utf-8-sig") as f:
        cfg = yaml.safe_load(f)

    rep = parse_representation_config(cfg)
    assert rep.enabled is True
    assert rep.input_mode == "obsm"
    assert rep.obsm_key == "X_rep"
    assert rep.final_dim == 50
    assert rep.reducer_fit_scope == "train_only"
    assert rep.representation_id in REPRESENTATION_ARMS

    assert cfg.get("formal_benchmark") is False
    assert cfg.get("result_class") == "representation_dynamics_supplementary"
    assert cfg["evaluation"]["mode"] == "representation"
    assert "eval_representation_forecast" in cfg["evaluation"]["forecast_eval_script"]


@pytest.mark.parametrize("path", _config_files(), ids=lambda p: p.name)
def test_no_forbidden_representations(path):
    # Word-boundary match so legitimate substrings (e.g. the "uce" inside
    # "reducer_fit_scope") don't trip the scope guard.
    text = path.read_text(encoding="utf-8-sig").lower()
    for token in FORBIDDEN_TOKENS:
        assert not re.search(rf"\b{re.escape(token)}\b", text), (
            f"{path.name} references forbidden arm {token!r}"
        )
