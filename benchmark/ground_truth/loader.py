"""
Ground-truth provider loader.

The benchmark treats cell-state labels and reference lineage graphs as a
replaceable provider layer. Current reportable runs resolve to registered
ground-truth providers in ``benchmark/ground_truth/registry.yaml`` or to a
``ground_truth`` block directly in a method config.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Optional


@dataclass
class GroundTruthSpec:
    provider_id: str
    state_key: Optional[str] = None
    label_path: Optional[str] = None
    metadata_path: Optional[str] = None
    reference_graph_path: Optional[str] = None
    reference_edges_path: Optional[str] = None
    confidence_mode: str = "all"
    exclude_uncertain_states: bool = False
    status: Optional[str] = None
    n_states: Optional[int] = None
    n_graph_edges: Optional[int] = None
    source_h5ad: Optional[str] = None
    annotation_method: Optional[str] = None
    notes: Optional[str] = None
    label_mode: Optional[str] = None
    analysis_role: Optional[str] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def _project_root() -> Path:
    here = Path(__file__).resolve()
    for candidate in [here.parent, *here.parents]:
        if (candidate / "benchmark").exists() and (
            candidate / "benchmark" / "evaluation"
        ).exists():
            return candidate
    return here.parents[2]


def _read_yaml(path: Path) -> Dict[str, Any]:
    import yaml

    with open(path, encoding="utf-8-sig") as f:
        data = yaml.safe_load(f) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Expected mapping in YAML file: {path}")
    return data


def _registry_path(project_root: Path) -> Path:
    return project_root / "benchmark" / "ground_truth" / "registry.yaml"


def _provider_from_registry(provider_id: str, project_root: Path) -> Dict[str, Any]:
    registry = _read_yaml(_registry_path(project_root))
    providers = registry.get("providers") or {}
    if provider_id not in providers:
        raise KeyError(
            f"Unknown ground-truth provider {provider_id!r}. "
            f"Registered providers: {sorted(providers)}"
        )
    cfg = providers[provider_id] or {}
    if not isinstance(cfg, dict):
        raise ValueError(f"Provider entry {provider_id!r} must be a mapping.")
    return cfg


def _config_ground_truth_config(method_config: Dict[str, Any]) -> Dict[str, Any]:
    """Build a provider-like config from pre-provider lineage/state_system keys."""
    lineage = method_config.get("lineage") or {}
    state_system = method_config.get("state_system") or {}
    version = state_system.get("version") or lineage.get("state_system_version")
    return {
        "provider_id": version or "config_lineage_provider",
        "state_key": lineage.get("cell_state_key"),
        "label_path": state_system.get("label_tsv"),
        "metadata_path": state_system.get("metadata_tsv"),
        "reference_graph_path": lineage.get("reference_graph_path")
        or state_system.get("graph_json"),
        "reference_edges_path": state_system.get("graph_csv"),
        "confidence_mode": lineage.get("edge_confidence_mode", "all"),
        "exclude_uncertain_states": bool(lineage.get("exclude_uncertain_states", False)),
        "status": state_system.get("status"),
        "n_states": state_system.get("n_states"),
        "n_graph_edges": state_system.get("n_graph_edges"),
        "source_h5ad": state_system.get("source_h5ad"),
        "notes": state_system.get("note"),
    }


def _merge_provider_config(
    provider_cfg: Dict[str, Any],
    override_cfg: Dict[str, Any],
) -> Dict[str, Any]:
    merged = dict(provider_cfg)
    for key, value in override_cfg.items():
        if value is not None:
            merged[key] = value
    return merged


def load_ground_truth(
    method_config: Optional[Dict[str, Any]] = None,
    provider_id: Optional[str] = None,
    project_root: Optional[Path | str] = None,
) -> GroundTruthSpec:
    """
    Resolve a ground-truth provider into a uniform spec.

    Resolution order:
      1. ``method_config["ground_truth"]`` if present.
      2. Registry entry named by ``provider_id`` or ``ground_truth.provider_id``.
      3. Config-local ``lineage`` + ``state_system`` fields.

    Config-local fields override registry defaults, which lets a run vary
    confidence mode or state key without creating a new provider entry.
    """
    # Convenience shorthand: load_ground_truth("some_provider_id")
    if isinstance(method_config, str):
        provider_id = provider_id or method_config
        method_config = {}
    method_config = method_config or {}
    root = Path(project_root) if project_root is not None else _project_root()
    gt_cfg = method_config.get("ground_truth") or {}
    if gt_cfg is None:
        gt_cfg = {}
    if not isinstance(gt_cfg, dict):
        raise ValueError("method_config['ground_truth'] must be a mapping.")

    resolved_provider_id = provider_id or gt_cfg.get("provider_id")
    provider_cfg: Dict[str, Any] = {}
    if resolved_provider_id:
        provider_cfg = _provider_from_registry(str(resolved_provider_id), root)

    if provider_cfg:
        cfg = _merge_provider_config(provider_cfg, gt_cfg)
    elif gt_cfg:
        cfg = dict(gt_cfg)
    else:
        config_cfg = _config_ground_truth_config(method_config)
        config_provider_id = config_cfg.get("provider_id")
        if config_provider_id:
            try:
                provider_cfg = _provider_from_registry(str(config_provider_id), root)
            except (KeyError, FileNotFoundError):
                provider_cfg = {}
        cfg = _merge_provider_config(provider_cfg, config_cfg)

    if "provider_id" not in cfg or not cfg.get("provider_id"):
        cfg["provider_id"] = resolved_provider_id or "config_lineage_provider"

    # Accept both provider vocabulary and config-local lineage vocabulary.
    if not cfg.get("state_key"):
        cfg["state_key"] = cfg.get("cell_state_key")
    if not cfg.get("confidence_mode"):
        cfg["confidence_mode"] = cfg.get("edge_confidence_mode", "all")
    if not cfg.get("reference_graph_path"):
        cfg["reference_graph_path"] = cfg.get("graph_json")
    if not cfg.get("reference_edges_path"):
        cfg["reference_edges_path"] = cfg.get("graph_csv")

    return GroundTruthSpec(
        provider_id=str(cfg.get("provider_id")),
        state_key=cfg.get("state_key"),
        label_path=cfg.get("label_path"),
        metadata_path=cfg.get("metadata_path"),
        reference_graph_path=cfg.get("reference_graph_path"),
        reference_edges_path=cfg.get("reference_edges_path"),
        confidence_mode=str(cfg.get("confidence_mode") or "all"),
        exclude_uncertain_states=bool(cfg.get("exclude_uncertain_states", False)),
        status=cfg.get("status"),
        n_states=cfg.get("n_states"),
        n_graph_edges=cfg.get("n_graph_edges"),
        source_h5ad=cfg.get("source_h5ad"),
        annotation_method=cfg.get("annotation_method"),
        notes=cfg.get("notes") or cfg.get("note"),
        label_mode=cfg.get("label_mode") or None,
        analysis_role=cfg.get("analysis_role") or None,
    )
