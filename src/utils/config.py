"""Config loading and path resolution.

Reads config.yaml, selects the active profile (test|server), and returns a
single flat-ish dict so downstream code never has to know which mode it is in.
"""
from __future__ import annotations

import os
import hashlib
import json
from typing import Any, Dict

import yaml


def load_config(path: str = "config.yaml") -> Dict[str, Any]:
    """Load config.yaml and merge in the active mode profile under ``active``.

    Returns a dict with the original top-level keys plus:
        cfg["mode"]    -> "test" | "server"
        cfg["active"]  -> the resolved profile block for that mode
        cfg["_config_hash"] -> short hash of the raw config (provenance)
    """
    with open(path, "r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)

    mode = raw.get("mode", "test")
    if mode not in raw:
        raise KeyError(f"config mode '{mode}' has no matching profile block")

    raw["active"] = raw[mode]
    raw["_config_hash"] = hashlib.sha1(
        json.dumps(raw, sort_keys=True, default=str).encode()
    ).hexdigest()[:10]
    return raw


def resolve_path(cfg: Dict[str, Any], path: str) -> str:
    """Resolve a possibly-relative path against the project root in cfg."""
    root = cfg.get("project", {}).get("root", ".")
    if os.path.isabs(path):
        return path
    return os.path.normpath(os.path.join(root, path))
