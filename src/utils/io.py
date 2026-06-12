"""Small IO helpers: directory creation and table/JSON writers."""
from __future__ import annotations

import json
import os
from typing import Any


def ensure_dir(path: str) -> str:
    os.makedirs(path, exist_ok=True)
    return path


def save_json(obj: Any, path: str) -> str:
    ensure_dir(os.path.dirname(path) or ".")
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, default=str)
    return path


def save_csv(df, path: str, index: bool = True) -> str:
    """Save a pandas DataFrame to CSV, creating parent dirs."""
    ensure_dir(os.path.dirname(path) or ".")
    df.to_csv(path, index=index)
    return path
