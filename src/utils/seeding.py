"""Deterministic seeding across numpy / random / (optional) torch / scanpy."""
from __future__ import annotations

import os
import random


def set_global_seed(seed: int = 0) -> None:
    """Set seeds everywhere we can, including libraries that may be absent."""
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    try:
        import numpy as np

        np.random.seed(seed)
    except Exception:
        pass
    try:  # optional
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except Exception:
        pass
    try:  # optional
        import scanpy as sc

        sc.settings.seed = seed
    except Exception:
        pass
