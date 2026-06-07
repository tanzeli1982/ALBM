"""TOML configuration loading helpers for isoALBM runs."""

from __future__ import annotations

import os
import tomllib
from copy import deepcopy
from typing import Any, Dict, Tuple


def load_config(path: str) -> Tuple[Dict[str, Any], str]:
    """Load a TOML config and return (config, base_dir)."""
    config_path = os.path.abspath(path)
    with open(config_path, "rb") as f:
        config = tomllib.load(f)

    config_dir = os.path.dirname(config_path)
    configured_base = config.get("base_dir", ".")
    base_dir = resolve_path(configured_base, config_dir)
    config["_config_path"] = config_path
    config["_base_dir"] = base_dir
    return config, base_dir


def resolve_path(path: str, base_dir: str) -> str:
    """Resolve a config path relative to base_dir unless it is already absolute."""
    if path is None or path == "":
        return path
    expanded = os.path.expanduser(os.path.expandvars(str(path)))
    if os.path.isabs(expanded):
        return os.path.normpath(expanded)
    return os.path.normpath(os.path.join(base_dir, expanded))


def public_config(config: Dict[str, Any]) -> Dict[str, Any]:
    """Return a copy suitable for JSON metadata without private helper keys."""
    cleaned = deepcopy(config)
    cleaned.pop("_config_path", None)
    cleaned.pop("_base_dir", None)
    return cleaned
