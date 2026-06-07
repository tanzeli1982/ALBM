"""CSV and metadata writers for configurable isoALBM analyses."""

from __future__ import annotations

import json
import os
import shutil
from datetime import datetime, timezone
from typing import Any, Dict, Iterable, List, Optional

import numpy as np
import pandas as pd


def ensure_output_dirs(output_dir: str) -> str:
    """Create output_dir and its figures subdirectory. Return figures path."""
    os.makedirs(output_dir, exist_ok=True)
    figures_dir = os.path.join(output_dir, "figures")
    os.makedirs(figures_dir, exist_ok=True)
    return figures_dir


def write_run_metadata(output_dir: str, metadata: Dict[str, Any]) -> str:
    """Write run metadata JSON."""
    path = os.path.join(output_dir, "run_metadata.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(_jsonable(metadata), f, indent=2, sort_keys=True)
    return path


def write_resolved_config(output_dir: str, config_path: str, config: Dict[str, Any]) -> None:
    """Copy the source TOML and write a JSON representation of the resolved config."""
    if config_path and os.path.exists(config_path):
        shutil.copy2(config_path, os.path.join(output_dir, "resolved_config.toml"))
    path = os.path.join(output_dir, "resolved_config.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(_jsonable(config), f, indent=2, sort_keys=True)


def write_targets(output_dir: str, targets: Iterable[Dict[str, Any]]) -> str:
    """Write normalized isotope targets."""
    rows = []
    for target in targets:
        rows.append({k: v for k, v in target.items() if k != "mask"})
    path = os.path.join(output_dir, "isotope_observation_targets.csv")
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def write_optimization_parameters(
    output_dir: str,
    condition_results: Dict[str, Dict[str, Any]],
) -> str:
    """Write one row per condition with parameters and optimizer metadata."""
    rows: List[Dict[str, Any]] = []
    for name, result in condition_results.items():
        row = {
            "condition": name,
            "model": result.get("model", ""),
            "optimized": result.get("optimized", ""),
        }
        row.update(result.get("optimization", {}))
        row.update(result.get("params", {}))
        rows.append(row)
    path = os.path.join(output_dir, "optimization_parameters.csv")
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def write_isotope_timeseries(
    output_dir: str,
    time_array: Any,
    condition_results: Dict[str, Dict[str, Any]],
) -> str:
    """Write modeled isotope time series for each condition."""
    df = pd.DataFrame({"time": pd.to_datetime(time_array)})
    preferred_keys = [
        "delta_dissolved",
        "delta_emission",
        "delta_bubble",
        "C13_sed_prod",
        "ch4_conc_mb",
    ]
    for condition, result in condition_results.items():
        safe = _safe_name(condition)
        ts = result.get("timeseries", {})
        for key in preferred_keys:
            if key in ts:
                values = np.asarray(ts[key])
                if values.ndim == 1 and len(values) == len(df):
                    df[f"{safe}.{key}"] = values
    path = os.path.join(output_dir, "isotope_timeseries.csv")
    df.to_csv(path, index=False)
    return path


def write_flux_timeseries(output_dir: str, albm_data: Any) -> str:
    """Write ALBM methane flux and related time series."""
    sediment_temp_avg = getattr(albm_data, "sediment_temp_avg", None)
    if sediment_temp_avg is None:
        sediment_temp_avg = np.full(len(albm_data.time_array), np.nan)
    df = pd.DataFrame({
        "time": pd.to_datetime(albm_data.time_array),
        "total_flux_kg_m2_s": albm_data.total_flux,
        "diffusion_flux_kg_m2_s": albm_data.dflux_data,
        "ebullition_flux_kg_m2_s": albm_data.eflux_data,
        "upward_bubbling_flux_kg_m2_s": albm_data.ch4upb_data,
        "sediment_temperature_avg_c": sediment_temp_avg,
    })
    path = os.path.join(output_dir, "flux_timeseries.csv")
    df.to_csv(path, index=False)
    return path


def write_flux_comparison(output_dir: str, albm_data: Any, eddy_data: Optional[Any]) -> Optional[str]:
    """Write daily ALBM vs methane flux observation comparison if observations exist."""
    if eddy_data is None:
        return None
    model = pd.DataFrame({
        "date": pd.to_datetime(albm_data.time_array).date,
        "albm_total_flux_kg_m2_s": albm_data.total_flux,
    })
    obs = eddy_data.ch4_flux_daily_avg.copy()
    obs["date"] = pd.to_datetime(obs["datetime"]).dt.date
    obs = obs.rename(columns={"ch4_flux_daily_avg_kg_m2_s": "observed_flux_kg_m2_s"})
    out = pd.merge(model, obs[["date", "observed_flux_kg_m2_s"]], on="date", how="inner")
    out["albm_minus_observed_kg_m2_s"] = out["albm_total_flux_kg_m2_s"] - out["observed_flux_kg_m2_s"]
    path = os.path.join(output_dir, "methane_flux_comparison.csv")
    out.to_csv(path, index=False)
    return path


def write_target_evaluation(output_dir: str, rows: Iterable[Dict[str, Any]]) -> str:
    """Write per-condition model-vs-target evaluations."""
    path = os.path.join(output_dir, "observation_overlap_statistics.csv")
    pd.DataFrame(list(rows)).to_csv(path, index=False)
    return path


def _safe_name(value: str) -> str:
    return "".join(c.lower() if c.isalnum() else "_" for c in value).strip("_")


def _jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, (datetime, pd.Timestamp)):
        return value.isoformat()
    if value is timezone.utc:
        return "UTC"
    return value
