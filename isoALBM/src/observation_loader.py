"""Observation loading and target normalization for configurable isoALBM runs."""

from __future__ import annotations

import os
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

from config_loader import resolve_path
from data_loader import EddyFluxData, load_eddy_flux_data


STANDARD_ISOTOPE_COLUMNS = {
    "date": "date",
    "value": "delta13c_ch4",
    "component": "component",
    "season": "season",
    "depth": "depth_m",
    "site": "site",
}


COMPONENT_ALIASES = {
    "dissolved": "dissolved",
    "dissolved_ch4": "dissolved",
    "water": "dissolved",
    "emission": "emission",
    "emissions": "emission",
    "emitted": "emission",
    "surface_emission": "emission",
    "surface_emissions": "emission",
    "bubble": "bubble",
    "ebullition": "bubble",
}


LEGACY_TARGET_RULES = {
    "dissolved_summer": ("dissolved", "summer"),
    "dissolved_winter": ("dissolved", "winter"),
    "emissions_summer": ("emission", "summer"),
    "emissions_winter": ("emission", "winter"),
}


def normalize_component(value: Any) -> str:
    """Normalize component names used in configs and raw tables."""
    text = str(value).strip().lower().replace(" ", "_").replace("-", "_")
    return COMPONENT_ALIASES.get(text, text)


def load_isotope_observation_tables(
    table_configs: Iterable[Dict[str, Any]],
    base_dir: str,
) -> Dict[str, pd.DataFrame]:
    """Load configured isotope observation CSV/XLSX tables by id."""
    tables: Dict[str, pd.DataFrame] = {}
    for idx, table_cfg in enumerate(table_configs or []):
        table_id = table_cfg.get("id", f"table_{idx + 1}")
        path = resolve_path(table_cfg["path"], base_dir)
        sheet = table_cfg.get("sheet")
        ext = os.path.splitext(path)[1].lower()
        if ext in {".xlsx", ".xls"}:
            df = pd.read_excel(path, sheet_name=sheet or 0)
        else:
            df = pd.read_csv(path)
        tables[table_id] = df
    return tables


def build_isotope_targets(
    run_config: Dict[str, Any],
    base_dir: str,
    default_targets: Optional[Dict[str, Dict[str, float]]] = None,
) -> Tuple[List[Dict[str, Any]], Dict[str, Dict[str, Any]]]:
    """Build normalized isotope targets from summary config and raw tables."""
    table_configs = run_config.get("isotope_observation_tables", [])
    tables = load_isotope_observation_tables(table_configs, base_dir)
    table_config_by_id = {
        cfg.get("id", f"table_{idx + 1}"): cfg
        for idx, cfg in enumerate(table_configs or [])
    }

    targets: List[Dict[str, Any]] = []
    configured_targets = run_config.get("isotope_targets", [])

    if not configured_targets and default_targets:
        for name, values in default_targets.items():
            component, season = LEGACY_TARGET_RULES.get(name, ("dissolved", "all"))
            targets.append(_summary_target(name, component, season, values, {}))
    else:
        for target_cfg in configured_targets:
            source = target_cfg.get("source", "summary")
            if source == "table":
                target = _target_from_table(target_cfg, tables, table_config_by_id)
            else:
                defaults = default_targets.get(target_cfg["name"], {}) if default_targets else {}
                merged = {**defaults, **target_cfg}
                target = _summary_target(
                    name=target_cfg["name"],
                    component=target_cfg.get("component", _default_component(target_cfg["name"])),
                    season=target_cfg.get("season"),
                    values=merged,
                    target_cfg=target_cfg,
                )
            targets.append(target)

    target_map = {target["name"]: target for target in targets}
    return targets, target_map


def _summary_target(
    name: str,
    component: Any,
    season: Optional[str],
    values: Dict[str, Any],
    target_cfg: Dict[str, Any],
) -> Dict[str, Any]:
    target = {
        "name": name,
        "component": normalize_component(component),
        "season": _normalize_optional(values.get("season", season)),
        "start_date": _normalize_optional(values.get("start_date")),
        "end_date": _normalize_optional(values.get("end_date")),
        "mean": float(values["mean"]),
        "std": float(values.get("std", 0.0)),
        "weight": float(values.get("weight", 1.0)),
        "enabled": bool(values.get("enabled", True)),
        "source": values.get("source", "summary"),
        "source_file": values.get("source_file", ""),
        "n_observations": values.get("n_observations", ""),
    }
    for key in (
        "dissolved_mean",
        "dissolved_std",
        "ebullition_mean",
        "ebullition_std",
        "diffusion_fraction",
        "notes",
    ):
        if key in values:
            target[key] = values[key]
    target["legacy_compatible"] = _is_legacy_compatible(target)
    return target


def _target_from_table(
    target_cfg: Dict[str, Any],
    tables: Dict[str, pd.DataFrame],
    table_config_by_id: Dict[str, Dict[str, Any]],
) -> Dict[str, Any]:
    table_id = target_cfg.get("table")
    if table_id not in tables:
        raise KeyError(f"Isotope target {target_cfg.get('name')} references unknown table {table_id!r}")
    table_cfg = table_config_by_id.get(table_id, {})
    df = tables[table_id].copy()
    columns = _table_columns(table_cfg)
    value_col = target_cfg.get("value_column", columns["value"])
    date_col = target_cfg.get("date_column", columns["date"])
    component_col = target_cfg.get("component_column", columns["component"])
    season_col = target_cfg.get("season_column", columns["season"])
    site_col = target_cfg.get("site_column", columns["site"])

    if value_col not in df.columns:
        raise KeyError(f"Value column {value_col!r} not found in isotope table {table_id!r}")

    if date_col in df.columns:
        df["_isoalbm_date"] = pd.to_datetime(df[date_col], errors="coerce")
    else:
        df["_isoalbm_date"] = pd.NaT

    if component_col in df.columns:
        df["_isoalbm_component"] = df[component_col].map(normalize_component)
    else:
        df["_isoalbm_component"] = ""

    if season_col in df.columns:
        df["_isoalbm_season"] = df[season_col].astype(str).str.lower().str.strip()
    else:
        df["_isoalbm_season"] = df["_isoalbm_date"].map(_season_for_timestamp)

    component = normalize_component(target_cfg.get("component", target_cfg.get("filter_component", "")))
    if component and component_col in df.columns:
        df = df[df["_isoalbm_component"] == component]

    season = _normalize_optional(target_cfg.get("season", target_cfg.get("filter_season")))
    if season:
        df = df[df["_isoalbm_season"] == season]

    start_date = _normalize_optional(target_cfg.get("start_date"))
    end_date = _normalize_optional(target_cfg.get("end_date"))
    if start_date:
        df = df[df["_isoalbm_date"] >= pd.Timestamp(start_date)]
    if end_date:
        df = df[df["_isoalbm_date"] <= pd.Timestamp(end_date)]

    filter_site = target_cfg.get("filter_site")
    if filter_site and site_col in df.columns:
        df = df[df[site_col].astype(str) == str(filter_site)]

    values = pd.to_numeric(df[value_col], errors="coerce").dropna()
    if values.empty:
        raise ValueError(f"No isotope observations remain for target {target_cfg.get('name')!r}")

    mean = float(values.mean())
    std = float(values.std(ddof=1)) if len(values) > 1 else 0.0
    target = {
        "name": target_cfg["name"],
        "component": component or _default_component(target_cfg["name"]),
        "season": season,
        "start_date": start_date,
        "end_date": end_date,
        "mean": mean,
        "std": float(target_cfg.get("std", std)),
        "weight": float(target_cfg.get("weight", 1.0)),
        "enabled": bool(target_cfg.get("enabled", True)),
        "source": "table",
        "source_file": table_cfg.get("path", ""),
        "table": table_id,
        "n_observations": int(len(values)),
    }
    target["legacy_compatible"] = _is_legacy_compatible(target)
    return target


def load_flux_observations(
    flux_config: Dict[str, Any],
    base_dir: str,
) -> Optional[EddyFluxData]:
    """Load optional methane flux observations for comparison plots/statistics."""
    if not flux_config or not bool(flux_config.get("enabled", False)):
        return None

    path = resolve_path(flux_config.get("path", ""), base_dir)
    required = bool(flux_config.get("required", False))
    if not path or not os.path.exists(path):
        message = f"Methane flux observation file not found: {path}"
        if required:
            raise FileNotFoundError(message)
        print(f"Warning: {message}; skipping flux comparison.")
        return None

    fmt = flux_config.get("format", "auto")
    if fmt == "auto":
        return load_eddy_flux_data(path)

    df = _read_table(path, flux_config.get("sheet"))
    date_col = flux_config.get("date_column", "date")
    flux_col = flux_config.get("flux_column", "ch4_flux")
    units = flux_config.get("units", "kg_m2_s")
    if date_col not in df.columns or flux_col not in df.columns:
        raise KeyError(f"Flux table must contain {date_col!r} and {flux_col!r}")

    daily = pd.DataFrame({
        "datetime": pd.to_datetime(df[date_col], errors="coerce"),
        "flux": pd.to_numeric(df[flux_col], errors="coerce"),
    }).dropna()
    daily = daily.groupby(daily["datetime"].dt.date)["flux"].mean().reset_index()
    daily["datetime"] = pd.to_datetime(daily["datetime"])
    daily["ch4_flux_daily_avg_kg_m2_s"] = _convert_flux_units(daily["flux"], units)
    daily = daily[["datetime", "ch4_flux_daily_avg_kg_m2_s"]]
    return EddyFluxData(
        datetime=pd.DatetimeIndex(daily["datetime"]),
        ch4_flux=daily["ch4_flux_daily_avg_kg_m2_s"].values,
        ch4_flux_daily_avg=daily,
    )


def targets_to_obs_dict(targets: Iterable[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """Convert normalized targets to the obs_targets dict expected by plot helpers."""
    out: Dict[str, Dict[str, Any]] = {}
    for target in targets:
        public = {k: v for k, v in target.items() if k != "mask"}
        out[target["name"]] = public
    return out


def _read_table(path: str, sheet: Optional[str] = None) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    if ext in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=sheet or 0)
    return pd.read_csv(path)


def _table_columns(table_cfg: Dict[str, Any]) -> Dict[str, str]:
    columns = dict(STANDARD_ISOTOPE_COLUMNS)
    mapping = table_cfg.get("columns", {})
    columns.update({k: v for k, v in mapping.items() if v})
    for key in STANDARD_ISOTOPE_COLUMNS:
        config_key = f"{key}_column"
        if config_key in table_cfg:
            columns[key] = table_cfg[config_key]
    return columns


def _convert_flux_units(values: pd.Series, units: str) -> pd.Series:
    units_norm = units.lower().replace(" ", "").replace("-", "_")
    if units_norm in {"kg_m2_s", "kg/m2/s", "kgm-2s-1"}:
        return values.astype(float)
    if units_norm in {"mg_m2_day", "mg/m2/day", "mgm-2day-1"}:
        return values.astype(float) * 1e-6 / 86400.0
    if units_norm in {"umol_m2_s", "umol/m2/s", "umolm-2s-1"}:
        return values.astype(float) * 1e-6 * 0.01604
    raise ValueError(f"Unsupported methane flux units: {units}")


def _default_component(target_name: str) -> str:
    return LEGACY_TARGET_RULES.get(target_name, ("dissolved", ""))[0]


def _normalize_optional(value: Any) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    return text.lower() if text else None


def _season_for_timestamp(value: Any) -> str:
    if pd.isna(value):
        return ""
    ts = pd.Timestamp(value)
    if (ts.month == 6 and ts.day >= 21) or ts.month in {7, 8} or (ts.month == 9 and ts.day <= 22):
        return "summer"
    if (ts.month == 12 and ts.day >= 21) or ts.month in {1, 2} or (ts.month == 3 and ts.day <= 20):
        return "winter"
    return "shoulder"


def _is_legacy_compatible(target: Dict[str, Any]) -> bool:
    rule = LEGACY_TARGET_RULES.get(target["name"])
    if not rule:
        return False
    component, season = rule
    return (
        target.get("component") == component
        and target.get("season") == season
        and not target.get("start_date")
        and not target.get("end_date")
        and float(target.get("weight", 1.0)) == 1.0
    )
