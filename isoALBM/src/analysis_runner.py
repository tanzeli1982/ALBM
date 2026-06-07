"""Configurable isoALBM analysis runner."""

from __future__ import annotations

import os
import platform
import time
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution

from config_loader import load_config, public_config, resolve_path
from data_loader import DEFAULT_CACHE_PATH, DEFAULT_EDDY_FLUX_FILE, get_flux_data_dict, get_processed_data
from isotope_model import (
    compute_isotope_timeseries,
    compute_isotope_timeseries_temp,
    get_default_bounds,
    get_default_bounds_temp,
    get_default_obs_targets,
    get_default_params,
    get_default_params_temp,
    get_seasonal_mean,
    run_optimization,
    run_optimization_temp,
)
from observation_loader import (
    LEGACY_TARGET_RULES,
    build_isotope_targets,
    load_flux_observations,
    targets_to_obs_dict,
)
from output_writer import (
    ensure_output_dirs,
    write_flux_comparison,
    write_flux_timeseries,
    write_isotope_timeseries,
    write_optimization_parameters,
    write_resolved_config,
    write_run_metadata,
    write_target_evaluation,
    write_targets,
)


STATIC_PARAM_ORDER = ["alpha_am", "alpha_hm", "alpha_mo", "alpha_op", "f_am", "C13_POM"]
TEMP_PARAM_ORDER = ["m", "b", "alpha_mo", "alpha_op"]
COMPONENT_TIMESERIES_KEYS = {
    "dissolved": "delta_dissolved",
    "emission": "delta_emission",
    "emissions": "delta_emission",
    "emitted": "delta_emission",
    "bubble": "delta_bubble",
    "sediment_production": "C13_sed_prod",
}
COST_PENALTY = 1e12


def run_config_file(
    config_path: str,
    *,
    condition_names: Optional[Sequence[str]] = None,
    output_dir_override: Optional[str] = None,
    maxiter_override: Optional[int] = None,
    workers_override: Optional[int] = None,
    popsize_override: Optional[int] = None,
    no_plots: bool = False,
) -> List[Dict[str, Any]]:
    """Run all configured analyses in a TOML file."""
    config, base_dir = load_config(config_path)
    runs = config.get("runs", [])
    if not runs:
        raise ValueError("Config must contain at least one [[runs]] table")

    results = []
    for run_cfg in runs:
        result = run_single_config(
            run_cfg=run_cfg,
            global_config=config,
            base_dir=base_dir,
            condition_names=condition_names,
            output_dir_override=output_dir_override,
            maxiter_override=maxiter_override,
            workers_override=workers_override,
            popsize_override=popsize_override,
            no_plots=no_plots,
        )
        results.append(result)
    return results


def run_single_config(
    *,
    run_cfg: Dict[str, Any],
    global_config: Dict[str, Any],
    base_dir: str,
    condition_names: Optional[Sequence[str]] = None,
    output_dir_override: Optional[str] = None,
    maxiter_override: Optional[int] = None,
    workers_override: Optional[int] = None,
    popsize_override: Optional[int] = None,
    no_plots: bool = False,
) -> Dict[str, Any]:
    """Run one configured site/ALBM analysis."""
    run_id = run_cfg.get("id", "isoalbm_run")
    output_dir = output_dir_override or run_cfg.get("output_dir", os.path.join("outputs", run_id))
    output_dir = resolve_path(output_dir, base_dir)
    figures_dir = ensure_output_dirs(output_dir)

    print("\n" + "=" * 70)
    print(f"[{time.strftime('%H:%M:%S')}] isoALBM configurable run: {run_id}")
    print("=" * 70)

    data_cfg = run_cfg.get("data", {})
    cache_path = resolve_path(data_cfg.get("cache_path", DEFAULT_CACHE_PATH), base_dir)
    results_dir = resolve_path(data_cfg.get("results_dir", "ALBM data"), base_dir)
    date_range = data_cfg.get("date_range", "20210101_20250101")
    lake_area = float(data_cfg.get("lake_area", 40000.0))
    use_cache = bool(data_cfg.get("use_cache", True))
    save_cache = bool(data_cfg.get("save_cache", False))
    data_eddy_file = resolve_path(data_cfg.get("eddy_flux_file", DEFAULT_EDDY_FLUX_FILE), base_dir)

    albm_data, _cached_eddy = get_processed_data(
        cache_path=cache_path,
        results_dir=results_dir,
        date_range=date_range,
        eddy_flux_file=data_eddy_file,
        use_cache=use_cache,
        save_cache=save_cache,
        lake_area=lake_area,
    )
    flux_data = get_flux_data_dict(albm_data)

    eddy_data = load_flux_observations(run_cfg.get("flux_observations", {}), base_dir)

    default_obs_targets = get_default_obs_targets()
    targets, target_map = build_isotope_targets(run_cfg, base_dir, default_obs_targets)
    obs_targets = targets_to_obs_dict(targets)
    prepared_targets = prepare_targets_for_model(targets, flux_data["time_array"])

    optimizer_cfg = dict(run_cfg.get("optimizer", {}))
    if maxiter_override is not None:
        optimizer_cfg["maxiter"] = maxiter_override
    if workers_override is not None:
        optimizer_cfg["workers"] = workers_override
    if popsize_override is not None:
        optimizer_cfg["popsize"] = popsize_override

    conditions = list(run_cfg.get("conditions", []))
    if condition_names:
        selected = set(condition_names)
        conditions = [condition for condition in conditions if condition.get("name") in selected]
        missing = selected - {condition.get("name") for condition in conditions}
        if missing:
            raise KeyError(f"Configured conditions not found: {sorted(missing)}")
    if not conditions:
        raise ValueError("No conditions selected for this run")

    condition_results = run_conditions(
        conditions=conditions,
        flux_data=flux_data,
        obs_targets=obs_targets,
        prepared_targets=prepared_targets,
        optimizer_cfg=optimizer_cfg,
    )

    evaluation_rows = evaluate_condition_results(condition_results, prepared_targets)

    write_resolved_config(output_dir, global_config.get("_config_path", ""), public_config(global_config))
    write_run_metadata(
        output_dir,
        {
            "run_id": run_id,
            "created_at": datetime.now().isoformat(timespec="seconds"),
            "python": platform.python_version(),
            "platform": platform.platform(),
            "cache_path": cache_path,
            "results_dir": results_dir,
            "date_range": date_range,
            "albm_start_date": str(albm_data.start_date),
            "albm_end_date": str(albm_data.end_date),
            "n_timesteps": int(albm_data.n_timesteps),
            "sediment_temp_shape": list(np.shape(albm_data.sediment_temp_depth)),
            "flux_observations_loaded": eddy_data is not None,
            "conditions": list(condition_results.keys()),
        },
    )
    write_targets(output_dir, targets)
    write_optimization_parameters(output_dir, condition_results)
    write_isotope_timeseries(output_dir, flux_data["time_array"], condition_results)
    write_flux_timeseries(output_dir, albm_data)
    write_flux_comparison(output_dir, albm_data, eddy_data)
    write_target_evaluation(output_dir, evaluation_rows)

    run_plots(
        run_cfg=run_cfg,
        figures_dir=figures_dir,
        albm_data=albm_data,
        eddy_data=eddy_data,
        condition_results=condition_results,
        flux_data=flux_data,
        obs_targets=obs_targets,
        no_plots=no_plots,
    )

    print(f"[{time.strftime('%H:%M:%S')}] Run complete: {output_dir}")
    return {
        "run_id": run_id,
        "output_dir": output_dir,
        "figures_dir": figures_dir,
        "condition_results": condition_results,
        "targets": targets,
    }


def run_conditions(
    *,
    conditions: Iterable[Dict[str, Any]],
    flux_data: Dict[str, Any],
    obs_targets: Dict[str, Dict[str, Any]],
    prepared_targets: List[Dict[str, Any]],
    optimizer_cfg: Dict[str, Any],
) -> Dict[str, Dict[str, Any]]:
    """Run each configured optimization condition."""
    results: Dict[str, Dict[str, Any]] = {}
    condition_list = list(conditions)
    total = len(condition_list)
    for idx, condition in enumerate(condition_list, start=1):
        name = condition["name"]
        model = condition.get("model", "static")
        optimize = bool(condition.get("optimize", True))
        target_names = condition.get("target_names")
        selected_targets = select_targets(prepared_targets, target_names)

        print(f"[{time.strftime('%H:%M:%S')}] Condition {idx}/{total}: {name}")
        if optimize:
            params, opt_meta = run_condition_optimization(
                model=model,
                flux_data=flux_data,
                obs_targets=obs_targets,
                selected_targets=selected_targets,
                optimizer_cfg=optimizer_cfg,
            )
        else:
            params = get_default_params_temp() if model in {"temp", "temperature"} else get_default_params()
            opt_meta = {
                "backend": "initial",
                "cost": np.nan,
                "success": "",
                "n_evaluations": 0,
            }

        if model in {"temp", "temperature"}:
            timeseries = compute_isotope_timeseries_temp(params, flux_data)
        else:
            timeseries = compute_isotope_timeseries(params, flux_data)

        eval_rows = evaluate_timeseries_against_targets(timeseries, selected_targets)
        finite_costs = [row["weighted_error"] for row in eval_rows if np.isfinite(row["weighted_error"])]
        if not optimize:
            opt_meta["cost"] = float(np.sum(finite_costs)) if finite_costs else np.nan

        results[name] = {
            "params": params,
            "timeseries": timeseries,
            "color": condition.get("color", "#555555"),
            "linestyle": condition.get("linestyle", "-"),
            "model": model,
            "optimized": optimize,
            "target_names": [target["name"] for target in selected_targets],
            "optimization": opt_meta,
        }
    return results


def run_condition_optimization(
    *,
    model: str,
    flux_data: Dict[str, Any],
    obs_targets: Dict[str, Dict[str, Any]],
    selected_targets: List[Dict[str, Any]],
    optimizer_cfg: Dict[str, Any],
) -> Tuple[Dict[str, float], Dict[str, Any]]:
    """Run either the legacy paper optimizer or the generic target optimizer."""
    backend = optimizer_cfg.get("backend", "legacy_when_possible")
    use_legacy = backend in {"legacy", "legacy_when_possible"} and _legacy_target_set(selected_targets)
    if backend == "legacy" and not use_legacy:
        raise ValueError("optimizer.backend='legacy' requires legacy-compatible summer/winter targets")

    normalize_by_std = bool(optimizer_cfg.get("normalize_by_std", False))
    maxiter = int(optimizer_cfg.get("maxiter", 200))
    tol = float(optimizer_cfg.get("tol", 0.01))
    popsize = int(optimizer_cfg.get("popsize", 15))
    workers = int(optimizer_cfg.get("workers", -1))
    seed = optimizer_cfg.get("seed", 42)
    verbose = bool(optimizer_cfg.get("verbose", False))
    polish = bool(optimizer_cfg.get("polish", True))

    if use_legacy:
        toggles = {key: False for key in LEGACY_TARGET_RULES}
        for target in selected_targets:
            toggles[target["name"]] = True
        if model in {"temp", "temperature"}:
            result = run_optimization_temp(
                flux_data=flux_data,
                obs_targets=obs_targets,
                target_toggles=toggles,
                normalize_by_std=normalize_by_std,
                maxiter=maxiter,
                tol=tol,
                popsize=popsize,
                workers=workers,
                verbose=verbose,
                seed=seed,
                polish=polish,
            )
            params = {
                "m": float(result.m),
                "b": float(result.b),
                "alpha_mo": float(result.alpha_mo),
                "alpha_op": float(result.alpha_op),
            }
        else:
            result = run_optimization(
                flux_data=flux_data,
                obs_targets=obs_targets,
                target_toggles=toggles,
                normalize_by_std=normalize_by_std,
                maxiter=maxiter,
                tol=tol,
                popsize=popsize,
                workers=workers,
                verbose=verbose,
                seed=seed,
                polish=polish,
            )
            params = {
                "alpha_am": float(result.alpha_am),
                "alpha_hm": float(result.alpha_hm),
                "alpha_mo": float(result.alpha_mo),
                "alpha_op": float(result.alpha_op),
                "f_am": float(result.f_am),
                "C13_POM": float(result.C13_POM),
            }
        return params, {
            "backend": "legacy",
            "cost": float(result.cost),
            "success": bool(result.success),
            "n_evaluations": int(result.n_evaluations),
        }

    params, meta = run_generic_optimization(
        model=model,
        flux_data=flux_data,
        selected_targets=selected_targets,
        optimizer_cfg=optimizer_cfg,
    )
    return params, meta


def run_generic_optimization(
    *,
    model: str,
    flux_data: Dict[str, Any],
    selected_targets: List[Dict[str, Any]],
    optimizer_cfg: Dict[str, Any],
) -> Tuple[Dict[str, float], Dict[str, Any]]:
    """Run generic differential evolution against arbitrary isotope targets."""
    is_temp = model in {"temp", "temperature"}
    param_order = TEMP_PARAM_ORDER if is_temp else STATIC_PARAM_ORDER
    default_bounds = get_default_bounds_temp() if is_temp else get_default_bounds()
    bounds_cfg = optimizer_cfg.get("bounds", {})
    bounds = []
    for name in param_order:
        value = bounds_cfg.get(name, default_bounds[name])
        bounds.append(tuple(value))

    maxiter = int(optimizer_cfg.get("maxiter", 200))
    tol = float(optimizer_cfg.get("tol", 0.01))
    popsize = int(optimizer_cfg.get("popsize", 15))
    workers = int(optimizer_cfg.get("workers", -1))
    seed = optimizer_cfg.get("seed", 42)
    polish = bool(optimizer_cfg.get("polish", True))
    verbose = bool(optimizer_cfg.get("verbose", False))
    normalize_by_std = bool(optimizer_cfg.get("normalize_by_std", False))

    result = differential_evolution(
        _generic_objective,
        bounds=bounds,
        args=(model, param_order, flux_data, selected_targets, normalize_by_std),
        maxiter=maxiter,
        tol=tol,
        popsize=popsize,
        workers=workers,
        updating="deferred" if workers != 1 else "immediate",
        polish=polish,
        disp=verbose,
        seed=seed,
    )
    params = _vector_to_params(result.x, param_order)
    if not is_temp and "f_am" in params:
        params["f_hm"] = 1.0 - params["f_am"]
    return params, {
        "backend": "generic",
        "cost": float(result.fun),
        "success": bool(result.success),
        "n_evaluations": int(result.nfev),
    }


def _generic_objective(
    x: np.ndarray,
    model: str,
    param_order: Sequence[str],
    flux_data: Dict[str, Any],
    selected_targets: List[Dict[str, Any]],
    normalize_by_std: bool,
) -> float:
    params = _vector_to_params(x, param_order)
    if "f_am" in params:
        params["f_hm"] = 1.0 - params["f_am"]
    if model in {"temp", "temperature"}:
        timeseries = compute_isotope_timeseries_temp(params, flux_data)
    else:
        timeseries = compute_isotope_timeseries(params, flux_data)
    rows = evaluate_timeseries_against_targets(timeseries, selected_targets, normalize_by_std)
    costs = [row["weighted_error"] for row in rows if np.isfinite(row["weighted_error"])]
    return float(np.sum(costs)) if costs else COST_PENALTY


def prepare_targets_for_model(
    targets: Iterable[Dict[str, Any]],
    time_array: Any,
) -> List[Dict[str, Any]]:
    """Add model-time masks to normalized targets."""
    time_index = pd.DatetimeIndex(time_array)
    prepared = []
    for target in targets:
        item = dict(target)
        item["mask"] = _target_mask(item, time_index)
        prepared.append(item)
    return prepared


def select_targets(
    prepared_targets: List[Dict[str, Any]],
    target_names: Optional[Sequence[str]],
) -> List[Dict[str, Any]]:
    """Select enabled targets by name, or all enabled targets if no names are supplied."""
    if target_names:
        by_name = {target["name"]: target for target in prepared_targets}
        missing = [name for name in target_names if name not in by_name]
        if missing:
            raise KeyError(f"Unknown isotope target(s): {missing}")
        return [by_name[name] for name in target_names if by_name[name].get("enabled", True)]
    return [target for target in prepared_targets if target.get("enabled", True)]


def evaluate_condition_results(
    condition_results: Dict[str, Dict[str, Any]],
    prepared_targets: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Evaluate every condition against its configured targets."""
    rows: List[Dict[str, Any]] = []
    for condition, result in condition_results.items():
        selected = select_targets(prepared_targets, result.get("target_names"))
        for row in evaluate_timeseries_against_targets(result["timeseries"], selected):
            row["condition"] = condition
            rows.append(row)
    return rows


def evaluate_timeseries_against_targets(
    timeseries: Dict[str, Any],
    selected_targets: List[Dict[str, Any]],
    normalize_by_std: bool = False,
) -> List[Dict[str, Any]]:
    """Return model-vs-observation rows for target means."""
    rows: List[Dict[str, Any]] = []
    for target in selected_targets:
        key = COMPONENT_TIMESERIES_KEYS.get(target["component"], target["component"])
        values = np.asarray(timeseries.get(key, []), dtype=float)
        mask = np.asarray(target.get("mask", np.ones(len(values), dtype=bool)), dtype=bool)
        if len(values) != len(mask):
            model_mean = np.nan
            model_std = np.nan
        else:
            vals = values[mask]
            vals = vals[np.isfinite(vals)]
            model_mean = float(np.mean(vals)) if vals.size else np.nan
            model_std = float(np.std(vals, ddof=1)) if vals.size > 1 else (0.0 if vals.size == 1 else np.nan)

        obs_mean = float(target["mean"])
        obs_std = float(target.get("std", 0.0))
        raw_error = (model_mean - obs_mean) ** 2 if np.isfinite(model_mean) else np.nan
        if normalize_by_std:
            scaled_error = raw_error / (obs_std ** 2) if obs_std > 0 and np.isfinite(raw_error) else np.nan
        else:
            scaled_error = raw_error
        weighted_error = scaled_error * float(target.get("weight", 1.0)) if np.isfinite(scaled_error) else np.nan
        rows.append({
            "target": target["name"],
            "component": target["component"],
            "season": target.get("season") or "",
            "start_date": target.get("start_date") or "",
            "end_date": target.get("end_date") or "",
            "model_mean": model_mean,
            "model_std": model_std,
            "observed_mean": obs_mean,
            "observed_std": obs_std,
            "raw_error": raw_error,
            "weighted_error": weighted_error,
            "weight": float(target.get("weight", 1.0)),
            "n_observations": target.get("n_observations", ""),
            "source": target.get("source", ""),
            "source_file": target.get("source_file", ""),
        })
    return rows


def run_plots(
    *,
    run_cfg: Dict[str, Any],
    figures_dir: str,
    albm_data: Any,
    eddy_data: Any,
    condition_results: Dict[str, Any],
    flux_data: Dict[str, Any],
    obs_targets: Dict[str, Any],
    no_plots: bool = False,
) -> None:
    """Run selected plot functions from flux_plots.py."""
    if no_plots:
        return
    plot_cfg = run_cfg.get("plots", {})
    if not plot_cfg.get("enabled", True):
        return
    show = bool(plot_cfg.get("show", False))

    from flux_plots import (
        plot_ch4_budget,
        plot_dissolved_pool_timeseries,
        plot_eddy_comparison,
        plot_flux_components,
        plot_flux_components_linear,
        plot_optimization_comparison,
        plot_overlap_violins,
        plot_oxic_and_exchange_heatmaps,
        plot_pool_sources_sinks,
        plot_sediment_temperature,
        plot_water_column_doy_climatology_heatmaps,
        plot_water_column_heatmaps,
    )

    plot_calls = {
        "water_column_heatmaps": lambda: plot_water_column_heatmaps(albm_data, figures_dir, show),
        "water_column_doy_climatology_heatmaps": lambda: plot_water_column_doy_climatology_heatmaps(albm_data, figures_dir, show),
        "oxic_and_exchange_heatmaps": lambda: plot_oxic_and_exchange_heatmaps(albm_data, figures_dir, show),
        "eddy_comparison": lambda: plot_eddy_comparison(albm_data, eddy_data, figures_dir, show),
        "flux_components": lambda: plot_flux_components(albm_data, figures_dir, show),
        "flux_components_linear": lambda: plot_flux_components_linear(albm_data, figures_dir, show),
        "ch4_budget": lambda: plot_ch4_budget(albm_data, figures_dir, show),
        "pool_sources_sinks": lambda: plot_pool_sources_sinks(albm_data, figures_dir, show),
        "sediment_temperature": lambda: plot_sediment_temperature(albm_data, figures_dir, show),
        "dissolved_pool_timeseries": lambda: plot_dissolved_pool_timeseries(condition_results, flux_data, figures_dir, show=show),
    }

    for name, call in plot_calls.items():
        if bool(plot_cfg.get(name, False)):
            if name == "eddy_comparison" and eddy_data is None:
                print("Skipping eddy_comparison plot because no methane flux observations were loaded.")
                continue
            call()

    if _has_paper_obs_targets(obs_targets):
        if bool(plot_cfg.get("optimization_comparison", False)):
            plot_optimization_comparison(condition_results, flux_data, obs_targets, figures_dir, show=show)
        if bool(plot_cfg.get("overlap_violins", False)):
            plot_overlap_violins(condition_results, flux_data, obs_targets, figures_dir, show=show)
    elif plot_cfg.get("optimization_comparison", False) or plot_cfg.get("overlap_violins", False):
        print("Skipping paper-style isotope comparison plots because summer/winter paper targets are not all present.")


def _target_mask(target: Dict[str, Any], time_index: pd.DatetimeIndex) -> np.ndarray:
    if target.get("start_date") or target.get("end_date"):
        mask = np.ones(len(time_index), dtype=bool)
        if target.get("start_date"):
            mask &= time_index >= pd.Timestamp(target["start_date"])
        if target.get("end_date"):
            mask &= time_index <= pd.Timestamp(target["end_date"])
        return mask
    season = target.get("season")
    if season == "summer":
        return np.asarray(
            ((time_index.month == 6) & (time_index.day >= 21))
            | (time_index.month == 7)
            | (time_index.month == 8)
            | ((time_index.month == 9) & (time_index.day <= 22)),
            dtype=bool,
        )
    if season == "winter":
        return np.asarray(
            ((time_index.month == 12) & (time_index.day >= 21))
            | (time_index.month == 1)
            | (time_index.month == 2)
            | ((time_index.month == 3) & (time_index.day <= 20)),
            dtype=bool,
        )
    return np.ones(len(time_index), dtype=bool)


def _legacy_target_set(selected_targets: List[Dict[str, Any]]) -> bool:
    if not selected_targets:
        return False
    return all(target.get("legacy_compatible", False) for target in selected_targets)


def _has_paper_obs_targets(obs_targets: Dict[str, Any]) -> bool:
    return all(key in obs_targets for key in LEGACY_TARGET_RULES)


def _vector_to_params(x: Sequence[float], param_order: Sequence[str]) -> Dict[str, float]:
    return {name: float(value) for name, value in zip(param_order, x)}
