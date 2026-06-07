"""
Multi-Condition Isotope Optimization Script

Runs optimizations under multiple observation conditions (e.g. all obs, no winter obs)
and plots the comparison. Uses data_loader, isotope_model, and flux_plots.

Usage:
    python src/manuscript_scripts/multi_condition_optimization.py

    # Optional: show matplotlib windows while also saving figures
    $env:ISOALBM_SHOW_PLOTS = "1"
    python src/manuscript_scripts/multi_condition_optimization.py
"""

import os
import sys
import time
import numpy as np
import pandas as pd

def _find_src_dir():
    """Find src/ when run as a script or from an interactive cell."""
    script_file = globals().get("__file__", "")
    script_dir = os.path.dirname(os.path.abspath(script_file)) if script_file else ""
    cwd = os.path.abspath(os.getcwd())
    candidates = [
        os.path.join(script_dir, os.pardir),
        os.path.join(script_dir, "src"),
        script_dir,
        os.path.join(cwd, "src"),
        cwd,
        os.path.join(cwd, os.pardir),
        os.path.join(cwd, os.pardir, "src"),
    ]
    for candidate in candidates:
        candidate = os.path.normpath(os.path.abspath(candidate))
        if os.path.isfile(os.path.join(candidate, "data_loader.py")):
            return candidate
    raise ModuleNotFoundError(
        "Could not locate src/data_loader.py. Run from the repository root, "
        "or add the repository's src directory to sys.path before importing."
    )


SRC_DIR = _find_src_dir()
PROJECT_ROOT = os.path.abspath(os.path.join(SRC_DIR, os.pardir))
if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from data_loader import (
    get_processed_data,
    get_flux_data_dict,
    DEFAULT_CACHE_PATH,
    DEFAULT_EDDY_FLUX_FILE,
)
from isotope_model import (
    run_optimization,
    run_optimization_temp,
    compute_isotope_timeseries,
    compute_isotope_timeseries_temp,
    get_default_obs_targets,
    get_default_params,
    get_default_params_temp,
)
from flux_plots import (
    plot_optimization_comparison,
    plot_dissolved_pool_timeseries,
    CB_ORANGE,
    CB_GREEN,
    CB_BLUE,
    CB_VERMILLION,
)
try:
    from flux_plots import plot_overlap_violins
except ImportError:
    plot_overlap_violins = None  # optional: upgrade flux_plots or restart kernel
try:
    from flux_plots import compute_and_print_paper_statistics
except ImportError:
    compute_and_print_paper_statistics = None  # optional: run paper stats when available

# =============================================================================
# Configuration
# =============================================================================

RESULTS_DIR = 'ALBM data'
DATE_RANGE = '20210101_20250101'
EDDY_FLUX_FILE = DEFAULT_EDDY_FLUX_FILE
FIGS_DIR = os.path.join(PROJECT_ROOT, 'outputs')
CACHE_PATH = DEFAULT_CACHE_PATH
USE_CACHE = True  # Use cached data if present (for switching machines)

# DE stops when std(population cost) <= atol + tol*|mean cost|
def _env_int(name, default):
    value = os.environ.get(name)
    if value is None or value.strip() == "":
        return default
    try:
        return int(value)
    except ValueError:
        print(f"Warning: ignoring invalid {name}={value!r}; using {default}.")
        return default


def _env_bool(name, default=False):
    value = os.environ.get(name)
    if value is None or value.strip() == "":
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


OPT_MAXITER = 200
OPT_TOL = 0.01
OPT_POPSIZE = 15
OPT_WORKERS = _env_int('ISOALBM_WORKERS', os.cpu_count() or 1)
OPT_VERBOSE = True
OPT_SEED = 42       # Fixed seed for reproducible optimization (set None for random)
NORMALIZE_BY_STD = False
SHOW_PLOTS = _env_bool('ISOALBM_SHOW_PLOTS', False)

LAST_ALBM_DATA = None
LAST_EDDY_DATA = None
LAST_FLUX_DATA = None
LAST_OBS_TARGETS = None
LAST_CONDITION_RESULTS = None

# Conditions to run (name -> config with target_toggles, color, linestyle, marker)
CONDITIONS = {
    'Initial guess': {
        'run_optimization': False,
        'target_toggles': None,
        'color': CB_GREEN,
        'linestyle': '-',
        'marker': 'o'
    },
    'All Observations': {
        'run_optimization': True,
        'target_toggles': {
            'dissolved_summer': True,
            'dissolved_winter': True,
            'emissions_summer': True,
            'emissions_winter': True
        },
        'color': CB_ORANGE,
        'linestyle': '-',
        'marker': 's'
    },
    'No Winter Obs.': {
        'run_optimization': True,
        'target_toggles': {
            'dissolved_summer': True,
            'dissolved_winter': False,
            'emissions_summer': True,
            'emissions_winter': False
        },
        'color': CB_GREEN,
        'linestyle': '-',
        'marker': '^'
    },
    'No Winter Emission Obs.': {
        'run_optimization': True,
        'target_toggles': {
            'dissolved_summer': True,
            'dissolved_winter': True,
            'emissions_summer': True,
            'emissions_winter': False
        },
        'color': CB_BLUE,
        'linestyle': '-',
        'marker': 'D'
    },
    'Temperature-Based': {
        'run_optimization': True,
        'use_temp_model': True,  # Flag to use temperature-based model
        'target_toggles': {
            'dissolved_summer': True,
            'dissolved_winter': True,
            'emissions_summer': True,
            'emissions_winter': True
        },
        'color': CB_VERMILLION,
        'linestyle': '-.',
        'marker': 'P'
    }
}


def _observation_period_masks(time_array):
    """
    Return list of (period_name, mask, season) for each observation period, plus
    (all_summer_mask, all_winter_mask, all_overlap_mask).
    season is 'summer' or 'winter'. Summer: Jun 21–Sep 22; Winter: Dec 21–Mar 20 (next year).
    """
    time_pd = pd.DatetimeIndex(time_array)
    n = len(time_pd)
    all_summer_mask = np.zeros(n, dtype=bool)
    all_winter_mask = np.zeros(n, dtype=bool)
    periods = []
    min_year = int(time_pd.year.min())
    max_year = int(time_pd.year.max())
    for year in range(min_year, max_year + 1):
        summer_start = pd.Timestamp(f'{year}-06-21')
        summer_end = pd.Timestamp(f'{year}-09-22')
        summer_mask = (time_pd >= summer_start) & (time_pd <= summer_end)
        if np.any(summer_mask):
            periods.append((f'Summer {year}', summer_mask, 'summer'))
            all_summer_mask |= summer_mask
        winter_dec_start = pd.Timestamp(f'{year}-12-21')
        winter_dec_end = pd.Timestamp(f'{year}-12-31')
        winter_jan_start = pd.Timestamp(f'{year + 1}-01-01')
        winter_jan_end = pd.Timestamp(f'{year + 1}-03-20')
        winter_mask = (
            ((time_pd >= winter_dec_start) & (time_pd <= winter_dec_end)) |
            ((time_pd >= winter_jan_start) & (time_pd <= winter_jan_end))
        )
        if year < max_year and np.any(winter_mask):
            periods.append((f'Winter {year}-{year + 1}', winter_mask, 'winter'))
            all_winter_mask |= winter_mask
    all_overlap = all_summer_mask | all_winter_mask
    return periods, all_summer_mask, all_winter_mask, all_overlap


def _mean_std(arr):
    """Return (mean, std) for an array; std uses ddof=1. Handles empty/single-point."""
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan, np.nan
    m = float(np.mean(arr))
    s = float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0
    return m, s


def print_overlap_statistics(condition_results, flux_data, obs_targets):
    """
    Report model mean ± std per observation period (Summer 2021, Summer 2022, …),
    then All summers, All winters, All observation periods; compare to observations.
    """
    time_array = flux_data['time_array']
    periods, all_summer_mask, all_winter_mask, all_mask = _observation_period_masks(time_array)

    obs_diss_s = obs_targets['dissolved_summer']
    obs_diss_w = obs_targets['dissolved_winter']
    obs_emis_s = obs_targets['emissions_summer']
    obs_emis_w = obs_targets['emissions_winter']

    print("\n" + "="*100)
    print("MODEL vs OBSERVATIONS: δ¹³C in observation overlap periods (mean ± std, diff from obs mean and std)")
    print("="*100)

    for condition_name, result in condition_results.items():
        ts = result['timeseries']
        d = ts['delta_dissolved']
        e = ts['delta_emission']

        print(f"\n  {condition_name}")
        print("  " + "-"*90)
        print("  Dissolved CH₄ δ¹³C (‰)")

        for period_name, mask, season in periods:
            vals = d[mask]
            m, s = _mean_std(vals)
            obs = obs_diss_s if season == 'summer' else obs_diss_w
            diff_mean = m - obs['mean'] if not np.isnan(m) else np.nan
            diff_std = s - obs['std'] if not np.isnan(s) else np.nan
            diff_str = ""
            if not np.isnan(diff_mean) or not np.isnan(diff_std):
                parts = []
                if not np.isnan(diff_mean):
                    parts.append(f"diff mean: {diff_mean:+.2f} ‰")
                if not np.isnan(diff_std):
                    parts.append(f"diff std: {diff_std:+.2f} ‰")
                diff_str = "   |   " + "  ".join(parts)
            print(f"    {period_name}:  Model {m:.2f} ± {s:.2f}{diff_str}  (obs: {obs['mean']:.2f} ± {obs['std']:.2f})")

        m_s, s_s = _mean_std(d[all_summer_mask])
        m_w, s_w = _mean_std(d[all_winter_mask])
        m_all, s_all = _mean_std(d[all_mask])
        print(f"    All summers:   Model {m_s:.2f} ± {s_s:.2f}   |   diff mean: {m_s - obs_diss_s['mean']:+.2f} ‰  diff std: {s_s - obs_diss_s['std']:+.2f} ‰  (obs: {obs_diss_s['mean']:.2f} ± {obs_diss_s['std']:.2f})")
        print(f"    All winters:   Model {m_w:.2f} ± {s_w:.2f}   |   diff mean: {m_w - obs_diss_w['mean']:+.2f} ‰  diff std: {s_w - obs_diss_w['std']:+.2f} ‰  (obs: {obs_diss_w['mean']:.2f} ± {obs_diss_w['std']:.2f})")
        print(f"    All periods:   Model {m_all:.2f} ± {s_all:.2f}")

        print("  Surface emitted CH₄ δ¹³C (‰)")
        for period_name, mask, season in periods:
            vals = e[mask]
            m, s = _mean_std(vals)
            obs = obs_emis_s if season == 'summer' else obs_emis_w
            diff_mean = m - obs['mean'] if not np.isnan(m) else np.nan
            diff_std = s - obs['std'] if not np.isnan(s) else np.nan
            diff_str = ""
            if not np.isnan(diff_mean) or not np.isnan(diff_std):
                parts = []
                if not np.isnan(diff_mean):
                    parts.append(f"diff mean: {diff_mean:+.2f} ‰")
                if not np.isnan(diff_std):
                    parts.append(f"diff std: {diff_std:+.2f} ‰")
                diff_str = "   |   " + "  ".join(parts)
            print(f"    {period_name}:  Model {m:.2f} ± {s:.2f}{diff_str}  (obs: {obs['mean']:.2f} ± {obs['std']:.2f})")

        m_s, s_s = _mean_std(e[all_summer_mask])
        m_w, s_w = _mean_std(e[all_winter_mask])
        m_all, s_all = _mean_std(e[all_mask])
        print(f"    All summers:   Model {m_s:.2f} ± {s_s:.2f}   |   diff mean: {m_s - obs_emis_s['mean']:+.2f} ‰  diff std: {s_s - obs_emis_s['std']:+.2f} ‰  (obs: {obs_emis_s['mean']:.2f} ± {obs_emis_s['std']:.2f})")
        print(f"    All winters:   Model {m_w:.2f} ± {s_w:.2f}   |   diff mean: {m_w - obs_emis_w['mean']:+.2f} ‰  diff std: {s_w - obs_emis_w['std']:+.2f} ‰  (obs: {obs_emis_w['mean']:.2f} ± {obs_emis_w['std']:.2f})")
        print(f"    All periods:   Model {m_all:.2f} ± {s_all:.2f}")

    print("="*100 + "\n")


def rerun_paper_statistics(
    condition_results=None,
    flux_data=None,
    obs_targets=None,
    albm_data=None,
    eddy_data=None,
    figs_dir=None,
):
    """
    Print/write paper statistics for the most recent run.

    In an interactive console, run this after main() if you want to print the
    paper-stat section again without rerunning the optimizations.
    """
    condition_results = condition_results if condition_results is not None else LAST_CONDITION_RESULTS
    flux_data = flux_data if flux_data is not None else LAST_FLUX_DATA
    obs_targets = obs_targets if obs_targets is not None else LAST_OBS_TARGETS
    albm_data = albm_data if albm_data is not None else LAST_ALBM_DATA
    eddy_data = eddy_data if eddy_data is not None else LAST_EDDY_DATA
    figs_dir = figs_dir or FIGS_DIR

    missing = []
    if condition_results is None:
        missing.append("condition_results")
    if flux_data is None:
        missing.append("flux_data")
    if obs_targets is None:
        missing.append("obs_targets")
    if missing:
        raise RuntimeError(
            "Cannot print paper statistics yet; missing "
            + ", ".join(missing)
            + ". Run main() first, or pass those objects explicitly."
        )

    if albm_data is None:
        print(f"[{time.strftime('%H:%M:%S')}] Skipping paper statistics: albm_data was not provided.", flush=True)
        return None
    if compute_and_print_paper_statistics is None:
        print(f"[{time.strftime('%H:%M:%S')}] Skipping paper statistics: compute_and_print_paper_statistics "
              "could not be imported from flux_plots.", flush=True)
        return None

    os.makedirs(figs_dir, exist_ok=True)
    paper_stats_path = os.path.join(figs_dir, 'paper_statistics_flux_isotope.csv')
    print(f"[{time.strftime('%H:%M:%S')}] Printing paper statistics; CSV -> {paper_stats_path}", flush=True)
    return compute_and_print_paper_statistics(
        albm_data=albm_data,
        eddy_data=eddy_data,
        condition_results=condition_results,
        flux_data=flux_data,
        obs_targets=obs_targets,
        csv_path=paper_stats_path,
    )


def run_multi_condition_comparison(flux_data, obs_targets, figs_dir='figs', show=True, albm_data=None, eddy_data=None):
    """
    Run optimization under each condition and create comparison plot.

    If albm_data and eddy_data are provided, paper statistics (B–J) are printed
    after the comparison. Section A (ERA5) is produced by climateChangeAnalysis.py.

    Returns:
    --------
    dict : condition_results (name -> params, timeseries, color, linestyle, marker)
    """
    global LAST_ALBM_DATA, LAST_EDDY_DATA, LAST_FLUX_DATA, LAST_OBS_TARGETS, LAST_CONDITION_RESULTS

    mc_start = time.time()
    print("\n" + "="*70)
    print(f"[{time.strftime('%H:%M:%S')}] MULTI-CONDITION OPTIMIZATION COMPARISON")
    print(f"[{time.strftime('%H:%M:%S')}] Optimizer workers: {OPT_WORKERS}")
    print("="*70)

    condition_results = {}
    total_conditions = len(CONDITIONS)

    for idx, (condition_name, config) in enumerate(CONDITIONS.items()):
        progress = (idx + 1) / total_conditions
        bar_length = 40
        filled = int(bar_length * progress)
        bar = '█' * filled + '░' * (bar_length - filled)
        print(f"\r[{time.strftime('%H:%M:%S')}] [{bar}] {progress*100:.0f}% - {condition_name}", end='', flush=True)

        if config['run_optimization']:
            # Check if this is a temperature-based optimization
            if config.get('use_temp_model', False):
                result = run_optimization_temp(
                    flux_data=flux_data,
                    obs_targets=obs_targets,
                    target_toggles=config['target_toggles'],
                    normalize_by_std=NORMALIZE_BY_STD,
                    maxiter=OPT_MAXITER,
                    tol=OPT_TOL,
                    popsize=OPT_POPSIZE,
                    workers=OPT_WORKERS,
                    verbose=OPT_VERBOSE,
                    seed=OPT_SEED
                )
                params = {
                    'm': result.m,
                    'b': result.b,
                    'alpha_mo': result.alpha_mo,
                    'alpha_op': result.alpha_op
                }
                timeseries = compute_isotope_timeseries_temp(params, flux_data)
            else:
                result = run_optimization(
                    flux_data=flux_data,
                    obs_targets=obs_targets,
                    target_toggles=config['target_toggles'],
                    normalize_by_std=NORMALIZE_BY_STD,
                    maxiter=OPT_MAXITER,
                    tol=OPT_TOL,
                    popsize=OPT_POPSIZE,
                    workers=OPT_WORKERS,
                    verbose=OPT_VERBOSE,
                    seed=OPT_SEED
                )
                params = {
                    'alpha_am': result.alpha_am,
                    'alpha_hm': result.alpha_hm,
                    'alpha_mo': result.alpha_mo,
                    'alpha_op': result.alpha_op,
                    'f_am': result.f_am,
                    'C13_POM': result.C13_POM
                }
                timeseries = compute_isotope_timeseries(params, flux_data)
        else:
            params = get_default_params()
            timeseries = compute_isotope_timeseries(params, flux_data)

        condition_results[condition_name] = {
            'params': params,
            'timeseries': timeseries,
            'color': config['color'],
            'linestyle': config['linestyle'],
            'marker': config['marker']
        }

    opt_elapsed = time.time() - mc_start
    print(f"\n[{time.strftime('%H:%M:%S')}] All optimizations complete ({opt_elapsed:.1f}s)")
    LAST_CONDITION_RESULTS = condition_results
    LAST_FLUX_DATA = flux_data
    LAST_OBS_TARGETS = obs_targets
    if albm_data is not None:
        LAST_ALBM_DATA = albm_data
    if eddy_data is not None:
        LAST_EDDY_DATA = eddy_data

    # Report resulting values (per condition)
    print("\n" + "="*100)
    print("OPTIMIZATION RESULTS (by condition)")
    print("="*100)
    for condition_name, result in condition_results.items():
        p = result['params']
        config = CONDITIONS[condition_name]
        
        if config.get('use_temp_model', False):
            # Temperature-based model
            sed_mean = result['timeseries'].get('C13_sed_prod_mean', np.mean(result['timeseries']['C13_sed_prod']))
            print(f"\n  {condition_name}:")
            print(f"    m={p['m']:.4f} °C⁻¹  b={p['b']:.2f} ‰  α_mo={p['alpha_mo']:.4f}  α_op={p['alpha_op']:.4f}")
            print(f"    δ¹³C_sed (mean)={sed_mean:.2f}‰")
        else:
            # Standard fractionation model
            sed = result['timeseries'].get('C13_sed_prod', result['timeseries'].get('C13_sed_prod_mean', -25.0))
            if 'f_am' in p:
                f_hm = 1.0 - p['f_am']
                print(f"\n  {condition_name}:")
                print(f"    α_am={p['alpha_am']:.4f}  α_hm={p['alpha_hm']:.4f}  α_mo={p['alpha_mo']:.4f}  α_op={p['alpha_op']:.4f}")
                print(f"    f_am={p['f_am']:.4f}  f_hm={f_hm:.4f}  C13_POM={p['C13_POM']:.2f}‰  δ¹³C_sed={sed:.2f}‰")
            else:
                print(f"\n  {condition_name}:")
                print(f"    (Initial guess - using default parameters)")
    print()

    # Ensure output directory exists before writing figures and statistics.
    os.makedirs(figs_dir, exist_ok=True)
    # Statistics: model mean ± std in summer, winter, and all overlap; comparison to observations
    print_overlap_statistics(condition_results, flux_data, obs_targets)

    # Paper statistics (numeric only, no plots). Print before any optional plot display.
    rerun_paper_statistics(
        condition_results=condition_results,
        flux_data=flux_data,
        obs_targets=obs_targets,
        albm_data=albm_data,
        eddy_data=eddy_data,
        figs_dir=figs_dir,
    )

    # Create comparison figure (plot in flux_plots)
    plot_optimization_comparison(
        condition_results=condition_results,
        flux_data=flux_data,
        obs_targets=obs_targets,
        figs_dir=figs_dir,
        show=show
    )
    plot_dissolved_pool_timeseries(
        condition_results=condition_results,
        flux_data=flux_data,
        figs_dir=figs_dir,
        show=show
    )

    # Violin plots of overlap periods: Initial guess vs Temperature-Based (dissolved and emitted)
    if plot_overlap_violins is not None:
        plot_overlap_violins(
            condition_results=condition_results,
            flux_data=flux_data,
            obs_targets=obs_targets,
            figs_dir=figs_dir,
            show=show
        )

    # Print summary table
    total_elapsed = time.time() - mc_start
    print("\n" + "="*100)
    print(f"[{time.strftime('%H:%M:%S')}] OPTIMIZATION COMPARISON SUMMARY (total: {total_elapsed:.1f}s)")
    print("="*100)
    # Print header - check if we have temperature-based models
    has_temp = any(CONDITIONS[name].get('use_temp_model', False) for name in condition_results.keys())
    if has_temp:
        print(f"{'Condition':<25} {'m (°C⁻¹)':<12} {'b (‰)':<10} {'α_mo':<8} {'α_op':<8} {'δ¹³C_sed':<10}")
        print("-"*100)
        for condition_name, result in condition_results.items():
            p = result['params']
            config = CONDITIONS[condition_name]
            if config.get('use_temp_model', False):
                sed_mean = result['timeseries'].get('C13_sed_prod_mean', np.mean(result['timeseries']['C13_sed_prod']))
                print(f"{condition_name:<25} {p['m']:<12.4f} {p['b']:<10.2f} {p['alpha_mo']:<8.4f} "
                      f"{p['alpha_op']:<8.4f} {sed_mean:<10.2f}")
            else:
                # For non-temp models, show N/A or standard params
                sed = result['timeseries'].get('C13_sed_prod', result['timeseries'].get('C13_sed_prod_mean', -25.0))
                if 'f_am' in p:
                    print(f"{condition_name:<25} {'N/A':<12} {'N/A':<10} {p.get('alpha_mo', 0):<8.4f} "
                          f"{p.get('alpha_op', 0):<8.4f} {sed:<10.2f}")
                else:
                    print(f"{condition_name:<25} {'N/A':<12} {'N/A':<10} {'N/A':<8} {'N/A':<8} {'N/A':<10}")
    else:
        print(f"{'Condition':<25} {'α_am':<8} {'α_hm':<8} {'α_mo':<8} {'α_op':<8} {'f_am':<8} {'f_hm':<8} {'C13_POM':<10} {'δ¹³C_sed':<10}")
        print("-"*100)
        for condition_name, result in condition_results.items():
            p = result['params']
            sed = result['timeseries'].get('C13_sed_prod', result['timeseries'].get('C13_sed_prod_mean', -25.0))
            if 'f_am' in p:
                f_hm = 1.0 - p['f_am']
                print(f"{condition_name:<25} {p['alpha_am']:<8.4f} {p['alpha_hm']:<8.4f} {p['alpha_mo']:<8.4f} "
                      f"{p['alpha_op']:<8.4f} {p['f_am']:<8.4f} {f_hm:<8.4f} {p['C13_POM']:<10.2f} {sed:<10.2f}")
            else:
                print(f"{condition_name:<25} {'N/A':<8} {'N/A':<8} {'N/A':<8} {'N/A':<8} {'N/A':<8} {'N/A':<8} {'N/A':<10} {sed:<10.2f}")
    print("="*100 + "\n")

    return condition_results


def main():
    global LAST_ALBM_DATA, LAST_EDDY_DATA, LAST_FLUX_DATA, LAST_OBS_TARGETS

    script_start_time = time.time()
    source_file = globals().get("__file__", "<interactive input>")
    paper_stats_path = os.path.join(FIGS_DIR, 'paper_statistics_flux_isotope.csv')

    print("\n" + "="*70)
    print(f"[{time.strftime('%H:%M:%S')}] MULTI-CONDITION ISOTOPE OPTIMIZATION")
    print(f"[{time.strftime('%H:%M:%S')}] Source: {source_file}", flush=True)
    print(f"[{time.strftime('%H:%M:%S')}] Paper statistics: "
          f"{'enabled' if compute_and_print_paper_statistics is not None else 'disabled'}", flush=True)
    print(f"[{time.strftime('%H:%M:%S')}] Paper statistics print after all optimizations finish.", flush=True)
    print(f"[{time.strftime('%H:%M:%S')}] Paper statistics CSV: {paper_stats_path}", flush=True)
    print(f"[{time.strftime('%H:%M:%S')}] Show plots: {SHOW_PLOTS}", flush=True)
    print("="*70)

    albm_data, eddy_data = get_processed_data(
        cache_path=CACHE_PATH,
        results_dir=RESULTS_DIR,
        date_range=DATE_RANGE,
        eddy_flux_file=EDDY_FLUX_FILE,
        use_cache=USE_CACHE,
        save_cache=False
    )
    flux_data = get_flux_data_dict(albm_data)
    obs_targets = get_default_obs_targets()
    LAST_ALBM_DATA = albm_data
    LAST_EDDY_DATA = eddy_data
    LAST_FLUX_DATA = flux_data
    LAST_OBS_TARGETS = obs_targets

    os.makedirs(FIGS_DIR, exist_ok=True)

    condition_results = run_multi_condition_comparison(
        flux_data=flux_data,
        obs_targets=obs_targets,
        figs_dir=FIGS_DIR,
        show=SHOW_PLOTS,
        albm_data=albm_data,
        eddy_data=eddy_data,
    )

    total_runtime = time.time() - script_start_time
    hours, remainder = divmod(total_runtime, 3600)
    minutes, seconds = divmod(remainder, 60)
    print(f"\n[{time.strftime('%H:%M:%S')}] Total runtime: {int(hours):02d}:{int(minutes):02d}:{seconds:05.2f} (hh:mm:ss)")
    print(f"               {total_runtime:.2f} seconds")

    return condition_results


if __name__ == '__main__':
    main()
