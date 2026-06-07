"""
Multi-start optimization (uniqueness / identifiability check).

Runs differential evolution many times with independent random seeds for the static
fractionation model and/or the temperature-based model (default: both), then
summarizes cost and parameter spread and saves tables plus distribution figures.

Usage (from repo root):
    python src/manuscript_scripts/multi_start_optimization.py
    python src/manuscript_scripts/multi_start_optimization.py --n-runs 10 --model static

From Jupyter, call main() (uses defaults) or main(['--n-runs', '2', ...]) so arguments are
explicit; the kernel's sys.argv is not used.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, *args, **kwargs):  # type: ignore[misc]
        return iterable

from data_loader import (
    get_processed_data,
    get_flux_data_dict,
    DEFAULT_CACHE_PATH,
    DEFAULT_EDDY_FLUX_FILE,
)
from isotope_model import (
    run_optimization,
    run_optimization_temp,
    get_default_obs_targets,
    OptimizationResult,
    OptimizationResultTemp,
    delta_to_ratio,
    ratio_to_delta,
)

# Match multi_condition_optimization.py
RESULTS_DIR = 'ALBM data'
DATE_RANGE = '20210101_20250101'
EDDY_FLUX_FILE = DEFAULT_EDDY_FLUX_FILE
CACHE_PATH = DEFAULT_CACHE_PATH
USE_CACHE = True

OPT_MAXITER = 200
OPT_TOL = 0.01
OPT_POPSIZE = 15
OPT_WORKERS = -1
NORMALIZE_BY_STD = False

# Same toggles as "All Observations" / "Temperature-Based" in CONDITIONS
ALL_OBS_TARGET_TOGGLES: Dict[str, bool] = {
    'dissolved_summer': True,
    'dissolved_winter': True,
    'emissions_summer': True,
    'emissions_winter': True,
}


def _make_seed_list(n_runs: int, master_seed: int) -> List[int]:
    """Reproducible distinct integer seeds for each DE run."""
    rng = np.random.default_rng(master_seed)
    # Large positive ints; scipy's DE accepts these as seed
    return [int(x) for x in rng.integers(1, 2**31 - 1, size=n_runs, dtype=np.int64)]


def _set_run_pbar_postfix(pbar: Any, seed: int, row: Dict[str, Any]) -> None:
    """Update tqdm postfix with last cost or error (no-op if tqdm fallback)."""
    if not hasattr(pbar, 'set_postfix'):
        return
    cost = row.get('cost')
    if isinstance(cost, (int, float)) and np.isfinite(cost):
        pbar.set_postfix(seed=seed, cost=f'{float(cost):.4g}')
    elif row.get('error'):
        pbar.set_postfix(seed=seed, err=str(row['error'])[:28])
    else:
        pbar.set_postfix(seed=seed)


def _run_static_batch(
    flux_data: Dict[str, Any],
    obs_targets: Dict[str, Dict[str, float]],
    seeds: Sequence[int],
    *,
    maxiter: int,
    tol: float,
    popsize: int,
    workers: int,
    normalize_by_std: bool,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    run_iter = tqdm(
        enumerate(seeds),
        total=len(seeds),
        desc='[static] DE runs',
        unit='run',
    )
    for run_idx, seed in run_iter:
        t0 = time.time()
        try:
            r: OptimizationResult = run_optimization(
                flux_data=flux_data,
                obs_targets=obs_targets,
                target_toggles=ALL_OBS_TARGET_TOGGLES,
                normalize_by_std=normalize_by_std,
                maxiter=maxiter,
                tol=tol,
                popsize=popsize,
                workers=workers,
                verbose=False,
                seed=int(seed),
            )
            R_POM = delta_to_ratio(r.C13_POM)
            R_sed = r.f_hm * (R_POM / r.alpha_hm) + r.f_am * (R_POM / r.alpha_am)
            c13_sed = ratio_to_delta(R_sed)
            row = {
                'model': 'static',
                'run': run_idx,
                'seed': int(seed),
                'cost': r.cost,
                'success': r.success,
                'n_evaluations': r.n_evaluations,
                'alpha_am': r.alpha_am,
                'alpha_hm': r.alpha_hm,
                'alpha_mo': r.alpha_mo,
                'alpha_op': r.alpha_op,
                'f_am': r.f_am,
                'f_hm': r.f_hm,
                'C13_POM': r.C13_POM,
                'C13_sed_prod': c13_sed,
                'elapsed_s': time.time() - t0,
                'error': '',
            }
        except Exception as e:
            row = {
                'model': 'static',
                'run': run_idx,
                'seed': int(seed),
                'cost': np.nan,
                'success': False,
                'n_evaluations': 0,
                'alpha_am': np.nan,
                'alpha_hm': np.nan,
                'alpha_mo': np.nan,
                'alpha_op': np.nan,
                'f_am': np.nan,
                'f_hm': np.nan,
                'C13_POM': np.nan,
                'C13_sed_prod': np.nan,
                'elapsed_s': time.time() - t0,
                'error': str(e),
            }
        rows.append(row)
        _set_run_pbar_postfix(run_iter, int(seed), row)
    return rows


def _run_temp_batch(
    flux_data: Dict[str, Any],
    obs_targets: Dict[str, Dict[str, float]],
    seeds: Sequence[int],
    *,
    maxiter: int,
    tol: float,
    popsize: int,
    workers: int,
    normalize_by_std: bool,
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    run_iter = tqdm(
        enumerate(seeds),
        total=len(seeds),
        desc='[temp] DE runs',
        unit='run',
    )
    for run_idx, seed in run_iter:
        t0 = time.time()
        try:
            r: OptimizationResultTemp = run_optimization_temp(
                flux_data=flux_data,
                obs_targets=obs_targets,
                target_toggles=ALL_OBS_TARGET_TOGGLES,
                normalize_by_std=normalize_by_std,
                maxiter=maxiter,
                tol=tol,
                popsize=popsize,
                workers=workers,
                verbose=False,
                seed=int(seed),
            )
            row = {
                'model': 'temp',
                'run': run_idx,
                'seed': int(seed),
                'cost': r.cost,
                'success': r.success,
                'n_evaluations': r.n_evaluations,
                'm': r.m,
                'b': r.b,
                'alpha_mo': r.alpha_mo,
                'alpha_op': r.alpha_op,
                'elapsed_s': time.time() - t0,
                'error': '',
            }
        except Exception as e:
            row = {
                'model': 'temp',
                'run': run_idx,
                'seed': int(seed),
                'cost': np.nan,
                'success': False,
                'n_evaluations': 0,
                'm': np.nan,
                'b': np.nan,
                'alpha_mo': np.nan,
                'alpha_op': np.nan,
                'elapsed_s': time.time() - t0,
                'error': str(e),
            }
        rows.append(row)
        _set_run_pbar_postfix(run_iter, int(seed), row)
    return rows


def _print_summary(df: pd.DataFrame, label: str, cost_rtol: float, param_atol: float) -> None:
    print('\n' + '=' * 80)
    print(f'SUMMARY — {label}')
    print('=' * 80)
    cost_col = pd.to_numeric(df['cost'], errors='coerce')
    ok = df['success'] & np.isfinite(cost_col)
    finite = np.isfinite(cost_col)
    usable = ok if ok.any() else finite
    costs = cost_col[usable].astype(float)
    if costs.empty:
        print('No finite-cost runs.')
        return
    print(f"Successful runs: {ok.sum()} / {len(df)}")
    if not ok.any() and finite.any():
        print(f"Finite-cost runs used for summary despite success=False: {int(finite.sum())} / {len(df)}")
    print(f"Cost: min={costs.min():.6g}  median={costs.median():.6g}  max={costs.max():.6g}")
    cmin = costs.min()
    tol_c = max(cost_rtol * max(abs(cmin), 1e-12), 1e-15)
    near_best = usable & ((cost_col - cmin).abs() <= tol_c)
    print(f"Runs within cost tol (|c - c_min| <= {tol_c:.3g}): {int(near_best.sum())}")

    model = df['model'].iloc[0]
    if model == 'static':
        pcols = ['alpha_am', 'alpha_hm', 'alpha_mo', 'alpha_op', 'f_am', 'f_hm', 'C13_POM', 'C13_sed_prod']
    else:
        pcols = ['m', 'b', 'alpha_mo', 'alpha_op']

    best_idx = cost_col[usable].idxmin()
    ref = df.loc[best_idx, pcols].astype(float)
    print(f"Reference row (best cost): run index {df.loc[best_idx, 'run']}")
    for col in pcols:
        v = df.loc[usable, col].astype(float)
        dev = (v - ref[col]).abs().max()
        print(f"  {col}: max |x - x_best| over summarized runs = {dev:.6g}")
    within = usable & near_best
    for col in pcols:
        within &= (pd.to_numeric(df[col], errors='coerce') - float(ref[col])).abs() <= param_atol
    print(f"Runs within cost tol AND all params within atol={param_atol}: {int(within.sum())}")
    print('=' * 80 + '\n')


def _summary_rows(df: pd.DataFrame, label: str, cost_rtol: float, param_atol: float) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []

    def add(statistic_id: str, description: str, value: Any, units: str = "", notes: str = "") -> None:
        if isinstance(value, (np.floating, np.integer)):
            value = value.item()
        rows.append({
            'source_file': 'multi_start_optimization.py',
            'model': label,
            'statistic_id': statistic_id,
            'description': description,
            'value': value,
            'units': units,
            'manuscript_line': '447',
            'notes': notes,
        })

    add('run_count', 'Number of multi-start optimization runs', len(df), 'count')
    if df.empty or 'success' not in df or 'cost' not in df:
        add('successful_finite_runs', 'Number of successful finite-cost runs', 0, 'count')
        return rows

    cost_col = pd.to_numeric(df['cost'], errors='coerce')
    ok = df['success'].astype(bool) & np.isfinite(cost_col)
    finite = np.isfinite(cost_col)
    usable = ok if ok.any() else finite
    costs = cost_col[usable].astype(float)
    add('successful_finite_runs', 'Number of successful finite-cost runs', int(ok.sum()), 'count')
    add('finite_cost_runs', 'Number of finite-cost runs', int(finite.sum()), 'count')
    if not ok.any() and finite.any():
        add('summary_basis', 'Rows used for cost and parameter summaries', 'finite_cost_runs', '',
            'All success flags were false, so finite-cost rows were summarized.')
    else:
        add('summary_basis', 'Rows used for cost and parameter summaries', 'successful_finite_runs')
    if costs.empty:
        return rows

    cmin = float(costs.min())
    tol_c = max(cost_rtol * max(abs(cmin), 1e-12), 1e-15)
    near_best = usable & ((cost_col - cmin).abs() <= tol_c)
    add('cost_min', 'Minimum optimization cost', cmin)
    add('cost_median', 'Median optimization cost among summarized runs', float(costs.median()))
    add('cost_max', 'Maximum optimization cost among summarized runs', float(costs.max()))
    add('cost_tolerance', 'Cost tolerance used to define near-best runs', tol_c)
    add('near_best_runs', 'Runs within the cost tolerance of the best run', int(near_best.sum()), 'count')

    model = str(df['model'].iloc[0]) if 'model' in df else label
    if model == 'static':
        pcols = ['alpha_am', 'alpha_hm', 'alpha_mo', 'alpha_op', 'f_am', 'f_hm', 'C13_POM', 'C13_sed_prod']
    else:
        pcols = ['m', 'b', 'alpha_mo', 'alpha_op']
    pcols = [col for col in pcols if col in df.columns]

    best_idx = cost_col[usable].idxmin()
    best_run = df.loc[best_idx, 'run'] if 'run' in df.columns else best_idx
    add('best_run', 'Run index with the minimum cost', int(best_run) if pd.notna(best_run) else best_idx)
    ref = df.loc[best_idx, pcols].astype(float)
    within = usable & near_best
    for col in pcols:
        vals = pd.to_numeric(df.loc[usable, col], errors='coerce').astype(float)
        dev = (vals - ref[col]).abs().max()
        add(f'{col}_mean', f'Mean {col} among summarized runs', float(vals.mean()))
        add(f'{col}_std', f'Standard deviation of {col} among summarized runs', float(vals.std(ddof=1)))
        add(f'{col}_min', f'Minimum {col} among summarized runs', float(vals.min()))
        add(f'{col}_max', f'Maximum {col} among summarized runs', float(vals.max()))
        add(f'{col}_max_abs_dev_from_best',
            f'Maximum absolute deviation of {col} from best-cost run among summarized runs',
            float(dev))
        within &= (pd.to_numeric(df[col], errors='coerce') - float(ref[col])).abs() <= param_atol

    add('near_best_and_param_close_runs',
        f'Runs within cost tolerance and all parameters within atol={param_atol}',
        int(within.sum()), 'count')
    return rows


def _save_summary_statistics(
    df_static: Optional[pd.DataFrame],
    df_temp: Optional[pd.DataFrame],
    output_dir: str,
    cost_rtol: float,
    param_atol: float,
    filename: str = 'multi_start_summary_statistics.csv',
) -> Optional[str]:
    rows: List[Dict[str, Any]] = []
    if df_static is not None and len(df_static):
        rows.extend(_summary_rows(df_static, 'static', cost_rtol, param_atol))
    if df_temp is not None and len(df_temp):
        rows.extend(_summary_rows(df_temp, 'temp', cost_rtol, param_atol))
    if not rows:
        return None
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, filename)
    pd.DataFrame(rows).to_csv(out_path, index=False, encoding='utf-8-sig')
    print(f"Saved: {out_path}")
    return out_path


def _save_tables(rows: List[Dict[str, Any]], output_dir: str, base_name: str) -> Tuple[str, str]:
    os.makedirs(output_dir, exist_ok=True)
    df = pd.DataFrame(rows)
    csv_path = os.path.join(output_dir, f'{base_name}.csv')
    json_path = os.path.join(output_dir, f'{base_name}.json')
    df.to_csv(csv_path, index=False)
    # JSON: records, numpy-friendly
    df.to_json(json_path, orient='records', indent=2)
    print(f"Saved: {csv_path}")
    print(f"Saved: {json_path}")
    return csv_path, json_path


def _hist_panel(ax: plt.Axes, values: np.ndarray, title: str, xlabel: str) -> None:
    v = values[np.isfinite(values)]
    if v.size == 0:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title)
        return
    ax.hist(v, bins=max(5, min(20, len(v))), color='steelblue', edgecolor='black', alpha=0.85)
    if v.size >= 15:
        try:
            from scipy.stats import gaussian_kde
            xs = np.linspace(v.min(), v.max(), 200)
            kde = gaussian_kde(v)
            ax2 = ax.twinx()
            ax2.plot(xs, kde(xs), color='darkorange', linewidth=1.5)
            ax2.set_ylabel('KDE', fontsize=8)
            ax2.tick_params(labelsize=7)
        except Exception:
            pass
    ax.set_title(title, fontsize=10)
    ax.set_xlabel(xlabel, fontsize=9)


def plot_static_distributions(df: pd.DataFrame, out_path: str, show: bool) -> None:
    df = df[np.isfinite(df['cost'])].copy()
    params = [
        ('alpha_am', r'$\alpha_{\mathrm{AM}}$'),
        ('alpha_hm', r'$\alpha_{\mathrm{HM}}$'),
        ('alpha_mo', r'$\alpha_{\mathrm{MO}}$'),
        ('alpha_op', r'$\alpha_{\mathrm{OP}}$'),
        ('f_am', r'$f_{\mathrm{AM}}$'),
        ('f_hm', r'$f_{\mathrm{HM}}$'),
        ('C13_POM', r'$\delta^{13}$C POM (‰)'),
        ('C13_sed_prod', r'$\delta^{13}$C sed. prod. (‰)'),
    ]
    fig, axes = plt.subplots(3, 3, figsize=(11, 9))
    axes = axes.ravel()
    for i, (col, title) in enumerate(params):
        _hist_panel(axes[i], df[col].values, title, col)
    for j in range(len(params), len(axes)):
        axes[j].axis('off')
    fig.suptitle('Static model: parameter distributions across seeds', fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_temp_distributions(df: pd.DataFrame, out_path: str, show: bool) -> None:
    df = df[np.isfinite(df['cost'])].copy()
    params = [
        ('m', r'$m$ (‰ °C$^{-1}$)'),
        ('b', r'$b$ (‰)'),
        ('alpha_mo', r'$\alpha_{\mathrm{MO}}$'),
        ('alpha_op', r'$\alpha_{\mathrm{OP}}$'),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    axes = axes.ravel()
    for i, (col, title) in enumerate(params):
        _hist_panel(axes[i], df[col].values, title, col)
    fig.suptitle('Temperature model: parameter distributions across seeds', fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_cost_distributions(
    df_static: Optional[pd.DataFrame],
    df_temp: Optional[pd.DataFrame],
    out_path: str,
    show: bool,
) -> None:
    n = int(df_static is not None) + int(df_temp is not None)
    if n == 0:
        return
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4))
    if n == 1:
        axes = [axes]
    idx = 0
    if df_static is not None and len(df_static):
        c = df_static['cost'].values
        c = c[np.isfinite(c)]
        axes[idx].hist(c, bins=max(5, min(15, len(c))), color='C0', edgecolor='black', alpha=0.85)
        axes[idx].set_title('Static model: cost')
        axes[idx].set_xlabel('Cost')
        idx += 1
    if df_temp is not None and len(df_temp):
        c = df_temp['cost'].values
        c = c[np.isfinite(c)]
        axes[idx].hist(c, bins=max(5, min(15, len(c))), color='C1', edgecolor='black', alpha=0.85)
        axes[idx].set_title('Temperature model: cost')
        axes[idx].set_xlabel('Cost')
    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(fig)


def _print_run_table(df: pd.DataFrame) -> None:
    cols = [c for c in df.columns if c not in ('error',)]
    print(df[cols].to_string(index=False))


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """
    Parse CLI arguments. When ``argv`` is None, uses ``sys.argv[1:]``.

    Uses ``parse_known_args`` so that unrecognised tokens (Jupyter kernel flags,
    ``-X``, ``-f kernel.json``, VS Code args, etc.) are silently ignored rather
    than causing argparse to exit with code 2.

    From a notebook call ``main()`` (defaults) or ``main(['--n-runs', '5', ...])``.
    """
    p = argparse.ArgumentParser(description='Multi-start DE optimization for uniqueness check')
    p.add_argument('--n-runs', type=int, default=100, help='Number of independent runs per model')
    p.add_argument('--master-seed', type=int, default=42, help='Seed used to generate per-run seeds')
    p.add_argument('--model', choices=('both', 'static', 'temp'), default='both',
                   help='Which optimizer(s) to run (default: both)')
    p.add_argument('--maxiter', type=int, default=OPT_MAXITER)
    p.add_argument('--tol', type=float, default=OPT_TOL)
    p.add_argument('--popsize', type=int, default=OPT_POPSIZE)
    p.add_argument('--workers', type=int, default=OPT_WORKERS)
    p.add_argument('--normalize-by-std', action='store_true', help='Normalize cost by obs std')
    p.add_argument('--output-dir', type=str, default='figs/multi_start',
                   help='Directory for CSV, JSON, and figures')
    p.add_argument('--cost-rtol', type=float, default=1e-4, help='Relative cost tolerance for "near best"')
    p.add_argument('--param-atol', type=float, default=1e-3, help='Absolute param tolerance vs best run')
    p.add_argument('--show', action='store_true', help='Show matplotlib windows')
    args, _ = p.parse_known_args(argv)
    return args


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    repo_root = PROJECT_ROOT
    cwd = os.getcwd()
    # Project root: ./ALBM data ; from this script: ../../ALBM data
    cand = os.path.join(cwd, RESULTS_DIR)
    if os.path.isdir(cand):
        results_dir = cand
    else:
        cand2 = os.path.join(repo_root, RESULTS_DIR)
        results_dir = cand2 if os.path.isdir(cand2) else cand

    eddy_path = EDDY_FLUX_FILE
    if not os.path.isfile(eddy_path):
        alt = os.path.join(repo_root, EDDY_FLUX_FILE)
        if os.path.isfile(alt):
            eddy_path = alt

    cache_path = CACHE_PATH
    if not os.path.isfile(cache_path):
        alt = os.path.join(repo_root, CACHE_PATH)
        if os.path.isfile(alt):
            cache_path = alt

    print(f"[{time.strftime('%H:%M:%S')}] Loading data (results_dir={results_dir})...")
    albm_data, _eddy = get_processed_data(
        cache_path=cache_path,
        results_dir=results_dir,
        date_range=DATE_RANGE,
        eddy_flux_file=eddy_path,
        use_cache=USE_CACHE,
        save_cache=False,
    )
    flux_data = get_flux_data_dict(albm_data)
    obs_targets = get_default_obs_targets()

    seeds = _make_seed_list(args.n_runs, args.master_seed)
    print(f"[{time.strftime('%H:%M:%S')}] {args.n_runs} seeds generated (master={args.master_seed})")

    all_rows: List[Dict[str, Any]] = []
    df_static: Optional[pd.DataFrame] = None
    df_temp: Optional[pd.DataFrame] = None

    if args.model in ('both', 'static'):
        print(f"\n[{time.strftime('%H:%M:%S')}] --- Static (fractionation) model: {args.n_runs} runs ---\n")
        static_rows = _run_static_batch(
            flux_data, obs_targets, seeds,
            maxiter=args.maxiter, tol=args.tol, popsize=args.popsize,
            workers=args.workers, normalize_by_std=args.normalize_by_std,
        )
        all_rows.extend(static_rows)
        df_static = pd.DataFrame(static_rows)
        _print_run_table(df_static)
        _print_summary(df_static, 'static', args.cost_rtol, args.param_atol)

    if args.model in ('both', 'temp'):
        print(f"\n[{time.strftime('%H:%M:%S')}] --- Temperature-based model: {args.n_runs} runs ---\n")
        temp_rows = _run_temp_batch(
            flux_data, obs_targets, seeds,
            maxiter=args.maxiter, tol=args.tol, popsize=args.popsize,
            workers=args.workers, normalize_by_std=args.normalize_by_std,
        )
        all_rows.extend(temp_rows)
        df_temp = pd.DataFrame(temp_rows)
        _print_run_table(df_temp)
        _print_summary(df_temp, 'temp', args.cost_rtol, args.param_atol)

    out_dir = args.output_dir
    if not os.path.isabs(out_dir):
        out_dir = os.path.join(repo_root, out_dir)
    os.makedirs(out_dir, exist_ok=True)
    ts = time.strftime('%Y%m%d_%H%M%S')
    base = f'multi_start_{ts}'
    _save_tables(all_rows, out_dir, base)

    # Figures
    if df_static is not None and len(df_static):
        plot_static_distributions(
            df_static,
            os.path.join(out_dir, 'multi_start_static_params.png'),
            args.show,
        )
        print(f"Saved: {os.path.join(out_dir, 'multi_start_static_params.png')}")
    if df_temp is not None and len(df_temp):
        plot_temp_distributions(
            df_temp,
            os.path.join(out_dir, 'multi_start_temp_params.png'),
            args.show,
        )
        print(f"Saved: {os.path.join(out_dir, 'multi_start_temp_params.png')}")
    plot_cost_distributions(
        df_static,
        df_temp,
        os.path.join(out_dir, 'multi_start_costs.png'),
        args.show,
    )
    print(f"Saved: {os.path.join(out_dir, 'multi_start_costs.png')}")
    _save_summary_statistics(
        df_static,
        df_temp,
        out_dir,
        args.cost_rtol,
        args.param_atol,
    )

    print(f"\n[{time.strftime('%H:%M:%S')}] Done. Outputs in: {out_dir}")
    return 0


if __name__ == '__main__':
    # Do not pass sys.argv[1:] here: in Jupyter, that slice is still the kernel's -m/-f
    # tail. main() uses _argv_for_argparse(None) which maps real ``python script.py`` to
    # sys.argv[1:] and notebooks to defaults.
    sys.exit(main())
