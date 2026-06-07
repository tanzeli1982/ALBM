"""
ALBM plotting module: flux comparison, water-column and sediment heatmaps,
isotope results, and multi-condition optimization comparison.

All plot functions can be called individually. Use plot_all() to generate
the standard set. Control log output via logging (e.g. logging.getLogger('flux_plots').setLevel(logging.WARNING)).

Usage:
    from data_loader import load_albm_data, load_eddy_flux_data
    from flux_plots import plot_all, plot_eddy_comparison, plot_isotope_results

    albm_data = load_albm_data()
    eddy_data = load_eddy_flux_data()
    plot_all(albm_data, eddy_data, figs_dir='figs')
"""

import logging
import os
import re
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import BoundaryNorm, LogNorm
from matplotlib.ticker import LogFormatterSciNotation, ScalarFormatter
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerLine2D
from typing import Optional, Dict, Any, List, Tuple, Union

from data_loader import ALBMData, EddyFluxData

logger = logging.getLogger(__name__)

# =============================================================================
# Colorblind-friendly palette (Okabe-Ito)
# =============================================================================
CB_ORANGE = '#E69F00'
CB_SKY_BLUE = '#56B4E9'
CB_GREEN = '#009E73'
CB_YELLOW = '#F0E442'
CB_BLUE = '#0072B2'
CB_VERMILLION = '#D55E00'
CB_PURPLE = '#CC79A7'
CB_BLACK = '#000000'
CB_GRAY = '#999999'

# Observation colors (summer reuses CB_PURPLE for consistency)
OBS_SUMMER_COLOR = CB_PURPLE
OBS_WINTER_COLOR = '#661100'

# Sediment depth: model uses 25 m total, 41 Z-nodes (0.625 m spacing)
HSED_M = 25.0
SED_DEPTH_NODES_TYPICAL = 41
SED_DEPTH_SPACING_M = HSED_M / (SED_DEPTH_NODES_TYPICAL - 1)  # 0.625

ISOTOPE_Y_MIN = -85
ISOTOPE_Y_MAX = -10

# Legend convention: LEGEND_UPPER_RIGHT for most plots; use LEGEND_LOWER_RIGHT for dense panels
LEGEND_UPPER_RIGHT = 'upper right'
LEGEND_LOWER_RIGHT = 'lower right'


def _apply_figure_style(fig: plt.Figure, axes: Union[plt.Axes, List[plt.Axes], None] = None) -> None:
    """Set transparent figure background and transparent axis face for all axes."""
    fig.patch.set_alpha(0.0)
    if axes is None:
        return
    ax_list = np.atleast_1d(axes).flatten().tolist()
    for ax in ax_list:
        ax.set_facecolor('none')


def _sediment_depth_array(n_sed_depths: int) -> np.ndarray:
    """Return depth array (m) for sediment layers. Uses 0.625 m spacing for Z=41 (25 m), else 0.1 m."""
    if n_sed_depths <= 0:
        return np.array([0.0])
    if n_sed_depths > 5:
        return np.linspace(0.0, SED_DEPTH_SPACING_M * (n_sed_depths - 1), n_sed_depths)
    return np.linspace(0.0, 0.1 * (n_sed_depths - 1), n_sed_depths)


class _LegendSeparatorLine(Line2D):
    """Marker type so only the legend separator gets the full-width handler."""


class _FullWidthLineHandler(HandlerLine2D):
    """Draw a horizontal line across the full legend width; long line in handlebox coords, clipped to legend frame."""
    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        big = 10000.0
        line = Line2D(
            [-big, width + big], [height / 2.0, height / 2.0],
            color=orig_handle.get_color(),
            linewidth=max(1.0, orig_handle.get_linewidth()),
            transform=trans,
            solid_capstyle='butt',
            clip_on=True,
        )
        # Clip to legend frame so the line doesn't extend across the figure
        line.set_clip_path(legend.get_frame())
        return [line]


def _subplot_label(i: int) -> str:
    """Return subplot label 'a)', 'b)', ... for index i (0-based)."""
    return f'{chr(97 + i)})'


def _observation_season_ranges(year: int, max_year: int) -> List[Tuple[pd.Timestamp, pd.Timestamp, str]]:
    """
    Yield (start, end, season_key) for each observation period in the year.
    season_key: 'summer' | 'winter_dec' | 'winter_jan'
    Summer: Jun 21–Sep 22; Winter: Dec 21–Dec 31 and Jan 1–Mar 20 (next year).
    """
    out = [
        (pd.Timestamp(f'{year}-06-21'), pd.Timestamp(f'{year}-09-22'), 'summer'),
        (pd.Timestamp(f'{year}-12-21'), pd.Timestamp(f'{year}-12-31'), 'winter_dec'),
    ]
    if year < max_year:
        out.append((pd.Timestamp(f'{year+1}-01-01'), pd.Timestamp(f'{year+1}-03-20'), 'winter_jan'))
    return out


def _format_axes_dates(axes: Union[plt.Axes, List[plt.Axes]], rotation: int = 0) -> None:
    """Apply year-based date formatting and optional tick rotation to axes."""
    ax_list = np.atleast_1d(axes).flatten().tolist()
    for ax in ax_list:
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(mdates.MonthLocator())
        plt.setp(ax.xaxis.get_majorticklabels(), fontsize=10, rotation=rotation)
        plt.setp(ax.xaxis.get_minorticklabels(), rotation=rotation)


def _no_data_figure(message: str = 'No data') -> plt.Figure:
    """Return a minimal figure with a single axis showing message (for consistent API instead of None)."""
    fig, ax = plt.subplots(figsize=(4, 2))
    _apply_figure_style(fig, ax)
    ax.text(0.5, 0.5, message, ha='center', va='center', fontsize=14, transform=ax.transAxes)
    ax.axis('off')
    return fig


def setup_axis_style(ax: plt.Axes, linewidth: int = 2) -> None:
    """Apply standard axis spine and tick styling. Call after other axis setup to avoid overrides."""
    ax.spines['bottom'].set_linewidth(linewidth)
    ax.spines['left'].set_linewidth(linewidth)
    ax.spines['right'].set_linewidth(linewidth)
    ax.spines['top'].set_linewidth(linewidth)
    ax.tick_params(axis='both', which='major', labelsize=12, width=2)


def plot_water_column_heatmaps(
    albm_data: ALBMData,
    figs_dir: str = 'figs',
    show: bool = True
) -> plt.Figure:
    """
    Plot consolidated water column heatmaps (CH4 concentration, oxidation, production, exchange).
    
    Parameters:
    -----------
    albm_data : ALBMData
        ALBM data container
    figs_dir : str
        Directory to save figures
    show : bool
        Whether to display the plot
    
    Returns:
    --------
    fig : matplotlib Figure
    """
    logger.info("Creating consolidated water column heatmap subplot...")
    MW_CH4 = 16.04e-3
    time_numeric = mdates.date2num(albm_data.time_array)
    has_sed_temp = (hasattr(albm_data, 'sediment_temp_depth') and
                    albm_data.sediment_temp_depth is not None and
                    albm_data.sediment_temp_depth.ndim >= 2)
    n_rows = 4 if has_sed_temp else 3
    fig, axes = plt.subplots(n_rows, 1, figsize=(14, 4 * n_rows), sharex=True)
    axes = np.atleast_1d(axes).flatten().tolist()
    _apply_figure_style(fig, axes)

    consolidated_plots = [
        {'data': albm_data.ch4conc * MW_CH4, 'title': 'Dissolved CH₄ Concentration',
         'cbar_label': 'CH₄ Concentration (kg m⁻³)', 'cmap': 'plasma'},
        {'data': albm_data.ch4oxid * MW_CH4, 'title': 'CH₄ Oxidation Rate',
         'cbar_label': 'CH₄ Oxidation (kg m⁻³ s⁻¹)', 'cmap': 'RdPu'},
        {'data': albm_data.och4prod * MW_CH4, 'title': 'Oxic CH₄ Production',
         'cbar_label': 'Oxic CH₄ Production (kg m⁻³ s⁻¹)', 'cmap': 'cividis'},
    ]

    for i, (ax, plot_info) in enumerate(zip(axes, consolidated_plots)):
        plot_data = plot_info['data'].T
        im = ax.pcolormesh(time_numeric, albm_data.water_depth_array, plot_data,
                          shading='auto', cmap=plot_info['cmap'])
        ax.invert_yaxis()
        ax.set_ylim(4.1, 0)
        ax.xaxis_date()
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=3))
        ax.set_ylabel('Depth (m)', fontsize=14)
        ax.set_title(plot_info['title'], fontsize=16)
        ax.text(-0.08, 1.02, _subplot_label(i), transform=ax.transAxes, fontsize=16, fontweight='bold',
                verticalalignment='bottom', horizontalalignment='left')
        cbar = plt.colorbar(im, ax=ax, label=plot_info['cbar_label'])
        cbar.ax.tick_params(labelsize=10)
        setup_axis_style(ax)

    if has_sed_temp:
        ax = axes[3]
        n_sed_depths = albm_data.sediment_temp_depth.shape[1]
        sed_depths = _sediment_depth_array(n_sed_depths)
        hsed = float(sed_depths[-1])
        im = ax.pcolormesh(time_numeric, sed_depths, albm_data.sediment_temp_depth.T,
                          shading='nearest', cmap='coolwarm')
        ax.set_ylim(0, hsed)
        ax.invert_yaxis()
        ax.xaxis_date()
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=3))
        ax.set_ylabel('Sediment Depth (m)', fontsize=14)
        ax.set_title('Sediment Temperature by Depth', fontsize=16)
        ax.text(-0.08, 1.02, _subplot_label(3), transform=ax.transAxes, fontsize=16, fontweight='bold',
                verticalalignment='bottom', horizontalalignment='left')
        cbar = plt.colorbar(im, ax=ax, label='Temperature (°C)')
        cbar.ax.tick_params(labelsize=10)
        setup_axis_style(ax)

    axes[-1].set_xlabel('Year', fontsize=14)
    plt.tight_layout()
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, 'isoALBM_consolidated_water_column_heatmaps.png'),
                transparent=True, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    return fig


def _doy_mean_climatology(
    time_array: np.ndarray,
    data_tx: np.ndarray,
) -> np.ndarray:
    """
    Mean over calendar years by day-of-year (1–366). For each DOY, average all
    timesteps in the model record that fall on that DOY.

    Parameters
    ----------
    time_array : (n_time,)
    data_tx : (n_time, n_depth) depth-resolved series

    Returns
    -------
    mean_doy : (n_depth, 366) with NaN where that DOY never occurs
    """
    if data_tx.ndim != 2:
        raise ValueError("data_tx must be 2D (n_time, n_depth)")
    n_time, n_dep = data_tx.shape
    if n_time == 0:
        return np.full((n_dep, 366), np.nan)
    t = pd.DatetimeIndex(pd.to_datetime(np.asarray(time_array).ravel()))
    doy = t.dayofyear.to_numpy()
    out = np.full((n_dep, 366), np.nan, dtype=float)
    for d in range(1, 367):
        mask = doy == d
        if np.any(mask):
            out[:, d - 1] = np.nanmean(data_tx[mask, :], axis=0)
    return out


def _centers_to_edges_1d(centers: np.ndarray) -> np.ndarray:
    """
    Build monotonic bin edges for pcolormesh (shading='flat') from 1D cell centers.
    Length = len(centers) + 1.
    """
    c = np.asarray(centers, dtype=float).ravel()
    if c.size == 0:
        return np.array([0.0, 1.0])
    if c.size == 1:
        return np.array([c[0] - 0.5, c[0] + 0.5])
    d = np.diff(c)
    left = c[0] - d[0] / 2.0
    right = c[-1] + d[-1] / 2.0
    inner = (c[:-1] + c[1:]) / 2.0
    return np.concatenate([[left], inner, [right]])


def _resolve_data_path(path: Optional[str]) -> Optional[str]:
    """Resolve a data path from either the current directory or the ALBM module directory."""
    if path is None:
        return None
    path_str = os.fspath(path)
    candidates = [path_str]
    if not os.path.isabs(path_str):
        candidates.append(os.path.join(os.path.dirname(__file__), path_str))
    for candidate in candidates:
        if os.path.exists(candidate):
            return candidate
    return path_str


def _numeric_or_nan(value: Any) -> float:
    parsed = pd.to_numeric(pd.Series([value]), errors='coerce').iloc[0]
    return float(parsed) if pd.notna(parsed) else np.nan


def _processed_mtt_date(filepath: str, sheet_name: str) -> pd.Timestamp:
    """Use the filename date, except for multi-date sheets that clearly identify another nearby date."""
    filename_match = re.search(r'Processed_(\d{8})', os.path.basename(filepath))
    file_date = (
        pd.to_datetime(filename_match.group(1), format='%Y%m%d')
        if filename_match else pd.NaT
    )
    sheet_match = re.search(r'(\d{8})', str(sheet_name))
    if sheet_match:
        sheet_date = pd.to_datetime(sheet_match.group(1), format='%Y%m%d')
        if pd.isna(file_date) or abs((sheet_date - file_date).days) <= 45:
            return sheet_date
    return file_date


def _load_mtt_oxidation_observations(
    mtt_dir: str = 'data/BTL_MTT',
) -> pd.DataFrame:
    """
    Load processed Big Trail Lake MTT CH4 oxidation observations.

    Source rates are CH4 net rates in mmol m-3 d-1. Negative net CH4 rates are
    methane consumption, so the plotted oxidation rate is -net_rate converted to
    kg m-3 s-1.
    """
    root = _resolve_data_path(mtt_dir)
    if root is None or not os.path.isdir(root):
        logger.info("MTT oxidation directory not found: %s", mtt_dir)
        return pd.DataFrame()

    rows: List[Dict[str, Any]] = []
    pattern = os.path.join(root, '**', 'BTL_MTT_Processed*.xlsx')
    for filepath in sorted(glob_path for glob_path in _recursive_glob(pattern)
                           if not os.path.basename(glob_path).startswith('~$')):
        try:
            workbook = pd.ExcelFile(filepath)
        except Exception as exc:
            logger.warning("Skipping MTT workbook %s (%s)", filepath, exc)
            continue

        for sheet_name in workbook.sheet_names:
            try:
                sheet = pd.read_excel(filepath, sheet_name=sheet_name, header=None, dtype=object)
            except Exception as exc:
                logger.warning("Skipping MTT sheet %s:%s (%s)", filepath, sheet_name, exc)
                continue

            rate_cell = None
            for row_idx in range(sheet.shape[0]):
                for col_idx in range(sheet.shape[1]):
                    value = sheet.iat[row_idx, col_idx]
                    if isinstance(value, str) and 'CH4 net rate' in value:
                        rate_cell = (row_idx, col_idx)
                        break
                if rate_cell is not None:
                    break
            if rate_cell is None:
                continue

            header_row, rate_col = rate_cell
            header = [
                str(value).strip().lower() if pd.notna(value) else ''
                for value in sheet.iloc[header_row].tolist()
            ]
            depth_col = next((i for i, label in enumerate(header) if 'depth' in label), None)
            site_col = next((i for i, label in enumerate(header) if label == 'site'), None)
            if depth_col is None:
                continue

            measurement_date = _processed_mtt_date(filepath, sheet_name)
            if pd.isna(measurement_date):
                continue

            saw_data = False
            for row_idx in range(header_row + 1, sheet.shape[0]):
                row = sheet.iloc[row_idx]
                row_text = ' '.join(str(value) for value in row.tolist() if isinstance(value, str))
                if saw_data and ('CO2' in row_text or 'Sample ID' in row_text):
                    break

                depth_m = _numeric_or_nan(row.iloc[depth_col])
                net_rate = _numeric_or_nan(row.iloc[rate_col])
                if np.isfinite(depth_m) and np.isfinite(net_rate):
                    saw_data = True
                    rows.append({
                        'datetime': measurement_date,
                        'doy': int(measurement_date.dayofyear),
                        'depth_m': depth_m,
                        'site': (
                            str(row.iloc[site_col]).strip()
                            if site_col is not None and pd.notna(row.iloc[site_col])
                            else ''
                        ),
                        'ch4_net_rate_mmol_m3_d': net_rate,
                    })
                elif saw_data and row.isna().all():
                    break

    if not rows:
        return pd.DataFrame()

    obs = pd.DataFrame(rows)
    grouped = (
        obs.groupby(['datetime', 'doy', 'depth_m'], as_index=False)
        .agg(
            ch4_net_rate_mmol_m3_d=('ch4_net_rate_mmol_m3_d', 'mean'),
            n=('ch4_net_rate_mmol_m3_d', 'size'),
        )
    )
    grouped['oxidation_kg_m3_s'] = (
        -grouped['ch4_net_rate_mmol_m3_d'] * 16.04e-6 / 86400.0
    )
    return grouped


def _recursive_glob(pattern: str) -> List[str]:
    """Small wrapper so recursive file matching stays easy to test and mock."""
    import glob
    return glob.glob(pattern, recursive=True)


def _load_dissolved_ch4_observations(
    filepath: str = 'data/2022_11_21_CH4CO2_DeepShallow.xlsx',
) -> pd.DataFrame:
    """Load dissolved CH4 observations and convert CH4mgL to kg m-3."""
    resolved = _resolve_data_path(filepath)
    if resolved is None or not os.path.exists(resolved):
        logger.info("Dissolved CH4 workbook not found: %s", filepath)
        return pd.DataFrame()

    rows = []
    try:
        workbook = pd.ExcelFile(resolved)
    except Exception as exc:
        logger.warning("Could not open dissolved CH4 workbook %s (%s)", resolved, exc)
        return pd.DataFrame()

    for sheet_name in workbook.sheet_names:
        try:
            sheet = pd.read_excel(resolved, sheet_name=sheet_name)
        except Exception as exc:
            logger.warning("Skipping dissolved CH4 sheet %s:%s (%s)", resolved, sheet_name, exc)
            continue

        date_col = next((col for col in sheet.columns if str(col).strip().lower() == 'date'), None)
        depth_col = next((col for col in sheet.columns if str(col).strip().lower() == 'depth'), None)
        ch4_col = next((col for col in sheet.columns if str(col).strip().lower() == 'ch4mgl'), None)
        site_col = next((col for col in sheet.columns if str(col).strip().lower() == 'site'), None)
        if date_col is None or depth_col is None or ch4_col is None:
            continue

        obs = pd.DataFrame({
            'datetime': pd.to_datetime(sheet[date_col], errors='coerce'),
            'depth_m': pd.to_numeric(sheet[depth_col], errors='coerce'),
            'ch4_mg_l': pd.to_numeric(sheet[ch4_col], errors='coerce'),
            'site': sheet[site_col].astype(str) if site_col is not None else '',
        }).dropna(subset=['datetime', 'depth_m', 'ch4_mg_l'])
        rows.append(obs)

    if not rows:
        return pd.DataFrame()

    combined = pd.concat(rows, ignore_index=True)
    grouped = (
        combined.groupby(['datetime', 'depth_m'], as_index=False)
        .agg(ch4_mg_l=('ch4_mg_l', 'mean'), n=('ch4_mg_l', 'size'))
    )
    grouped['doy'] = grouped['datetime'].dt.dayofyear
    grouped['ch4_kg_m3'] = grouped['ch4_mg_l'] * 1e-3
    return grouped


def _overlay_observations(
    ax: plt.Axes,
    observations: pd.DataFrame,
    value_col: str,
    cmap: str,
    norm: Any,
    label: str,
    marker: str,
    size: float = 26.0,
    linewidth: float = 0.3,
    positive_only: bool = False,
) -> None:
    """Overlay depth-vs-DOY observations colored in the same units as the heatmap."""
    if observations is None or observations.empty:
        return
    valid = observations[
        np.isfinite(observations['doy'])
        & np.isfinite(observations['depth_m'])
        & np.isfinite(observations[value_col])
    ]
    if positive_only or isinstance(norm, LogNorm):
        valid = valid[valid[value_col] > 0]
    if valid.empty:
        return

    ax.scatter(
        valid['doy'],
        valid['depth_m'],
        c=valid[value_col],
        cmap=cmap,
        norm=norm,
        s=size,
        marker=marker,
        edgecolors='black',
        linewidths=linewidth,
        alpha=0.92,
        zorder=5,
        label=label,
    )
    ax.legend(
        loc='lower center',
        bbox_to_anchor=(0.5, 0.02),
        frameon=True,
        framealpha=0.85,
        fontsize=8,
    )


def plot_water_column_doy_climatology_heatmaps(
    albm_data: ALBMData,
    figs_dir: str = 'figs',
    show: bool = True,
    mtt_dir: str = 'data/BTL_MTT',
    dissolved_profile_file: str = 'data/2022_11_21_CH4CO2_DeepShallow.xlsx',
) -> plt.Figure:
    """
    Four-panel (2×2) heatmaps of the same variables as ``plot_water_column_heatmaps``:
    dissolved CH₄, oxidation, oxic production, and sediment temperature (if available)
    or CH₄ exchange. Each panel shows the multi-year mean by day-of-year (DOY 1–366);
    x-axis is day of year, y-axis is depth. Subplots use a portrait aspect ratio.
    """
    logger.info("Creating water column DOY climatology heatmaps (2×2)...")
    MW_CH4 = 16.04e-3
    has_sed_temp = (
        hasattr(albm_data, 'sediment_temp_depth')
        and albm_data.sediment_temp_depth is not None
        and albm_data.sediment_temp_depth.ndim >= 2
    )
    time_array = albm_data.time_array

    # 2×2 figure; height moderate — set_box_aspect is not used so rows can pack tightly
    fig, axes = plt.subplots(2, 2, figsize=(10, 12), constrained_layout=False)
    axes_flat = axes.flatten().tolist()
    _apply_figure_style(fig, axes_flat)

    panels = [
        {
            'data': albm_data.ch4conc * MW_CH4,
            'title': 'Dissolved CH₄ Concentration',
            'cbar_label': 'CH₄ Concentration (kg m⁻³)',
            'cmap': 'plasma',
            'y': albm_data.water_depth_array,
            'ylim': (4.1, 0),  # 0 m at top of axis (surface)
            'ylabel': 'Water Depth (m)',
            'cbar_scientific': True,
        },
        {
            'data': albm_data.ch4oxid * MW_CH4,
            'title': 'CH₄ Oxidation Rate',
            'cbar_label': 'CH₄ Oxidation (kg m⁻³ s⁻¹)',
            'cmap': 'RdPu',
            'y': albm_data.water_depth_array,
            'ylim': (4.1, 0),
            'ylabel': 'Water Depth (m)',
            'cbar_scientific': True,
            'log_scale': True,
            'log_bounds': [1e-11, 5e-11, 1e-10, 5e-10, 1e-9, 5e-9, 1e-8],
            'cbar_extend': 'both',
        },
        {
            'data': albm_data.och4prod * MW_CH4,
            'title': 'Oxic CH₄ Production',
            'cbar_label': 'Oxic CH₄ Production (kg m⁻³ s⁻¹)',
            'cmap': 'cividis',
            'y': albm_data.water_depth_array,
            'ylim': (4.1, 0),
            'ylabel': 'Water Depth (m)',
            'cbar_scientific': True,
        },
    ]

    if has_sed_temp:
        n_sed = albm_data.sediment_temp_depth.shape[1]
        sed_depths = _sediment_depth_array(n_sed)
        hsed = float(sed_depths[-1])
        panels.append({
            'data': albm_data.sediment_temp_depth,
            'title': 'Sediment Temperature by Depth',
            'cbar_label': 'Temperature (°C)',
            'cmap': 'coolwarm',
            'y': sed_depths,
            'ylim': (hsed, 0),  # 0 m sediment top at top of axis
            'ylabel': 'Sediment Depth (m)',
        })
    else:
        panels.append({
            'data': albm_data.ch4exc * MW_CH4,
            'title': 'CH₄ Exchange Rate',
            'cbar_label': 'CH₄ Exchange (kg m⁻³ s⁻¹)',
            'cmap': 'viridis',
            'y': albm_data.water_depth_array,
            'ylim': (4.1, 0),
            'ylabel': 'Water Depth (m)',
        })

    dissolved_obs = _load_dissolved_ch4_observations(dissolved_profile_file)
    oxidation_obs = _load_mtt_oxidation_observations(mtt_dir)
    if not dissolved_obs.empty:
        logger.info("Loaded %d averaged dissolved CH4 profile observations.", len(dissolved_obs))
    if not oxidation_obs.empty:
        logger.info("Loaded %d averaged MTT oxidation observations.", len(oxidation_obs))

    # 366 DOY bins: edges 0.5 … 366.5 (367 edges)
    doy_edges = np.linspace(0.5, 366.5, 367)

    # Pass 1: heatmaps only (colorbars after layout so bar height matches axes bbox)
    pending_cbars = []

    for i, (ax, plot_info) in enumerate(zip(axes_flat, panels)):
        raw = plot_info['data']
        if raw.ndim != 2:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            continue
        data_tx = np.asarray(raw)
        if data_tx.shape[0] != len(time_array):
            # Align on minimum length if needed
            n = min(data_tx.shape[0], len(time_array))
            data_tx = data_tx[:n, :]
            t_use = time_array[:n]
        else:
            t_use = time_array

        mean_doy = _doy_mean_climatology(t_use, data_tx)
        plot_data = mean_doy
        y_centers = np.asarray(plot_info['y'], dtype=float)
        y_edges = _centers_to_edges_1d(y_centers)

        norm = None
        cmap = plot_info['cmap']
        if plot_info.get('log_scale', False):
            positive = plot_data[np.isfinite(plot_data) & (plot_data > 0)]
            if positive.size > 0:
                cmap = plt.get_cmap(cmap).copy()
                log_bounds = plot_info.get('log_bounds')
                if log_bounds is not None:
                    log_bounds = np.asarray(log_bounds, dtype=float)
                    norm = BoundaryNorm(
                        log_bounds,
                        ncolors=cmap.N,
                        extend=plot_info.get('cbar_extend', 'neither'),
                    )
                else:
                    vmin = float(plot_info.get('log_vmin', np.nanmin(positive)))
                    vmax = float(plot_info.get('log_vmax', np.nanmax(positive)))
                    if vmax <= vmin:
                        vmax = vmin * 10.0
                    norm = LogNorm(vmin=vmin, vmax=vmax)
                plot_data = np.ma.masked_where(plot_data <= 0, plot_data)
                cmap.set_bad('white')

        im = ax.pcolormesh(
            doy_edges,
            y_edges,
            plot_data,
            shading='flat',
            cmap=cmap,
            norm=norm,
        )
        # ylim (deeper, 0): surface depth 0 m at top of axis (no invert_yaxis needed)
        ax.set_ylim(plot_info['ylim'])
        # Right column (b,d): extra label padding so y-axis titles sit clear of left colorbars
        ax.set_ylabel(
            plot_info['ylabel'],
            fontsize=12,
            labelpad=12 if i in (1, 3) else 6,
        )
        ax.set_title(plot_info['title'], fontsize=14)
        ax.set_xlim(0.5, 366.5)  # DOY 1–366
        # Day of year on every panel; tighter labelpad on top row to limit vertical gap
        ax.set_xlabel(
            'Day of year',
            fontsize=12,
            labelpad=3 if i in (0, 1) else 5,
        )
        ax.text(
            -0.12,
            1.02,
            _subplot_label(i),
            transform=ax.transAxes,
            fontsize=14,
            fontweight='bold',
            verticalalignment='bottom',
            horizontalalignment='left',
        )
        # No set_box_aspect here: it shrinks the axes inside each grid cell and leaves
        # large empty bands between rows; heatmaps fill cells for tighter vertical packing.
        setup_axis_style(ax)
        if i == 0:
            _overlay_observations(
                ax=ax,
                observations=dissolved_obs,
                value_col='ch4_kg_m3',
                cmap=plot_info['cmap'],
                norm=im.norm,
                label='Observed dissolved CH4',
                marker='o',
                size=12,
                linewidth=0.18,
            )
        elif i == 1:
            _overlay_observations(
                ax=ax,
                observations=oxidation_obs,
                value_col='oxidation_kg_m3_s',
                cmap=plot_info['cmap'],
                norm=im.norm,
                label='Observed MTT oxidation',
                marker='^',
                size=32,
                linewidth=0.35,
                positive_only=True,
            )
        pending_cbars.append((i, ax, im, plot_info))

    # w_pad: horizontal gap for colorbars vs right-column y-labels.
    # h_pad + positive hspace: comfortable vertical gap between top and bottom rows
    plt.tight_layout(pad=0.5, w_pad=6.5, h_pad=1.5)
    fig.subplots_adjust(hspace=0.2)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    # Pass 2: skinny colorbars exactly the height of each heatmap (drawn bbox)
    for i, ax, im, plot_info in pending_cbars:
        bbox = ax.get_window_extent(renderer=renderer).transformed(
            fig.transFigure.inverted()
        )
        width_frac = 0.012
        # Left column (a,c): smaller pad keeps the colorbar closer to its heatmap so the
        # bar + tick labels intrude less on the right panel's y-axis label (w_pad also helps)
        pad_frac = 0.005 if i in (0, 2) else 0.008
        cax = fig.add_axes(
            [bbox.x1 + pad_frac, bbox.y0, width_frac, bbox.height],
        )
        cbar = plt.colorbar(
            im,
            cax=cax,
            label=plot_info['cbar_label'],
            extend=plot_info.get('cbar_extend', 'neither'),
        )
        cbar.ax.tick_params(labelsize=9)
        if isinstance(im.norm, BoundaryNorm) and plot_info.get('log_bounds') is not None:
            ticks = np.asarray(plot_info['log_bounds'], dtype=float)
            cbar.set_ticks(ticks)
            cbar.ax.set_yticklabels([f'{tick:.0e}' for tick in ticks])
        elif isinstance(im.norm, LogNorm):
            cbar.ax.yaxis.set_major_formatter(LogFormatterSciNotation())
        elif plot_info.get('cbar_scientific', False):
            sci_fmt = ScalarFormatter(useMathText=True)
            sci_fmt.set_scientific(True)
            sci_fmt.set_powerlimits((0, 0))
            cbar.ax.yaxis.set_major_formatter(sci_fmt)
    os.makedirs(figs_dir, exist_ok=True)
    out_path = os.path.join(figs_dir, 'isoALBM_water_column_doy_climatology_heatmaps.png')
    plt.savefig(out_path, transparent=True, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    return fig


def plot_oxic_and_exchange_heatmaps(
    albm_data: ALBMData,
    figs_dir: str = 'figs',
    show: bool = True
) -> plt.Figure:
    """
    Plot Oxic CH₄ production and CH₄ exchange rate heatmaps (depth vs time) in a single
    figure with two subplots. Both use colorblind-friendly colormap and log scale.
    """
    logger.info("Creating oxic production and exchange rate heatmaps...")
    MW_CH4 = 16.04e-3
    time_numeric = mdates.date2num(albm_data.time_array)
    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharex=True)
    _apply_figure_style(fig, axes)
    plots = [
        {'data': albm_data.och4prod * MW_CH4, 'title': 'Oxic CH₄ Production',
         'cbar_label': 'Oxic CH₄ Production (kg m⁻³ s⁻¹)', 'cmap': 'cividis'},
        {'data': np.abs(albm_data.ch4exc * MW_CH4), 'title': 'CH₄ Exchange Rate',
         'cbar_label': 'CH₄ Exchange (kg m⁻³ s⁻¹)', 'cmap': 'cividis'},
    ]
    for i, (ax, plot_info) in enumerate(zip(axes, plots)):
        plot_data = plot_info['data'].T
        positive = plot_data[plot_data > 0]
        small = (np.nanmin(positive) * 1e-6) if len(positive) > 0 else 1e-20
        plot_data_log = np.where(plot_data > 0, plot_data, small)
        vmin = max(float(np.nanmin(plot_data_log)), 1e-20)
        vmax = max(float(np.nanmax(plot_data_log)), vmin * 10)
        im = ax.pcolormesh(time_numeric, albm_data.water_depth_array, plot_data_log,
                          shading='auto', cmap=plot_info['cmap'], norm=LogNorm(vmin=vmin, vmax=vmax))
        ax.invert_yaxis()
        ax.set_ylim(4.1, 0)
        ax.xaxis_date()
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=3))
        ax.set_ylabel('Depth (m)', fontsize=14)
        ax.set_title(plot_info['title'], fontsize=16)
        ax.text(-0.08, 1.02, _subplot_label(i), transform=ax.transAxes, fontsize=16, fontweight='bold',
                verticalalignment='bottom', horizontalalignment='left')
        cbar = plt.colorbar(im, ax=ax, label=plot_info['cbar_label'])
        cbar.ax.tick_params(labelsize=10)
        setup_axis_style(ax)
    axes[-1].set_xlabel('Year', fontsize=14)
    plt.tight_layout()
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, 'isoALBM_oxic_and_exchange_heatmaps.png'),
                transparent=True, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    return fig


def plot_eddy_comparison(
    albm_data: ALBMData,
    eddy_data: Optional[EddyFluxData],
    figs_dir: str = 'figs',
    show: bool = True
) -> plt.Figure:
    """
    Plot comparison between ALBM and eddy covariance fluxes.
    
    Parameters:
    -----------
    albm_data : ALBMData
        ALBM data container
    eddy_data : EddyFluxData
        Eddy covariance data container
    figs_dir : str
        Directory to save figures
    show : bool
        Whether to display the plot
    
    Returns
    -------
    matplotlib.figure.Figure
        The figure (or a minimal no-data figure if no eddy data).
    """
    if eddy_data is None:
        logger.info("No eddy flux data available for comparison plot")
        return _no_data_figure("No eddy flux data")
    logger.info("Creating eddy covariance comparison plot...")
    eddy_df = eddy_data.ch4_flux_daily_avg
    eddy_start = eddy_df['datetime'].min()
    eddy_end = eddy_df['datetime'].max()
    model_start = albm_data.time_array[0]
    model_end = albm_data.time_array[-1]
    overlap_start = max(eddy_start, model_start)
    overlap_end = min(eddy_end, model_end)
    eddy_mask = (eddy_df['datetime'] >= overlap_start) & (eddy_df['datetime'] <= overlap_end)
    model_mask = (albm_data.time_array >= overlap_start) & (albm_data.time_array <= overlap_end)
    eddy_dates = eddy_df.loc[eddy_mask, 'datetime']
    n_model = np.sum(model_mask)
    n_eddy = np.sum(eddy_mask)
    if n_model < 2 or n_eddy == 0:
        model_interp = np.full(n_eddy, np.nan) if n_eddy else np.array([])
        eddy_values = eddy_df.loc[eddy_mask, 'ch4_flux_daily_avg_kg_m2_s'].values if n_eddy else np.array([])
    else:
        model_interp = np.interp(
            mdates.date2num(eddy_dates),
            mdates.date2num(albm_data.time_array[model_mask]),
            albm_data.total_flux[model_mask]
        )
        eddy_values = eddy_df.loc[eddy_mask, 'ch4_flux_daily_avg_kg_m2_s'].values
    valid = np.isfinite(model_interp) & np.isfinite(eddy_values) if len(model_interp) else np.array([], dtype=bool)
    rmse = np.sqrt(np.mean((model_interp[valid] - eddy_values[valid]) ** 2)) if np.sum(valid) > 0 else np.nan
    if np.isfinite(rmse):
        logger.info("RMSE (isoALBM Surface Flux vs Eddy Covariance): %s kg m^-2 s^-1", f"{rmse:.2e}")
    else:
        logger.info("RMSE (isoALBM Surface Flux vs Eddy Covariance): N/A (insufficient overlap or invalid data)")
    
    # Create figure with 3 subplots: fluxes, anomaly, flux components (log)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(14, 11))
    _apply_figure_style(fig, [ax1, ax2, ax3])

    # Subplot a: Flux comparison
    ax1.plot(eddy_df.loc[eddy_mask, 'datetime'],
             eddy_df.loc[eddy_mask, 'ch4_flux_daily_avg_kg_m2_s'],
             label='Eddy Covariance Flux', alpha=0.7, color=CB_BLACK)
    ax1.plot(albm_data.time_array[model_mask], albm_data.total_flux[model_mask],
             label='isoALBM Surface Flux', alpha=0.7, color=CB_ORANGE)
    ax1.set_ylabel('CH₄ Flux (kg m⁻² s⁻¹)', fontsize=16)
    ax1.grid(True)
    ax1.legend(loc=LEGEND_UPPER_RIGHT, framealpha=1, edgecolor='black')
    setup_axis_style(ax1)
    ax1.text(-0.08, 1.02, _subplot_label(0), transform=ax1.transAxes, fontsize=16, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='left')
    
    # Subplot b: Flux anomaly
    model_difference = model_interp - eddy_values
    ax2.scatter(eddy_dates, model_difference, alpha=0.5, color=CB_ORANGE, s=10, label='isoALBM - Eddy')
    ax2.set_ylabel('CH₄ Flux Anomaly\n(kg m⁻² s⁻¹)', fontsize=16)
    ax2.grid(True)
    ax2.legend(loc=LEGEND_UPPER_RIGHT, framealpha=1, edgecolor='black')
    setup_axis_style(ax2)
    ax2.text(-0.08, 1.02, _subplot_label(1), transform=ax2.transAxes, fontsize=16, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='left')
    # Subplot c: Flux components (log scale) over overlap period
    t_overlap = albm_data.time_array[model_mask]
    ax3.plot(t_overlap, albm_data.dflux_data[model_mask], color=CB_GREEN,
             label='Diffusion CH₄ Flux', linewidth=2)
    ax3.plot(t_overlap, albm_data.eflux_data[model_mask] + albm_data.ch4upb_data[model_mask],
             color=CB_ORANGE, label='Ebullition CH₄ Flux', linewidth=2)
    ax3.set_yscale('log')
    ax3.set_ylabel('CH₄ Flux (kg m⁻² s⁻¹)', fontsize=16)
    ax3.set_xlabel('Year', fontsize=16)
    ax3.grid(True, which='both', linestyle='--', alpha=0.7)
    ax3.legend(loc=LEGEND_UPPER_RIGHT, framealpha=1, edgecolor='black')
    setup_axis_style(ax3)
    ax3.text(-0.08, 1.02, _subplot_label(2), transform=ax3.transAxes, fontsize=16, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='left')
    
    # Format x-axis for all
    for ax in [ax1, ax2, ax3]:
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter(r'$\mathbf{%Y}$'))
        ax.xaxis.set_minor_locator(mdates.MonthLocator())
        plt.setp(ax.xaxis.get_majorticklabels(), fontsize=10)
        plt.setp(ax.xaxis.get_minorticklabels(), rotation=0)
    
    plt.tight_layout()
    
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, 'ALBM_vs_eddy_comparison.png'), 
                dpi=300, transparent=True)
    
    if show:
        plt.show()
    
    return fig


def plot_flux_components(
    albm_data: ALBMData,
    figs_dir: str = 'figs',
    show: bool = True
) -> plt.Figure:
    """
    Plot diffusion and ebullition flux components.
    
    Parameters:
    -----------
    albm_data : ALBMData
        ALBM data container
    figs_dir : str
        Directory to save figures
    show : bool
        Whether to display the plot
    
    Returns:
    --------
    fig : matplotlib Figure
    """
    logger.info("Creating flux components plot (log scale)...")
    fig, ax = plt.subplots(figsize=(12, 6))
    _apply_figure_style(fig, ax)

    # Plot diffusion and ebullition fluxes
    ax.plot(albm_data.time_array[:len(albm_data.dflux_data)], albm_data.dflux_data, 
            color=CB_GREEN, label='Diffusion CH₄ Flux', linewidth=2)
    ax.plot(albm_data.time_array[:len(albm_data.eflux_data)], 
            albm_data.eflux_data + albm_data.ch4upb_data, 
            color=CB_ORANGE, label='Ebullition CH₄ Flux', linewidth=2)
    
    ax.set_ylabel('CH₄ Flux (kg m⁻² s⁻¹)', fontsize=16)
    ax.set_xlabel('Year', fontsize=16)
    ax.set_yscale('log')
    ax.grid(True, which='both', linestyle='--', alpha=0.7)
    ax.legend(fontsize=14, loc=LEGEND_LOWER_RIGHT, framealpha=1, edgecolor='black')
    ax.tick_params(axis='both', which='major', labelsize=16, width=2)
    setup_axis_style(ax)
    ax.set_xlim(albm_data.start_date, albm_data.end_date)
    ax.xaxis.set_major_locator(mdates.YearLocator(1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    
    plt.tight_layout()
    
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, 'ALBM_separate_diffusion_ebullition_fluxes.png'), 
                transparent=True, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    
    return fig


def plot_flux_components_linear(
    albm_data: ALBMData,
    figs_dir: str = 'figs',
    show: bool = True
) -> plt.Figure:
    """
    Plot diffusion and ebullition flux components (linear scale).
    
    Parameters:
    -----------
    albm_data : ALBMData
        ALBM data container
    figs_dir : str
        Directory to save figures
    show : bool
        Whether to display the plot
    
    Returns:
    --------
    fig : matplotlib Figure
    """
    logger.info("Creating flux components plot (linear scale)...")
    fig, ax = plt.subplots(figsize=(12, 6))
    _apply_figure_style(fig, ax)

    # Plot diffusion and ebullition fluxes
    ax.plot(albm_data.time_array[:len(albm_data.dflux_data)], albm_data.dflux_data, 
            color=CB_GREEN, label='Diffusion CH₄ Flux', linewidth=2)
    ax.plot(albm_data.time_array[:len(albm_data.eflux_data)], 
            albm_data.eflux_data + albm_data.ch4upb_data, 
            color=CB_ORANGE, label='Ebullition CH₄ Flux', linewidth=2)
    
    ax.set_ylabel('CH₄ Flux (kg m⁻² s⁻¹)', fontsize=16)
    ax.set_xlabel('Year', fontsize=16)
    ax.set_title('isoALBM CH₄ Diffusion and Ebullition Flux', fontsize=18)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend(fontsize=14, loc=LEGEND_UPPER_RIGHT, framealpha=1, edgecolor='black')
    ax.tick_params(axis='both', which='major', labelsize=16, width=2)
    setup_axis_style(ax)
    ax.set_xlim(albm_data.start_date, albm_data.end_date)
    ax.xaxis.set_major_locator(mdates.YearLocator(1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    
    plt.tight_layout()
    
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, 'ALBM_diffusion_ebullition_linear.png'), 
                transparent=True, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    
    return fig


def plot_sediment_temperature_heatmap(
    albm_data: ALBMData,
    figs_dir: str = 'figs',
    show: bool = True
) -> plt.Figure:
    """
    Plot sediment temperature heatmap (depth vs time). The single time series
    used in the isotope model is the CH4-flux-weighted mean over depth.

    Parameters
    ----------
    albm_data : ALBMData
        ALBM data container (must have sediment_temp_depth).
    figs_dir : str
        Directory to save figures.
    show : bool
        Whether to display the plot.

    Returns
    -------
    matplotlib.figure.Figure
        The figure (or a minimal no-data figure if sediment data is missing).
    """
    if not hasattr(albm_data, 'sediment_temp_depth') or albm_data.sediment_temp_depth is None:
        logger.warning("No sediment temperature data available")
        return _no_data_figure("No sediment temperature data")
    fig, ax = plt.subplots(figsize=(14, 6))
    _apply_figure_style(fig, ax)
    n_sed_depths = albm_data.sediment_temp_depth.shape[1] if albm_data.sediment_temp_depth.ndim > 1 else 1
    sed_depths = _sediment_depth_array(n_sed_depths)
    hsed = float(sed_depths[-1])
    time_numeric = mdates.date2num(albm_data.time_array)

    im = ax.pcolormesh(time_numeric, sed_depths, albm_data.sediment_temp_depth.T,
                       shading='nearest', cmap='coolwarm')
    ax.set_ylim(0, hsed)
    ax.invert_yaxis()
    ax.xaxis_date()
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=3))
    ax.set_xlabel('Year', fontsize=16)
    ax.set_ylabel('Sediment Depth (m)', fontsize=16)
    ax.set_title('Sediment Temperature by Depth', fontsize=18)
    cbar = plt.colorbar(im, ax=ax, label='Temperature (°C)')
    cbar.ax.tick_params(labelsize=12)
    setup_axis_style(ax)
    plt.tight_layout()
    
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, 'isoALBM_sediment_temperature_heatmap.png'),
                transparent=True, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    return fig


def plot_sediment_temperature_avg(
    albm_data: ALBMData,
    figs_dir: str = 'figs',
    show: bool = True
) -> plt.Figure:
    """
    Plot flux-weighted sediment temperature over time (mean weighted by CH4
    diffusion + ebullition per depth, i.e. temperature experienced by methane).

    Parameters
    ----------
    albm_data : ALBMData
        ALBM data container (must have sediment_temp_depth or sediment_temp_avg).
    figs_dir : str
        Directory to save figures.
    show : bool
        Whether to display the plot.

    Returns
    -------
    matplotlib.figure.Figure
        The figure (or a minimal no-data figure if sediment data is missing).
    """
    if not hasattr(albm_data, 'sediment_temp_depth') or albm_data.sediment_temp_depth is None:
        logger.warning("No sediment temperature data available")
        return _no_data_figure("No sediment temperature data")
    if hasattr(albm_data, 'sediment_temp_avg') and albm_data.sediment_temp_avg is not None:
        sediment_temp_avg = albm_data.sediment_temp_avg
    else:
        sediment_temp_avg = np.mean(albm_data.sediment_temp_depth, axis=1)
    
    fig, ax = plt.subplots(figsize=(14, 5))
    _apply_figure_style(fig, ax)
    time_numeric = mdates.date2num(albm_data.time_array)
    ax.plot(time_numeric, sediment_temp_avg, color=CB_BLUE, linewidth=2)
    ax.xaxis_date()
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=3))
    ax.set_xlabel('Year', fontsize=16)
    ax.set_ylabel('Temperature (°C)', fontsize=16)
    ax.set_title('Flux-Weighted Sediment Temperature', fontsize=18)
    ax.grid(True, alpha=0.3)
    setup_axis_style(ax)
    plt.tight_layout()
    
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, 'isoALBM_sediment_temperature_avg.png'),
                transparent=True, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    return fig


def plot_sediment_temperature_by_depth(
    albm_data: ALBMData,
    figs_dir: str = 'figs',
    show: bool = True
) -> plt.Figure:
    """
    Line plot of sediment temperature over time with one line per depth.
    The flux-weighted mean of these depths is used in the isotope model.

    Parameters
    ----------
    albm_data : ALBMData
        ALBM data container (must have sediment_temp_depth).
    figs_dir : str
        Directory to save figures.
    show : bool
        Whether to display the plot.

    Returns
    -------
    matplotlib.figure.Figure
        The figure (or a minimal no-data figure if sediment data is missing).
    """
    if not hasattr(albm_data, 'sediment_temp_depth') or albm_data.sediment_temp_depth is None:
        logger.warning("No sediment temperature data available")
        return _no_data_figure("No sediment temperature data")
    n_sed_depths = albm_data.sediment_temp_depth.shape[1] if albm_data.sediment_temp_depth.ndim > 1 else 1
    sed_depths = _sediment_depth_array(n_sed_depths)
    time_numeric = mdates.date2num(albm_data.time_array)
    max_legend_lines = 10
    if n_sed_depths > max_legend_lines:
        label_indices = np.round(np.linspace(0, n_sed_depths - 1, max_legend_lines)).astype(int)
    else:
        label_indices = np.arange(n_sed_depths)
    _palette = [CB_VERMILLION, CB_ORANGE, CB_GREEN, CB_SKY_BLUE, CB_BLUE, CB_PURPLE]
    fig, ax = plt.subplots(figsize=(14, 5))
    _apply_figure_style(fig, ax)

    cmap_depth = plt.cm.coolwarm_r
    for j in range(n_sed_depths):
        depth_m = sed_depths[j]
        color = cmap_depth(j / max(n_sed_depths - 1, 1))
        label = f'{depth_m:.2f} m' if j in label_indices else None
        ax.plot(
            time_numeric,
            albm_data.sediment_temp_depth[:, j],
            color=color,
            linewidth=0.8 if n_sed_depths > 10 else 1.2,
            label=label
        )

    ax.xaxis_date()
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=3))
    ax.set_xlabel('Year', fontsize=16)
    ax.set_ylabel('Temperature (°C)', fontsize=16)
    ax.set_title('Sediment Temperature by Depth', fontsize=18)
    ax.legend(loc='upper right', fontsize=10, ncol=2 if n_sed_depths > 3 else 1)
    ax.grid(True, alpha=0.3)
    setup_axis_style(ax)
    plt.tight_layout()

    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, 'isoALBM_sediment_temperature_by_depth.png'),
                transparent=True, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    return fig


def plot_sediment_temperature(
    albm_data: ALBMData,
    figs_dir: str = 'figs',
    show: bool = True
) -> plt.Figure:
    """
    Plot sediment temperature: heatmap (depth vs time), flux-weighted mean time series,
    and line plot by depth as separate figures.

    Parameters
    ----------
    albm_data : ALBMData
        ALBM data container.
    figs_dir : str
        Directory to save figures.
    show : bool
        Whether to display the plots.

    Returns
    -------
    matplotlib.figure.Figure
        The flux-weighted mean time series figure (or no-data figure if no sediment data).
    """
    logger.info("Creating sediment temperature plots (heatmap, flux-weighted mean, by-depth)...")
    fig_heat = plot_sediment_temperature_heatmap(albm_data, figs_dir=figs_dir, show=show)
    fig_avg = plot_sediment_temperature_avg(albm_data, figs_dir=figs_dir, show=show)
    fig_lines = plot_sediment_temperature_by_depth(albm_data, figs_dir=figs_dir, show=show)
    return fig_avg


def plot_ch4_budget(
    albm_data: ALBMData,
    figs_dir: str = 'figs',
    show: bool = True
) -> plt.Figure:
    """
    Plot CH4 budget components time series.
    
    Parameters:
    -----------
    albm_data : ALBMData
        ALBM data container
    figs_dir : str
        Directory to save figures
    show : bool
        Whether to display the plot
    
    Returns:
    --------
    fig : matplotlib Figure
    """
    logger.info("Creating CH4 budget components plot...")
    fig, ax = plt.subplots(figsize=(12, 6))
    _apply_figure_style(fig, ax)
    # Plot budget components
    ax.plot(albm_data.time_array[:len(albm_data.sedch4df_data)], albm_data.sedch4df_data, 
            color=CB_VERMILLION, label='Sediment diffusion', linewidth=2)
    ax.plot(albm_data.time_array[:len(albm_data.och4prod_data)], albm_data.och4prod_data, 
            color=CB_GREEN, label='Oxic production', linewidth=2)
    ax.plot(albm_data.time_array[:len(albm_data.ch4exc_data)], albm_data.ch4exc_data, 
            color=CB_SKY_BLUE, label='Exchange', linewidth=2)
    ax.plot(albm_data.time_array[:len(albm_data.och4_data)], -albm_data.och4_data, 
            color=CB_ORANGE, label='-Oxidation', linewidth=2)
    ax.plot(albm_data.time_array[:len(albm_data.dflux_data)], -albm_data.dflux_data, 
            color=CB_PURPLE, label='-Diffusion flux', linewidth=2)
    
    ax.set_ylabel('CH₄ Flux (kg m⁻³ s⁻¹)', fontsize=16)
    ax.set_xlabel('Year', fontsize=16)
    ax.set_title('ALBM CH₄ Budget Components (2021-2024)', fontsize=18)
    ax.grid(True)
    ax.legend(fontsize=12, loc='best')
    ax.tick_params(axis='both', which='major', labelsize=14)
    setup_axis_style(ax)
    ax.set_xlim(albm_data.start_date, albm_data.end_date)
    ax.xaxis.set_major_locator(mdates.YearLocator(1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    
    plt.tight_layout()
    
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, 'ALBM_ch4_budget_components.png'), 
                transparent=True, dpi=300, bbox_inches='tight')
    
    if show:
        plt.show()
    
    return fig


def plot_pool_sources_sinks(
    albm_data: ALBMData,
    figs_dir: str = 'figs',
    show: bool = True
) -> plt.Figure:
    """
    Two subplots: sources and sinks for the dissolved CH₄ pool (top) and for the
    bubble pool (bottom), as used in the isotope mass balance. All fluxes in
    mass rate (kg m⁻³ s⁻¹ or kg m⁻² s⁻¹ as provided by ALBM).
    """
    time_arr = albm_data.time_array
    n = min(
        len(time_arr),
        len(albm_data.sedch4df_data),
        len(albm_data.sedch4eb_data),
        len(albm_data.och4prod_data),
        len(albm_data.ch4exc_data),
        len(albm_data.och4_data),
        len(albm_data.dflux_data),
        len(albm_data.eflux_data),
        len(albm_data.ch4upb_data),
    )
    t = time_arr[:n]
    sedch4df = albm_data.sedch4df_data[:n]
    sedch4eb = albm_data.sedch4eb_data[:n]
    och4prod = albm_data.och4prod_data[:n]
    ch4exc = albm_data.ch4exc_data[:n]
    och4 = albm_data.och4_data[:n]
    dflux = albm_data.dflux_data[:n]
    eflux = albm_data.eflux_data[:n]
    ch4upb = albm_data.ch4upb_data[:n]

    exc_bub_to_diss = np.maximum(ch4exc, 0.0)
    exc_diss_to_bub = np.maximum(-ch4exc, 0.0)

    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharex=True)
    _apply_figure_style(fig, axes)
    # --- Dissolved pool: sources and sinks ---
    ax1 = axes[0]
    ax1.plot(t, sedch4df, color=CB_BLUE, linewidth=1.5, label='Sediment diffusion (source)')
    ax1.plot(t, och4prod, color=CB_GREEN, linewidth=1.5, label='Oxic production (source)')
    ax1.plot(t, exc_bub_to_diss, color=CB_PURPLE, linewidth=1.5, label='Bubble→dissolved (source)')
    ax1.plot(t, -och4, color=CB_ORANGE, linewidth=1.5, label='Oxidation (sink)')
    ax1.plot(t, -dflux, color=CB_VERMILLION, linewidth=1.5, label='Surface diffusion (sink)')
    ax1.plot(t, -exc_diss_to_bub, color=CB_SKY_BLUE, linewidth=1.5, label='Dissolved→bubble (sink)')
    ax1.axhline(0, color=CB_BLACK, linewidth=0.8, linestyle='-', alpha=0.5)
    ax1.set_ylabel('Flux (kg m⁻³ s⁻¹)', fontsize=14)
    ax1.set_title('Dissolved CH₄ pool: sources and sinks', fontsize=16)
    ax1.legend(loc=LEGEND_LOWER_RIGHT, fontsize=9, ncol=2)
    ax1.grid(True, alpha=0.3)
    setup_axis_style(ax1)
    ax1.text(-0.08, 1.02, _subplot_label(0), transform=ax1.transAxes, fontsize=16, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='left')
    ax1.set_xlim(t[0], t[-1])
    ax1.xaxis.set_major_locator(mdates.YearLocator(1))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

    # --- Bubble pool: sources and sinks ---
    ax2 = axes[1]
    ax2.plot(t, sedch4eb, color=CB_BLUE, linewidth=1.5, label='Sediment ebullition (source)')
    ax2.plot(t, exc_diss_to_bub, color=CB_SKY_BLUE, linewidth=1.5, label='Dissolved→bubble (source)')
    ax2.plot(t, -exc_bub_to_diss, color=CB_PURPLE, linewidth=1.5, label='Bubble→dissolved (sink)')
    ax2.plot(t, -eflux, color=CB_VERMILLION, linewidth=1.5, label='Surface ebullition (sink)')
    ax2.plot(t, -ch4upb, color=CB_ORANGE, linewidth=1.5, label='Upward bubbling (sink)')
    ax2.axhline(0, color=CB_BLACK, linewidth=0.8, linestyle='-', alpha=0.5)
    ax2.set_ylabel('Flux (kg m⁻³ s⁻¹ or kg m⁻² s⁻¹)', fontsize=14)
    ax2.set_xlabel('Year', fontsize=14)
    ax2.set_title('Bubble pool: sources and sinks', fontsize=16)
    ax2.legend(loc=LEGEND_LOWER_RIGHT, fontsize=9, ncol=2)
    ax2.text(-0.08, 1.02, _subplot_label(1), transform=ax2.transAxes, fontsize=16, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='left')
    ax2.grid(True, alpha=0.3)
    setup_axis_style(ax2)
    ax2.set_xlim(t[0], t[-1])
    ax2.xaxis.set_major_locator(mdates.YearLocator(1))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

    plt.tight_layout()
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, 'ALBM_pool_sources_sinks.png'),
                transparent=True, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    return fig


# =============================================================================
# Isotope plots
# =============================================================================

def plot_isotope_results(
    time_array: np.ndarray,
    delta_dissolved: np.ndarray,
    delta_emission: np.ndarray,
    C13_sed_prod: float,
    obs_targets: Dict[str, Dict[str, float]],
    params: Dict[str, float],
    figs_dir: str = 'figs',
    filename: str = 'ALBM_optimized_isotopes.png',
    show: bool = True
) -> plt.Figure:
    """
    Plot optimized isotope results with observations (single optimization).
    
    Parameters:
    -----------
    time_array : array
        Time array for x-axis
    delta_dissolved : array
        Dissolved CH4 δ13C values
    delta_emission : array
        Surface emission δ13C values
    C13_sed_prod : float
        Sediment production δ13C
    obs_targets : dict
        Observational targets (dissolved_summer, dissolved_winter, etc.)
    params : dict
        Model parameters used (alpha_am, alpha_hm, etc.)
    figs_dir : str
        Output directory
    filename : str
        Output filename
    show : bool
        Whether to display plot
    
    Returns:
    --------
    fig : matplotlib Figure
    """
    plot_start = time.time()
    logger.info("Generating isotope results plot: %s", filename)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    _apply_figure_style(fig, [ax1, ax2])
    time_pd = pd.DatetimeIndex(time_array)
    start_date = time_array[0]
    end_date = time_array[-1]
    max_year = time_pd.year.max()
    # --- Subplot 1: Dissolved CH4 ---
    ax1.plot(time_array, delta_dissolved, color=CB_ORANGE, linewidth=2,
             label=f'Model (α_am={params["alpha_am"]:.3f}, α_hm={params["alpha_hm"]:.3f})')
    for year in range(time_pd.year.min(), max_year + 1):
        for start_ts, end_ts, season_key in _observation_season_ranges(year, max_year):
            obs_key = 'dissolved_summer' if season_key == 'summer' else 'dissolved_winter'
            obs = obs_targets[obs_key]
            color = OBS_SUMMER_COLOR if season_key == 'summer' else OBS_WINTER_COLOR
            ax1.hlines(y=obs['mean'], xmin=start_ts, xmax=end_ts, color=color, linestyle='-', linewidth=2)
            ax1.fill_between([start_ts, end_ts], obs['mean'] - obs['std'], obs['mean'] + obs['std'], color=color, alpha=0.2)
    ax1.axhline(y=C13_sed_prod, color=CB_GREEN, linestyle='--', linewidth=2,
               label=f'Sediment Production δ¹³C = {C13_sed_prod:.1f}‰')
    
    ax1.set_ylabel('δ¹³C-CH₄ (‰)', fontsize=16)
    ax1.set_title('Optimized Dissolved CH₄ δ¹³C', fontsize=18)
    ax1.set_ylim(ISOTOPE_Y_MIN, ISOTOPE_Y_MAX)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.tick_params(axis='both', which='major', labelsize=14, width=2)
    setup_axis_style(ax1)
    ax1.set_xlim(start_date, end_date)
    ax1.xaxis.set_major_locator(mdates.YearLocator(1))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax1.text(-0.08, 1.02, _subplot_label(0), transform=ax1.transAxes, fontsize=16, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='left')
    legend_elements = [
        Line2D([0], [0], color=CB_ORANGE, linewidth=2, label='Optimized Model'),
        Line2D([0], [0], color=CB_GREEN, linestyle='--', linewidth=2,
               label=f'Sediment Prod. ({C13_sed_prod:.1f}‰)'),
        Line2D([0], [0], color=OBS_SUMMER_COLOR, linewidth=2, label='Obs. Summer Mean'),
        Patch(facecolor=OBS_SUMMER_COLOR, alpha=0.2, label='Obs. Summer ±1σ'),
        Line2D([0], [0], color=OBS_WINTER_COLOR, linewidth=2, label='Obs. Winter Mean'),
        Patch(facecolor=OBS_WINTER_COLOR, alpha=0.2, label='Obs. Winter ±1σ')
    ]
    ax1.legend(handles=legend_elements, fontsize=10, framealpha=0, loc='lower left')
    # --- Subplot 2: Surface Emissions ---
    ax2.plot(time_array, delta_emission, color=CB_ORANGE, linewidth=2,
             label=f'Model (f_am={params["f_am"]:.3f}, C13_POM={params["C13_POM"]:.1f}‰)')
    for year in range(time_pd.year.min(), max_year + 1):
        for start_ts, end_ts, season_key in _observation_season_ranges(year, max_year):
            obs_key = 'emissions_summer' if season_key == 'summer' else 'emissions_winter'
            obs = obs_targets[obs_key]
            color = OBS_SUMMER_COLOR if season_key == 'summer' else OBS_WINTER_COLOR
            ax2.hlines(y=obs['mean'], xmin=start_ts, xmax=end_ts, color=color, linestyle='-', linewidth=2)
            ax2.fill_between([start_ts, end_ts], obs['mean'] - obs['std'], obs['mean'] + obs['std'], color=color, alpha=0.2)
    ax2.axhline(y=C13_sed_prod, color=CB_GREEN, linestyle='--', linewidth=2,
               label=f'Sediment Production δ¹³C = {C13_sed_prod:.1f}‰')
    ax2.set_ylabel('δ¹³C-CH₄ (‰)', fontsize=16)
    ax2.set_xlabel('Year', fontsize=16)
    ax2.set_title('Optimized Surface Emissions δ¹³C', fontsize=18)
    ax2.set_ylim(ISOTOPE_Y_MIN, -42)
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.tick_params(axis='both', which='major', labelsize=14, width=2)
    setup_axis_style(ax2)
    ax2.set_xlim(start_date, end_date)
    ax2.xaxis.set_major_locator(mdates.YearLocator(1))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax2.text(-0.08, 1.02, _subplot_label(1), transform=ax2.transAxes, fontsize=16, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='left')
    ax2.legend(handles=legend_elements, fontsize=10, framealpha=0, loc='lower left')
    
    plt.tight_layout()
    
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, filename), transparent=True, dpi=300, bbox_inches='tight')
    
    plot_elapsed = time.time() - plot_start
    logger.info("Plot saved: %s (%.1fs)", filename, plot_elapsed)
    
    if show:
        plt.show()
    
    return fig


def plot_optimization_comparison(
    condition_results: Dict[str, Any],
    flux_data: Dict[str, np.ndarray],
    obs_targets: Dict[str, Dict[str, float]],
    figs_dir: str = 'figs',
    filename: str = 'isoALBM_optimization_comparison.png',
    show: bool = True
) -> plt.Figure:
    """
    Plot multi-condition optimization comparison (dissolved and surface emitted δ¹³C).
    
    Parameters:
    -----------
    condition_results : dict
        Per-condition dict with 'params', 'timeseries', 'color', 'linestyle'
    flux_data : dict
        Flux data (must contain 'time_array')
    obs_targets : dict
        Observational targets
    figs_dir : str
        Output directory
    filename : str
        Output filename
    show : bool
        Whether to display plot
    
    Returns:
    --------
    fig : matplotlib Figure
    """
    plot_start = time.time()
    logger.info("Generating multi-condition comparison figure")
    fig, axes = plt.subplots(2, 1, figsize=(16, 10), sharex=True)
    _apply_figure_style(fig, axes)
    fig.tight_layout(rect=[0, 0, 0.85, 1])
    time_array = flux_data['time_array']
    time_pd = pd.DatetimeIndex(time_array)
    start_date = time_array[0]
    end_date = time_array[-1]
    ax1 = axes[0]
    excluded = {'No Winter Obs.', 'No Winter Emission Obs.'}
    def _display_name(name):
        if name == 'All Observations':
            return 'Optimization: Obs.'
        if name == 'Temperature-Based':
            return 'Optimization: Temp.'
        return name
    def _trace_color(condition_name, result):
        return CB_BLUE if condition_name == 'All Observations' else result['color']
    plot_conditions = {k: v for k, v in condition_results.items() if k not in excluded}
    for condition_name, result in plot_conditions.items():
        ax1.plot(time_array, result['timeseries']['delta_dissolved'],
                color=_trace_color(condition_name, result), linestyle=result['linestyle'], linewidth=2,
                label=_display_name(condition_name))
    max_year = time_pd.year.max()
    for year in range(time_pd.year.min(), max_year + 1):
        for start_ts, end_ts, season_key in _observation_season_ranges(year, max_year):
            obs_key = 'dissolved_summer' if season_key == 'summer' else 'dissolved_winter'
            obs = obs_targets[obs_key]
            color = OBS_SUMMER_COLOR if season_key == 'summer' else OBS_WINTER_COLOR
            ax1.fill_between([start_ts, end_ts], obs['mean'] - obs['std'], obs['mean'] + obs['std'],
                            color=color, alpha=0.2, zorder=5)
            ax1.hlines(y=obs['mean'], xmin=start_ts, xmax=end_ts, color=color, linestyle='-', linewidth=3, zorder=10)
    ax1.set_ylabel('Dissolved δ¹³C-CH₄ (‰)', fontsize=16)
    # y-axis autoscales based on data
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.tick_params(axis='both', which='major', labelsize=14, width=2)
    setup_axis_style(ax1)
    ax1.set_xlim(start_date, end_date)
    ax1.text(-0.08, 1.02, _subplot_label(0), transform=ax1.transAxes, fontsize=16, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='left')
    model_legend_elements = []
    for condition_name, result in plot_conditions.items():
        model_legend_elements.append(
            Line2D([0], [0], color=_trace_color(condition_name, result), linestyle=result['linestyle'],
                   linewidth=2, label=_display_name(condition_name))
        )
    obs_legend_elements = [
        Line2D([0], [0], color=OBS_SUMMER_COLOR, linewidth=3, label='Obs. Summer Mean'),
        Patch(facecolor=OBS_SUMMER_COLOR, alpha=0.2, label='Obs. Summer ±1σ'),
        Line2D([0], [0], color=OBS_WINTER_COLOR, linewidth=3, label='Obs. Winter Mean'),
        Patch(facecolor=OBS_WINTER_COLOR, alpha=0.2, label='Obs. Winter ±1σ')
    ]
    sep_isoalbm = Line2D([0], [0], color='none', label='isoALBM')
    hline = _LegendSeparatorLine([0, 1], [0, 0], color=CB_BLACK, linewidth=1, label=' ')
    sep_obs = Line2D([0], [0], color='none', label='Observations')
    all_legend_elements = [sep_isoalbm] + model_legend_elements + [hline] + [sep_obs] + obs_legend_elements
    leg1 = ax1.legend(handles=all_legend_elements, loc='center left', bbox_to_anchor=(1.02, 0.5),
                      fontsize=10, framealpha=0.9, edgecolor='black', frameon=True,
                      handler_map={_LegendSeparatorLine: _FullWidthLineHandler()})
    
    # --- Subplot 2: Surface Emissions ---
    ax2 = axes[1]
    for condition_name, result in plot_conditions.items():
        ax2.plot(time_array, result['timeseries']['delta_emission'],
                color=_trace_color(condition_name, result), linestyle=result['linestyle'], linewidth=2,
                label=_display_name(condition_name))
    for year in range(time_pd.year.min(), max_year + 1):
        for start_ts, end_ts, season_key in _observation_season_ranges(year, max_year):
            obs_key = 'emissions_summer' if season_key == 'summer' else 'emissions_winter'
            obs = obs_targets[obs_key]
            color = OBS_SUMMER_COLOR if season_key == 'summer' else OBS_WINTER_COLOR
            ax2.fill_between([start_ts, end_ts], obs['mean'] - obs['std'], obs['mean'] + obs['std'],
                            color=color, alpha=0.2, zorder=5)
            ax2.hlines(y=obs['mean'], xmin=start_ts, xmax=end_ts, color=color, linestyle='-', linewidth=3, zorder=10)
    ax2.set_ylabel('Surface Emitted δ¹³C-CH₄ (‰)', fontsize=16)
    ax2.set_xlabel('Year', fontsize=16)
    # y-axis autoscales based on data
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.tick_params(axis='both', which='major', labelsize=14, width=2)
    setup_axis_style(ax2)
    ax2.xaxis.set_major_locator(mdates.YearLocator(1))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax2.text(-0.08, 1.02, _subplot_label(1), transform=ax2.transAxes, fontsize=16, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='left')
    model_legend_elements_2 = []
    for condition_name, result in plot_conditions.items():
        model_legend_elements_2.append(
            Line2D([0], [0], color=_trace_color(condition_name, result), linestyle=result['linestyle'],
                   linewidth=2, label=_display_name(condition_name))
        )
    obs_legend_elements_2 = [
        Line2D([0], [0], color=OBS_SUMMER_COLOR, linewidth=3, label='Obs. Summer Mean'),
        Patch(facecolor=OBS_SUMMER_COLOR, alpha=0.2, label='Obs. Summer ±1σ'),
        Line2D([0], [0], color=OBS_WINTER_COLOR, linewidth=3, label='Obs. Winter Mean'),
        Patch(facecolor=OBS_WINTER_COLOR, alpha=0.2, label='Obs. Winter ±1σ')
    ]
    sep_isoalbm_2 = Line2D([0], [0], color='none', label='isoALBM')
    hline_2 = _LegendSeparatorLine([0, 1], [0, 0], color=CB_BLACK, linewidth=1, label=' ')
    sep_obs_2 = Line2D([0], [0], color='none', label='Observations')
    all_legend_elements_2 = [sep_isoalbm_2] + model_legend_elements_2 + [hline_2] + [sep_obs_2] + obs_legend_elements_2
    leg2 = ax2.legend(handles=all_legend_elements_2, loc='center left', bbox_to_anchor=(1.02, 0.5),
                      fontsize=10, framealpha=0.9, edgecolor='black', frameon=True,
                      handler_map={_LegendSeparatorLine: _FullWidthLineHandler()})
    # Note: tight_layout() does not accept bbox_extra_artists (that is for savefig only)
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, filename), transparent=True, dpi=300, bbox_inches='tight')
    
    plot_elapsed = time.time() - plot_start
    logger.info("Figure saved: %s (%.1fs)", filename, plot_elapsed)
    if show:
        plt.show()
    return fig


def plot_dissolved_pool_timeseries(
    condition_results: Dict[str, Any],
    flux_data: Dict[str, np.ndarray],
    figs_dir: str = 'figs',
    filename: str = 'isoALBM_dissolved_ch4_pool.png',
    show: bool = True
) -> plt.Figure:
    """
    Plot the internally modeled dissolved CH4 pool above the isotope threshold.
    """
    if not condition_results:
        return _no_data_figure("No isotope condition results")

    preferred_name = 'Initial guess' if 'Initial guess' in condition_results else next(iter(condition_results))
    timeseries = condition_results[preferred_name].get('timeseries', {})
    if 'ch4_conc_mb' not in timeseries:
        return _no_data_figure("No dissolved CH4 pool in isotope results")

    time_array = np.asarray(flux_data['time_array'])
    ch4_pool = np.asarray(timeseries['ch4_conc_mb'], dtype=float)
    n = min(len(time_array), len(ch4_pool))
    if n == 0:
        return _no_data_figure("No dissolved CH4 pool data")

    time_array = time_array[:n]
    ch4_pool = ch4_pool[:n]
    positive_pool = np.where(ch4_pool > 0.0, ch4_pool, np.nan)
    valid = np.asarray(
        timeseries.get('dissolved_concentration_valid', np.isfinite(positive_pool)),
        dtype=bool
    )[:n]
    threshold = timeseries.get('dissolved_concentration_filter_min', None)
    if threshold is not None and np.isfinite(threshold):
        threshold = float(threshold)
    else:
        threshold = None
    plotted_pool = positive_pool.copy()
    if threshold is not None:
        plotted_pool[plotted_pool < threshold] = np.nan

    fig, ax = plt.subplots(figsize=(14, 5))
    _apply_figure_style(fig, ax)
    ax.plot(
        time_array,
        plotted_pool,
        color=CB_BLUE,
        linewidth=2,
        label='ALBM dissolved CH$_4$ pool'
    )
    if threshold is not None and threshold > 0.0:
        ax.axhline(
            threshold,
            color=CB_VERMILLION,
            linestyle='--',
            linewidth=2,
            label=f'Dissolved concentration threshold, {threshold:.0e} kg m$^{{-3}}$'
        )

    if np.any(np.isfinite(plotted_pool)):
        ax.set_yscale('log')
    ax.set_ylabel('Dissolved CH$_4$ pool (kg m$^{-3}$)', fontsize=16)
    ax.set_xlabel('Year', fontsize=16)
    ax.grid(True, which='both', linestyle='--', alpha=0.35)
    ax.tick_params(axis='both', which='major', labelsize=14, width=2)
    setup_axis_style(ax)
    ax.set_xlim(time_array[0], time_array[-1])
    ax.xaxis.set_major_locator(mdates.YearLocator(1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.legend(loc='upper right', fontsize=11, framealpha=0.9, edgecolor='black', frameon=True)

    plt.tight_layout()
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, filename), transparent=True, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    return fig


def _observation_overlap_mask(time_array: np.ndarray) -> np.ndarray:
    """
    Boolean mask True for timesteps that fall in observation periods (summer and winter).
    Summer: Jun 21–Sep 22; Winter: Dec 21–Dec 31 and Jan 1–Mar 20 (next year).
    """
    time_pd = pd.DatetimeIndex(time_array)
    mask = np.zeros(len(time_pd), dtype=bool)
    min_year = int(time_pd.year.min())
    max_year = int(time_pd.year.max())
    for year in range(min_year, max_year + 1):
        for start_ts, end_ts, _ in _observation_season_ranges(year, max_year):
            mask |= (time_pd >= start_ts) & (time_pd <= end_ts)
    return mask


def _observation_overlap_periods(time_array: np.ndarray) -> List[Tuple[str, np.ndarray]]:
    """
    List of (period_name, mask) for each measurement overlap period.
    Summer: Jun 21–Sep 22; Winter: Dec 21–Dec 31 and Jan 1–Mar 20 (next year).
    """
    time_pd = pd.DatetimeIndex(time_array)
    periods: List[Tuple[str, np.ndarray]] = []
    min_year = int(time_pd.year.min())
    max_year = int(time_pd.year.max())
    for year in range(min_year, max_year + 1):
        winter_masks = []
        for start_ts, end_ts, season_key in _observation_season_ranges(year, max_year):
            mask = (time_pd >= start_ts) & (time_pd <= end_ts)
            if season_key == 'summer' and np.any(mask):
                periods.append((f'Summer {year}', mask))
            elif season_key in ('winter_dec', 'winter_jan'):
                winter_masks.append(mask)
        if winter_masks:
            combined = winter_masks[0].copy()
            for m in winter_masks[1:]:
                combined |= m
            if np.any(combined):
                periods.append((f'Winter {year}-{year + 1}', combined))
    return periods


def plot_overlap_violins(
    condition_results: Dict[str, Any],
    flux_data: Dict[str, Any],
    obs_targets: Optional[Dict[str, Dict[str, float]]] = None,
    figs_dir: str = 'figs',
    filename: str = 'isoALBM_overlap_violins.png',
    show: bool = True
) -> plt.Figure:
    """
    Violin plots of modeled δ¹³C per observation-overlap period: Initial guess vs Temperature-Based.
    Two stacked subplots: dissolved CH4 (top), emitted CH4 (bottom). Each period has two violins.
    If obs_targets is provided, draws a flat bar (horizontal band) for observation mean ± std per period.

    Parameters
    ----------
    condition_results : dict
        Must contain 'Initial guess' and/or 'Temperature-Based' with timeseries and color/linestyle.
    flux_data : dict
        Must contain 'time_array'.
    obs_targets : dict, optional
        Observation targets (dissolved_summer, dissolved_winter, emissions_summer, emissions_winter).
    figs_dir : str
        Directory to save figures.
    show : bool
        Whether to display the plot.

    Returns
    -------
    matplotlib.figure.Figure
        The figure (or a minimal no-data figure if no conditions or no observation periods).
    """
    names = ['Initial guess', 'Temperature-Based']
    subset = {k: condition_results[k] for k in names if k in condition_results}
    if len(subset) == 0:
        logger.warning("No 'Initial guess' or 'Temperature-Based' in condition_results; skipping overlap violins.")
        return _no_data_figure("No conditions to plot")
    active_names = [n for n in names if n in subset]
    show_obs = obs_targets is not None and 'dissolved_summer' in obs_targets
    if show_obs:
        active_names = active_names + ['Observations']
    OBS_VIOLIN_COLOR = CB_GRAY
    OBS_VIOLIN_MARKER = 'D'
    DEFAULT_VIOLIN_MARKERS = {
        'Initial guess': 'o',
        'Temperature-Based': 'P',
        'Observations': OBS_VIOLIN_MARKER,
    }

    def _marker_for_name(name):
        if name == 'Observations':
            return OBS_VIOLIN_MARKER
        return subset[name].get('marker', DEFAULT_VIOLIN_MARKERS.get(name, 'o'))

    time_array = flux_data['time_array']
    periods = _observation_overlap_periods(time_array)
    if not periods:
        return _no_data_figure("No observation periods")

    def _obs_for_period(period_name, obs):
        is_summer = period_name.strip().lower().startswith('summer')
        if is_summer:
            return obs['dissolved_summer'], obs['emissions_summer']
        return obs['dissolved_winter'], obs['emissions_winter']

    n_periods = len(periods)
    n_conditions = len(active_names)
    dissolved_flat = []
    emitted_flat = []
    colors_flat = []
    markers_flat = []
    for period_name, period_mask in periods:
        for name in active_names:
            if name == 'Observations':
                obs_diss, obs_emis = _obs_for_period(period_name, obs_targets)
                # Small array so violin shows observation mean ± std
                dissolved_flat.append(np.array([obs_diss['mean'] - obs_diss['std'], obs_diss['mean'], obs_diss['mean'] + obs_diss['std']]))
                emitted_flat.append(np.array([obs_emis['mean'] - obs_emis['std'], obs_emis['mean'], obs_emis['mean'] + obs_emis['std']]))
                colors_flat.append(OBS_VIOLIN_COLOR)
                markers_flat.append(_marker_for_name(name))
            else:
                result = subset[name]
                ts = result['timeseries']
                dissolved_flat.append(ts['delta_dissolved'][period_mask])
                emitted_flat.append(ts['delta_emission'][period_mask])
                colors_flat.append(result['color'])
                markers_flat.append(_marker_for_name(name))

    def _draw_violinplot(ax, values_flat, values_positions, values_colors, values_markers):
        finite_items = []
        for idx, values in enumerate(values_flat):
            finite = np.asarray(values, dtype=float)
            finite = finite[np.isfinite(finite)]
            if finite.size > 0:
                finite_items.append((idx, finite))
        if not finite_items:
            return None
        item_indices = [idx for idx, _ in finite_items]
        parts = ax.violinplot(
            [values for _, values in finite_items],
            positions=[values_positions[idx] for idx in item_indices],
            showmeans=True,
            showmedians=True
        )
        for pc, idx in zip(parts['bodies'], item_indices):
            pc.set_facecolor(values_colors[idx])
            pc.set_alpha(0.7)
        for idx, finite in finite_items:
            ax.scatter(
                values_positions[idx],
                float(np.mean(finite)),
                marker=values_markers[idx],
                s=90,
                facecolors='white',
                edgecolors=values_colors[idx],
                linewidths=1.8,
                alpha=0.95,
                zorder=6,
            )
        return parts

    # Positions: two violins per period with a small gap between period groups
    positions = []
    tick_positions = []
    period_x_ranges = []  # (x_left, x_right) for each period for obs bar
    for i in range(n_periods):
        base = i * (n_conditions + 0.5)
        for j in range(n_conditions):
            positions.append(base + j)
        tick_positions.append(base + (n_conditions - 1) / 2.0)
        period_x_ranges.append((base - 0.25, base + n_conditions - 1 + 0.25))
    period_labels = [p[0] for p in periods]

    fig, axes = plt.subplots(2, 1, figsize=(max(10, n_periods * 1.8), 10), sharex=False)
    _apply_figure_style(fig, axes)
    # --- Dissolved CH4 (top) ---
    ax1 = axes[0]
    _draw_violinplot(ax1, dissolved_flat, positions, colors_flat, markers_flat)
    if show_obs:
        for i, (period_name, _) in enumerate(periods):
            x_left, x_right = period_x_ranges[i]
            obs_diss, _ = _obs_for_period(period_name, obs_targets)
            mean, std = obs_diss['mean'], obs_diss['std']
            ax1.fill_between([x_left, x_right], mean - std, mean + std,
                             color=OBS_SUMMER_COLOR if 'summer' in period_name.lower() else OBS_WINTER_COLOR,
                             alpha=0.35, zorder=0)
            ax1.hlines(mean, x_left, x_right, colors='k', linewidths=1.5, linestyles='solid', zorder=1)
    ax1.set_xticks(tick_positions)
    ax1.set_xticklabels(period_labels, rotation=25, ha='right')
    ax1.set_ylabel('Dissolved δ¹³C-CH₄ (‰)', fontsize=14)
    ax1.grid(True, axis='y', alpha=0.3)
    setup_axis_style(ax1)
    ax1.text(-0.08, 1.02, _subplot_label(0), transform=ax1.transAxes, fontsize=16, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='left')
    def _violin_legend_label(name):
        if name == 'Observations':
            return 'Observations'
        if name == 'Temperature-Based':
            return 'Optimization: Temp.'
        return name
    legend_handles = [
        Line2D(
            [0], [0],
            marker=_marker_for_name(n),
            linestyle='none',
            markerfacecolor='white',
            markeredgecolor=subset[n]['color'] if n in subset else OBS_VIOLIN_COLOR,
            markeredgewidth=1.2,
            markersize=7,
            label=_violin_legend_label(n),
        )
        for n in active_names
    ]
    if show_obs:
        legend_handles.append(Patch(facecolor=OBS_SUMMER_COLOR, alpha=0.35, label='Obs. (summer ±1σ)'))
        legend_handles.append(Patch(facecolor=OBS_WINTER_COLOR, alpha=0.35, label='Obs. (winter ±1σ)'))
    ax1.legend(handles=legend_handles, loc='upper right', fontsize=10)

    # --- Emitted CH4 (bottom) ---
    ax2 = axes[1]
    _draw_violinplot(ax2, emitted_flat, positions, colors_flat, markers_flat)
    if show_obs:
        for i, (period_name, _) in enumerate(periods):
            x_left, x_right = period_x_ranges[i]
            _, obs_emis = _obs_for_period(period_name, obs_targets)
            mean, std = obs_emis['mean'], obs_emis['std']
            ax2.fill_between([x_left, x_right], mean - std, mean + std,
                             color=OBS_SUMMER_COLOR if 'summer' in period_name.lower() else OBS_WINTER_COLOR,
                             alpha=0.35, zorder=0)
            ax2.hlines(mean, x_left, x_right, colors='k', linewidths=1.5, linestyles='solid', zorder=1)
    ax2.set_xticks(tick_positions)
    ax2.set_xticklabels(period_labels, rotation=25, ha='right')
    ax2.set_ylabel('Surface Emitted δ¹³C-CH₄ (‰)', fontsize=14)
    ax2.set_xlabel('Observation period', fontsize=14)
    ax2.grid(True, axis='y', alpha=0.3)
    setup_axis_style(ax2)
    ax2.text(-0.08, 1.02, _subplot_label(1), transform=ax2.transAxes, fontsize=16, fontweight='bold',
             verticalalignment='bottom', horizontalalignment='left')
    ax2.legend(handles=legend_handles, loc='upper right', fontsize=10)

    plt.tight_layout()
    os.makedirs(figs_dir, exist_ok=True)
    plt.savefig(os.path.join(figs_dir, filename), transparent=True, dpi=300, bbox_inches='tight')
    if show:
        plt.show()
    return fig


def plot_all(
    albm_data: ALBMData,
    eddy_data: Optional[EddyFluxData] = None,
    figs_dir: str = 'figs',
    show: bool = True
):
    """
    Generate all flux plots (water column, eddy comparison, flux components,
    CH4 budget, and sediment). Isotope comparison plots are generated by the
    configurable runner or manuscript script when condition results are in memory.
    
    Parameters:
    -----------
    albm_data : ALBMData
        ALBM data container
    eddy_data : EddyFluxData, optional
        Eddy covariance data container
    figs_dir : str
        Directory to save figures
    show : bool
        Whether to display plots
    """
    total_start = time.time()
    logger.info("GENERATING FLUX PLOTS")
    os.makedirs(figs_dir, exist_ok=True)
    t0 = time.time()
    plot_water_column_heatmaps(albm_data, figs_dir, show)
    logger.info("Water column heatmaps: %.1fs", time.time() - t0)
    t0 = time.time()
    plot_water_column_doy_climatology_heatmaps(albm_data, figs_dir, show)
    logger.info("Water column DOY climatology heatmaps: %.1fs", time.time() - t0)
    t0 = time.time()
    plot_oxic_and_exchange_heatmaps(albm_data, figs_dir, show)
    logger.info("Oxic and exchange heatmaps: %.1fs", time.time() - t0)
    t0 = time.time()
    plot_eddy_comparison(albm_data, eddy_data, figs_dir, show)
    logger.info("Eddy comparison: %.1fs", time.time() - t0)
    t0 = time.time()
    plot_flux_components(albm_data, figs_dir, show)
    logger.info("Flux components (log): %.1fs", time.time() - t0)
    t0 = time.time()
    plot_flux_components_linear(albm_data, figs_dir, show)
    logger.info("Flux components (linear): %.1fs", time.time() - t0)
    t0 = time.time()
    plot_ch4_budget(albm_data, figs_dir, show)
    logger.info("CH4 budget: %.1fs", time.time() - t0)
    t0 = time.time()
    plot_pool_sources_sinks(albm_data, figs_dir, show)
    logger.info("Pool sources/sinks: %.1fs", time.time() - t0)
    t0 = time.time()
    plot_sediment_temperature(albm_data, figs_dir, show)
    logger.info("Sediment temperature: %.1fs", time.time() - t0)

    logger.info("Skipping isotope comparison plots in plot_all; run an optimization workflow to generate them.")

    total_elapsed = time.time() - total_start
    logger.info("All flux plots saved to %s/ (total: %.1fs)", figs_dir, total_elapsed)


# =============================================================================
# Paper statistics (numeric only, no plots)
# =============================================================================

# Observation date ranges at BTL (update from paper/supplement if needed)
OBS_DISSOLVED_DATE_RANGE = "Summer 2021 – Summer 2024 (Jun 21–Sep 22; Dec 21–Mar 20)"
OBS_EMISSION_DATE_RANGE = "Summer 2021 – Winter 2023-24 (Jun 21–Sep 22; Dec 21–Mar 20)"
# Literature fractionation ranges (Conrad 2005, Hornibrook 2000, Whiticar 1999; optimization bounds)
PAPER_ALPHA_AM_RANGE = (1.000, 1.040)
PAPER_ALPHA_HM_RANGE = (1.030, 1.080)
# Oh et al. (2022) / Whiticar (1999) adopted values used as initial guess (boreal forested)
PAPER_ALPHA_AM_OH2022 = 1.0238
PAPER_ALPHA_HM_OH2022 = 1.0456
# δ¹³C-CH₄ analytical precision (‰), from methods
PAPER_ISOTOPE_PRECISION_PERCML = 0.2
# Lake overturning: temperature gradient threshold (°C) — from ALBM config if available
PAPER_OVERTURN_TEMP_GRADIENT_THRESHOLD = None  # set if known, e.g. 0.01


def _summer_winter_masks(time_array: np.ndarray):
    """Return (summer_mask, winter_mask) for time_array. Summer: Jun 21–Sep 22; Winter: Dec 21–Mar 20."""
    t = pd.DatetimeIndex(time_array)
    n = len(t)
    summer_mask = np.zeros(n, dtype=bool)
    winter_mask = np.zeros(n, dtype=bool)
    for year in range(int(t.year.min()), int(t.year.max()) + 1):
        summer_start = pd.Timestamp(f'{year}-06-21')
        summer_end = pd.Timestamp(f'{year}-09-22')
        summer_mask |= (t >= summer_start) & (t <= summer_end)
        winter_dec_start = pd.Timestamp(f'{year}-12-21')
        winter_jan_end = pd.Timestamp(f'{year + 1}-03-20')
        winter_mask |= ((t >= winter_dec_start) & (t.year == year)) | (
            (t <= winter_jan_end) & (t.year == year + 1))
    return summer_mask, winter_mask


def _dissolved_ch4_observation_point_comparison(
    albm_data: ALBMData,
    dissolved_profile_file: str = 'data/2022_11_21_CH4CO2_DeepShallow.xlsx',
) -> Dict[str, Any]:
    """
    Compare modeled dissolved CH4 against profile observations at matching dates
    and nearest model water-column depths.

    Observations are loaded in kg m-3. Model ch4conc is converted from model molar
    units to kg m-3 using CH4 molecular mass, then reported here as mg L-1 for
    manuscript readability.
    """
    obs = _load_dissolved_ch4_observations(dissolved_profile_file)
    empty = {
        'n_observations': len(obs),
        'n_matched': 0,
        'model_lower_count': np.nan,
        'model_lower_fraction': np.nan,
        'model_lower_percent': np.nan,
        'model_median_mg_l': np.nan,
        'observed_median_mg_l': np.nan,
        'model_q05_mg_l': np.nan,
        'model_q95_mg_l': np.nan,
        'observed_q05_mg_l': np.nan,
        'observed_q95_mg_l': np.nan,
    }
    if obs.empty or not hasattr(albm_data, 'ch4conc') or albm_data.ch4conc is None:
        return empty

    ch4conc = np.asarray(albm_data.ch4conc, dtype=float)
    if ch4conc.ndim != 2 or ch4conc.size == 0:
        return empty

    time_index = pd.DatetimeIndex(albm_data.time_array).normalize()
    n_time = min(len(time_index), ch4conc.shape[0])
    time_index = time_index[:n_time]
    ch4conc = ch4conc[:n_time, :]
    if n_time == 0:
        return empty

    date_to_index = {date: idx for idx, date in enumerate(time_index)}
    depths = np.asarray(albm_data.water_depth_array, dtype=float)
    if depths.size != ch4conc.shape[1]:
        depths = np.arange(ch4conc.shape[1], dtype=float) * 0.1

    mw_ch4_kg_per_mol = 16.04e-3
    rows: List[Tuple[float, float]] = []
    for _, row in obs.iterrows():
        obs_date = pd.Timestamp(row['datetime']).normalize()
        model_idx = date_to_index.get(obs_date)
        if model_idx is None:
            continue
        depth_m = float(row['depth_m'])
        obs_kg_m3 = float(row['ch4_kg_m3'])
        if not np.isfinite(depth_m) or not np.isfinite(obs_kg_m3):
            continue
        depth_idx = int(np.nanargmin(np.abs(depths - depth_m)))
        model_kg_m3 = float(ch4conc[model_idx, depth_idx] * mw_ch4_kg_per_mol)
        if np.isfinite(model_kg_m3):
            rows.append((model_kg_m3, obs_kg_m3))

    if not rows:
        return empty

    matched = pd.DataFrame(rows, columns=['model_kg_m3', 'observed_kg_m3'])
    model_mg_l = matched['model_kg_m3'].to_numpy() * 1000.0
    observed_mg_l = matched['observed_kg_m3'].to_numpy() * 1000.0
    lower = matched['model_kg_m3'] < matched['observed_kg_m3']
    lower_count = int(lower.sum())
    lower_fraction = float(lower.mean())
    return {
        'n_observations': len(obs),
        'n_matched': int(len(matched)),
        'model_lower_count': lower_count,
        'model_lower_fraction': lower_fraction,
        'model_lower_percent': lower_fraction * 100.0,
        'model_median_mg_l': float(np.nanmedian(model_mg_l)),
        'observed_median_mg_l': float(np.nanmedian(observed_mg_l)),
        'model_q05_mg_l': float(np.nanquantile(model_mg_l, 0.05)),
        'model_q95_mg_l': float(np.nanquantile(model_mg_l, 0.95)),
        'observed_q05_mg_l': float(np.nanquantile(observed_mg_l, 0.05)),
        'observed_q95_mg_l': float(np.nanquantile(observed_mg_l, 0.95)),
    }


def compute_and_print_paper_statistics(
    albm_data: ALBMData,
    eddy_data: Optional[EddyFluxData],
    condition_results: Dict[str, Any],
    flux_data: Dict[str, Any],
    obs_targets: Dict[str, Dict[str, float]],
    csv_path: Optional[str] = None,
) -> pd.DataFrame:
    """
    Compute and print all paper statistics (A–J) to the terminal. No plots.
    Run after multi_condition_optimization to fill condition_results. Section A
    (ERA5 trends) is printed by climateChangeAnalysis.py when run separately.

    Returns a tidy DataFrame with one statistic per row. If csv_path is provided,
    the same table is written to disk.
    """
    time_arr = flux_data['time_array']
    summer_mask, winter_mask = _summer_winter_masks(time_arr)
    rows: List[Dict[str, Any]] = []

    def _clean_scalar(value: Any) -> Any:
        if isinstance(value, (np.floating, np.integer)):
            return value.item()
        if isinstance(value, np.ndarray):
            if value.size == 1:
                return _clean_scalar(value.item())
            return value.tolist()
        return value

    def _add_stat(
        section: str,
        statistic_id: str,
        description: str,
        value: Any,
        units: str = "",
        manuscript_line: str = "",
        notes: str = "",
    ) -> None:
        rows.append({
            'source_file': 'flux_plots.py',
            'section': section,
            'statistic_id': statistic_id,
            'description': description,
            'value': _clean_scalar(value),
            'units': units,
            'manuscript_line': manuscript_line,
            'notes': notes,
        })

    def _mean_std(arr):
        if arr.size == 0:
            return np.nan, np.nan
        m = float(np.mean(arr))
        s = float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0
        return m, s

    # ----- B. Isotope mechanics (literature / Oh et al. 2022) -----
    print("\n" + "="*70)
    print("PAPER STATISTICS — B. Isotope Mechanics (Methods)")
    print("="*70)
    print("B.1 AM fractionation factor range (literature):  α_AM ∈ [%.3f, %.3f]" % PAPER_ALPHA_AM_RANGE)
    print("B.2 HM fractionation factor range (literature):  α_HM ∈ [%.3f, %.3f]" % PAPER_ALPHA_HM_RANGE)
    print("B.3 α_AM adopted (Oh et al. 2022, boreal):       %.4f" % PAPER_ALPHA_AM_OH2022)
    print("B.4 α_HM adopted (Oh et al. 2022, boreal):       %.4f" % PAPER_ALPHA_HM_OH2022)
    _add_stat('B. Isotope Mechanics', 'B.1_alpha_am_range_min',
              'AM fractionation factor literature range minimum',
              PAPER_ALPHA_AM_RANGE[0], 'unitless alpha', '299')
    _add_stat('B. Isotope Mechanics', 'B.1_alpha_am_range_max',
              'AM fractionation factor literature range maximum',
              PAPER_ALPHA_AM_RANGE[1], 'unitless alpha', '299')
    _add_stat('B. Isotope Mechanics', 'B.2_alpha_hm_range_min',
              'HM fractionation factor literature range minimum',
              PAPER_ALPHA_HM_RANGE[0], 'unitless alpha', '299')
    _add_stat('B. Isotope Mechanics', 'B.2_alpha_hm_range_max',
              'HM fractionation factor literature range maximum',
              PAPER_ALPHA_HM_RANGE[1], 'unitless alpha', '299')
    _add_stat('B. Isotope Mechanics', 'B.3_alpha_am_oh2022',
              'Adopted Oh et al. boreal AM fractionation factor',
              PAPER_ALPHA_AM_OH2022, 'unitless alpha', '299')
    _add_stat('B. Isotope Mechanics', 'B.4_alpha_hm_oh2022',
              'Adopted Oh et al. boreal HM fractionation factor',
              PAPER_ALPHA_HM_OH2022, 'unitless alpha', '299')
    print("="*70)

    # ----- C. Observation date ranges -----
    print("\n" + "="*70)
    print("PAPER STATISTICS — C. Model Optimization (Observation Date Ranges)")
    print("="*70)
    print("C.1 Dissolved surface measurement date range:    %s" % OBS_DISSOLVED_DATE_RANGE)
    print("C.2 Surface emission measurement date range:    %s" % OBS_EMISSION_DATE_RANGE)
    _add_stat('C. Model Optimization', 'C.1_dissolved_date_range',
              'Dissolved surface isotope observation date range',
              OBS_DISSOLVED_DATE_RANGE, '', '413')
    _add_stat('C. Model Optimization', 'C.2_emission_date_range',
              'Surface-emission isotope observation date range',
              OBS_EMISSION_DATE_RANGE, '', '413')
    print("="*70)

    # ----- D. Methane flux (EC + ALBM) -----
    print("\n" + "="*70)
    print("PAPER STATISTICS — D. Results: Methane Flux at Big Trail Lake")
    print("Model period: March 2021 – December 2024. Summer = Jun 21–Sep 22, Winter = Dec 21–Mar 20.")
    print("="*70)
    kg_m2_s_to_mg_m2_d = 1e6 * 86400  # kg m-2 s-1 -> mg m-2 d-1
    if eddy_data is not None:
        eddy_df = eddy_data.ch4_flux_daily_avg
        model_start, model_end = albm_data.time_array[0], albm_data.time_array[-1]
        overlap_start = max(eddy_df['datetime'].min(), pd.Timestamp(model_start))
        overlap_end = min(eddy_df['datetime'].max(), pd.Timestamp(model_end))
        eddy_mask_overlap = (eddy_df['datetime'] >= overlap_start) & (eddy_df['datetime'] <= overlap_end)
        eddy_dates = eddy_df.loc[eddy_mask_overlap, 'datetime'].values
        eddy_flux = eddy_df.loc[eddy_mask_overlap, 'ch4_flux_daily_avg_kg_m2_s'].values
        model_interp = np.interp(
            mdates.date2num(pd.DatetimeIndex(eddy_dates)),
            mdates.date2num(albm_data.time_array),
            albm_data.total_flux
        )
        eddy_summer, eddy_winter = _summer_winter_masks(pd.DatetimeIndex(eddy_dates).to_numpy())
        mean_ec_summer = np.nanmean(eddy_flux[eddy_summer]) if np.any(eddy_summer) else np.nan
        mean_ec_winter = np.nanmean(eddy_flux[eddy_winter]) if np.any(eddy_winter) else np.nan
        daily_std_ec_summer = np.nanstd(eddy_flux[eddy_summer], ddof=1) if np.sum(eddy_summer) > 1 else np.nan
        daily_std_ec_winter = np.nanstd(eddy_flux[eddy_winter], ddof=1) if np.sum(eddy_winter) > 1 else np.nan
        valid = np.isfinite(model_interp) & np.isfinite(eddy_flux)
        rmse = np.sqrt(np.mean((model_interp[valid] - eddy_flux[valid]) ** 2)) if np.sum(valid) > 0 else np.nan
        res = model_interp - eddy_flux
        mean_res_summer = np.nanmean(res[eddy_summer]) if np.any(eddy_summer) else np.nan
        mean_res_winter = np.nanmean(res[eddy_winter]) if np.any(eddy_winter) else np.nan
        comp = pd.DataFrame({
            'datetime': pd.DatetimeIndex(eddy_dates),
            'eddy_flux': eddy_flux,
            'model_flux': model_interp,
        }).dropna(subset=['eddy_flux', 'model_flux'])
        weekly = comp.set_index('datetime').resample('7D', origin=overlap_start).mean()
        weekly = weekly.dropna(subset=['eddy_flux', 'model_flux'])
        weekly_summer, weekly_winter = _summer_winter_masks(weekly.index.to_numpy())
        weekly_ec_summer = weekly.loc[weekly_summer, 'eddy_flux'].to_numpy()
        weekly_ec_winter = weekly.loc[weekly_winter, 'eddy_flux'].to_numpy()
        weekly_model = weekly['model_flux'].to_numpy()
        weekly_ec = weekly['eddy_flux'].to_numpy()
        std_ec_summer = np.nanstd(weekly_ec_summer, ddof=1) if weekly_ec_summer.size > 1 else np.nan
        std_ec_winter = np.nanstd(weekly_ec_winter, ddof=1) if weekly_ec_winter.size > 1 else np.nan
        weekly_valid = np.isfinite(weekly_model) & np.isfinite(weekly_ec)
        std_model = np.nanstd(weekly_model[weekly_valid], ddof=1) if np.sum(weekly_valid) > 1 else np.nan
        std_ec = np.nanstd(weekly_ec[weekly_valid], ddof=1) if np.sum(weekly_valid) > 1 else np.nan
        var_underestimate_pct = (1 - std_model / std_ec) * 100 if std_ec and std_ec > 0 else np.nan
        pct_resid = np.full(weekly_ec.shape, np.nan, dtype=float)
        pct_ok = weekly_valid & (np.abs(weekly_ec) > 1e-15)
        pct_resid[pct_ok] = (weekly_model[pct_ok] - weekly_ec[pct_ok]) / weekly_ec[pct_ok] * 100.0
        weekly_flux_bias_pct = np.nanmean(pct_resid[pct_ok]) if np.any(pct_ok) else np.nan
        weekly_flux_bias_pct_std = np.nanstd(pct_resid[pct_ok], ddof=1) if np.sum(pct_ok) > 1 else np.nan
        print("D.1  Mean summer EC flux:        %.4e kg m⁻² s⁻¹  (%.2f mg CH₄ m⁻² d⁻¹)" % (
            mean_ec_summer, mean_ec_summer * kg_m2_s_to_mg_m2_d))
        print("D.2  Mean winter EC flux:        %.4e kg m⁻² s⁻¹  (%.2f mg CH₄ m⁻² d⁻¹)" % (
            mean_ec_winter, mean_ec_winter * kg_m2_s_to_mg_m2_d))
        print("D.3  Summer EC flux variability (weekly std): %.4e kg m⁻² s⁻¹" % std_ec_summer)
        print("D.4  Winter EC flux variability (weekly std): %.4e kg m⁻² s⁻¹" % std_ec_winter)
        print("D.5  ALBM vs EC RMSE:            %.4e kg m⁻² s⁻¹" % rmse)
        print("D.6  Mean model–obs residual (summer): %.4e kg m⁻² s⁻¹" % mean_res_summer)
        print("D.7  Mean model–obs residual (winter): %.4e kg m⁻² s⁻¹" % mean_res_winter)
        print("D.8  Weekly flux percent bias (model vs EC): %.1f%% ± %.1f%%" % (
            weekly_flux_bias_pct, weekly_flux_bias_pct_std))
        print("D.8b Weekly variability underestimate (%%): %.1f%%" % var_underestimate_pct)
        print("D.9  std(EC)=%.4e  std(ALBM)=%.4e kg m⁻² s⁻¹" % (std_ec, std_model))
        _add_stat('D. Methane Flux', 'D.1_mean_summer_ec_flux',
                  'Mean summer EC methane flux', mean_ec_summer, 'kg m-2 s-1', '464',
                  'Summer = Jun 21-Sep 22; daily EC means.')
        _add_stat('D. Methane Flux', 'D.1_mean_summer_ec_flux_mg_m2_d',
                  'Mean summer EC methane flux', mean_ec_summer * kg_m2_s_to_mg_m2_d,
                  'mg CH4 m-2 d-1', '464')
        _add_stat('D. Methane Flux', 'D.2_mean_winter_ec_flux',
                  'Mean winter EC methane flux', mean_ec_winter, 'kg m-2 s-1', '464',
                  'Winter = Dec 21-Mar 20; daily EC means.')
        _add_stat('D. Methane Flux', 'D.2_mean_winter_ec_flux_mg_m2_d',
                  'Mean winter EC methane flux', mean_ec_winter * kg_m2_s_to_mg_m2_d,
                  'mg CH4 m-2 d-1', '464')
        _add_stat('D. Methane Flux', 'D.3_weekly_summer_ec_std',
                  'Summer EC methane flux variability from consecutive non-overlapping weekly bins',
                  std_ec_summer, 'kg m-2 s-1', '464')
        _add_stat('D. Methane Flux', 'D.4_weekly_winter_ec_std',
                  'Winter EC methane flux variability from consecutive non-overlapping weekly bins',
                  std_ec_winter, 'kg m-2 s-1', '464')
        _add_stat('D. Methane Flux', 'D.3a_daily_summer_ec_std',
                  'Summer EC methane flux variability from daily bins',
                  daily_std_ec_summer, 'kg m-2 s-1', '464',
                  'Reported for traceability; manuscript language asks for weekly bins.')
        _add_stat('D. Methane Flux', 'D.4a_daily_winter_ec_std',
                  'Winter EC methane flux variability from daily bins',
                  daily_std_ec_winter, 'kg m-2 s-1', '464',
                  'Reported for traceability; manuscript language asks for weekly bins.')
        _add_stat('D. Methane Flux', 'D.5_daily_model_ec_rmse',
                  'RMSE between daily ALBM total flux and daily EC flux',
                  rmse, 'kg m-2 s-1', '472')
        _add_stat('D. Methane Flux', 'D.6_summer_mean_residual',
                  'Mean daily ALBM minus EC residual in summer',
                  mean_res_summer, 'kg m-2 s-1', '472')
        _add_stat('D. Methane Flux', 'D.7_winter_mean_residual',
                  'Mean daily ALBM minus EC residual in winter',
                  mean_res_winter, 'kg m-2 s-1', '472')
        _add_stat('D. Methane Flux', 'D.8_weekly_flux_bias_pct',
                  'Mean weekly percent bias, ALBM relative to EC',
                  weekly_flux_bias_pct, 'percent', '472',
                  'Computed as mean((ALBM - EC) / EC * 100) for weekly bins with nonzero EC.')
        _add_stat('D. Methane Flux', 'D.8_weekly_flux_bias_pct_std',
                  'Standard deviation of weekly percent bias, ALBM relative to EC',
                  weekly_flux_bias_pct_std, 'percentage points', '472')
        _add_stat('D. Methane Flux', 'D.8b_weekly_variability_underestimate_pct',
                  'Underestimate of weekly flux variability by ALBM relative to EC',
                  var_underestimate_pct, 'percent', '472',
                  'Computed as (1 - std(ALBM_weekly) / std(EC_weekly)) * 100.')
        _add_stat('D. Methane Flux', 'D.9_weekly_ec_std_all',
                  'Standard deviation of all weekly EC flux bins',
                  std_ec, 'kg m-2 s-1', '472')
        _add_stat('D. Methane Flux', 'D.9_weekly_albm_std_all',
                  'Standard deviation of all weekly ALBM flux bins',
                  std_model, 'kg m-2 s-1', '472')
    else:
        print("D.1–D.9  (No eddy data — skip flux comparison)")
    # Ebullition / diffusion ratio (summer)
    eflux = albm_data.eflux_data + albm_data.ch4upb_data
    dflux = albm_data.dflux_data
    e_summer = np.sum(eflux[summer_mask])
    d_summer = np.sum(dflux[summer_mask])
    ebull_diff_ratio = e_summer / d_summer if d_summer > 0 else np.nan
    obs_summer_diff_frac = obs_targets.get('emissions_summer', {}).get('diffusion_fraction', np.nan)
    obs_winter_diff_frac = obs_targets.get('emissions_winter', {}).get('diffusion_fraction', np.nan)
    obs_summer_ebull_diff_ratio = ((1.0 - obs_summer_diff_frac) / obs_summer_diff_frac
                                   if obs_summer_diff_frac and obs_summer_diff_frac > 0 else np.nan)
    obs_winter_ebull_diff_ratio = ((1.0 - obs_winter_diff_frac) / obs_winter_diff_frac
                                   if obs_winter_diff_frac and obs_winter_diff_frac > 0 else np.nan)
    print("D.10 Ebullition-to-diffusion ratio (summer): %.2f" % ebull_diff_ratio)
    print("D.11 Lake overturning temp. gradient threshold: %s" % (
        str(PAPER_OVERTURN_TEMP_GRADIENT_THRESHOLD) if PAPER_OVERTURN_TEMP_GRADIENT_THRESHOLD is not None else "N/A (see ALBM config)"))
    # Observed vs modeled diffusion:ebullition — we don't have observed component split; report ALBM ratio
    print("D.12 ALBM-modeled diffusion:ebullition ratio (summer): 1 : %.2f" % ebull_diff_ratio)
    print("D.13 Observed summer ebullition-to-diffusion ratio implied by isotope target: %.2f" %
          obs_summer_ebull_diff_ratio)
    _add_stat('D. Methane Flux', 'D.10_model_summer_ebullition_diffusion_ratio',
              'ALBM summer surface ebullition plus ice-leakage flux divided by diffusive flux',
              ebull_diff_ratio, 'ratio', '167;472;532;575')
    _add_stat('D. Methane Flux', 'D.11_overturn_threshold',
              'Lake overturning temperature-gradient threshold',
              PAPER_OVERTURN_TEMP_GRADIENT_THRESHOLD, 'deg C', '480',
              'No fixed temperature-gradient threshold is reported by this Python workflow; the ALBM source thermal_mod.f90 uses kinetic/potential energy mixing logic instead.')
    _add_stat('D. Methane Flux', 'D.12_model_summer_diffusion_ebullition_ratio_text',
              'Text version of ALBM modeled summer diffusion:ebullition ratio',
              f'1:{ebull_diff_ratio:.2f}', 'ratio', '472')
    _add_stat('D. Methane Flux', 'D.13_observed_summer_ebullition_diffusion_ratio',
              'Observed summer ebullition-to-diffusion ratio implied by diffusion fraction',
              obs_summer_ebull_diff_ratio, 'ratio', '532',
              'Computed as (1 - diffusion_fraction) / diffusion_fraction from isotope target metadata.')
    _add_stat('D. Methane Flux', 'D.14_observed_winter_ebullition_diffusion_ratio',
              'Observed winter ebullition-to-diffusion ratio implied by diffusion fraction',
              obs_winter_ebull_diff_ratio, 'ratio', '',
              'Computed as (1 - diffusion_fraction) / diffusion_fraction from isotope target metadata.')

    dissolved_ch4_cmp = _dissolved_ch4_observation_point_comparison(albm_data)
    print("D.15 Dissolved CH4 model below observations: %.1f%% (%s/%s matched profile points)" % (
        dissolved_ch4_cmp['model_lower_percent'],
        dissolved_ch4_cmp['model_lower_count'],
        dissolved_ch4_cmp['n_matched'],
    ))
    print("D.16 Dissolved CH4 median concentration: model %.2f mg L-1 vs observed %.2f mg L-1" % (
        dissolved_ch4_cmp['model_median_mg_l'],
        dissolved_ch4_cmp['observed_median_mg_l'],
    ))
    dissolved_notes = (
        'Dissolved CH4 observations from 2022_11_21_CH4CO2_DeepShallow.xlsx; '
        'model values use the same calendar day and nearest 0.1 m water-column layer; '
        'model ch4conc converted to mg L-1 using CH4 molecular mass.'
    )
    _add_stat('D. Methane Flux', 'D.15_dissolved_ch4_observation_count',
              'Number of dissolved CH4 profile observations loaded',
              dissolved_ch4_cmp['n_observations'], 'count', '', dissolved_notes)
    _add_stat('D. Methane Flux', 'D.16_dissolved_ch4_matched_count',
              'Number of dissolved CH4 observations matched to model date-depth points',
              dissolved_ch4_cmp['n_matched'], 'count', '', dissolved_notes)
    _add_stat('D. Methane Flux', 'D.17_dissolved_ch4_model_below_obs_count',
              'Count of matched dissolved CH4 observations where model is lower than observed',
              dissolved_ch4_cmp['model_lower_count'], 'count', '', dissolved_notes)
    _add_stat('D. Methane Flux', 'D.18_dissolved_ch4_model_below_obs_percent',
              'Percent of matched dissolved CH4 observations where model is lower than observed',
              dissolved_ch4_cmp['model_lower_percent'], 'percent', '', dissolved_notes)
    _add_stat('D. Methane Flux', 'D.19_dissolved_ch4_model_median',
              'Median modeled dissolved CH4 concentration at observation points',
              dissolved_ch4_cmp['model_median_mg_l'], 'mg CH4 L-1', '', dissolved_notes)
    _add_stat('D. Methane Flux', 'D.20_dissolved_ch4_observed_median',
              'Median observed dissolved CH4 concentration at matched observation points',
              dissolved_ch4_cmp['observed_median_mg_l'], 'mg CH4 L-1', '', dissolved_notes)
    _add_stat('D. Methane Flux', 'D.21_dissolved_ch4_model_q05',
              'Fifth percentile modeled dissolved CH4 concentration at observation points',
              dissolved_ch4_cmp['model_q05_mg_l'], 'mg CH4 L-1', '', dissolved_notes)
    _add_stat('D. Methane Flux', 'D.22_dissolved_ch4_model_q95',
              'Ninety-fifth percentile modeled dissolved CH4 concentration at observation points',
              dissolved_ch4_cmp['model_q95_mg_l'], 'mg CH4 L-1', '', dissolved_notes)
    _add_stat('D. Methane Flux', 'D.23_dissolved_ch4_observed_q05',
              'Fifth percentile observed dissolved CH4 concentration at matched observation points',
              dissolved_ch4_cmp['observed_q05_mg_l'], 'mg CH4 L-1', '', dissolved_notes)
    _add_stat('D. Methane Flux', 'D.24_dissolved_ch4_observed_q95',
              'Ninety-fifth percentile observed dissolved CH4 concentration at matched observation points',
              dissolved_ch4_cmp['observed_q95_mg_l'], 'mg CH4 L-1', '', dissolved_notes)
    print("="*70)

    # ----- E. Isotope observational summary -----
    print("\n" + "="*70)
    print("PAPER STATISTICS — E. Results: Methane Isotopes — Observational Summary")
    print("="*70)
    ds, de = obs_targets['dissolved_summer'], obs_targets['emissions_summer']
    dw, ew = obs_targets['dissolved_winter'], obs_targets['emissions_winter']
    print("E.1 Mean ± std dissolved δ¹³C-CH₄, summer: %.2f ± %.2f ‰" % (ds['mean'], ds['std']))
    print("E.2 Mean ± std dissolved δ¹³C-CH₄, winter: %.2f ± %.2f ‰" % (dw['mean'], dw['std']))
    print("E.3 Mean ± std emitted δ¹³C-CH₄, summer:   %.2f ± %.2f ‰" % (de['mean'], de['std']))
    print("E.4 Mean ± std emitted δ¹³C-CH₄, winter:   %.2f ± %.2f ‰" % (ew['mean'], ew['std']))
    print("E.5 δ¹³C-CH₄ measurement precision:       ±%.1f ‰" % PAPER_ISOTOPE_PRECISION_PERCML)
    dissolved_snr = abs(ds['mean'] - dw['mean']) / PAPER_ISOTOPE_PRECISION_PERCML
    emitted_snr = abs(de['mean'] - ew['mean']) / PAPER_ISOTOPE_PRECISION_PERCML
    print("E.6 Seasonal signal-to-noise ratio, dissolved: %.1f" % dissolved_snr)
    print("E.7 Seasonal signal-to-noise ratio, emitted:   %.1f" % emitted_snr)
    _add_stat('E. Isotope Observations', 'E.1_dissolved_summer_mean',
              'Observed dissolved methane isotope mean in summer',
              ds['mean'], 'per mil VPDB', '413;485')
    _add_stat('E. Isotope Observations', 'E.1_dissolved_summer_std',
              'Observed dissolved methane isotope standard deviation in summer',
              ds['std'], 'per mil', '413;485')
    _add_stat('E. Isotope Observations', 'E.2_dissolved_winter_mean',
              'Observed dissolved methane isotope mean in winter',
              dw['mean'], 'per mil VPDB', '413;485')
    _add_stat('E. Isotope Observations', 'E.2_dissolved_winter_std',
              'Observed dissolved methane isotope standard deviation in winter',
              dw['std'], 'per mil', '413;485')
    _add_stat('E. Isotope Observations', 'E.3_emitted_summer_mean',
              'Observed surface-emitted methane isotope mean in summer',
              de['mean'], 'per mil VPDB', '167;413;485;577')
    _add_stat('E. Isotope Observations', 'E.3_emitted_summer_std',
              'Observed surface-emitted methane isotope propagated standard deviation in summer',
              de['std'], 'per mil', '167;413;485;577')
    _add_stat('E. Isotope Observations', 'E.4_emitted_winter_mean',
              'Observed surface-emitted methane isotope mean in winter',
              ew['mean'], 'per mil VPDB', '167;413;485;577')
    _add_stat('E. Isotope Observations', 'E.4_emitted_winter_std',
              'Observed surface-emitted methane isotope propagated standard deviation in winter',
              ew['std'], 'per mil', '167;413;485;577')
    _add_stat('E. Isotope Observations', 'E.5_measurement_precision',
              'Carbon-isotope analytical precision',
              PAPER_ISOTOPE_PRECISION_PERCML, 'per mil 1 sigma', '485')
    _add_stat('E. Isotope Observations', 'E.6_dissolved_seasonal_snr',
              'Dissolved isotope seasonal signal-to-noise ratio',
              dissolved_snr, 'ratio', '485',
              'Computed as abs(summer mean - winter mean) / analytical precision.')
    _add_stat('E. Isotope Observations', 'E.7_emitted_seasonal_snr',
              'Surface-emission isotope seasonal signal-to-noise ratio',
              emitted_snr, 'ratio', '485',
              'Computed as abs(summer mean - winter mean) / analytical precision.')
    print("="*70)

    # ----- F. Case 1 (Initial guess) -----
    case1 = condition_results.get('Initial guess')
    if case1 is not None:
        ts = case1['timeseries']
        d = ts['delta_dissolved']
        print("\n" + "="*70)
        print("PAPER STATISTICS — F. Case 1 (Initial Parameterization, oh2022improved)")
        print("="*70)
        m_s, s_s = _mean_std(d[summer_mask])
        m_w, s_w = _mean_std(d[winter_mask])
        print("F.1 Case 1 dissolved δ¹³C model vs obs, summer: model %.2f ‰  obs %.2f ‰" % (m_s, ds['mean']))
        print("F.2 Case 1 dissolved δ¹³C model vs obs, winter: model %.2f ‰  obs %.2f ‰" % (m_w, dw['mean']))
        print("F.3 Case 1 dissolved variability, summer: model std %.2f  obs std %.2f" % (s_s, ds['std']))
        print("F.4 Case 1 dissolved variability, winter: model std %.2f  obs std %.2f" % (s_w, dw['std']))
        sed_prod = ts.get('C13_sed_prod', ts.get('C13_sed_prod_mean', np.nan))
        sed_mean = np.nanmean(sed_prod) if np.size(sed_prod) else np.nan
        print("F.5 Case 1 sediment production δ¹³C-CH₄,prod: %.2f ‰" % sed_mean)
        _add_stat('F. Case 1', 'F.1_case1_dissolved_summer_model_mean',
                  'Case 1 modeled dissolved methane isotope mean in summer',
                  m_s, 'per mil VPDB', '495')
        _add_stat('F. Case 1', 'F.1_case1_dissolved_summer_obs_mean',
                  'Observed dissolved methane isotope mean in summer used for Case 1 comparison',
                  ds['mean'], 'per mil VPDB', '495')
        _add_stat('F. Case 1', 'F.2_case1_dissolved_winter_model_mean',
                  'Case 1 modeled dissolved methane isotope mean in winter',
                  m_w, 'per mil VPDB', '495')
        _add_stat('F. Case 1', 'F.2_case1_dissolved_winter_obs_mean',
                  'Observed dissolved methane isotope mean in winter used for Case 1 comparison',
                  dw['mean'], 'per mil VPDB', '495')
        _add_stat('F. Case 1', 'F.3_case1_dissolved_summer_model_std',
                  'Case 1 modeled dissolved methane isotope standard deviation in summer',
                  s_s, 'per mil', '495')
        _add_stat('F. Case 1', 'F.3_case1_dissolved_summer_obs_std',
                  'Observed dissolved methane isotope standard deviation in summer',
                  ds['std'], 'per mil', '495')
        _add_stat('F. Case 1', 'F.4_case1_dissolved_winter_model_std',
                  'Case 1 modeled dissolved methane isotope standard deviation in winter',
                  s_w, 'per mil', '495')
        _add_stat('F. Case 1', 'F.4_case1_dissolved_winter_obs_std',
                  'Observed dissolved methane isotope standard deviation in winter',
                  dw['std'], 'per mil', '495')
        _add_stat('F. Case 1', 'F.5_case1_sediment_production_mean',
                  'Case 1 mean sediment methane production isotope value',
                  sed_mean, 'per mil VPDB', '495')
        print("="*70)

    # ----- G. Case 2 (Optimized static) -----
    case2 = condition_results.get('All Observations')
    if case2 is not None:
        ts = case2['timeseries']
        p = case2['params']
        d = ts['delta_dissolved']
        m_s, s_s = _mean_std(d[summer_mask])
        m_w, s_w = _mean_std(d[winter_mask])
        print("\n" + "="*70)
        print("PAPER STATISTICS — G. Case 2 (Optimized Fractionation, Static)")
        print("="*70)
        print("G.1 Case 2 dissolved δ¹³C model vs obs, summer: model %.2f ‰  obs %.2f ‰" % (m_s, ds['mean']))
        print("G.2 Case 2 dissolved δ¹³C model vs obs, winter: model %.2f ‰  obs %.2f ‰" % (m_w, dw['mean']))
        print("G.3 Case 2 dissolved variability, summer: model std %.2f  obs std %.2f" % (s_s, ds['std']))
        print("G.4 Case 2 dissolved variability, winter: model std %.2f  obs std %.2f" % (s_w, dw['std']))
        alpha_mo = p.get('alpha_mo', np.nan)
        print("G.5 Case 2 inferred oxidation fractionation α_MO: %.4f  (ε_MO = %.2f ‰)" % (alpha_mo, (alpha_mo - 1) * 1000))
        alpha_op = p.get('alpha_op', np.nan)
        print("G.6 Case 2 oxic production α_OP (at bounds): [1.000, 1.080]  optimized: %.4f  (ε_OP = %.2f ‰)" % (alpha_op, (alpha_op - 1) * 1000))
        sed_prod = ts.get('C13_sed_prod', ts.get('C13_sed_prod_mean', np.nan))
        sed_mean = np.nanmean(sed_prod) if np.size(sed_prod) else np.nan
        print("G.7 Case 2 sediment production δ¹³C-CH₄,prod: %.2f ‰" % sed_mean)
        _add_stat('G. Case 2', 'G.1_case2_dissolved_summer_model_mean',
                  'Case 2 modeled dissolved methane isotope mean in summer',
                  m_s, 'per mil VPDB', '503')
        _add_stat('G. Case 2', 'G.1_case2_dissolved_summer_obs_mean',
                  'Observed dissolved methane isotope mean in summer used for Case 2 comparison',
                  ds['mean'], 'per mil VPDB', '503')
        _add_stat('G. Case 2', 'G.2_case2_dissolved_winter_model_mean',
                  'Case 2 modeled dissolved methane isotope mean in winter',
                  m_w, 'per mil VPDB', '503')
        _add_stat('G. Case 2', 'G.2_case2_dissolved_winter_obs_mean',
                  'Observed dissolved methane isotope mean in winter used for Case 2 comparison',
                  dw['mean'], 'per mil VPDB', '503')
        _add_stat('G. Case 2', 'G.3_case2_dissolved_summer_model_std',
                  'Case 2 modeled dissolved methane isotope standard deviation in summer',
                  s_s, 'per mil', '503')
        _add_stat('G. Case 2', 'G.3_case2_dissolved_summer_obs_std',
                  'Observed dissolved methane isotope standard deviation in summer',
                  ds['std'], 'per mil', '503')
        _add_stat('G. Case 2', 'G.4_case2_dissolved_winter_model_std',
                  'Case 2 modeled dissolved methane isotope standard deviation in winter',
                  s_w, 'per mil', '503')
        _add_stat('G. Case 2', 'G.4_case2_dissolved_winter_obs_std',
                  'Observed dissolved methane isotope standard deviation in winter',
                  dw['std'], 'per mil', '503')
        _add_stat('G. Case 2', 'G.5_case2_alpha_mo',
                  'Case 2 optimized methane oxidation fractionation factor',
                  alpha_mo, 'unitless alpha', '503')
        _add_stat('G. Case 2', 'G.5_case2_epsilon_mo',
                  'Case 2 optimized methane oxidation enrichment factor',
                  (alpha_mo - 1) * 1000, 'per mil', '503')
        _add_stat('G. Case 2', 'G.6_case2_alpha_op',
                  'Case 2 optimized oxic production fractionation factor',
                  alpha_op, 'unitless alpha', '503')
        _add_stat('G. Case 2', 'G.6_case2_epsilon_op',
                  'Case 2 optimized oxic production enrichment factor',
                  (alpha_op - 1) * 1000, 'per mil', '503;567')
        _add_stat('G. Case 2', 'G.7_case2_sediment_production_mean',
                  'Case 2 mean sediment methane production isotope value',
                  sed_mean, 'per mil VPDB', '503')
        print("="*70)

    # ----- H. Case 3 (Temperature-dependent) -----
    case3 = condition_results.get('Temperature-Based')
    if case3 is not None:
        ts = case3['timeseries']
        p = case3['params']
        d = ts['delta_dissolved']
        e = ts['delta_emission']
        m_s, s_s = _mean_std(d[summer_mask])
        m_w, s_w = _mean_std(d[winter_mask])
        me_s, se_s = _mean_std(e[summer_mask])
        me_w, se_w = _mean_std(e[winter_mask])
        print("\n" + "="*70)
        print("PAPER STATISTICS — H. Case 3 (Temperature-Dependent Production)")
        print("="*70)
        print("H.1 Case 3 dissolved δ¹³C model vs obs, summer: model %.2f ‰  obs %.2f ‰" % (m_s, ds['mean']))
        print("H.2 Case 3 dissolved δ¹³C model vs obs, winter: model %.2f ‰  obs %.2f ‰" % (m_w, dw['mean']))
        print("H.3 Case 3 dissolved variability, summer: model std %.2f  obs std %.2f" % (s_s, ds['std']))
        print("H.4 Case 3 dissolved variability, winter: model std %.2f  obs std %.2f" % (s_w, dw['std']))
        print("H.5 Case 3 oxidation fractionation α_MO: %.4f" % p.get('alpha_mo', np.nan))
        sed_prod = ts.get('C13_sed_prod', ts.get('C13_sed_prod_mean', np.nan))
        sed_mean = np.nanmean(sed_prod) if np.size(sed_prod) else np.nan
        print("H.6 Case 3 mean sediment production δ¹³C-CH₄,prod: %.2f ‰" % sed_mean)
        print("H.7 Case 3 emitted δ¹³C model vs obs, summer: model %.2f ‰  obs %.2f ‰" % (me_s, de['mean']))
        print("H.8 Case 3 emitted δ¹³C model vs obs, winter: model %.2f ‰  obs %.2f ‰" % (me_w, ew['mean']))
        print("H.9 Case 3 optimized temperature slope m: %.4f ‰ °C⁻¹" % p.get('m', np.nan))
        print("H.10 Case 3 emitted variability, summer: model std %.2f  obs std %.2f" % (se_s, de['std']))
        print("H.11 Case 3 emitted variability, winter: model std %.2f  obs std %.2f" % (se_w, ew['std']))
        temp_eq = "delta13C_prod = %.4f * T_C + %.2f" % (p.get('m', np.nan), p.get('b', np.nan))
        print("H.12 Case 3 temperature-response fit: %s" % temp_eq)
        _add_stat('H. Case 3', 'H.1_case3_dissolved_summer_model_mean',
                  'Case 3 modeled dissolved methane isotope mean in summer',
                  m_s, 'per mil VPDB', '521')
        _add_stat('H. Case 3', 'H.1_case3_dissolved_summer_obs_mean',
                  'Observed dissolved methane isotope mean in summer used for Case 3 comparison',
                  ds['mean'], 'per mil VPDB', '521')
        _add_stat('H. Case 3', 'H.2_case3_dissolved_winter_model_mean',
                  'Case 3 modeled dissolved methane isotope mean in winter',
                  m_w, 'per mil VPDB', '521')
        _add_stat('H. Case 3', 'H.2_case3_dissolved_winter_obs_mean',
                  'Observed dissolved methane isotope mean in winter used for Case 3 comparison',
                  dw['mean'], 'per mil VPDB', '521')
        _add_stat('H. Case 3', 'H.3_case3_dissolved_summer_model_std',
                  'Case 3 modeled dissolved methane isotope standard deviation in summer',
                  s_s, 'per mil', '521')
        _add_stat('H. Case 3', 'H.3_case3_dissolved_summer_obs_std',
                  'Observed dissolved methane isotope standard deviation in summer',
                  ds['std'], 'per mil', '521')
        _add_stat('H. Case 3', 'H.4_case3_dissolved_winter_model_std',
                  'Case 3 modeled dissolved methane isotope standard deviation in winter',
                  s_w, 'per mil', '521')
        _add_stat('H. Case 3', 'H.4_case3_dissolved_winter_obs_std',
                  'Observed dissolved methane isotope standard deviation in winter',
                  dw['std'], 'per mil', '521')
        _add_stat('H. Case 3', 'H.5_case3_alpha_mo',
                  'Case 3 optimized methane oxidation fractionation factor',
                  p.get('alpha_mo', np.nan), 'unitless alpha', '521')
        _add_stat('H. Case 3', 'H.6_case3_sediment_production_mean',
                  'Case 3 mean sediment methane production isotope value',
                  sed_mean, 'per mil VPDB', '521')
        _add_stat('H. Case 3', 'H.7_case3_emitted_summer_model_mean',
                  'Case 3 modeled surface-emitted methane isotope mean in summer',
                  me_s, 'per mil VPDB', '521')
        _add_stat('H. Case 3', 'H.7_case3_emitted_summer_obs_mean',
                  'Observed surface-emitted methane isotope mean in summer used for Case 3 comparison',
                  de['mean'], 'per mil VPDB', '521')
        _add_stat('H. Case 3', 'H.8_case3_emitted_winter_model_mean',
                  'Case 3 modeled surface-emitted methane isotope mean in winter',
                  me_w, 'per mil VPDB', '521')
        _add_stat('H. Case 3', 'H.8_case3_emitted_winter_obs_mean',
                  'Observed surface-emitted methane isotope mean in winter used for Case 3 comparison',
                  ew['mean'], 'per mil VPDB', '521')
        _add_stat('H. Case 3', 'H.9_case3_temperature_slope',
                  'Case 3 optimized sediment production isotope temperature-response slope',
                  p.get('m', np.nan), 'per mil deg C-1', '521;563;577')
        _add_stat('H. Case 3', 'H.10_case3_emitted_summer_model_std',
                  'Case 3 modeled surface-emitted methane isotope standard deviation in summer',
                  se_s, 'per mil', '521')
        _add_stat('H. Case 3', 'H.10_case3_emitted_summer_obs_std',
                  'Observed surface-emitted methane isotope standard deviation in summer',
                  de['std'], 'per mil', '521')
        _add_stat('H. Case 3', 'H.11_case3_emitted_winter_model_std',
                  'Case 3 modeled surface-emitted methane isotope standard deviation in winter',
                  se_w, 'per mil', '521')
        _add_stat('H. Case 3', 'H.11_case3_emitted_winter_obs_std',
                  'Observed surface-emitted methane isotope standard deviation in winter',
                  ew['std'], 'per mil', '521')
        _add_stat('H. Case 3', 'H.12_case3_temperature_response_equation',
                  'Case 3 optimized linear temperature-response equation',
                  temp_eq, 'per mil VPDB', '521',
                  'T_C is flux-weighted mean sediment temperature in deg C.')
        print("="*70)

    # ----- I. Discussion -----
    print("\n" + "="*70)
    print("PAPER STATISTICS — I. Discussion")
    print("="*70)
    emitted_seasonal_difference = de['mean'] - ew['mean']
    summer_mean_model_flux = float(np.nanmean(albm_data.total_flux[summer_mask])) if np.any(summer_mask) else np.nan
    winter_mean_model_flux = float(np.nanmean(albm_data.total_flux[winter_mask])) if np.any(winter_mask) else np.nan
    winter_flux_pct_of_summer = (
        winter_mean_model_flux / summer_mean_model_flux * 100.0
        if summer_mean_model_flux and summer_mean_model_flux > 0 else np.nan
    )
    print("I.1 Summer–winter difference in emitted δ¹³C (‰): %.2f" % emitted_seasonal_difference)
    print("I.1b Winter flux as percent of summer flux in isoALBM: %.1f%%" % winter_flux_pct_of_summer)
    _add_stat('I. Discussion', 'I.1_emitted_summer_winter_difference',
              'Observed summer minus winter difference in surface-emitted methane isotope means',
              emitted_seasonal_difference, 'per mil', '552')
    _add_stat('I. Discussion', 'I.1b_winter_flux_pct_of_summer',
              'Winter mean isoALBM total methane flux as percent of summer mean isoALBM total methane flux',
              winter_flux_pct_of_summer, 'percent', '565')
    case3 = condition_results.get('Temperature-Based')
    if case3 is not None:
        m_slope = case3['params'].get('m', np.nan)
        b_intercept = case3['params'].get('b', np.nan)
        print("I.2 Temperature sensitivity: Case 3 slope m = %.4f ‰ °C⁻¹ (Rozmiarek et al. 2025: report value for comparison)" % m_slope)
        print("I.3 Case 3 optimized intercept b (δ¹³C at T=0°C): %.2f ‰" % b_intercept)
        alpha_op = case3['params'].get('alpha_op', np.nan)
        print("I.4 Inferred oxic production α_OP (Case 3): %.4f  (ε_OP = %.2f ‰)" % (alpha_op, (alpha_op - 1) * 1000))
        _add_stat('I. Discussion', 'I.2_case3_temperature_slope',
                  'Case 3 optimized sediment production isotope temperature-response slope',
                  m_slope, 'per mil deg C-1', '563;577',
                  'Rozmiarek et al. comparison value must be added from the cited paper.')
        _add_stat('I. Discussion', 'I.3_case3_temperature_intercept',
                  'Case 3 optimized sediment production isotope temperature-response intercept',
                  b_intercept, 'per mil VPDB', '563')
        _add_stat('I. Discussion', 'I.4_case3_alpha_op',
                  'Case 3 optimized oxic production fractionation factor',
                  alpha_op, 'unitless alpha', '567')
        _add_stat('I. Discussion', 'I.4_case3_epsilon_op',
                  'Case 3 optimized oxic production enrichment factor',
                  (alpha_op - 1) * 1000, 'per mil', '567')
        ts = case3['timeseries']
        if 'sediment_temp_avg' in flux_data and 'C13_sed_prod' in ts:
            T_sed = flux_data['sediment_temp_avg']
            T_summer = T_sed[summer_mask]
            prod_summer = ts['C13_sed_prod'][summer_mask]
            emission_summer = ts.get('delta_emission', np.full_like(T_sed, np.nan))[summer_mask]
            if np.sum(summer_mask) > 0:
                print("I.5 Summer sediment temp range: %.2f – %.2f °C; δ¹³C-CH₄,prod range: %.2f – %.2f ‰" % (
                    float(np.nanmin(T_summer)), float(np.nanmax(T_summer)),
                    float(np.nanmin(prod_summer)), float(np.nanmax(prod_summer))))
                print("I.6 Summer emitted δ¹³C-CH₄ range: %.2f – %.2f ‰" % (
                    float(np.nanmin(emission_summer)), float(np.nanmax(emission_summer))))
                _add_stat('I. Discussion', 'I.5_summer_sediment_temp_min',
                          'Minimum summer flux-weighted mean sediment temperature',
                          float(np.nanmin(T_summer)), 'deg C', '565')
                _add_stat('I. Discussion', 'I.5_summer_sediment_temp_max',
                          'Maximum summer flux-weighted mean sediment temperature',
                          float(np.nanmax(T_summer)), 'deg C', '565')
                _add_stat('I. Discussion', 'I.5_summer_delta13c_prod_min',
                          'Minimum summer Case 3 sediment production isotope value',
                          float(np.nanmin(prod_summer)), 'per mil VPDB', '565')
                _add_stat('I. Discussion', 'I.5_summer_delta13c_prod_max',
                          'Maximum summer Case 3 sediment production isotope value',
                          float(np.nanmax(prod_summer)), 'per mil VPDB', '565')
                _add_stat('I. Discussion', 'I.6_summer_delta13c_emission_min',
                          'Minimum summer Case 3 surface-emitted isotope value',
                          float(np.nanmin(emission_summer)), 'per mil VPDB', '565')
                _add_stat('I. Discussion', 'I.6_summer_delta13c_emission_max',
                          'Maximum summer Case 3 surface-emitted isotope value',
                          float(np.nanmax(emission_summer)), 'per mil VPDB', '565')
            else:
                print("I.5 Summer sediment temp / δ¹³C,prod range: N/A")
        else:
            print("I.5 Summer sediment temp and δ¹³C,prod range: N/A (no temp-based timeseries)")
    else:
        print("I.2–I.5 (Case 3 not available)")
    print("="*70)

    # ----- J. Conclusion -----
    print("\n" + "="*70)
    print("PAPER STATISTICS — J. Conclusion (Confirmatory)")
    print("="*70)
    print("J.1 Ebullition/diffusion ratio: %.2f (same as D.10)" % ebull_diff_ratio)
    _add_stat('J. Conclusion', 'J.1_ebullition_diffusion_ratio',
              'Confirmatory ALBM summer ebullition-to-diffusion ratio',
              ebull_diff_ratio, 'ratio', '575')
    if case3 is not None:
        temp_slope = case3['params'].get('m', np.nan)
        print("J.2 Temperature response slope: %.4f ‰ °C⁻¹ (same as H.9 / I.2)" % temp_slope)
        _add_stat('J. Conclusion', 'J.2_temperature_response_slope',
                  'Confirmatory Case 3 temperature-response slope',
                  temp_slope, 'per mil deg C-1', '577')
    else:
        print("J.2 Temperature response slope: N/A")
        _add_stat('J. Conclusion', 'J.2_temperature_response_slope',
                  'Confirmatory Case 3 temperature-response slope',
                  np.nan, 'per mil deg C-1', '577')
    print("="*70 + "\n")

    stats_df = pd.DataFrame(rows)
    if csv_path is not None:
        out_dir = os.path.dirname(csv_path)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
        stats_df.to_csv(csv_path, index=False, encoding='utf-8-sig')
        print(f"Saved paper statistics CSV: {csv_path}")
    return stats_df


# =============================================================================
# Main execution
# =============================================================================
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    from data_loader import get_processed_data, DEFAULT_CACHE_PATH, DEFAULT_EDDY_FLUX_FILE

    # Load data (from cache if present, else from raw and save cache)
    albm_data, eddy_data = get_processed_data(
        cache_path=DEFAULT_CACHE_PATH,
        results_dir='ALBM data',
        date_range='20210101_20250101',
        eddy_flux_file=DEFAULT_EDDY_FLUX_FILE,
        use_cache=True,
        save_cache=False
    )
    
    # Generate all plots
    plot_all(albm_data, eddy_data, figs_dir='figs', show=True)
