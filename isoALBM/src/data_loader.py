"""
ALBM Data Loader Module

Loads ALBM model output from NetCDF files and eddy covariance flux data.
Creates standardized data structures for use by analysis modules.

Can save/load processed data to a single cache file so other machines (e.g. after
cloning from GitHub) can run analysis without the raw NetCDF/eddy files.

Usage:
    from data_loader import load_albm_data, load_eddy_flux_data
    
    albm_data = load_albm_data(results_dir='ALBM data', date_range='20210101_20250101')
    eddy_data = load_eddy_flux_data(filepath='data/2026_05_08_EC_CH4.xlsx')

    # Prefer: use cache when available (for switching machines)
    from data_loader import get_processed_data
    albm_data, eddy_data = get_processed_data()
    flux_data = get_flux_data_dict(albm_data)
"""

import os
import sys
import time
import pickle
import argparse
import numpy as np
import pandas as pd
import netCDF4 as nc
from dataclasses import dataclass, fields
from typing import Optional, Dict, Any, Tuple, Sequence

for _stream in (sys.stdout, sys.stderr):
    if hasattr(_stream, "reconfigure"):
        _stream.reconfigure(encoding="utf-8")


def _running_in_interactive_kernel() -> bool:
    """Return True when code is being executed by IPython/Jupyter."""
    return hasattr(sys, "ps1") or "ipykernel" in sys.modules


def _candidate_roots():
    """Yield likely repository roots for script, import, and interactive use."""
    seen = set()
    script_file = globals().get("__file__", "")
    starts = []
    if script_file:
        starts.append(os.path.dirname(os.path.abspath(script_file)))
    starts.append(os.path.abspath(os.getcwd()))
    for start in starts:
        current = os.path.abspath(start)
        while current and current not in seen:
            seen.add(current)
            yield current
            parent = os.path.dirname(current)
            if parent == current:
                break
            current = parent


def _find_project_root() -> str:
    for root in _candidate_roots():
        if (
            os.path.isfile(os.path.join(root, "src", "data_loader.py"))
            and os.path.isdir(os.path.join(root, "ALBM data"))
        ):
            return root
    return os.path.abspath(os.getcwd())


PROJECT_ROOT = _find_project_root()


def _repo_path(path: str) -> str:
    """Resolve relative project paths from the repository root."""
    expanded = os.path.expanduser(os.path.expandvars(str(path)))
    if os.path.isabs(expanded):
        return os.path.normpath(expanded)
    cwd_path = os.path.abspath(expanded)
    if os.path.exists(cwd_path):
        return os.path.normpath(cwd_path)
    return os.path.normpath(os.path.join(PROJECT_ROOT, expanded))


# Default path for cached processed data (can be committed to GitHub)
DEFAULT_CACHE_PATH = _repo_path('data/processed_albm_data.pkl')
DEFAULT_ALBM_RESULTS_DIR = _repo_path('ALBM data')
DEFAULT_EDDY_FLUX_FILE = _repo_path('data/2026_05_08_EC_CH4.xlsx')
LEGACY_EDDY_FLUX_FILE = _repo_path('data/summary_report_2021-01-01_2025-01-27.txt')


@dataclass
class ALBMData:
    """Container for ALBM model output data."""
    # Time information
    time_array: np.ndarray
    start_date: np.datetime64
    end_date: np.datetime64
    n_timesteps: int
    
    # Raw data arrays (in model units: mol/m³ or mol/m³/s)
    ch4eb: np.ndarray        # Surface ebullition flux
    ch4df: np.ndarray        # Surface diffusion flux
    ch4conc: np.ndarray      # Depth-resolved CH4 concentration
    ch4oxid: np.ndarray      # Depth-resolved CH4 oxidation
    sedch4df: np.ndarray     # Sediment CH4 diffusion
    sedch4eb: np.ndarray     # Sediment CH4 ebullition
    och4prod: np.ndarray     # Oxic CH4 production
    ch4exc: np.ndarray       # CH4 exchange
    ch4upb: np.ndarray       # Upward bubbling/ice leakage
    icebch4: np.ndarray      # Ice-trapped bubble CH4
    icebleak: np.ndarray     # Ice bubble leakage
    Az: np.ndarray           # Lake area by depth
    
    # Processed flux data (in kg/m²/s or kg/m³)
    eflux_data: np.ndarray       # Ebullition flux (kg m⁻² s⁻¹)
    dflux_data: np.ndarray       # Diffusion flux (kg m⁻² s⁻¹)
    och4_data: np.ndarray        # Oxidation rate (kg m⁻³ s⁻¹)
    ch4prod_data: np.ndarray     # Sediment CH4 production (kg m⁻³ s⁻¹)
    och4prod_data: np.ndarray    # Oxic CH4 production (kg m⁻³ s⁻¹)
    sedch4df_data: np.ndarray    # Mean sediment diffusion (kg m⁻³ s⁻¹)
    sedch4eb_data: np.ndarray    # Mean sediment ebullition (kg m⁻³ s⁻¹)
    ch4exc_data: np.ndarray      # CH4 exchange (kg m⁻³ s⁻¹)
    ch4conc_data: np.ndarray     # CH4 concentration (kg m⁻³)
    ch4upb_data: np.ndarray      # Upward bubbling (kg m⁻² s⁻¹)
    exbch4_data: np.ndarray      # Ice bubble exchange (kg m⁻³ s⁻¹)
    icebleak_data: np.ndarray    # Ice leakage (kg m⁻³ s⁻¹)
    cch4_data: np.ndarray        # CH4 concentration change rate
    
    # Depth arrays
    water_depth_array: np.ndarray
    sed_depth_array: np.ndarray
    n_water_layers: int
    n_sed_layers: int
    
    # Total flux for comparison
    total_flux: np.ndarray  # Total CH4 flux (kg m⁻² s⁻¹)
    
    # Sediment temperature data (from isobLakeOut.sedtemp.*.nc, dims Lake/Time/Z/COL)
    # sediment_temp_depth: flux-weighted mean over COL=5 sub-columns → (n_timesteps, Z=41) in °C
    #   Depths: 0-25 m with 0.625 m spacing (hsed=25 m, Z=41 nodes).
    # sediment_temp_avg: mean over Z=41 depth layers → (n_timesteps,) in °C
    sediment_temp_depth: Optional[np.ndarray]  # (n_timesteps, Z=41) depth-resolved sediment temperature (°C)
    sediment_temp_avg: Optional[np.ndarray]    # (n_timesteps,) depth-mean sediment temperature (°C)


@dataclass
class EddyFluxData:
    """Container for eddy covariance flux data."""
    datetime: pd.DatetimeIndex
    ch4_flux: np.ndarray           # CH4 flux (kg m⁻² s⁻¹)
    ch4_flux_daily_avg: pd.DataFrame  # Daily averaged flux data
    

if __name__ == "__main__":
    sys.modules.setdefault("data_loader", sys.modules[__name__])
    ALBMData.__module__ = "data_loader"
    EddyFluxData.__module__ = "data_loader"


class _ProcessedDataUnpickler(pickle.Unpickler):
    """Read caches written when data_loader.py was run as __main__."""

    def find_class(self, module: str, name: str) -> Any:
        if module in {"__main__", "data_loader"} and name in {"ALBMData", "EddyFluxData"}:
            return globals()[name]
        return super().find_class(module, name)


def _dataclass_payload(instance: Any) -> Optional[Dict[str, Any]]:
    if instance is None:
        return None
    return {field.name: getattr(instance, field.name) for field in fields(instance)}


def load_netcdf_variable(results_dir: str, filename: str, variable_name: str) -> np.ndarray:
    """Load a variable from a NetCDF file."""
    filepath = os.path.join(results_dir, filename)
    if not os.path.exists(filepath):
        raise FileNotFoundError(
            f"Required ALBM NetCDF file not found: {filepath}. "
            "Check results_dir or run from the repository root."
        )
    with nc.Dataset(filepath, 'r') as ds:
        if variable_name not in ds.variables:
            raise KeyError(f"Variable {variable_name!r} not found in {filepath}")
        return ds.variables[variable_name][:]


def _sediment_temp_weighted_mean(
    sediment_temp_depth: np.ndarray,
    sedch4df: np.ndarray,
    sedch4eb: np.ndarray
) -> np.ndarray:
    """
    Mean sediment temperature weighted by methane output (diffusion + ebullition)
    per depth, so depths with more CH4 flux have higher weight (temperature
    experienced by methane).
    """
    T, D = sediment_temp_depth.shape[0], sediment_temp_depth.shape[1]
    # Flux arrays may be (T,) or (T, D); ensure 2D for per-depth weights
    if sedch4df.ndim == 1:
        sedch4df = np.broadcast_to(sedch4df[:, np.newaxis], (T, D))
    elif sedch4df.shape[0] != T or sedch4df.shape[1] != D:
        return np.mean(sediment_temp_depth, axis=1)
    if sedch4eb.ndim == 1:
        sedch4eb = np.broadcast_to(sedch4eb[:, np.newaxis], (T, D))
    elif sedch4eb.shape[0] != T or sedch4eb.shape[1] != D:
        return np.mean(sediment_temp_depth, axis=1)
    weights = np.maximum(sedch4df + sedch4eb, 0.0)
    out = np.zeros(T)
    for i in range(T):
        w = weights[i, :]
        s = w.sum()
        if s > 0:
            out[i] = np.average(sediment_temp_depth[i, :], weights=w)
        else:
            out[i] = np.mean(sediment_temp_depth[i, :])
    return out


def load_albm_data(
    results_dir: str = DEFAULT_ALBM_RESULTS_DIR,
    date_range: str = '20210101_20250101',
    lake_area: float = 40000.0
) -> ALBMData:
    """
    Load ALBM model output from NetCDF files.
    
    Parameters:
    -----------
    results_dir : str
        Directory containing ALBM NetCDF output files
    date_range : str
        Date range string in format 'YYYYMMDD_YYYYMMDD'
    lake_area : float
        Lake surface area in m² (default: 40000 for BTL)
    
    Returns:
    --------
    ALBMData : Dataclass containing all ALBM data arrays
    """
    start_time = time.time()
    results_dir = _repo_path(results_dir)
    print(f"[{time.strftime('%H:%M:%S')}] Loading ALBM data from {results_dir}...")
    
    # Molecular weight of CH4 for unit conversion (mol to kg)
    MW_CH4 = 16.04e-3  # kg/mol
    
    # Load raw data from NetCDF files
    Az = load_netcdf_variable(results_dir, f'isobLakeOut.Az.{date_range}.nc', 'Az')
    ch4oxid = load_netcdf_variable(results_dir, f'isobLakeOut.ch4oxid.{date_range}.nc', 'ch4oxid')
    ch4df = load_netcdf_variable(results_dir, f'isobLakeOut.ch4df.{date_range}.nc', 'ch4df')
    sedch4df = load_netcdf_variable(results_dir, f'isobLakeOut.sedch4df.{date_range}.nc', 'sedch4df')
    ch4conc = load_netcdf_variable(results_dir, f'isobLakeOut.ch4conc.{date_range}.nc', 'ch4conc')
    sedch4eb = load_netcdf_variable(results_dir, f'isobLakeOut.sedch4eb.{date_range}.nc', 'sedch4eb')
    och4prod = load_netcdf_variable(results_dir, f'isobLakeOut.och4prod.{date_range}.nc', 'och4prod')
    ch4exc = load_netcdf_variable(results_dir, f'isobLakeOut.ch4exc.{date_range}.nc', 'ch4exc')
    ch4eb = load_netcdf_variable(results_dir, f'isobLakeOut.ch4eb.{date_range}.nc', 'ch4eb')
    ch4exb = load_netcdf_variable(results_dir, f'isobLakeOut.exbch4.{date_range}.nc', 'exbch4')
    icebch4 = load_netcdf_variable(results_dir, f'isobLakeOut.icebch4.{date_range}.nc', 'icebch4')
    icebleak = load_netcdf_variable(results_dir, f'isobLakeOut.icebleak.{date_range}.nc', 'icebleak')
    ch4upb = load_netcdf_variable(results_dir, f'isobLakeOut.ch4upb.{date_range}.nc', 'ch4upb')
    sedtemp = load_netcdf_variable(results_dir, f'isobLakeOut.sedtemp.{date_range}.nc', 'sedtemp')
    
    print(f"[{time.strftime('%H:%M:%S')}] Raw data loaded successfully ({time.time() - start_time:.1f}s)")
    
    # Remove singleton dimensions
    Az = np.squeeze(Az)
    ch4oxid = np.squeeze(ch4oxid)
    ch4df = np.squeeze(ch4df)
    sedch4df = np.squeeze(sedch4df)
    ch4conc = np.squeeze(ch4conc)
    sedch4eb = np.squeeze(sedch4eb)
    och4prod = np.squeeze(och4prod)
    ch4exc = np.squeeze(ch4exc)
    ch4eb = np.squeeze(ch4eb)
    ch4exb = np.squeeze(ch4exb)
    icebch4 = np.squeeze(icebch4)
    icebleak = np.squeeze(icebleak)
    ch4upb = np.squeeze(ch4upb)
    
    # Calculate derived quantities
    # Oxidized CH4: sum over depth with layer thickness (0.1 m) and area weighting
    och4 = np.sum(ch4oxid.T * 0.1 * Az[:, np.newaxis], axis=0) / lake_area
    
    # Mean sediment CH4 production (ebullition + diffusion)
    ch4p = np.mean(sedch4eb + sedch4df, axis=1)
    
    # Oxidized CH4 production: sum over depth with layer thickness and area weighting
    och4p = np.sum(och4prod.T * 0.1 * Az[:, np.newaxis], axis=0) / lake_area
    
    # CH4 exchange: sum over depth with layer thickness and area weighting
    exch4 = np.sum(ch4exc.T * 0.1 * Az[:, np.newaxis], axis=0) / lake_area
    
    # Ice-trapped bubble CH4 exchange
    exbch4 = ch4exb / lake_area
    
    # CH4 concentration: sum over depth with layer thickness and area weighting
    ch4con = np.sum(ch4conc.T * 0.1 * Az[:, np.newaxis], axis=0) / lake_area
    
    # CH4 concentration change rate using finite difference
    ch4con_padded_start = np.concatenate([[0], ch4con])
    ch4con_padded_end = np.concatenate([ch4con, [0]])
    cch4dif = (ch4con_padded_start - ch4con_padded_end) / 86400
    cch4 = cch4dif[:-1]
    
    # Ice bubble leakage rate using finite difference
    icebch4_padded_start = np.concatenate([[0], icebch4])
    icebch4_padded_end = np.concatenate([icebch4, [0]])
    icebleakdif = icebch4_padded_start - icebch4_padded_end
    leak = icebleakdif[:-1]
    
    # Set ch4upb to 0 when surface ebullition (ch4eb) isn't 0
    ch4upb_processed = ch4upb.copy()
    ch4upb_processed[np.abs(ch4eb) > 0] = 0
    
    # Convert to kg units for flux calculations
    eflux_data = ch4eb * MW_CH4
    dflux_data = ch4df * MW_CH4
    och4_data = och4 * MW_CH4
    ch4prod_data = ch4p * MW_CH4
    och4prod_data = och4p * MW_CH4
    sedch4df_data = np.mean(sedch4df, axis=1) * MW_CH4
    sedch4eb_data = np.mean(sedch4eb, axis=1) * MW_CH4
    ch4exc_data = exch4 * MW_CH4
    ch4conc_data = ch4con * MW_CH4
    ch4upb_data = ch4upb_processed * MW_CH4
    exbch4_data = exbch4 * MW_CH4
    icebleak_data = icebleak * MW_CH4
    cch4_data = cch4 * MW_CH4
    
    # Total flux
    total_flux = dflux_data + eflux_data + ch4upb_data
    
    # Time series length and calendar
    data_length = len(ch4eb) if len(ch4eb) > 0 else 1461
    start_date_str = date_range.split('_')[0]
    start_date = np.datetime64(f'{start_date_str[:4]}-{start_date_str[4:6]}-{start_date_str[6:8]}')
    
    # Sediment temperature: load from NetCDF.
    # File dims: (Lake=1, Time, Z=41, COL=5). After squeeze: (Time, Z, COL).
    # Processing:
    #   1. Average over COL weighted by per-column (sedch4df + sedch4eb) flux.
    #   2. Convert Kelvin -> Celsius (model outputs in K).
    #   3. sediment_temp_depth: (Time, Z=41); sediment_temp_avg: mean over Z.
    KELVIN_OFFSET = 273.15
    sediment_temp_depth = None
    sediment_temp_avg = None
    sediment_temp_error = None
    sedtemp = np.squeeze(sedtemp) if len(sedtemp) > 0 else np.array([])
    if len(sedtemp) > 0:
        sedtemp = sedtemp.astype(float)
        if sedtemp.ndim == 3:
            # Expected layout after squeeze: (Time, Z, COL)
            T_dim, Z_dim, COL_dim = sedtemp.shape
            if T_dim == data_length:
                # Flux-weighted average over COL (per-column sedch4df + sedch4eb weights)
                if (sedch4df.ndim == 2 and sedch4df.shape == (T_dim, COL_dim) and
                        sedch4eb.ndim == 2 and sedch4eb.shape == (T_dim, COL_dim)):
                    col_weights = np.maximum(sedch4df + sedch4eb, 0.0)       # (T, COL)
                    col_weight_sum = col_weights.sum(axis=1)                  # (T,)
                    col_weights_bc = col_weights[:, np.newaxis, :]            # (T, 1, COL)
                    weighted_col_sum = (sedtemp * col_weights_bc).sum(axis=2) # (T, Z)
                    w_sum_bc = np.maximum(col_weight_sum[:, np.newaxis], 1e-300)
                    no_flux = (col_weight_sum == 0)[:, np.newaxis]
                    col_mean = sedtemp.mean(axis=2)                           # (T, Z) zero-flux substitute
                    temp_tz = np.where(no_flux, col_mean, weighted_col_sum / w_sum_bc)
                    weighting_desc = f"flux-weighted over COL={COL_dim}"
                else:
                    temp_tz = sedtemp.mean(axis=2)                            # (T, Z)
                    weighting_desc = f"uniform mean over COL={COL_dim}"
                # Convert K -> °C (model outputs in Kelvin)
                if temp_tz.mean() > 200:
                    temp_tz -= KELVIN_OFFSET
                sediment_temp_depth = temp_tz                                 # (T, Z=41)
                sediment_temp_avg = sediment_temp_depth.mean(axis=1)          # (T,) mean over Z
                hsed_m = 0.625 * (Z_dim - 1)  # 25 m for Z=41 with 0.625 m spacing
                print(f"[{time.strftime('%H:%M:%S')}] Sediment temperature: 3D (Time={T_dim}, "
                      f"Z={Z_dim}, COL={COL_dim}) -> {weighting_desc} -> ({T_dim}, {Z_dim}); "
                      f"depth 0-{hsed_m:.2f} m; "
                      f"mean range {sediment_temp_avg.min():.1f}-{sediment_temp_avg.max():.1f}°C")
            else:
                sediment_temp_error = (
                    f"sedtemp 3D shape {sedtemp.shape} has time axis {T_dim}, "
                    f"expected {data_length}"
                )
        elif sedtemp.ndim == 2:
            # (Time, Z) or (Z, Time)
            if sedtemp.shape[0] == data_length:
                sediment_temp_depth = sedtemp
            elif sedtemp.shape[1] == data_length:
                sediment_temp_depth = sedtemp.T
            if sediment_temp_depth is not None:
                if sediment_temp_depth.mean() > 200:
                    sediment_temp_depth = sediment_temp_depth - KELVIN_OFFSET
                sediment_temp_avg = _sediment_temp_weighted_mean(
                    sediment_temp_depth, sedch4df, sedch4eb
                )
            else:
                sediment_temp_error = (
                    f"sedtemp 2D shape {sedtemp.shape} does not include "
                    f"the expected time length {data_length}"
                )
        elif sedtemp.ndim == 1 and sedtemp.shape[0] == data_length:
            sediment_temp_avg = sedtemp.copy()
            if sediment_temp_avg.mean() > 200:
                sediment_temp_avg -= KELVIN_OFFSET
            sediment_temp_depth = np.atleast_2d(sediment_temp_avg).T
        else:
            sediment_temp_error = (
                f"unhandled sedtemp shape {sedtemp.shape}; expected "
                f"(Time, Z, COL), (Time, Z), (Z, Time), or (Time,)"
            )
    else:
        sediment_temp_error = (
            f"missing or empty sedtemp variable in "
            f"isobLakeOut.sedtemp.{date_range}.nc"
        )
    if sediment_temp_avg is None:
        detail = sediment_temp_error or "sediment temperature could not be parsed"
        print(f"[{time.strftime('%H:%M:%S')}] ERROR: Sediment temperature unavailable: {detail}. "
              "No synthetic sediment temperature will be generated.")
    
    # Create time array (start_date already computed above for day-of-year)
    time_array = np.arange(start_date, start_date + np.timedelta64(data_length, 'D'), np.timedelta64(1, 'D'))
    end_date = time_array[-1]
    
    # Define depth arrays.
    # sed_depth_array reflects the sediment CH4 flux grid (COL columns, no Z depth).
    # sediment_temp_depth has Z=41 depth layers over hsed=25 m (0.625 m spacing).
    n_water_layers = ch4oxid.shape[1] if ch4oxid.ndim > 1 else 1
    n_sed_layers = sedch4df.shape[1] if sedch4df.ndim > 1 else 1
    water_depth_array = np.arange(0, n_water_layers * 0.1, 0.1)
    # sed_depth_array: placeholder using COL count; physical depths are in sediment_temp_depth (Z)
    sed_depth_array = np.arange(0, n_sed_layers * 0.1, 0.1)
    
    elapsed = time.time() - start_time
    print(f"[{time.strftime('%H:%M:%S')}] ALBM data loaded: {data_length} timesteps from {start_date} to {end_date} ({elapsed:.1f}s)")
    
    return ALBMData(
        time_array=time_array,
        start_date=start_date,
        end_date=end_date,
        n_timesteps=data_length,
        ch4eb=ch4eb,
        ch4df=ch4df,
        ch4conc=ch4conc,
        ch4oxid=ch4oxid,
        sedch4df=sedch4df,
        sedch4eb=sedch4eb,
        och4prod=och4prod,
        ch4exc=ch4exc,
        ch4upb=ch4upb_processed,
        icebch4=icebch4,
        icebleak=icebleak,
        Az=Az,
        eflux_data=eflux_data,
        dflux_data=dflux_data,
        och4_data=och4_data,
        ch4prod_data=ch4prod_data,
        och4prod_data=och4prod_data,
        sedch4df_data=sedch4df_data,
        sedch4eb_data=sedch4eb_data,
        ch4exc_data=ch4exc_data,
        ch4conc_data=ch4conc_data,
        ch4upb_data=ch4upb_data,
        exbch4_data=exbch4_data,
        icebleak_data=icebleak_data,
        cch4_data=cch4_data,
        water_depth_array=water_depth_array,
        sed_depth_array=sed_depth_array,
        n_water_layers=n_water_layers,
        n_sed_layers=n_sed_layers,
        total_flux=total_flux,
        sediment_temp_depth=sediment_temp_depth,
        sediment_temp_avg=sediment_temp_avg
    )


def _canonical_path(path: Optional[str]) -> Optional[str]:
    """Return a stable absolute path string for cache source comparisons."""
    if path is None:
        return None
    return os.path.normcase(os.path.abspath(os.fspath(path)))


def _load_summary_report_daily_flux(filepath: str) -> pd.DataFrame:
    """
    Load an EddyPro summary report.

    The report's ch4_flux column is in umol CH4 m-2 s-1. We average the
    half-hourly values by day and convert to kg CH4 m-2 s-1.
    """
    df = pd.read_csv(filepath, sep='\t', skiprows=[1])
    required_cols = {'date', 'time', 'ch4_flux'}
    missing = sorted(required_cols - set(df.columns))
    if missing:
        raise ValueError(f"Missing required summary-report columns: {missing}")

    df['datetime'] = pd.to_datetime(
        df['date'].astype(str) + ' ' + df['time'].astype(str),
        errors='coerce'
    )
    df['ch4_flux'] = pd.to_numeric(df['ch4_flux'], errors='coerce')
    df = df.dropna(subset=['datetime'])

    daily = df.groupby(df['datetime'].dt.date)['ch4_flux'].mean().reset_index()
    daily['ch4_flux_daily_avg_kg_m2_s'] = daily['ch4_flux'] * 16.04e-9
    daily = daily[['datetime', 'ch4_flux_daily_avg_kg_m2_s']]
    daily['datetime'] = pd.to_datetime(daily['datetime'])
    return daily


def _load_gapfilled_excel_daily_flux(filepath: str) -> pd.DataFrame:
    """
    Load the daily MDS gap-filled CH4 flux workbook.

    The workbook stores daily CH4 flux in mg CH4 m-2 day-1. Convert it to
    kg CH4 m-2 s-1 to match ALBM and the existing downstream code.
    """
    df = pd.read_excel(filepath)
    df.columns = [str(col).strip() for col in df.columns]

    date_col = 'Date' if 'Date' in df.columns else None
    if date_col is None:
        date_matches = [col for col in df.columns if col.lower().startswith('date')]
        date_col = date_matches[0] if date_matches else None

    preferred_flux_col = 'mgCH4/m2/day from MDS gapfill'
    flux_col = preferred_flux_col if preferred_flux_col in df.columns else None
    if flux_col is None:
        flux_matches = [
            col for col in df.columns
            if 'mgch4' in col.lower().replace(' ', '') and 'day' in col.lower()
        ]
        flux_col = flux_matches[0] if flux_matches else None

    if date_col is None or flux_col is None:
        raise ValueError(
            "Expected an Excel workbook with Date and mgCH4/m2/day flux columns; "
            f"found columns: {list(df.columns)}"
        )

    daily = pd.DataFrame({
        'datetime': pd.to_datetime(df[date_col], errors='coerce'),
        'mg_ch4_m2_day': pd.to_numeric(df[flux_col], errors='coerce')
    }).dropna(subset=['datetime', 'mg_ch4_m2_day'])

    daily = daily.groupby(daily['datetime'].dt.date)['mg_ch4_m2_day'].mean().reset_index()
    daily['datetime'] = pd.to_datetime(daily['datetime'])
    daily['ch4_flux_daily_avg_kg_m2_s'] = daily['mg_ch4_m2_day'] * 1e-6 / 86400.0
    return daily[['datetime', 'ch4_flux_daily_avg_kg_m2_s']]


def load_eddy_flux_data(filepath: str = DEFAULT_EDDY_FLUX_FILE) -> Optional[EddyFluxData]:
    """
    Load eddy covariance flux data from a supported source file.
    
    Parameters:
    -----------
    filepath : str
        Path to an eddy flux source. Supports the legacy tab-delimited EddyPro
        summary report and the daily MDS gap-filled Excel workbook.
    
    Returns:
    --------
    EddyFluxData : Dataclass containing eddy flux data, or None if file not found
    """
    filepath = _repo_path(filepath)
    if not os.path.exists(filepath):
        print(f"[{time.strftime('%H:%M:%S')}] Warning: Eddy flux file not found: {filepath}")
        return None
    
    start_time = time.time()
    print(f"[{time.strftime('%H:%M:%S')}] Loading eddy flux data from {filepath}...")
    
    ext = os.path.splitext(filepath)[1].lower()
    if ext in {'.xlsx', '.xls'}:
        eddy_flux_daily_avg = _load_gapfilled_excel_daily_flux(filepath)
    else:
        eddy_flux_daily_avg = _load_summary_report_daily_flux(filepath)
    
    eddy_flux_daily_avg = eddy_flux_daily_avg.sort_values('datetime').reset_index(drop=True)

    elapsed = time.time() - start_time
    print(f"[{time.strftime('%H:%M:%S')}] Eddy flux loaded: {len(eddy_flux_daily_avg)} days ({eddy_flux_daily_avg['datetime'].min()} to {eddy_flux_daily_avg['datetime'].max()}) ({elapsed:.1f}s)")
    
    return EddyFluxData(
        datetime=pd.DatetimeIndex(eddy_flux_daily_avg['datetime']),
        ch4_flux=eddy_flux_daily_avg['ch4_flux_daily_avg_kg_m2_s'].values,
        ch4_flux_daily_avg=eddy_flux_daily_avg
    )


def save_processed_data(
    albm_data: ALBMData,
    eddy_data: Optional[EddyFluxData],
    path: str = DEFAULT_CACHE_PATH,
    metadata: Optional[Dict[str, str]] = None
) -> str:
    """
    Save processed ALBM and eddy data to a cache file (e.g. for GitHub / other machines).
    
    Parameters:
    -----------
    albm_data : ALBMData
        Loaded ALBM data
    eddy_data : EddyFluxData or None
        Loaded eddy flux data (may be None if file missing)
    path : str
        Output path for the cache file
    metadata : dict, optional
        Optional metadata (e.g. results_dir, date_range, eddy_file) for reference
    
    Returns:
    --------
    str : Path to the saved file
    """
    path = _repo_path(path)
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    payload = {
        'cache_schema_version': 2,
        'albm_data': _dataclass_payload(albm_data),
        'eddy_data': _dataclass_payload(eddy_data),
        'metadata': metadata or {}
    }
    with open(path, 'wb') as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"[{time.strftime('%H:%M:%S')}] Processed data saved to {path}")
    return path


def load_processed_data(
    path: str = DEFAULT_CACHE_PATH,
    return_metadata: bool = False
) -> Optional[Tuple[Any, ...]]:
    """
    Load processed ALBM and eddy data from a cache file.
    
    Parameters:
    -----------
    path : str
        Path to the cache file
    return_metadata : bool
        If True, include cache metadata as the third returned item
    
    Returns:
    --------
    tuple (albm_data, eddy_data), optionally with metadata, or None if file missing
    Raises RuntimeError if a present cache cannot be loaded.
    """
    path = _repo_path(path)
    if not os.path.exists(path):
        return None
    try:
        with open(path, 'rb') as f:
            payload = _ProcessedDataUnpickler(f).load()
        if payload.get('cache_schema_version') == 2:
            albm_data = ALBMData(**payload['albm_data'])
            eddy_payload = payload.get('eddy_data')
            eddy_data = EddyFluxData(**eddy_payload) if eddy_payload is not None else None
        else:
            albm_data = payload['albm_data']
            eddy_data = payload.get('eddy_data')
        meta = payload.get('metadata', {})
        if meta:
            print(f"[{time.strftime('%H:%M:%S')}] Loaded from cache: {path} (metadata: {meta})")
        else:
            print(f"[{time.strftime('%H:%M:%S')}] Loaded from cache: {path}")
        if return_metadata:
            return albm_data, eddy_data, meta
        return albm_data, eddy_data
    except Exception as e:
        raise RuntimeError(f"Cache exists but could not be loaded: {path}") from e


def get_processed_data(
    cache_path: str = DEFAULT_CACHE_PATH,
    results_dir: str = DEFAULT_ALBM_RESULTS_DIR,
    date_range: str = '20210101_20250101',
    eddy_flux_file: str = DEFAULT_EDDY_FLUX_FILE,
    use_cache: bool = True,
    save_cache: bool = False,
    lake_area: float = 40000.0
) -> Tuple[ALBMData, Optional[EddyFluxData]]:
    """
    Get processed data: load from cache if present and use_cache is True,
    otherwise load from raw sources and optionally save to cache.
    
    Use this in scripts so that after cloning the repo (e.g. on another machine),
    running with the cached file avoids needing the raw NetCDF and eddy CSV.
    
    Parameters:
    -----------
    cache_path : str
        Path to the cache file
    results_dir : str
        ALBM data directory (used only if loading from raw)
    date_range : str
        Date range string (used only if loading from raw)
    eddy_flux_file : str
        Eddy flux source path
    use_cache : bool
        If True and cache_path exists, load from cache
    save_cache : bool
        If True and we load from raw, save to cache_path
    lake_area : float
        Lake area in m² (used only if loading from raw)
    
    Returns:
    --------
    tuple : (albm_data, eddy_data)
    """
    cache_path = _repo_path(cache_path)
    results_dir = _repo_path(results_dir)
    eddy_flux_file = _repo_path(eddy_flux_file)
    if use_cache:
        cached = load_processed_data(cache_path, return_metadata=True)
        if cached is not None:
            albm_data, eddy_data, meta = cached
            cached_eddy_file = meta.get('eddy_flux_file') if isinstance(meta, dict) else None
            _eddy_source_stale = (
                _canonical_path(cached_eddy_file) != _canonical_path(eddy_flux_file)
            )
            # Check if temperature data exists and has the expected Big Trail depth shape.
            temp_depth = getattr(albm_data, 'sediment_temp_depth', None)
            temp_avg = getattr(albm_data, 'sediment_temp_avg', None)
            _temp_stale = (
                temp_depth is None
                or temp_avg is None
                or np.ndim(temp_depth) < 2
                or np.shape(temp_depth)[0] != getattr(albm_data, 'n_timesteps', 0)
                or np.shape(temp_depth)[1] <= 5
                or np.ndim(temp_avg) != 1
                or np.shape(temp_avg)[0] != getattr(albm_data, 'n_timesteps', 0)
            )
            if _temp_stale:
                print(f"[{time.strftime('%H:%M:%S')}] Stale/missing sediment temperature in cache; "
                      f"reloading from NetCDF...")
                # Reload from raw using the same logic as load_albm_data.
                # The easiest way is to call load_albm_data directly and take only the temp fields.
                try:
                    _fresh = load_albm_data(
                        results_dir=results_dir, date_range=date_range,
                        lake_area=lake_area
                    )
                    albm_data.sediment_temp_depth = _fresh.sediment_temp_depth
                    albm_data.sediment_temp_avg = _fresh.sediment_temp_avg
                    if albm_data.sediment_temp_depth is None or albm_data.sediment_temp_avg is None:
                        print(f"[{time.strftime('%H:%M:%S')}] ERROR: Sediment temperature could not be "
                              "reloaded from NetCDF; cached temperature fields were cleared.")
                    else:
                        print(f"[{time.strftime('%H:%M:%S')}] Sediment temperature regenerated: "
                              f"shape {albm_data.sediment_temp_depth.shape}, "
                              f"range {albm_data.sediment_temp_avg.min():.1f}-"
                              f"{albm_data.sediment_temp_avg.max():.1f}°C")
                    if save_cache:
                        save_processed_data(
                            albm_data, eddy_data, path=cache_path,
                            metadata={'results_dir': results_dir, 'date_range': date_range,
                                      'eddy_flux_file': eddy_flux_file}
                        )
                except Exception as _e:
                    albm_data.sediment_temp_depth = None
                    albm_data.sediment_temp_avg = None
                    print(f"[{time.strftime('%H:%M:%S')}] ERROR: could not reload sedtemp ({_e}); "
                          "cached temperature fields were cleared.")
                    if save_cache:
                        save_processed_data(
                            albm_data, eddy_data, path=cache_path,
                            metadata={'results_dir': results_dir, 'date_range': date_range,
                                      'eddy_flux_file': eddy_flux_file}
                        )
            if _eddy_source_stale:
                print(f"[{time.strftime('%H:%M:%S')}] Cached eddy source differs from requested source; "
                      f"reloading eddy flux from {eddy_flux_file}...")
                eddy_data = load_eddy_flux_data(eddy_flux_file)
                if save_cache:
                    save_processed_data(
                        albm_data, eddy_data, path=cache_path,
                        metadata={'results_dir': results_dir, 'date_range': date_range,
                                  'eddy_flux_file': eddy_flux_file}
                    )
            return albm_data, eddy_data

    # Load from raw sources
    albm_data = load_albm_data(results_dir=results_dir, date_range=date_range, lake_area=lake_area)
    eddy_data = load_eddy_flux_data(eddy_flux_file)

    if save_cache:
        save_processed_data(
            albm_data, eddy_data, path=cache_path,
            metadata={'results_dir': results_dir, 'date_range': date_range, 'eddy_flux_file': eddy_flux_file}
        )

    return albm_data, eddy_data


def get_flux_data_dict(albm_data: ALBMData) -> Dict[str, Any]:
    """
    Create flux data dictionary for isotope model functions.
    
    Parameters:
    -----------
    albm_data : ALBMData
        ALBM data container
    
    Returns:
    --------
    dict : Dictionary with flux data arrays for isotope calculations
    """
    return {
        'sedch4df_data': albm_data.sedch4df_data,
        'sedch4eb_data': albm_data.sedch4eb_data,   # sediment ebullition → bubble pool
        'och4prod_data': albm_data.och4prod_data,
        'ch4exc_data': albm_data.ch4exc_data,
        'och4_data': albm_data.och4_data,
        'dflux_data': albm_data.dflux_data,
        'ch4conc_data': albm_data.ch4conc_data,
        'eflux_data': albm_data.eflux_data,
        'ch4upb_data': albm_data.ch4upb_data,
        'time_array': albm_data.time_array,
        'sediment_temp_avg': albm_data.sediment_temp_avg,
        'sediment_temp_depth': albm_data.sediment_temp_depth
    }


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Inspect or refresh the processed ALBM cache.")
    parser.add_argument("--cache-path", default=DEFAULT_CACHE_PATH, help="Processed cache path.")
    parser.add_argument("--results-dir", default=DEFAULT_ALBM_RESULTS_DIR, help="ALBM NetCDF directory.")
    parser.add_argument("--date-range", default="20210101_20250101", help="ALBM filename date range.")
    parser.add_argument("--eddy-flux-file", default=DEFAULT_EDDY_FLUX_FILE, help="Optional eddy flux source.")
    parser.add_argument("--lake-area", type=float, default=40000.0, help="Lake area in m^2.")
    parser.add_argument(
        "--refresh-cache",
        action="store_true",
        help="Reload raw NetCDF/eddy inputs and write the processed cache.",
    )
    parse_argv = [] if argv is None and _running_in_interactive_kernel() else argv
    return parser.parse_args(parse_argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    albm_data, eddy_data = get_processed_data(
        cache_path=args.cache_path,
        results_dir=args.results_dir,
        date_range=args.date_range,
        eddy_flux_file=args.eddy_flux_file,
        use_cache=not args.refresh_cache,
        save_cache=args.refresh_cache,
        lake_area=args.lake_area,
    )
    
    print(f"\nALBM Data Summary:")
    print(f"  Time range: {albm_data.start_date} to {albm_data.end_date}")
    print(f"  Number of timesteps: {albm_data.n_timesteps}")
    print(f"  Water layers: {albm_data.n_water_layers}")
    print(f"  Sediment layers: {albm_data.n_sed_layers}")
    print(f"  Mean total flux: {np.mean(albm_data.total_flux):.2e} kg m⁻² s⁻¹")
    
    if eddy_data is not None:
        print(f"\nEddy Flux Data Summary:")
        print(f"  Date range: {eddy_data.datetime.min()} to {eddy_data.datetime.max()}")
        print(f"  Number of observations: {len(eddy_data.ch4_flux)}")
        print(f"  Mean flux: {np.mean(eddy_data.ch4_flux):.2e} kg m⁻² s⁻¹")
    else:
        print(f"\nEddy Flux Data: not loaded (file missing or not in cache)")
    return 0


# =============================================================================
# Main execution - inspect by default; use --refresh-cache to rewrite cache
# =============================================================================
if __name__ == '__main__':
    _exit_code = main()
    if not _running_in_interactive_kernel():
        sys.exit(_exit_code)
