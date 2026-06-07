"""
Isotope Model Library

Contains functions for calculating CH4 isotope mass balance and optimization.
This module provides the core isotope modeling and optimization capabilities.

Usage:
    from isotope_model import (
        calculate_dissolved_ch4_isotopes,
        calculate_surface_emission_isotopes,
        calculate_cost,
        run_optimization,
        get_default_obs_targets
    )
"""

import time
import numpy as np
import pandas as pd
from typing import Dict, Any, Tuple, Optional
from scipy.optimize import differential_evolution
from dataclasses import dataclass
import threading
from collections import OrderedDict
try:
    from tqdm import tqdm
except ImportError:
    # Fallback if tqdm not available
    def tqdm(iterable, *args, **kwargs):
        return iterable

# Try to import GPU libraries
GPU_AVAILABLE = False
GPU_INFO = {}
try:
    import cupy as cp
    GPU_AVAILABLE = True
    
    # Get GPU information
    try:
        device = cp.cuda.Device(0)
        GPU_INFO = {
            'name': 'NVIDIA GPU',  # Default name
            'compute_capability': device.compute_capability,
            'memory_total': device.mem_info[1] / (1024**3),  # GB
            'memory_free': device.mem_info[0] / (1024**3),   # GB
        }
        # Try to get device name properly
        try:
            props = cp.cuda.runtime.getDeviceProperties(0)
            if 'name' in props:
                GPU_INFO['name'] = props['name'].decode() if isinstance(props['name'], bytes) else props['name']
        except Exception:
            pass
    except Exception:
        pass
except ImportError:
    pass

# =============================================================================
# Constants
# =============================================================================

# VPDB standard 13C/12C ratio
R_STD = 0.0112372

# Default fractionation factors
DEFAULT_ALPHA_AM = 1.0238   # Acetoclastic methanogenesis
DEFAULT_ALPHA_HM = 1.0456   # Hydrogenotrophic methanogenesis
DEFAULT_ALPHA_MO = 1.0151   # Methanotrophy
DEFAULT_ALPHA_OP = 1.0238   # Oxic production
DEFAULT_ALPHA_TE = 1.000    # Ebullition (no fractionation)
DEFAULT_ALPHA_TD = 1.005    # Diffusion
DEFAULT_ALPHA_EX = 1.000    # Exchange (no fractionation)

# Default source parameters
DEFAULT_F_AM = 0.5          # Fraction acetoclastic methanogenesis
DEFAULT_C13_POM = -25.6     # δ13C of particulate organic matter


# Default observational isotope targets. Surface-emission targets are calculated
# as a flux-weighted mix of dissolved diffusive emissions and observed ebullition.
DEFAULT_OBS_DISSOLVED_SUMMER = -39.08
DEFAULT_OBS_DISSOLVED_WINTER = -63.56
DEFAULT_OBS_EBULLITION_SUMMER = -56.93
DEFAULT_OBS_EBULLITION_WINTER = -64.28
DEFAULT_OBS_DIFFUSION_FRACTION_SUMMER = 0.13
DEFAULT_OBS_DIFFUSION_FRACTION_WINTER = 0.24

DEFAULT_OBS_DISSOLVED_SUMMER_STD = 10.80
DEFAULT_OBS_DISSOLVED_WINTER_STD = 0.0738
DEFAULT_OBS_EBULLITION_SUMMER_STD = 0.223
DEFAULT_OBS_EBULLITION_WINTER_STD = 6.19

# Minimum internally modeled dissolved CH4 concentration required for treating
# dissolved isotope ratios as resolved model predictions (kg m-3).
DEFAULT_DISSOLVED_CH4_MIN_CONC_FOR_ISOTOPE = 1e-5

# =============================================================================
# Isotope Conversion Functions
# =============================================================================

def delta_to_ratio(delta: float) -> float:
    """Convert delta notation (‰) to isotope ratio."""
    return R_STD * (1 + delta / 1000)


def ratio_to_delta(ratio: float) -> float:
    """Convert isotope ratio to delta notation (‰)."""
    return (ratio / R_STD - 1) * 1000


def mix_observed_surface_emission_delta(
    dissolved_delta: float,
    ebullition_delta: float,
    diffusion_fraction: float,
) -> float:
    """
    Calculate observed surface-emission delta from dissolved and ebullition means.

    The diffusion fraction is the fraction of total surface emission represented
    by diffusive flux. The remaining fraction is assigned to ebullition.
    """
    if not 0.0 <= diffusion_fraction <= 1.0:
        raise ValueError("diffusion_fraction must be between 0 and 1")
    return (
        diffusion_fraction * dissolved_delta
        + (1.0 - diffusion_fraction) * ebullition_delta
    )


def propagate_observed_surface_emission_std(
    dissolved_std: float,
    ebullition_std: float,
    diffusion_fraction: float,
) -> float:
    """
    Propagate uncertainty for a weighted dissolved/ebullition emission mixture.

    Assumes independent dissolved and ebullition uncertainties:
    sigma_total = sqrt((f_diff*sigma_diss)^2 + (f_ebull*sigma_ebull)^2).
    """
    if not 0.0 <= diffusion_fraction <= 1.0:
        raise ValueError("diffusion_fraction must be between 0 and 1")
    ebullition_fraction = 1.0 - diffusion_fraction
    return float(np.sqrt(
        (diffusion_fraction * dissolved_std) ** 2
        + (ebullition_fraction * ebullition_std) ** 2
    ))


def filter_delta_by_concentration(
    delta_values: np.ndarray,
    concentrations: np.ndarray,
    min_concentration: float = DEFAULT_DISSOLVED_CH4_MIN_CONC_FOR_ISOTOPE,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Mask isotope values when the modeled dissolved CH4 pool is too small to
    resolve a meaningful isotope ratio.
    """
    delta_values = np.asarray(delta_values, dtype=float)
    concentrations = np.asarray(concentrations, dtype=float)
    valid = (
        np.isfinite(delta_values)
        & np.isfinite(concentrations)
        & (concentrations >= min_concentration)
    )
    filtered = delta_values.copy()
    filtered[~valid] = np.nan
    return filtered, valid


# =============================================================================
# Isotope Calculation Functions
# =============================================================================

def calculate_dissolved_ch4_isotopes(
    sedch4df_data: np.ndarray,
    och4prod_data: np.ndarray,
    ch4exc_data: np.ndarray,
    och4_data: np.ndarray,
    dflux_data: np.ndarray,
    ch4conc_data: np.ndarray,
    alpha_am: float,
    alpha_hm: float,
    alpha_mo: float,
    alpha_op: float,
    alpha_te: float,
    alpha_td: float,
    alpha_ex: float,
    C13_POM: float,
    f_am: float,
    f_hm: float
) -> Dict[str, Any]:
    """
    Calculate dissolved CH4 isotope mass balance.
    
    Parameters:
    -----------
    sedch4df_data : array - Sediment CH4 diffusion flux (kg m⁻³ s⁻¹)
    och4prod_data : array - Oxic CH4 production flux (kg m⁻³ s⁻¹)
    ch4exc_data : array - CH4 exchange flux (kg m⁻³ s⁻¹)
    och4_data : array - CH4 oxidation consumption flux (kg m⁻³ s⁻¹)
    dflux_data : array - Diffusion flux at surface (kg m⁻³ s⁻¹)
    ch4conc_data : array - CH4 concentration (kg m⁻³)
    alpha_am : float - Fractionation factor for acetoclastic methanogenesis
    alpha_hm : float - Fractionation factor for hydrogenotrophic methanogenesis
    alpha_mo : float - Fractionation factor for methanotrophy
    alpha_op : float - Fractionation factor for oxic production
    alpha_te : float - Fractionation factor for ebullition (typically 1.0)
    alpha_td : float - Fractionation factor for diffusion
    alpha_ex : float - Fractionation factor for exchange (typically 1.0)
    C13_POM : float - δ13C of particulate organic matter (‰)
    f_am : float - Fraction of acetoclastic methanogenesis
    f_hm : float - Fraction of hydrogenotrophic methanogenesis
    
    Returns:
    --------
    dict containing:
        - ch4_conc_mb : array - CH4 concentration from mass balance
        - delta_conc_mb : array - δ13C of dissolved CH4 from mass balance
        - isotope_mass_mb : array - Isotope mass for stability tracking
        - C13_sed_prod : float - δ13C of sediment-produced methane
        - R_sed_prod : float - Isotope ratio of sediment-produced methane
    """
    # Calculate source ratios
    R_POM = delta_to_ratio(C13_POM)
    R_HM = R_POM / alpha_hm
    R_AM = R_POM / alpha_am
    
    # Weighted average in ratio space
    R_sed_prod = f_hm * R_HM + f_am * R_AM
    
    # Convert back to delta notation
    C13_sed_prod = ratio_to_delta(R_sed_prod)
    
    # Initialize arrays for mass balance
    n_steps = len(sedch4df_data)
    ch4_conc_mb = np.zeros(n_steps)
    delta_conc_mb = np.zeros(n_steps)
    isotope_mass_mb = np.zeros(n_steps)
    
    # Set initial conditions
    ch4_conc_mb[0] = ch4conc_data[0] if ch4conc_data[0] > 0 else 1e-10
    delta_conc_mb[0] = C13_sed_prod
    isotope_mass_mb[0] = ch4_conc_mb[0] * delta_to_ratio(delta_conc_mb[0])
    
    # Time stepping loop
    for i in range(1, n_steps):
        dt = 86400  # time step in seconds (1 day)
        
        # Previous values
        C_prev = ch4_conc_mb[i-1]
        R_prev = delta_to_ratio(delta_conc_mb[i-1])
        
        # Fluxes for this time step
        sedch4df_flux = sedch4df_data[i]
        och4prod_flux = och4prod_data[i]
        ch4exc_flux = ch4exc_data[i]
        och4_flux = -och4_data[i]
        dflux_val = -dflux_data[i]
        
        # Calculate isotope ratios for each flux
        R_sedch4df = R_sed_prod
        R_och4prod = R_POM / alpha_op
        R_ch4exc = R_prev if not np.isnan(R_prev) else R_sed_prod
        R_och4 = (R_prev / alpha_mo) if not np.isnan(R_prev) else (R_sed_prod / alpha_mo)
        R_dflux = (R_prev / alpha_td) if not np.isnan(R_prev) else (R_sed_prod / alpha_td)
        
        # Mass balance
        dC_dt = sedch4df_flux + och4prod_flux + ch4exc_flux + och4_flux + dflux_val
        
        # Isotope mass balance
        isotope_flux = (sedch4df_flux * R_sedch4df + 
                       och4prod_flux * R_och4prod +
                       ch4exc_flux * R_ch4exc +
                       och4_flux * R_och4 +
                       dflux_val * R_dflux)
        
        if np.isnan(isotope_flux):
            isotope_flux = 0.0
        
        # Update concentration
        C_new = C_prev + dC_dt * dt
        if C_new < 0:
            C_new = 0.0
        ch4_conc_mb[i] = C_new
        
        # Update isotope mass and ratio
        if C_new == 0:
            delta_conc_mb[i] = np.nan
            isotope_mass_mb[i] = 0.0
        elif C_new > 0:
            if C_prev > 0 and not np.isnan(delta_conc_mb[i-1]):
                isotope_mass_prev = isotope_mass_mb[i-1]
            else:
                isotope_mass_prev = 0.0
            
            isotope_mass_new = isotope_mass_prev + isotope_flux * dt
            isotope_mass_mb[i] = isotope_mass_new
            
            if C_new > 1e-12 and isotope_mass_new != 0:
                R_new = isotope_mass_new / C_new
                if R_new > 0:
                    delta_new = ratio_to_delta(R_new)
                    delta_conc_mb[i] = np.clip(delta_new, -100, 10)
                else:
                    delta_conc_mb[i] = ratio_to_delta(R_sed_prod)
            else:
                delta_conc_mb[i] = ratio_to_delta(R_sed_prod)
        else:
            delta_conc_mb[i] = ratio_to_delta(R_sed_prod)
            isotope_mass_mb[i] = 0.0
    
    return {
        'ch4_conc_mb': ch4_conc_mb,
        'delta_conc_mb': delta_conc_mb,
        'isotope_mass_mb': isotope_mass_mb,
        'C13_sed_prod': C13_sed_prod,
        'R_sed_prod': R_sed_prod
    }


def calculate_dissolved_and_bubble_pool_isotopes(
    sedch4df_data: np.ndarray,
    sedch4eb_data: np.ndarray,
    och4prod_data: np.ndarray,
    ch4exc_data: np.ndarray,
    och4_data: np.ndarray,
    dflux_data: np.ndarray,
    eflux_data: np.ndarray,
    ch4upb_data: np.ndarray,
    ch4conc_data: np.ndarray,
    R_sed_prod,             # scalar float or 1-D array (time-varying)
    alpha_mo: float,
    alpha_op: float,
    alpha_td: float = DEFAULT_ALPHA_TD,
    alpha_te: float = DEFAULT_ALPHA_TE,
) -> Dict[str, Any]:
    """
    Co-evolve dissolved CH4 and bubble pool isotopes in a single time loop.

    Physical description
    --------------------
    Dissolved pool sources/sinks (kg m⁻³ s⁻¹):
      + sedch4df   : sediment diffusion into dissolved pool           → at R_sed_prod
      + och4prod   : oxic production                                  → at R_sed_prod/alpha_op
      + ch4exc > 0 : bubble pool dissolving into water               → at R_bubble (prev step)
      - ch4exc < 0 : dissolved absorbing into bubbles (rare)         → at R_diss (current)
      - och4       : oxidation (fractionated)                        → R_diss / alpha_mo
      - dflux      : surface diffusion loss (fractionated)           → R_diss / alpha_td

    Bubble pool sources/sinks (same flux units):
      + sedch4eb   : sediment ebullition into bubble pool            → at R_sed_prod
      - ch4exc > 0 : bubbles dissolving into water                   → at R_bubble (prev step)
      + ch4exc < 0 : dissolved absorbing into bubbles                → at R_diss (current)
      - eflux      : surface ebullition to atmosphere                → at R_bubble (prev step)
      - ch4upb     : upward bubbling to ice/surface                  → at R_bubble (prev step)

    Coupling is handled sequentially at daily timesteps:
      dissolved update uses R_bubble from the *previous* step (explicit coupling),
      bubble update uses R_diss after the daily dissolved update.

    Parameters
    ----------
    sedch4df_data  : sediment diffusion flux entering dissolved pool (kg m⁻³ s⁻¹)
    sedch4eb_data  : sediment ebullition flux entering bubble pool   (kg m⁻³ s⁻¹)
    och4prod_data  : oxic CH4 production                             (kg m⁻³ s⁻¹)
    ch4exc_data    : net bubble→dissolved exchange (>0 = bubbles → dissolved)
    och4_data      : oxidation consumption                           (kg m⁻³ s⁻¹, positive)
    dflux_data     : surface diffusion flux (loss for dissolved)     (kg m⁻³ or m⁻² s⁻¹)
    eflux_data     : surface ebullition flux (loss for bubble pool)  (same units as dflux)
    ch4upb_data    : upward bubbling flux (loss for bubble pool)     (same units)
    ch4conc_data   : measured dissolved CH4 concentration for init   (kg m⁻³)
    R_sed_prod     : sediment production isotope ratio; scalar or array (time-varying)
    alpha_mo       : oxidation fractionation factor
    alpha_op       : oxic production fractionation factor
    alpha_td       : diffusion fractionation factor
    alpha_te       : ebullition fractionation factor (default 1.0)

    Returns
    -------
    dict with:
        delta_dissolved  : δ¹³C of dissolved CH4 (‰)
        delta_bubble     : δ¹³C of bubble pool (‰)
        R_bubble         : isotope ratio of bubble pool
        ch4_conc_mb      : dissolved CH4 concentration from mass balance
        iso_mass_diss    : dissolved isotope mass (diagnostic)
        iso_mass_bub     : bubble isotope mass (diagnostic)
        B_pool           : bubble pool concentration (diagnostic)
    """
    n = len(sedch4df_data)
    dt = 86400.0  # 1 day in seconds

    # Broadcast R_sed_prod to array for uniform indexing
    R_sp = np.broadcast_to(
        np.asarray(R_sed_prod, dtype=float),
        (n,)
    ).copy()

    # --- Dissolved pool arrays ---
    ch4_conc_mb  = np.zeros(n)
    R_diss       = np.zeros(n)
    delta_diss   = np.zeros(n)
    iso_mass_d   = np.zeros(n)

    # --- Bubble pool arrays ---
    B_pool       = np.zeros(n)
    R_bub        = np.zeros(n)
    delta_bub    = np.zeros(n)
    iso_mass_b   = np.zeros(n)
    dissolved_sink_scale = np.ones(n)
    bubble_sink_scale    = np.ones(n)

    delta_min = -100.0
    delta_max = 10.0
    dissolved_floor = 1e-12
    bubble_floor = 1e-12

    def _finite_flux(value: float) -> float:
        """Convert non-finite fluxes to zero before mass-balance arithmetic."""
        value = float(value)
        return value if np.isfinite(value) else 0.0

    def _limit_sink_fluxes(pool_mass: float, source_flux: float, sink_fluxes, step_dt: float):
        """
        Proportionally limit same-step sinks so they cannot remove more mass
        than the pool plus current-step sources can supply.
        """
        positive_sinks = [max(_finite_flux(sink), 0.0) for sink in sink_fluxes]
        total_sink_flux = sum(positive_sinks)
        if total_sink_flux <= 0.0:
            return positive_sinks, 1.0

        available_mass = max(_finite_flux(pool_mass), 0.0) + max(_finite_flux(source_flux), 0.0) * step_dt
        requested_sink_mass = total_sink_flux * step_dt
        if requested_sink_mass <= available_mass:
            return positive_sinks, 1.0

        scale = max(available_mass / requested_sink_mass, 0.0) if requested_sink_mass > 0.0 else 1.0
        return [sink * scale for sink in positive_sinks], scale

    def _consistent_pool_state(
        pool_mass: float,
        isotope_mass: float,
        fallback_ratio: float,
        floor: float,
        depleted_mass: float,
    ):
        """
        Return concentration, isotope mass, ratio, and delta as one consistent
        state. If the pool is depleted or isotope mass becomes invalid, reset
        to the fallback ratio instead of carrying a hidden invalid isotope mass.
        """
        pool_mass = max(_finite_flux(pool_mass), 0.0)
        fallback_ratio = fallback_ratio if np.isfinite(fallback_ratio) and fallback_ratio > 0.0 else R_STD

        if pool_mass <= floor:
            pool_mass = depleted_mass
            ratio = fallback_ratio
        elif not np.isfinite(isotope_mass) or isotope_mass <= 0.0:
            ratio = fallback_ratio
        else:
            ratio = isotope_mass / pool_mass
            if not np.isfinite(ratio) or ratio <= 0.0:
                ratio = fallback_ratio

        delta = float(np.clip(ratio_to_delta(ratio), delta_min, delta_max))
        ratio = delta_to_ratio(delta)
        isotope_mass = pool_mass * ratio
        return pool_mass, isotope_mass, ratio, delta

    # ---- Initial conditions ----
    R0 = R_sp[0]
    ch4_conc_mb[0] = max(float(ch4conc_data[0]), 1e-10)
    R_diss[0]      = R0
    delta_diss[0]  = ratio_to_delta(R0)
    iso_mass_d[0]  = ch4_conc_mb[0] * R0

    # Bubble pool: seed with one day's sediment ebullition
    B_pool[0]     = max(float(sedch4eb_data[0]) * dt, 1e-12)
    R_bub[0]      = R0
    delta_bub[0]  = ratio_to_delta(R0)
    iso_mass_b[0] = B_pool[0] * R0

    for i in range(1, n):
        R_sed = R_sp[i]

        # Previous dissolved state
        C_d    = ch4_conc_mb[i - 1]
        Rd     = R_diss[i - 1]
        iso_d  = iso_mass_d[i - 1]

        # Previous bubble state
        B_b    = B_pool[i - 1]
        Rb     = R_bub[i - 1]
        iso_b  = iso_mass_b[i - 1]

        # --- Exchange splitting ---
        exc      = float(ch4exc_data[i])
        exc_bub_to_diss = max(exc,  0.0)   # bubbles dissolving → source for dissolved
        exc_diss_to_bub = max(-exc, 0.0)   # dissolved absorbing into bubbles → sink for dissolved

        sedch4df_flux = _finite_flux(sedch4df_data[i])
        sedch4eb_flux = _finite_flux(sedch4eb_data[i])
        och4prod_flux = _finite_flux(och4prod_data[i])
        och4_flux = max(_finite_flux(och4_data[i]), 0.0)
        dflux_val = max(_finite_flux(dflux_data[i]), 0.0)
        eflux_val = max(_finite_flux(eflux_data[i]), 0.0)
        ch4upb_val = max(_finite_flux(ch4upb_data[i]), 0.0)

        C_d, iso_d, Rd, _ = _consistent_pool_state(
            C_d,
            iso_d,
            R_sed,
            dissolved_floor,
            0.0
        )
        B_b, iso_b, Rb, _ = _consistent_pool_state(
            B_b,
            iso_b,
            R_sed,
            bubble_floor,
            bubble_floor
        )

        # Limit same-day sinks directly, without substeps or isotope solves.
        (exc_bub_to_diss_limited,
         eflux_limited,
         ch4upb_limited), bubble_sink_scale[i] = _limit_sink_fluxes(
            B_b,
            sedch4eb_flux,
            (exc_bub_to_diss, eflux_val, ch4upb_val),
            dt
        )

        dissolved_source_flux = sedch4df_flux + och4prod_flux + exc_bub_to_diss_limited
        (exc_diss_to_bub_limited,
         och4_flux_limited,
         dflux_limited), dissolved_sink_scale[i] = _limit_sink_fluxes(
            C_d,
            dissolved_source_flux,
            (exc_diss_to_bub, och4_flux, dflux_val),
            dt
        )

        # --- Dissolved pool update ---
        R_dflux  = Rd / alpha_td
        R_och4   = Rd / alpha_mo
        R_op     = R_sed / alpha_op

        dC_d_dt = (sedch4df_flux
                   + och4prod_flux
                   + exc_bub_to_diss_limited
                   - exc_diss_to_bub_limited
                   - och4_flux_limited
                   - dflux_limited)

        iso_flux_d = (sedch4df_flux * R_sed
                      + och4prod_flux * R_op
                      + exc_bub_to_diss_limited * Rb
                      - exc_diss_to_bub_limited * Rd
                      - och4_flux_limited * R_och4
                      - dflux_limited * R_dflux)

        if np.isnan(iso_flux_d):
            iso_flux_d = 0.0

        C_d_new  = C_d + dC_d_dt * dt
        iso_d_new = iso_d + iso_flux_d * dt
        C_d_new, iso_d_new, Rd_new, delta_diss[i] = _consistent_pool_state(
            C_d_new,
            iso_d_new,
            R_sed,
            dissolved_floor,
            0.0
        )

        ch4_conc_mb[i] = C_d_new
        iso_mass_d[i]  = iso_d_new
        R_diss[i] = Rd_new

        # --- Bubble pool update (uses freshly computed dissolved R for gain) ---
        bubble_sed_source_mass = max(sedch4eb_flux, 0.0) * dt
        bubble_exchange_source_mass = exc_diss_to_bub_limited * dt
        bubble_source_mass = bubble_sed_source_mass + bubble_exchange_source_mass
        bubble_source_iso = (
            bubble_sed_source_mass * R_sed
            + bubble_exchange_source_mass * Rd_new
        )
        if bubble_source_mass > bubble_floor and np.isfinite(bubble_source_iso):
            R_bubble_source = bubble_source_iso / bubble_source_mass
        else:
            R_bubble_source = R_sed

        bubble_exchange_loss_mass = exc_bub_to_diss_limited * dt
        bubble_eflux_loss_mass = eflux_limited * dt
        bubble_upb_loss_mass = ch4upb_limited * dt
        bubble_loss_mass = (
            bubble_exchange_loss_mass
            + bubble_eflux_loss_mass
            + bubble_upb_loss_mass
        )

        # Isotope-aware sink split: old storage can only supply B_b of loss.
        # Any excess same-day loss is assigned to same-day bubble sources.
        if bubble_loss_mass > 0.0:
            old_pool_loss_fraction = min(max(B_b, 0.0), bubble_loss_mass) / bubble_loss_mass
            source_loss_fraction = 1.0 - old_pool_loss_fraction
        else:
            old_pool_loss_fraction = 0.0
            source_loss_fraction = 0.0

        old_pool_iso_loss = old_pool_loss_fraction * (
            bubble_exchange_loss_mass * Rb
            + bubble_eflux_loss_mass * Rb * alpha_te
            + bubble_upb_loss_mass * Rb
        )
        source_iso_loss = source_loss_fraction * (
            bubble_exchange_loss_mass * R_bubble_source
            + bubble_eflux_loss_mass * R_bubble_source * alpha_te
            + bubble_upb_loss_mass * R_bubble_source
        )

        B_new = (
            B_b
            + (sedch4eb_flux + exc_diss_to_bub_limited) * dt
            - bubble_loss_mass
        )
        iso_b_new = iso_b + bubble_source_iso - old_pool_iso_loss - source_iso_loss
        B_new, iso_b_new, Rb_new, delta_bub[i] = _consistent_pool_state(
            B_new,
            iso_b_new,
            R_sed,
            bubble_floor,
            bubble_floor
        )

        B_pool[i]     = B_new
        iso_mass_b[i] = iso_b_new
        R_bub[i] = Rb_new

    return {
        'delta_dissolved': delta_diss,
        'delta_bubble':    delta_bub,
        'R_bubble':        R_bub,
        'ch4_conc_mb':     ch4_conc_mb,
        'iso_mass_diss':   iso_mass_d,
        'iso_mass_bub':    iso_mass_b,
        'B_pool':          B_pool,
        'dissolved_sink_scale': dissolved_sink_scale,
        'bubble_sink_scale':    bubble_sink_scale,
    }


def calculate_surface_emission_isotopes(
    dflux_data: np.ndarray,
    eflux_data: np.ndarray,
    ch4upb_data: np.ndarray,
    delta_conc_mb: np.ndarray,
    R_sed_prod: float,
    alpha_td: float = DEFAULT_ALPHA_TD,
    delta_bubble: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Calculate surface emission isotope signatures.

    When delta_bubble is provided (bubble pool tracking active):
      - Diffusion uses the dissolved pool delta (fractionated by alpha_td).
      - Ebullition and upward bubbling use the bubble pool delta.
    Otherwise falls back to R_sed_prod for ebullition/upward-bubbling (legacy behaviour).

    Parameters
    ----------
    delta_bubble : optional array - δ¹³C of bubble pool; if None uses R_sed_prod for eflux/ch4upb.
    """
    n_steps = len(delta_conc_mb)
    total_surface_flux = dflux_data + eflux_data + ch4upb_data

    isotope_flux_dflux  = np.zeros(n_steps)
    isotope_flux_eflux  = np.zeros(n_steps)
    isotope_flux_ch4upb = np.zeros(n_steps)
    delta_surface_emission = np.zeros(n_steps)

    use_bubble = (delta_bubble is not None)

    for i in range(n_steps):
        # ---- Diffusion: always from dissolved pool (fractionated) ----
        if dflux_data[i] > 0:
            if not np.isnan(delta_conc_mb[i]):
                isotope_flux_dflux[i] = dflux_data[i] * delta_to_ratio(delta_conc_mb[i])
            else:
                isotope_flux_dflux[i] = dflux_data[i] * R_sed_prod

        # ---- Ebullition & upward bubbling: from bubble pool if tracked ----
        if use_bubble and not np.isnan(delta_bubble[i]):
            R_bub_i = delta_to_ratio(delta_bubble[i])
        else:
            R_bub_i = R_sed_prod

        if eflux_data[i] > 0:
            isotope_flux_eflux[i]  = eflux_data[i]  * R_bub_i
        if ch4upb_data[i] > 0:
            isotope_flux_ch4upb[i] = ch4upb_data[i] * R_bub_i

        # ---- Weighted emission delta ----
        if total_surface_flux[i] > 0:
            dflux_contrib  = max(0.0, dflux_data[i])  / total_surface_flux[i]
            eflux_contrib  = max(0.0, eflux_data[i])  / total_surface_flux[i]
            ch4upb_contrib = max(0.0, ch4upb_data[i]) / total_surface_flux[i]

            delta_dflux  = delta_conc_mb[i] if not np.isnan(delta_conc_mb[i]) else ratio_to_delta(R_sed_prod)
            delta_eflux  = ratio_to_delta(R_bub_i)
            delta_ch4upb = ratio_to_delta(R_bub_i)

            delta_surface_emission[i] = (dflux_contrib  * delta_dflux +
                                         eflux_contrib  * delta_eflux +
                                         ch4upb_contrib * delta_ch4upb)
        else:
            delta_surface_emission[i] = (delta_conc_mb[i] if not np.isnan(delta_conc_mb[i])
                                         else ratio_to_delta(R_sed_prod))

    return {
        'delta_surface_emission': delta_surface_emission,
        'isotope_flux_dflux':     isotope_flux_dflux,
        'isotope_flux_eflux':     isotope_flux_eflux,
        'total_surface_flux':     total_surface_flux,
    }


def compute_isotope_timeseries_temp(
    params: Dict[str, float],
    flux_data: Dict[str, np.ndarray]
) -> Dict[str, Any]:
    """
    Compute isotope time series using temperature-based sediment production.
    C13_sed_prod(t) = m * T_flux_weighted(t) + b
    Tracks both the dissolved pool and the bubble pool so surface emissions
    use the bubble pool isotopic composition, not the raw sediment production.

    Parameters
    ----------
    params : dict  - keys: m, b, alpha_mo, alpha_op
    flux_data : dict - must include 'sediment_temp_avg' and 'sedch4eb_data'
    """
    temp_avg = flux_data.get('sediment_temp_avg')
    if temp_avg is None:
        raise ValueError(
            "flux_data must include valid 'sediment_temp_avg' data for "
            "temperature-based calculation"
        )
    temp_avg = np.asarray(temp_avg, dtype=float)
    if temp_avg.ndim != 1 or temp_avg.size == 0 or not np.isfinite(temp_avg).all():
        raise ValueError(
            "flux_data must include a non-empty fully finite 'sediment_temp_avg' array "
            "for temperature-based calculation"
        )
    m = params['m']
    b = params['b']

    C13_sed_prod_time = m * temp_avg + b
    R_sed_prod_time   = delta_to_ratio(C13_sed_prod_time)

    alpha_mo = params.get('alpha_mo', DEFAULT_ALPHA_MO)
    alpha_op = params.get('alpha_op', DEFAULT_ALPHA_OP)

    # Co-evolve dissolved and bubble pools with time-varying R_sed_prod
    pools = calculate_dissolved_and_bubble_pool_isotopes(
        sedch4df_data=flux_data['sedch4df_data'],
        sedch4eb_data=flux_data['sedch4eb_data'],
        och4prod_data=flux_data['och4prod_data'],
        ch4exc_data=flux_data['ch4exc_data'],
        och4_data=flux_data['och4_data'],
        dflux_data=flux_data['dflux_data'],
        eflux_data=flux_data['eflux_data'],
        ch4upb_data=flux_data['ch4upb_data'],
        ch4conc_data=flux_data['ch4conc_data'],
        R_sed_prod=R_sed_prod_time,
        alpha_mo=alpha_mo,
        alpha_op=alpha_op,
        alpha_td=DEFAULT_ALPHA_TD,
        alpha_te=DEFAULT_ALPHA_TE,
    )

    delta_dissolved_raw = pools['delta_dissolved']
    delta_dissolved, dissolved_concentration_valid = filter_delta_by_concentration(
        delta_dissolved_raw,
        pools['ch4_conc_mb']
    )

    emissions = calculate_surface_emission_isotopes_temp(
        dflux_data=flux_data['dflux_data'],
        eflux_data=flux_data['eflux_data'],
        ch4upb_data=flux_data['ch4upb_data'],
        delta_conc_mb=delta_dissolved,
        R_sed_prod_time=R_sed_prod_time,
        alpha_td=DEFAULT_ALPHA_TD,
        delta_bubble=pools['delta_bubble'],
    )

    return {
        'delta_dissolved':    delta_dissolved,
        'delta_dissolved_raw': delta_dissolved_raw,
        'dissolved_concentration_valid': dissolved_concentration_valid,
        'dissolved_concentration_filter_min': DEFAULT_DISSOLVED_CH4_MIN_CONC_FOR_ISOTOPE,
        'delta_emission':     emissions['delta_surface_emission'],
        'delta_bubble':       pools['delta_bubble'],
        'C13_sed_prod':       C13_sed_prod_time,
        'C13_sed_prod_mean':  float(np.mean(C13_sed_prod_time)),
        'ch4_conc_mb':        pools['ch4_conc_mb'],
        'iso_mass_diss':      pools['iso_mass_diss'],
        'iso_mass_bub':       pools['iso_mass_bub'],
        'B_pool':             pools['B_pool'],
        'dissolved_sink_scale': pools['dissolved_sink_scale'],
        'bubble_sink_scale':    pools['bubble_sink_scale'],
    }


def calculate_surface_emission_isotopes_temp(
    dflux_data: np.ndarray,
    eflux_data: np.ndarray,
    ch4upb_data: np.ndarray,
    delta_conc_mb: np.ndarray,
    R_sed_prod_time: np.ndarray,
    alpha_td: float,
    delta_bubble: Optional[np.ndarray] = None,
) -> Dict[str, Any]:
    """
    Calculate surface emission isotopes using time-varying sediment production.
    When delta_bubble is provided, ebullition and upward bubbling use the bubble pool delta.
    """
    n_steps = len(delta_conc_mb)
    total_surface_flux = dflux_data + eflux_data + ch4upb_data

    isotope_flux_dflux  = np.zeros(n_steps)
    isotope_flux_eflux  = np.zeros(n_steps)
    isotope_flux_ch4upb = np.zeros(n_steps)
    delta_surface_emission = np.zeros(n_steps)

    use_bubble = (delta_bubble is not None)

    for i in range(n_steps):
        R_sed_prod = R_sed_prod_time[i]

        if dflux_data[i] > 0:
            if not np.isnan(delta_conc_mb[i]):
                isotope_flux_dflux[i] = dflux_data[i] * delta_to_ratio(delta_conc_mb[i])
            else:
                isotope_flux_dflux[i] = dflux_data[i] * R_sed_prod

        if use_bubble and not np.isnan(delta_bubble[i]):
            R_bub_i = delta_to_ratio(delta_bubble[i])
        else:
            R_bub_i = R_sed_prod

        if eflux_data[i] > 0:
            isotope_flux_eflux[i]  = eflux_data[i]  * R_bub_i
        if ch4upb_data[i] > 0:
            isotope_flux_ch4upb[i] = ch4upb_data[i] * R_bub_i

        if total_surface_flux[i] > 0:
            dflux_contrib  = max(0.0, dflux_data[i])  / total_surface_flux[i]
            eflux_contrib  = max(0.0, eflux_data[i])  / total_surface_flux[i]
            ch4upb_contrib = max(0.0, ch4upb_data[i]) / total_surface_flux[i]

            delta_dflux  = delta_conc_mb[i] if not np.isnan(delta_conc_mb[i]) else ratio_to_delta(R_sed_prod)
            delta_eflux  = ratio_to_delta(R_bub_i)
            delta_ch4upb = ratio_to_delta(R_bub_i)

            delta_surface_emission[i] = (dflux_contrib  * delta_dflux +
                                         eflux_contrib  * delta_eflux +
                                         ch4upb_contrib * delta_ch4upb)
        else:
            delta_surface_emission[i] = (delta_conc_mb[i] if not np.isnan(delta_conc_mb[i])
                                         else ratio_to_delta(R_sed_prod))

    return {
        'delta_surface_emission': delta_surface_emission,
        'isotope_flux_dflux':     isotope_flux_dflux,
        'isotope_flux_eflux':     isotope_flux_eflux,
        'total_surface_flux':     total_surface_flux,
    }


def calculate_dissolved_ch4_isotopes_temp(
    sedch4df_data: np.ndarray,
    och4prod_data: np.ndarray,
    ch4exc_data: np.ndarray,
    och4_data: np.ndarray,
    dflux_data: np.ndarray,
    ch4conc_data: np.ndarray,
    R_sed_prod_time: np.ndarray,  # Time-varying sediment production ratio
    alpha_mo: float,
    alpha_op: float,
    alpha_te: float,
    alpha_td: float,
    alpha_ex: float
) -> Dict[str, Any]:
    """
    Calculate dissolved CH4 isotopes with time-varying sediment production.
    Similar to calculate_dissolved_ch4_isotopes but uses time-varying R_sed_prod.
    """
    n_steps = len(sedch4df_data)
    
    # Initialize arrays
    ch4_conc_mb = np.zeros(n_steps)
    delta_conc_mb = np.zeros(n_steps)
    isotope_mass_mb = np.zeros(n_steps)
    
    dt = 86400  # time step in seconds (1 day)
    
    # Initial conditions
    ch4_conc_mb[0] = max(ch4conc_data[0], 1e-10)
    # Use first timestep's R_sed_prod for initial delta
    R_initial = R_sed_prod_time[0]
    delta_conc_mb[0] = ratio_to_delta(R_initial)
    isotope_mass_mb[0] = ch4_conc_mb[0] * R_initial
    
    # Time stepping
    for i in range(1, n_steps):
        C_prev = ch4_conc_mb[i-1]
        delta_prev = delta_conc_mb[i-1]
        isotope_mass_prev = isotope_mass_mb[i-1]
        
        # Get time-varying sediment production ratio
        R_sed_prod = R_sed_prod_time[i]
        
        # Calculate R_prev
        R_prev = delta_to_ratio(delta_prev) if not np.isnan(delta_prev) else R_sed_prod
        
        # Isotope ratios for other processes
        R_och4prod = R_sed_prod / alpha_op  # Using sediment production as base
        R_ch4exc = R_prev if not np.isnan(R_prev) else R_sed_prod
        R_och4 = (R_prev / alpha_mo) if not np.isnan(R_prev) else (R_sed_prod / alpha_mo)
        R_dflux = (R_prev / alpha_td) if not np.isnan(R_prev) else (R_sed_prod / alpha_td)
        
        # Mass balance
        dC_dt = sedch4df_data[i] + och4prod_data[i] + ch4exc_data[i] - och4_data[i] - dflux_data[i]
        
        # Isotope flux
        isotope_flux = (sedch4df_data[i] * R_sed_prod + 
                       och4prod_data[i] * R_och4prod +
                       ch4exc_data[i] * R_ch4exc -
                       och4_data[i] * R_och4 -
                       dflux_data[i] * R_dflux)
        
        if np.isnan(isotope_flux):
            isotope_flux = 0.0
        
        C_new = C_prev + dC_dt * dt
        if C_new < 0:
            C_new = 0
        
        ch4_conc_mb[i] = C_new
        
        # Isotope mass balance
        if C_prev > 0 and not np.isnan(delta_conc_mb[i-1]):
            isotope_mass_prev = isotope_mass_mb[i-1]
            isotope_mass_new = isotope_mass_prev + isotope_flux * dt
            isotope_mass_mb[i] = isotope_mass_new
            
            # Calculate delta
            if C_new > 1e-12 and isotope_mass_new != 0:
                R_new = isotope_mass_new / C_new
                delta_new = ratio_to_delta(R_new)
                delta_new = max(min(delta_new, 10.0), -100.0)
            else:
                delta_new = ratio_to_delta(R_sed_prod)
        else:
            isotope_mass_mb[i] = 0.0
            delta_new = ratio_to_delta(R_sed_prod)
        
        delta_conc_mb[i] = delta_new
    
    return {
        'ch4_conc_mb': ch4_conc_mb,
        'delta_conc_mb': delta_conc_mb,
        'isotope_mass_mb': isotope_mass_mb
    }


def compute_isotope_timeseries(
    params: Dict[str, float],
    flux_data: Dict[str, np.ndarray]
) -> Dict[str, Any]:
    """
    Compute isotope time series for given parameters.
    Tracks both the dissolved pool and the bubble pool so that surface ebullition
    and upward-bubbling carry the isotopic composition of the bubble pool rather
    than the raw sediment production value.

    Parameters
    ----------
    params : dict  - keys: alpha_am, alpha_hm, alpha_mo, alpha_op, f_am, C13_POM
    flux_data : dict - flux data arrays from ALBM (must include 'sedch4eb_data')
    """
    # Derive sediment production isotope ratio (constant)
    R_POM    = delta_to_ratio(params['C13_POM'])
    f_am     = params['f_am']
    f_hm     = 1.0 - f_am
    R_sed    = f_hm * (R_POM / params['alpha_hm']) + f_am * (R_POM / params['alpha_am'])
    C13_sed  = ratio_to_delta(R_sed)

    # Co-evolve dissolved and bubble pools
    pools = calculate_dissolved_and_bubble_pool_isotopes(
        sedch4df_data=flux_data['sedch4df_data'],
        sedch4eb_data=flux_data['sedch4eb_data'],
        och4prod_data=flux_data['och4prod_data'],
        ch4exc_data=flux_data['ch4exc_data'],
        och4_data=flux_data['och4_data'],
        dflux_data=flux_data['dflux_data'],
        eflux_data=flux_data['eflux_data'],
        ch4upb_data=flux_data['ch4upb_data'],
        ch4conc_data=flux_data['ch4conc_data'],
        R_sed_prod=R_sed,
        alpha_mo=params['alpha_mo'],
        alpha_op=params['alpha_op'],
        alpha_td=DEFAULT_ALPHA_TD,
        alpha_te=DEFAULT_ALPHA_TE,
    )

    delta_dissolved_raw = pools['delta_dissolved']
    delta_dissolved, dissolved_concentration_valid = filter_delta_by_concentration(
        delta_dissolved_raw,
        pools['ch4_conc_mb']
    )

    # Surface emissions: diffusion from dissolved, ebullition/upwelling from bubble pool
    emissions = calculate_surface_emission_isotopes(
        dflux_data=flux_data['dflux_data'],
        eflux_data=flux_data['eflux_data'],
        ch4upb_data=flux_data['ch4upb_data'],
        delta_conc_mb=delta_dissolved,
        R_sed_prod=R_sed,
        alpha_td=DEFAULT_ALPHA_TD,
        delta_bubble=pools['delta_bubble'],
    )

    return {
        'delta_dissolved':  delta_dissolved,
        'delta_dissolved_raw': delta_dissolved_raw,
        'dissolved_concentration_valid': dissolved_concentration_valid,
        'dissolved_concentration_filter_min': DEFAULT_DISSOLVED_CH4_MIN_CONC_FOR_ISOTOPE,
        'delta_emission':   emissions['delta_surface_emission'],
        'delta_bubble':     pools['delta_bubble'],
        'C13_sed_prod':     C13_sed,
        'R_sed_prod':       R_sed,
        'ch4_conc_mb':      pools['ch4_conc_mb'],
        'iso_mass_diss':    pools['iso_mass_diss'],
        'iso_mass_bub':     pools['iso_mass_bub'],
        'B_pool':           pools['B_pool'],
        'dissolved_sink_scale': pools['dissolved_sink_scale'],
        'bubble_sink_scale':    pools['bubble_sink_scale'],
    }


# =============================================================================
# Seasonal Mean Calculation
# =============================================================================

def get_seasonal_mean(
    data: np.ndarray,
    time_array: np.ndarray,
    season: str
) -> float:
    """
    Calculate mean value for a given season across all years.
    
    Parameters:
    -----------
    data : array - Data array to average
    time_array : array - Corresponding time array
    season : str - 'summer' or 'winter'
    
    Returns:
    --------
    float : Seasonal mean value
    """
    time_pd = pd.DatetimeIndex(time_array)
    
    if season == 'summer':
        # Summer: June 21 - September 22
        mask = ((time_pd.month > 6) | ((time_pd.month == 6) & (time_pd.day >= 21))) & \
               ((time_pd.month < 9) | ((time_pd.month == 9) & (time_pd.day <= 22)))
    elif season == 'winter':
        # Winter: December 21 - March 20
        mask = ((time_pd.month == 12) & (time_pd.day >= 21)) | \
               (time_pd.month <= 2) | \
               ((time_pd.month == 3) & (time_pd.day <= 20))
    else:
        raise ValueError(f"Unknown season: {season}")
    
    valid_data = data[mask]
    valid_data = valid_data[~np.isnan(valid_data)]

    return np.mean(valid_data) if len(valid_data) > 0 else np.nan


# =============================================================================
# Observational Targets
# =============================================================================

def get_default_obs_targets() -> Dict[str, Dict[str, float]]:
    """
    Get default observational targets for optimization.
    
    Returns:
    --------
    dict : Observational targets with mean and std for each target
    """
    emissions_summer = mix_observed_surface_emission_delta(
        DEFAULT_OBS_DISSOLVED_SUMMER,
        DEFAULT_OBS_EBULLITION_SUMMER,
        DEFAULT_OBS_DIFFUSION_FRACTION_SUMMER,
    )
    emissions_winter = mix_observed_surface_emission_delta(
        DEFAULT_OBS_DISSOLVED_WINTER,
        DEFAULT_OBS_EBULLITION_WINTER,
        DEFAULT_OBS_DIFFUSION_FRACTION_WINTER,
    )
    emissions_summer_std = propagate_observed_surface_emission_std(
        DEFAULT_OBS_DISSOLVED_SUMMER_STD,
        DEFAULT_OBS_EBULLITION_SUMMER_STD,
        DEFAULT_OBS_DIFFUSION_FRACTION_SUMMER,
    )
    emissions_winter_std = propagate_observed_surface_emission_std(
        DEFAULT_OBS_DISSOLVED_WINTER_STD,
        DEFAULT_OBS_EBULLITION_WINTER_STD,
        DEFAULT_OBS_DIFFUSION_FRACTION_WINTER,
    )
    return {
        'dissolved_summer': {
            'mean': DEFAULT_OBS_DISSOLVED_SUMMER,
            'std': DEFAULT_OBS_DISSOLVED_SUMMER_STD,
        },
        'dissolved_winter': {
            'mean': DEFAULT_OBS_DISSOLVED_WINTER,
            'std': DEFAULT_OBS_DISSOLVED_WINTER_STD,
        },
        'emissions_summer': {
            'mean': emissions_summer,
            'std': emissions_summer_std,
            'dissolved_mean': DEFAULT_OBS_DISSOLVED_SUMMER,
            'dissolved_std': DEFAULT_OBS_DISSOLVED_SUMMER_STD,
            'ebullition_mean': DEFAULT_OBS_EBULLITION_SUMMER,
            'ebullition_std': DEFAULT_OBS_EBULLITION_SUMMER_STD,
            'diffusion_fraction': DEFAULT_OBS_DIFFUSION_FRACTION_SUMMER,
        },
        'emissions_winter': {
            'mean': emissions_winter,
            'std': emissions_winter_std,
            'dissolved_mean': DEFAULT_OBS_DISSOLVED_WINTER,
            'dissolved_std': DEFAULT_OBS_DISSOLVED_WINTER_STD,
            'ebullition_mean': DEFAULT_OBS_EBULLITION_WINTER,
            'ebullition_std': DEFAULT_OBS_EBULLITION_WINTER_STD,
            'diffusion_fraction': DEFAULT_OBS_DIFFUSION_FRACTION_WINTER,
        },
    }


# =============================================================================
# Cost Function
# =============================================================================

def _target_error_or_nan(
    model_val: float,
    obs_val: float,
    obs_std: float,
    normalize_by_std: bool,
) -> float:
    """Return target error, or NaN when the model has no valid prediction."""
    if not np.isfinite(model_val):
        return np.nan
    if normalize_by_std:
        return ((model_val - obs_val) / obs_std) ** 2
    return (model_val - obs_val) ** 2


def calculate_cost(
    params: Dict[str, float],
    flux_data: Dict[str, np.ndarray],
    obs_targets: Dict[str, Dict[str, float]],
    target_toggles: Dict[str, bool],
    normalize_by_std: bool = False
) -> Tuple[float, Dict[str, float], Dict[str, float]]:
    """
    Calculate cost function for optimization.
    
    Parameters:
    -----------
    params : dict
        Parameter dictionary with keys: alpha_am, alpha_hm, alpha_mo, alpha_op, f_am, C13_POM
    flux_data : dict
        Flux data dictionary from ALBM
    obs_targets : dict
        Observational targets dictionary
    target_toggles : dict
        Boolean toggles for which targets to include
    normalize_by_std : bool
        Whether to normalize errors by standard deviation
    
    Returns:
    --------
    tuple : (total_cost, individual_costs, model_predictions)
    """
    # Run isotope calculation
    timeseries = compute_isotope_timeseries(params, flux_data)
    
    delta_conc = timeseries['delta_dissolved']
    delta_emission = timeseries['delta_emission']
    time_arr = flux_data['time_array']
    
    model_predictions = {}
    individual_costs = {}
    total_cost = 0.0
    
    # Dissolved summer
    if target_toggles.get('dissolved_summer', False):
        model_val = get_seasonal_mean(delta_conc, time_arr, 'summer')
        obs_val = obs_targets['dissolved_summer']['mean']
        obs_std = obs_targets['dissolved_summer']['std']
        error = _target_error_or_nan(model_val, obs_val, obs_std, normalize_by_std)
        model_predictions['dissolved_summer'] = model_val
        individual_costs['dissolved_summer'] = error
        if np.isfinite(error):
            total_cost += error
    
    # Dissolved winter
    if target_toggles.get('dissolved_winter', False):
        model_val = get_seasonal_mean(delta_conc, time_arr, 'winter')
        obs_val = obs_targets['dissolved_winter']['mean']
        obs_std = obs_targets['dissolved_winter']['std']
        error = _target_error_or_nan(model_val, obs_val, obs_std, normalize_by_std)
        model_predictions['dissolved_winter'] = model_val
        individual_costs['dissolved_winter'] = error
        if np.isfinite(error):
            total_cost += error
    
    # Emissions summer
    if target_toggles.get('emissions_summer', False):
        model_val = get_seasonal_mean(delta_emission, time_arr, 'summer')
        obs_val = obs_targets['emissions_summer']['mean']
        obs_std = obs_targets['emissions_summer']['std']
        error = _target_error_or_nan(model_val, obs_val, obs_std, normalize_by_std)
        model_predictions['emissions_summer'] = model_val
        individual_costs['emissions_summer'] = error
        if np.isfinite(error):
            total_cost += error
    
    # Emissions winter
    if target_toggles.get('emissions_winter', False):
        model_val = get_seasonal_mean(delta_emission, time_arr, 'winter')
        obs_val = obs_targets['emissions_winter']['mean']
        obs_std = obs_targets['emissions_winter']['std']
        error = _target_error_or_nan(model_val, obs_val, obs_std, normalize_by_std)
        model_predictions['emissions_winter'] = model_val
        individual_costs['emissions_winter'] = error
        if np.isfinite(error):
            total_cost += error
    
    return total_cost, individual_costs, model_predictions


def calculate_cost_temp(
    params: Dict[str, float],
    flux_data: Dict[str, np.ndarray],
    obs_targets: Dict[str, Dict[str, float]],
    target_toggles: Dict[str, bool],
    normalize_by_std: bool = False
) -> Tuple[float, Dict[str, float], Dict[str, float]]:
    """
    Calculate cost function for temperature-based optimization.
    
    Parameters:
    -----------
    params : dict
        Parameter dictionary with keys: m (slope), b (intercept), alpha_mo, alpha_op
    flux_data : dict
        Flux data dictionary from ALBM (must include 'sediment_temp_avg')
    obs_targets : dict
        Observational targets dictionary
    target_toggles : dict
        Boolean toggles for which targets to include
    normalize_by_std : bool
        Whether to normalize errors by standard deviation
    
    Returns:
    --------
    tuple : (total_cost, individual_costs, model_predictions)
    """
    # Run temperature-based isotope calculation
    timeseries = compute_isotope_timeseries_temp(params, flux_data)
    
    delta_conc = timeseries['delta_dissolved']
    delta_emis = timeseries['delta_emission']
    
    # Calculate seasonal means
    time_array = flux_data['time_array']
    time_pd = pd.DatetimeIndex(time_array)
    
    # Summer mask: June 21 - September 22
    summer_mask = ((time_pd.month > 6) | ((time_pd.month == 6) & (time_pd.day >= 21))) & \
                  ((time_pd.month < 9) | ((time_pd.month == 9) & (time_pd.day <= 22)))
    
    # Winter mask: December 21 - March 20
    winter_mask = ((time_pd.month == 12) & (time_pd.day >= 21)) | \
                  (time_pd.month <= 2) | \
                  ((time_pd.month == 3) & (time_pd.day <= 20))
    
    # Calculate model predictions
    model_predictions = {}
    individual_costs = {}
    total_cost = 0.0
    
    # Dissolved summer
    if target_toggles.get('dissolved_summer', False):
        valid_data = delta_conc[summer_mask]
        valid_data = valid_data[~np.isnan(valid_data)]
        model_val = np.mean(valid_data) if len(valid_data) > 0 else np.nan
        obs_mean = obs_targets['dissolved_summer']['mean']
        obs_std = obs_targets['dissolved_summer']['std']
        error = _target_error_or_nan(model_val, obs_mean, obs_std, normalize_by_std)
        individual_costs['dissolved_summer'] = error
        model_predictions['dissolved_summer'] = model_val
        if np.isfinite(error):
            total_cost += error
    
    # Dissolved winter
    if target_toggles.get('dissolved_winter', False):
        valid_data = delta_conc[winter_mask]
        valid_data = valid_data[~np.isnan(valid_data)]
        model_val = np.mean(valid_data) if len(valid_data) > 0 else np.nan
        obs_mean = obs_targets['dissolved_winter']['mean']
        obs_std = obs_targets['dissolved_winter']['std']
        error = _target_error_or_nan(model_val, obs_mean, obs_std, normalize_by_std)
        individual_costs['dissolved_winter'] = error
        model_predictions['dissolved_winter'] = model_val
        if np.isfinite(error):
            total_cost += error
    
    # Emissions summer
    if target_toggles.get('emissions_summer', False):
        valid_data = delta_emis[summer_mask]
        valid_data = valid_data[~np.isnan(valid_data)]
        model_val = np.mean(valid_data) if len(valid_data) > 0 else np.nan
        obs_mean = obs_targets['emissions_summer']['mean']
        obs_std = obs_targets['emissions_summer']['std']
        error = _target_error_or_nan(model_val, obs_mean, obs_std, normalize_by_std)
        individual_costs['emissions_summer'] = error
        model_predictions['emissions_summer'] = model_val
        if np.isfinite(error):
            total_cost += error
    
    # Emissions winter
    if target_toggles.get('emissions_winter', False):
        valid_data = delta_emis[winter_mask]
        valid_data = valid_data[~np.isnan(valid_data)]
        model_val = np.mean(valid_data) if len(valid_data) > 0 else np.nan
        obs_mean = obs_targets['emissions_winter']['mean']
        obs_std = obs_targets['emissions_winter']['std']
        error = _target_error_or_nan(model_val, obs_mean, obs_std, normalize_by_std)
        individual_costs['emissions_winter'] = error
        model_predictions['emissions_winter'] = model_val
        if np.isfinite(error):
            total_cost += error
    
    return total_cost, individual_costs, model_predictions


# =============================================================================
# GPU-Accelerated Batch Cost Calculation
# =============================================================================

def calculate_optimal_gpu_batch_size(
    flux_data: Dict[str, np.ndarray],
    safety_factor: float = 0.7
) -> int:
    """
    Calculate optimal GPU batch size based on available memory and data size.
    
    Parameters:
    -----------
    flux_data : dict
        Flux data dictionary to estimate memory requirements
    safety_factor : float
        Fraction of free GPU memory to use (default: 0.7 = 70%)
    
    Returns:
    --------
    int : Recommended batch size (power of 2)
    """
    if not GPU_AVAILABLE:
        return 1
    
    try:
        device = cp.cuda.Device(0)
        free_memory = device.mem_info[0]  # bytes
        
        # Estimate number of timesteps from flux data
        n_timesteps = len(flux_data.get('sedch4df_data', []))
        if n_timesteps == 0:
            n_timesteps = 1461  # Default for 4-year run
        
        # Estimate memory per batch element (more accurate calculation):
        # - 8 float64 flux arrays: 8 * n_timesteps * 8 bytes
        # - 3 float64 state arrays (ch4_conc_mb, delta_conc_mb, isotope_mass_mb): 3 * n_timesteps * 8 bytes
        # - 2 float64 emission arrays (delta_emission, R_emission): 2 * n_timesteps * 8 bytes
        # - Intermediate arrays (R_dflux_arr, masks, etc.): ~2 * n_timesteps * 8 bytes
        # - Total per element: (8 + 3 + 2 + 2) * n_timesteps * 8 bytes
        bytes_per_element = 15 * n_timesteps * 8
        
        # Add 20% overhead for temporary allocations and kernel workspace
        bytes_per_element = int(bytes_per_element * 1.2)
        
        # Calculate batch size using safety factor
        usable_memory = int(free_memory * safety_factor)
        
        if bytes_per_element <= 0:
            return 256  # Safe default if calculation fails
        
        optimal_batch = usable_memory // bytes_per_element
        
        # Clamp to reasonable range
        # - Minimum: 32 (enough to hide latency)
        # - Maximum: 2048 (diminishing returns beyond this)
        optimal_batch = max(32, min(optimal_batch, 2048))
        
        # Round to power of 2 for optimal GPU performance
        import math
        if optimal_batch <= 0:
            return 256  # Safe default
        optimal_batch = 2 ** int(math.log2(optimal_batch))
        
        # Ensure we have at least 32
        return max(32, optimal_batch)
        
    except (RuntimeError, AttributeError, ValueError, ZeroDivisionError) as e:
        # More specific exception handling
        return 256  # Safe default

# CUDA kernel for parallel batch processing
# Each thread processes one batch element (one optimization run)
_batch_isotope_kernel = None

def _get_batch_isotope_kernel():
    """Get or create the CUDA kernel for parallel batch processing."""
    global _batch_isotope_kernel
    # Force recompilation by clearing cache (in case kernel code changed)
    _batch_isotope_kernel = None
    if _batch_isotope_kernel is None:
        kernel_code = '''
        // CUDA kernel - no includes needed, uses built-in math functions
        
        // Constants
        #define DT_SECONDS 86400.0
        #define MIN_CONC 1e-10
        #define MIN_SURFACE_FLUX 1e-20
        #define MIN_VALID_CONC 1e-12
        #define DELTA_MIN -100.0
        #define DELTA_MAX 10.0
        #define C13_TO_RATIO_DIVISOR 1000.0
            
        // Helper function to calculate seasonal mean
        __device__ double calc_seasonal_mean(
            const double* data, const bool* mask, int b_offset, int n_steps)
        {
            double sum = 0.0;
            int count = 0;
            for (int i = 0; i < n_steps; i++) {
                int idx = b_offset + i;
                if (mask[i] && !isnan(data[idx])) {
                    sum += data[idx];
                    count++;
                }
            }
            // Return NaN as 0.0/0.0 if no valid data (CUDA doesn't have NAN constant)
            return (count > 0) ? (sum / count) : (0.0/0.0);
        }
        
        // Helper function to calculate cost contribution
        __device__ double calc_cost_contribution(
            double model_val, double obs_mean, double obs_std, bool normalize_by_std)
        {
            if (isnan(model_val)) return 0.0;
            double diff = model_val - obs_mean;
            if (normalize_by_std) {
                double normalized_diff = diff / obs_std;
                return normalized_diff * normalized_diff;  // x*x instead of pow(x, 2.0)
            } else {
                return diff * diff;  // x*x instead of pow(x, 2.0)
            }
        }
        
        extern "C" __global__
        void batch_isotope_kernel(
            const double* params,           // [batch_size, 6] parameters
            const double* sedch4df,         // Flux arrays (shared across all batches)
            const double* och4prod,
            const double* ch4exc,
            const double* och4,
            const double* dflux,
            const double* ch4conc,
            const double* eflux,
            const double* ch4upb,
            double* ch4_conc_mb,            // [batch_size, n_steps] output arrays
            double* delta_conc_mb,
            double* isotope_mass_mb,
            double* delta_emission,         // [batch_size, n_steps]
            const bool* summer_mask,         // Seasonal masks (shared)
            const bool* winter_mask,
            double* costs,                  // Output costs [batch_size]
            double R_STD,
            double alpha_td,
            double obs_diss_summer_mean, double obs_diss_summer_std,
            double obs_diss_winter_mean, double obs_diss_winter_std,
            double obs_emis_summer_mean, double obs_emis_summer_std,
            double obs_emis_winter_mean, double obs_emis_winter_std,
            bool use_diss_summer, bool use_diss_winter,
            bool use_emis_summer, bool use_emis_winter,
            bool normalize_by_std,
            int n_steps, int batch_size)
        {
            int b = blockIdx.x * blockDim.x + threadIdx.x;
            if (b >= batch_size) return;
            
            // Bounds checking
            if (b < 0 || b >= batch_size) {
                costs[b] = 1e10;  // Error value
                return;
            }
            
            // Calculate array offsets for this batch element
            int b_offset = b * n_steps;
            
            // Bounds checking for array access
            if (b_offset < 0 || b_offset + n_steps > batch_size * n_steps) {
                costs[b] = 1e10;  // Error value
                return;
            }
            
            // Get parameters for this batch element with bounds checking
            int param_base = b * 6;
            if (param_base + 5 >= batch_size * 6) {
                costs[b] = 1e10;  // Error value
                return;
            }
            
            double alpha_am = params[param_base + 0];
            double alpha_hm = params[param_base + 1];
            double alpha_mo = params[param_base + 2];
            double alpha_op = params[param_base + 3];
            double f_am = params[param_base + 4];
            double C13_POM = params[param_base + 5];
            
            // Validate parameters
            if (alpha_hm <= 0 || alpha_am <= 0 || alpha_mo <= 0 || alpha_op <= 0 || 
                f_am < 0 || f_am > 1) {
                costs[b] = 1e10;  // Error value
                return;
            }
            
            double f_hm = 1.0 - f_am;
            
            // Calculate source ratios
            double R_POM = R_STD * (1.0 + C13_POM / C13_TO_RATIO_DIVISOR);
            double R_HM = R_POM / alpha_hm;
            double R_AM = R_POM / alpha_am;
            double R_sed_prod = f_hm * R_HM + f_am * R_AM;
            double C13_sed_prod = (R_sed_prod / R_STD - 1.0) * C13_TO_RATIO_DIVISOR;
            
            // Initial conditions
            ch4_conc_mb[b_offset] = fmax(ch4conc[0], MIN_CONC);
            delta_conc_mb[b_offset] = C13_sed_prod;
            isotope_mass_mb[b_offset] = ch4_conc_mb[b_offset] * R_STD * 
                                       (1.0 + delta_conc_mb[b_offset] / C13_TO_RATIO_DIVISOR);
            
            // Time stepping loop (sequential within each thread)
            // Also calculates surface emission isotopes during the loop to avoid separate pass
            for (int i = 1; i < n_steps; i++) {
                int idx = b_offset + i;
                int prev_idx = b_offset + i - 1;
                
                // Bounds checking
                if (idx >= batch_size * n_steps || prev_idx < 0) break;
                
                // Get previous values
                double C_prev = ch4_conc_mb[prev_idx];
                double delta_prev = delta_conc_mb[prev_idx];
                double isotope_mass_prev = isotope_mass_mb[prev_idx];
                
                // Calculate R_prev
                double R_prev;
                if (isnan(delta_prev)) {
                    R_prev = R_sed_prod;
                } else {
                    R_prev = R_STD * (1.0 + delta_prev / C13_TO_RATIO_DIVISOR);
                }
                
                // Fluxes with bounds checking
                if (i >= n_steps) break;
                double sedch4df_flux = sedch4df[i];
                double och4prod_flux = och4prod[i];
                double ch4exc_flux = ch4exc[i];
                double och4_flux = -och4[i];
                double dflux_val = -dflux[i];
                
                // Isotope ratios
                double R_och4prod = R_POM / alpha_op;
                double R_och4 = R_prev / alpha_mo;
                double R_dflux = R_prev / alpha_td;
                
                // Mass balance
                double dC_dt = sedch4df_flux + och4prod_flux + ch4exc_flux + och4_flux + dflux_val;
                double isotope_flux = (sedch4df_flux * R_sed_prod + 
                                     och4prod_flux * R_och4prod +
                                     ch4exc_flux * R_prev +
                                     och4_flux * R_och4 +
                                     dflux_val * R_dflux);
                
                double C_new = C_prev + dC_dt * DT_SECONDS;
                C_new = fmax(C_new, 0.0);
                ch4_conc_mb[idx] = C_new;
                
                // Isotope mass balance
                double isotope_mass_new = (C_prev > 0 ? isotope_mass_prev : 0.0) + isotope_flux * DT_SECONDS;
                isotope_mass_mb[idx] = isotope_mass_new;
                
                // Calculate delta
                bool C_new_valid = (C_new > MIN_VALID_CONC) && (isotope_mass_new != 0);
                double R_new = C_new_valid ? (isotope_mass_new / C_new) : R_sed_prod;
                double delta_new;
                if (R_new > 0 && C_new_valid) {
                    delta_new = (R_new / R_STD - 1.0) * C13_TO_RATIO_DIVISOR;
                    delta_new = fmax(fmin(delta_new, DELTA_MAX), DELTA_MIN);
                } else {
                    delta_new = C13_sed_prod;
                }
                if (C_new <= 0) {
                    isotope_mass_mb[idx] = 0.0;
                }
                
                delta_conc_mb[idx] = delta_new;
                
                // Calculate surface emission isotopes (merged into time-stepping loop)
                double total_surface_flux = dflux[i] + eflux[i] + ch4upb[i];
                
                if (total_surface_flux > MIN_SURFACE_FLUX) {
                    double R_dflux_val = isnan(delta_conc_mb[idx]) ? R_sed_prod : 
                                        (R_STD * (1.0 + delta_conc_mb[idx] / C13_TO_RATIO_DIVISOR) / alpha_td);
                    double isotope_flux_dflux = dflux[i] * R_dflux_val;
                    double isotope_flux_eflux = eflux[i] * R_sed_prod;
                    double isotope_flux_upb = ch4upb[i] * R_sed_prod;
                    double total_isotope_flux = isotope_flux_dflux + isotope_flux_eflux + isotope_flux_upb;
                    // Division by zero protection
                    if (total_surface_flux > 0) {
                        double R_emission = total_isotope_flux / total_surface_flux;
                        delta_emission[idx] = (R_emission / R_STD - 1.0) * C13_TO_RATIO_DIVISOR;
                    } else {
                        delta_emission[idx] = C13_sed_prod;
                    }
                } else {
                    delta_emission[idx] = C13_sed_prod;
                }
            }
            
            // Handle initial timestep emission calculation
            int idx_0 = b_offset;
            double total_surface_flux_0 = dflux[0] + eflux[0] + ch4upb[0];
            if (total_surface_flux_0 > MIN_SURFACE_FLUX) {
                double R_dflux_val_0 = isnan(delta_conc_mb[idx_0]) ? R_sed_prod : 
                                      (R_STD * (1.0 + delta_conc_mb[idx_0] / C13_TO_RATIO_DIVISOR) / alpha_td);
                double isotope_flux_dflux_0 = dflux[0] * R_dflux_val_0;
                double isotope_flux_eflux_0 = eflux[0] * R_sed_prod;
                double isotope_flux_upb_0 = ch4upb[0] * R_sed_prod;
                double total_isotope_flux_0 = isotope_flux_dflux_0 + isotope_flux_eflux_0 + isotope_flux_upb_0;
                if (total_surface_flux_0 > 0) {
                    double R_emission_0 = total_isotope_flux_0 / total_surface_flux_0;
                    delta_emission[idx_0] = (R_emission_0 / R_STD - 1.0) * C13_TO_RATIO_DIVISOR;
                } else {
                    delta_emission[idx_0] = C13_sed_prod;
                }
            } else {
                delta_emission[idx_0] = C13_sed_prod;
            }
            
            // Calculate seasonal means and cost using helper functions
            double cost = 0.0;
            
            // Dissolved summer
            if (use_diss_summer) {
                double model_val = calc_seasonal_mean(delta_conc_mb, summer_mask, b_offset, n_steps);
                cost += calc_cost_contribution(model_val, obs_diss_summer_mean, obs_diss_summer_std, normalize_by_std);
            }
            
            // Dissolved winter
            if (use_diss_winter) {
                double model_val = calc_seasonal_mean(delta_conc_mb, winter_mask, b_offset, n_steps);
                cost += calc_cost_contribution(model_val, obs_diss_winter_mean, obs_diss_winter_std, normalize_by_std);
            }
            
            // Emissions summer
            if (use_emis_summer) {
                double model_val = calc_seasonal_mean(delta_emission, summer_mask, b_offset, n_steps);
                cost += calc_cost_contribution(model_val, obs_emis_summer_mean, obs_emis_summer_std, normalize_by_std);
            }
            
            // Emissions winter
            if (use_emis_winter) {
                double model_val = calc_seasonal_mean(delta_emission, winter_mask, b_offset, n_steps);
                cost += calc_cost_contribution(model_val, obs_emis_winter_mean, obs_emis_winter_std, normalize_by_std);
            }
            
            costs[b] = cost;
        }
        '''
        _batch_isotope_kernel = cp.RawKernel(kernel_code, 'batch_isotope_kernel')
    return _batch_isotope_kernel


def _batch_calculate_cost_gpu(
    params_batch: np.ndarray,
    flux_data: Dict[str, np.ndarray],
    obs_targets: Dict[str, Dict[str, float]],
    target_toggles: Dict[str, bool],
    normalize_by_std: bool,
    batch_size: int
) -> np.ndarray:
    """
    Calculate cost for a batch of parameter sets using GPU-accelerated isotope computation.
    Uses parallel CUDA kernel where each GPU thread processes one batch element.
    
    Parameters:
    -----------
    params_batch : np.ndarray
        Array of shape (batch_size, 6) with parameter sets
        [alpha_am, alpha_hm, alpha_mo, alpha_op, f_am, C13_POM]
    flux_data : dict
        Flux data dictionary from ALBM
    obs_targets : dict
        Observational targets dictionary
    target_toggles : dict
        Boolean toggles for which targets to include
    normalize_by_std : bool
        Whether to normalize errors by standard deviation
    batch_size : int
        Number of parameter sets in batch
    
    Returns:
    --------
    np.ndarray : Array of shape (batch_size,) with cost values
    """
    if not GPU_AVAILABLE:
        raise RuntimeError("GPU not available")
    
    # Verify GPU is accessible
    try:
        device = cp.cuda.Device(0)
        device.use()
    except Exception as e:
        raise RuntimeError(f"Failed to access GPU device: {e}")
    
    # Transfer flux data to GPU once
    sedch4df = cp.asarray(flux_data['sedch4df_data'])
    och4prod = cp.asarray(flux_data['och4prod_data'])
    ch4exc = cp.asarray(flux_data['ch4exc_data'])
    och4 = cp.asarray(flux_data['och4_data'])
    dflux = cp.asarray(flux_data['dflux_data'])
    ch4conc = cp.asarray(flux_data['ch4conc_data'])
    eflux = cp.asarray(flux_data['eflux_data'])
    ch4upb = cp.asarray(flux_data['ch4upb_data'])
    time_array = flux_data['time_array']
    
    n_steps = len(sedch4df)
    params_gpu = cp.asarray(params_batch.flatten(), dtype=cp.float64)  # Flatten for kernel
    
    # Pre-compute seasonal masks on CPU once (pandas operation)
    time_pd = pd.DatetimeIndex(time_array)
    summer_mask = ((time_pd.month > 6) | ((time_pd.month == 6) & (time_pd.day >= 21))) & \
                  ((time_pd.month < 9) | ((time_pd.month == 9) & (time_pd.day <= 22)))
    winter_mask = ((time_pd.month == 12) & (time_pd.day >= 21)) | \
                  (time_pd.month <= 2) | \
                  ((time_pd.month == 3) & (time_pd.day <= 20))
    summer_mask_gpu = cp.asarray(summer_mask, dtype=cp.bool_)
    winter_mask_gpu = cp.asarray(winter_mask, dtype=cp.bool_)
    
    # Pre-allocate arrays for all batch elements [batch_size, n_steps]
    # Flatten to 1D for kernel (kernel uses b_offset = b * n_steps)
    # Use empty() instead of zeros() since kernel initializes all values
    ch4_conc_mb = cp.empty(batch_size * n_steps, dtype=cp.float64)
    delta_conc_mb = cp.empty(batch_size * n_steps, dtype=cp.float64)
    isotope_mass_mb = cp.empty(batch_size * n_steps, dtype=cp.float64)
    delta_emission = cp.empty(batch_size * n_steps, dtype=cp.float64)
    costs_gpu = cp.zeros(batch_size, dtype=cp.float64)  # Keep zeros for costs (error value)
    
    # Pre-compute observation values
    obs_diss_summer_mean = float(obs_targets.get('dissolved_summer', {}).get('mean', 0.0))
    obs_diss_summer_std = float(obs_targets.get('dissolved_summer', {}).get('std', 1.0))
    obs_diss_winter_mean = float(obs_targets.get('dissolved_winter', {}).get('mean', 0.0))
    obs_diss_winter_std = float(obs_targets.get('dissolved_winter', {}).get('std', 1.0))
    obs_emis_summer_mean = float(obs_targets.get('emissions_summer', {}).get('mean', 0.0))
    obs_emis_summer_std = float(obs_targets.get('emissions_summer', {}).get('std', 1.0))
    obs_emis_winter_mean = float(obs_targets.get('emissions_winter', {}).get('mean', 0.0))
    obs_emis_winter_std = float(obs_targets.get('emissions_winter', {}).get('std', 1.0))
    
    # Thread-safe state tracking
    _gpu_batch_call_count = getattr(_batch_calculate_cost_gpu, '_call_count', 0)
    _batch_calculate_cost_gpu._call_count = _gpu_batch_call_count + 1
    is_first_call = (_gpu_batch_call_count == 0)
    
    if is_first_call:
        print(f"[{time.strftime('%H:%M:%S')}] GPU parallel batch kernel: batch_size={batch_size}")
        print(f"[{time.strftime('%H:%M:%S')}] Processing all batch elements in parallel on GPU...")
    
    # Get kernel
    kernel = _get_batch_isotope_kernel()
    
    # Launch kernel: one thread per batch element
    # Thread block size configurable (optimal depends on GPU architecture)
    # Default to 256, but can be tuned for specific GPUs
    threads_per_block = getattr(_batch_calculate_cost_gpu, '_threads_per_block', 256)
    # Auto-tune for RTX 4000 Ada (compute capability 8.9) - prefer 256 or 512
    if not hasattr(_batch_calculate_cost_gpu, '_threads_per_block'):
        # Try to get optimal block size from GPU
        try:
            device = cp.cuda.Device(0)
            # For compute capability 8.x, 256-512 threads per block is optimal
            threads_per_block = 256
            _batch_calculate_cost_gpu._threads_per_block = threads_per_block
        except Exception:
            threads_per_block = 256
    
    blocks = (batch_size + threads_per_block - 1) // threads_per_block
    
    # Progress indication: record start time
    import time as time_module
    start_time = time_module.time() if is_first_call else None
    
    kernel(
        (blocks,), (threads_per_block,),
        (params_gpu, sedch4df, och4prod, ch4exc, och4, dflux, ch4conc, eflux, ch4upb,
         ch4_conc_mb, delta_conc_mb, isotope_mass_mb, delta_emission,
         summer_mask_gpu, winter_mask_gpu, costs_gpu,
         float(R_STD), float(DEFAULT_ALPHA_TD),
         obs_diss_summer_mean, obs_diss_summer_std,
         obs_diss_winter_mean, obs_diss_winter_std,
         obs_emis_summer_mean, obs_emis_summer_std,
         obs_emis_winter_mean, obs_emis_winter_std,
         bool(target_toggles.get('dissolved_summer', False)),
         bool(target_toggles.get('dissolved_winter', False)),
         bool(target_toggles.get('emissions_summer', False)),
         bool(target_toggles.get('emissions_winter', False)),
         bool(normalize_by_std),
         n_steps, batch_size)
    )
    
    # Synchronize GPU to ensure all computations complete
    cp.cuda.Stream.null.synchronize()
    
    # Progress indication: report timing
    if is_first_call and start_time is not None:
        elapsed = time_module.time() - start_time
        print(f"[{time.strftime('%H:%M:%S')}] GPU parallel batch computation complete in {elapsed:.2f}s, transferring results to CPU...")
        print(f"[{time.strftime('%H:%M:%S')}] Average time per batch element: {elapsed/batch_size*1000:.2f}ms")
    
    # Transfer results back to CPU
    costs_cpu = cp.asnumpy(costs_gpu)
    
    return costs_cpu


# =============================================================================
# Batched Objective Function Wrapper for GPU
# =============================================================================

class BatchedObjectiveWrapper:
    """
    Wrapper that batches objective function calls for GPU evaluation.
    Thread-safe for use with multiprocessing.
    """
    def __init__(self, base_func, flux_data, obs_targets, target_toggles, 
                 normalize_by_std, batch_size=32, use_gpu=False):
        self.base_func = base_func
        self.flux_data = flux_data
        self.obs_targets = obs_targets
        self.target_toggles = target_toggles
        self.normalize_by_std = normalize_by_std
        self.batch_size = batch_size
        self.use_gpu = use_gpu and GPU_AVAILABLE
        
        # Thread-safe batching
        self.lock = threading.Lock()
        self.pending_params = []
        self.pending_indices = []
        self.results_cache = OrderedDict()
        self.cache_max_size = 10000  # Limit cache size
        self.gpu_call_count = 0  # Track GPU usage
        self.cpu_fallback_count = 0  # Track CPU fallbacks
        
    def _flush_batch(self):
        """Evaluate pending batch and cache results."""
        if not self.pending_params:
            return
        
        params_array = np.array(self.pending_params)
        n_pending = len(self.pending_params)
        
        if self.use_gpu:
            try:
                # Force GPU synchronization to ensure computation happens
                costs = _batch_calculate_cost_gpu(
                    params_array, self.flux_data, self.obs_targets,
                    self.target_toggles, self.normalize_by_std, n_pending
                )
                # Synchronize GPU to ensure computation completes
                if GPU_AVAILABLE:
                    cp.cuda.Stream.null.synchronize()
                self.gpu_call_count += n_pending
            except Exception as e:
                # Fallback to CPU if GPU fails
                if not hasattr(self, '_gpu_batch_error_logged'):
                    import traceback
                    print(f"[{time.strftime('%H:%M:%S')}] GPU batch evaluation failed: {e}")
                    traceback.print_exc()
                    self._gpu_batch_error_logged = True
                self.cpu_fallback_count += n_pending
                costs = np.array([self.base_func(p) for p in self.pending_params])
        else:
            costs = np.array([self.base_func(p) for p in self.pending_params])
        
        # Cache results - ensure all results are cached before clearing pending
        for param_tuple, cost in zip(self.pending_indices, costs):
            self.results_cache[param_tuple] = float(cost)  # Ensure float type
            if len(self.results_cache) > self.cache_max_size:
                # Remove oldest entry
                self.results_cache.popitem(last=False)
        
        # Clear pending AFTER caching (important for thread safety)
        pending_copy = list(self.pending_indices)  # Keep copy for debugging if needed
        self.pending_params = []
        self.pending_indices = []
    
    def __call__(self, x):
        """Evaluate objective function, batching if using GPU."""
        param_tuple = tuple(x)
        
        # Check cache first
        with self.lock:
            if param_tuple in self.results_cache:
                return self.results_cache[param_tuple]
        
        # For GPU, accumulate calls and batch them
        if self.use_gpu:
            with self.lock:
                # Add to pending batch
                self.pending_params.append(x)
                self.pending_indices.append(param_tuple)
                n_pending = len(self.pending_params)
                
                # Print first GPU call only once (use a flag to ensure it only prints once)
                if not hasattr(self, '_batching_message_printed'):
                    print(f"[{time.strftime('%H:%M:%S')}] GPU batching enabled: accumulating calls for batch_size={self.batch_size}")
                    self._batching_message_printed = True
                
                # Flush batch if it's full
                if n_pending >= self.batch_size:
                    # Make a copy of pending_indices to check after flush
                    pending_before_flush = list(self.pending_indices)
                    self._flush_batch()
                    # Result should definitely be in cache now (it was in the batch we just flushed)
                    if param_tuple in self.results_cache:
                        return self.results_cache[param_tuple]
                    # If not in cache, something went wrong - but continue to wait/retry logic below
            
            # Since differential_evolution calls one at a time, flush when we have a reasonable number
            # This balances batching efficiency with responsiveness
            import time as time_module
            min_flush_size = max(1, min(16, self.batch_size // 128))  # Flush with 1-16 calls
            start_wait = time_module.time()
            max_wait = 0.05  # Maximum 50ms wait (reduced for responsiveness)
            
            # Wait and periodically check/flush
            while time_module.time() - start_wait < max_wait:
                with self.lock:
                    # Check if result is already cached (from a previous flush)
                    if param_tuple in self.results_cache:
                        return self.results_cache[param_tuple]
                    n_pending = len(self.pending_params)
                    # Flush if we have minimum batch size
                    if n_pending >= min_flush_size:
                        self._flush_batch()
                        # Result should definitely be in cache after flush
                        if param_tuple in self.results_cache:
                            return self.results_cache[param_tuple]
                time_module.sleep(0.005)  # Small sleep to avoid busy-wait
            
            # Final flush - force flush whatever we have
            with self.lock:
                # One more check before flushing
                if param_tuple in self.results_cache:
                    return self.results_cache[param_tuple]
                # Force flush any remaining pending
                if len(self.pending_params) > 0:
                    self._flush_batch()
                    # Result should be in cache now
                    if param_tuple in self.results_cache:
                        return self.results_cache[param_tuple]
                
                # If still not in cache, something went wrong - evaluate on GPU individually
                # This should be very rare
                try:
                    params_array = np.array([x])
                    costs = _batch_calculate_cost_gpu(
                        params_array, self.flux_data, self.obs_targets,
                        self.target_toggles, self.normalize_by_std, 1
                    )
                    result = float(costs[0])
                    if GPU_AVAILABLE:
                        cp.cuda.Stream.null.synchronize()
                    self.gpu_call_count += 1
                    self.results_cache[param_tuple] = result
                    if len(self.results_cache) > self.cache_max_size:
                        self.results_cache.popitem(last=False)
                    return result
                except Exception as e:
                    # Last resort: CPU fallback
                    if not hasattr(self, '_cpu_fallback_logged'):
                        print(f"[{time.strftime('%H:%M:%S')}] WARNING: GPU evaluation failed, using CPU fallback: {e}")
                        self._cpu_fallback_logged = True
                    self.cpu_fallback_count += 1
                    result = self.base_func(x)
                    self.results_cache[param_tuple] = result
                    if len(self.results_cache) > self.cache_max_size:
                        self.results_cache.popitem(last=False)
                    return result
        else:
            # CPU mode: evaluate directly
            result = self.base_func(x)
            with self.lock:
                self.results_cache[param_tuple] = result
                if len(self.results_cache) > self.cache_max_size:
                    self.results_cache.popitem(last=False)
            return result
    
    def finalize(self):
        """Flush any remaining pending evaluations."""
        with self.lock:
            self._flush_batch()
        # Print GPU usage statistics (always print when GPU is enabled)
        if self.use_gpu:
            total = self.gpu_call_count + self.cpu_fallback_count
            if total > 0:
                gpu_pct = (self.gpu_call_count / total * 100)
                print(f"[{time.strftime('%H:%M:%S')}] GPU usage stats: {self.gpu_call_count} GPU calls, {self.cpu_fallback_count} CPU fallbacks ({gpu_pct:.1f}% GPU)")
            else:
                print(f"[{time.strftime('%H:%M:%S')}] WARNING: No GPU calls made - check if GPU is actually being used!")


# =============================================================================
# Objective functions for optimization (picklable for multiprocessing)
# =============================================================================

# Legacy globals for single-worker / GPU path
_opt_flux_data = None
_opt_obs_targets = None
_opt_target_toggles = None
_opt_normalize_by_std = False

# Penalty cost when evaluation fails or is non-finite
_COST_PENALTY = 1e10


def _objective_with_args(
    x: np.ndarray,
    flux_data: Dict[str, np.ndarray],
    obs_targets: Dict[str, Dict[str, float]],
    target_toggles: Dict[str, bool],
    normalize_by_std: bool
) -> float:
    """
    Module-level objective for differential_evolution. Pass data via args=()
    so it is picklable when using workers > 1.
    """
    params = {
        'alpha_am': float(x[0]),
        'alpha_hm': float(x[1]),
        'alpha_mo': float(x[2]),
        'alpha_op': float(x[3]),
        'f_am': float(x[4]),
        'C13_POM': float(x[5])
    }
    cost, _, _ = calculate_cost(
        params, flux_data, obs_targets, target_toggles, normalize_by_std
    )
    return float(cost) if np.isfinite(cost) else _COST_PENALTY


def _objective_function_global(x: np.ndarray) -> float:
    """
    Objective using module globals (only works in main process; use
    _make_objective_with_data when workers > 1).
    """
    global _opt_flux_data, _opt_obs_targets, _opt_target_toggles, _opt_normalize_by_std
    params = {
        'alpha_am': x[0], 'alpha_hm': x[1], 'alpha_mo': x[2],
        'alpha_op': x[3], 'f_am': x[4], 'C13_POM': x[5]
    }
    cost, _, _ = calculate_cost(
        params, _opt_flux_data, _opt_obs_targets,
        _opt_target_toggles, _opt_normalize_by_std
    )
    return float(cost) if np.isfinite(cost) else _COST_PENALTY


# =============================================================================
# Optimization Function
# =============================================================================

@dataclass
class OptimizationResult:
    """Container for optimization results."""
    alpha_am: float
    alpha_hm: float
    alpha_mo: float
    alpha_op: float
    f_am: float
    f_hm: float
    C13_POM: float
    cost: float
    success: bool
    n_evaluations: int


@dataclass
class OptimizationResultTemp:
    """Container for temperature-based optimization results."""
    m: float  # Slope (°C⁻¹)
    b: float  # Intercept (‰)
    alpha_mo: float
    alpha_op: float
    cost: float
    success: bool
    n_evaluations: int


def get_default_bounds() -> Dict[str, Tuple[float, float]]:
    """Get default parameter bounds for optimization."""
    return {
        'alpha_am': (1.000, 1.040),
        'alpha_hm': (1.030, 1.080),
        'alpha_mo': (1.015, 1.035),
        'alpha_op': (1.000, 1.080),
        'f_am': (0.0, 1.0),
        'C13_POM': (-28.0, -22.0)
    }


def get_default_bounds_temp() -> Dict[str, Tuple[float, float]]:
    """Get default parameter bounds for temperature-based optimization."""
    return {
        'm': (-10.0, 10.0),  # Slope (°C⁻¹)
        'b': (-100.0, -20.0),  # Intercept (‰)
        'alpha_mo': (1.015, 1.035),
        'alpha_op': (1.000, 1.080)
    }


def run_optimization(
    flux_data: Dict[str, np.ndarray],
    obs_targets: Dict[str, Dict[str, float]],
    target_toggles: Dict[str, bool],
    bounds: Optional[Dict[str, Tuple[float, float]]] = None,
    normalize_by_std: bool = False,
    maxiter: int = 200,
    tol: float = 0.01,
    popsize: int = 15,
    mutation: Tuple[float, float] = (0.5, 1.0),
    recombination: float = 0.7,
    polish: bool = True,
    workers: int = -1,
    verbose: bool = True,
    seed: Optional[int] = None,
    use_gpu: bool = False,
    gpu_batch_size: Optional[int] = None
) -> OptimizationResult:
    """
    Run differential evolution optimization for isotope fractionation factors.
    
    Parameters:
    -----------
    flux_data : dict
        Flux data dictionary from ALBM
    obs_targets : dict
        Observational targets dictionary
    target_toggles : dict
        Boolean toggles for which targets to include
    bounds : dict, optional
        Parameter bounds dictionary
    normalize_by_std : bool
        Whether to normalize errors by standard deviation
    maxiter : int
        Maximum number of iterations
    tol : float
        Tolerance for convergence
    popsize : int
        Population size multiplier
    mutation : tuple
        Mutation constant range
    recombination : float
        Recombination constant
    polish : bool
        Whether to polish with L-BFGS-B
    workers : int
        Number of parallel workers (-1 for all cores)
    verbose : bool
        Whether to print progress
    seed : int, optional
        Random seed for reproducibility
    use_gpu : bool
        Whether to use GPU acceleration (if available)
    gpu_batch_size : int, optional
        Batch size for GPU evaluation (only used if use_gpu=True).
        If None or 0, automatically calculates optimal batch size based on GPU memory.
    
    Returns:
    --------
    OptimizationResult : Optimization results
    """
    global _opt_flux_data, _opt_obs_targets, _opt_target_toggles, _opt_normalize_by_std
    
    # Check GPU availability
    using_gpu = use_gpu and GPU_AVAILABLE
    if use_gpu:
        if not GPU_AVAILABLE:
            print(f"[{time.strftime('%H:%M:%S')}] Warning: GPU requested but not available. Using CPU.")
            using_gpu = False
        else:
            # Always print GPU info when GPU is enabled (even if verbose=False)
            try:
                device = cp.cuda.Device(0)
                mem_info = device.mem_info
                print(f"[{time.strftime('%H:%M:%S')}] GPU enabled: {GPU_INFO.get('name', 'NVIDIA GPU')}")
                print(f"[{time.strftime('%H:%M:%S')}] GPU memory: {mem_info[0]/(1024**3):.1f} GB free / {mem_info[1]/(1024**3):.1f} GB total")
            except Exception as e:
                print(f"[{time.strftime('%H:%M:%S')}] Warning: Could not access GPU device: {e}")
                using_gpu = False
    
    # Calculate optimal GPU batch size if not specified
    if using_gpu:
        if gpu_batch_size is None or gpu_batch_size == 0:
            gpu_batch_size = calculate_optimal_gpu_batch_size(flux_data)
            print(f"[{time.strftime('%H:%M:%S')}] Auto-calculated optimal GPU batch size: {gpu_batch_size}")
    
    # Module-level objective + args is picklable for multiprocessing; data passed via args=()
    opt_args = (flux_data, obs_targets, target_toggles, normalize_by_std)
    global _opt_flux_data, _opt_obs_targets, _opt_target_toggles, _opt_normalize_by_std
    _opt_flux_data = flux_data
    _opt_obs_targets = obs_targets
    _opt_target_toggles = target_toggles
    _opt_normalize_by_std = normalize_by_std
    
    # Get bounds
    if bounds is None:
        bounds = get_default_bounds()
    
    bounds_list = [
        bounds['alpha_am'],
        bounds['alpha_hm'],
        bounds['alpha_mo'],
        bounds['alpha_op'],
        bounds['f_am'],
        bounds['C13_POM']
    ]
    
    # Determine number of workers (disable multiprocessing if using GPU)
    if using_gpu:
        n_workers = 1  # GPU doesn't work well with multiprocessing
        print(f"[{time.strftime('%H:%M:%S')}] Running optimization with GPU acceleration (batch_size={gpu_batch_size})...")
        # Test GPU with a simple operation to verify it's working
        try:
            test_array = cp.array([1.0, 2.0, 3.0])
            test_result = cp.sum(test_array ** 2)
            _ = float(test_result)  # Force computation
            cp.cuda.Stream.null.synchronize()  # Ensure computation completes
            print(f"[{time.strftime('%H:%M:%S')}] GPU test passed - GPU is ready for computation")
        except Exception as e:
            print(f"[{time.strftime('%H:%M:%S')}] WARNING: GPU test failed: {e}")
            import traceback
            traceback.print_exc()
            using_gpu = False
    else:
        if workers == -1:
            from multiprocessing import cpu_count
            n_workers = cpu_count()
        else:
            n_workers = workers
        if verbose:
            print(f"[{time.strftime('%H:%M:%S')}] Running optimization with {n_workers} CPU workers...")
    
    opt_start = time.time()
    
    # Objective: module-level func + args (picklable for workers); GPU wrapper needs (x)-only callable
    if using_gpu:
        base_for_gpu = lambda x: _objective_with_args(x, flux_data, obs_targets, target_toggles, normalize_by_std)
        batched_obj = BatchedObjectiveWrapper(
            base_for_gpu,
            flux_data,
            obs_targets,
            target_toggles,
            normalize_by_std,
            batch_size=gpu_batch_size,
            use_gpu=True
        )
        objective_func = batched_obj
        opt_args_use = None  # not used
    else:
        objective_func = _objective_with_args
        opt_args_use = opt_args
    
    # Run optimization (args= passes data to workers so objective is picklable)
    try:
        result = differential_evolution(
            objective_func,
            bounds=bounds_list,
            args=opt_args_use if opt_args_use is not None else (),
            maxiter=maxiter,
            tol=tol,
            popsize=popsize,
            mutation=mutation,
            recombination=recombination,
            polish=polish,
            workers=n_workers,
            updating='deferred' if n_workers != 1 else 'immediate',
            disp=verbose,
            seed=seed
        )
    finally:
        # Finalize batched wrapper if used
        if using_gpu:
            batched_obj.finalize()
            # Clean up GPU memory
            try:
                if GPU_AVAILABLE:
                    cp.get_default_memory_pool().free_all_blocks()
            except Exception:
                pass
    
    opt_elapsed = time.time() - opt_start
    if verbose:
        print(f"[{time.strftime('%H:%M:%S')}] Optimization complete: cost={result.fun:.4f}, {result.nfev} evaluations ({opt_elapsed:.1f}s)")
    
    return OptimizationResult(
        alpha_am=result.x[0],
        alpha_hm=result.x[1],
        alpha_mo=result.x[2],
        alpha_op=result.x[3],
        f_am=result.x[4],
        f_hm=1.0 - result.x[4],
        C13_POM=result.x[5],
        cost=result.fun,
        success=result.success,
        n_evaluations=result.nfev
    )


def get_default_params() -> Dict[str, float]:
    """Get default (initial guess) parameter values."""
    return {
        'alpha_am': DEFAULT_ALPHA_AM,
        'alpha_hm': DEFAULT_ALPHA_HM,
        'alpha_mo': DEFAULT_ALPHA_MO,
        'alpha_op': DEFAULT_ALPHA_OP,
        'f_am': DEFAULT_F_AM,
        'C13_POM': DEFAULT_C13_POM
    }


def get_default_params_temp() -> Dict[str, float]:
    """Get default (initial guess) parameter values for temperature-based optimization."""
    return {
        'm': 0.0,  # Slope (°C⁻¹)
        'b': -25.0,  # Intercept (‰)
        'alpha_mo': DEFAULT_ALPHA_MO,
        'alpha_op': DEFAULT_ALPHA_OP
    }


# Global variables for temperature-based optimization
_opt_temp_flux_data = None
_opt_temp_obs_targets = None
_opt_temp_target_toggles = None
_opt_temp_normalize_by_std = False


def _objective_temp_with_args(
    x: np.ndarray,
    flux_data: Dict[str, np.ndarray],
    obs_targets: Dict[str, Dict[str, float]],
    target_toggles: Dict[str, bool],
    normalize_by_std: bool
) -> float:
    """Module-level objective for temperature-based DE; pass data via args=() for pickling."""
    params = {
        'm': float(x[0]), 'b': float(x[1]),
        'alpha_mo': float(x[2]), 'alpha_op': float(x[3])
    }
    cost, _, _ = calculate_cost_temp(
        params, flux_data, obs_targets, target_toggles, normalize_by_std
    )
    return float(cost) if np.isfinite(cost) else _COST_PENALTY


def _objective_function_temp_global(x: np.ndarray) -> float:
    """Objective using module globals (use _make_objective_temp_with_data when workers > 1)."""
    global _opt_temp_flux_data, _opt_temp_obs_targets, _opt_temp_target_toggles, _opt_temp_normalize_by_std
    params = {'m': x[0], 'b': x[1], 'alpha_mo': x[2], 'alpha_op': x[3]}
    cost, _, _ = calculate_cost_temp(
        params, _opt_temp_flux_data, _opt_temp_obs_targets,
        _opt_temp_target_toggles, _opt_temp_normalize_by_std
    )
    return float(cost) if np.isfinite(cost) else _COST_PENALTY


def run_optimization_temp(
    flux_data: Dict[str, np.ndarray],
    obs_targets: Dict[str, Dict[str, float]],
    target_toggles: Dict[str, bool],
    bounds: Optional[Dict[str, Tuple[float, float]]] = None,
    normalize_by_std: bool = False,
    maxiter: int = 200,
    tol: float = 0.01,
    popsize: int = 15,
    mutation: Tuple[float, float] = (0.5, 1.0),
    recombination: float = 0.7,
    polish: bool = True,
    workers: int = -1,
    verbose: bool = True,
    seed: Optional[int] = None
) -> OptimizationResultTemp:
    """
    Run differential evolution optimization for temperature-based sediment production.
    Fits linear relationship: C13_sed_prod = m * temp_avg + b
    
    Parameters:
    -----------
    flux_data : dict
        Flux data dictionary from ALBM (must include 'sediment_temp_avg')
    obs_targets : dict
        Observational targets dictionary
    target_toggles : dict
        Boolean toggles for which targets to include
    bounds : dict, optional
        Parameter bounds dictionary (default: m, b, alpha_mo, alpha_op)
    normalize_by_std : bool
        Whether to normalize errors by standard deviation
    maxiter : int
        Maximum number of iterations
    tol : float
        Tolerance for convergence
    popsize : int
        Population size multiplier
    mutation : tuple
        Mutation constant range
    recombination : float
        Recombination constant
    polish : bool
        Whether to polish with L-BFGS-B
    workers : int
        Number of parallel workers (-1 for all cores)
    verbose : bool
        Whether to print progress
    seed : int, optional
        Random seed for reproducibility
    
    Returns:
    --------
    OptimizationResultTemp : Optimization results
    """
    global _opt_temp_flux_data, _opt_temp_obs_targets, _opt_temp_target_toggles, _opt_temp_normalize_by_std
    
    # Check that temperature data is available
    temp_avg = flux_data.get('sediment_temp_avg')
    if temp_avg is None:
        raise ValueError(
            "flux_data must include valid 'sediment_temp_avg' data for "
            "temperature-based optimization"
        )
    temp_avg = np.asarray(temp_avg, dtype=float)
    if temp_avg.ndim != 1 or temp_avg.size == 0 or not np.isfinite(temp_avg).all():
        raise ValueError(
            "flux_data must include a non-empty fully finite 'sediment_temp_avg' array "
            "for temperature-based optimization"
        )
    
    # Module-level objective + args is picklable for workers
    opt_args_temp = (flux_data, obs_targets, target_toggles, normalize_by_std)
    _opt_temp_flux_data = flux_data
    _opt_temp_obs_targets = obs_targets
    _opt_temp_target_toggles = target_toggles
    _opt_temp_normalize_by_std = normalize_by_std
    
    # Get bounds
    if bounds is None:
        bounds = get_default_bounds_temp()
    
    bounds_list = [
        bounds['m'],
        bounds['b'],
        bounds['alpha_mo'],
        bounds['alpha_op']
    ]
    
    # Determine number of workers
    if workers == -1:
        from multiprocessing import cpu_count
        n_workers = cpu_count()
    else:
        n_workers = workers
    if verbose:
        print(f"[{time.strftime('%H:%M:%S')}] Running temperature-based optimization with {n_workers} workers...")
        print(f"[{time.strftime('%H:%M:%S')}] Fitting: C13_sed_prod = m * temp_avg + b")
    
    opt_start = time.time()
    
    # Run optimization (args= passes data so module-level objective is picklable)
    result = differential_evolution(
        _objective_temp_with_args,
        bounds=bounds_list,
        args=opt_args_temp,
        maxiter=maxiter,
        tol=tol,
        popsize=popsize,
        mutation=mutation,
        recombination=recombination,
        polish=polish,
        workers=n_workers,
        updating='deferred' if n_workers != 1 else 'immediate',
        disp=verbose,
        seed=seed
    )
    
    opt_elapsed = time.time() - opt_start
    if verbose:
        print(f"[{time.strftime('%H:%M:%S')}] Temperature-based optimization complete: cost={result.fun:.4f}, {result.nfev} evaluations ({opt_elapsed:.1f}s)")
        print(f"[{time.strftime('%H:%M:%S')}] Fitted parameters: m={result.x[0]:.4f} °C⁻¹, b={result.x[1]:.2f} ‰")
    
    return OptimizationResultTemp(
        m=result.x[0],
        b=result.x[1],
        alpha_mo=result.x[2],
        alpha_op=result.x[3],
        cost=result.fun,
        success=result.success,
        n_evaluations=result.nfev
    )


# =============================================================================
# Main execution - demonstrates usage
# =============================================================================

if __name__ == '__main__':
    print("Isotope Model Library")
    print("=" * 50)
    print("\nThis module provides isotope modeling functions.")
    print("Import and use in other scripts:")
    print("  from isotope_model import run_optimization, compute_isotope_timeseries")
