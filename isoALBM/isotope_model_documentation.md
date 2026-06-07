# isoALBM Isotope Model Implementation

This document describes the current implementation of the methane stable-carbon
isotope workflow in this repository. It covers the data objects, isotope state
updates, optimization targets, configurable runner, and generated outputs.

The canonical command is:

```bash
python src/run_isoalbm.py --config configs/paper.toml
```

The legacy paper scripts remain available:

```bash
python src/manuscript_scripts/multi_condition_optimization.py
python src/manuscript_scripts/multi_start_optimization.py
```

The Python modules live in `src/`. The processed cache is a schema-versioned
payload that stores ALBM and eddy data fields without depending on dataclass
import paths; legacy class-pickle caches are still read for compatibility.

---

## 1. Main Modules

| File | Role |
| --- | --- |
| `src/run_isoalbm.py` | Command-line entry point for TOML-configured runs. |
| `src/config_loader.py` | Loads TOML config files and resolves paths. |
| `src/observation_loader.py` | Loads optional flux observations and isotope observation tables. |
| `src/analysis_runner.py` | Runs configured conditions, optimizers, plots, and output writers. |
| `src/output_writer.py` | Writes CSV/JSON output products. |
| `src/data_loader.py` | Loads ALBM NetCDF files or the processed cache. |
| `src/isotope_model.py` | Contains isotope mass balance functions and optimizers. |
| `src/flux_plots.py` | Contains plotting and paper-statistics functions. |
| `src/manuscript_scripts/multi_condition_optimization.py` | Legacy Big Trail Lake paper workflow. |
| `src/manuscript_scripts/multi_start_optimization.py` | Legacy multi-start optimization check. |

---

## 2. Data Loading

### 2.1 ALBM NetCDF Inputs

`data_loader.load_albm_data()` loads precomputed ALBM NetCDF files from a
results directory. The default Big Trail Lake date range is
`20210101_20250101`.

| NetCDF variable | Processed field | Description |
| --- | --- | --- |
| `ch4eb` | `eflux_data` | Surface ebullition flux. |
| `ch4df` | `dflux_data` | Surface diffusive flux. |
| `ch4conc` | `ch4conc_data` | Depth-integrated dissolved CH4 concentration. |
| `ch4oxid` | `och4_data` | Depth-integrated CH4 oxidation. |
| `sedch4df` | `sedch4df_data` | Mean sediment diffusive CH4 input. |
| `sedch4eb` | `sedch4eb_data` | Mean sediment ebullitive CH4 input. |
| `och4prod` | `och4prod_data` | Depth-integrated oxic CH4 production. |
| `ch4exc` | `ch4exc_data` | Depth-integrated bubble/dissolved exchange. |
| `ch4upb` | `ch4upb_data` | Upward bubbling, with values set to zero when `ch4eb` is nonzero. |
| `exbch4` | `exbch4_data` | Ice bubble exchange, normalized by lake area. |
| `icebleak` | `icebleak_data` | Ice leakage. |
| `sedtemp` | `sediment_temp_depth`, `sediment_temp_avg` | Sediment temperature. |
| `Az` | `Az` | Lake area by water-column depth. |

Molar CH4 variables are converted to mass units with:

\[
M_{\text{CH}_4} = 16.04 \times 10^{-3}\;\text{kg mol}^{-1}
\]

Surface total flux is:

\[
F_{\text{total}} = F_{\text{diff}} + F_{\text{eb}} + F_{\text{upb}}
\]

### 2.2 Depth Integration

Water-column variables are area-weighted over water depth using 0.1 m layer
thickness:

\[
\bar X(t) =
\frac{1}{A_{\text{lake}}}
\sum_k X(z_k,t)\,\Delta z\,A(z_k)
\]

The default lake area is 40,000 m2.

Sediment CH4 flux arrays are averaged along their sediment/COL dimension before
conversion to kg units:

\[
\bar F_{\text{sed,diff}}(t) =
\frac{1}{N}\sum_k F_{\text{sed,diff}}(t,k)
\]

and similarly for sediment ebullition.

### 2.3 Concentration Difference Diagnostic

`cch4_data` is computed from the depth-integrated concentration as:

\[
cch4_i = \frac{C_{i-1} - C_i}{86400}
\]

with `C_{-1}` represented by zero for the first timestep.

### 2.4 Sediment Temperature

For the Big Trail Lake NetCDF layout, `sedtemp` is expected after squeezing as:

```text
(Time, Z, COL)
```

The loader handles that case as follows:

1. Average across `COL` using nonnegative `sedch4df + sedch4eb` as the per-time
   weights.
2. Use an arithmetic mean over `COL` for timesteps with zero total sediment CH4
   flux.
3. Convert Kelvin to Celsius when values are above 200.
4. Store the resulting `(Time, Z)` array as `sediment_temp_depth`.
5. Store the arithmetic mean over `Z` as `sediment_temp_avg`.

For 2-D sediment temperature arrays, `_sediment_temp_weighted_mean` uses
available sediment CH4 flux columns as weights where dimensions allow it.
If sediment temperature is missing or has an unhandled shape, the loader prints
an `ERROR:` message and leaves `sediment_temp_depth` and `sediment_temp_avg`
unset. No synthetic sediment temperature is generated.

### 2.5 Processed Cache

`get_processed_data()` loads `data/processed_albm_data.pkl` when `use_cache` is
true and the cache exists. If the cache is stale or missing sediment temperature
with the expected depth shape, the loader tries to reload sediment temperature
from the NetCDF files. Analysis scripts keep `save_cache = false` by default;
use `python src/data_loader.py --refresh-cache` to intentionally rewrite the
processed cache.
Temperature-based isotope calculations require valid `sediment_temp_avg` data
and raise an error if it is unavailable.

The current Big Trail Lake cache spans 2021-01-01 through 2024-12-31 and has
`sediment_temp_depth.shape == (1461, 41)`.

### 2.6 Optional Methane Flux Observations

Methane flux observations are comparison data only. They are not used in the
isotope optimization objective.

The loader supports:

- the paper Excel workbook with daily mg CH4 m-2 day-1 values, converted to
  kg CH4 m-2 s-1;
- a legacy EddyPro summary report with umol CH4 m-2 s-1 values, converted to
  kg CH4 m-2 s-1;
- generic configured flux tables through `src/observation_loader.py`.

If flux observations are disabled or unavailable with `required = false`, the
configured run skips flux-comparison outputs and continues.

---

## 3. Isotope Notation and Constants

The isotope model uses VPDB isotope ratios:

\[
R_{\text{VPDB}} = 0.0112372
\]

Delta notation is:

\[
\delta^{13}\mathrm{C} =
\left(\frac{R_{\text{sample}}}{R_{\text{VPDB}}} - 1\right)1000
\]

The conversion functions are:

\[
R = R_{\text{VPDB}}\left(1 + \frac{\delta}{1000}\right)
\]

\[
\delta =
\left(\frac{R}{R_{\text{VPDB}}} - 1\right)1000
\]

implemented as `delta_to_ratio()` and `ratio_to_delta()`.

Default isotope constants in `src/isotope_model.py` are:

| Constant | Value |
| --- | ---: |
| `DEFAULT_ALPHA_AM` | 1.0238 |
| `DEFAULT_ALPHA_HM` | 1.0456 |
| `DEFAULT_ALPHA_MO` | 1.0151 |
| `DEFAULT_ALPHA_OP` | 1.0238 |
| `DEFAULT_ALPHA_TE` | 1.000 |
| `DEFAULT_ALPHA_TD` | 1.005 |
| `DEFAULT_ALPHA_EX` | 1.000 |
| `DEFAULT_F_AM` | 0.5 |
| `DEFAULT_C13_POM` | -25.6 per mil |
| `DEFAULT_DISSOLVED_CH4_MIN_CONC_FOR_ISOTOPE` | 1e-5 kg m-3 |

---

## 4. Standard Fractionation Model

The standard model is computed by:

```python
compute_isotope_timeseries(params, flux_data)
```

Required parameters:

| Parameter | Default | Bounds from `get_default_bounds()` |
| --- | ---: | --- |
| `alpha_am` | 1.0238 | 1.000 to 1.040 |
| `alpha_hm` | 1.0456 | 1.030 to 1.080 |
| `alpha_mo` | 1.0151 | 1.015 to 1.035 |
| `alpha_op` | 1.0238 | 1.000 to 1.080 |
| `f_am` | 0.5 | 0.0 to 1.0 |
| `C13_POM` | -25.6 | -28.0 to -22.0 |

The code computes:

\[
R_{\text{POM}} =
R_{\text{VPDB}}\left(1 + \frac{C13_{\text{POM}}}{1000}\right)
\]

\[
f_{\text{HM}} = 1 - f_{\text{AM}}
\]

\[
R_{\text{sed}} =
f_{\text{HM}}\frac{R_{\text{POM}}}{\alpha_{\text{HM}}}
+ f_{\text{AM}}\frac{R_{\text{POM}}}{\alpha_{\text{AM}}}
\]

`C13_sed_prod` is the scalar delta value corresponding to `R_sed`.

---

## 5. Temperature-Based Model

The temperature-based model is computed by:

```python
compute_isotope_timeseries_temp(params, flux_data)
```

Required parameters:

| Parameter | Default | Bounds from `get_default_bounds_temp()` |
| --- | ---: | --- |
| `m` | 0.0 | -10.0 to 10.0 |
| `b` | -25.0 | -100.0 to -20.0 |
| `alpha_mo` | 1.0151 | 1.015 to 1.035 |
| `alpha_op` | 1.0238 | 1.000 to 1.080 |

The sediment production delta time series is:

\[
C13_{\text{sed,prod}}(t) = m\,T_{\text{sed,avg}}(t) + b
\]

and the corresponding isotope ratio is:

\[
R_{\text{sed}}(t) =
R_{\text{VPDB}}
\left(1 + \frac{C13_{\text{sed,prod}}(t)}{1000}\right)
\]

`flux_data` must include `sediment_temp_avg`.

---

## 6. Coupled Dissolved and Bubble Pool Model

Both `compute_isotope_timeseries()` and `compute_isotope_timeseries_temp()` call:

```python
calculate_dissolved_and_bubble_pool_isotopes(...)
```

The function co-evolves:

- dissolved CH4 concentration, isotope ratio, delta, and isotope mass;
- bubble pool concentration, isotope ratio, delta, and isotope mass.

The timestep is 86,400 s.

### 6.1 Initial Conditions

At timestep zero:

\[
C_0 = \max(C_{\text{ALBM},0}, 10^{-10})
\]

\[
B_0 = \max(F_{\text{sed,eb},0}\Delta t, 10^{-12})
\]

Both isotope ratios are initialized to the first sediment production ratio.

### 6.2 Exchange Split

The exchange flux is split into:

\[
F_{\text{exc}}^+ = \max(F_{\text{exc}}, 0)
\]

\[
F_{\text{exc}}^- = \max(-F_{\text{exc}}, 0)
\]

`F_exc+` is bubble-to-dissolved exchange. `F_exc-` is dissolved-to-bubble
exchange.

### 6.3 Sink Limiting

For each pool, same-step sink fluxes are proportionally limited when requested
sink mass exceeds:

\[
\text{available mass} = \text{pool mass} + \text{current-step source flux}\Delta t
\]

The scale factors are returned as:

- `dissolved_sink_scale`
- `bubble_sink_scale`

### 6.4 Dissolved Pool Update

The dissolved mass update is:

\[
\frac{dC}{dt} =
F_{\text{sed,diff}}
+ F_{\text{ox,prod}}
+ F_{\text{exc}}^+
- F_{\text{exc}}^-
- F_{\text{oxid}}
- F_{\text{diff}}
\]

using the limited sink terms.

The isotope flux is:

\[
\frac{dM_C}{dt} =
F_{\text{sed,diff}}R_{\text{sed}}
+ F_{\text{ox,prod}}\frac{R_{\text{sed}}}{\alpha_{\text{OP}}}
+ F_{\text{exc}}^+R_{\text{bubble,prev}}
- F_{\text{exc}}^-R_{\text{diss,prev}}
- F_{\text{oxid}}\frac{R_{\text{diss,prev}}}{\alpha_{\text{MO}}}
- F_{\text{diff}}\frac{R_{\text{diss,prev}}}{\alpha_{\text{TD}}}
\]

The updated concentration and isotope mass are:

\[
C_{i} = C_{i-1} + \frac{dC}{dt}\Delta t
\]

\[
M_{C,i} = M_{C,i-1} + \frac{dM_C}{dt}\Delta t
\]

The helper `_consistent_pool_state()` clips delta to -100 to +10 per mil and
resets invalid or depleted states to the current sediment production ratio. The
dissolved depleted mass is zero.

### 6.5 Bubble Pool Update

The bubble source masses are:

\[
S_{\text{sed}} = F_{\text{sed,eb}}\Delta t
\]

\[
S_{\text{exc}} = F_{\text{exc}}^-\Delta t
\]

The bubble loss masses are:

\[
L_{\text{exc}} = F_{\text{exc}}^+\Delta t
\]

\[
L_{\text{eb}} = F_{\text{eb}}\Delta t
\]

\[
L_{\text{upb}} = F_{\text{upb}}\Delta t
\]

The source isotope mass is:

\[
M_{\text{source}} =
S_{\text{sed}}R_{\text{sed}}
+ S_{\text{exc}}R_{\text{diss,new}}
\]

The source isotope ratio is `M_source / (S_sed + S_exc)` when source mass is
positive, otherwise `R_sed`.

Bubble losses are split between isotope mass from the previous bubble pool and
isotope mass from same-day bubble sources when loss mass exceeds old pool mass.
The ebullition isotope-loss term uses `alpha_te`; its default value is 1.0.

The bubble pool concentration and isotope mass are then updated and passed
through `_consistent_pool_state()`. The bubble depleted mass is `1e-12`.

---

## 7. Dissolved Delta Filtering

After the coupled pool update, the returned dissolved delta is filtered by:

```python
filter_delta_by_concentration(delta_values, concentrations)
```

Values are set to `NaN` where modeled dissolved CH4 concentration is below
`DEFAULT_DISSOLVED_CH4_MIN_CONC_FOR_ISOTOPE`, currently `1e-5 kg m-3`.

The raw unfiltered series is returned as `delta_dissolved_raw`. The filtered
series is returned as `delta_dissolved` and is used by the cost function and
surface-emission calculation.

---

## 8. Surface Emission Delta

Surface emission delta is computed by:

```python
calculate_surface_emission_isotopes(...)
calculate_surface_emission_isotopes_temp(...)
```

The total surface flux is:

\[
F_{\text{total}} =
F_{\text{diff}} + F_{\text{eb}} + F_{\text{upb}}
\]

When `delta_bubble` is available:

\[
\delta_{\text{diff}} = \delta_{\text{dissolved}}
\]

\[
\delta_{\text{eb}} = \delta_{\text{bubble}}
\]

\[
\delta_{\text{upb}} = \delta_{\text{bubble}}
\]

If `delta_dissolved` is `NaN`, the diffusion component falls back to the
sediment production ratio for that timestep. If `delta_bubble` is unavailable
or `NaN`, ebullition and upward bubbling fall back to the sediment production
ratio.

The flux-weighted surface-emission delta is:

\[
\delta_{\text{emission}} =
\frac{F_{\text{diff}}}{F_{\text{total}}}\delta_{\text{diff}}
+ \frac{F_{\text{eb}}}{F_{\text{total}}}\delta_{\text{eb}}
+ \frac{F_{\text{upb}}}{F_{\text{total}}}\delta_{\text{upb}}
\]

The `alpha_td` argument remains in the function signature for compatibility;
the weighted surface-emission delta uses the dissolved-pool delta directly.

---

## 9. Return Values from Isotope Time-Series Functions

`compute_isotope_timeseries()` returns:

| Key | Description |
| --- | --- |
| `delta_dissolved` | Filtered dissolved CH4 delta. |
| `delta_dissolved_raw` | Unfiltered dissolved CH4 delta. |
| `dissolved_concentration_valid` | Boolean validity mask from the concentration filter. |
| `dissolved_concentration_filter_min` | Concentration threshold used by the filter. |
| `delta_emission` | Surface-emission delta. |
| `delta_bubble` | Bubble-pool delta. |
| `C13_sed_prod` | Scalar sediment production delta. |
| `R_sed_prod` | Scalar sediment production isotope ratio. |
| `ch4_conc_mb` | Dissolved CH4 concentration from the mass balance. |
| `iso_mass_diss` | Dissolved isotope mass diagnostic. |
| `iso_mass_bub` | Bubble isotope mass diagnostic. |
| `B_pool` | Bubble pool concentration diagnostic. |
| `dissolved_sink_scale` | Dissolved sink scaling diagnostic. |
| `bubble_sink_scale` | Bubble sink scaling diagnostic. |

`compute_isotope_timeseries_temp()` returns the same keys except `R_sed_prod`;
`C13_sed_prod` is a time series and `C13_sed_prod_mean` is also returned.

---

## 10. Default Observation Targets

`get_default_obs_targets()` returns four Big Trail Lake paper targets:

| Target key | Season | Mean | Std |
| --- | --- | ---: | ---: |
| `dissolved_summer` | summer | -39.08 | 10.80 |
| `dissolved_winter` | winter | -63.56 | 0.0738 |
| `emissions_summer` | summer | -54.6095 | 1.417341130462247 |
| `emissions_winter` | winter | -64.1072 | 4.704433342597597 |

The emission targets are computed from dissolved and ebullition isotope means:

\[
\delta_{\text{emission,obs}} =
f_{\text{diff}}\delta_{\text{dissolved,obs}}
+ (1 - f_{\text{diff}})\delta_{\text{ebullition,obs}}
\]

with propagated standard deviation:

\[
\sigma_{\text{emission}} =
\sqrt{
(f_{\text{diff}}\sigma_{\text{dissolved}})^2
+ ((1 - f_{\text{diff}})\sigma_{\text{ebullition}})^2
}
\]

Default diffusion fractions are:

| Season | Diffusion fraction |
| --- | ---: |
| summer | 0.13 |
| winter | 0.24 |

---

## 11. Seasonal Masks and Legacy Cost Functions

The legacy cost functions are:

```python
calculate_cost(...)
calculate_cost_temp(...)
```

They evaluate only the four default target keys. The seasonal masks are:

| Season | Dates |
| --- | --- |
| summer | June 21 through September 22 |
| winter | December 21 through March 20 |

The unnormalized cost is:

\[
J = \sum_k
\left(\bar\delta_{\text{model},k}
- \bar\delta_{\text{obs},k}\right)^2
\]

With `normalize_by_std=True`, each target error is:

\[
\left(
\frac{\bar\delta_{\text{model},k}
- \bar\delta_{\text{obs},k}}
{\sigma_{\text{obs},k}}
\right)^2
\]

If a model target mean is not finite, that target contributes no finite error
to the total cost.

---

## 12. Optimization

The legacy optimizers are:

```python
run_optimization(...)
run_optimization_temp(...)
```

Both use `scipy.optimize.differential_evolution`.

Default optimizer settings:

| Setting | Default |
| --- | ---: |
| `maxiter` | 200 |
| `tol` | 0.01 |
| `popsize` | 15 |
| `mutation` | (0.5, 1.0) |
| `recombination` | 0.7 |
| `polish` | True |
| `workers` | -1 |
| `seed` | None in the function signature; 42 in `configs/paper.toml` and legacy paper scripts |

`run_optimization()` also accepts `use_gpu` and `gpu_batch_size`. GPU evaluation
is optional, applies to the standard model path, and requires CuPy.

The configured runner uses `optimizer.backend = "legacy_when_possible"` in
`configs/paper.toml`. With that setting, conditions whose targets are the four
paper-compatible summer/winter targets call the legacy optimizers directly.
Other target sets use the generic optimizer in `src/analysis_runner.py`.

---

## 13. Configured Targets and Generic Optimization

The canonical runner supports summary isotope targets and raw observation table
targets.

Summary target example:

```toml
[[runs.isotope_targets]]
name = "dissolved_summer"
source = "summary"
component = "dissolved"
season = "summer"
mean = -39.08
std = 10.80
```

Raw-table target example:

```toml
[[runs.isotope_observation_tables]]
id = "field_isotopes"
path = "data/isotope_observations.csv"

[[runs.isotope_targets]]
name = "field_dissolved_summer"
source = "table"
table = "field_isotopes"
component = "dissolved"
season = "summer"
```

For raw tables, `src/observation_loader.py` can use standard columns:

```text
date,component,delta13c_ch4,season,depth_m,site,notes
```

or configured column mappings.

Target masks are prepared from either:

- `season = "summer"` or `season = "winter"`;
- explicit `start_date` and `end_date`;
- the full model period if neither season nor dates are supplied.

Targets may also define:

| Field | Use |
| --- | --- |
| `enabled` | Disabled targets are not selected when a condition omits `target_names`. |
| `weight` | In the generic optimizer, multiplies the target error after optional standard-deviation normalization. |
| `std` | Used for standard-deviation normalization when `normalize_by_std = true`. |

Generic target components map to isotope time-series keys as follows:

| Component | Time-series key |
| --- | --- |
| `dissolved` | `delta_dissolved` |
| `emission` | `delta_emission` |
| `emissions` | `delta_emission` |
| `emitted` | `delta_emission` |
| `bubble` | `delta_bubble` |
| `sediment_production` | `C13_sed_prod` |

`sediment_production` is a time series for the temperature model. In the static
model, `C13_sed_prod` is scalar, so it is not written as a per-timestep CSV
column and is not a time-window target for generic static optimization.

The generic optimizer uses the same parameter bounds as the legacy optimizers
unless bounds are overridden in the config.

---

## 14. Big Trail Lake Paper Conditions

`configs/paper.toml` defines five conditions:

| Condition | Model | Optimized | Targets |
| --- | --- | --- | --- |
| `Initial guess` | static | no | all four paper targets |
| `All Observations` | static | yes | all four paper targets |
| `No Winter Obs.` | static | yes | `dissolved_summer`, `emissions_summer` |
| `No Winter Emission Obs.` | static | yes | `dissolved_summer`, `dissolved_winter`, `emissions_summer` |
| `Temperature-Based` | temperature | yes | all four paper targets |

These conditions use `seed = 42`, `maxiter = 200`, `popsize = 15`, `tol = 0.01`,
and `normalize_by_std = false`.

---

## 15. Plotting

The configured runner calls plotting functions from `src/flux_plots.py` according
to the `[runs.plots]` flags in the TOML config.

Available plot flags include:

| Config flag | Plot function |
| --- | --- |
| `water_column_heatmaps` | `plot_water_column_heatmaps` |
| `water_column_doy_climatology_heatmaps` | `plot_water_column_doy_climatology_heatmaps` |
| `oxic_and_exchange_heatmaps` | `plot_oxic_and_exchange_heatmaps` |
| `eddy_comparison` | `plot_eddy_comparison` |
| `flux_components` | `plot_flux_components` |
| `flux_components_linear` | `plot_flux_components_linear` |
| `ch4_budget` | `plot_ch4_budget` |
| `pool_sources_sinks` | `plot_pool_sources_sinks` |
| `sediment_temperature` | `plot_sediment_temperature` |
| `optimization_comparison` | `plot_optimization_comparison` |
| `dissolved_pool_timeseries` | `plot_dissolved_pool_timeseries` |
| `overlap_violins` | `plot_overlap_violins` |

The paper-style isotope comparison and overlap violin plots require the four
paper target keys.

---

## 16. Output Files

The canonical runner writes one output directory per configured run:

```text
outputs/<run_id>/
  figures/
  resolved_config.toml
  resolved_config.json
  run_metadata.json
  isotope_observation_targets.csv
  optimization_parameters.csv
  isotope_timeseries.csv
  flux_timeseries.csv
  methane_flux_comparison.csv
  observation_overlap_statistics.csv
```

`methane_flux_comparison.csv` is written only when methane flux observations are
loaded. Manuscript paper-statistics output is written by
`src/manuscript_scripts/multi_condition_optimization.py`.

The legacy `src/manuscript_scripts/multi_condition_optimization.py` script writes to `outputs/`.

---

## 17. Paper Statistics

`flux_plots.compute_and_print_paper_statistics()` writes:

```text
paper_statistics_flux_isotope.csv
```

The function expects the four paper target keys and uses named paper conditions
when present:

- `Initial guess`
- `All Observations`
- `Temperature-Based`

If a condition is absent, its section is skipped or reported as not available.
