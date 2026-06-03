# Output Formats

The runner writes one folder per configured run:

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

`methane_flux_comparison.csv` is written only when methane flux observations
are enabled and loaded. Manuscript paper-statistics output is written by
`src/manuscript_scripts/multi_condition_optimization.py`.

## `run_metadata.json`

Machine-readable run metadata:

- run id
- config and data paths
- ALBM date range
- number of timesteps
- sediment temperature array shape
- whether methane flux observations were loaded
- condition names

## `isotope_observation_targets.csv`

One row per target after config and raw-table processing. Important columns:

- `name`
- `component`
- `season`
- `start_date`
- `end_date`
- `mean`
- `std`
- `weight`
- `source`
- `source_file`
- `n_observations`

This file records the exact targets used by the optimizer.

## `optimization_parameters.csv`

One row per condition. Columns include:

- `condition`
- `model`
- `optimized`
- `backend`
- `cost`
- `success`
- `n_evaluations`
- parameter columns such as `alpha_am`, `alpha_hm`, `alpha_mo`, `alpha_op`,
  `f_am`, `C13_POM`, `m`, and `b` depending on model type.

For very small smoke tests, `success` can be `False` because the optimizer is
stopped by `maxiter`, but the finite cost and parameters are still useful for
regression checks.

## `isotope_timeseries.csv`

Daily modeled isotope and pool time series. Column names use:

```text
<condition>.<variable>
```

Typical variables:

- `delta_dissolved`
- `delta_emission`
- `delta_bubble`
- `C13_sed_prod`
- `ch4_conc_mb`

## `flux_timeseries.csv`

Daily ALBM flux and temperature data:

- `time`
- `total_flux_kg_m2_s`
- `diffusion_flux_kg_m2_s`
- `ebullition_flux_kg_m2_s`
- `upward_bubbling_flux_kg_m2_s`
- `sediment_temperature_avg_c`

## `methane_flux_comparison.csv`

Written only when flux observations are loaded. Columns:

- `date`
- `albm_total_flux_kg_m2_s`
- `observed_flux_kg_m2_s`
- `albm_minus_observed_kg_m2_s`

## `observation_overlap_statistics.csv`

One row per condition and target:

- `condition`
- `target`
- `component`
- `model_mean`
- `model_std`
- `observed_mean`
- `observed_std`
- `raw_error`
- `weighted_error`
- target source metadata

## Figures

Figures are written under `outputs/<run_id>/figures/`. The available plot flags
are documented in `docs/configuration.md` and map to functions in
`src/flux_plots.py`.
