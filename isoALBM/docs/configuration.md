# Configuration

To run isoALBM:

```bash
python src/run_isoalbm.py --config configs/paper.toml
```

Configs are TOML files. Relative paths are resolved from `base_dir`, which is
itself relative to the config file. In `configs/paper.toml`, `base_dir = ".."`
means paths such as `data/processed_albm_data.pkl` resolve from the repository
root.

## Top Level

```toml
version = 1
base_dir = ".."

[[runs]]
id = "example_site"
output_dir = "outputs/example_site"
```

Use multiple `[[runs]]` blocks to run multiple sites or ALBM experiments. Each
run is independent and writes to its own `output_dir`.

## Data

```toml
[runs.data]
cache_path = "data/processed_albm_data.pkl"
results_dir = "ALBM data"
date_range = "20210101_20250101"
use_cache = true
save_cache = false
lake_area = 40000.0
```

If `use_cache = true` and `cache_path` exists, the runner loads the processed
ALBM cache. If not, it loads the NetCDF files in `results_dir` for `date_range`.

## Optional Methane Flux Observations

```toml
[runs.flux_observations]
enabled = true
required = false
format = "auto"
path = "data/2026_05_08_EC_CH4.xlsx"
```

When `enabled = false`, or when `required = false` and the file is missing, the
runner skips methane flux comparison outputs and continues with ALBM/isotope
analysis. This is useful for public releases where flux observations are
embargoed.

For generic flux tables:

```toml
[runs.flux_observations]
enabled = true
format = "generic_table"
path = "data/flux_observations.csv"
date_column = "date"
flux_column = "ch4_flux"
units = "kg_m2_s"
```

Supported generic units are `kg_m2_s`, `mg_m2_day`, and `umol_m2_s`.

## Optimizer

```toml
[runs.optimizer]
backend = "legacy_when_possible"
maxiter = 200
tol = 0.01
popsize = 15
workers = -1
seed = 42
normalize_by_std = false
polish = true
verbose = false
```

`legacy_when_possible` uses the original paper optimizer whenever the selected
targets are the four paper-compatible summer/winter targets. Other target sets
use the generic optimizer, which can optimize against arbitrary configured
target windows.

## Isotope Targets

Summary targets:

```toml
[[runs.isotope_targets]]
name = "dissolved_summer"
source = "summary"
component = "dissolved"
season = "summer"
mean = -39.08
std = 10.80
```

Raw-table targets:

```toml
[[runs.isotope_observation_tables]]
id = "field_isotopes"
path = "data/isotope_observations.csv"

[[runs.isotope_targets]]
name = "custom_dissolved_summer"
source = "table"
table = "field_isotopes"
component = "dissolved"
season = "summer"
```

Targets can use `season = "summer"` or `season = "winter"`, or explicit
`start_date` and `end_date` fields. Components currently mapped to model output
are `dissolved`, `emission`, `bubble`, and `sediment_production`.

## Conditions

Each condition chooses a model and a target subset:

```toml
[[runs.conditions]]
name = "All Observations"
model = "static"
optimize = true
target_names = ["dissolved_summer", "dissolved_winter", "emissions_summer", "emissions_winter"]
color = "#E69F00"
linestyle = "-"
```

Use `model = "temperature"` for the temperature-dependent sediment isotope
production model. Use `optimize = false` for an initial/default parameter run.

## Plots

Plot flags correspond to functions in `src/flux_plots.py`:

```toml
[runs.plots]
enabled = true
show = false
eddy_comparison = true
flux_components = true
flux_components_linear = true
ch4_budget = true
pool_sources_sinks = true
sediment_temperature = true
optimization_comparison = true
dissolved_pool_timeseries = true
overlap_violins = true
```

The paper-style isotope comparison and overlap violin plots require the four
paper target keys. Other configured analyses still write CSV outputs even when
those plots are skipped.

## Command Overrides

Useful smoke-test overrides:

```bash
python src/run_isoalbm.py --config configs/paper.toml --condition Temperature-Based --maxiter 1 --workers 1 --no-plots --output-dir outputs/smoke_config_runner
```
