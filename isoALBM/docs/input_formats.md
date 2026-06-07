# Input Formats

## ALBM Inputs

The runner accepts either:

- A processed pickle cache, normally `data/processed_albm_data.pkl`.
- A directory of precomputed ALBM NetCDF outputs named like
  `isobLakeOut.<variable>.20210101_20250101.nc`.

## Methane Flux Observations

Methane flux observations are optional and used only for comparison plots and
statistics. They are not used in the isotope optimization objective.

The paper config uses:

```toml
[runs.flux_observations]
enabled = true
required = false
format = "auto"
path = "data/2026_05_08_EC_CH4.xlsx"
```

If the file is absent and `required = false`, the runner continues without
`methane_flux_comparison.csv` or the flux-comparison figure.

For generic CSV/XLSX flux tables:

```toml
[runs.flux_observations]
enabled = true
format = "generic_table"
path = "data/flux_observations.csv"
date_column = "date"
flux_column = "ch4_flux"
units = "mg_m2_day"
```

Supported units:

- `kg_m2_s`
- `mg_m2_day`
- `umol_m2_s`

## Isotope Observation Summary Targets

Summary targets provide mean and standard deviation directly:

```toml
[[runs.isotope_targets]]
name = "dissolved_summer"
source = "summary"
component = "dissolved"
season = "summer"
mean = -39.08
std = 10.80
```

The model period can be a season:

- `summer`: June 21 through September 22.
- `winter`: December 21 through March 20.

Or an explicit date window:

```toml
start_date = "2022-07-01"
end_date = "2022-08-15"
```

## Raw Isotope Observation Tables

The recommended standard table columns are:

```csv
date,component,delta13c_ch4,season,depth_m,site,notes
2022-07-12,dissolved,-41.2,summer,1.5,Big Trail Lake,
2022-07-12,emission,-55.1,summer,,Big Trail Lake,
```

With standard columns, the table config can be minimal:

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

For real-world files with different headings, map the columns:

```toml
[[runs.isotope_observation_tables]]
id = "field_isotopes"
path = "data/isotope_observations.xlsx"
sheet = "Sheet1"
date_column = "Sample Date"
value_column = "d13C_CH4"
component_column = "Sample Type"
season_column = "Season"
depth_column = "Depth_m"
site_column = "Lake"
```

The target loader filters by `component`, `season`, `start_date`, `end_date`,
and optional `filter_site`, then computes the target mean and standard
deviation from the selected rows.

## Component Names

Common aliases are normalized:

- `dissolved`, `dissolved_ch4`, `water` -> dissolved model isotope.
- `emission`, `emissions`, `emitted`, `surface_emission` -> surface emitted isotope.
- `bubble`, `ebullition` -> modeled bubble pool isotope.

The generic optimizer can use any configured target set. The original paper
optimizer is used only when selected targets match the four legacy-compatible
summer/winter targets.
