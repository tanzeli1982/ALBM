# isoALBM

This repository is the code and data bundle for isoALBM analysis from precomputed ALBM output developed by Kevin Rozmiarek. Please refer to the Documentation folder and `isotope_model_documentation` for a full description. The default example case is for Big Trail Lake, Alaska and corresponds to the manuscript "Controls on δ13C − CH4 from Arctic Thermokarst Lakes".

## Maintainer
Kevin Rozmiarek - kevin.rozmiarek@colorado.edu

## Contents

Configurable workflow:

- `src/run_isoalbm.py` - command-line entry point for TOML-configured isoALBM runs.
- `configs/paper.toml` - Example configuration for Big Trail Lake
- `src/config_loader.py`, `src/observation_loader.py`, `src/analysis_runner.py`,
  and `src/output_writer.py` - support modules for config loading, raw observation
  tables, configurable optimization, and CSV output.

Core model:

- `src/data_loader.py` - loads ALBM NetCDF, observation, and cached ALBM data.
- `src/isotope_model.py` - isotope mass balance, cost functions, and optimizers.
- `src/flux_plots.py` - flux, isotope, overlap, and paper-statistic utilities.
- `isotope_model_documentation.md` - detailed model equations and reference notes.

Manuscript scripts:

- `src/manuscript_scripts/multi_condition_optimization.py` - legacy manuscript-specific isotope comparison run.
- `src/manuscript_scripts/multi_start_optimization.py` - multi-start uniqueness and robustness checks.

Existing inputs for manuscript:

- `data/processed_albm_data.pkl` - paper-period cache used by the default runs.
- `data/2026_05_08_EC_CH4.xlsx` - eddy covariance methane flux observations.
- `data/2022_11_21_CH4CO2_DeepShallow.xlsx` - dissolved CH4/CO2 profile observations.
- `ALBM data/isobLakeOut.*.20210101_20250101.nc` - precomputed Big Trail Lake ALBM NetCDF outputs for transparency and cache refreshes.

## Quickstart

Create and activate a local Python environment, then install the runtime
dependencies:

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

Verify the curated cache:

```bash
python -c "import sys; sys.path.insert(0, 'src'); from data_loader import load_processed_data; a,e=load_processed_data('data/processed_albm_data.pkl'); print(a.start_date, a.end_date, a.sediment_temp_depth.shape)"
```

Run the canonical configurable manuscript preset:

```bash
python src/run_isoalbm.py --config configs/paper.toml
```

Run a very small canonical smoke test:

```bash
python src/run_isoalbm.py --config configs/paper.toml --condition Temperature-Based --maxiter 1 --workers 1 --no-plots --output-dir outputs/smoke_config_runner
```

Run the legacy paper script, which is kept as a direct record of the manuscript
workflow before the configurable runner:

```bash
python src/manuscript_scripts/multi_condition_optimization.py
```

The full multi-condition optimization can take substantially longer than the
smoke test. The canonical runner writes structured outputs under
`outputs/<run_id>/`, including CSV analysis tables and `figures/`.

## Data

The paper workflow reads `data/processed_albm_data.pkl`, a compact cache
containing the paper-period ALBM fields for 2021-01-01 through 2024-12-31.
Optional methane flux observations are loaded separately for comparison plots
and statistics when configured and present. If the flux file is absent or
disabled, isotope optimization and ALBM-only outputs still run.

The cache can be regenerated from the included NetCDF set and configured flux
workbook with:

```bash
python src/data_loader.py --refresh-cache
```

The included NetCDF files are precomputed ALBM output for Big Trail Lake over
the date range encoded in the filenames, `20210101_20250101`. This repository
does not include the upstream climate forcing preparation or GIS/DSM layers

## Manuscript Figure and Statistic Mapping

| Manuscript item | Canonical command | Primary outputs |
| --- | --- | --- |
| ALBM flux diagnostics and optional model-data flux comparison | `python src/run_isoalbm.py --config configs/paper.toml` | `outputs/big_trail_lake_paper/figures/ALBM_*.png`, `outputs/big_trail_lake_paper/flux_timeseries.csv`, optional `methane_flux_comparison.csv` |
| Water-column and sediment-temperature diagnostics | enable the relevant plot flags in `configs/paper.toml` | `outputs/big_trail_lake_paper/figures/isoALBM_*heatmaps.png`, `isoALBM_sediment_temperature_*.png` |
| Main isotope condition comparison | `python src/run_isoalbm.py --config configs/paper.toml` | `outputs/big_trail_lake_paper/figures/isoALBM_optimization_comparison.png`, `isotope_timeseries.csv`, `optimization_parameters.csv` |
| Observation-overlap isotope distributions | `python src/run_isoalbm.py --config configs/paper.toml` | `outputs/big_trail_lake_paper/figures/isoALBM_overlap_violins.png`, `observation_overlap_statistics.csv` |
| Manuscript numeric statistics, sections B-J | `python src/manuscript_scripts/multi_condition_optimization.py` | `outputs/paper_statistics_flux_isotope.csv` |
| Multi-start optimizer robustness check | `python src/manuscript_scripts/multi_start_optimization.py --n-runs 100 --model both` | `figs/multi_start/*.csv`, `*.json`, `multi_start_*_params.png`, `multi_start_costs.png`, `multi_start_summary_statistics.csv` |

## Documentation

- `docs/configuration.md` - TOML schema and multi-run examples.
- `docs/input_formats.md` - ALBM, methane flux, and isotope observation inputs.
- `docs/output_formats.md` - CSV and figure outputs from the canonical runner.
- `docs/paper_reproduction.md` - reproduction of the work in the manuscript ""Controls on δ13C − CH4 from Arctic Thermokarst Lakes""

## Requirements

The runtime dependencies are listed in `requirements.txt`:

- `numpy`
- `pandas`
- `matplotlib`
- `scipy`
- `tqdm`
- `netCDF4`
- `openpyxl`

## License and Citation

Software source code is licensed under the BSD 3-Clause License; see `LICENSE`.
Unless otherwise noted, curated data files and generated CSV outputs are
licensed under Creative Commons Attribution 4.0 International (CC BY 4.0); see
`DATA_LICENSE.md`.

See `CITATION.cff` for a placeholder citation record. Will be updated with the final
manuscript DOI, author list, and release metadata when available.
