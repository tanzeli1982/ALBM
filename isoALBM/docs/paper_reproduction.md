# Paper Reproduction

The Big Trail Lake manuscript case is preserved in `configs/paper.toml` and in
the legacy scripts that were used before the configurable runner was added.

## Canonical Command

```bash
python src/run_isoalbm.py --config configs/paper.toml
```

This writes structured outputs to:

```text
outputs/big_trail_lake_paper/
```

The paper config uses the same:

- processed ALBM cache: `data/processed_albm_data.pkl`
- ALBM NetCDF directory: `ALBM data/`
- date range: `20210101_20250101`
- isotope targets
- optimization conditions
- optimizer seed and settings

The optimizer backend is `legacy_when_possible`, so the four paper-compatible
summer/winter target conditions use the same optimizer functions as
`src/manuscript_scripts/multi_condition_optimization.py`. This is intended to keep manuscript results
within numerical tolerance after the refactor.

## Legacy Command

```bash
python src/manuscript_scripts/multi_condition_optimization.py
```

This script is kept as the paper-specific workflow. It writes to
`figs/` and remains useful for direct comparison with pre-refactor behavior.

## Regression Smoke Test

A fast smoke test for the configurable runner is:

```bash
python src/run_isoalbm.py --config configs/paper.toml --condition Temperature-Based --maxiter 1 --workers 1 --no-plots --output-dir outputs/smoke_config_runner
```

With `maxiter = 1`, optimizer `success` is expected to be false. The regression
signal is that the canonical runner uses the legacy backend and produces the
same finite cost and parameter values as a direct call to
`run_optimization_temp` with the same settings.

## Full Results

For configurable model-output reproduction, run the canonical command without
smoke-test overrides. The most important CSV outputs are:

- `optimization_parameters.csv`
- `isotope_timeseries.csv`
- `flux_timeseries.csv`
- `observation_overlap_statistics.csv`

For manuscript numeric statistics, run:

```bash
python src/manuscript_scripts/multi_condition_optimization.py
```

That script writes:

- `outputs/paper_statistics_flux_isotope.csv`

If methane flux observations are available which is currently not true due to embargoed data, the run also writes:

- `methane_flux_comparison.csv`
- the eddy comparison figure
