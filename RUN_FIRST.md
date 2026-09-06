# Run first

The GitHub package contains code and environment files only. First obtain the companion `NeuroTRACE_Reproducibility_Archive_FINAL.zip` from Zenodo (https://doi.org/10.5281/zenodo.22525118) and extract its directories into this repository.

From the repository root, prepare the data directories:

```bash
mkdir -p data/processed data/figure_source_data data/intermediate_results
ARCHIVE_ROOT=../NeuroTRACE_Reproducibility_Archive
cp -R "$ARCHIVE_ROOT/processed_data/." data/processed/
cp -R "$ARCHIVE_ROOT/figure_source_data/." data/figure_source_data/
cp -R "$ARCHIVE_ROOT/intermediate_results/." data/intermediate_results/
```

Install the environment:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r environment/requirements.txt
```

Run the final workflow in this order:

```bash
python scripts/analysis/run_dual_head_ppr.py
python scripts/simulation/run_matched_random_calibration.py
python scripts/simulation/run_confirmatory_simulation.py
python scripts/simulation/summarize_simulation_tradeoffs.py
python scripts/validation/run_external_asd_validation.py
python scripts/figures/render_final_figures.py
```

The first five commands require the prepared upstream inputs described in `data/README.md`. The figure renderer consumes the extracted source tables and intermediate plot inputs. These commands are a reproducibility workflow; they do not alter the manuscript or download biological data.
