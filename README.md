# NeuroTRACE public reproducibility release

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22525118.svg)](https://doi.org/10.5281/zenodo.22525118)

This repository contains the NeuroTRACE analysis code, execution environment, and usage documentation. Reproducibility data and figure source tables are distributed in the companion archive available from Zenodo (https://doi.org/10.5281/zenodo.22525118); publication-formatted figures are not included in this repository.

## Final method

The final NeuroTRACE method is:

signed native modules  
→ positive/negative restart channels  
→ gene-only PPR  
→ `q_pos` / `q_neg`  
→ `q_magnitude` for developmental localization  
→ `q_signed` for direction-aware gene scoring

Positive and negative native weights are propagated independently on a gene-only graph. `q_magnitude = q_pos + q_neg` is used for developmental localization against six processed developmental profiles. `q_signed = q_pos - q_neg` is used for direction-aware gene scoring; gene priority is `abs(q_signed)` and direction is `sign(q_signed)`.

The final graph uses gene-gene embedding k-nearest-neighbour edges, `alpha = 0.35`, tolerance `1e-10`, and at most 120 iterations. No heterogeneous stage-node or risk-node graph is used as the final method.

## Repository structure

* `scripts/analysis`: dual-head gene-only PPR analysis.
* `scripts/simulation`: matched-random calibration and confirmatory simulation.
* `scripts/validation`: external ASD validation.
* `scripts/figures`: deterministic rendering from the extracted source tables.
* `environment`: Python and conda environment specifications.
* `data/README.md`: data acquisition and archive-extraction instructions.
* `metadata`: public release metadata and code/file inventories.

This GitHub package intentionally contains no rendered figures, supplementary files, figure source data, processed data, or intermediate result tables.

## Data and accessions

The reproducibility archive contains derived, redistributable tables associated with public resources including GEO GSE102741 and GSE64018, BrainSpan developmental profiles, PsychENCODE resources, and the Velmeshev single-cell reference. Raw expression matrices and controlled-access records are not redistributed. See `data/README.md`.

## Environment setup

Python 3.10 or newer is supported:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -r environment/requirements.txt
```

Alternatively, create the conda environment from `environment/environment.yml`.

## Reproducibility workflow

1. Download the companion NeuroTRACE reproducibility archive from Zenodo: https://doi.org/10.5281/zenodo.22525118
2. Extract `processed_data/`, `figure_source_data/`, and `intermediate_results/` into the repository as described in `RUN_FIRST.md`.
3. Run the analysis and validation scripts in the declared order.
4. Run the figure renderer after the source tables have been extracted.

The scripts use repository-relative paths and do not download biological data automatically.

## Citation and license

Please cite the associated NeuroTRACE manuscript, this repository, the reproducibility archive, and the original data providers. `CITATION.cff` contains the software citation metadata without a DOI placeholder. Code is distributed under the MIT License; source datasets remain subject to their original providers' terms.
