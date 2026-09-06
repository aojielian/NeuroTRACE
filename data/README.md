# Data availability and archive extraction

This GitHub package contains no processed data or figure source data. Obtain `NeuroTRACE_Reproducibility_Archive_FINAL.zip` from Zenodo (https://doi.org/10.5281/zenodo.22525118) and extract the following directories into the repository root:

* `processed_data/` → `data/processed/`
* `figure_source_data/` → `data/figure_source_data/`
* `intermediate_results/` → `data/intermediate_results/`

The archive contains derived, redistributable tables only. It does not contain raw expression matrices, subject-level records, or controlled-access resources.

Public resources used by the final workflow include:

* GEO GSE102741;
* GEO GSE64018;
* BrainSpan developmental reference profiles;
* PsychENCODE resources;
* the Velmeshev single-cell reference.

Obtain raw or controlled-access inputs from their original providers using these identifiers and prepare the upstream graph, native module weights, developmental profiles, cohort matrices, metadata, matched-random draws, and simulation helper inputs expected by the scripts. The extracted `figure_source_data/FIGURE_SOURCE_DATA_INDEX.tsv` maps each final figure panel to its source table.

The `processed_data/README_PROCESSED_DATA.txt` file documents the derived tables, source resources, processing role, consuming scripts, and redistribution scope. The archive metadata directory contains inventories and the reproducibility manifest.
