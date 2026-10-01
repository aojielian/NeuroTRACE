# NeuroTRACE reproducibility release

Developmental mapping resolves opposing expression programs in autism cortex.

This repository contains the code, processed reference inputs, final source tables and methods corresponding to the submission snapshot dated 2026-10-01. No raw expression records, individual-level scores, or controlled-access records are redistributed.

## Developmental reference

The six stages, in order, are Early prenatal, Mid prenatal, Late prenatal, Early postnatal, Childhood–adolescence and Adulthood. These are descriptive labels for the processed NeuroTRACE reference. Their observed ranges are 8–12, 13–21, 24–37, 56–80, 92–820 and 976–2120 post-conception weeks, respectively. The nomenclature does not change membership or numeric profiles. Figure abbreviations are EPr, MPr, LPr, EPost, C/A and Adult.

## Final calibration and intervals

The final calibration uses one standardized developmental-profile mean/variance matcher without replacement for NTM1–NTM3 at Top 200 and Top 500, 1,000 null draws each, 6,000 total, seed 202609301. Propagation settings are restart probability 0.35, tolerance 1e-10 and maximum 120 iterations. `scripts/calibration/unified_matcher.R` is the final calibration script. Frozen full replicate and summary tables are in `processed_results/`; Table S6 is in `supplementary_tables/`.

External ASD intervals use Student-t residual degrees of freedom 18 (GSE64018) and 47 (GSE102741), with coefficients, SE, P and FDR retained. All six GSE53987 MDD individual models use df_resid=30 and Student-t 95% confidence intervals from saved beta and SE; saved t/P independently verify df. Table S10 includes all six models. S7 cross-disorder means and direction counts are descriptive, with no pooled inference.

## Structure and use

- `processed_inputs/`: aggregate six-stage BrainSpan profiles, native module weights and frozen graph.
- `processed_results/`: complete null replicates, calibration, external intervals and six MDD models.
- `figure_source_data/`: final Fig1–Fig5 and S1–S7 source tables.
- `supplementary_tables/`: supplementary table mirrors.
- `scripts/figures/`: final deterministic plotting code, including the sequential Figure5/S2 quantitative palette.
- `documentation/`: supplementary methods, processed-input guidance and source index.
- `metadata/`: release file and submission snapshot hashes.

See `RUN_FIRST.md` for environment and commands. Historical analysis/simulation comparator scripts require unbundled upstream inputs and helpers. The runnable final entry points are the frozen-table renderer, saved MDD CI verification and unified matcher documented in RUN_FIRST.md. No analysis was rerun to prepare this snapshot.

## Archive availability

The reproducibility archive for this release is available at
https://doi.org/10.5281/zenodo.23083241.

The concept DOI for the NeuroTRACE archive is
https://doi.org/10.5281/zenodo.20159203.

This version archives the final processed reference inputs, saved calibration replicates, source tables, plotting scripts and Supplementary Methods.

## Citation and license

Please cite the associated manuscript, this repository and the original data providers. Software citation metadata are in `CITATION.cff`. Code remains under the repository MIT license; original resource reuse conditions apply.
