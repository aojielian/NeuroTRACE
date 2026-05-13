# NeuroTRACE

This repository contains the **final analysis scripts**, **main-figure scripts**, and **supplementary-figure scripts** used for the manuscript:

**NeuroTRACE links adult cortical disease modules to developmental timing, genetic convergence, and cross-disorder psychiatric axes**

NeuroTRACE is a graph-informed disease-to-developmental mapping framework that derives signed native transcriptomic modules from adult cortical disease expression data, embeds them in a developmental reference space, and applies personalized PageRank-based heterogeneous graph diffusion to estimate developmental stage affinity and graph-prioritized genes.

## Repository scope

This repository provides the code required to inspect and reproduce the final NeuroTRACE analysis workflow. It includes:

- final analysis scripts used for the manuscript;
- final scripts for main Figures 2–6;
- documentation for the manually prepared Figure 1 schematic;
- final scripts for Supplementary Figures S1–S5;
- code and figure manifests;
- environment information and basic reproducibility notes.

Large raw public datasets and processed analysis outputs are stored outside this code repository. Processed analysis outputs, figure source data, supplementary-table source data, and data manifests are available from Zenodo:

**https://doi.org/10.5281/zenodo.20159204**

Raw public datasets should be obtained from their original sources as listed in Supplementary Table S1 of the manuscript and in the data manifest provided with the Zenodo archive.

## Repository layout

```text
analysis/
  step02_simulation_benchmark/
  step03_feature_embedding/
  step05_optimal_transport_alignment/
  step06_cross_cohort_validation/
  step08_algorithm_strengthening/
  step10_risk_free_ablation/
  step10_transport_null/
  step10_zhou2022_genesets_fixed/
  step10_independent_genetic_risk/
  step10_composition_adjustment/
  step10_single_cell_validation/
  step11_robustness_sensitivity/
  step12_revision_strengthening/

figures/
  main_figures/
    Figure1/
    Figure2/
    Figure3/
    Figure4/
    Figure5/
    Figure6/
  supplementary_figures/
    Supplementary_Figure_S1/
    Supplementary_Figure_S2/
    Supplementary_Figure_S3/
    Supplementary_Figure_S4/
    Supplementary_Figure_S5/

docs/
  code_manifest.tsv
  figure_script_manifest.tsv
  input_data_manifest_template.tsv
  zenodo_file_manifest_template.tsv

tables/
  supplementary_table_notes.md
```

## Included analysis code

### 1. Simulation benchmark
Located in `analysis/step02_simulation_benchmark/`

These scripts implement the simulation benchmark used to evaluate NeuroTRACE and baseline methods under weak-signal, low-overlap, composition-confounded, high-dropout, false-prior-stress, and combined-hard settings. They also generate ablation summaries and final simulation freeze outputs.

### 2. Developmental embedding and native module derivation
Located in `analysis/step03_feature_embedding/`

These scripts construct the BrainSpan developmental embedding, derive native Gandal adult cortical disease modules, perform gene-symbol harmonization, calculate NTM-to-BrainSpan stage alignment, and generate native module summaries.

### 3. Graph-informed stage affinity and graph priority
Located in `analysis/step05_optimal_transport_alignment/`

These scripts construct graph-informed module-to-stage affinity outputs, personalized PageRank-based transport probability tables, top transport stage summaries, graph-prioritized gene lists, and graph-priority enrichment results.

### 4. External ASD validation
Located in `analysis/step06_cross_cohort_validation/`

These scripts identify external adult ASD cohort inputs, project NTM scores into GSE102741 and GSE64018, fit external validation models, and perform matched-random specificity analyses.

### 5. Algorithm strengthening and scoring baselines
Located in `analysis/step08_algorithm_strengthening/`

These scripts perform expression/variance-matched random specificity tests, signed versus unsigned scoring comparisons, cross-disease resource preparation, microarray probe-to-symbol collapsing, and cross-disease NTM scoring.

### 6. Risk-prior ablation
Located in `analysis/step10_risk_free_ablation/`

This script reruns graph diffusion after removing genetic-risk prior edges and evaluates no-risk graph outputs, including NTM stage affinity and downstream genetic enrichment.

### 7. Matched-random developmental transport nulls
Located in `analysis/step10_transport_null/`

This script generates matched-random developmental transport nulls used to calibrate observed stage-affinity behavior.

### 8. Zhou 2022 gene-set enrichment
Located in `analysis/step10_zhou2022_genesets_fixed/`

These scripts extract Zhou 2022 gene sets and test graph-prioritized and native module genes against Zhou 2022 genetic-resource categories.

### 9. Independent Satterstrom ASD exome-risk enrichment
Located in `analysis/step10_independent_genetic_risk/`

This script evaluates enrichment of NTM-associated gene sets against the Satterstrom ASD exome-risk gene set.

### 10. Bulk composition adjustment
Located in `analysis/step10_composition_adjustment/`

This script calculates bulk NTM scores, marker-based cell-type proxy scores, marker principal components, and diagnosis models before and after marker-based composition adjustment.

### 11. Single-cell and single-nucleus contextualization
Located in `analysis/step10_single_cell_validation/`

These scripts audit single-cell/single-nucleus objects and calculate broad cell-type NTM localization summaries in PsychENCODE and Velmeshev resources.

### 12. Robustness and sensitivity analyses
Located in `analysis/step11_robustness_sensitivity/`

These scripts implement hyperparameter sensitivity, graph versus baseline comparisons, gene-perturbation/hub sensitivity, NTM2 transport nulls, cross-disease multiple-testing summaries, alternative-weight sensitivity, BrainSpan-universe sensitivity, NNLS deconvolution, and robustness freeze summaries.

### 13. Revision-stage strengthening analyses
Located in `analysis/step12_revision_strengthening/`

These scripts implement discovery-module bootstrap stability, degree-matched graph nulls, stage-alignment weighting sensitivity, and held-out gnomAD constraint validation.

## Figure-generation code

### Main figures

Final main-figure scripts are located in `figures/main_figures/`.

- `Figure1/`: conceptual NeuroTRACE schematic documentation. Figure 1 is manually prepared and does not have quantitative plot data.
- `Figure2/figure2_simulation_signed_scoring_v3.R`: simulation benchmark and signed-scoring validation.
- `Figure3/figure3_native_modules_developmental_transport_FINAL_v1.R`: native module stage affinity and developmental transport.
- `Figure4/figure4_genetic_convergence_graph_prioritization_v3.R`: genetic convergence and graph-prioritized gene validation.
- `Figure5/figure5_external_validation_composition_v4.R`: external ASD validation and composition sensitivity.
- `Figure6/figure6_cross_disease_robustness_v4.R`: cross-disease generalization and MDD boundary.

### Supplementary figures

Final supplementary-figure scripts are located in `figures/supplementary_figures/`.

- `Supplementary_Figure_S1/supplementary_figure_S1_extended_simulation_v3.R`: extended simulation and scoring-baseline analyses.
- `Supplementary_Figure_S2/supplementary_figure_S2_developmental_transport_v8.R`: developmental transport calibration and stage-alignment sensitivity.
- `Supplementary_Figure_S3/supplementary_figure_S3_genetic_functional_convergence_v5.R`: genetic and functional convergence analyses.
- `Supplementary_Figure_S4/supplementary_figure_S4_composition_singlecell_v3.R`: composition and single-nucleus contextualization analyses.
- `Supplementary_Figure_S5/supplementary_figure_S5_robustness_boundary_v3.R`: robustness and cross-disease boundary analyses.

The corresponding figure source data are available in the Zenodo data/results release:

**https://doi.org/10.5281/zenodo.20159204**

## Processed data and Zenodo release

This code repository is designed to be used together with the Zenodo data/results archive. The Zenodo archive contains:

- figure source data for main Figures 2–6 and Supplementary Figures S1–S5;
- processed outputs from NeuroTRACE analysis steps;
- supplementary-table source/evidence files;
- file-level and input-data manifests;
- environment and session information.

In general:

1. Use this GitHub repository for code and workflow structure.
2. Use the Zenodo archive for processed data inputs and figure source tables.
3. Use `docs/figure_script_manifest.tsv` to identify the final script for each figure.
4. Use the Zenodo `figure_plotdata/` directory as the direct data source for final figure reproduction.

## Expected usage

This repository is intended for readers and reviewers who want to:

- inspect the final analytical workflow;
- understand how NeuroTRACE modules were derived;
- review the graph-informed developmental mapping procedure;
- reproduce the final figures from processed plot data;
- trace simulation, validation, enrichment, composition, cross-disease, and robustness analyses;
- adapt NeuroTRACE scripts to related disease-to-developmental mapping problems.

A typical reproduction workflow is:

```text
1. Clone this GitHub repository.
2. Download and unpack the Zenodo processed data/results archive.
3. Install the R/Python environment described in environment.yml and sessionInfo.txt.
4. Update project-root paths in the relevant scripts if needed.
5. Run final figure scripts under figures/main_figures/ and figures/supplementary_figures/ using the Zenodo figure_plotdata inputs.
6. Inspect processed_results/ in the Zenodo archive for intermediate analysis outputs.
```

## Environment

The project was developed in an HPC R/Python environment. The repository includes:

- `environment.yml`: conda-style environment specification;
- `sessionInfo.txt`: R session information from the analysis environment.

Some scripts may require path adjustment outside the original project directory. Slurm submission scripts are not included because they are cluster-specific and not required to inspect or reproduce the analytical workflow from scripts and processed data.

## Manifests

The `docs/` directory contains:

- `code_manifest.tsv`: list of files included in this code release;
- `figure_script_manifest.tsv`: final script associated with each main and supplementary figure;
- `input_data_manifest_template.tsv`: public input-data resource template;
- `zenodo_file_manifest_template.tsv`: expected Zenodo archive structure.

The Zenodo archive provides the processed data file manifest and MD5 checksums for archived data/result files.

## Repository philosophy

This repository prioritizes **final manuscript-matching code** rather than complete script-history preservation. Development versions, failed intermediate figure attempts, cluster-specific Slurm files, logs, temporary files, and large raw data objects are not part of the public code release. The final scripts included here correspond to the manuscript figures, supplementary figures, and final analysis outputs described in the manuscript and supplementary tables.

## Citation

If you use this code or processed results, please cite the associated manuscript and the Zenodo data/results archive:

**Processed data/results DOI:** https://doi.org/10.5281/zenodo.20159204

A separate Zenodo code DOI can be added after archiving the GitHub release.
