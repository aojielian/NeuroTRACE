# Reproduction commands

Work from the repository root. The repository does not automatically download biological data. No raw records or individual-level scores are included.

Install the existing Python environment from `environment/requirements.txt`. Figure rendering additionally uses matplotlib, pandas, NumPy and Pillow. R plotting uses ggplot2, patchwork, cowplot, dplyr, tidyr, readr, scales, ggrepel, ragg and systemfonts. Calibration uses Matrix and data.table. Render from the frozen source tables:

```sh
export NEUROTRACE_PROJECT_ROOT="$PWD"
Rscript scripts/figures/render_submission_figures.R
python scripts/figures/render_structure_figures.py
Rscript scripts/calibration/update_MDD_CI.R
```

To reproduce the identical final matched nulls, select a separate output directory without overwriting frozen tables:

```sh
export NEUROTRACE_OUTPUT_RUN="$PWD/reproduced_calibration"
Rscript scripts/calibration/unified_matcher.R
```

Compare outputs with `processed_results/unified_matching_sensitivity_results.tsv` and `processed_results/unified_matching_null_replicates.tsv.gz`. Do not interpret null median/percentile ranges as observed-effect uncertainty. The six-stage profiles and membership are fixed; only the two postnatal stage names differ from historical identifiers.
