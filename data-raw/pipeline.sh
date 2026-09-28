#!/usr/bin/env bash
# Full simulation and report pipeline for recurSurvTests.
# Run from the repository root on a Slurm cluster with R and package
# dependencies installed. Arrays run synchronously (--wait) so aggregation
# starts only after every task in that array has completed successfully.
set -euo pipefail

if [[ ! -f DESCRIPTION || ! -d data-raw ]]; then
  echo "Run this script from the recurSurvTests repository root." >&2
  exit 1
fi
if ! command -v sbatch >/dev/null 2>&1; then
  echo "Slurm sbatch is required. This pipeline submits the production arrays." >&2
  exit 1
fi
if ! command -v Rscript >/dev/null 2>&1; then
  echo "Rscript is required." >&2
  exit 1
fi

mkdir -p data-raw/slurm-logs

submit_wait() {
  local label="$1"
  shift
  echo "Submitting ${label}: $*"
  sbatch --wait --parsable "$@"
  echo "Finished ${label}."
}

# PSG-shaped two-arm/paired scenarios (6,000 chunks; 10 scenario-replicates
# per chunk) and the original three-group K validation.
submit_wait "PSG validation" --array=1-6000%100 data-raw/psg-validation-simulations.sbatch
Rscript data-raw/psg-validation-aggregate.R

submit_wait "K-group validation" --array=1-20%10 data-raw/wc-logrank-k-validation-simulation.sbatch
Rscript data-raw/wc-logrank-k-validation-aggregate.R

# Zhao table study, using identical task seeds for the two variance options.
submit_wait "Zhao pooled-risk table simulation" \
  --export=ALL,ZHAO_VARIANCE_METHOD=pooled_risk --array=1-120%20 \
  data-raw/zhao-2020-table-simulations.sbatch
Rscript data-raw/zhao-2020-validation-report.R

submit_wait "Zhao printed Equation 6 table simulation" \
  --export=ALL,ZHAO_VARIANCE_METHOD=zhao_eq6 --array=1-120%20 \
  data-raw/zhao-2020-table-simulations.sbatch
Rscript data-raw/zhao-2020-validation-report.R \
  data-raw/zhao-2020-table-simulations-zhao_eq6 \
  docs/zhao-validation-zhao-eq6
Rscript data-raw/zhao-variance-comparison.R

# Newer diagnostic/power designs, kept in their own data-raw subfolders.
submit_wait "K-group null calibration diagnostics" \
  --array=1-4%4 data-raw/wc-logrank-k-diagnostics/simulation.sbatch
Rscript data-raw/wc-logrank-k-diagnostics/aggregate.R

submit_wait "Two-group wc_logrank power grid" \
  --array=1-20%10 data-raw/wc-logrank-power/simulation.sbatch
Rscript data-raw/wc-logrank-power/aggregate.R

# Curve confidence-interval coverage is a separate estimator validation. Each
# array task writes task-specific RDS/CSV outputs; there is currently no
# dedicated aggregate script for these files.
submit_wait "Curve confidence-interval coverage" \
  --array=1-20%10 data-raw/curve-ci-coverage-simulation.sbatch

echo "All simulation arrays and available aggregators completed."
echo "Review the reports and task outputs before treating the simulations as validated."

# Optional downstream reproducibility checks; enable explicitly, e.g.
# RUN_PACKAGE_CHECK=1 RUN_MANUSCRIPT_RENDER=1 bash data-raw/pipeline.sh
if [[ "${RUN_PACKAGE_CHECK:-0}" == "1" ]]; then
  Rscript -e 'roxygen2::roxygenise(); devtools::check(args = c("--no-manual", "--no-build-vignettes"))'
fi
if [[ "${RUN_MANUSCRIPT_RENDER:-0}" == "1" ]]; then
  quarto render vignettes/jss-manuscript.qmd
fi
