#!/usr/bin/env bash
#SBATCH --job-name=recur-ci
#SBATCH --array=1-20%10
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --output=data-raw/slurm-logs/curve-ci-%A_%a.out
#SBATCH --error=data-raw/slurm-logs/curve-ci-%A_%a.err
set +u
source ~/.bash_profile
set -euo pipefail

cd "${SLURM_SUBMIT_DIR}"; mkdir -p data-raw/slurm-logs
export RUN_CURVE_CI_COVERAGE=true CURVE_CI_N_SIM=1000 CURVE_CI_B=999 CURVE_CI_TRUTH_N=25000
export CURVE_CI_TASK_ID="${SLURM_ARRAY_TASK_ID}" CURVE_CI_N_TASKS="${SLURM_ARRAY_TASK_COUNT}"
Rscript data-raw/curve-ci-coverage-simulation.R
