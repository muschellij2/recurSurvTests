#!/usr/bin/env bash
#SBATCH --job-name=wrs-kgroup
#SBATCH --array=1-20%10
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=data-raw/slurm-logs/wrs-k-%A_%a.out
#SBATCH --error=data-raw/slurm-logs/wrs-k-%A_%a.err
set +u
source ~/.bash_profile
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"; mkdir -p data-raw/slurm-logs
export RUN_WC_LOGRANK_K_VALIDATION=true WC_LOGRANK_K_N_SIM=20000
export WC_LOGRANK_K_TASK_ID="${SLURM_ARRAY_TASK_ID}" WC_LOGRANK_K_N_TASKS="${SLURM_ARRAY_TASK_COUNT}"
Rscript data-raw/wc-logrank-k-validation-simulation.R
