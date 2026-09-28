#!/usr/bin/env bash
#SBATCH --job-name=wc-power
#SBATCH --array=1-20%10
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=data-raw/slurm-logs/wc-power-%A_%a.out
#SBATCH --error=data-raw/slurm-logs/wc-power-%A_%a.err

set +u
source ~/.bash_profile
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p data-raw/wc-logrank-power/slurm-logs
export RUN_WC_POWER_SIMULATION=true
export WC_POWER_N_SIM=5000
export WC_POWER_N_TASKS="${SLURM_ARRAY_TASK_COUNT}"
export WC_POWER_TASK_ID="${SLURM_ARRAY_TASK_ID}"
Rscript data-raw/wc-logrank-power/simulate.R
