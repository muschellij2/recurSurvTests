#!/usr/bin/env bash
#SBATCH --job-name=ss-rcif-tables
#SBATCH --array=1-24%4
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=48:00:00
#SBATCH --output=data-raw/slurm-logs/ss-rcif-%A_%a.out
#SBATCH --error=data-raw/slurm-logs/ss-rcif-%A_%a.err

set +u
source ~/.bash_profile
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p data-raw/slurm-logs
export RUN_SS_RCIF_TABLE_SIMULATION=true
export SS_RCIF_TASK_ID="${SLURM_ARRAY_TASK_ID}"
export SS_RCIF_REPS="${SS_RCIF_REPS:-1000}"
Rscript data-raw/ss-rcif-table-simulation.R
