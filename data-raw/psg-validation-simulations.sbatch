#!/usr/bin/env bash
#SBATCH --job-name=recur-psg
# 6 scenarios x 10,000 Monte-Carlo replicates are split into 6,000 chunks of
# 10 scenario-replicates.  Adjust the concurrency cap to cluster policy.
#SBATCH --array=1-6000%100
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=data-raw/slurm-logs/psg-%A_%a.out
#SBATCH --error=data-raw/slurm-logs/psg-%A_%a.err

set +u
source ~/.bash_profile
set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"; mkdir -p data-raw/slurm-logs
export RUN_PSG_VALIDATION=true PSG_N_SIM=10000 PSG_B=999
export PSG_TASK_ID="${SLURM_ARRAY_TASK_ID}" PSG_N_TASKS="${SLURM_ARRAY_TASK_COUNT}"
Rscript data-raw/psg-validation-simulations.R
