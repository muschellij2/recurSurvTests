#!/usr/bin/env bash
#SBATCH --job-name=zhao2020
#SBATCH --array=1-120%20
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=data-raw/slurm-logs/zhao-%A_%a.out
#SBATCH --error=data-raw/slurm-logs/zhao-%A_%a.err

# Zhao et al. (2020) Table 1--3 simulation array.
#
# Submit from the package root:
#   sbatch data-raw/zhao-2020-table-simulations.sbatch
#
# The 120 array tasks correspond to the 12 Table 1, 36 Table 2, and 72 Table 3
# cells. `%20` limits the scheduler to 20 concurrent array tasks; change it to
# suit the allocation. Every task uses four local R workers and writes one
# task-specific RDS result. Increase --time if your cluster is slower than the
# benchmarked host, or lower --cpus-per-task if the partition limits cores.

set -euo pipefail

source ~/.bash_profile
cd "${SLURM_SUBMIT_DIR}"
mkdir -p data-raw/slurm-logs
mkdir -p data-raw/zhao-2020-table-simulations

export RUN_ZHAO_100K=true
export ZHAO_PARALLEL=multicore
export ZHAO_N_CORES="${SLURM_CPUS_PER_TASK}"
export ZHAO_TASK_ID="${SLURM_ARRAY_TASK_ID}"
export ZHAO_N_TASKS="${SLURM_ARRAY_TASK_COUNT}"
export ZHAO_CHUNK_SIZE=100

# Use the R installation configured on the cluster. If your site requires a
# module, add its `module load R/...` command above this line.
Rscript data-raw/zhao-2020-table-simulations.R
