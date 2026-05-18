#!/bin/bash
# =============================================================================
# submit.sh
# SLURM array job submission script for CASPOC benchmarking
#
# Usage:
#   cd benchmarking_caspoc
#   sbatch cluster/submit.sh
#
# Adjust CHUNK_SIZE and --array range to control granularity.
#
# Total jobs depend on cluster/config.R — datasets can opt out of permutation
# testing (run_perm = FALSE) and signal datasets sweep over SIGNAL_STRENGTHS.
# Get the exact number after editing config with:
#   Rscript -e 'source("cluster/config.R"); print_grid_summary(200)'
#
# Defaults (sim_null without perms, sim_signal with N_PERM=100 and one
# signal strength) give: 100*4*1 + 100*4*101*1 = 40,800 jobs
# -> 204 array tasks with CHUNK_SIZE=200.
# =============================================================================

#SBATCH --account=nn9114k
#SBATCH --job-name=caspoc_bench
#SBATCH --array=1-204
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --output=cluster/logs/slurm_%A_%a.out
#SBATCH --error=cluster/logs/slurm_%A_%a.err

# --- Configuration ---
CHUNK_SIZE=200

# --- Setup ---
module load R/4.5.2-gfbf-2025b

# Create log directory if needed
mkdir -p cluster/logs

# --- Run ---
echo "Starting array task ${SLURM_ARRAY_TASK_ID} at $(date)"
echo "Hostname: $(hostname)"

Rscript cluster/run_chunk.R ${SLURM_ARRAY_TASK_ID} ${CHUNK_SIZE}

echo "Finished array task ${SLURM_ARRAY_TASK_ID} at $(date)"
