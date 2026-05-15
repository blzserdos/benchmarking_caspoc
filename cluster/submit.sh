#!/bin/bash
# =============================================================================
# submit.sh
# SLURM array job submission script for CASPOC benchmarking
#
# Usage:
#   cd benchmarking_caspoc
#   sbatch cluster/submit.sh
#
# Adjust CHUNK_SIZE and --array range to control granularity:
#   Total jobs = N_ITERATIONS * N_APPROACHES * (N_PERM + 1) * N_DATASETS
#             = 100 * 4 * 101 * 2 = 80,800
#   Array tasks = ceil(80800 / CHUNK_SIZE)
#
# With CHUNK_SIZE=400: ceil(80800/400) = 202 array tasks
# =============================================================================

#SBATCH --account=nn9114k
#SBATCH --job-name=caspoc_bench
#SBATCH --array=1-202
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=cluster/logs/slurm_%A_%a.out
#SBATCH --error=cluster/logs/slurm_%A_%a.err

# --- Configuration ---
CHUNK_SIZE=400

# --- Setup ---
module load R/4.5.2-gfbf-2025b

# Create log directory if needed
mkdir -p cluster/logs

# --- Run ---
echo "Starting array task ${SLURM_ARRAY_TASK_ID} at $(date)"
echo "Hostname: $(hostname)"

Rscript cluster/run_chunk.R ${SLURM_ARRAY_TASK_ID} ${CHUNK_SIZE}

echo "Finished array task ${SLURM_ARRAY_TASK_ID} at $(date)"
