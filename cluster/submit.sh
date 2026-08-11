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
#   Rscript -e 'source("cluster/config.R"); print_grid_summary(100)'
#
# Current config (sim_null without perms, sim_signal with N_PERM=100 across 7
# signal strengths) gives: 100*4*1 + 100*4*101*7 = 283,200 jobs
# -> 2,832 array tasks at CHUNK_SIZE=100.
#
# CHUNK_SIZE sizing. Measured job runtimes from the previous run: naive_cv 4.3s,
# repeated_cv 46.5s, caspoc 110.3s, nested_cv 222.7s -> mean 95.6s, SD ~84s.
# The grid is shuffled, so a chunk is a random mix and its walltime is roughly
# CHUNK_SIZE * 95.6s +/- sqrt(CHUNK_SIZE) * 84s:
#
#   CHUNK_SIZE=200 -> 5.31h mean, 5.97h at +2SD   <- exceeded the old 6h limit;
#                                                    cost 10 whole chunks
#   CHUNK_SIZE=100 -> 2.66h mean, 3.60h at +4SD   <- comfortable headroom
#
# A timed-out task loses its entire chunk (run_chunk.R only writes at the end),
# so keep the headroom.
#
# NOTE: 2,832 tasks exceeds the MaxArraySize on many clusters (often 1001).
# Check with:  scontrol show config | grep MaxArraySize
#
# If it is too low, submit one array per signal strength — config.R honours a
# 0-based STRENGTH_INDEX environment variable (408 tasks for index 0, which also
# carries sim_null; 404 for the rest):
#
#   sbatch --array=1-408 --export=ALL,STRENGTH_INDEX=0 cluster/submit.sh
#   for i in $(seq 1 6); do
#     sbatch --array=1-404 --export=ALL,STRENGTH_INDEX=$i cluster/submit.sh
#   done
#
# Leave STRENGTH_INDEX unset when running collect_results.R — it needs the full
# grid to verify completeness.
# =============================================================================

#SBATCH --account=nn9114k
#SBATCH --job-name=caspoc_bench
#SBATCH --array=1-2832
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=06:00:00
#SBATCH --output=cluster/logs/slurm_%A_%a.out
#SBATCH --error=cluster/logs/slurm_%A_%a.err

# --- Configuration ---
CHUNK_SIZE=100

# --- Setup ---
module load R/4.5.2-gfbf-2025b

# Create log directory if needed
mkdir -p cluster/logs

# --- Run ---
echo "Starting array task ${SLURM_ARRAY_TASK_ID} at $(date)"
echo "Hostname: $(hostname)"

Rscript cluster/run_chunk.R ${SLURM_ARRAY_TASK_ID} ${CHUNK_SIZE}

echo "Finished array task ${SLURM_ARRAY_TASK_ID} at $(date)"
