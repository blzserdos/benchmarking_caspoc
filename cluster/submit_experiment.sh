#!/bin/bash
# =============================================================================
# submit_experiment.sh
# Runs the pooled-null design (run_experiment.R) as ONE multi-core job.
#
#   cd benchmarking_caspoc
#   sbatch cluster/submit_experiment.sh
#
# NOT an array job. cluster/submit.sh needs 2,832 array tasks because the
# permutation design is 283,200 jobs (~7,500 CPU-hours). This design is 3,200
# jobs -- roughly 85 CPU-hours at the runtimes measured on Saga (naive_cv 4.3s,
# repeated_cv 46.5s, caspoc 110.3s, nested_cv 222.7s; mean ~96s). On 32 cores
# that is under 3 hours of wall clock, which fits comfortably in one job.
#
# run_experiment.R checkpoints to results/experiment_raw.rds every 200 jobs and
# skips completed work on restart, so a timeout is not a disaster: resubmit the
# same script and it picks up where it stopped.
#
# IMPORTANT: leave STRENGTH_INDEX unset. cluster/config.R honours it and would
# silently reduce the sweep to a single signal strength.
# =============================================================================

#SBATCH --account=nn9114k
#SBATCH --job-name=caspoc_experiment
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=cluster/logs/experiment_%j.out
#SBATCH --error=cluster/logs/experiment_%j.err

set -euo pipefail

module load R/4.5.2-gfbf-2025b

mkdir -p cluster/logs results

# run_experiment.R reads SLURM_CPUS_PER_TASK for its worker count. detectCores()
# would report every core on the node rather than the 32 allocated here.
echo "Job ${SLURM_JOB_ID} on $(hostname), ${SLURM_CPUS_PER_TASK} cores, $(date)"

unset STRENGTH_INDEX

Rscript run_experiment.R

echo "Finished at $(date)"
