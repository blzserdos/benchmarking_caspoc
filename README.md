## Methods Compared

| Approach | Description | Bias |
|----------|-------------|------|
| Naïve CV | Single k-fold pass; selects hyperparameters and evaluates on the same folds | Biased (optimistic) |
| Repeated CV | Multiple repeats of k-fold; median across repeats, but still uses same data for tuning and evaluation | Biased (optimistic) |
| Nested CV | Outer loop for evaluation, inner loop for tuning (repeated for stability) | Unbiased |
| CASPOC | Double-split k-fold with circular fold assignment (tune → test → train) and consensus over odd repeats | Unbiased |

## Significance Testing

All methods use a uniform permutation-based approach: the full CV pipeline is re-run on `N_PERM` permuted datasets (Y rows shuffled) to build a null distribution. The empirical p-value uses the conservative formula from Phipson & Smyth (2010): `(n_extreme + 1) / (n_valid + 1)`.

## Project Structure

```
benchmarking_caspoc/
├── run_benchmarks.R          # Local analysis script (uses future/furrr)
├── R/
│   ├── generate_data.R       # Simulated and real data generators
│   ├── cv_approaches.R       # Uniform wrappers for all 4 CV methods
│   ├── permutation_test.R    # Permutation-based significance testing
│   └── evaluate_results.R    # Summarisation and reporting utilities
├── cluster/                  # SLURM cluster execution
│   ├── config.R              # Shared configuration (grid, datasets, HP settings)
│   ├── submit.sh             # SLURM array job submission script
│   ├── run_chunk.R           # Worker script (processes a chunk of the job grid)
│   ├── collect_results.R     # Post-processing (combines chunks, computes p-values)
│   ├── chunks/               # Partial results from each array task
│   └── logs/                 # SLURM stdout/stderr logs
├── results/                  # Final output .rds files
├── data/                     # Real datasets (if applicable)
└── figures/                  # Generated plots
```

## Usage

### Local (laptop/workstation)

```r
# From the benchmarking_caspoc/ directory:
Rscript run_benchmarks.R
```

Or interactively:

```r
source("run_benchmarks.R")
run_all_benchmarks()
```

Uses the `future`/`furrr` framework for multi-core parallelism. Set `N_CORES <- 1` for sequential execution.

### Cluster (SLURM)

```bash
cd benchmarking_caspoc

# 1. Submit array job (404 tasks, each processing 200 jobs)
sbatch cluster/submit.sh

# 2. Monitor progress
squeue -u $USER
ls cluster/chunks/ | wc -l   # completed chunks

# 3. After all tasks finish, collect and summarise
Rscript cluster/collect_results.R
```

Edit `cluster/config.R` to change simulation parameters. Edit `cluster/submit.sh` to adjust SLURM resources or chunk size.

## Configuration

Key parameters in `run_benchmarks.R`:

- `N_ITERATIONS = 100` — simulation replicates per dataset (increase for final paper)
- `N_PERM = 100` — permutations per iteration
- `N_CORES` — number of parallel workers (default: all cores minus one)
- `CV_CONFIG` — shared settings: number of folds (10), repeats (11), inner folds (5)
- `HP_GRID` — keepX/keepY sparsity options to search over

## Datasets

- **sim_null**: Independent X and Y (n=100, p=200, q=50). Used to assess Type I error (false positive rate).
- **sim_signal**: Shared latent structure with sparse loadings (20 relevant X, 10 relevant Y). Used to assess power.
- Real datasets (breast TCGA, microbiome?).

## Dependencies

```r
install.packages(c("mixOmics", "caret", "dplyr", "tibble", "MASS",
                   "future", "furrr", "parallelly"))

# CASPOC package (from GitHub):
# devtools::install_github("jonathanth/caspoc", ref = "CASPOC-v1.1")
```

## Output

- `results/<dataset>_flat_results.rds` — one row per job (iteration × approach × perm_id), the full raw data
- `results/benchmark_results.rds` — one row per (iteration × approach) with observed statistic and permutation p-value
- `results/benchmark_summary.rds` — aggregated: mean/SD of test statistics, rejection rates at α = 0.05 (labeled "FPR" for null data, "Power" for signal data)
