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
├── run_benchmarks.R          # Main analysis script (entry point)
├── R/
│   ├── generate_data.R       # Simulated and real data generators
│   ├── cv_approaches.R       # Uniform wrappers for all 4 CV methods
│   ├── permutation_test.R    # Permutation-based significance testing
│   └── evaluate_results.R    # Summarisation and reporting utilities
├── results/                  # Output .rds files (raw + summary)
├── data/                     # Real datasets (if applicable)
└── figures/                  # Generated plots
```

## Usage

```r
# From the benchmarking_caspoc/ directory:
Rscript run_benchmarks.R
```

Or interactively:

```r
source("run_benchmarks.R")
run_all_benchmarks()
```

Iterations are parallelized across cores using the `future`/`furrr` framework. Set `N_CORES <- 1` in the config section to run sequentially.

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

- `results/benchmark_raw_results.rds` — one row per approach × iteration with observed statistic and permutation p-value
- `results/benchmark_summary.rds` — aggregated: mean/SD of test statistics, rejection rates at α = 0.05 (labeled "FPR" for null data, "Power" for signal data)
