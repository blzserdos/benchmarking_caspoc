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

# 0. Confirm the grid size after editing config
Rscript -e 'source("cluster/config.R"); print_grid_summary(100)'

# 1. Submit array job (2,832 tasks, each processing 100 jobs)
#    Check your cluster's cap first: scontrol show config | grep MaxArraySize
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
- `BLOCK_STRUCTURE` — within-block correlation (see below); `NULL` for independent features
- `SIGNAL_STRENGTHS` — signal strengths swept on signal datasets

## Datasets

- **sim_null**: No association between X and Y (n=100, p=200, q=50). Used to assess Type I error (false positive rate).
- **sim_signal**: Shared latent structure with sparse loadings (20 relevant X, 10 relevant Y), swept over `SIGNAL_STRENGTHS`. Used to assess power.
- Real datasets (breast TCGA, microbiome?).

### Within-block correlation

Both simulated datasets apply **block-diagonal compound symmetry with
heterogeneous module strengths**. Module `k` has size ∝ `k^-alpha` and
within-module correlation `rho * k^-beta`; correlation across modules is zero.
It is realised as one latent factor per module, so every feature keeps unit
marginal variance and `signal_strength` means the same thing regardless of
structure.

Defaults (`n_modules=8, alpha=0.5, rho=0.8, beta=0.5`) were fitted to the
`breast.TCGA` mRNA block by minimising RMSE over the first 20 eigenvalue
fractions at matched n and p:

| PC | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|---|
| real (%) | 18.9 | 13.1 | 5.5 | 5.4 | 4.2 | 3.3 | 2.4 | 2.1 |
| sim (%) | 19.2 | 10.0 | 6.6 | 5.2 | 4.1 | 3.6 | 3.1 | 2.6 |

`beta > 0` matters more than it looks. Equal-strength modules (`beta = 0`)
produce a few equal eigenvalue spikes and then a cliff — 3 modules gives
19.5/14.7/12.6/**1.2**/1.2/… — which leaves sPLS facing several equally
attractive spurious directions. That idealisation is measurably pessimistic:
power at `signal_strength = 4` was 0.49 under equal modules versus 0.85 under
the fitted decaying spectrum.

The fit prioritises the eigenvalue spectrum, which is what governs sPLS's
behaviour, and pays for it elsewhere: simulated mean |correlation| is 0.135
against 0.200 in the real block. The structure also has no negative
correlations, which is a real gap for compositional (CLR-transformed)
microbiome data.

None of this is cosmetic. Module factors are independent *across* blocks but
high variance *within* one, so at n=100 sPLS is readily distracted by spurious
factor-to-factor alignments. Adding the structure shifts the detectability
threshold right by roughly 3–4× relative to independent noise.

### Interpreting `signal_strength`

With unit-norm loadings and independent noise, the oracle cross-block
correlation is exactly `s/(s+1)`, and each relevant feature carries
`0.05s/(1+0.05s)` of its variance on the shared axis (20 relevant X features).
Under `BLOCK_STRUCTURE` this becomes `s/(s + l'Σl)`, so it depends on
`SIGNAL_ALIGNMENT`:

| | oracle cor at s=1 | power at s=4 |
|---|---|---|
| `spread` (default) | 0.50 ± 0.03 | 0.38 |
| `aligned` | 0.58 ± 0.15 | 0.72 |

`spread` is the default precisely because it preserves the identity, so a power
curve indexed by `s` is interpretable.

For scale: applying the same CV statistic to real `breast.TCGA` block pairs
(subsampled to n=100) calibrates them to `s ≈ 28–51` — far above the range
where the four methods differ in power (the transition spans roughly `s = 2` to
`s = 10`). The sweep therefore characterises the **hard regime**, not the regime
real omics data occupies, where all four methods saturate at power 1.

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
