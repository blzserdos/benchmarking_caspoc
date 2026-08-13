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

Real omics blocks are not collections of independent features: they carry a few
dominant co-expression axes. Both simulated datasets therefore apply an
**overlapping factor structure**, where each feature loads on several latent
factors:

```
N_j = Σ_k a_jk F_k + √ψ_j ε_j        F_k, ε_j ~ N(0,1) iid
```

so `cor(i,j) = Σ_k a_ik a_jk` takes a continuum of values. Every feature is
rescaled to its own target communality, so marginal variance stays exactly 1 and
`signal_strength` means the same thing regardless of structure.

| parameter | default | controls |
|---|---|---|
| `K` | 20 | number of latent factors |
| `gamma` | 1.5 | spectrum slope — factor weights decay as `k^-gamma` |
| `pi_join` | 0.5 | membership density — how many factors a feature loads on |
| `h2_mean` | 0.6 | communality — share of a feature's variance that is shared |

**This is not cosmetic.** Nuisance factors are independent *across* blocks but
high variance *within* one, so at n=100 sPLS is readily distracted by spurious
factor-to-factor alignments. Structured noise shifts the detectability threshold
right substantially relative to independent features. Set `BLOCK_STRUCTURE` to
`NULL` to recover the original independent-feature simulation.

The goal is **omics-like** data, not a replica of any dataset. The values above
are deliberately round; a real block is used only to confirm the right ballpark:

| | mean\|c\| | RMS | median \|c\| | % negative | spectrum RMSE |
|---|---|---|---|---|---|
| real `breast.TCGA` mRNA | 0.200 | 0.249 | 0.173 | 44.0% | — |
| **overlapping (default)** | 0.198 | 0.248 | 0.162 | 50.0% | 0.67 |
| disjoint modules (previous) | 0.134 | 0.235 | 0.066 | 43.2% | 0.91 |

Tighter fitting would be false precision — draw-to-draw variation in the
covariance is as large as the effect of the parameters themselves.

#### Why this replaced disjoint modules

The previous structure partitioned features into modules, each driven by one
shared factor. It fitted the eigenvalue spectrum but had three artifacts it
could not represent its way out of:

- **86.3% of feature pairs were exactly zero**, so a wrongly-selected feature was
  unrelated to everything else; in real data it is still entangled with the rest
- **negative correlations were impossible** — all loadings on a module factor
  share a sign — which matters for compositional (CLR-transformed) microbiome data
- **only `n_modules` directions of variation existed**, so the spectrum fell off
  a cliff where real blocks have a long tail of moderate directions

Since the difficulty here is sPLS being distracted by directions competing with
the true signal, all three plausibly affect how the CV methods behave.

`pi_join` is a continuum rather than a switch: at `pi_join ≈ 1/K` each feature
loads on roughly one factor, which *is* the disjoint structure. The two are
endpoints of one dial, which makes robustness across them a sweep rather than a
two-point comparison.

### Interpreting `signal_strength`

With unit-norm loadings and independent noise, the oracle cross-block
correlation is exactly `s/(s+1)`, and each relevant feature carries
`0.05s/(1+0.05s)` of its variance on the shared axis (20 relevant X features).
Under `BLOCK_STRUCTURE` this becomes `s/(s + l'Σl)`, since the noise projected
onto the loading direction no longer has unit variance. Measured under the
overlapping structure, the identity holds closely:

| s | `s/(s+1)` | observed |
|---|---|---|
| 1 | 0.500 | 0.508 ± 0.051 |
| 2 | 0.667 | 0.671 ± 0.045 |
| 4 | 0.800 | 0.802 ± 0.033 |

so a power curve indexed by `s` means what it says.

The measured power transition (200 replicates, proxy statistic):

| s | 1 | 2 | 3 | 4 | 6 | 8 | 12 | 20 |
|---|---|---|---|---|---|---|---|---|
| power | 0.06 | 0.07 | 0.28 | 0.53 | 0.75 | 0.89 | 0.97 | 1.00 |

For scale: applying the same CV statistic to real `breast.TCGA` block pairs
(subsampled to n=100) calibrates them to **`s ≈ 37–48`** — far above the range
where the four methods differ in power. The sweep therefore characterises the
**hard regime**, not the regime real omics data occupies, where all four methods
saturate at power 1. This is a legitimate claim but must be stated that way.

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
