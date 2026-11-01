# What does a ridge network estimate? Reproducibility archive

Code, results and figures for the article

> Arashi, M. and Ebrahim, E. K. *What does a ridge network estimate? Effective penalty and leakage in ridge precision matrix estimation.*

The archive also contains the R package **ridgegraph** 0.1.1, which implements the ridge partial correlations, the estimated effective penalty and its limit as the penalty tends to zero, the two-point extrapolated estimator and the multi-penalty permutation selection.

## Contents

| Folder / file | What it holds |
|---|---|
| `code/sim_ridge_support.R` | Main simulation study (p = 100, four graphs, n = 50, 100, 200, 400, 100 replicates), the effective-penalty grid (band model, p = 200) and the timings |
| `code/sim_extensions.R` | Further simulations: the main design rerun with the same random numbers and extra quantities stored (signed partial correlations, population margins, further competitors), the stored Table 1 replicates, and the p = 400 block |
| `code/sim_summaries.R` | Turns the stored replicates into the tables of the extensions (Monte Carlo standard errors, family-wise error rates, margins, competitors, p = 400) |
| `code/table1_fixed_point.R` | The two deterministic columns of Table 1: the fixed-point effective penalty and the Ky Fan floor |
| `code/data_stock_network.R` | Stock-return analysis (yearly windows and full sample) |
| `code/data_extensions.R` | Further data analyses: seed variability of the permutation selection, factor adjustment, block permutations, leakage rates, rank correlations, cleaning check |
| `code/sim_main2.R`, `code/sim_main2_summary.R` | The main design extended to 500 replicates, with GeneNet, beam, parcor (replicates 1 to 100), the default-target van Wieringen-Peeters rule, matched-size and matched-FDR comparisons, penalty choice, the leakage cost, and the extrapolation checks |
| `code/sim_designs2.R`, `code/sim_designs2_summary.R`, `code/run_designs2.sh` | Further designs: heterogeneous mixed-sign weights at p = 100, large n, a factor model and heavy-tailed data |
| `code/sim_rescore.R` | Rescores the lasso paths of the further designs as in the main design, the factor-adjusted lasso comparison, the tail-area rules (Benjamini-Hochberg and tail-area Fdr) on the ridge partial correlations, and three constants of the text; rebuilds each replicate from its stored random numbers |
| `code/data_checks.R` | Further stock-data analyses: matched effective penalties, subsamples, a leakage diagnostic, retuned lasso methods, GeneNet and beam, GARCH standardization and a semi-synthetic check |
| `code/fig_headline.R` | Figure 1 (effective penalty against nominal penalty, and leakage by graph distance); deterministic, called by `make_figures.R` |
| `code/make_figures.R` | Draws all figures from the saved results |
| `results/` | All outputs (CSV and RDS) and the session information of every run (`sessionInfo_*.txt`) |
| `results/ext_work/`, `results/ext_work2/`, `results/ext_work3/`, `results/ext_work4/`, `results/ext_cache/`, `results/checks_cache/` | Stored replicates and work units; with them the tables can be rebuilt without rerunning the simulations |
| `figures/` | The figures of the article and the supplement (PDF, fonts embedded) |
| `ridgegraph/`, `ridgegraph_0.1.1.tar.gz` | Source and built package |

## Which script gives which result

| Result in the article | Script | Output files |
|---|---|---|
| Figure 1 (headline figure) | `fig_headline.R` | `sim_kappa_grid.csv` and deterministic computations |
| Table 2 (effective penalty, band model, p = 200) | `sim_ridge_support.R`, `sim_extensions.R kgrid`, `table1_fixed_point.R`, `sim_main2.R` | `sim_kappa_grid.csv`, `sim_kappa_grid_reps.csv`, `sim_kappa_fixedpoint.csv`, `sim_main2_kgrid_relerr.csv` |
| Table 3 (leakage and noise; main design and p = 400 block) | `sim_ridge_support.R`, `sim_extensions.R main highdim`, `sim_summaries.R`, `sim_main2.R` | `sim_kappa.csv`, `sim_margin.csv`, `sim_highdim.csv`, `sim_mcse.csv`, `sim_main2_leakage_oracle*.csv` |
| Table 4 (selection, main design, 500 replicates, with further competitors and a tail-area rule) | `sim_main2.R`, `sim_main2_summary.R`, `sim_rescore.R` | `sim_main2_selection.csv`, `sim_main2_fdr.csv`, `sim_main2_multiplicity*.csv`, `sim_main2_tailrules.csv` |
| Table 5 (rankings, main design, 500 replicates) | `sim_main2.R`, `sim_main2_summary.R` | `sim_main2_ranking.csv`, `sim_main2_tuning.csv` |
| Table 6 (further designs: heterogeneous weights, large n, factor model, heavy tails) | `sim_designs2.R`, `sim_rescore.R`, `sim_designs2_summary.R` (or `run_designs2.sh`) | `sim_designs2_*.csv` (lasso paths: `sim_designs2_paths.csv`, factor-adjusted comparison: `sim_designs2_factor_fair.csv`) |
| Tables 7 and 8 (stock returns) | `data_stock_network.R`, `data_checks.R` (some sentences of Section 6 use `data_extensions.R`) | `data_kappa.csv`, `data_methods.csv`, `data_matched.csv`, `data_checks_*.csv` |
| Supplement, simulation sections | `sim_ridge_support.R`, `sim_summaries.R`, `sim_main2_summary.R`, `sim_designs2_summary.R` | `sim_mcse.csv`, `sim_claims.csv`, `sim_fwer.csv`, `sim_competitors.csv`, `sim_highdim.csv`, `sim_leakage.csv`, `sim_pairs_*.csv`, `sim_sensitivity.csv`, `sim_main2_*.csv`, `sim_designs2_*.csv`, `sim_timing.csv` (machine dependent) |
| Supplement, data sections | `data_stock_network.R`, `data_extensions.R`, `data_checks.R` | `data_kappa_curve.csv`, `data_marginal.csv`, `data_components_W1.csv`, `data_mps_B999*.csv`, `data_factor*.csv`, `data_blockperm*.csv`, `data_leakage_rates.csv`, `data_spearman.csv`, `data_cleaning*.csv`, `data_checks_*.csv` |
| All figures | `make_figures.R` | `figures/*.pdf` |

## How to run

R 4.1 or later with the CRAN packages huge, glmnet, fdrtool, igraph, GGMridge, SILGGM, GeneNet, beam and parcor (and rugarch for the GARCH check). Install the package first:

```
R CMD INSTALL ridgegraph_0.1.1.tar.gz
```

Then run the scripts from the `code/` folder, in this order:

```
Rscript sim_ridge_support.R
Rscript sim_extensions.R
Rscript sim_summaries.R
Rscript table1_fixed_point.R
Rscript data_stock_network.R
Rscript data_extensions.R
Rscript sim_main2.R
Rscript sim_main2_summary.R
sh run_designs2.sh
Rscript sim_rescore.R
Rscript sim_designs2_summary.R
Rscript data_checks.R
Rscript make_figures.R
```

If the packages sit in a separate library folder, give it in the environment variable `RIDGE_RLIB`.

**Random numbers.** Every script sets its seeds. The main simulation uses two L'Ecuyer random-number streams per setting, built explicitly, so its results do not depend on the number of cores or on the operating system. The extensions reuse the same streams and keep any extra random draws apart, so they reproduce the published replicates exactly. The data analysis fixes the seed of the permutation selection (999 permutations) and of neighbourhood selection in every window. `sim_summaries.R` needs the stored replicates in `results/ext_work/`.

**Package versions.** The results of neighbourhood selection depend on the version of huge; the runs reported in the article used huge 1.6. The exact versions of every run are in `results/sessionInfo_*.txt`.

**Running time.** On a 12-core desktop (22 parallel workers) the rerun of the main design in `sim_extensions.R` took about 2.5 hours and the p = 400 block about one hour (mostly the graphical lasso path); the data analysis with 999 permutations takes about ten minutes per window on one core. The 500-replicate main design in `sim_main2.R` took about 12 hours on 11 workers, the further designs about 10.5 hours on 11 workers, and `data_checks.R` about 2.5 hours on 3 to 4 cores. With `results/ext_work/` present, `sim_summaries.R` rebuilds the tables without rerunning the simulations.

## Data

The stock data are the `stockdata` object of the CRAN package huge (daily closing prices of 452 S&P 500 stocks, 2003 to 2008, with their GICS sectors). No other data are used.

## Licence

Code: MIT licence (see `LICENSE`). Results and figures: CC BY 4.0.
