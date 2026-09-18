# R code for *High-Dimensional Mean Testing under a Controlled Dense Missing Pattern*

This repository contains the R code accompanying the article **“High-Dimensional Mean Testing under a Controlled Dense Missing Pattern.”** The scripts implement simulation studies and an EEG-data illustration for comparing the proposed high-dimensional mean test with the test of Yin under missing observations.

## Repository contents

| File | Description |
|---|---|
| [`lighttail_revised.R`](lighttail_revised.R) | Simulates a light-tailed Gaussian linear process, applies a dense MCAR missingness pattern, and compares the empirical rejection probabilities of the proposed and Yin tests. |
| [`linear-example2-simulation-revised.R`](linear-example2-simulation-revised.R) | Runs the heavy-tailed linear-process simulation from Example 2. Innovations follow a Gaussian–symmetric-Pareto mixture, and missingness varies across coordinates. |
| [`garch-simulation-simple.R`](garch-simulation-simple.R) | Runs a GARCH-type heavy-tailed simulation and compares the proposed test with Yin's test under coordinate-dependent missingness. |
| [`R code to compare power between our and Yin model-2.R`](R%20code%20to%20compare%20power%20between%20our%20and%20Yin%20model-2-2.R) | Applies both tests to the EEG example for a selected signal strength and reports their rejection indicators. |
| [`R code to check columns are not identical.R`](R%20code%20to%20check%20columns%20are%20not%20identical.R) | Performs pairwise two-sample Kolmogorov–Smirnov tests on selected blocks of the EEG data and applies the Benjamini–Hochberg adjustment. |

## Statistical setting

Let \(X_i=(X_{i1},\ldots,X_{ip})^\top\), \(i=1,\ldots,n\), be high-dimensional observations with potentially missing coordinates. The observed data are represented using indicators \(\varepsilon_{ik}\), where \(\varepsilon_{ik}=1\) means that coordinate \(k\) is observed for subject \(i\).

The scripts study the hypothesis

\[
H_0:\mu=\mu_0
\qquad\text{versus}\qquad
H_1:\mu\ne\mu_0,
\]

using a quadratic statistic. The proposed procedure employs truncation in the null-variance estimator to improve robustness to heavy-tailed observations. The simulation scripts select the truncation threshold from an independent null pilot sample; Yin's estimator remains untruncated.

The reported rejection probability at zero signal estimates empirical size, while rejection probabilities at nonzero signal values estimate empirical power.

## Requirements

- R version 4.0 or later is recommended.
- The three simulation scripts use base R only.
- The EEG comparison script loads the `readr` package, which can be installed with:

```r
install.packages("readr")
```

The simulation settings can be computationally demanding because the default number of Monte Carlo replications is `B = 1000` and the dimension can be proportional to the sample size.

## Running the simulations

Clone the repository and open its directory in R or RStudio:

```bash
git clone https://github.com/shivkky/R-code-for-article-High-Dimensional-Mean-Testing-under-Controlled-dense-Missing-Pattern.git
cd R-code-for-article-High-Dimensional-Mean-Testing-under-Controlled-dense-Missing-Pattern
```

Run a script from the R console, for example:

```r
source("lighttail_revised.R")
source("linear-example2-simulation-revised.R")
source("garch-simulation-simple.R")
```

Before running, edit the **VALUES TO CHANGE** section near the beginning of the selected script. The main parameters include:

- `n`: sample size;
- `p`: dimension;
- `B`: number of Monte Carlo replications;
- `test_level`: nominal significance level;
- `master_seed`: seed for reproducibility;
- `theta` or `b`: dependence parameter;
- `q_observed` or `beta_values`: missingness parameters; and
- `show_pilot_boxplots`: whether to display pilot-sample boxplots.

For a quick code check, first set `B` to a small value such as `10`. Use `B = 1000` (or a larger value) for the final Monte Carlo study.

### Simulation outputs

The scripts print empirical rejection probabilities to the R console. Depending on the script, they also generate:

- LaTeX tables for the proposed and Yin tests;
- a CSV file containing numerical results;
- an RDS file containing settings, the pilot threshold, signal direction, and results; and
- a PDF showing the pilot-sample threshold selection.

Output files are written to the current R working directory. Existing files with the same names may be overwritten.

## EEG-data scripts

The EEG data are not included in this repository. Before running the EEG scripts, load the required objects into the R environment.

### Power comparison

`R code to compare power between our and Yin model-2-2.R` expects:

- `Final_data_for_control_dence_paper`: an EEG data matrix with coordinates in rows and individuals in columns; and
- `EEG_control_group_null_mean_for_ALEA`: the null mean vector.

Then run:

```r
source("R code to compare power between our and Yin model-2-2.R")
```

The script calculates the standardized statistics and rejection indicators for the proposed and Yin procedures. The signal parameter `alpha` and tuning quantities can be changed near the beginning of the file.

### Pairwise distributional comparisons

`R code to check columns are not identical.R` expects:

- `Final_data_for_control_dence_paper`; and
- `Final_data_for_ALEA_paper`.

It first compares all column pairs using rows 1–16. Pairs for which no significant difference is detected are then examined using rows 17–32. P-values are adjusted using the Benjamini–Hochberg method.

```r
source("R code to check columns are not identical.R")
```

These pairwise KS comparisons describe differences between the selected empirical coordinate blocks; they should not be interpreted as a direct test that the individual-level multivariate observations are identically distributed.

## Reproducibility notes

- Run each script from the beginning after changing its settings.
- Keep the supplied random seeds to reproduce the corresponding experiment.
- Each pilot threshold is estimated independently of the Monte Carlo samples and then held fixed across replications for the relevant setting.
- Check the working directory with `getwd()` before running a script if you want outputs saved in a particular location.
- Numerical results can vary across R versions and computing platforms because of differences in random-number generation and numerical libraries.

## Citation

If you use this code, please cite the accompanying article:

> Shiv Kumar Yadav and Monika Bhattacharjee. *High-Dimensional Mean Testing under a Controlled Dense Missing Pattern*. Manuscript.

Please update this entry with the journal, year, volume, pages, and DOI when the final bibliographic information becomes available.

## Contact

For questions about the code, please open an issue in this GitHub repository or contact the repository owner through the [GitHub profile](https://github.com/shivkky).
