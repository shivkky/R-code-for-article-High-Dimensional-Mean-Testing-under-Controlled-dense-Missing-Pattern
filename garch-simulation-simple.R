# Simple GARCH-type simulation comparing the proposed test and Yin's test
#
# Run the complete script from the beginning after changing any value below.
# The same truncation threshold, selected from an independent null pilot
# sample, is used in every Monte Carlo replication of the proposed test.
# Yin's variance estimator is untruncated and therefore does not use tau.


# ===========================================================================
# VALUES TO CHANGE
# ===========================================================================

n <- 500                       # sample size
p <- 3*n/2                       # dimension
B <- 1000                      # number of Monte Carlo replications
b <- 0.5                       # parameter in the GARCH-type recursion
test_level <- 0.05             # level of the two-sided test
master_seed <- 20260827        # new seed for the complete experiment
show_pilot_boxplots <- TRUE    # change to FALSE to suppress the plots
proposed_latex_output_file <- "garch_proposed_power_table_n100.tex"
yin_latex_output_file <- "garch_yin_power_table_n100.tex"

# The nonzero rows estimate empirical power.  The row alpha = 0 estimates
# empirical size.  In every case, Delta[k] = alpha * n^(-1/8) * v[k].


# ===========================================================================
# FUNCTIONS USED TO GENERATE AND ANALYSE ONE DATASET
# ===========================================================================

# Generate a symmetric Pareto variable with
# P(|Z| > x) = (2/x)^3.5, x >= 2.
generate_signed_pareto <- function(number) {
  if (number == 0L) {
    return(numeric(0L))
  }

  magnitude <- 2 * runif(number)^(-1 / 3.5)
  random_sign <- sample(c(-1, 1), number, replace = TRUE)
  random_sign * magnitude
}


# STEP 1: Generate the complete GARCH-type data.
generate_complete_data <- function(
    n,
    p,
    alpha_value,
    v,
    b,
    seed) {
  set.seed(seed)

  # Mixture innovations Z have variance 10.  The small correction
  # (1 + i*k)^(-2) avoids the undefined boundary case at (i,k) = (1,1).
  index_product <- outer(seq_len(n), seq_len(p), `*`)
  pareto_probability <- 1 / (1 + index_product)^2

  pareto_second_moment <- 28 / 3
  gaussian_variance <-
    (10 - pareto_probability * pareto_second_moment) /
    (1 - pareto_probability)

  is_pareto <-
    matrix(runif(n * p), nrow = n, ncol = p) < pareto_probability

  Z <- matrix(rnorm(n * p), nrow = n, ncol = p) *
    sqrt(gaussian_variance)

  number_pareto <- sum(is_pareto)
  Z[is_pareto] <- generate_signed_pareto(number_pareto)

  # The simulation design assumes that 10*u has the distribution of Z.
  u <- Z / 10

  # GARCH-type recursion across the coordinate index k.
  W <- matrix(0, nrow = n, ncol = p)
  previous_W <- rep(0, n)

  for (k in seq_len(p)) {
    conditional_scale <- sqrt(1 - b^2) + b * abs(previous_W)
    W[, k] <- conditional_scale * u[, k]
    previous_W <- W[, k]
  }

  # Add the mean vector.  Under H0, alpha_value = 0.
  Delta <- alpha_value * n^(-1 / 8) * v
  X <- sweep(W, MARGIN = 2, STATS = Delta, FUN = "+")

  list(X = X, W = W, Delta = Delta)
}


# STEP 2: Generate the missingness indicators and mask the observations.
give_missingness <- function(complete_data, beta_missing, seed) {
  set.seed(seed)

  X <- complete_data$X
  n <- nrow(X)
  p <- ncol(X)
  k <- seq_len(p)

  missing_probability <-
    1 / (log(n) * n^0.02) + beta_missing / (1 + k)

  if (any(missing_probability < 0) ||
      any(missing_probability > 1)) {
    stop("The missingness probabilities are outside [0,1].")
  }

  observed_probability <- 1 - missing_probability
  epsilon <- matrix(
    rbinom(
      n * p,
      size = 1,
      prob = rep(observed_probability, each = n)
    ),
    nrow = n,
    ncol = p
  )

  X_observed <- X
  X_observed[epsilon == 0] <- NA_real_

  list(
    X_observed = X_observed,
    epsilon = epsilon,
    Delta = complete_data$Delta,
    observed_proportion = mean(epsilon)
  )
}


# STEPS 4 AND 5A: Calculate the proposed statistic and perform its test.
calculate_proposed_test <- function(data, tau, test_level = 0.05, mu0 = 0) {
  X_observed <- data$X_observed
  n <- nrow(X_observed)
  p <- ncol(X_observed)

  if (length(mu0) == 1L) {
    mu0 <- rep(mu0, p)
  }

  # Missing entries are replaced by zero after centering.
  centered <- sweep(X_observed, MARGIN = 2, STATS = mu0, FUN = "-")
  A <- centered
  A[!is.finite(A)] <- 0

  # U_n = {n(n-1)}^(-1) sum_{i != j} sum_k A[i,k] A[j,k].
  U_n <- sum(colSums(A)^2 - colSums(A^2)) / (n * (n - 1))

  # Truncation is applied only in the variance estimator.
  A_truncated <- A
  A_truncated[abs(A_truncated) > tau] <- 0

  G <- tcrossprod(A_truncated)
  diag(G) <- 0
  variance_estimate <-
    (2 / (n^2 * (n - 1)^2)) * sum(G^2)

  if (!is.finite(variance_estimate) || variance_estimate <= 0) {
    stop("The estimated null variance is not positive.")
  }

  test_statistic <- U_n / sqrt(variance_estimate)
  critical_value <- qnorm(1 - test_level / 2)
  reject <- abs(test_statistic) > critical_value
  p_value <- 2 * pnorm(-abs(test_statistic))

  c(
    U_n = U_n,
    variance_estimate = variance_estimate,
    test_statistic = test_statistic,
    p_value = p_value,
    reject = as.integer(reject)
  )
}


# STEPS 4 AND 5B: Calculate Yin's statistic and perform its test.
#
# Yin's test has the same quadratic numerator U_n.  Its variance estimator
# uses the observed centered residuals without truncation.
calculate_yin_test <- function(data, test_level = 0.05, mu0 = 0) {
  X_observed <- data$X_observed
  n <- nrow(X_observed)
  p <- ncol(X_observed)

  if (length(mu0) == 1L) {
    mu0 <- rep(mu0, p)
  }

  # Missing entries are replaced by zero after centering.
  centered <- sweep(X_observed, MARGIN = 2, STATS = mu0, FUN = "-")
  A <- centered
  A[!is.finite(A)] <- 0

  # The numerator is identical to that of the proposed test.
  U_n <- sum(colSums(A)^2 - colSums(A^2)) / (n * (n - 1))

  # Yin's pairwise squared-kernel variance estimator is untruncated.
  G <- tcrossprod(A)
  diag(G) <- 0
  variance_estimate <-
    (8 / (n^2 * (n - 1)^2)) * sum(G^2)

  if (!is.finite(variance_estimate) || variance_estimate <= 0) {
    stop("Yin's estimated null variance is not positive.")
  }

  test_statistic <- U_n / sqrt(variance_estimate)
  critical_value <- qnorm(1 - test_level / 2)
  reject <- abs(test_statistic) > critical_value
  p_value <- 2 * pnorm(-abs(test_statistic))

  c(
    U_n = U_n,
    variance_estimate = variance_estimate,
    test_statistic = test_statistic,
    p_value = p_value,
    reject = as.integer(reject)
  )
}


# ===========================================================================
# STEP 3: SET THE PARAMETER GRID AND GENERATE v ONCE
# ===========================================================================

alpha_values <- c(
  n^(-7 / 8) / 100,
  n^(-5 / 8) / 100,
  0,
  0.004,
  0.006,
  0.008,
  0.01,
  0.012
)

alpha_labels <- c(
  "$n^{-7/8}/100$",
  "$n^{-5/8}/100$",
  "0",
  "0.004",
  "0.006",
  "0.008",
  "0.01",
  "0.012"
)

beta_values <- c(0, 0.1, 0.2, 0.4)
beta_labels <- c("0", "0.1", "0.2", "0.4")

set.seed(master_seed)
v <- runif(p, min = 2, max = 7)

proposed_rejection_table <- matrix(
  NA_real_,
  nrow = length(alpha_values),
  ncol = length(beta_values),
  dimnames = list(alpha_labels, beta_labels)
)

yin_rejection_table <- matrix(
  NA_real_,
  nrow = length(alpha_values),
  ncol = length(beta_values),
  dimnames = list(alpha_labels, beta_labels)
)

thresholds <- numeric(length(beta_values))

if (show_pilot_boxplots) {
  old_graphical_parameters <- par(no.readonly = TRUE)
  par(mfrow = c(2, 2))
}


# ===========================================================================
# STEPS 3--6: PILOT THRESHOLD AND MONTE CARLO CALCULATION FOR EACH CELL
# ===========================================================================

for (beta_index in seq_along(beta_values)) {
  beta_current <- beta_values[beta_index]

  # Independent null pilot sample for the current value of beta.
  pilot_complete <- generate_complete_data(
    n = n,
    p = p,
    alpha_value = 0,
    v = v,
    b = b,
    seed = master_seed + 1000 * beta_index
  )

  pilot_data <- give_missingness(
    complete_data = pilot_complete,
    beta_missing = beta_current,
    seed = master_seed + 5000 + beta_index
  )

  pilot_absolute_values <- abs(as.vector(pilot_data$X_observed))
  pilot_absolute_values <-
    pilot_absolute_values[is.finite(pilot_absolute_values)]

  # The threshold is twice the upper whisker of the corresponding boxplot.
  pilot_boxplot <- boxplot.stats(pilot_absolute_values)
  tau <- 2 * pilot_boxplot$stats[5L]
  thresholds[beta_index] <- tau

  if (show_pilot_boxplots) {
    boxplot(
      pilot_absolute_values,
      horizontal = TRUE,
      xlab = expression(abs(X[ik] - mu[0*k])),
      main = paste("beta =", beta_current)
    )
    abline(v = tau, col = "red", lwd = 2, lty = 2)
  }

  cat(
    "beta =", beta_current,
    ": threshold =", tau,
    "; observed proportion =", pilot_data$observed_proportion,
    "\n"
  )

  for (alpha_index in seq_along(alpha_values)) {
    alpha_current <- alpha_values[alpha_index]
    proposed_rejections <- numeric(B)
    yin_rejections <- numeric(B)

    for (replication in seq_len(B)) {
      complete_data <- generate_complete_data(
        n = n,
        p = p,
        alpha_value = alpha_current,
        v = v,
        b = b,
        seed = master_seed +
          100000 * beta_index + 10000 * alpha_index + replication
      )

      observed_data <- give_missingness(
        complete_data = complete_data,
        beta_missing = beta_current,
        seed = master_seed + 1000000 +
          100000 * beta_index + 10000 * alpha_index + replication
      )

      proposed_result <- calculate_proposed_test(
        data = observed_data,
        tau = tau,
        test_level = test_level,
        mu0 = 0
      )

      yin_result <- calculate_yin_test(
        data = observed_data,
        test_level = test_level,
        mu0 = 0
      )

      proposed_rejections[replication] <- proposed_result["reject"]
      yin_rejections[replication] <- yin_result["reject"]
    }

    proposed_rejection_table[alpha_index, beta_index] <-
      mean(proposed_rejections)

    yin_rejection_table[alpha_index, beta_index] <-
      mean(yin_rejections)

    cat(
      "  alpha =", alpha_labels[alpha_index],
      ": proposed =",
      proposed_rejection_table[alpha_index, beta_index],
      "; Yin =",
      yin_rejection_table[alpha_index, beta_index],
      "\n"
    )
  }
}

if (show_pilot_boxplots) {
  par(old_graphical_parameters)
}


# ===========================================================================
# NUMERICAL AND LATEX OUTPUTS
# ===========================================================================

proposed_numerical_table <- data.frame(
  alpha_prime = alpha_labels,
  proposed_rejection_table,
  check.names = FALSE
)

yin_numerical_table <- data.frame(
  alpha_prime = alpha_labels,
  yin_rejection_table,
  check.names = FALSE
)

cat("\nProposed-test Monte Carlo rejection probabilities:\n")
print(proposed_numerical_table, row.names = FALSE)

cat("\nYin-test Monte Carlo rejection probabilities:\n")
print(yin_numerical_table, row.names = FALSE)

make_latex_table <- function(rejection_table, method_name) {
  latex_rows <- vapply(
    seq_along(alpha_labels),
    function(row_index) {
      paste0(
        alpha_labels[row_index],
        " & ",
        paste(
          sprintf("%.3f", rejection_table[row_index, ]),
          collapse = " & "
        ),
        " \\\\"
      )
    },
    character(1L)
  )

  c(
    "\\begin{minipage}{0.32\\textwidth}\\centering",
    paste0("\\textbf{", method_name, "}\\\\[2pt]"),
    "\\begin{tabular}{ccccc}",
    "\\toprule",
    paste0("\\multicolumn{5}{c}{$n=", n, "$}\\\\"),
    "\\midrule",
    "$\\alpha^\\prime$ & \\multicolumn{4}{c}{$\\beta$}\\\\",
    "\\cmidrule(lr){2-5}",
    paste0("& ", paste(beta_labels, collapse = " & "), "\\\\"),
    "\\midrule",
    latex_rows,
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{minipage}\\hfill"
  )
}

proposed_latex_table <- make_latex_table(
  proposed_rejection_table,
  "Proposed test"
)

yin_latex_table <- make_latex_table(
  yin_rejection_table,
  "Yin's test"
)

cat("\nProposed-test LaTeX table:\n")
cat(paste(proposed_latex_table, collapse = "\n"), "\n")
writeLines(proposed_latex_table, con = proposed_latex_output_file)

cat("\nYin-test LaTeX table:\n")
cat(paste(yin_latex_table, collapse = "\n"), "\n")
writeLines(yin_latex_table, con = yin_latex_output_file)

cat(
  "\nLaTeX tables written to:",
  proposed_latex_output_file,
  "and",
  yin_latex_output_file,
  "\n"
)
