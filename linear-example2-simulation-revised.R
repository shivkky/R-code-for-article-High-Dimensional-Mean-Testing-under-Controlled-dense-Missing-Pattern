# Linear-process simulation with the innovations from Example 2
#
# The script compares the proposed test with Yin's test.  It follows these
# six steps:
#   1. Generate complete linear-process observations.
#   2. Generate missingness indicators and mask the observations.
#   3. Select one threshold for each beta from an independent null pilot sample.
#   4. Calculate the two test statistics.
#   5. Apply two-sided tests at the chosen significance level.
#   6. Repeat the experiment B times to estimate rejection probabilities.
#
# The same pilot-selected threshold is used for every Monte Carlo replication
# having the same value of beta.  Yin's variance estimator is untruncated and
# therefore does not use this threshold.


# ===========================================================================
# VALUES TO CHANGE
# ===========================================================================

n <- 200                  # sample size
p <- 3*n/2 #nteger(3 * n / 2)        # dimension; also try n or n / 2
B <- 1000                        # number of Monte Carlo replications
theta <- 0.5                      # linear-process coefficient, 0 <= theta < 1
test_level <- 0.05                # level of the two-sided test
master_seed <- 20260910         # seed for the complete experiment
show_pilot_boxplots <- TRUE       # FALSE suppresses the pilot plots

proposed_latex_output_file <- "linear_proposed_rejection_table.tex"
yin_latex_output_file <- "linear_yin_rejection_table.tex"

# The nonzero alpha values estimate empirical power.  The row alpha = 0
# estimates empirical size.  In every case,
#
#     Delta[k] = alpha * n^(-1/8) * v[k].


# ===========================================================================
# FUNCTIONS USED TO GENERATE ONE DATASET
# ===========================================================================

# Generate the symmetric Pareto component in Example 2:
#
#     P(|Y| > x) = (2/x)^3.5,  x >= 2.
#
# Its density is 14*sqrt(2)*|x|^(-4.5) on {|x| > 2}, its variance is
# 28/3, its third moment is finite, and its fourth moment is infinite.
generate_signed_pareto <- function(number) {
  if (number == 0L) {
    return(numeric(0L))
  }

  magnitude <- 2 * runif(number)^(-1 / 3.5)
  random_sign <- sample(c(-1, 1), number, replace = TRUE)
  random_sign * magnitude
}


# Generate the independent innovations used in Example 2.
#
# The printed probability (i*k)^(-2) is one at (i,k)=(1,1), at which point
# its displayed Gaussian variance is undefined.  As in the accompanying
# GARCH script, the code uses the boundary-corrected probability
#
#     pi[i,k] = (1 + i*k)^(-2).
#
# This gives every innovation mean zero and variance 10.
generate_example2_innovations <- function(n, p) {
  row_index <- seq_len(n)
  coordinate_index <- seq_len(p)
  index_product <- outer(row_index, coordinate_index, `*`)
  pareto_probability <- 1 / (1 + index_product)^2

  pareto_second_moment <- 28 / 3
  gaussian_variance <-
    (10 - pareto_probability * pareto_second_moment) /
    (1 - pareto_probability)

  is_pareto <-
    matrix(runif(n * p), nrow = n, ncol = p) < pareto_probability

  innovations <- matrix(rnorm(n * p), nrow = n, ncol = p) *
    sqrt(gaussian_variance)

  number_pareto <- sum(is_pareto)
  innovations[is_pareto] <- generate_signed_pareto(number_pareto)

  list(
    innovations = innovations,
    is_pareto = is_pareto,
    number_pareto = number_pareto
  )
}


# STEP 1: Generate the complete linear-process data.
#
# The manuscript writes
#
#     W[i,k] = zeta[i,k] + sum_{j>=1} theta^j zeta[i,k-j].
#
# Because Example 2 specifies the heterogeneous innovations only for k >= 1,
# the simulation initializes W[i,0] at zero and uses the equivalent recursion
#
#     W[i,k] = zeta[i,k] + theta * W[i,k-1].
#
# The omitted infinite-past boundary term decays geometrically as theta^k.
generate_complete_data <- function(
    n,
    p,
    alpha_value,
    v,
    theta,
    seed) {
  if (!is.finite(theta) || theta < 0 || theta >= 1) {
    stop("theta must satisfy 0 <= theta < 1.")
  }

  set.seed(seed)
  innovation_result <- generate_example2_innovations(n = n, p = p)
  zeta <- innovation_result$innovations

  W <- matrix(0, nrow = n, ncol = p)
  previous_W <- numeric(n)

  for (k in seq_len(p)) {
    W[, k] <- zeta[, k] + theta * previous_W
    previous_W <- W[, k]
  }

  # Add the mean vector.  Under H0, alpha_value = 0.
  Delta <- alpha_value * n^(-1 / 8) * v
  X <- sweep(W, MARGIN = 2, STATS = Delta, FUN = "+")

  list(
    X = X,
    W = W,
    zeta = zeta,
    Delta = Delta,
    number_pareto = innovation_result$number_pareto
  )
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


# ===========================================================================
# FUNCTIONS USED TO CALCULATE THE TESTS
# ===========================================================================

# STEPS 4 AND 5A: Calculate the proposed statistic and perform its test.
calculate_proposed_test <- function(data, tau, test_level = 0.05, mu0 = 0) {
  X_observed <- data$X_observed
  n <- nrow(X_observed)
  p <- ncol(X_observed)

  if (length(mu0) == 1L) {
    mu0 <- rep(mu0, p)
  }
  if (length(mu0) != p) {
    stop("mu0 must be a scalar or a vector of length p.")
  }

  # Missing entries are replaced by zero after centering.  Thus A[i,k]
  # equals epsilon[i,k] * (X[i,k] - mu0[k]).
  centered <- sweep(X_observed, MARGIN = 2, STATS = mu0, FUN = "-")
  A <- centered
  A[!is.finite(A)] <- 0

  # U_n = {n(n-1)}^(-1) sum_{i != j} sum_k A[i,k] A[j,k].
  U_n <- sum(colSums(A)^2 - colSums(A^2)) / (n * (n - 1))

  # The proposed variance estimator uses hard deletion above tau.
  # The numerator U_n itself remains untruncated, as in the manuscript.
  A_truncated <- A
  A_truncated[abs(A_truncated) > tau] <- 0

  G <- tcrossprod(A_truncated)
  diag(G) <- 0
  variance_estimate <-
    (2 / (n^2 * (n - 1)^2)) * sum(G^2)

  if (!is.finite(variance_estimate) || variance_estimate <= 0) {
    stop("The proposed estimated null variance is not positive.")
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
# Yin's test uses the same untruncated quadratic numerator.  Its variance
# estimator uses the observed centered values without truncation.
calculate_yin_test <- function(data, test_level = 0.05, mu0 = 0) {
  X_observed <- data$X_observed
  n <- nrow(X_observed)
  p <- ncol(X_observed)

  if (length(mu0) == 1L) {
    mu0 <- rep(mu0, p)
  }
  if (length(mu0) != p) {
    stop("mu0 must be a scalar or a vector of length p.")
  }

  centered <- sweep(X_observed, MARGIN = 2, STATS = mu0, FUN = "-")
  A <- centered
  A[!is.finite(A)] <- 0

  U_n <- sum(colSums(A)^2 - colSums(A^2)) / (n * (n - 1))

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
  0.04,
  0.06,
  0.08,
  0.10,
  0.12
)

alpha_labels <- c(
  "{\\fontsize{4pt}{5pt}\\selectfont $n^{-7/8}/100$}",
  "{\\fontsize{4pt}{5pt}\\selectfont $n^{-5/8}/100$}",
  "0",
  "0.04",
  "0.06",
  "0.08",
  "0.1",
  "0.12"
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
    theta = theta,
    seed = master_seed + 1000L * beta_index
  )

  pilot_data <- give_missingness(
    complete_data = pilot_complete,
    beta_missing = beta_current,
    seed = master_seed + 5000L + beta_index
  )

  pilot_absolute_values <- abs(as.vector(pilot_data$X_observed))
  pilot_absolute_values <-
    pilot_absolute_values[is.finite(pilot_absolute_values)]

  # The threshold is 1.5 times the upper whisker of the corresponding boxplot.
  pilot_boxplot <- boxplot.stats(pilot_absolute_values)
  tau <-  (3/2)*pilot_boxplot$stats[5L]
  thresholds[beta_index] <- tau

  if (show_pilot_boxplots) {
    boxplot(
      pilot_absolute_values,
      horizontal = TRUE,
      xlab = expression(abs(X[ik] - mu[0 * k])),
      main = paste("beta =", beta_current)
    )
    abline(v = tau, col = "red", lwd = 2, lty = 2)
  }

  cat(
    "beta =", beta_current,
    ": threshold =", tau,
    "; observed proportion =", pilot_data$observed_proportion,
    "; pilot Pareto entries =", pilot_complete$number_pareto,
    "\n"
  )

  for (alpha_index in seq_along(alpha_values)) {
    alpha_current <- alpha_values[alpha_index]
    proposed_rejections <- numeric(B)
    yin_rejections <- numeric(B)

    for (replication in seq_len(B)) {
      data_seed <- master_seed +
        100000L * beta_index + 10000L * alpha_index + replication

      missingness_seed <- master_seed + 1000000L +
        100000L * beta_index + 10000L * alpha_index + replication

      complete_data <- generate_complete_data(
        n = n,
        p = p,
        alpha_value = alpha_current,
        v = v,
        theta = theta,
        seed = data_seed
      )

      observed_data <- give_missingness(
        complete_data = complete_data,
        beta_missing = beta_current,
        seed = missingness_seed
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

threshold_table <- data.frame(
  beta = beta_labels,
  threshold = thresholds,
  check.names = FALSE
)

cat("\nPilot truncation thresholds:\n")
print(threshold_table, row.names = FALSE)

cat("\nProposed-test Monte Carlo rejection probabilities:\n")
print(proposed_numerical_table, row.names = FALSE)

cat("\nYin-test Monte Carlo rejection probabilities:\n")
print(yin_numerical_table, row.names = FALSE)


make_latex_table <- function(rejection_table, method_name) {
  threshold_row <- paste0(
    "$\\mathrm{tth}$ & ",
    paste(sprintf("%.3f", thresholds), collapse = " & "),
    " \\\\"
  )

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
    paste0("\\multicolumn{5}{c}{$n=", n, ",\\ p=", p, "$}\\\\"),
    "\\midrule",
    "$\\alpha^\\prime$ & \\multicolumn{4}{c}{$\\beta$}\\\\",
    "\\cmidrule(lr){2-5}",
    paste0("& ", paste(beta_labels, collapse = " & "), "\\\\"),
    "\\midrule",
    threshold_row,
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

