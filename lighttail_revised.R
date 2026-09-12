# ============================================================================
# LIGHT-TAILED LINEAR-PROCESS SIMULATION: ONE (n,p) PANEL
#
# Change only the values in the first section. The code:
#   1. selects the truncation threshold from an independent pilot sample;
#   2. calculates our test and Yin's test;
#   3. reports empirical rejection probabilities;
#   4. prints a copy-ready LaTeX panel in the R console.
# ============================================================================


# ============================================================================
# VALUES TO CHANGE
# ============================================================================

n <- 200
p <- n
# Other choices:
# p <- n / 2
# p <- 3 * n / 2
p <- as.integer(p)

B <- 1000L                 # number of Monte Carlo replications
theta <- 0.5
q_observed <- 0.5
null_mean <- 0
test_level <- 0.05
threshold_multiplier <- 2
master_seed <- 20261001L

output_prefix <- paste0("light_tailed_n", n, "_p", p)


# ============================================================================
# SIGNAL VALUES
# ============================================================================

alpha_values_for_n <- function(n) {
  c(
    n^(-7 / 8) / 10,
    n^(-5 / 8) / 10,
    0,
    0.02,
    0.04,
    0.06,
    0.08,
    0.10
  )
}

alpha_labels <- c(
  "$n^{-7/8}/10$",
  "$n^{-5/8}/10$",
  "0",
  "0.02",
  "0.04",
  "0.06",
  "0.08",
  "0.10"
)


# ============================================================================
# STEP 1: GENERATE THE LINEAR PROCESS
# ============================================================================

generate_linear_process <- function(n, p, theta) {
  if (!is.finite(theta) || abs(theta) >= 1) {
    stop("theta must satisfy |theta| < 1.")
  }

  zeta <- matrix(
    rnorm(n * p),
    nrow = n,
    ncol = p
  )

  W <- matrix(
    0,
    nrow = n,
    ncol = p
  )

  previous_W <- numeric(n)

  for (k in seq_len(p)) {
    W[, k] <- zeta[, k] + theta * previous_W
    previous_W <- W[, k]
  }

  W
}


# ============================================================================
# STEP 2: GENERATE THE MISSINGNESS INDICATORS
# ============================================================================

generate_missing_indicators <- function(n, p, q_observed) {
  if (!is.finite(q_observed) ||
      q_observed <= 0 ||
      q_observed > 1) {
    stop("q_observed must lie in (0,1].")
  }

  matrix(
    rbinom(
      n * p,
      size = 1,
      prob = q_observed
    ),
    nrow = n,
    ncol = p
  )
}


expand_null_mean <- function(null_mean, p) {
  if (length(null_mean) == 1L) {
    return(rep(null_mean, p))
  }

  if (length(null_mean) != p) {
    stop("null_mean must be a scalar or a vector of length p.")
  }

  null_mean
}


# ============================================================================
# STEP 3: SELECT THE TRUNCATION THRESHOLD FROM A PILOT SAMPLE
# ============================================================================

select_pilot_threshold <- function(settings, make_plot = TRUE) {
  set.seed(settings$seed + 1L)

  W <- generate_linear_process(
    n = settings$n,
    p = settings$p,
    theta = settings$theta
  )

  epsilon <- generate_missing_indicators(
    n = settings$n,
    p = settings$p,
    q_observed = settings$q_observed
  )

  mu0 <- expand_null_mean(
    settings$null_mean,
    settings$p
  )

  X <- sweep(
    W,
    MARGIN = 2,
    STATS = mu0,
    FUN = "+"
  )

  centered <- sweep(
    X,
    MARGIN = 2,
    STATS = mu0,
    FUN = "-"
  )

  absolute_observed_values <- abs(
    centered[epsilon == 1]
  )

  upper_whisker <- boxplot.stats(
    absolute_observed_values
  )$stats[5L]

  tau <- settings$threshold_multiplier * upper_whisker

  if (make_plot) {
    boxplot(
      absolute_observed_values,
      horizontal = TRUE,
      outline = TRUE,
      xlab = expression(abs(X[ik] - mu[0 * k])),
      main = paste0(
        "Pilot sample: n = ",
        settings$n,
        ", p = ",
        settings$p,
        ", tth = ",
        sprintf("%.3f", tau)
      )
    )

    abline(
      v = tau,
      col = "red",
      lty = 2,
      lwd = 2
    )
  }

  list(
    threshold = tau,
    upper_whisker = upper_whisker,
    pilot_observed_proportion = mean(epsilon)
  )
}


# ============================================================================
# STEP 4: CALCULATE THE TEST STATISTIC
# ============================================================================

# For a matrix A with
#
# A[i,k] = epsilon[i,k] * (X[i,k] - mu0[k]),
#
# this function calculates
#
# U_n = {n(n-1)}^(-1)
#       sum_{i != j} sum_k A[i,k]A[j,k],
#
# together with its null-variance estimator.

calculate_quadratic_test <- function(A, test_level) {
  sample_size <- nrow(A)

  column_sums <- colSums(A)
  column_square_sums <- colSums(A^2)

  U_n <- sum(
    column_sums^2 - column_square_sums
  ) / (sample_size * (sample_size - 1))

  inner_products <- tcrossprod(A)
  diag(inner_products) <- 0

  variance_estimate <- 2 * sum(inner_products^2) /
    (sample_size^2 * (sample_size - 1)^2)

  if (!is.finite(variance_estimate) ||
      variance_estimate <= 0) {
    return(
      c(
        U_n = U_n,
        variance_estimate = NA_real_,
        test_statistic = NA_real_,
        reject = NA_real_
      )
    )
  }

  test_statistic <- U_n / sqrt(variance_estimate)

  critical_value <- qnorm(
    1 - test_level / 2
  )

  c(
    U_n = U_n,
    variance_estimate = variance_estimate,
    test_statistic = test_statistic,
    reject = as.integer(
      abs(test_statistic) > critical_value
    )
  )
}


# Our proposed test truncates the data in both U_n and the
# variance estimator.

calculate_proposed_test <- function(A, tau, test_level) {
  A_truncated <- A

  A_truncated[
    abs(A_truncated) > tau
  ] <- 0

  calculate_quadratic_test(
    A = A_truncated,
    test_level = test_level
  )
}


# Yin's test uses the untruncated data.

calculate_yin_test <- function(A, test_level) {
  calculate_quadratic_test(
    A = A,
    test_level = test_level
  )
}


# ============================================================================
# LATEX OUTPUT
# ============================================================================

format_probability <- function(x) {
  if (!is.finite(x)) {
    return("NA")
  }

  x <- round(x, 3)

  if (isTRUE(all.equal(x, 1))) {
    return("1")
  }

  sprintf("%.3f", x)
}


make_latex_panel <- function(results) {
  required_columns <- c(
    "n",
    "threshold",
    "alpha_label",
    "proposed_rejection_probability",
    "yin_rejection_probability"
  )

  missing_columns <- setdiff(
    required_columns,
    names(results)
  )

  if (length(missing_columns) > 0L) {
    stop(
      "The following result columns are missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  latex_rows <- vapply(
    seq_len(nrow(results)),
    function(j) {
      paste0(
        results$alpha_label[j],
        " & ",
        format_probability(
          results$proposed_rejection_probability[j]
        ),
        " & ",
        format_probability(
          results$yin_rejection_probability[j]
        ),
        "\\\\"
      )
    },
    character(1L)
  )

  c(
    paste0(
      "\\LightPanel{",
      results$n[1L],
      "}{",
      sprintf("%.3f", results$threshold[1L]),
      "}{%"
    ),
    latex_rows,
    "}"
  )
}


# ============================================================================
# CHECK THE INPUT VALUES
# ============================================================================

validate_settings <- function(settings) {
  required_names <- c(
    "n",
    "p",
    "B",
    "theta",
    "q_observed",
    "null_mean",
    "test_level",
    "threshold_multiplier",
    "seed",
    "output_prefix"
  )

  missing_names <- setdiff(
    required_names,
    names(settings)
  )

  if (length(missing_names) > 0L) {
    stop(
      "Missing settings: ",
      paste(missing_names, collapse = ", "),
      "."
    )
  }

  if (settings$n < 2L) {
    stop("n must be at least 2.")
  }

  if (settings$p < 1L) {
    stop("p must be positive.")
  }

  if (settings$B < 1L) {
    stop("B must be positive.")
  }

  if (settings$test_level <= 0 ||
      settings$test_level >= 1) {
    stop("test_level must lie in (0,1).")
  }

  invisible(TRUE)
}


# ============================================================================
# STEPS 5 AND 6: TWO-SIDED TEST AND MONTE CARLO REPLICATIONS
# ============================================================================

run_light_panel <- function(settings) {
  validate_settings(settings)

  csv_file <- paste0(
    settings$output_prefix,
    "_results.csv"
  )

  rds_file <- paste0(
    settings$output_prefix,
    "_results.rds"
  )

  tex_file <- paste0(
    settings$output_prefix,
    "_panel.tex"
  )

  plot_file <- paste0(
    settings$output_prefix,
    "_pilot_boxplot.pdf"
  )

  pdf(
    plot_file,
    width = 7,
    height = 4
  )

  pilot <- tryCatch(
    select_pilot_threshold(
      settings,
      make_plot = TRUE
    ),
    finally = dev.off()
  )

  tau <- pilot$threshold

  alpha_values <- alpha_values_for_n(
    settings$n
  )

  mu0 <- expand_null_mean(
    settings$null_mean,
    settings$p
  )

  proposed_rejections <- integer(
    length(alpha_values)
  )

  yin_rejections <- integer(
    length(alpha_values)
  )

  proposed_valid_replications <- integer(
    length(alpha_values)
  )

  yin_valid_replications <- integer(
    length(alpha_values)
  )

  cat(
    "Running n =",
    settings$n,
    ", p =",
    settings$p,
    ", B =",
    settings$B,
    ", tth =",
    sprintf("%.3f", tau),
    "\n"
  )

  # Generate the signal direction once and keep it fixed
  # throughout this panel.

  set.seed(settings$seed + 2L)

  v <- runif(
    settings$p,
    min = 2,
    max = 7
  )

  for (replication in seq_len(settings$B)) {
    set.seed(
      settings$seed + 1000L + replication
    )

    W <- generate_linear_process(
      n = settings$n,
      p = settings$p,
      theta = settings$theta
    )

    epsilon <- generate_missing_indicators(
      n = settings$n,
      p = settings$p,
      q_observed = settings$q_observed
    )

    # The same W and epsilon are used for all signal levels
    # within a given Monte Carlo replication.

    for (j in seq_along(alpha_values)) {
      alpha_prime <- alpha_values[j]

      Delta <- alpha_prime *
        settings$n^(-1 / 8) *
        v

      X <- sweep(
        W,
        MARGIN = 2,
        STATS = mu0 + Delta,
        FUN = "+"
      )

      centered <- sweep(
        X,
        MARGIN = 2,
        STATS = mu0,
        FUN = "-"
      )

      A <- centered * epsilon

      proposed_result <- calculate_proposed_test(
        A = A,
        tau = tau,
        test_level = settings$test_level
      )

      yin_result <- calculate_yin_test(
        A = A,
        test_level = settings$test_level
      )

      if (!is.na(proposed_result["reject"])) {
        proposed_rejections[j] <-
          proposed_rejections[j] +
          as.integer(proposed_result["reject"])

        proposed_valid_replications[j] <-
          proposed_valid_replications[j] + 1L
      }

      if (!is.na(yin_result["reject"])) {
        yin_rejections[j] <-
          yin_rejections[j] +
          as.integer(yin_result["reject"])

        yin_valid_replications[j] <-
          yin_valid_replications[j] + 1L
      }
    }

    progress_interval <- max(
      1L,
      settings$B %/% 10L
    )

    if (replication %% progress_interval == 0L ||
        replication == settings$B) {
      cat(
        "Completed",
        replication,
        "of",
        settings$B,
        "replications\n"
      )
    }
  }

  proposed_probabilities <- ifelse(
    proposed_valid_replications > 0,
    proposed_rejections /
      proposed_valid_replications,
    NA_real_
  )

  yin_probabilities <- ifelse(
    yin_valid_replications > 0,
    yin_rejections /
      yin_valid_replications,
    NA_real_
  )

  results <- data.frame(
    n = rep(settings$n, length(alpha_values)),
    p = rep(settings$p, length(alpha_values)),
    threshold = rep(tau, length(alpha_values)),
    upper_whisker = rep(
      pilot$upper_whisker,
      length(alpha_values)
    ),
    pilot_observed_proportion = rep(
      pilot$pilot_observed_proportion,
      length(alpha_values)
    ),
    alpha_label = alpha_labels,
    alpha_prime = alpha_values,
    proposed_rejection_probability =
      proposed_probabilities,
    yin_rejection_probability =
      yin_probabilities,
    proposed_valid_replications =
      proposed_valid_replications,
    yin_valid_replications =
      yin_valid_replications,
    stringsAsFactors = FALSE
  )

  write.csv(
    results,
    csv_file,
    row.names = FALSE
  )

  saveRDS(
    list(
      settings = settings,
      pilot = pilot,
      signal_direction = v,
      results = results
    ),
    rds_file
  )

  latex_output <- make_latex_panel(results)

  writeLines(
    latex_output,
    con = tex_file
  )

  cat("\nNumerical results:\n")

  print(
    data.frame(
      alpha = results$alpha_label,
      proposed =
        results$proposed_rejection_probability,
      yin =
        results$yin_rejection_probability,
      check.names = FALSE
    ),
    row.names = FALSE
  )

  cat(
    "\nPilot upper whisker:",
    sprintf("%.3f", pilot$upper_whisker),
    "\nTruncation threshold:",
    sprintf("%.3f", tau),
    "\nPilot observed proportion:",
    sprintf("%.3f", pilot$pilot_observed_proportion),
    "\n"
  )

  cat(
    "\n============================================================\n",
    "COPY THE FOLLOWING OUTPUT DIRECTLY INTO THE LATEX TABLE\n",
    "============================================================\n\n",
    sep = ""
  )

  writeLines(latex_output)

  cat(
    "\n============================================================\n"
  )

  invisible(results)
}


# ============================================================================
# RUN THE SELECTED PANEL
# ============================================================================

settings <- list(
  n = n,
  p = p,
  B = B,
  theta = theta,
  q_observed = q_observed,
  null_mean = null_mean,
  test_level = test_level,
  threshold_multiplier = threshold_multiplier,
  seed = master_seed,
  output_prefix = output_prefix
)

results <- run_light_panel(settings)