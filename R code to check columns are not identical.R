# Example matrix
set.seed(123)

X = as.data.frame(Final_data_for_control_dence_paper)

X = X[1:16, , drop = FALSE]
# All pairs of columns
column_pairs = combn(ncol(X), 2)

# Store results
results = data.frame()

for (i in seq_len(ncol(column_pairs))) {
  
  j = column_pairs[1, i]
  k = column_pairs[2, i]
  
  x = X[, j]
  y = X[, k]
  
  # KS test
  test_result = ks.test(x, y, exact = FALSE)
  
  results = rbind(
    results,
    data.frame(
      Column_1   = colnames(X)[j],
      Column_2   = colnames(X)[k],
      KS_distance = unname(test_result$statistic),
      P_value     = test_result$p.value
    )
  )
}

# Adjust p-values because several tests are performed
results$Adjusted_P_value = p.adjust(
  results$P_value,
  method = "BH"
)

results$Conclusion = ifelse(
  results$Adjusted_P_value < 0.05,
  "Different distributions",
  "No significant difference detected"
)

print(results)

s = sum(results$Conclusion == "Different distributions", na.rm = TRUE)

diff_dist = s/length(results$Column_2)
diff_dist


### 
# Use the full dataset
X_full = as.data.frame(Final_data_for_ALEA_paper)

stopifnot(nrow(X_full) >= 32)

# Select pairs with no significant difference in rows 1–16
selected_pairs = results[
  !is.na(results$Conclusion) &
    results$Conclusion == "No significant difference detected",
  c("Column_1", "Column_2"),
  drop = FALSE
]

# Store results for rows 17–32
results_17_32 = data.frame(
  Column_1 = character(),
  Column_2 = character(),
  KS_distance = numeric(),
  P_value = numeric()
)

for (i in seq_len(nrow(selected_pairs))) {
  
  j = selected_pairs$Column_1[i]
  k = selected_pairs$Column_2[i]
  
  # Extract rows 17–32 of the selected columns
  x = X_full[17:32, j]
  y = X_full[17:32, k]
  
  x = x[!is.na(x)]
  y = y[!is.na(y)]
  
  test_result = ks.test(x, y, exact = FALSE)
  
  results_17_32 = rbind(
    results_17_32,
    data.frame(
      Column_1 = j,
      Column_2 = k,
      KS_distance = unname(test_result$statistic),
      P_value = test_result$p.value
    )
  )
}

# Adjust across the selected comparisons
results_17_32$Adjusted_P_value = p.adjust(
  results_17_32$P_value,
  method = "BH"
)

results_17_32$Indicator = as.integer(
  results_17_32$Adjusted_P_value < 0.05
)

results_17_32$Conclusion = ifelse(
  results_17_32$Indicator == 1,
  "Different distributions",
  "No significant difference detected"
)

print(results_17_32)

# Proportion of selected pairs now showing a significant difference
s = sum(results_17_32$Indicator, na.rm = TRUE)

diff_dist1 = s/2415
diff_dist1


