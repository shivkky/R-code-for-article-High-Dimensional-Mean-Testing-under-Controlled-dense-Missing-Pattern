rm(list = ls())
library(readr)
#set.seed(123)
data_matrix = data.matrix(Final_data_for_control_dence_paper)
n= ncol(data_matrix)    # Sample size
p= nrow(data_matrix)    # Dimensions

W_val= t(data_matrix)

eta_n = 1/(n^(7/24))
Upsilon_n =  sqrt(2.437725)
delta = 1
alpha = 0.309

data_mu_0 = data.matrix(EEG_control_group_null_mean_for_ALEA)
mu_0 = data_mu_0
v = c(runif(p,0,1))

#control_group_mean_data = read.csv("EEG_data_to_calculate_control_group_mean.csv")
#col_mean_nonzero = apply(control_group_mean_data, 1, function(x) {
#  nz = x[x != 0]           # keep only non-zero elements
#  if(length(nz) == 0) NA else mean(nz)  # handle case with all zeros
#})

#col_mean_nonzero
mu = alpha*mu_0

M = (W_val!=0)*1

Mu_matrix = matrix(mu_0, nrow = n, ncol = p, byrow = TRUE)
Z2 = (W_val - Mu_matrix)*M
Z1 = W_val + matrix(mu, nrow = n, ncol = p, byrow = TRUE)


T1 = T2 = T3 = 0
for (i in 1:n) {
  for (j in 1:n) {
    if (i != j) {
      T1 = T1 + sum(Z1[i, ] * Z1[j, ])
      T2 = T2 + sum(Z1[i, ] * (M[j, ] * mu_0))
      T3 = T3 + sum(mu_0 * (M[i, ] * M[j, ] * mu_0))
    }
  }
}
T1 = T1 / (n * (n - 1))
T2 = T2 / (n * (n - 1))
T3 = T3 / (n * (n - 1))
T_n = T1 - 2 * T2 + T3
T_n

W = Z2
threshold =  eta_n * n^(1 / (2 + delta)) * Upsilon_n
Wf = W * (abs(W) < threshold)
Wf_masked = Wf
G = Wf_masked %*% t(Wf_masked)
diag(G) = 0  # Set diagonal to 0 to exclude i == j
sigma2_hat = 2*sum(G^2) / (n^4)
sigma2_hat


T_q = qnorm(0.975, 0,1, lower.tail = TRUE)
T_n_hat = T_n / sqrt(sigma2_hat)
T_n_hat
indicator = as.numeric(abs(T_n_hat) >= T_q)
print(indicator)


## R code to calculate power in Yin Model

## Estimator to estimate sigma_n,0

W_Yin = Z2
G_Yin = W_Yin%*%t(W_Yin)
diag(G_Yin)= 0
sigma2_hat_Yin = 2*sum(G_Yin^2) / ((n*(n-1))^2)
sigma2_hat_Yin

## Test Statistics and decision

T_n_hat_Yin = T_n / sqrt(sigma2_hat_Yin)
T_n_hat_Yin
indicator_Yin = as.numeric(abs(T_n_hat_Yin) >= T_q)
indicator_Yin

