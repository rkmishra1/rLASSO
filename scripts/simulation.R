# Simulation Study for Reciprocal Bayesian LASSO and Competing Methods
# Scenarios evaluated:
# 1. Scenario 1 (Case I): n=50, p=20, Isotropic (rho=0.0), sigma=3.0, beta0 = (5, 0, ..., 0)
# 2. Scenario 2 (Case III): n=100, p=50, AR (rho=0.95), sigma=1.5, beta0 = (5, 0, ..., 0)
# 3. Scenario 3 (Case VIII): n=50, p=20, CS (rho=0.5), sigma=3.0, beta0 = (3, 1.5, 0, 0, 2, 0, ..., 0)

source("R/rLASSO.R")
library(glmnet)

# Function to generate design matrix X
generate_X <- function(n, p, rho, design) {
  if (design == "IS") {
    Sigma <- diag(1, p)
  } else if (design == "CS") {
    Sigma <- matrix(rho, p, p)
    diag(Sigma) <- 1
  } else if (design == "AR") {
    Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
  }
  # Chol decomposition
  ch <- chol(Sigma)
  X <- matrix(rnorm(n * p), n, p) %*% ch
  return(X)
}

# Function to calculate Balanced Accuracy Rate (BAR)
calculate_BAR <- function(beta0, beta_est) {
  # OLS fit doesn't give exact zeros, so we only compute BAR for sparse/regularized estimators
  tp <- sum(beta0 != 0 & beta_est != 0)
  tn <- sum(beta0 == 0 & beta_est == 0)
  fp <- sum(beta0 == 0 & beta_est != 0)
  fn <- sum(beta0 != 0 & beta_est == 0)
  
  sensitivity <- if (tp + fn > 0) tp / (tp + fn) else 0
  specificity <- if (tn + fp > 0) tn / (tn + fp) else 0
  
  return(0.5 * (sensitivity + specificity))
}

# Define Scenarios
scenarios <- list(
  list(name = "Scenario 1 (Case I)", n = 50, p = 20, rho = 0.0, sigma = 3.0, design = "IS", 
       beta0 = c(5, rep(0, 19))),
  list(name = "Scenario 2 (Case III)", n = 100, p = 50, rho = 0.95, sigma = 1.5, design = "AR", 
       beta0 = c(5, rep(0, 49))),
  list(name = "Scenario 3 (Case VIII)", n = 50, p = 20, rho = 0.5, sigma = 3.0, design = "CS", 
       beta0 = c(3, 1.5, 0, 0, 2, rep(0, 15)))
)

n_rep <- 5  # Number of replications (can be increased by user)
max_steps <- 500
burn_in <- 100

methods <- c("Lasso", "ALasso", "BLasso", "Horseshoe", "BayesA", "BayesB", "BayesC", "rLASSO (S5)")

# Results structure: list of matrices, one for each scenario
results_mse <- list()
results_bar <- list()

set.seed(42)

for (s_idx in seq_along(scenarios)) {
  scen <- scenarios[[s_idx]]
  cat("\nRunning", scen$name, "...\n")
  
  mse_mat <- matrix(0, nrow = n_rep, ncol = length(methods), dimnames = list(NULL, methods))
  bar_mat <- matrix(0, nrow = n_rep, ncol = length(methods), dimnames = list(NULL, methods))
  
  for (r in 1:n_rep) {
    cat("Replication", r, "/", n_rep, "\r")
    
    # 1. Generate Data
    X_train <- generate_X(scen$n, scen$p, scen$rho, scen$design)
    y_train <- X_train %*% scen$beta0 + rnorm(scen$n, mean = 0, sd = scen$sigma)
    
    # Center y and scale X
    X_train_std <- scale(X_train)
    y_train_ctr <- as.vector(y_train - mean(y_train))
    
    # Generate Test Data
    X_test <- generate_X(200, scen$p, scen$rho, scen$design)
    y_test <- X_test %*% scen$beta0 + rnorm(200, mean = 0, sd = scen$sigma)
    X_test_std <- scale(X_test, center = colMeans(X_train), scale = apply(X_train, 2, sd))
    y_test_ctr <- as.vector(y_test - mean(y_train))
    
    # --- 1. Lasso ---
    lasso_cv <- cv.glmnet(X_train_std, y_train_ctr, alpha = 1)
    beta_lasso <- as.vector(predict(lasso_cv, type = "coefficients", s = "lambda.min"))[-1]
    y_pred_lasso <- X_test_std %*% beta_lasso
    mse_mat[r, "Lasso"] <- mean((y_test_ctr - y_pred_lasso)^2)
    bar_mat[r, "Lasso"] <- calculate_BAR(scen$beta0, beta_lasso)
    
    # --- 2. Adaptive Lasso ---
    # Ridge regression for weights
    ridge_cv <- cv.glmnet(X_train_std, y_train_ctr, alpha = 0)
    beta_ridge <- as.vector(predict(ridge_cv, type = "coefficients", s = "lambda.min"))[-1]
    w <- 1 / (abs(ridge_cv$glmnet.fit$beta[, which(ridge_cv$lambda == ridge_cv$lambda.min)]) + 1e-5)
    alasso_cv <- cv.glmnet(X_train_std, y_train_ctr, alpha = 1, penalty.factor = w)
    beta_alasso <- as.vector(predict(alasso_cv, type = "coefficients", s = "lambda.min"))[-1]
    y_pred_alasso <- X_test_std %*% beta_alasso
    mse_mat[r, "ALasso"] <- mean((y_test_ctr - y_pred_alasso)^2)
    bar_mat[r, "ALasso"] <- calculate_BAR(scen$beta0, beta_alasso)
    
    # --- 3. Bayesian Lasso ---
    bl_fit <- BayesLasso(X_train_std, y_train_ctr, max.steps = max_steps, n.burn = burn_in, seed = r)
    # Variable selection by 95% credible interval
    beta_bl_sel <- variable_selection_CI(bl_fit$beta.post)
    y_pred_bl <- X_test_std %*% bl_fit$beta
    mse_mat[r, "BLasso"] <- mean((y_test_ctr - y_pred_bl)^2)
    bar_mat[r, "BLasso"] <- calculate_BAR(scen$beta0, beta_bl_sel)
    
    # --- 4. Horseshoe ---
    hs_fit <- BayesHorseshoe(X_train_std, y_train_ctr, max.steps = max_steps, n.burn = burn_in, seed = r)
    # Variable selection by 95% credible interval
    beta_hs_sel <- variable_selection_CI(hs_fit$beta.post)
    y_pred_hs <- X_test_std %*% hs_fit$beta
    mse_mat[r, "Horseshoe"] <- mean((y_test_ctr - y_pred_hs)^2)
    bar_mat[r, "Horseshoe"] <- calculate_BAR(scen$beta0, beta_hs_sel)
    
    # --- 5. BayesA ---
    ba_fit <- BayesRLasso(X_train_std, y_train_ctr, lambda.estimate = "AP", 
                         max.steps = max_steps, n.burn = burn_in, seed = r)
    beta_ba_sel <- variable_selection_FBP(X_train_std, y_train_ctr, ba_fit$beta.post, ba_fit$lambda.post)
    y_pred_ba <- X_test_std %*% ba_fit$beta
    mse_mat[r, "BayesA"] <- mean((y_test_ctr - y_pred_ba)^2)
    bar_mat[r, "BayesA"] <- calculate_BAR(scen$beta0, beta_ba_sel)
    
    # --- 6. BayesB ---
    bb_fit <- BayesRLasso(X_train_std, y_train_ctr, lambda.estimate = "EB", 
                         max.steps = max_steps, n.burn = burn_in, seed = r)
    beta_bb_sel <- variable_selection_FBP(X_train_std, y_train_ctr, bb_fit$beta.post, bb_fit$lambda.post)
    y_pred_bb <- X_test_std %*% bb_fit$beta
    mse_mat[r, "BayesB"] <- mean((y_test_ctr - y_pred_bb)^2)
    bar_mat[r, "BayesB"] <- calculate_BAR(scen$beta0, beta_bb_sel)
    
    # --- 7. BayesC ---
    bc_fit <- BayesRLasso(X_train_std, y_train_ctr, lambda.estimate = "MCMC", 
                         max.steps = max_steps, n.burn = burn_in, seed = r)
    beta_bc_sel <- variable_selection_FBP(X_train_std, y_train_ctr, bc_fit$beta.post, bc_fit$lambda.post)
    y_pred_bc <- X_test_std %*% bc_fit$beta
    mse_mat[r, "BayesC"] <- mean((y_test_ctr - y_pred_bc)^2)
    bar_mat[r, "BayesC"] <- calculate_BAR(scen$beta0, beta_bc_sel)
    
    # --- 8. rLASSO (S5) ---
    s5_fit <- rrLASSO.S5(X_train_std, y_train_ctr, lam = 1.0, intercept = FALSE, seed = r, IT = 3, ITER = 3)
    beta_s5 <- s5_fit$beta
    y_pred_s5 <- X_test_std %*% beta_s5
    mse_mat[r, "rLASSO (S5)"] <- mean((y_test_ctr - y_pred_s5)^2)
    bar_mat[r, "rLASSO (S5)"] <- calculate_BAR(scen$beta0, beta_s5)
  }
  
  results_mse[[scen$name]] <- mse_mat
  results_bar[[scen$name]] <- bar_mat
  
  cat("\nScenario results (Mean MSE):\n")
  print(colMeans(mse_mat))
  cat("Scenario results (Mean BAR):\n")
  print(colMeans(bar_mat))
}

# Save results to a file for table generation
save(results_mse, results_bar, scenarios, file = "results/simulation_results.RData")
cat("\nSimulation study completed! Results saved to results/simulation_results.RData\n")
