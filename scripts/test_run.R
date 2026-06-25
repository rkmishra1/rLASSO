# Test run script to verify all core algorithms work without errors

source("R/rLASSO.R")

cat("--- Initializing test run ---\n")
set.seed(123)

n <- 30
p <- 5
beta0 <- c(3, 1.5, 0, 0, 0)
X <- matrix(rnorm(n * p), n, p)
y <- X %*% beta0 + rnorm(n)

X_std <- scale(X)
y_ctr <- as.vector(y - mean(y))

steps <- 100
burn <- 20

# 1. Test BayesRLasso (AP)
cat("Testing BayesRLasso (AP)... ")
fit_ap <- BayesRLasso(X_std, y_ctr, lambda.estimate = "AP", max.steps = steps, n.burn = burn)
cat("OK\n")

# 2. Test BayesRLasso (EB)
cat("Testing BayesRLasso (EB)... ")
fit_eb <- BayesRLasso(X_std, y_ctr, lambda.estimate = "EB", max.steps = steps, n.burn = burn)
cat("OK\n")

# 3. Test BayesRLasso (MCMC)
cat("Testing BayesRLasso (MCMC)... ")
fit_mcmc <- BayesRLasso(X_std, y_ctr, lambda.estimate = "MCMC", max.steps = steps, n.burn = burn)
cat("OK\n")

# 4. Test rrLASSO.S5
cat("Testing rrLASSO.S5... ")
fit_s5 <- rrLASSO.S5(X_std, y_ctr, S = 3, lam = 1.0, intercept = FALSE)
cat("OK\n")

# 5. Test BayesLasso
cat("Testing BayesLasso... ")
fit_bl <- BayesLasso(X_std, y_ctr, max.steps = steps, n.burn = burn)
cat("OK\n")

# 6. Test BayesHorseshoe
cat("Testing BayesHorseshoe... ")
fit_hs <- BayesHorseshoe(X_std, y_ctr, max.steps = steps, n.burn = burn)
cat("OK\n")

# 7. Test Variable Selection CI
cat("Testing Variable Selection (CI)... ")
beta_ci <- variable_selection_CI(fit_mcmc$beta.post)
cat("OK\n")

# 8. Test Variable Selection FBP
cat("Testing Variable Selection (FBP)... ")
beta_fbp <- variable_selection_FBP(X_std, y_ctr, fit_mcmc$beta.post, fit_mcmc$lambda.post)
cat("OK\n")

cat("\n--- All tests completed successfully! ---\n")
