# Real Data Analysis on the Stamey Prostate Cancer Dataset
# Fits BayesA, BayesB, BayesC, Horseshoe, BayesLasso, Lasso, Adaptive Lasso, and Frequentist rLASSO
# Computes MSPE (out-of-sample) and Model Size (MS) using the train/test split
# Fits on the full dataset to plot posterior means, 95% CIs, and frequentist rLASSO estimates

source("R/rLASSO.R")
library(glmnet)
library(ggplot2)
library(gridExtra)

# 1. Load Data
df <- read.table("data/prostate.data", header = TRUE, row.names = 1)
y <- df$lpsa
X <- as.matrix(df[, 1:8])

# Split by the 'train' indicator
train_idx <- as.logical(df$train)
X_train <- X[train_idx, ]
y_train <- y[train_idx]
X_test <- X[!train_idx, ]
y_test <- y[!train_idx]

# Center response and standardize predictors
y_train_ctr <- as.vector(y_train - mean(y_train))
y_test_ctr <- as.vector(y_test - mean(y_train)) # center test response using training mean

mean_X_train <- colMeans(X_train)
sd_X_train <- apply(X_train, 2, sd)

X_train_std <- scale(X_train, center = mean_X_train, scale = sd_X_train)
X_test_std <- scale(X_test, center = mean_X_train, scale = sd_X_train)

# Methods list
methods <- c("Lasso", "ALasso", "BLasso", "Horseshoe", "BayesA", "BayesB", "BayesC", "rLASSO (S5)")
mspe <- rep(0, length(methods))
model_size <- rep(0, length(methods))
names(mspe) <- names(model_size) <- methods

max_steps <- 11000
burn_in <- 1000

cat("\n--- Running Real Data Analysis (Train/Test Split) ---\n")

# --- 1. Lasso ---
lasso_cv <- cv.glmnet(X_train_std, y_train_ctr, alpha = 1)
beta_lasso <- as.vector(predict(lasso_cv, type = "coefficients", s = "lambda.min"))[-1]
y_pred_lasso <- X_test_std %*% beta_lasso
mspe["Lasso"] <- mean((y_test_ctr - y_pred_lasso)^2)
model_size["Lasso"] <- sum(beta_lasso != 0)

# --- 2. Adaptive Lasso ---
ridge_cv <- cv.glmnet(X_train_std, y_train_ctr, alpha = 0)
w <- 1 / (abs(as.vector(predict(ridge_cv, type = "coefficients", s = "lambda.min"))[-1]) + 1e-5)
alasso_cv <- cv.glmnet(X_train_std, y_train_ctr, alpha = 1, penalty.factor = w)
beta_alasso <- as.vector(predict(alasso_cv, type = "coefficients", s = "lambda.min"))[-1]
y_pred_alasso <- X_test_std %*% beta_alasso
mspe["ALasso"] <- mean((y_test_ctr - y_pred_alasso)^2)
model_size["ALasso"] <- sum(beta_alasso != 0)

# --- 3. Bayesian Lasso ---
bl_fit <- BayesLasso(X_train_std, y_train_ctr, max.steps = max_steps, n.burn = burn_in, seed = 1234)
beta_bl_sel <- variable_selection_CI(bl_fit$beta.post)
y_pred_bl <- X_test_std %*% bl_fit$beta
mspe["BLasso"] <- mean((y_test_ctr - y_pred_bl)^2)
model_size["BLasso"] <- sum(beta_bl_sel != 0)

# --- 4. Horseshoe ---
hs_fit <- BayesHorseshoe(X_train_std, y_train_ctr, max.steps = max_steps, n.burn = burn_in, seed = 1234)
beta_hs_sel <- variable_selection_CI(hs_fit$beta.post)
y_pred_hs <- X_test_std %*% hs_fit$beta
mspe["Horseshoe"] <- mean((y_test_ctr - y_pred_hs)^2)
model_size["Horseshoe"] <- sum(beta_hs_sel != 0)

# --- 5. BayesA ---
ba_fit <- BayesRLasso(X_train_std, y_train_ctr, lambda.estimate = "AP", 
                     max.steps = max_steps, n.burn = burn_in, seed = 1234)
beta_ba_sel <- variable_selection_FBP(X_train_std, y_train_ctr, ba_fit$beta.post, ba_fit$lambda.post)
y_pred_ba <- X_test_std %*% ba_fit$beta
mspe["BayesA"] <- mean((y_test_ctr - y_pred_ba)^2)
model_size["BayesA"] <- sum(beta_ba_sel != 0)

# --- 6. BayesB ---
bb_fit <- BayesRLasso(X_train_std, y_train_ctr, lambda.estimate = "EB", 
                     max.steps = max_steps, n.burn = burn_in, seed = 1234)
beta_bb_sel <- variable_selection_FBP(X_train_std, y_train_ctr, bb_fit$beta.post, bb_fit$lambda.post)
y_pred_bb <- X_test_std %*% bb_fit$beta
mspe["BayesB"] <- mean((y_test_ctr - y_pred_bb)^2)
model_size["BayesB"] <- sum(beta_bb_sel != 0)

# --- 7. BayesC ---
bc_fit <- BayesRLasso(X_train_std, y_train_ctr, lambda.estimate = "MCMC", 
                     max.steps = max_steps, n.burn = burn_in, seed = 1234)
beta_bc_sel <- variable_selection_FBP(X_train_std, y_train_ctr, bc_fit$beta.post, bc_fit$lambda.post)
y_pred_bc <- X_test_std %*% bc_fit$beta
mspe["BayesC"] <- mean((y_test_ctr - y_pred_bc)^2)
model_size["BayesC"] <- sum(beta_bc_sel != 0)

# --- 8. rLASSO (S5) ---
s5_fit <- rrLASSO.S5(X_train_std, y_train_ctr, lam = 1.0, intercept = FALSE, seed = 1234)
beta_s5 <- s5_fit$beta
y_pred_s5 <- X_test_std %*% beta_s5
mspe["rLASSO (S5)"] <- mean((y_test_ctr - y_pred_s5)^2)
model_size["rLASSO (S5)"] <- sum(beta_s5 != 0)

# Show performance metrics
cat("\nPerformance on Prostate Cancer Dataset (MSPE and Model Size):\n")
res_table <- data.frame(MSPE = mspe, ModelSize = model_size)
print(res_table)
write.csv(res_table, "results/prostate_performance_table.csv")

# -------------------------------------------------------------
# Fit on Full Dataset for Parameter Estimation & Figure 2
# -------------------------------------------------------------

cat("\n--- Fitting on Full Dataset for Parameter Estimation (Figure 2) ---\n")

y_full_ctr <- as.vector(y - mean(y))
X_full_std <- scale(X)

# Fit BayesA
fit_A <- BayesRLasso(X_full_std, y_full_ctr, lambda.estimate = "AP", max.steps = max_steps, n.burn = burn_in, seed = 1234)
# Fit BayesB
fit_B <- BayesRLasso(X_full_std, y_full_ctr, lambda.estimate = "EB", max.steps = max_steps, n.burn = burn_in, seed = 1234)
# Fit BayesC
fit_C <- BayesRLasso(X_full_std, y_full_ctr, lambda.estimate = "MCMC", max.steps = max_steps, n.burn = burn_in, seed = 1234)
# Fit Horseshoe
fit_HS <- BayesHorseshoe(X_full_std, y_full_ctr, max.steps = max_steps, n.burn = burn_in, seed = 1234)
# Fit Frequentist rLASSO S5 (using BayesC's median lambda for FBP consistency)
fit_rLASSO <- rrLASSO.S5(X_full_std, y_full_ctr, lam = fit_C$lambda, intercept = FALSE, seed = 1234)

# Create coefficient comparison data frame
covariates <- colnames(X)
n_covs <- length(covariates)

plot_data <- data.frame(
  Covariate = rep(covariates, 5),
  Estimate = c(fit_A$beta, fit_B$beta, fit_C$beta, fit_HS$beta, fit_rLASSO$beta),
  Lower = c(fit_A$lowerbeta, fit_B$lowerbeta, fit_C$lowerbeta, 
            apply(fit_HS$beta.post, 2, quantile, probs = 0.025), 
            fit_rLASSO$beta),
  Upper = c(fit_A$upperbeta, fit_B$upperbeta, fit_C$upperbeta, 
            apply(fit_HS$beta.post, 2, quantile, probs = 0.975), 
            fit_rLASSO$beta),
  Method = rep(c("BayesA", "BayesB", "BayesC", "Horseshoe", "rLASSO"), each = n_covs)
)

# Plot Figure 2 (posterior mean and 95% CIs)
p2 <- ggplot(plot_data, aes(x = Covariate, y = Estimate, color = Method, shape = Method)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position = position_dodge(width = 0.6), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  theme_minimal() +
  labs(
    title = "Posterior Mean and 95% Credible Intervals for Prostate Cancer Covariates",
    x = "Covariates (Standardized)",
    y = "Coefficient Estimates",
    color = "Method",
    shape = "Method"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 12, face = "bold")
  )

ggsave("results/prostate_coefficients_comparison.png", plot = p2, width = 10, height = 6, dpi = 300)

# Save parameter estimation table
param_table <- data.frame(
  Covariate = covariates,
  BayesA = sprintf("%.3f (%.3f, %.3f)", fit_A$beta, fit_A$lowerbeta, fit_A$upperbeta),
  BayesB = sprintf("%.3f (%.3f, %.3f)", fit_B$beta, fit_B$lowerbeta, fit_B$upperbeta),
  BayesC = sprintf("%.3f (%.3f, %.3f)", fit_C$beta, fit_C$lowerbeta, fit_C$upperbeta),
  Horseshoe = sprintf("%.3f (%.3f, %.3f)", fit_HS$beta, 
                      apply(fit_HS$beta.post, 2, quantile, probs = 0.025), 
                      apply(fit_HS$beta.post, 2, quantile, probs = 0.975)),
  rLASSO_S5 = sprintf("%.3f", fit_rLASSO$beta)
)
print(param_table)
write.csv(param_table, "results/prostate_parameter_estimates.csv", row.names = FALSE)

# Generate a trace plot for BayesC's beta.post
png("results/prostate_bayesC_traceplot.png", width = 800, height = 600)
par(mfrow = c(3, 3))
for (j in 1:ncol(X)) {
  plot(fit_C$beta.post[, j], type = "l", main = covariates[j], xlab = "Iteration", ylab = "Beta")
}
plot(fit_C$lambda.post, type = "l", main = "lambda", xlab = "Iteration", ylab = "Lambda")
dev.off()

cat("\nReal data analysis completed! Outputs saved to results/\n")
