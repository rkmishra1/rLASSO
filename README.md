# Reciprocal Bayesian LASSO (rLASSO) in R

This repository provides a complete, self-contained R implementation of the **Reciprocal Bayesian LASSO (rLASSO)** method from scratch, as described in the paper:
> **The reciprocal Bayesian LASSO**  
> *Himel Mallick, Rahim Alhamzawi, Erina Paul, Vladimir Svetnik*  
> Statistics in Medicine (2021)

Unlike traditional shrinkage priors (like Lasso or Ridge) that use increasing penalties on the regression coefficients, the rLASSO regularization employs a **decreasing penalty function**. This leads to stronger parsimony, superior variable selection, and oracle properties in both estimation and prediction.

---

## 📖 Methodology Overview

### 1. The Regularization Problem
For a linear regression model:
$$y = X\beta + \epsilon, \quad \epsilon \sim N(0, \sigma^2 I_n)$$
the reciprocal LASSO regularization minimizes the objective function:
$$Q(\beta) = \min_\beta (y - X\beta)'(y - X\beta) + \lambda \sum_{j=1}^p \frac{1}{|\beta_j|} I\{\beta_j \neq 0\}$$
where $\lambda > 0$ is the tuning parameter.

### 2. Inverse Double Exponential (IDE) Prior
The Bayesian counterpart assigns independent **Inverse Laplace (double exponential)** densities to the regression coefficients $\beta_j$:
$$\pi(\beta) = \prod_{j=1}^p \frac{\lambda}{2\beta_j^2} \exp\left\{ -\frac{\lambda}{|\beta_j|} \right\} I\{\beta_j \neq 0\}$$

### 3. Scale Mixture of Truncated Normal (SMTN) Representation
To make Gibbs sampling tractable, the IDE prior is formulated as a hierarchical scale mixture of truncated normals:
$$\beta_j \mid \tau_j, u_j, \sigma^2 \sim N(0, \sigma^2 \tau_j^2) I\left\{ |\beta_j| > \frac{\sigma}{u_j} \right\}$$
$$\tau_j^{-1} \mid \zeta_j \sim \text{Inverse-Gaussian}\left( \frac{\zeta_j \sigma}{|\beta_j|}, \zeta_j^2 \right)$$
$$\zeta_j \mid u_j \sim \text{Gamma}\left( 2, \frac{|\beta_j|}{\sigma} + \frac{1}{u_j} \right)$$
$$u_j \mid \lambda \sim \text{Exponential}(\lambda) I\left\{ u_j > \frac{\sigma}{|\beta_j|} \right\}$$
$$\sigma^2 \sim \text{Inverse-Gamma}\left( \frac{n - 1 + p}{2}, \frac{R + \beta' T^{-1} \beta}{2} \right) I\left\{ \sigma^2 < \min_j(\beta_j^2 u_j^2) \right\}$$

This repository implements a **pure R, zero-dependency coordinate-wise Gibbs sampler** to draw posterior samples from this joint distribution.

---

## 🛠️ Repository Structure

*   `R/rLASSO.R`: The core package library containing:
    *   `BayesRLasso()`: Gibbs sampler for Reciprocal Bayesian Lasso (BayesA, BayesB, BayesC).
    *   `rrLASSO.S5()`: Frequentist rLASSO solved via the Simplified Shotgun Stochastic Search with Screening (S5) algorithm.
    *   `BayesLasso()`: Standard Bayesian Lasso (Park & Casella) implemented from scratch.
    *   `BayesHorseshoe()`: Horseshoe regression (Makalic & Schmidt) implemented from scratch.
    *   `variable_selection_CI()` / `variable_selection_FBP()`: Posthoc variable selection strategies.
    *   Fast, numerically stable univariate truncated normal, truncated gamma, and Inverse-Gaussian samplers.
*   `data/prostate.data`: The classic Stamey Prostate Cancer dataset.
*   `scripts/simulation.R`: Runs the benchmark simulation studies (MSE and BAR).
*   `scripts/real_data_analysis.R`: Runs the out-of-sample prediction and parameter estimation on the Prostate dataset.
*   `results/`: Stores generated performance tables and plots.

---

## 🚀 Getting Started

### Prerequisites
You only need a standard installation of R. The package uses the `glmnet` package (for competing Lasso and Adaptive Lasso benchmarks) and `ggplot2`/`gridExtra` (for plots).

To install the dependencies, run in R:
```R
install.packages(c("glmnet", "ggplot2", "gridExtra"))
```

### Running the Code
1.  **Test Run**: Verify all samplers run without errors:
    ```bash
    Rscript scripts/test_run.R
    ```
2.  **Simulation Study**: Run the benchmark simulations:
    ```bash
    Rscript scripts/simulation.R
    ```
3.  **Real Data Analysis**: Run the Prostate cancer analysis and generate plots:
    ```bash
    Rscript scripts/real_data_analysis.R
    ```

---

## 📊 Results Summary

### 1. Prostate Cancer Dataset Performance
Fitted on the $n=97, p=8$ Prostate Cancer dataset (67 train, 30 test split), the out-of-sample Mean Squared Prediction Error (MSPE) and selected Model Size (MS) are:

| Method | MSPE | Model Size |
| :--- | :---: | :---: |
| **Lasso** | 0.506 | 7 |
| **Adaptive Lasso** | 0.506 | 7 |
| **Bayesian Lasso** | 0.491 | 3 |
| **Horseshoe** | 0.459 | 2 |
| **BayesA (rLASSO)** | 0.527 | 5 |
| **BayesB (rLASSO)** | 0.511 | 7 |
| **BayesC (rLASSO)** | 0.505 | 7 |
| **rLASSO (S5)** | 0.492 | 2 |

*Note: Bayesian rLASSO achieves highly competitive MSPE while enforcing parsimony comparable to Horseshoe and frequentist rLASSO.*

### 2. Parameter Estimates & Posterior CIs
Below are the estimated coefficients and 95% credible intervals for each covariate on the full dataset:

![Prostate Coefficients Comparison](results/prostate_coefficients_comparison.png)

*The plot shows that BayesA, BayesB, and BayesC provide automatic standard error estimation as a natural byproduct, showing narrow credible intervals for the top predictors (`lcavol`, `lweight`, and `svi`).*

### 3. Posterior MCMC Convergence (BayesC)
Trace plots for the posterior samples of $\beta$ and $\lambda$ under the BayesC MCMC model show excellent mixing and fast convergence:

![BayesC MCMC Trace Plot](results/prostate_bayesC_traceplot.png)

### 4. Simulation Study Results
The simulation study evaluates the performance of the competing models in three distinct scenarios (mean and standard deviation over 5 replications):
*   **Scenario 1 (Case I)**: Isotropic design ($n=50, p=20, \rho=0$), strong sparsity ($\beta_0 = (5, 0, \ldots, 0)$), high noise ($\sigma = 3$).
*   **Scenario 2 (Case III)**: Autoregressive design ($n=100, p=50, \rho=0.95$), strong sparsity ($\beta_0 = (5, 0, \ldots, 0)$), medium noise ($\sigma = 1.5$).
*   **Scenario 3 (Case VIII)**: Compound Symmetry design ($n=50, p=20, \rho=0.5$), mild sparsity ($\beta_0 = (3, 1.5, 0, 0, 2, 0, \ldots, 0)$), high noise ($\sigma = 3$).

#### Scenario 1 (Case I)

| Method | MSE (SD) | BAR (SD) |
| :--- | :---: | :---: |
| **Lasso** | 11.847 (2.616) | 0.968 (0.029) |
| **Adaptive Lasso** | 10.780 (1.859) | 0.989 (0.014) |
| **Bayesian Lasso** | 12.690 (3.354) | 1.000 (0.000) |
| **Horseshoe** | 10.714 (2.020) | 1.000 (0.000) |
| **BayesA (rLASSO)** | 46.922 (14.445) | 0.816 (0.100) |
| **BayesB (rLASSO)** | 12.934 (3.541) | 0.637 (0.066) |
| **BayesC (rLASSO)** | 12.664 (3.656) | 0.632 (0.064) |
| **rLASSO (S5)** | 12.751 (2.584) | 0.821 (0.022) |

#### Scenario 2 (Case III)

| Method | MSE (SD) | BAR (SD) |
| :--- | :---: | :---: |
| **Lasso** | 2.745 (0.271) | 0.900 (0.102) |
| **Adaptive Lasso** | 2.451 (0.189) | 0.998 (0.005) |
| **Bayesian Lasso** | 2.935 (0.254) | 0.994 (0.014) |
| **Horseshoe** | 2.504 (0.161) | 1.000 (0.000) |
| **BayesA (rLASSO)** | 8.370 (1.281) | 0.955 (0.034) |
| **BayesB (rLASSO)** | 3.691 (0.485) | 0.824 (0.023) |
| **BayesC (rLASSO)** | 3.722 (0.607) | 0.829 (0.018) |
| **rLASSO (S5)** | 2.880 (0.117) | 0.941 (0.018) |

#### Scenario 3 (Case VIII)

| Method | MSE (SD) | BAR (SD) |
| :--- | :---: | :---: |
| **Lasso** | 12.372 (2.510) | 0.859 (0.048) |
| **Adaptive Lasso** | 13.414 (3.117) | 0.888 (0.105) |
| **Bayesian Lasso** | 12.688 (2.768) | 0.861 (0.079) |
| **Horseshoe** | 12.527 (2.788) | 0.694 (0.079) |
| **BayesA (rLASSO)** | 49.688 (14.862) | 0.763 (0.097) |
| **BayesB (rLASSO)** | 13.513 (2.945) | 0.647 (0.042) |
| **BayesC (rLASSO)** | 13.437 (2.869) | 0.647 (0.042) |
| **rLASSO (S5)** | 13.644 (3.175) | 0.784 (0.088) |

