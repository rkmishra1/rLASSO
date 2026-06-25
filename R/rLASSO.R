# Reciprocal Bayesian LASSO and Competing Methods in R
# Implements:
# 1. Truncated Normal and Gamma Samplers (pure R)
# 2. Inverse Gaussian Sampler (pure R)
# 3. Reciprocal Bayesian LASSO (BayesA, BayesB, BayesC) using coordinate-wise Gibbs
# 4. Frequentist rLASSO solved via the S5 algorithm
# 5. Bayesian LASSO (Park & Casella)
# 6. Horseshoe Regression (Makalic & Schmidt)
# 7. Posthoc Variable Selection (CI, FBP)

# -------------------------------------------------------------
# 1. Truncated Normal and Gamma Samplers
# -------------------------------------------------------------

# Sample from standard normal truncated below at a: Z ~ N(0, 1) I(Z > a)
rtruncnorm_left <- function(a) {
  if (is.na(a) || is.nan(a)) return(rnorm(1))
  if (a <= 1) {
    while (TRUE) {
      z <- rnorm(1)
      if (z > a) return(z)
    }
  } else {
    # Shifted exponential proposal (Robert, 1995)
    # Optimal rate is alpha = a
    while (TRUE) {
      e <- rexp(1, rate = a)
      z <- a + e
      if (runif(1) <= exp(-0.5 * e^2)) {
        return(z)
      }
    }
  }
}

# Sample from N(mu, s^2) truncated to |x| > c for c > 0
rtruncnorm_mid <- function(mu, s, c) {
  if (c <= 0) return(rnorm(1, mean = mu, sd = s))
  
  z_left <- (-c - mu) / s
  z_right <- (c - mu) / s
  
  log_p_left <- pnorm(z_left, log.p = TRUE)
  log_p_right <- pnorm(z_right, lower.tail = FALSE, log.p = TRUE)
  
  max_log <- max(log_p_left, log_p_right)
  p_l <- exp(log_p_left - max_log)
  p_r <- exp(log_p_right - max_log)
  p_sum <- p_l + p_r
  
  if (is.nan(p_sum) || p_sum <= 0) {
    # Fallback if both tails are extremely far out
    if (z_left > -z_right) {
      z <- -rtruncnorm_left(-z_left)
    } else {
      z <- rtruncnorm_left(z_right)
    }
  } else {
    if (runif(1) * p_sum < p_l) {
      z <- -rtruncnorm_left(-z_left)
    } else {
      z <- rtruncnorm_left(z_right)
    }
  }
  return(mu + s * z)
}

rtgamma <- function(shape, rate, lower) {
  if (lower <= 0) {
    return(rgamma(1, shape = shape, rate = rate))
  }
  p_lower <- pgamma(lower, shape = shape, rate = rate)
  if (p_lower < 1 - 1e-10) {
    u <- runif(1, min = p_lower, max = 1)
    val <- qgamma(u, shape = shape, rate = rate)
    if (is.finite(val) && val >= lower) {
      return(val)
    }
  }
  
  # Fallback:
  lambda <- rate - (shape - 1) / lower
  if (lambda <= 0) {
    # Simple rejection sampling
    while (TRUE) {
      w <- rgamma(1, shape = shape, rate = rate)
      if (w > lower) return(w)
    }
  }
  
  # Shifted exponential proposal: w = lower + e, where e ~ Exp(lambda)
  log_ratio_max <- (shape - 1) * log(lower) - (rate - lambda) * lower
  max_attempts <- 1000
  for (i in 1:max_attempts) {
    w <- lower + rexp(1, rate = lambda)
    log_ratio <- (shape - 1) * log(w) - (rate - lambda) * w
    if (log(runif(1)) <= log_ratio - log_ratio_max) {
      return(w)
    }
  }
  return(lower + rexp(1, rate = rate))
}

# -------------------------------------------------------------
# 2. Inverse Gaussian Sampler
# -------------------------------------------------------------

# Sample from Inverse Gaussian IG(mu, lambda)
rinvgauss_single <- function(mu, lambda) {
  if (is.infinite(mu) || is.na(mu) || is.nan(mu) || mu <= 0 || lambda <= 0) {
    return(rexp(1))
  }
  y <- rnorm(1)^2
  mu2 <- mu^2
  x <- mu + (mu2 * y) / (2 * lambda) - (mu / (2 * lambda)) * sqrt(4 * mu * lambda * y + mu2 * y^2)
  if (is.na(x) || is.nan(x) || !is.finite(x) || x <= 0) {
    return(mu)
  }
  if (runif(1) <= mu / (mu + x)) {
    return(x)
  } else {
    return(mu2 / x)
  }
}

# -------------------------------------------------------------
# 3. Multivariate Normal Sampler
# -------------------------------------------------------------

# Single draw from MVN(mean, sigma)
rmvnorm_single <- function(mean, sigma) {
  p <- length(mean)
  z <- rnorm(p)
  ch <- tryCatch({
    chol(sigma)
  }, error = function(e) {
    eg <- eigen(sigma, symmetric = TRUE)
    eg$vectors %*% diag(sqrt(pmax(eg$values, 0)), p)
  })
  if (is.matrix(ch)) {
    return(mean + t(ch) %*% z)
  } else {
    return(mean + ch %*% z)
  }
}

# -------------------------------------------------------------
# 4. Helper for Apriori Estimation of Lambda
# -------------------------------------------------------------

hyper_par_BayesRLasso <- function(x, y, threshold = NULL) {
  n <- nrow(x)
  p <- ncol(x)
  threshold <- ifelse(is.null(threshold), p^-0.5, threshold)
  
  # Generate null distribution of beta
  betas <- matrix(0, 3, 1000) # Use 1000 instead of 10000 for massive speedup
  for (k in 1:1000) {
    sam <- sample(1:p, 3)
    ind <- sample(1:n, n, replace = TRUE)
    XtX_sub <- crossprod(x[ind, sam])
    fit <- solve(XtX_sub + diag(1e-5, 3)) %*% crossprod(x[ind, sam], y[ind])
    betas[, k] <- as.vector(fit)
  }
  betas <- as.vector(betas)
  
  # Calculate candidate lambda
  lambda.cand <- seq(0.01, (sd(y) + 0.1), length.out = 30)^2
  pro <- rep(0, 30)
  den.null1 <- density(betas)
  den.null <- approxfun(den.null1$x, den.null1$y, rule = 1, method = "linear")
  
  for (k in 1:30) {
    lambda <- lambda.cand[k]
    den <- function(val) { lambda * val^-2 * exp(-1 * lambda / abs(val)) / 2 }
    f <- function(val) { 
      d_val <- den(val)
      n_val <- den.null(val)
      n_val[is.na(n_val)] <- 0
      d_val - n_val 
    }
    tryCatch({
      a <- uniroot(f, interval = c(0.001, max(betas)))
      th <- a$root
      loc <- integrate(den.null, lower = th, upper = max(betas) - 0.001)$value
      nonloc <- integrate(den, lower = 0, upper = th)$value
      pro[k] <- loc + nonloc
    }, error = function(e) {
      pro[k] <- 1.0 # fallback
    })
  }
  B <- lambda.cand[which.min((pro - threshold)^2)]
  return(B)
}

# -------------------------------------------------------------
# 5. Reciprocal Bayesian LASSO (BayesA, BayesB, BayesC)
# -------------------------------------------------------------

BayesRLasso <- function(X, y,
                        lambda.estimate = c("AP", "EB", "MCMC"),
                        lambda.fixed = NULL, # used if lambda.estimate = "AP" and not NULL
                        update.sigma2 = TRUE,
                        max.steps = 11000,
                        n.burn = 1000,
                        n.thin = 1,
                        a = 0.001,
                        b = 0.001,
                        posterior.summary.beta = 'mean',
                        posterior.summary.lambda = 'median',
                        beta.ci.level = 0.95,
                        lambda.ci.level = 0.95,
                        seed = 1234,
                        verbose = FALSE) {
  set.seed(seed)
  lambda.estimate <- match.arg(lambda.estimate)
  
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  # Standardize X and center y (assume already done, but enforce matrix type)
  XtX <- t(X) %*% X
  X_t_X_diag <- colSums(X^2)
  xy <- t(X) %*% y
  
  # Initialize coefficients
  # OLS or ridge
  v <- try(solve(XtX), silent = TRUE)
  if (inherits(v, "try-error")) {
    eps <- 1e-05
    beta <- as.vector(solve(XtX + diag(eps, p)) %*% xy)
  } else {
    beta <- as.vector(solve(XtX + diag(1e-5, p)) %*% xy)
  }
  
  # Avoid exact zeros for initializations
  beta[abs(beta) < 1e-5] <- 1e-5 * sign(beta[abs(beta) < 1e-5])
  beta[beta == 0] <- 1e-5
  
  residue <- as.vector(y - X %*% beta)
  sigma2 <- sum(residue^2) / n
  if (sigma2 == 0) sigma2 <- 1.0
  
  # Initialize latent variables
  d <- 1 / (beta^2) # T^-1 diagonal
  u <- 1 / abs(beta)
  zeta <- u
  
  # Initialize lambda
  if (lambda.estimate == "AP") {
    if (!is.null(lambda.fixed)) {
      lambda <- lambda.fixed
    } else {
      if (verbose) cat("Estimating lambda apriori...\n")
      lambda <- hyper_par_BayesRLasso(X, y)
    }
  } else {
    lambda <- 1.0
  }
  
  # Save MCMC samples
  beta_samples <- matrix(0, nrow = max.steps, ncol = p)
  sigma2_samples <- rep(0, max.steps)
  lambda_samples <- rep(0, max.steps)
  
  # Precompute residuals
  r <- residue
  
  # Gibbs Sampler
  for (step in 1:max.steps) {
    if (verbose && step %% 1000 == 0) {
      cat("Gibbs Iteration:", step, "\n")
    }
    
    # 1. Update beta component-by-component (coordinate-wise)
    for (j in 1:p) {
      r_j <- r + X[, j] * beta[j]
      V_j <- X_t_X_diag[j] + d[j]
      s_j <- sqrt(sigma2 / V_j)
      mu_j <- sum(X[, j] * r_j) / V_j
      c_j <- sqrt(sigma2) / u[j]
      
      # Sample beta_j from mid-truncated normal
      beta[j] <- rtruncnorm_mid(mu_j, s_j, c_j)
      r <- r_j - X[, j] * beta[j]
    }
    beta_samples[step, ] <- beta
    
    # 2. Update sigma^2
    if (update.sigma2) {
      shape_sig <- (n + p - 1) / 2
      resid_ss <- sum(r^2)
      prior_ss <- sum(beta^2 * d)
      scale_sig <- (resid_ss + prior_ss) / 2
      
      # Truncation: sigma2 < min_j(beta_j^2 * u_j^2)
      # So W = 1/sigma^2 > 1 / min_j(beta_j^2 * u_j^2)
      U <- min(beta^2 * u^2)
      lower_W <- 1 / U
      
      W <- rtgamma(shape = shape_sig, rate = scale_sig, lower = lower_W)
      sigma2 <- 1 / W
    }
    sigma2_samples[step] <- sigma2
    
    # 3. Update u
    T_vec <- sqrt(sigma2) / abs(beta)
    for (j in 1:p) {
      u[j] <- T_vec[j] + rexp(1, rate = lambda)
    }
    
    # 4. Update zeta
    for (j in 1:p) {
      zeta[j] <- rgamma(1, shape = 2, rate = abs(beta[j]) / sqrt(sigma2) + 1 / u[j])
    }
    
    # 5. Update d (diagonal of T^-1)
    for (j in 1:p) {
      d[j] <- rinvgauss_single(mu = zeta[j] * sqrt(sigma2) / abs(beta[j]), lambda = zeta[j]^2)
    }
    
    # 6. Update lambda
    if (lambda.estimate == "EB") {
      # Empirical Bayes update at each step (online EM)
      lambda <- (2 * p) / sum(u)
    } else if (lambda.estimate == "MCMC") {
      # MCMC update with Gamma prior
      lambda <- rgamma(1, shape = a + 2 * p, rate = b + sum(u))
    }
    lambda_samples[step] <- lambda
  }
  
  # Post-burnin samples
  keep_indices <- seq(n.burn + 1, max.steps, by = n.thin)
  beta_post <- beta_samples[keep_indices, , drop = FALSE]
  sigma2_post <- sigma2_samples[keep_indices]
  lambda_post <- lambda_samples[keep_indices]
  
  # Summarize posterior
  if (posterior.summary.beta == 'mean') {
    beta_est <- colMeans(beta_post)
  } else if (posterior.summary.beta == 'median') {
    beta_est <- apply(beta_post, 2, median)
  } else {
    # Mode estimation via density
    beta_est <- apply(beta_post, 2, function(col) {
      d_est <- density(col)
      d_est$x[which.max(d_est$y)]
    })
  }
  
  if (posterior.summary.lambda == 'mean') {
    lambda_est <- mean(lambda_post)
  } else if (posterior.summary.lambda == 'median') {
    lambda_est <- median(lambda_post)
  } else {
    d_est <- density(lambda_post)
    lambda_est <- d_est$x[which.max(d_est$y)]
  }
  
  # Credible intervals
  lower_beta <- apply(beta_post, 2, quantile, probs = (1 - beta.ci.level)/2)
  upper_beta <- apply(beta_post, 2, quantile, probs = (1 + beta.ci.level)/2)
  lambda_ci <- quantile(lambda_post, probs = c((1 - lambda.ci.level)/2, (1 + lambda.ci.level)/2))
  
  names(beta_est) <- names(lower_beta) <- names(upper_beta) <- colnames(beta_post) <- colnames(X)
  
  return(list(
    beta = beta_est,
    lowerbeta = lower_beta,
    upperbeta = upper_beta,
    lambda = lambda_est,
    lambdaci = lambda_ci,
    beta.post = beta_post,
    sigma2.post = sigma2_post,
    lambda.post = lambda_post
  ))
}

# -------------------------------------------------------------
# 6. Frequentist rLASSO solved via S5 Algorithm
# -------------------------------------------------------------

# Model probability / objective evaluated at the MAP for a subset ind2
# Uses precomputed cross-products to avoid expensive matrix-vector operations inside optim
solve_cubic_nr <- function(z, C) {
  w <- max(z, 0) + C^(1/3)
  for (iter in 1:10) {
    w2 <- w^2
    f_val <- w2 * (w - z) - C
    f_prime <- w * (3 * w - 2 * z)
    if (f_prime == 0) break
    diff <- f_val / f_prime
    w <- w - diff
    if (abs(diff) < 1e-7) break
  }
  return(w)
}

solve_rLASSO_cd <- function(M, v, B, initial_x, max_iter = 50, tol = 1e-6) {
  p <- length(v)
  x <- initial_x
  M_diag <- diag(M)
  for (iter in 1:max_iter) {
    x_old <- x
    for (j in 1:p) {
      if (p > 1) {
        sum_Mjk_xk <- sum(M[j, ] * x) - M_diag[j] * x[j]
      } else {
        sum_Mjk_xk <- 0
      }
      z_j <- (v[j] - sum_Mjk_xk) / M_diag[j]
      if (abs(z_j) < 1e-10) {
        x[j] <- 0
      } else {
        C_j <- B / (2 * M_diag[j])
        x_abs <- solve_cubic_nr(abs(z_j), C_j)
        x[j] <- sign(z_j) * x_abs
      }
    }
    if (sum((x - x_old)^2) < tol) break
  }
  return(x)
}

# Model probability / objective evaluated at the MAP for a subset ind2
# Uses precomputed cross-products to avoid expensive matrix-vector operations inside optim
ind_fun_rLASSO <- function(ind2, B, y, X, p) {
  p.g <- length(ind2)
  if (p.g > 1) {
    X_sub <- X[, ind2, drop = FALSE]
    M <- crossprod(X_sub)
    v <- as.vector(crossprod(X_sub, y))
    yy <- sum(y^2)
    
    # Initial OLS fit (regularized)
    fit <- solve(M + diag(1e-5, p.g)) %*% v
    initial_x <- as.vector(fit)
    
    # Fast coordinate descent
    x_opt <- tryCatch({
      solve_rLASSO_cd(M, v, B, initial_x)
    }, error = function(e) {
      initial_x
    })
    
    # Objective value
    val <- yy - 2 * sum(x_opt * v) + as.numeric(t(x_opt) %*% M %*% x_opt) + sum(B / (abs(x_opt) + 1e-8))
    int <- -1 * val
  } else if (p.g == 1) {
    # 1D optimization
    X_sub <- X[, ind2]
    M <- sum(X_sub^2)
    v <- sum(X_sub * y)
    yy <- sum(y^2)
    
    initial_x <- v / (M + 1e-5)
    
    # Solve cubic for 1D case directly
    z_1 <- v / M
    C_1 <- B / (2 * M)
    x_opt <- sign(z_1) * solve_cubic_nr(abs(z_1), C_1)
    
    val <- yy - 2 * x_opt * v + x_opt^2 * M + B / (abs(x_opt) + 1e-8)
    int <- -1 * val
  } else {
    int <- -1 * sum(y^2)
  }
  return(int)
}

rrLASSO.S5 <- function(X, y,
                       S = 30, # Screening size
                       lam = 1.0, # Penalty parameter
                       intercept = TRUE,
                       verbose = FALSE,
                       seed = 1234,
                       IT = 5,
                       ITER = 5) {
  set.seed(seed)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  a0 <- 0.01
  b0 <- 0.01
  tau <- 1
  tem <- seq(0.4, 1, length.out = IT)^2
  IT.seq <- rep(ITER, IT)
  A3 <- S
  gam <- rep(0, p)
  p.g <- sum(gam)
  ind2 <- which(gam == 1)
  
  curr <- -1000000
  GAM <- matrix(1, nrow = p, ncol = 1)
  OBJ <- -1000000
  OBJ.id <- -1000000
  OBJ.m0 <- matrix(-1000000, nrow = p, ncol = 1)
  OBJ.p0 <- matrix(-1000000, nrow = p, ncol = 1)
  ID <- -100
  ID.obj <- -100
  
  res <- y
  corr <- as.vector(crossprod(res, X))
  ind.ix <- sort.int(abs(corr), decreasing = TRUE, index.return = TRUE)$ix
  s <- ind.ix
  
  size <- min(A3, p)
  IND <- s[1:size]
  p.ind <- length(IND)
  
  B <- lam
  it <- 1
  
  for (it in 1:IT) {
    IT0 <- IT.seq[it]
    pq <- 0
    for (iter in 1:IT0) {
      id <- sum(2^(2 * log(ind2 + 1)))
      id.ind <- which(id == ID)
      leng <- length(id.ind)
      
      if (leng == 0) {
        ID <- c(ID, id)
        C.p <- rep(-100000, p)
        for (i in seq_along(IND)) {
          j <- IND[i]
          if (gam[j] == 0) {
            gam.p <- gam; gam.p[j] <- 1; ind.p <- which(gam.p == 1)
            obj.p <- ind_fun_rLASSO(ind.p, B, y, X, p)
            C.p[j] <- obj.p
          }
        }
        
        C.m <- rep(-100000, p)
        if (p.g > 0) {
          for (i in 1:p.g) {
            j <- ind2[i]
            gam.m <- gam; gam.m[j] <- 0; ind.m <- which(gam.m == 1)
            obj.m <- ind_fun_rLASSO(ind.m, B, y, X, p)
            C.m[j] <- obj.m
          }
        }
        OBJ.p0 <- cb<-cbind(OBJ.p0, C.p)
        OBJ.m0 <- cb<-cbind(OBJ.m0, C.m)
      } else {
        pq <- pq + 1
        C.p <- OBJ.p0[, (id.ind[1])]; C.m <- OBJ.m0[, (id.ind[1])]
      }
      
      prop_p <- exp(tem[it] * (C.p - max(C.p)))
      # Add small perturbation to avoid zero prob
      prop_p <- prop_p / sum(prop_p)
      sample.p <- sample(1:length(prop_p), 1, prob = prop_p)
      obj.p <- C.p[sample.p]
      
      if (p.g > 0) {
        prop_m <- exp(tem[it] * (C.m - max(C.m)))
        prop_m <- prop_m / sum(prop_m)
        sample.m <- sample(1:length(prop_m), 1, prob = prop_m)
        obj.m <- C.m[sample.m]
      } else {
        sample.m <- sample(1:p, 1)
        obj.m <- -100000
      }
      
      l <- 1 / (1 + exp(tem[it] * obj.m - tem[it] * obj.p))
      if (is.na(l)) l <- 0.5
      if (l > runif(1)) {
        gam[sample.p] <- 1; obj <- obj.p; curr <- obj.p
      } else {
        gam[sample.m] <- 0; obj <- obj.m; curr <- obj.m
      }
      ind2 <- which(gam == 1)
      p.g <- sum(gam)
      
      if (p.g > 0) {
        fit <- solve(crossprod(X[, ind2, drop = FALSE]) + diag(1e-5, p.g)) %*% crossprod(X[, ind2, drop = FALSE], y)
        res <- y - X[, ind2, drop = FALSE] %*% fit
        corr <- as.vector(crossprod(res, X))
        ind.ix <- sort.int(abs(corr), decreasing = TRUE, index.return = TRUE)$ix
        s <- c(ind2, ind.ix)
      } else {
        res <- y
        corr <- as.vector(crossprod(res, X))
        ind.ix <- sort.int(abs(corr), decreasing = TRUE, index.return = TRUE)$ix
        s <- ind.ix
      }
      size <- min(A3, p)
      IND <- unique(s[1:size])
      p.ind <- length(IND)
      
      id <- sum(2^(2 * log(ind2 + 1)))
      id.ind <- which(id == ID.obj)
      
      leng <- length(id.ind)
      if (leng == 0) {
        ID.obj <- c(ID.obj, id)
        OBJ <- c(OBJ, obj)
        GAM <- cbind(GAM, gam)
      }
    }
  }
  
  # Select best model
  best_idx <- which.max(OBJ)
  gam0 <- GAM[, best_idx]
  hppm <- which(gam0 == 1)
  
  # Refit selected model using OLS (debiasing)
  beta <- rep(0, p)
  if (length(hppm) >= 1) {
    X_hppm <- X[, hppm, drop = FALSE]
    if (intercept) {
      ols_fit <- lm(y ~ X_hppm)
      beta[hppm] <- ols_fit$coefficients[-1]
      beta <- c(ols_fit$coefficients[1], beta)
      names(beta) <- c('(Intercept)', colnames(X))
    } else {
      ols_fit <- lm(y ~ X_hppm - 1)
      beta[hppm] <- ols_fit$coefficients
      names(beta) <- colnames(X)
    }
  } else {
    if (intercept) {
      beta <- c(mean(y), rep(0, p))
      names(beta) <- c('(Intercept)', colnames(X))
    } else {
      beta <- rep(0, p)
      names(beta) <- colnames(X)
    }
  }
  
  return(list(hppm = hppm, beta = beta))
}

# -------------------------------------------------------------
# 7. Bayesian LASSO (Park & Casella)
# -------------------------------------------------------------

BayesLasso <- function(X, y,
                       a = 1.0,
                       b = 1.0,
                       max.steps = 11000,
                       n.burn = 1000,
                       n.thin = 1,
                       posterior.summary.beta = 'mean',
                       seed = 1234) {
  set.seed(seed)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  XtX <- t(X) %*% X
  xy <- t(X) %*% y
  
  # Initialize
  lambda_sq <- rgamma(1, shape = a, rate = b)
  sigma2 <- runif(1, 0.1, 10)
  tau_sq <- rexp(p, rate = lambda_sq / 2)
  beta <- rep(0, p)
  
  beta_samples <- matrix(0, nrow = max.steps, ncol = p)
  sigma2_samples <- rep(0, max.steps)
  lambda_samples <- rep(0, max.steps)
  
  for (step in 1:max.steps) {
    # 1. Update beta (block update)
    invD <- diag(1 / tau_sq, p)
    invA <- solve(XtX + invD)
    mean_be <- as.vector(invA %*% xy)
    cov_be <- sigma2 * invA
    beta <- rmvnorm_single(mean_be, cov_be)
    beta_samples[step, ] <- beta
    
    # 2. Update tau_sq
    for (j in 1:p) {
      inv_tau_j <- rinvgauss_single(mu = sqrt(lambda_sq * sigma2 / beta[j]^2), lambda = lambda_sq)
      tau_sq[j] <- 1 / inv_tau_j
    }
    
    # 3. Update sigma^2
    shape_sig <- (n + p - 1) / 2
    resid_ss <- sum((y - X %*% beta)^2)
    prior_ss <- sum(beta^2 / tau_sq)
    scale_sig <- (resid_ss + prior_ss) / 2
    sigma2 <- 1 / rgamma(1, shape = shape_sig, rate = scale_sig)
    sigma2_samples[step] <- sigma2
    
    # 4. Update lambda
    shape_lam <- p + a
    scale_lam <- sum(tau_sq) / 2 + b
    lambda_sq <- rgamma(1, shape = shape_lam, rate = scale_lam)
    lambda_samples[step] <- sqrt(lambda_sq)
  }
  
  keep_indices <- seq(n.burn + 1, max.steps, by = n.thin)
  beta_post <- beta_samples[keep_indices, , drop = FALSE]
  sigma2_post <- sigma2_samples[keep_indices]
  lambda_post <- lambda_samples[keep_indices]
  
  if (posterior.summary.beta == 'mean') {
    beta_est <- colMeans(beta_post)
  } else {
    beta_est <- apply(beta_post, 2, median)
  }
  
  names(beta_est) <- colnames(X)
  
  return(list(
    beta = beta_est,
    beta.post = beta_post,
    sigma2.post = sigma2_post,
    lambda.post = lambda_post
  ))
}

# -------------------------------------------------------------
# 8. Horseshoe Regression (Makalic & Schmidt)
# -------------------------------------------------------------

BayesHorseshoe <- function(X, y,
                           max.steps = 11000,
                           n.burn = 1000,
                           n.thin = 1,
                           posterior.summary.beta = 'mean',
                           seed = 1234) {
  set.seed(seed)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  XtX <- t(X) %*% X
  xy <- t(X) %*% y
  
  # Initialize
  beta <- rep(0, p)
  sigma2 <- 1.0
  lambda_sq <- rep(1.0, p)
  nu <- rep(1.0, p)
  tau_sq <- 1.0
  xi <- 1.0
  
  beta_samples <- matrix(0, nrow = max.steps, ncol = p)
  sigma2_samples <- rep(0, max.steps)
  
  for (step in 1:max.steps) {
    # 1. Update beta (block update)
    D_diag <- tau_sq * lambda_sq
    invD <- diag(1 / D_diag, p)
    invA <- solve(XtX + invD)
    mean_be <- as.vector(invA %*% xy)
    cov_be <- sigma2 * invA
    beta <- rmvnorm_single(mean_be, cov_be)
    beta_samples[step, ] <- beta
    
    # 2. Update local shrinkage parameters lambda_sq
    for (j in 1:p) {
      lambda_sq[j] <- 1 / rgamma(1, shape = 1, rate = 1 / nu[j] + beta[j]^2 / (2 * sigma2 * tau_sq))
    }
    
    # 3. Update local auxiliary variables nu
    for (j in 1:p) {
      nu[j] <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / lambda_sq[j])
    }
    
    # 4. Update global variance tau_sq
    tau_sq <- 1 / rgamma(1, shape = (p + 1) / 2, rate = 1 / xi + sum(beta^2 / (lambda_sq)) / (2 * sigma2))
    
    # 5. Update global auxiliary variable xi
    xi <- 1 / rgamma(1, shape = 1, rate = 1 + 1 / tau_sq)
    
    # 6. Update sigma^2
    shape_sig <- (n + p) / 2
    resid_ss <- sum((y - X %*% beta)^2)
    prior_ss <- sum(beta^2 / lambda_sq) / tau_sq
    scale_sig <- (resid_ss + prior_ss) / 2
    sigma2 <- 1 / rgamma(1, shape = shape_sig, rate = scale_sig)
    sigma2_samples[step] <- sigma2
  }
  
  keep_indices <- seq(n.burn + 1, max.steps, by = n.thin)
  beta_post <- beta_samples[keep_indices, , drop = FALSE]
  sigma2_post <- sigma2_samples[keep_indices]
  
  if (posterior.summary.beta == 'mean') {
    beta_est <- colMeans(beta_post)
  } else {
    beta_est <- apply(beta_post, 2, median)
  }
  
  names(beta_est) <- colnames(X)
  
  return(list(
    beta = beta_est,
    beta.post = beta_post,
    sigma2.post = sigma2_post
  ))
}

# -------------------------------------------------------------
# 9. Posthoc Variable Selection
# -------------------------------------------------------------

# Variable selection by 95% equal-tailed credible intervals (Park & Casella)
variable_selection_CI <- function(beta_post, ci_level = 0.95) {
  lower <- apply(beta_post, 2, quantile, probs = (1 - ci_level)/2)
  upper <- apply(beta_post, 2, quantile, probs = (1 + ci_level)/2)
  
  beta_est <- colMeans(beta_post)
  # Set to 0 if the credible interval contains 0
  beta_est[lower <= 0 & upper >= 0] <- 0
  return(beta_est)
}

# Variable selection by Frequentist Backpropagation (FBP)
variable_selection_FBP <- function(X, y, beta_post, lambda_post, intercept = TRUE, IT = 5, ITER = 5) {
  lambda_median <- median(lambda_post)
  
  # Run the frequentist rLASSO S5 algorithm with the median lambda
  s5_fit <- rrLASSO.S5(X, y, lam = lambda_median, intercept = intercept, IT = IT, ITER = ITER)
  
  # Keep only the variables selected by S5 (active set)
  beta_est <- colMeans(beta_post)
  selected_vars <- s5_fit$hppm
  
  active_beta <- rep(0, ncol(X))
  active_beta[selected_vars] <- beta_est[selected_vars]
  names(active_beta) <- colnames(X)
  return(active_beta)
}
