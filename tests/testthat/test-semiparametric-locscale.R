# ==============================================================================
# Tests for Location-Scale Extension (Heteroskedastic Model)
# ==============================================================================

test_that("ord_locscale_fit converges correctly", {
  # Purpose: Verify that the location-scale model (simultaneous mean/scale estimation)
  #          converges and produces valid parameter estimates
  #
  # What is tested:
  # - Convergence: Moment conditions satisfied (small score norm)
  # - Parameter structure: Correct dimensions for beta (mean) and gamma (log-scale)
  # - Normalization: beta[fixed] = 1
  # - Positivity: All sigma (scale) values are positive
  #
  # Method: Generate heteroskedastic ordered response data where variance depends
  #         on covariates, then estimate location-scale model and verify properties
  
  set.seed(1010)
  n <- 200
  X <- cbind(1, rnorm(n), runif(n, -1, 1))
  
  beta_true <- c(1, 0.5, -0.3)
  gamma_true <- c(0, 0.2)
  alpha_true <- 1.3
  
  W <- X[, 1:2]
  mu <- X %*% beta_true
  log_sigma <- W %*% gamma_true
  sigma <- exp(log_sigma)
  
  latent <- mu + sigma * rnorm(n)
  Y <- cut(latent, breaks = c(-Inf, 0, alpha_true, Inf), labels = FALSE)
  
  # Get Fhat from homoskedastic first stage
  stage1 <- ord_two_stage_fit(Y, X, verbose = FALSE)
  Fhat <- stage1$Fhat
  alpha <- stage1$alpha
  
  # Fit heteroskedastic model
  fit <- ord_locscale_fit(Y, X, Fhat, alpha, verbose = FALSE)
  
  # Check convergence
  expect_true(fit$converged)
  expect_lt(fit$score_norm, 1e-3)
  
  # Check parameter dimensions
  expect_equal(length(fit$beta), 3)
  expect_equal(fit$beta[1], 1)
  expect_equal(length(fit$gamma), 2)
  expect_equal(length(fit$mu), n)
  expect_equal(length(fit$sigma), n)
  expect_equal(length(fit$S), n)
  
  # Check scale is positive
  expect_true(all(fit$sigma > 0))
})

test_that("ord_locscale_fit prediction function works", {
  # Purpose: Verify that the prediction function for location-scale model
  #          produces valid probability predictions
  #
  # What is tested:
  # - Function exists and is callable
  # - Correct dimensions (n x 3)
  # - Column names: p1, p2, p3
  # - Probability validity: non-negative, sum to 1
  # - Out-of-sample prediction: Works on new data
  #
  # Method: Fit location-scale model, then predict on both training and
  #         new test data, verifying prediction validity
  
  set.seed(1111)
  n <- 150
  X <- cbind(1, rnorm(n))
  
  beta <- c(1, 0.6)
  gamma <- c(0, 0.1)
  alpha <- 1.0
  
  W <- X
  mu <- X %*% beta
  sigma <- exp(W %*% gamma)
  latent <- mu + sigma * rnorm(n)
  Y <- cut(latent, breaks = c(-Inf, 0, alpha, Inf), labels = FALSE)
  
  Fhat <- stage1i_isotonic(Y, X, beta)
  fit <- ord_locscale_fit(Y, X, Fhat, alpha, verbose = FALSE)
  
  # Check prediction function exists
  expect_true(is.function(fit$predict_probs))
  
  # In-sample predictions
  probs <- fit$predict_probs(X, X)
  
  expect_equal(nrow(probs), n)
  expect_equal(ncol(probs), 3)
  expect_equal(colnames(probs), c("p1", "p2", "p3"))
  
  expect_true(all(probs >= 0))
  expect_true(all(probs <= 1))
  expect_equal(rowSums(probs), rep(1, n), tolerance = 1e-6)
  
  # Out-of-sample predictions
  X_new <- cbind(1, rnorm(50))
  probs_new <- fit$predict_probs(X_new, X_new)
  expect_equal(nrow(probs_new), 50)
  expect_true(all(rowSums(probs_new) >= 0.99))
})

test_that("ord_locscale_fit with initial values", {
  # Purpose: Verify that custom initial values are accepted and used correctly
  #
  # What is tested:
  # - Accepts user-provided beta_init and gamma_init
  # - Converges with custom starting values
  # - Returns correct parameter dimensions
  #
  # Method: Provide custom initial values and verify successful estimation
  
  set.seed(1212)
  n <- 120
  X <- cbind(1, rnorm(n))
  
  beta_init <- c(1, 0.5)
  gamma_init <- c(0, 0)
  alpha <- 1.2
  
  mu <- X %*% beta_init
  latent <- mu + rnorm(n)
  Y <- cut(latent, breaks = c(-Inf, 0, alpha, Inf), labels = FALSE)
  
  Fhat <- stage1i_isotonic(Y, X, beta_init)
  
  fit <- ord_locscale_fit(
    Y, X, Fhat, alpha,
    beta_init = beta_init,
    gamma_init = gamma_init,
    verbose = FALSE
  )
  
  expect_true(fit$converged)
  expect_equal(length(fit$beta), 2)
  expect_equal(length(fit$gamma), 2)
})

test_that("ord_locscale_fit handles weights", {
  # Purpose: Verify that observation weights are correctly incorporated
  #
  # What is tested:
  # - Equal weights match unweighted results
  # - Different weights produce different estimates
  # - No errors with weighted estimation
  #
  # Method: Compare weighted and unweighted estimation, verify weights matter
  
  set.seed(1313)
  n <- 150
  X <- cbind(1, rnorm(n))
  
  beta <- c(1, 0.4)
  alpha <- 1.0
  mu <- X %*% beta
  latent <- mu + rnorm(n)
  Y <- cut(latent, breaks = c(-Inf, 0, alpha, Inf), labels = FALSE)
  
  Fhat <- stage1i_isotonic(Y, X, beta)
  
  # Equal weights should match unweighted
  fit_unweighted <- ord_locscale_fit(Y, X, Fhat, alpha, verbose = FALSE)
  fit_weighted <- ord_locscale_fit(Y, X, Fhat, alpha, w = rep(1, n), verbose = FALSE)
  
  expect_equal(fit_unweighted$beta, fit_weighted$beta, tolerance = 1e-6)
  expect_equal(fit_unweighted$gamma, fit_weighted$gamma, tolerance = 1e-6)
  
  # Non-uniform weights should differ
  w <- runif(n, 0.5, 1.5)
  fit_w <- ord_locscale_fit(Y, X, Fhat, alpha, w = w, verbose = FALSE)
  
  expect_false(all(abs(fit_w$beta - fit_unweighted$beta) < 1e-4))
})

test_that("predict_probs_locscale convenience wrapper", {
  # Purpose: Verify that the convenience wrapper function for prediction
  #          produces identical results to the fit object's built-in predictor
  #
  # What is tested:
  # - Wrapper exists and is callable
  # - Produces identical results to fit$predict_probs()
  #
  # Method: Compare predictions from wrapper vs. built-in predictor
  
  set.seed(1414)
  n <- 100
  X <- cbind(1, rnorm(n))
  
  beta <- c(1, 0.5)
  gamma <- c(0, 0.15)
  alpha <- 1.1
  
  mu <- X %*% beta
  sigma <- exp(X %*% gamma)
  latent <- mu + sigma * rnorm(n)
  Y <- cut(latent, breaks = c(-Inf, 0, alpha, Inf), labels = FALSE)
  
  Fhat <- stage1i_isotonic(Y, X, beta)
  fit <- ord_locscale_fit(Y, X, Fhat, alpha, verbose = FALSE)
  
  # Compare wrapper vs. built-in predictor
  probs1 <- predict_probs_locscale(fit, X, X, Fhat, alpha)
  probs2 <- fit$predict_probs(X, X)
  
  expect_equal(probs1, probs2, tolerance = 1e-10)
})

test_that("ord_locscale_fit different solvers agree", {
  # Purpose: Verify that spectral and BB solvers produce similar results
  #          (tests robustness of estimation across different algorithms)
  #
  # What is tested:
  # - Both solvers converge
  # - Solutions are numerically close
  #
  # Method: Solve same problem with spectral and BB methods, compare results
  
  skip_if_not_installed("BB")
  
  set.seed(1515)
  n <- 120
  X <- cbind(1, rnorm(n))
  
  beta <- c(1, 0.4)
  gamma <- c(0, 0.1)
  alpha <- 1.0
  
  mu <- X %*% beta
  sigma <- exp(X %*% gamma)
  latent <- mu + sigma * rnorm(n)
  Y <- cut(latent, breaks = c(-Inf, 0, alpha, Inf), labels = FALSE)
  
  Fhat <- stage1i_isotonic(Y, X, beta)
  
  fit_spectral <- ord_locscale_fit(
    Y, X, Fhat, alpha,
    method = "spectral",
    verbose = FALSE,
    maxit = 150
  )
  
  fit_BB <- ord_locscale_fit(
    Y, X, Fhat, alpha,
    method = "BB",
    verbose = FALSE,
    maxit = 150
  )
  
  # Both should converge
  expect_true(fit_spectral$converged)
  expect_true(fit_BB$converged)
  
  # Solutions should be close (relaxed tolerance across solvers)
  expect_equal(fit_spectral$beta, fit_BB$beta, tolerance = 0.05)
  expect_equal(fit_spectral$gamma, fit_BB$gamma, tolerance = 0.05)
})

test_that("ord_locscale_fit moment conditions satisfied", {
  # Purpose: Verify that estimated parameters satisfy the moment conditions
  #          E[X_{-fixed} * (1{Y=1} - F(S))] = 0 and E[W * (1-1{Y=3} - F(S+α/σ))] = 0
  #
  # What is tested:
  # - Score vector is close to zero
  # - All score components are finite
  #
  # Method: Estimate parameters and check that moment conditions are satisfied
  #         to numerical tolerance
  
  set.seed(1616)
  n <- 180
  X <- cbind(1, rnorm(n), rnorm(n))
  
  beta <- c(1, 0.3, -0.2)
  gamma <- c(0, 0.1, 0)
  alpha <- 1.2
  
  mu <- X %*% beta
  sigma <- exp(X %*% gamma)
  latent <- mu + sigma * rnorm(n)
  Y <- cut(latent, breaks = c(-Inf, 0, alpha, Inf), labels = FALSE)
  
  Fhat <- stage1i_isotonic(Y, X, beta)
  
  fit <- ord_locscale_fit(Y, X, Fhat, alpha, verbose = FALSE, tol = 1e-5)
  
  # Check moment conditions
  expect_lt(max(abs(fit$score)), 1e-4)
  
  # Check all finite
  expect_true(all(is.finite(fit$score)))
})
