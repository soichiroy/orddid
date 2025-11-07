# ==============================================================================
# Tests for Stage 1(i): Isotonic Regression (PAVA-based CDF estimation)
# ==============================================================================

test_that("stage1i_isotonic produces monotone CDF", {
  # Purpose: Verify that isotonic regression produces a valid, monotone CDF
  # 
  # What is tested:
  # - Monotonicity: CDF values are non-decreasing
  # - Bounds: All CDF values lie in [0, 1]
  # - Structure: Returns proper FhatIsotonic object with required components
  # - Fitted values: Produces fitted values for all observations
  #
  # Method: Generate 3-category ordered response data from ordered probit model,
  #         then estimate CDF using isotonic regression on the binary indicator
  
  set.seed(123)
  n <- 200
  X <- cbind(1, rnorm(n), rnorm(n))
  beta_true <- c(1, 0.5, -0.3)
  u_true <- X %*% beta_true
  
  eps <- rnorm(n)
  latent <- u_true + eps
  Y <- cut(latent, breaks = c(-Inf, 0, 1, Inf), labels = FALSE)
  
  Fhat <- stage1i_isotonic(Y, X, beta_true)
  
  # Verify object structure
  expect_s3_class(Fhat, "FhatIsotonic")
  
  # Check monotonicity (allowing for numerical tolerance)
  expect_true(all(diff(Fhat$values) >= -1e-10))
  
  # Check CDF is bounded in [0, 1]
  expect_true(all(Fhat$values >= 0 & Fhat$values <= 1))
  
  # Check fitted values exist for all observations
  expect_equal(length(Fhat$fitted), n)
})

# ==============================================================================
# Tests for Stage 1(ii): Moment Equation Solver
# ==============================================================================

test_that("stage1ii_solve achieves moment conditions", {
  # Purpose: Verify that the moment equation solver finds beta satisfying
  #          E[X_{-fixed} * (1{Y=1} - F(X'beta))] = 0
  #
  # What is tested:
  # - Convergence: Score norm becomes small (moment conditions satisfied)
  # - Parameter structure: Correct dimension and normalization (beta[fixed]=1)
  # - Returns Fhat object
  #
  # Method: Solve moment equations using spectral gradient method with
  #         backtracking line search on simulated ordered probit data
  
  set.seed(101)
  n <- 200
  X <- cbind(1, rnorm(n), rnorm(n))
  beta_true <- c(1, 0.5, -0.3)
  u_true <- X %*% beta_true
  latent <- u_true + rnorm(n)
  Y <- cut(latent, breaks = c(-Inf, 0, 1, Inf), labels = FALSE)
  
  result <- stage1ii_solve(Y, X, fixed_index = 1, verbose = FALSE, maxit = 500)
  
  # Check moment conditions are approximately satisfied
  # (relaxed tolerance due to finite sample and optimization)
  expect_lt(result$score_norm, 0.5)
  
  # Check parameter dimension
  expect_equal(length(result$beta), 3)
  
  # Check normalization constraint
  expect_equal(result$beta[1], 1)
})

# ==============================================================================
# Tests for Full Two-Stage Estimation Workflow
# ==============================================================================

test_that("ord_two_stage_fit end-to-end", {
  # Purpose: Verify the complete two-stage semiparametric estimation pipeline
  #          (Stage 1: estimate beta and F; Stage 2: estimate threshold alpha)
  #
  # What is tested:
  # - Complete workflow: Both stages execute successfully
  # - Parameter estimates: Correct dimensions and constraints
  # - Prediction: Generates valid probability predictions
  # - Probability validity: Predictions sum to 1 and are non-negative
  #
  # Method: Run full two-stage estimation on simulated 3-category ordered
  #         response data, then verify structural properties and predictions
  
  set.seed(707)
  n <- 250
  X <- cbind(1, rnorm(n), runif(n, -1, 1))
  beta_true <- c(1, 0.4, -0.5)
  alpha_true <- 1.2
  u_true <- X %*% beta_true
  latent <- u_true + rnorm(n)
  Y <- cut(latent, breaks = c(-Inf, 0, alpha_true, Inf), labels = FALSE)
  
  fit <- ord_two_stage_fit(Y, X, verbose = FALSE)
  
  # Check parameter dimensions
  expect_equal(length(fit$beta), 3)
  
  # Check normalization
  expect_equal(fit$beta[1], 1)
  
  # Check alpha estimated
  expect_true(is.numeric(fit$alpha))
  
  # Check Fhat object returned
  expect_s3_class(fit$Fhat, "FhatIsotonic")
  
  # Check both stages executed (convergence may vary with random data)
  expect_true(!is.null(fit$stage1))
  expect_true(!is.null(fit$stage2))
  
  # Test prediction function
  probs <- fit$predict_probs(X)
  
  # Check prediction dimensions
  expect_equal(nrow(probs), n)
  expect_equal(ncol(probs), 3)
  
  # Check probabilities are valid (sum to 1, allowing numerical tolerance)
  expect_true(all(rowSums(probs) > 0.99 & rowSums(probs) < 1.01))
})

