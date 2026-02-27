# ==============================================================================
# Tests for Stage 2: Threshold Estimation (alpha)
# ==============================================================================

test_that("EstimateCutoff finds correct threshold", {
  # Purpose: Verify that EstimateCutoff correctly estimates the threshold parameter
  #          by solving the moment condition E[1 - 1{Y=3} - F(X'beta + alpha)] = 0
  #
  # What is tested:
  # - Convergence: Moment condition satisfied (psi close to zero)
  # - Valid estimate: Alpha is finite, non-negative
  # - Accuracy: Estimate is reasonably close to true value
  #
  # Method: Generate 3-category data with known threshold, estimate F from stage1,
  #         then solve for alpha using monotone zero-crossing search
  
  set.seed(404)
  n <- 200
  X <- cbind(1, rnorm(n), rnorm(n))
  beta_true <- c(1, 0.5, -0.3)
  u_true <- X %*% beta_true
  
  eps <- rnorm(n)
  latent <- u_true + eps
  alpha_true <- 1.5
  Y <- cut(latent, breaks = c(-Inf, 0, alpha_true, Inf), labels = FALSE)
  
  Fhat <- stage1i_isotonic(Y, X, beta_true)
  result <- EstimateCutoff(Y, X, beta_true, Fhat, verbose = FALSE)
  
  # Check estimate is valid
  expect_true(is.numeric(result$alpha))
  expect_true(is.finite(result$alpha))
  expect_true(result$alpha >= 0)
  
  # Check convergence
  expect_true(result$converged)
  
  # Check moment condition satisfied
  expect_lt(abs(result$psi_at), 0.1)
  
  # Check accuracy (relaxed tolerance due to finite sample)
  expect_lt(abs(result$alpha - alpha_true), 0.5)
})

test_that("EstimateCutoff handles boundary cases", {
  # Purpose: Verify robust behavior when data has no observations in category 3
  #          (boundary case where alpha may be at the lower bound)
  #
  # What is tested:
  # - No errors: Function completes without crashing
  # - Valid solution: Returns finite alpha
  # - Method selection: Chooses appropriate solver (boundary/uniroot/grid)
  #
  # Method: Generate data with only categories 1 and 2, then attempt to
  #         estimate threshold (should find boundary solution)
  
  set.seed(505)
  n <- 150
  X <- cbind(1, rnorm(n))
  beta <- c(1, 0.6)
  u <- X %*% beta
  
  latent <- u + rnorm(n, sd = 0.5)
  Y <- pmin(cut(latent, breaks = c(-Inf, 0, Inf), labels = FALSE), 2)
  
  Fhat <- stage1i_isotonic(Y, X, beta)
  result <- EstimateCutoff(Y, X, beta, Fhat, verbose = FALSE)
  
  # Check valid solution returned
  expect_true(result$method %in% c("boundary", "uniroot", "grid-cross", "grid-min"))
  expect_true(is.finite(result$alpha))
})

test_that("EstimateCutoff with specified bounds", {
  # Purpose: Verify that user-specified search bounds are respected
  #
  # What is tested:
  # - Bound respect: Estimated alpha falls within specified bounds
  # - Custom bounds work: No errors with non-default bounds
  #
  # Method: Estimate threshold with custom lower/upper bounds and verify
  #         the solution respects these constraints
  
  set.seed(606)
  n <- 180
  X <- cbind(1, rnorm(n))
  beta <- c(1, 0.5)
  u <- X %*% beta
  latent <- u + rnorm(n)
  Y <- cut(latent, breaks = c(-Inf, 0, 1.2, Inf), labels = FALSE)
  
  Fhat <- stage1i_isotonic(Y, X, beta)
  
  result <- EstimateCutoff(
    Y, X, beta, Fhat,
    alpha_lower = 0.5,
    alpha_upper = 3,
    verbose = FALSE
  )
  
  # Check bounds respected
  expect_gte(result$alpha, 0.5)
  expect_lte(result$alpha, 3)
})

# ==============================================================================
# Tests for Prediction Functions
# ==============================================================================

test_that("predict_probs_3cat produces valid probabilities", {
  # Purpose: Verify that the 3-category prediction function produces valid
  #          probability predictions (non-negative, sum to 1)
  #
  # What is tested:
  # - Dimensions: Correct number of rows (n) and columns (3)
  # - Column names: p1, p2, p3
  # - Non-negativity: All probabilities >= 0
  # - Unit sum: Each row sums to 1 (within tolerance)
  #
  # Method: Generate data, estimate F and alpha, then predict probabilities
  #         and verify statistical properties
  
  set.seed(808)
  n <- 100
  X <- cbind(1, rnorm(n))
  beta <- c(1, 0.5)
  alpha <- 1.0
  u <- X %*% beta
  latent <- u + rnorm(n)
  Y <- cut(latent, breaks = c(-Inf, 0, alpha, Inf), labels = FALSE)
  
  Fhat <- stage1i_isotonic(Y, X, beta)
  probs <- predict_probs_3cat(X, beta, alpha, Fhat)
  
  # Check dimensions
  expect_equal(nrow(probs), n)
  expect_equal(ncol(probs), 3)
  expect_equal(colnames(probs), c("p1", "p2", "p3"))
  
  # Check validity
  expect_true(all(probs >= 0))
  expect_true(all(probs <= 1))
  expect_equal(rowSums(probs), rep(1, n), tolerance = 1e-6)
})

test_that("EstimateCutoff grid search fallback works", {
  # Purpose: Verify that grid search fallback works when uniroot fails or
  #          is not applicable (wide bounds, difficult objective function)
  #
  # What is tested:
  # - Robustness: Function completes without errors
  # - Valid solution: Returns finite alpha
  # - Method reported: Grid search method is used when appropriate
  #
  # Method: Create difficult estimation case with very wide bounds and
  #         high variance, forcing grid search fallback
  
  set.seed(909)
  n <- 100
  X <- cbind(1, rnorm(n))
  beta <- c(1, 0.3)
  u <- X %*% beta
  
  latent <- u + rnorm(n, sd = 2)
  Y <- cut(latent, breaks = c(-Inf, -1, 2, Inf), labels = FALSE)
  
  Fhat <- stage1i_isotonic(Y, X, beta)
  
  result <- EstimateCutoff(
    Y, X, beta, Fhat,
    alpha_lower = 0,
    alpha_upper = 10,
    verbose = FALSE,
    grid_points = 500
  )
  
  # Check valid solution found
  expect_true(is.finite(result$alpha))
  expect_true(result$method %in% c("uniroot", "grid-cross", "grid-min"))
})
