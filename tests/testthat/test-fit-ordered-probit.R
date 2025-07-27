# Test for FitOrderedProbit Function

test_that("FitOrderedProbit works correctly", {
  # Simulated data
  set.seed(123)
  n <- 50000
  E <- rnorm(n, mean = 0.5, sd = 1)
  outcome <- cut(
    E,
    breaks = c(-Inf, 0, 0.3, 1, Inf),
    labels = c("Low", "Medium", "High", "Very High")
  )
  Y <- as.numeric(outcome)

  # Test FitOrderedProbit without cutoff values
  result1 <- .FitOrderedProbit(Y = Y)
  # Check that the cutoff values are correctly estimated
  expect_equal(result1$cutoff, c(0, 0.3, 1), tolerance = 1e-4)
  expect_equal(result1$mu, 0.5, tolerance = 1e-4)

  # Test FitOrderedProbit with provided cutoff values
  result2 <- .FitOrderedProbit(Y = Y, cutoff_val = c(0, 0.3, 1))
  expect_equal(result2$mu, 0.5, tolerance = 1e-4)
})
