# Test for Log-Likelihood Function

test_that("Log-likelihood function computes correctly", {
  # Simulated data
  set.seed(123)
  Y <- c(1, 2, 2, 3, 3, 3)
  par <- c(0, log(1))  # mu = 0, sd = 1
  cutoff_val <- c(0, 1)

  # Expected log-likelihood (manual calculation)
  expected_ll <- sum(
    c(
      log(pnorm((0 - 0) / 1) - pnorm((-Inf - 0) / 1)),
      log(pnorm((1 - 0) / 1) - pnorm((0 - 0) / 1)),
      log(pnorm((1 - 0) / 1) - pnorm((0 - 0) / 1)),
      log(pnorm((Inf - 0) / 1) - pnorm((1 - 0) / 1)),
      log(pnorm((Inf - 0) / 1) - pnorm((1 - 0) / 1)),
      log(pnorm((Inf - 0) / 1) - pnorm((1 - 0) / 1))
    )
  )

  # Compute log-likelihood using the function
  computed_ll <- -.LogLikelihoodProbit(par, Y, cutoff_val)

  # Test if the computed log-likelihood matches the expected value
  expect_equal(computed_ll, expected_ll, tolerance = 1e-6)
})
