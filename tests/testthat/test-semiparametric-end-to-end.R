.make_appendix_semiparametric_dgp <- function(
    seed, delta, n = 1200L, alpha = 1.2) {
  set.seed(seed)
  X <- cbind(
    intercept = 1,
    normalized = stats::rnorm(n),
    x2 = stats::rnorm(n)
  )
  beta <- c(0.5, -0.8, 0.3)
  kappa <- c(-Inf, 0, alpha, Inf)
  draw_outcome <- function(shift = 0) {
    latent <- as.numeric(X %*% beta) + shift + stats::rnorm(n)
    as.integer(cut(latent, breaks = kappa, labels = FALSE)) - 1L
  }

  Y00 <- draw_outcome()
  Y01 <- draw_outcome()
  Y10 <- draw_outcome()
  Y11 <- draw_outcome(delta)

  mu <- as.numeric(X %*% beta)
  true_zeta <- vapply(seq_len(3L), function(j) {
    observed <- stats::pnorm(kappa[j + 1L] - mu - delta) -
      stats::pnorm(kappa[j] - mu - delta)
    counterfactual <- stats::pnorm(kappa[j + 1L] - mu) -
      stats::pnorm(kappa[j] - mu)
    mean(observed - counterfactual)
  }, numeric(1))

  list(
    X = X,
    Y00 = Y00,
    Y01 = Y01,
    Y10 = Y10,
    Y11 = Y11,
    true_zeta = true_zeta
  )
}

.fit_appendix_semiparametric_dgp <- function(data) {
  .semiparametric_did_once(
    data$Y00, data$X,
    data$Y01, data$X,
    data$Y10, data$X,
    data$Y11, data$X,
    fixed_index = 2L,
    semiparametric_regularization = "ispline"
  )
}

test_that("canonical estimator recovers a known category-effect shift", {
  .skip_semiparametric_stack()
  data <- .make_appendix_semiparametric_dgp(seed = 8001L, delta = 0.5)
  fit <- .fit_appendix_semiparametric_dgp(data)
  estimation_error <- fit$zeta - data$true_zeta

  expect_equal(unname(fit$convergence), c(0L, 0L, 0L))
  expect_equal(fit$first_stage$fixed_index, 2L)
  expect_equal(unname(fit$first_stage$constant_index), 1L)
  expect_equal(fit$first_stage$b[2L], 1)
  expect_setequal(
    unique(fit$profile_fits$control_post$candidate_fits$optimizer),
    c("BFGS", "Nelder-Mead")
  )
  expect_setequal(
    unique(fit$profile_fits$treated_pre$candidate_fits$optimizer),
    c("BFGS", "Nelder-Mead")
  )

  expect_lt(max(abs(estimation_error)), 0.05)
  expect_lt(fit$zeta[1L], 0)
  expect_gt(fit$zeta[3L], 0)
  expect_equal(sum(fit$zeta), 0, tolerance = 1e-8)
})

test_that("canonical estimator returns near-zero effects under the null", {
  .skip_semiparametric_stack()
  data <- .make_appendix_semiparametric_dgp(seed = 8002L, delta = 0)
  fit <- .fit_appendix_semiparametric_dgp(data)

  expect_equal(unname(fit$convergence), c(0L, 0L, 0L))
  expect_lt(max(abs(fit$zeta)), 0.06)
  expect_equal(sum(fit$zeta), 0, tolerance = 1e-8)
})
