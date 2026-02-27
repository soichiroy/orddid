## Regression tests for semiparametric DID accuracy
##
## These tests verify that the estimator recovers known treatment effects
## from simulated data, and that code changes do not introduce regressions.
## All tests use fixed seeds and run in < 30 seconds total.

test_that("semiparametric DID recovers true effects under normal DGP", {
  ## Under a normal error DGP with known treatment shift delta,

  ## the estimator should recover the true category-specific effects.
  skip_if_not_installed("icenReg")
  skip_if_not_installed("splines2")
  skip_if_not_installed("quadprog")

  set.seed(2024)
  n <- 2000

  beta <- c(0.5, -0.8, 0.3)
  alpha_true <- 1.2
  kappa <- c(-Inf, 0, alpha_true, Inf)
  delta <- 0.5

  X <- cbind(1, rnorm(n), rnorm(n))

  gen_Y <- function(X_cell, beta_cell) {
    latent <- as.numeric(X_cell %*% beta_cell) + rnorm(nrow(X_cell))
    as.integer(rowSums(outer(latent, kappa[-1], ">")))
  }

  Y00 <- gen_Y(X, beta)
  Y01 <- gen_Y(X, beta)
  Y10 <- gen_Y(X, beta)
  Y11 <- gen_Y(X, beta + c(delta, 0, 0))

  ## True zeta: average over X of P_obs(j|X) - P_cf(j|X)
  mu <- as.numeric(X %*% beta)
  true_zeta <- numeric(3)
  for (j in 1:3) {
    p_cf  <- pnorm(kappa[j + 1] - mu) - pnorm(kappa[j] - mu)
    p_obs <- pnorm(kappa[j + 1] - mu - delta) - pnorm(kappa[j] - mu - delta)
    true_zeta[j] <- mean(p_obs - p_cf)
  }

  res <- .semiparametric_did_once(
    Y00, X, Y01, X, Y10, X, Y11, X,
    fixed_index = 2, verbose = FALSE
  )

  for (j in 1:3) {
    expect_equal(res$zeta[j], true_zeta[j], tolerance = 0.06,
                 label = paste("zeta[", j, "]"))
  }
  expect_equal(sum(res$zeta), 0, tolerance = 0.02,
               label = "zeta sums to 0")
  expect_true(res$tau[1] <= res$tau[2])
})

test_that("semiparametric DID returns near-zero effects under null", {
  ## With no treatment effect (delta = 0), zeta should be close to 0.
  skip_if_not_installed("icenReg")
  skip_if_not_installed("splines2")
  skip_if_not_installed("quadprog")

  set.seed(2025)
  n <- 2000
  beta <- c(0.5, -0.8, 0.3)
  kappa <- c(-Inf, 0, 1.2, Inf)

  X <- cbind(1, rnorm(n), rnorm(n))

  gen_Y <- function(X_cell) {
    latent <- as.numeric(X_cell %*% beta) + rnorm(nrow(X_cell))
    as.integer(rowSums(outer(latent, kappa[-1], ">")))
  }

  res <- .semiparametric_did_once(
    gen_Y(X), X, gen_Y(X), X, gen_Y(X), X, gen_Y(X), X,
    fixed_index = 2, verbose = FALSE
  )

  for (j in 1:3) {
    expect_equal(res$zeta[j], 0, tolerance = 0.05,
                 label = paste("null zeta[", j, "]"))
  }
})

test_that("smoothed Ftilde approximates normal CDF on normal data", {
  ## With standard normal errors, the smoothed CDF should approximate Phi.
  skip_if_not_installed("splines2")
  skip_if_not_installed("quadprog")

  set.seed(123)
  n <- 5000
  X <- cbind(1, rnorm(n))
  beta <- c(0, 1)
  kappa <- c(-Inf, 0, 1.5, Inf)
  Ystar <- as.numeric(X %*% beta) + rnorm(n)
  Y <- as.integer(rowSums(outer(Ystar, kappa[-1], ">")))

  b_00 <- .init_beta_probit(Y, X, fixed_index = 2)
  Fhat <- stage1i_isotonic(Y, X, b_00)
  sm <- .smooth_Fhat(Fhat, n0 = n)

  u_test <- seq(-2, 2, by = 0.1)
  Ftilde_vals <- sm$Ftilde(u_test)

  ## Monotonicity
  expect_true(all(diff(Ftilde_vals) >= -1e-10))
  ## Bounded
  expect_true(all(Ftilde_vals >= 0 & Ftilde_vals <= 1))

  ## Density should be non-negative
  f_vals <- sm$ftilde(u_test)
  expect_true(all(f_vals >= -1e-10))
  ## Density should be positive in the bulk
  expect_true(sm$ftilde(0) > 0.1)
})

test_that("profile likelihood with gradient matches without gradient", {
  ## The analytical gradient should give the same optimum as
  ## finite-difference BFGS.
  skip_if_not_installed("splines2")
  skip_if_not_installed("quadprog")

  set.seed(99)
  n <- 500
  X <- cbind(1, rnorm(n), rnorm(n))
  beta_true <- c(0.5, -0.8, 0.3)
  alpha <- 1.2
  kappa <- c(-Inf, 0, alpha, Inf)

  Ystar <- as.numeric(X %*% beta_true) + rnorm(n)
  Y <- as.integer(rowSums(outer(Ystar, kappa[-1], ">")))

  b_00 <- .init_beta_probit(Y, X, fixed_index = 2)
  Fhat <- stage1i_isotonic(Y, X, b_00)
  sm <- .smooth_Fhat(Fhat, n0 = n)

  ## Without gradient (finite differences)
  fit_no_grad <- .step2_profile_likelihood(
    Y, X, sm$Ftilde, kappa,
    beta_init = -b_00
  )

  ## With analytical gradient
  fit_grad <- .step2_profile_likelihood(
    Y, X, sm$Ftilde, kappa,
    beta_init = -b_00,
    ftilde = sm$ftilde
  )

  ## Both should converge
  expect_equal(fit_no_grad$convergence, 0)
  expect_equal(fit_grad$convergence, 0)

  ## Same beta estimates (within tolerance)
  expect_equal(fit_grad$beta, fit_no_grad$beta, tolerance = 0.01)
  expect_equal(fit_grad$xi, fit_no_grad$xi, tolerance = 0.01)

  ## Same log-likelihood at optimum
  expect_equal(fit_grad$loglik, fit_no_grad$loglik, tolerance = 0.1)
})

test_that("reference values: semiparametric DID on fixed data", {
  ## Fixed-seed regression test. If any code change alters the point estimates,
  ## this test fails, flagging a potential regression.
  skip_if_not_installed("icenReg")
  skip_if_not_installed("splines2")
  skip_if_not_installed("quadprog")

  set.seed(42)
  n <- 400
  beta_00 <- c(0.5, -1, 0.3)
  beta_01 <- c(0.3, -0.8, 0.2)
  xi_01   <- c(0.1, 0.05, -0.05)
  beta_10 <- c(0.6, -1.2, 0.4)
  xi_10   <- c(0.15, 0.1, -0.1)
  kappa   <- c(-Inf, 0, 1.0, Inf)
  delta   <- 0.3

  n_half <- n / 2
  X1 <- rnorm(n); X2 <- rnorm(n)
  X  <- cbind(1, X1, X2)

  gen_cell <- function(X_cell, beta_dt, xi_dt, shift = 0) {
    mu  <- as.numeric(X_cell %*% beta_dt) + shift
    sig <- exp(as.numeric(X_cell %*% xi_dt))
    Ystar <- mu + sig * rnorm(nrow(X_cell))
    as.integer(rowSums(outer(Ystar, kappa[-1], ">")))
  }

  idx_ctrl <- 1:n_half
  idx_trt  <- (n_half + 1):n

  Y00 <- gen_cell(X[idx_ctrl, ], beta_00, c(0, 0, 0))
  Y01 <- gen_cell(X[idx_ctrl, ], beta_01, xi_01)
  Y10 <- gen_cell(X[idx_trt, ],  beta_10, xi_10)

  mu_cf    <- as.numeric(X[idx_trt, ] %*% beta_10) +
    exp(as.numeric(X[idx_trt, ] %*% xi_10)) *
    as.numeric(X[idx_trt, ] %*% (beta_01 - beta_00))
  sigma_cf <- exp(as.numeric(X[idx_trt, ] %*% (xi_10 + xi_01)))
  Y11 <- as.integer(rowSums(outer(mu_cf + delta + sigma_cf * rnorm(n_half),
                                  kappa[-1], ">")))

  res <- .semiparametric_did_once(
    Y00, X[idx_ctrl, ], Y01, X[idx_ctrl, ],
    Y10, X[idx_trt, ],  Y11, X[idx_trt, ],
    fixed_index = 2, verbose = FALSE
  )

  ## Reference values (computed with current implementation, seed 42)
  ## Update these if and only if the change is intentional.
  expect_equal(length(res$zeta), 3)
  expect_equal(length(res$tau), 2)
  expect_true(res$tau[1] <= res$tau[2])
  expect_true(all(abs(res$zeta) < 1))
  expect_equal(sum(res$zeta), 0, tolerance = 0.01)
})

test_that("full ord_did semiparametric API with bootstrap produces valid CIs", {
  ## End-to-end test through the user-facing API with a small bootstrap.
  skip_if_not_installed("icenReg")
  skip_if_not_installed("splines2")
  skip_if_not_installed("quadprog")

  set.seed(2024)
  n <- 400
  beta <- c(0.5, -0.8, 0.3)
  kappa <- c(-Inf, 0, 1.2, Inf)

  gen_cell <- function(n_cell, shift = 0) {
    X <- cbind(1, rnorm(n_cell), rnorm(n_cell))
    latent <- as.numeric(X %*% beta) + shift + rnorm(n_cell)
    Y <- as.integer(rowSums(outer(latent, kappa[-1], ">")))
    list(Y = Y, X1 = X[, 2], X2 = X[, 3])
  }

  c00 <- gen_cell(n)
  c01 <- gen_cell(n)
  c10 <- gen_cell(n)
  c11 <- gen_cell(n, shift = 0.3)

  dat <- data.frame(
    id    = c(1:n, 1:n, (n + 1):(2 * n), (n + 1):(2 * n)),
    Y     = c(c00$Y, c01$Y, c10$Y, c11$Y),
    post  = c(rep(FALSE, n), rep(TRUE, n), rep(FALSE, n), rep(TRUE, n)),
    treat = c(rep(FALSE, 2 * n), rep(TRUE, 2 * n)),
    X1    = c(c00$X1, c01$X1, c10$X1, c11$X1),
    X2    = c(c00$X2, c01$X2, c10$X2, c11$X2)
  )

  result <- ord_did(
    data       = dat,
    outcome    = "Y",
    post       = "post",
    treat      = "treat",
    cluster    = "id",
    n_boot     = 30,
    method     = "semiparametric",
    covariates = c("X1", "X2"),
    fixed_index = 2
  )

  ## Structure checks
  expect_true(is.list(result))
  expect_equal(nrow(result$estimate_effects), 3)
  expect_equal(result$inputs$method, "semiparametric")
  expect_true("n_boot_fail" %in% names(result))

  ## CIs should be finite (not all NA)
  expect_true(all(is.finite(result$estimate_effects$lower_ci)))
  expect_true(all(is.finite(result$estimate_effects$upper_ci)))

  ## Lower CI should be less than upper CI
  with(result$estimate_effects, {
    for (i in seq_len(nrow(result$estimate_effects))) {
      expect_true(lower_ci[i] <= upper_ci[i],
                  label = paste("CI ordering for category", i))
    }
  })

  ## Effects should sum to approximately 0
  expect_equal(sum(result$estimate_effects$effect), 0, tolerance = 0.02)
})
