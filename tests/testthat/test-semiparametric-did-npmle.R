## Test semiparametric DID (NPMLE method) internals

test_that(".semiparametric_did_once runs without error on simulated data", {
  skip_if_not_installed("icenReg")
  skip_if_not_installed("splines2")
  skip_if_not_installed("quadprog")

  set.seed(42)
  n <- 400

  ## True parameters
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

  X00 <- X[idx_ctrl, ]; X01 <- X[idx_ctrl, ]
  X10 <- X[idx_trt, ];  X11 <- X[idx_trt, ]

  res <- .semiparametric_did_once(
    Y00, X00, Y01, X01, Y10, X10, Y11, X11,
    fixed_index = 2, verbose = FALSE
  )

  expect_true(is.list(res))
  expect_length(res$zeta, 3)
  expect_length(res$tau, 2)
  expect_true(res$tau[1] <= res$tau[2])
  expect_true(all(abs(res$zeta) < 1))
})

test_that(".smooth_Fhat produces valid monotone CDF", {
  skip_if_not_installed("splines2")
  skip_if_not_installed("quadprog")

  set.seed(123)
  n <- 500
  X <- cbind(1, rnorm(n), rnorm(n))
  beta <- c(0.5, -1, 0.3)
  kappa <- c(-Inf, 0, 1.0, Inf)
  Ystar <- as.numeric(X %*% beta) + rnorm(n)
  Y <- as.integer(rowSums(outer(Ystar, kappa[-1], ">")))

  b_00 <- .init_beta_probit(Y, X, fixed_index = 2)
  Fhat <- stage1i_isotonic(Y, X, b_00)
  sm <- .smooth_Fhat(Fhat, n0 = n)

  u_test <- seq(sm$u_range[1], sm$u_range[2], length.out = 100)
  F_vals <- sm$Ftilde(u_test)

  ## Monotone
  expect_true(all(diff(F_vals) >= -1e-10))
  ## Bounded in [0, 1]
  expect_true(all(F_vals >= 0 & F_vals <= 1))
})
