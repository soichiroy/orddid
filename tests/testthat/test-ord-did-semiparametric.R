## Test ord_did(method = "semiparametric") through user-facing API

test_that("ord_did with method='semiparametric' runs and returns correct structure", {
  skip_if_not_installed("icenReg")
  skip_if_not_installed("splines2")
  skip_if_not_installed("quadprog")

  set.seed(42)
  n <- 400
  n_half <- n / 2
  beta_00 <- c(0.5, -1, 0.3)
  kappa <- c(-Inf, 0, 1.0, Inf)

  gen_cell <- function(n_cell, beta, xi, shift = 0) {
    X <- cbind(1, rnorm(n_cell), rnorm(n_cell))
    mu <- as.numeric(X %*% beta) + shift
    sig <- exp(as.numeric(X %*% xi))
    Ystar <- mu + sig * rnorm(n_cell)
    Y <- as.integer(rowSums(outer(Ystar, kappa[-1], ">")))
    list(Y = Y, X1 = X[, 2], X2 = X[, 3])
  }

  c00 <- gen_cell(n_half, c(0.5, -1, 0.3), c(0, 0, 0))
  c01 <- gen_cell(n_half, c(0.3, -0.8, 0.2), c(0.1, 0.05, -0.05))
  c10 <- gen_cell(n_half, c(0.6, -1.2, 0.4), c(0.15, 0.1, -0.1))
  c11 <- gen_cell(n_half, c(0.6, -1.2, 0.4), c(0.15, 0.1, -0.1), shift = 0.3)

  dat <- data.frame(
    id   = c(1:n_half, 1:n_half, (n_half+1):n, (n_half+1):n),
    Y    = c(c00$Y, c01$Y, c10$Y, c11$Y),
    post = c(rep(FALSE, n_half), rep(TRUE, n_half),
             rep(FALSE, n_half), rep(TRUE, n_half)),
    treat = c(rep(FALSE, 2*n_half), rep(TRUE, 2*n_half)),
    X1   = c(c00$X1, c01$X1, c10$X1, c11$X1),
    X2   = c(c00$X2, c01$X2, c10$X2, c11$X2)
  )

  result <- ord_did(
    data = dat,
    outcome = "Y",
    post = "post",
    treat = "treat",
    cluster = "id",
    n_boot = 20,
    method = "semiparametric",
    covariates = c("X1", "X2"),
    fixed_index = 2
  )

  expect_true(is.list(result))
  expect_true("estimate_effects" %in% names(result))
  expect_true("relative_effects" %in% names(result))
  expect_equal(nrow(result$estimate_effects), 3)
  expect_equal(result$inputs$method, "semiparametric")
})

test_that("ord_did semiparametric errors without covariates", {
  dat <- data.frame(id = 1:10, Y = sample(0:2, 10, TRUE),
                    post = rep(c(FALSE, TRUE), 5),
                    treat = rep(c(FALSE, TRUE), each = 5))

  expect_error(
    ord_did(data = dat, outcome = "Y", post = "post", treat = "treat",
            cluster = "id", n_boot = 5, method = "semiparametric"),
    "covariates must be specified"
  )
})
