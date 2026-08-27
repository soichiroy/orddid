.skip_semiparametric_stack <- function() {
  testthat::skip_if_not_installed("icenReg")
  testthat::skip_if_not_installed("splines2")
  testthat::skip_if_not_installed("quadprog")
}

.draw_three_category_outcome <- function(X, beta, alpha = 1.1,
                                         shift = 0, xi = NULL) {
  if (is.null(xi)) xi <- rep(0, ncol(X))
  latent <- as.numeric(X %*% beta) + shift +
    exp(as.numeric(X %*% xi)) * stats::rnorm(nrow(X))
  as.integer(cut(
    latent,
    breaks = c(-Inf, 0, alpha, Inf),
    labels = FALSE
  )) - 1L
}

.make_semiparametric_first_stage_data <- function(
    n = 180L, seed = 541L, first_slope = -1) {
  set.seed(seed)
  X <- cbind(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  Y <- .draw_three_category_outcome(
    X,
    beta = c(first_slope, 0.3),
    alpha = 1.1
  )
  if (!identical(sort(unique(Y)), 0:2)) {
    Y[seq_len(3L)] <- 0:2
  }
  list(Y = Y, X = X)
}

.make_semiparametric_panel <- function(n_per_group = 80L, seed = 922L) {
  set.seed(seed)
  X_control <- cbind(
    x1 = stats::rnorm(n_per_group),
    x2 = stats::rnorm(n_per_group)
  )
  X_treated <- cbind(
    x1 = stats::rnorm(n_per_group),
    x2 = stats::rnorm(n_per_group)
  )

  beta_00 <- c(-1, 0.3)
  beta_01 <- c(-0.85, 0.2)
  beta_10 <- c(-1.1, 0.4)
  xi_01 <- c(0.05, -0.03)
  xi_10 <- c(-0.04, 0.04)

  Y00 <- .draw_three_category_outcome(X_control, beta_00)
  Y01 <- .draw_three_category_outcome(X_control, beta_01, xi = xi_01)
  Y10 <- .draw_three_category_outcome(X_treated, beta_10, xi = xi_10)
  Y11 <- .draw_three_category_outcome(
    X_treated,
    beta = beta_10,
    shift = 0.25,
    xi = xi_10 + xi_01
  )
  if (!identical(sort(unique(Y00)), 0:2)) {
    Y00[seq_len(3L)] <- 0:2
  }

  data.frame(
    id = c(
      seq_len(n_per_group),
      seq_len(n_per_group),
      n_per_group + seq_len(n_per_group),
      n_per_group + seq_len(n_per_group)
    ),
    Y = c(Y00, Y01, Y10, Y11),
    post = rep(c(FALSE, TRUE, FALSE, TRUE), each = n_per_group),
    treat = rep(c(FALSE, FALSE, TRUE, TRUE), each = n_per_group),
    X1 = c(X_control[, 1L], X_control[, 1L],
           X_treated[, 1L], X_treated[, 1L]),
    X2 = c(X_control[, 2L], X_control[, 2L],
           X_treated[, 2L], X_treated[, 2L])
  )
}

.fit_semiparametric_test_panel <- function(
    data, n_boot = 0L, semiparametric_regularization = "ispline") {
  ord_did(
    data = data,
    outcome = "Y",
    post = "post",
    treat = "treat",
    cluster = "id",
    n_boot = n_boot,
    method = "semiparametric",
    covariates = c("X1", "X2"),
    fixed_index = 1L,
    semiparametric_regularization = semiparametric_regularization
  )
}
