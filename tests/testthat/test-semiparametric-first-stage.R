test_that("ordered-probit initializer applies the paper normalization", {
  data <- .make_semiparametric_first_stage_data(n = 700L, seed = 110L)

  init <- .initialize_semiparametric_j3(
    data$Y,
    data$X,
    fixed_index = 1L
  )

  expect_equal(init$b[1L], 1)
  expect_equal(init$beta, -init$b)
  expect_gt(init$scale, 0)
  expect_gt(init$alpha, 0)
  expect_equal(
    init$alpha,
    init$ordered_probit$alpha / init$scale,
    tolerance = 1e-10
  )
  expect_equal(
    init$b,
    -init$ordered_probit$beta / init$scale,
    tolerance = 1e-10
  )
})

test_that("ordered-probit initializer rejects the wrong orientation", {
  data <- .make_semiparametric_first_stage_data(
    n = 900L,
    seed = 111L,
    first_slope = 1.5
  )

  expect_error(
    .initialize_semiparametric_j3(data$Y, data$X, fixed_index = 1L),
    "wrong or near-zero orientation"
  )
})

test_that("beta moment update freezes normalized and constant columns", {
  set.seed(112)
  n <- 1200L
  X <- cbind(
    normalized = stats::rnorm(n),
    intercept = 1,
    free = stats::rnorm(n)
  )
  b_true <- c(1, 0.25, -0.7)
  is_lowest <- stats::rbinom(
    n,
    size = 1L,
    prob = stats::pnorm(as.numeric(X %*% b_true))
  )
  Y <- ifelse(is_lowest == 1L, 0L, rep(c(1L, 2L), length.out = n))
  b_start <- c(1, 0.6, 0.9)

  update <- .update_beta_npmle_moment(
    Y = Y,
    X = X,
    b_start = b_start,
    Ftilde = stats::pnorm,
    ftilde = stats::dnorm,
    fixed_index = 1L
  )

  expect_equal(unname(update$constant_index), 2L)
  expect_equal(unname(update$free_index), 3L)
  expect_equal(update$b[c(1L, 2L)], b_start[c(1L, 2L)])
  expect_lt(update$objective, update$start_objective)
  expect_lt(update$score_norm, update$start_score_norm)
})

test_that("beta moment update lowers a misspecified starting moment", {
  set.seed(113)
  n <- 1500L
  X <- cbind(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  b_true <- c(1, -0.55)
  is_lowest <- stats::rbinom(
    n,
    size = 1L,
    prob = stats::pnorm(as.numeric(X %*% b_true))
  )
  Y <- ifelse(is_lowest == 1L, 0L, rep(c(1L, 2L), length.out = n))

  update <- .update_beta_npmle_moment(
    Y = Y,
    X = X,
    b_start = c(1, 1.25),
    Ftilde = stats::pnorm,
    ftilde = stats::dnorm,
    fixed_index = 1L
  )

  expect_equal(update$b[1L], 1)
  expect_true(update$updated)
  expect_lt(update$objective, update$start_objective)
  expect_lt(update$score_norm, update$start_score_norm)
})

test_that("alternating first stage records diagnostics and ends with a refit", {
  .skip_semiparametric_stack()
  data <- .make_semiparametric_first_stage_data(n = 150L, seed = 114L)

  fit <- .run_semiparametric_first_stage(
    Y00 = data$Y,
    X00 = data$X,
    fixed_index = 1L,
    min_iter = 2L,
    max_iter = 3L,
    stable_passes = 1L
  )

  expect_gte(fit$iterations, 2L)
  expect_lte(fit$iterations, 3L)
  expect_length(fit$trace, fit$iterations)
  expect_named(
    fit$trace[[1L]],
    c(
      "iteration", "b", "alpha", "parameter_change",
      "beta_score_norm", "alpha_score", "beta_update"
    )
  )
  expect_equal(fit$trace[[fit$iterations]]$b, fit$b)
  expect_equal(fit$trace[[fit$iterations]]$alpha, fit$alpha)
  expect_true(all(is.finite(c(
    fit$final_beta_score,
    fit$final_alpha_score
  ))))
  expect_type(fit$cycle_detected, "logical")
  if (fit$cycle_detected) expect_false(fit$fixed_point_stable)
  expect_true(
    fit$stop_reason %in% c(
      "stable", "max_iter", "stable_update_but_final_score_large"
    ) || grepl("^cycle_[2-5]$", fit$stop_reason)
  )

  refit <- .compute_joint_npmle(data$Y, data$X, fit$b, fit$alpha)
  expect_false(is.null(refit))
  expect_equal(fit$npmle$knots, refit$knots, tolerance = 1e-8)
  expect_equal(fit$npmle$values, refit$values, tolerance = 1e-8)
})
