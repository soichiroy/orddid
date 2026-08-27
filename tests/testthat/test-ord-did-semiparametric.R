test_that("public semiparametric API supports I-spline and ridge", {
  .skip_semiparametric_stack()
  data <- .make_semiparametric_panel(n_per_group = 70L, seed = 301L)

  fit_default <- .fit_semiparametric_test_panel(data, n_boot = 0L)
  fit_ispline <- .fit_semiparametric_test_panel(
    data,
    n_boot = 0L,
    semiparametric_regularization = "ispline"
  )
  fit_ridge <- .fit_semiparametric_test_panel(
    data,
    n_boot = 0L,
    semiparametric_regularization = "ridge"
  )
  fit_probit <- ord_did(
    data = data,
    outcome = "Y",
    post = "post",
    treat = "treat",
    cluster = "id",
    n_boot = 0L,
    method = "parametric_covariates",
    covariates = c("X1", "X2")
  )

  expect_equal(
    fit_default$estimate_effects,
    fit_ispline$estimate_effects
  )
  expect_equal(
    fit_default$relative_effects,
    fit_ispline$relative_effects
  )
  expect_identical(
    fit_default$inputs$semiparametric_regularization,
    "ispline"
  )
  expect_identical(
    fit_ridge$inputs$semiparametric_regularization,
    "ridge"
  )
  expect_equal(
    fit_ridge$point_diagnostics$regularization$epsilon,
    1 / 140
  )
  expect_equal(nrow(fit_ridge$estimate_effects), 3L)
  expect_equal(nrow(fit_probit$estimate_effects), 3L)
  expect_true(all(is.finite(fit_ridge$estimate_effects$effect)))
  expect_true(all(is.finite(fit_probit$estimate_effects$effect)))
  expect_identical(fit_probit$inputs$method, "parametric_covariates")
  expect_equal(
    fit_ispline$point_diagnostics$first_stage$beta,
    fit_ridge$point_diagnostics$first_stage$beta
  )
  expect_equal(
    fit_ispline$point_diagnostics$first_stage$alpha,
    fit_ridge$point_diagnostics$first_stage$alpha
  )
  expect_equal(
    fit_ispline$point_diagnostics$first_stage$npmle_values,
    fit_ridge$point_diagnostics$first_stage$npmle_values
  )
  expect_equal(fit_ispline$point_diagnostics$first_stage$iterations, 2L)
  expect_true(fit_ispline$point_diagnostics$first_stage$completed)
  expect_identical(
    fit_ispline$point_diagnostics$first_stage$stop_reason,
    "two_updates"
  )
  for (profile in fit_ispline$point_diagnostics$profile) {
    expect_setequal(
      unique(profile$candidate_fits$optimizer),
      c("BFGS", "Nelder-Mead")
    )
    expect_false(any(profile$candidate_fits$analytic_gradient))
  }
  ## Fixed-cell regression values for the consolidated estimator. These are
  ## intentionally tied to the canonical first-stage and multistart choices.
  expect_equal(
    fit_ispline$estimate_effects$effect,
    c(-0.0927329483900696, 0.113966265808133, -0.0212333174180639),
    tolerance = 1e-4
  )
  expect_equal(
    fit_ridge$estimate_effects$effect,
    c(-0.135733862397050, 0.110015567962128, 0.0257182944349219),
    tolerance = 1e-4
  )
  expect_equal(
    fit_ispline$point_diagnostics$first_stage$beta,
    c(-1, 0.2499660),
    tolerance = 1e-4
  )
  expect_equal(
    fit_ispline$point_diagnostics$first_stage$alpha,
    0.9779432,
    tolerance = 1e-4
  )
  expect_equal(dim(fit_default$bootstrap_effects), c(0L, 3L))
  expect_equal(dim(fit_ridge$bootstrap_effects), c(0L, 3L))
  expect_equal(fit_default$n_boot_fail, 0L)
  expect_equal(fit_ridge$n_boot_fail, 0L)
})

test_that("semiparametric bootstrap exposes one category-effect row per draw", {
  .skip_semiparametric_stack()
  data <- .make_semiparametric_panel(n_per_group = 70L, seed = 302L)
  set.seed(303)

  fit <- suppressWarnings(ord_did(
    data = data,
    outcome = "Y",
    post = "post",
    treat = "treat",
    cluster = "id",
    n_boot = 1L,
    method = "semiparametric",
    covariates = c("X1", "X2"),
    fixed_index = 1L,
    semiparametric_regularization = "ispline"
  ))

  expect_equal(dim(fit$bootstrap_effects), c(1L, 3L))
  expect_equal(nrow(fit$estimate_effects), 3L)
  expect_equal(
    sum(stats::complete.cases(fit$bootstrap_effects)) + fit$n_boot_fail,
    1L
  )
})

test_that("Appendix intercept design runs all three covariate estimators", {
  .skip_semiparametric_stack()
  data <- .make_semiparametric_panel(n_per_group = 60L, seed = 304L)
  data$intercept <- 1

  fit_method <- function(method, regularization = "ispline") {
    ord_did(
      data = data,
      outcome = "Y",
      post = "post",
      treat = "treat",
      cluster = "id",
      n_boot = 0L,
      method = method,
      covariates = c("intercept", "X1", "X2"),
      fixed_index = 2L,
      semiparametric_regularization = regularization
    )
  }

  fit_probit <- fit_method("parametric_covariates")
  fit_ispline <- fit_method("semiparametric", "ispline")
  fit_ridge <- fit_method("semiparametric", "ridge")

  expect_true(all(is.finite(fit_probit$estimate_effects$effect)))
  expect_true(all(is.finite(fit_ispline$estimate_effects$effect)))
  expect_true(all(is.finite(fit_ridge$estimate_effects$effect)))
  expect_equal(
    unname(fit_ispline$point_diagnostics$first_stage$constant_index),
    1L
  )
  expect_equal(fit_ispline$point_diagnostics$first_stage$beta[2L], -1)
  expect_equal(
    fit_ispline$point_diagnostics$first_stage$beta,
    fit_ridge$point_diagnostics$first_stage$beta
  )
})

test_that("public semiparametric API requires covariates", {
  data <- data.frame(
    id = seq_len(12L),
    Y = rep(0:2, 4L),
    post = rep(c(FALSE, TRUE), 6L),
    treat = rep(c(FALSE, TRUE), each = 6L)
  )

  expect_error(
    ord_did(
      data = data,
      outcome = "Y",
      post = "post",
      treat = "treat",
      cluster = "id",
      n_boot = 0L,
      method = "semiparametric"
    ),
    "covariates must be specified"
  )
})
