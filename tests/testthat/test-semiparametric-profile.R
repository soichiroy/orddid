test_that("profile likelihood selects the largest multistart likelihood", {
  set.seed(211)
  n <- 320L
  X <- cbind(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  beta <- c(-0.9, 0.35)
  kappa <- c(-Inf, 0, 1.15, Inf)
  Y <- .draw_three_category_outcome(X, beta = beta, alpha = 1.15)

  fit <- .step2_profile_likelihood(
    Y = Y,
    X = X,
    Ftilde = stats::pnorm,
    kappa = kappa,
    beta_init = c(2, -1.5),
    xi_init = c(0.4, -0.4),
    ftilde = stats::dnorm,
    multistart = TRUE
  )

  expect_gte(nrow(fit$candidate_fits), 2L)
  converged <- fit$candidate_fits$convergence == 0L
  eligible <- if (any(converged)) converged else rep(TRUE, nrow(fit$candidate_fits))
  expect_equal(fit$loglik, max(fit$candidate_fits$loglik[eligible]))
  selected_row <- fit$candidate_fits$start_id == fit$selected_start &
    fit$candidate_fits$optimizer == fit$optimizer
  expect_equal(sum(selected_row), 1L)
  expect_equal(fit$loglik, fit$candidate_fits$loglik[selected_row])
})

test_that("profile selection prioritizes a converged candidate", {
  fits <- list(
    list(value = 10, convergence = 0L),
    list(value = 8, convergence = 10L),
    list(value = 12, convergence = 0L)
  )
  expect_equal(.select_profile_fit(fits), 1L)

  fits[[1L]]$convergence <- 1L
  fits[[3L]]$convergence <- 1L
  expect_equal(.select_profile_fit(fits), 2L)
})

test_that("nonsmooth profile evaluates both numerical optimizers", {
  set.seed(214)
  n <- 180L
  X <- cbind(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  kappa <- c(-Inf, 0, 1.1, Inf)
  Y <- .draw_three_category_outcome(X, beta = c(-0.8, 0.2), alpha = 1.1)

  fit <- .step2_profile_likelihood(
    Y, X, stats::pnorm, kappa,
    beta_init = c(-0.8, 0.2),
    ftilde = NULL,
    optimizer = "both"
  )

  expect_setequal(unique(fit$candidate_fits$optimizer), c("BFGS", "Nelder-Mead"))
  converged <- fit$candidate_fits$convergence == 0L
  eligible <- if (any(converged)) converged else rep(TRUE, nrow(fit$candidate_fits))
  expect_equal(fit$loglik, max(fit$candidate_fits$loglik[eligible]))
})

test_that("deterministic multistart result is insensitive to supplied start", {
  set.seed(212)
  n <- 360L
  X <- cbind(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  beta <- c(-0.8, 0.25)
  kappa <- c(-Inf, 0, 1.1, Inf)
  Y <- .draw_three_category_outcome(X, beta = beta, alpha = 1.1)

  fit_near <- .step2_profile_likelihood(
    Y, X, stats::pnorm, kappa,
    beta_init = beta,
    xi_init = c(0, 0),
    ftilde = stats::dnorm,
    multistart = TRUE
  )
  fit_far <- .step2_profile_likelihood(
    Y, X, stats::pnorm, kappa,
    beta_init = c(2.5, -2),
    xi_init = c(0.7, -0.7),
    ftilde = stats::dnorm,
    multistart = TRUE
  )

  expect_equal(fit_near$loglik, fit_far$loglik, tolerance = 1e-5)
  expect_equal(unname(fit_near$beta), unname(fit_far$beta), tolerance = 2e-3)
  expect_equal(unname(fit_near$xi), unname(fit_far$xi), tolerance = 2e-3)
})

test_that("analytic and finite-difference profiles agree on the selected fit", {
  set.seed(213)
  n <- 280L
  X <- cbind(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  beta <- c(-0.9, 0.3)
  kappa <- c(-Inf, 0, 1.2, Inf)
  Y <- .draw_three_category_outcome(X, beta = beta, alpha = 1.2)

  fit_gradient <- .step2_profile_likelihood(
    Y, X, stats::pnorm, kappa,
    beta_init = beta,
    ftilde = stats::dnorm
  )
  fit_numeric <- .step2_profile_likelihood(
    Y, X, stats::pnorm, kappa,
    beta_init = beta,
    ftilde = NULL
  )

  expect_equal(fit_gradient$loglik, fit_numeric$loglik, tolerance = 1e-4)
  expect_equal(fit_gradient$beta, fit_numeric$beta, tolerance = 5e-3)
  expect_equal(fit_gradient$xi, fit_numeric$xi, tolerance = 5e-3)
})

test_that("profile likelihood reuses CDF evaluations at shared cutoffs", {
  set.seed(814)
  n <- 100L
  X <- cbind(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  Y <- .draw_three_category_outcome(
    X,
    beta = c(-0.7, 0.25),
    alpha = 1.05,
    xi = c(0.08, -0.04)
  )

  cdf_inputs <- list()
  counted_pnorm <- function(u) {
    cdf_inputs[[length(cdf_inputs) + 1L]] <<- u
    stats::pnorm(u)
  }
  density_inputs <- list()
  counted_dnorm <- function(u) {
    density_inputs[[length(density_inputs) + 1L]] <<- u
    stats::dnorm(u)
  }
  fit <- .step2_profile_likelihood(
    Y, X, counted_pnorm, c(-Inf, 0, 1.05, Inf),
    beta_init = c(-0.65, 0.2),
    xi_init = c(0, 0),
    ftilde = counted_dnorm,
    optimizer = "BFGS",
    multistart = FALSE
  )

  consecutive_duplicates <- vapply(
    seq.int(2L, length(cdf_inputs)),
    function(i) identical(cdf_inputs[[i - 1L]], cdf_inputs[[i]]),
    logical(1)
  )
  expect_gt(length(cdf_inputs), 0L)
  expect_false(any(consecutive_duplicates))
  expect_equal(density_inputs, cdf_inputs)

  # Frozen pre-cache result: reusing a shared boundary must not alter the fit.
  expect_equal(
    unname(c(fit$beta, fit$xi)),
    c(-0.601918677783507, 0.172026061937291,
      0.246298823487888, -0.052926651473940),
    tolerance = 1e-7
  )
  expect_equal(fit$loglik, -89.4314561977679, tolerance = 1e-8)
})
