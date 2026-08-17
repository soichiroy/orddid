test_that("stochastic_dominance helper tightens known margins", {
  p1 <- c(0.2, 0.3, 0.5)
  p0 <- c(0.4, 0.4, 0.2)

  unconstrained <- .ComputeRelativeEffect(p1, p0, tau_assumption = "none")
  constrained <- .ComputeRelativeEffect(
    p1, p0, tau_assumption = "stochastic_dominance"
  )

  expect_true(constrained$effect_lb[1] >= unconstrained$effect_lb[1] - 1e-10)
  expect_true(constrained$effect_ub[1] <= unconstrained$effect_ub[1] + 1e-10)
  expect_equal(constrained$effect_lb[1], 0.3125, tolerance = 1e-8)
  expect_equal(constrained$effect_ub[1], 0.5133333333, tolerance = 1e-8)
  expect_equal(constrained$status[1], "ok")
})

test_that("stochastic_dominance bound is narrower under a DGP satisfying the assumption", {
  Pi <- matrix(
    c(
      0.18, 0.06, 0.01,
      0.07, 0.18, 0.07,
      0.01, 0.08, 0.34
    ),
    nrow = 3,
    byrow = TRUE
  )

  check_prd <- function(Pi, tol = 1e-10) {
    J <- nrow(Pi)

    for (k_low in seq_len(J - 1L)) {
      k_high <- k_low + 1L
      for (threshold in seq_len(J - 1L)) {
        lhs <- sum(Pi[seq_len(threshold), k_low]) / sum(Pi[, k_low])
        rhs <- sum(Pi[seq_len(threshold), k_high]) / sum(Pi[, k_high])
        if (lhs < rhs - tol) return(FALSE)
      }
    }

    for (j_low in seq_len(J - 1L)) {
      j_high <- j_low + 1L
      for (threshold in seq_len(J - 1L)) {
        lhs <- sum(Pi[j_low, seq_len(threshold)]) / sum(Pi[j_low, ])
        rhs <- sum(Pi[j_high, seq_len(threshold)]) / sum(Pi[j_high, ])
        if (lhs < rhs - tol) return(FALSE)
      }
    }

    TRUE
  }

  expect_true(check_prd(Pi))

  p1 <- rowSums(Pi)
  p0 <- colSums(Pi)
  tau_true <- sum(Pi * sign(row(Pi) - col(Pi)))

  bounds <- .ComputeRelativeEffect(
    p1, p0, tau_assumption = c("none", "stochastic_dominance")
  )

  width_none <- with(bounds[bounds$tau_assumption == "none", ], effect_ub - effect_lb)
  width_sd <- with(bounds[bounds$tau_assumption == "stochastic_dominance", ], effect_ub - effect_lb)
  sd_row <- bounds[bounds$tau_assumption == "stochastic_dominance", , drop = FALSE]

  expect_equal(sd_row$status[1], "ok")
  expect_true(width_sd < width_none)
  expect_true(sd_row$effect_lb[1] <= tau_true)
  expect_true(tau_true <= sd_row$effect_ub[1])
})

test_that("stochastic_dominance helper handles zero marginal masses", {
  p1 <- c(0.5, 0.5, 0)
  p0 <- c(0.4, 0, 0.6)

  constrained <- .ComputeRelativeEffect(
    p1, p0, tau_assumption = "stochastic_dominance"
  )

  expect_equal(constrained$effect_lb[1], -0.6, tolerance = 1e-8)
  expect_equal(constrained$effect_ub[1], -0.4, tolerance = 1e-8)
  expect_equal(constrained$status[1], "ok")
})

test_that("helper accepts multiple tau assumptions", {
  res <- .ComputeRelativeEffect(
    c(0.2, 0.3, 0.5),
    c(0.4, 0.4, 0.2),
    tau_assumption = c("none", "stochastic_dominance")
  )

  expect_equal(
    res$tau_assumption,
    c("none", "stochastic_dominance")
  )
  expect_equal(nrow(res), 2)
})

test_that("ord_did accepts tau_assumption for the parametric path", {
  set.seed(321)

  n_per_group <- 150
  cutpoints <- c(-Inf, 0, 1, Inf)

  gen_outcome <- function(n, mean_shift) {
    y_star <- rnorm(n, mean = mean_shift, sd = 1)
    as.integer(cut(y_star, breaks = cutpoints, labels = FALSE)) - 1L
  }

  control_ids <- seq_len(n_per_group)
  treated_ids <- seq_len(n_per_group) + n_per_group

  dat <- rbind(
    data.frame(
      id = control_ids,
      Y = gen_outcome(n_per_group, mean_shift = 0.0),
      post = FALSE,
      treat = FALSE
    ),
    data.frame(
      id = control_ids,
      Y = gen_outcome(n_per_group, mean_shift = 0.2),
      post = TRUE,
      treat = FALSE
    ),
    data.frame(
      id = treated_ids,
      Y = gen_outcome(n_per_group, mean_shift = 0.1),
      post = FALSE,
      treat = TRUE
    ),
    data.frame(
      id = treated_ids,
      Y = gen_outcome(n_per_group, mean_shift = 0.6),
      post = TRUE,
      treat = TRUE
    )
  )

  result <- ord_did(
    data = dat,
    outcome = "Y",
    post = "post",
    treat = "treat",
    cluster = "id",
    n_boot = 20,
    tau_assumption = "stochastic_dominance"
  )

  expect_equal(result$inputs$tau_assumption, "stochastic_dominance")
  expect_equal(nrow(result$relative_effects), 1)
  expect_true(result$relative_effects$effect_lb[1] <= result$relative_effects$effect_ub[1])
  expect_true(is.finite(result$relative_effects$effect_lb[1]))
  expect_true(is.finite(result$relative_effects$effect_ub[1]))
})

test_that("ord_did returns both tau assumptions by default", {
  set.seed(654)

  n_per_group <- 100
  cutpoints <- c(-Inf, 0, 1, Inf)

  gen_outcome <- function(n, mean_shift) {
    y_star <- rnorm(n, mean = mean_shift, sd = 1)
    as.integer(cut(y_star, breaks = cutpoints, labels = FALSE)) - 1L
  }

  control_ids <- seq_len(n_per_group)
  treated_ids <- seq_len(n_per_group) + n_per_group

  dat <- rbind(
    data.frame(id = control_ids, Y = gen_outcome(n_per_group, 0.0), post = FALSE, treat = FALSE),
    data.frame(id = control_ids, Y = gen_outcome(n_per_group, 0.2), post = TRUE, treat = FALSE),
    data.frame(id = treated_ids, Y = gen_outcome(n_per_group, 0.1), post = FALSE, treat = TRUE),
    data.frame(id = treated_ids, Y = gen_outcome(n_per_group, 0.6), post = TRUE, treat = TRUE)
  )

  result <- ord_did(
    data = dat,
    outcome = "Y",
    post = "post",
    treat = "treat",
    cluster = "id",
    n_boot = 10
  )

  expect_equal(
    result$relative_effects$tau_assumption,
    c("none", "stochastic_dominance")
  )
  expect_equal(result$inputs$tau_assumption, c("none", "stochastic_dominance"))
})
