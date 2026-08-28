.make_sparse_support_panel <- function() {
  make_rows <- function(ids, post, treat, y) {
    data.frame(
      id = ids,
      Y = y,
      post = post,
      treat = treat,
      stringsAsFactors = FALSE
    )
  }

  control_ids <- paste0("c", seq_len(30L))
  treated_ids <- paste0("t", seq_len(30L))

  rbind(
    make_rows(control_ids, FALSE, FALSE, rep(1:3, each = 10L)),
    make_rows(control_ids, TRUE, FALSE, rep(c(2L, 3L, 1L), each = 10L)),
    make_rows(
      treated_ids, FALSE, TRUE,
      c(1L, rep(2L, 14L), rep(3L, 15L))
    ),
    make_rows(
      treated_ids, TRUE, TRUE,
      c(rep(1L, 5L), rep(2L, 10L), rep(3L, 15L))
    )
  )
}


test_that("fixed-support helpers distinguish legacy and covariate requirements", {
  cells <- list(
    Y00 = rep(0:2, 3L),
    Y01 = rep(0:2, 3L),
    Y10 = rep(0:2, 3L),
    Y11 = rep(0:2, 3L)
  )

  expect_equal(.orddid_pooled_category_count(cells), 3L)
  expect_true(.orddid_has_full_category_support(cells, 3L))
  expect_silent(.orddid_require_full_category_support(cells, 3L))

  cells$Y10 <- rep(1:2, 3L)
  expect_false(.orddid_has_full_category_support(cells, 3L))
  expect_error(
    .orddid_require_full_category_support(cells, 3L),
    "Y10 \\(2/3\\)"
  )

  # Covariate estimators fit latent distributions in Y00/Y01/Y10, so all
  # three need full support. Y11 is tabulated against fixed category labels.
  expect_false(.orddid_has_full_category_support(
    cells, 3L, required_cells = c("Y00", "Y01", "Y10")
  ))
  cells$Y10 <- rep(0:2, 3L)
  cells$Y11 <- rep(1:2, 3L)
  expect_false(.orddid_has_full_category_support(cells, 3L))
  expect_true(.orddid_has_full_category_support(
    cells, 3L, required_cells = c("Y00", "Y01", "Y10")
  ))
})


test_that("legacy point estimators reject incomplete cell support", {
  data <- .make_sparse_support_panel()
  data <- data[!(data$treat & !data$post & data$Y == 1L), ]

  expect_error(
    ord_did(
      data, "Y", "post", "treat", "id",
      n_boot = 1L, use_parallel = FALSE
    ),
    "Y10 \\(2/3\\)"
  )
  expect_error(
    run_equivalence(
      data, "Y", "post", "treat", "id",
      n_boot = 1L, use_parallel = FALSE,
      v_range = c(0.25, 0.5, 0.75)
    ),
    "Y10 \\(2/3\\)"
  )
  expect_error(
    run_equivalence_zeta(
      data, "Y", "post", "treat", "id",
      n_boot = 1L, use_parallel = FALSE
    ),
    "Y10 \\(2/3\\)"
  )
})


test_that("legacy bootstraps filter support-deficient draws", {
  data <- .make_sparse_support_panel()

  set.seed(412)
  fit <- ord_did(
    data, "Y", "post", "treat", "id",
    n_boot = 40L, use_parallel = FALSE,
    tau_assumption = "none"
  )
  expect_gt(fit$n_boot_fail, 0L)
  expect_lt(fit$n_boot_fail, 40L)
  expect_equal(fit$n_boot_valid + fit$n_boot_fail, 40L)
  expect_true(all(is.finite(fit$estimate_effects$lower_ci)))
  expect_true(all(is.finite(fit$estimate_effects$upper_ci)))

  set.seed(413)
  eq <- run_equivalence(
    data, "Y", "post", "treat", "id",
    n_boot = 30L, use_parallel = FALSE,
    v_range = c(0.25, 0.5, 0.75)
  )
  expect_gt(eq$n_boot_fail, 0L)
  expect_equal(eq$n_boot_valid + eq$n_boot_fail, 30L)
  expect_true(all(is.finite(eq$estimates$lower_ci)))
  expect_true(all(is.finite(eq$estimates$upper_ci)))
})


.make_staggered_support_panel <- function() {
  make_cell <- function(prefix, time, group, y) {
    data.frame(
      unit = paste0(prefix, seq_along(y)),
      time = time,
      group = group,
      Y = y,
      stringsAsFactors = FALSE
    )
  }

  rbind(
    make_cell("n", 1, Inf, rep(1:3, each = 4L)),
    make_cell("n", 2, Inf, rep(c(2L, 3L, 1L), each = 4L)),
    make_cell("g", 1, 2, rep(1:3, each = 4L)),
    make_cell("g", 2, 2, rep(1:2, each = 6L))
  )
}


test_that("staggered fits require model-cell support, not treated-post support", {
  data <- .make_staggered_support_panel()
  fit <- .staggered_estimate_all(
    data = data,
    outcome = "Y",
    unit = "unit",
    time = "time",
    group = "group",
    base_period = "previous",
    tau_assumption = "none"
  )

  expect_length(fit$gt_results, 1L)
  expect_equal(fit$gt_results[[1L]]$prob_obs, c(0.5, 0.5, 0))

  missing_baseline <- data[
    !(is.infinite(data$group) & data$time == 1L & data$Y == 1L),
  ]
  expect_error(
    .staggered_estimate_all(
      missing_baseline, "Y", "unit", "time", "group"
    ),
    "never_treated_baseline \\(2/3\\)"
  )

  missing_group_base <- data[
    !(data$group == 2 & data$time == 1L & data$Y == 1L),
  ]
  expect_error(
    .staggered_estimate_all(
      missing_group_base, "Y", "unit", "time", "group"
    ),
    "adoption_group_base \\(2/3\\)"
  )
})


test_that("covariate workhorses require support in fitted cells", {
  Y00 <- rep(1:2, 5L)
  Y01 <- rep(0:2, length.out = 10L)
  Y10 <- rep(0:2, length.out = 10L)
  Y11 <- rep(0:2, length.out = 10L)
  X <- cbind(x = seq_len(10L) / 10)

  expect_error(
    .parametric_covariates_did_once(
      Y00, X, Y01, X, Y10, X, Y11, X,
      tau_assumption = "none"
    ),
    "Y00 \\(2/3\\)"
  )
  expect_error(
    .semiparametric_did_once(
      Y00, X, Y01, X, Y10, X, Y11, X,
      tau_assumption = "none"
    ),
    "Y00 \\(2/3\\)"
  )

  expect_error(
    .parametric_covariates_did_once(
      rep(0:2, length.out = 10L), X,
      rep(1:2, length.out = 10L), X,
      Y10, X, Y11, X,
      tau_assumption = "none"
    ),
    "Y01 \\(2/3\\)"
  )
  expect_error(
    .semiparametric_did_once(
      rep(0:2, length.out = 10L), X,
      rep(1:2, length.out = 10L), X,
      Y10, X, Y11, X,
      tau_assumption = "none"
    ),
    "Y01 \\(2/3\\)"
  )
})


test_that("covariate public paths reject support loss in profile cells", {
  data <- .make_sparse_support_panel()
  data$x <- ave(seq_len(nrow(data)), interaction(data$post, data$treat),
                FUN = function(z) seq_along(z) / length(z))
  data <- data[!(!data$treat & data$post & data$Y == 1L), ]

  for (method in c("parametric_covariates", "semiparametric")) {
    expect_error(
      ord_did(
        data, "Y", "post", "treat", "id",
        n_boot = 0L, use_parallel = FALSE,
        method = method, covariates = "x"
      ),
      "Y01 \\(2/3\\)"
    )
  }
})


test_that("equivalence bootstraps retain matrix shape for one draw", {
  data <- .make_sparse_support_panel()
  data$one_cluster <- 1L

  eq <- run_equivalence(
    data, "Y", "post", "treat", "one_cluster",
    n_boot = 1L, use_parallel = FALSE,
    v_range = c(0.25, 0.5, 0.75)
  )
  expect_equal(eq$n_boot_valid, 1L)
  expect_equal(eq$n_boot_fail, 0L)

  eq_zeta <- run_equivalence_zeta(
    data, "Y", "post", "treat", "one_cluster",
    n_boot = 1L, use_parallel = FALSE
  )
  expect_equal(eq_zeta$n_boot_valid, 1L)
  expect_equal(eq_zeta$n_boot_fail, 0L)
  expect_equal(dim(eq_zeta$boot_estimates), c(1L, 3L))
})


test_that("covariate support checks run on parallel workers", {
  detected_cores <- parallel::detectCores()
  skip_if(is.na(detected_cores) || detected_cores < 2L)
  data <- .make_semiparametric_panel(n_per_group = 70L, seed = 719L)

  set.seed(720)
  fit <- ord_did(
    data = data,
    outcome = "Y",
    post = "post",
    treat = "treat",
    cluster = "id",
    n_boot = 2L,
    use_parallel = TRUE,
    n_cores = 2L,
    method = "parametric_covariates",
    covariates = c("X1", "X2"),
    tau_assumption = "none"
  )

  expect_equal(fit$n_boot_valid + fit$n_boot_fail, 2L)
  expect_equal(dim(fit$bootstrap_effects), c(2L, 3L))
})


test_that("legacy support checks run on parallel workers", {
  detected_cores <- parallel::detectCores()
  skip_if(is.na(detected_cores) || detected_cores < 2L)
  data <- .make_sparse_support_panel()
  data$one_cluster <- 1L

  fit <- ord_did(
    data, "Y", "post", "treat", "one_cluster",
    n_boot = 2L, use_parallel = TRUE, n_cores = 2L,
    tau_assumption = "none"
  )

  expect_equal(fit$n_boot_valid, 2L)
  expect_equal(fit$n_boot_fail, 0L)
})
