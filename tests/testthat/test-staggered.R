## Test staggered adoption design

test_that("ord_did_staggered runs on simulated staggered data", {
  set.seed(123)

  ## Simulate: T=3 periods (1,2,3), 2 groups + never-treated
  ## Group 2: adopts at t=2, Group 3: adopts at t=3, Inf: never-treated
  n_per_group <- 100
  kappa <- c(0, 1)

  gen_y <- function(n, mu, sd_y) {
    Ystar <- rnorm(n, mean = mu, sd = sd_y)
    cut_y <- c(-Inf, kappa, Inf)
    as.integer(cut(Ystar, breaks = cut_y, labels = FALSE))
  }

  rows <- list()
  for (gg in c(2, 3, Inf)) {
    for (tt in 1:3) {
      ids <- if (is.finite(gg)) {
        paste0("g", gg, "_", seq_len(n_per_group))
      } else {
        paste0("ginf_", seq_len(n_per_group))
      }

      ## mu varies by time for everyone (parallel trends on latent)
      mu_base <- 0.2 * (tt - 1)

      ## Treatment effect only if t >= g
      te <- if (is.finite(gg) && tt >= gg) 0.4 else 0

      Y <- gen_y(n_per_group, mu = mu_base + te, sd_y = 1)

      rows[[length(rows) + 1]] <- data.frame(
        unit  = ids,
        time  = tt,
        group = gg,
        Y     = Y,
        stringsAsFactors = FALSE
      )
    }
  }

  dat <- do.call(rbind, rows)

  result <- ord_did_staggered(
    data = dat,
    outcome = "Y",
    unit = "unit",
    time = "time",
    group = "group",
    cluster = "unit",
    n_boot = 30,
    base_period = "previous",
    weights = "uniform"
  )

  expect_s3_class(result, "orddid_staggered")
  expect_true(nrow(result$group_time_effects) > 0)
  expect_true(nrow(result$group_time_tau) > 0)
  expect_true(!is.null(result$aggregate_effects$agg_zeta))
})

test_that("staggered with no treatment effect gives near-zero zeta", {
  set.seed(456)

  n_per <- 200
  kappa <- c(0, 1)

  rows <- list()
  for (gg in c(2, Inf)) {
    for (tt in 1:3) {
      ids <- if (is.finite(gg)) {
        paste0("g", gg, "_", seq_len(n_per))
      } else {
        paste0("ginf_", seq_len(n_per))
      }

      mu <- 0.1 * (tt - 1)  # no treatment effect
      Ystar <- rnorm(n_per, mean = mu, sd = 1)
      Y <- as.integer(cut(Ystar, breaks = c(-Inf, kappa, Inf), labels = FALSE))

      rows[[length(rows) + 1]] <- data.frame(
        unit = ids, time = tt, group = gg, Y = Y,
        stringsAsFactors = FALSE
      )
    }
  }

  dat <- do.call(rbind, rows)

  result <- ord_did_staggered(
    data = dat,
    outcome = "Y",
    unit = "unit",
    time = "time",
    group = "group",
    n_boot = 20,
    base_period = "previous"
  )

  ## Under no treatment, zeta should be near zero
  zeta_vals <- result$group_time_effects$effect
  expect_true(all(abs(zeta_vals) < 0.15))
})
