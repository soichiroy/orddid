##
## test code for odered probit

context("ord_probit")



## --------------------------------------------------------------- ##
##                       Define functions                          ##
## --------------------------------------------------------------- ##

## define function for DGP
Y_gen <- function(n_obs, mu, sd, cutoffs) {
  n_cat <- length(cutoffs) + 1
  ## simulate latent variables
  Yutil <- rnorm(n_obs, mean = mu, sd = sd)

  ## create categrical outcomes
  kappaJ <- c(-Inf, cutoffs, Inf)
  Y <- as.numeric(cut(Yutil, breaks = kappaJ))
  return(Y)
}

## simulation function
run_sim_ord_probit <- function(n_sim, n_obs, mu, sd, cutoffs) {
  est_save <- list()
  count <- 1
  for (i in 1:n_sim) {
    Y <- Y_gen(n_obs = n_obs, mu = mu, sd = sd, cutoffs = cutoffs)
    if (length(table(Y)) == (length(cutoffs)+1)) {
      ## use data only when we observe full categories
      est <- fit_ord_probit(Y, cut = c(0, 1))
      est_save[[count]] <- est
      count <- count + 1
    }
  }
  return(est_save)
}


## --------------------------------------------------------------- ##
##                        Testing Functions                        ##
## --------------------------------------------------------------- ##

test_that("run check", {
  # gen outcome
  Y <- Y_gen(n_obs = 500, mu = 1, sd = 2, cutoffs = c(0, 1, 1.5))

  # fit the probit
  debugonce(fit_ord_probit)
  fit <- fit_ord_probit(Y, cut = c(0, 1))
  expect_length(fit, 4)
  expect_length(fit$cutoff, 3)

})


test_that("input check", {
  # missing one category
  Y <- c(rep(0, 10), rep(1, 10), rep(3, 10))
  expect_error(fit_ord_probit(Y, cut = c(0, 1)), "Outcome has a missing categories.")

  # missing two category
  Y <- c(rep(0, 10), rep(3, 10))
  expect_error(fit_ord_probit(Y, cut = c(0, 1)), "Outcome has a missing categories.")
})


test_that("accurate parameter estimation (J = 3)", {
  n_sim  <- 100
  n_obs  <- 5000
  set.seed(1234)
  fit <- run_sim_ord_probit(n_sim, n_obs,
    mu = -0.5, sd = 1.5,
    cutoffs = c(0, 1)
  )

  ## compute bias
  bias_mu <- mean(sapply(fit, function(x) x$mu)) - (-0.5)
  bias_sd <- mean(sapply(fit, function(x) x$sd)) - 1.5

  ## checks
  expect_lte(bias_mu, 0.01)
  expect_lte(bias_sd, 0.01)

})


test_that("accurate parameter estimation (J = 4)", {
  n_sim  <- 100
  n_obs  <- 5000
  set.seed(1234)
  fit <- run_sim_ord_probit(n_sim, n_obs,
    mu = -0.5, sd = 1.5,
    cutoffs = c(0, 1, 1.5)
  )

  ## compute bias
  bias_mu <- mean(sapply(fit, function(x) x$mu)) - (-0.5)
  bias_sd <- mean(sapply(fit, function(x) x$sd)) - 1.5
  bias_ct <- mean(sapply(fit, function(x) x$cutoff[3])) - 1.5

  ## checks
  expect_lte(abs(bias_mu), 0.01)
  expect_lte(abs(bias_sd), 0.01)
  expect_lte(abs(bias_ct), 0.01)

})


test_that("accurate parameter estimation (J = 5)", {
  n_sim  <- 100
  n_obs  <- 5000
  set.seed(1234)
  fit <- run_sim_ord_probit(n_sim, n_obs,
    mu = 0.5, sd = 2,
    cutoffs = c(0, 1, 1.5, 2)
  )

  ## compute bias
  bias_mu <- mean(sapply(fit, function(x) x$mu)) - (0.5)
  bias_sd <- mean(sapply(fit, function(x) x$sd)) - 2
  bias_ct1 <- mean(sapply(fit, function(x) x$cutoff[3])) - 1.5
  bias_ct2 <- mean(sapply(fit, function(x) x$cutoff[4])) - 2

  ## checks
  expect_lte(abs(bias_mu), 0.01)
  expect_lte(abs(bias_sd), 0.01)
  expect_lte(abs(bias_ct1), 0.01)
  expect_lte(abs(bias_ct2), 0.01)

})
