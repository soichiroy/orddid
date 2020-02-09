
##
## test code for ordered probit with group indicator
##

context("ord_probit_group")



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
run_sim_ord_probit_gr <- function(n_sim, n_obs, mu, sd, cutoffs) {
  est_save <- list()
  count <- 1
  for (i in 1:n_sim) {
    Y00 <- Y_gen(n_obs = n_obs, mu = mu[1], sd = sd[1], cutoffs = cutoffs)
    Y01 <- Y_gen(n_obs = n_obs, mu = mu[2], sd = sd[2], cutoffs = cutoffs)
    Y10 <- Y_gen(n_obs = n_obs, mu = mu[3], sd = sd[3], cutoffs = cutoffs)
    Y   <- c(Y00, Y01, Y10)
    gr  <- c(rep(1, n_obs), rep(2, n_obs), rep(3, n_obs))

    ## check conditions
    c00 <- length(table(Y00)) == (length(cutoffs)+1)
    c01 <- length(table(Y01)) == (length(cutoffs)+1)
    c10 <- length(table(Y10)) == (length(cutoffs)+1)
    if (isTRUE(all(c00, c01, c10))) {
      ## use data only when we observe full categories
      est <- fit_ord_probit_gr(Y, id_group = gr, cut = c(0, 1))
      est_save[[count]] <- est
      count <- count + 1
    }
  }
  return(est_save)
}


# debugonce(fit_ord_probit_gr)
# debugonce(log_like_probit_group)
# xx <- run_sim_ord_probit_gr(
#   n_sim = 1, n_obs = 1000,
#   mu = c(0.5, 0.5, 0.5),
#   sd = c(1.5, 1.5, 1.5),
#   cutoffs = c(0, 1, 1.5)
# )
#

test_that("accurate parameter estimation (J = 3)", {
  n_sim  <- 150
  n_obs  <- 2000
  set.seed(1234)
  fit <- run_sim_ord_probit_gr(
    n_sim = n_sim, n_obs = n_obs,
    mu = c(0.5, 0.5, 0.5),
    sd = c(1.5, 1.5, 1.5),
    cutoffs = c(0, 1, 1.5)
  )


  ## compute bias
  bias_mu00 <- mean(sapply(fit, function(x) x$mu[1])) - 0.5
  bias_mu01 <- mean(sapply(fit, function(x) x$mu[2])) - 0.5
  bias_mu10 <- mean(sapply(fit, function(x) x$mu[3])) - 0.5
  bias_sd00 <- mean(sapply(fit, function(x) x$sd[1])) - 1.5
  bias_sd01 <- mean(sapply(fit, function(x) x$sd[2])) - 1.5
  bias_sd10 <- mean(sapply(fit, function(x) x$sd[3])) - 1.5

  # cat("bias_sd00 = ", bias_sd00)

  ## checks
  expect_lte(abs(bias_mu00), 0.01)
  expect_lte(abs(bias_mu01), 0.01)
  expect_lte(abs(bias_mu10), 0.01)
  expect_lte(abs(bias_sd00), 0.02)
  expect_lte(abs(bias_sd01), 0.02)
  expect_lte(abs(bias_sd10), 0.02)
})



test_that("accurate parameter estimation, eq-param (J = 4)", {
  n_sim  <- 150
  n_obs  <- 2000
  set.seed(1234)
  fit <- run_sim_ord_probit_gr(
    n_sim = n_sim, n_obs = n_obs,
    mu = c(0.5, 0.5, 0.5),
    sd = c(1.5, 1.5, 1.5),
    cutoffs = c(0, 1, 1.5)
  )


  ## compute bias
  bias_mu00 <- mean(sapply(fit, function(x) x$mu[1])) - 0.5
  bias_mu01 <- mean(sapply(fit, function(x) x$mu[2])) - 0.5
  bias_mu10 <- mean(sapply(fit, function(x) x$mu[3])) - 0.5
  bias_sd00 <- mean(sapply(fit, function(x) x$sd[1])) - 1.5
  bias_sd01 <- mean(sapply(fit, function(x) x$sd[2])) - 1.5
  bias_sd10 <- mean(sapply(fit, function(x) x$sd[3])) - 1.5
  bias_ct   <- mean(sapply(fit, function(x) x$cutoff[3])) - 1.5

  # cat("bias_sd00 = ", bias_sd00)

  ## checks
  expect_lte(abs(bias_mu00), 0.01)
  expect_lte(abs(bias_mu01), 0.01)
  expect_lte(abs(bias_mu10), 0.01)
  expect_lte(abs(bias_sd00), 0.02)
  expect_lte(abs(bias_sd01), 0.02)
  expect_lte(abs(bias_sd10), 0.02)
  expect_lte(abs(bias_ct), 0.01)
})



test_that("accurate parameter estimation (J = 5)", {
  n_sim  <- 150
  n_obs  <- 2000
  set.seed(1234)
  fit <- run_sim_ord_probit_gr(
    n_sim = n_sim, n_obs = n_obs,
    mu = c(0.5, 1, -0.25),
    sd = c(1.5, 2, 1),
    cutoffs = c(0, 1, 1.25, 1.55)
  )


  ## compute bias
  bias_mu00 <- mean(sapply(fit, function(x) x$mu[1])) - 0.5
  bias_mu01 <- mean(sapply(fit, function(x) x$mu[2])) - 1
  bias_mu10 <- mean(sapply(fit, function(x) x$mu[3])) - (-0.25)
  bias_sd00 <- mean(sapply(fit, function(x) x$sd[1])) - 1.5
  bias_sd01 <- mean(sapply(fit, function(x) x$sd[2])) - 2
  bias_sd10 <- mean(sapply(fit, function(x) x$sd[3])) - 1
  bias_ct1  <- mean(sapply(fit, function(x) x$cutoff[3])) - 1.25
  bias_ct2  <- mean(sapply(fit, function(x) x$cutoff[4])) - 1.55

  # cat("bias_sd00 = ", bias_sd00)

  ## checks
  expect_lte(abs(bias_mu00), 0.01)
  expect_lte(abs(bias_mu01), 0.01)
  expect_lte(abs(bias_mu10), 0.01)
  expect_lte(abs(bias_sd00), 0.02)
  expect_lte(abs(bias_sd01), 0.02)
  expect_lte(abs(bias_sd10), 0.02)
  expect_lte(abs(bias_ct1), 0.01)
  expect_lte(abs(bias_ct2), 0.01)
})
