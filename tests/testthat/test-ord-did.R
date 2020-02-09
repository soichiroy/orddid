
##
## test code for ordered probit wi
##

context("ord_did")



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
run_sim_ord_did <- function(n_sim, n_obs, mu, sd, cutoffs) {
  est_save <- list()
  count <- 1
  for (i in 1:n_sim) {
    Y00 <- Y_gen(n_obs = n_obs, mu = mu[1], sd = sd[1], cutoffs = cutoffs)
    Y01 <- Y_gen(n_obs = n_obs, mu = mu[2], sd = sd[2], cutoffs = cutoffs)
    Y10 <- Y_gen(n_obs = n_obs, mu = mu[3], sd = sd[3], cutoffs = cutoffs)
    Y11 <- Y_gen(n_obs = n_obs, mu = mu[4], sd = sd[4], cutoffs = cutoffs)

    ## check conditions
    c00 <- length(table(Y00)) == (length(cutoffs)+1)
    c01 <- length(table(Y01)) == (length(cutoffs)+1)
    c10 <- length(table(Y10)) == (length(cutoffs)+1)
    c11 <- length(table(Y11)) == (length(cutoffs)+1)

    Ynew <- c(Y11, Y01)
    Yold <- c(Y10, Y00)
    treat <- c(rep(1, n_obs), rep(0, n_obs))

    if (isTRUE(all(c00, c01, c10, c11))) {
      ## use data only when we observe full categories
      est <- ord_did_run(Ynew = Ynew, Yold = Yold, treat = treat, cut = c(0, 1))
      est_save[[count]] <- est
      count <- count + 1
    }
  }
  return(est_save)
}


compute_effect <- function(mu, sd, cutoffs) {
  mu00 <- mu[1]; mu01 <- mu[2]; mu10 <- mu[3]
  sd00 <- sd[1]; sd01 <- sd[2]; sd10 <- sd[3]

  ## identification
  mu11 <- mu10 + (mu01 - mu00) / (sd00 / sd10)
  sd11 <- sd10 * sd01 / sd00

  Ypred <- rep(NA, length(cutoffs)+1)
  cutoffs_add <- c(-Inf, cutoffs, Inf)
  for (i in 1:(length(cutoffs)+1)) {
    Ypred[i] <- pnorm(cutoffs_add[i+1], mean = mu11, sd = sd11) -
                pnorm(cutoffs_add[i], mean = mu11, sd = sd11)
  }

  return(list(mu11 = mu11, sd11 = sd11, Y0 = Ypred))
}



## --------------------------------------------------------------- ##
##                            Testing                              ##
## --------------------------------------------------------------- ##


test_that("no effect (J = 4)", {
  n_sim <- 150
  n_obs <- 2000
  set.seed(1234)
  fit <- run_sim_ord_did(
    n_sim = n_sim,
    n_obs = n_obs,
    mu = c(0.5, 0.5, 0.5, 0),
    sd = c(1.5, 1.5, 1.5, 2),
    cutoffs = c(0, 1, 1.5)
  )

  ## bias
  bias_mu11 <- mean(sapply(fit, function(x) x$mu11)) - 0.5
  bias_sd11 <- mean(sapply(fit, function(x) x$ss11)) - 1.5

  ## checks
  expect_lte(abs(bias_mu11), 0.01)
  expect_lte(abs(bias_sd11), 0.01)
})


test_that("some effect (J = 4)", {
  n_sim <- 150
  n_obs <- 2000
  mu <- c(0.5, 1.5, -0.5, 0)
  sd <- c(2.5, 1.5, 1, 2)
  ct <- c(0, 1, 1.5)
  set.seed(1234)
  fit <- run_sim_ord_did(
    n_sim = n_sim,
    n_obs = n_obs,
    mu = mu, sd = sd,
    cutoffs = ct
  )

  ## truth
  th <- compute_effect(mu, sd, ct)

  ## bias
  bias_mu11 <- mean(sapply(fit, function(x) x$mu11)) - th$mu11
  bias_sd11 <- mean(sapply(fit, function(x) x$ss11)) - th$sd11
  bias_Y0   <- mean(sapply(fit, function(x) mean(x$Y0 - th$Y0)))

  ## checks
  expect_lte(abs(bias_mu11), 0.01)
  expect_lte(abs(bias_sd11), 0.01)
  expect_lte(abs(bias_Y0), 0.01)
})



test_that("orddid input check", {
  ## different length
  set.seed(1234)
  Y1 <- sample(1:3, size = 1001, replace = TRUE)
  Y0 <- sample(1:3, size = 1000, replace = TRUE)
  treat <- sample(0:1, size = 1000, replace = TRUE)

  expect_error(ord_did(Y1, Y0, treat, n_boot = 2))

  ## differnt length of id_cluster
  set.seed(1234)
  Y1 <- sample(1:3, size = 2000, replace = TRUE)
  Y0 <- sample(1:3, size = 2000, replace = TRUE)
  treat <- sample(0:1, size = 2000, replace = TRUE)
  cluster <- sample(1:4, size = 1000, replace = TRUE)
  expect_error(ord_did(Y1, Y0, treat, cluster, n_boot = 2))
})
