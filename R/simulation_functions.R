##
## Functions used for simulation studies 
##


latent_sim <- function(link, mu, sd, n_obs) {
  if (link == 'probit') {
    Y <- rnorm(n_obs, mean = mu, sd = sd)
  } else if (link == 'logit') {
    Y <- rlogis(n_obs, location = mu, scale = sd)
  } else {
    stop("not supported link function \n")
  }

  return(Y)
}


orddid_sim <- function(n_unit, p, mu, sd, mu_obs, sd_obs, n_choice,
link = 'probit') {
  n1 <- ceiling(n_unit * p)
  n0 <- n_unit - n1
  treat <- c(rep(1, n1), rep(0, n0))

  kappa <- seq(-0.5, 1, length.out = n_choice-1)
  kappaJ <- c(-Inf, kappa, Inf)

  # mu10 <- mu$mu10
  # mu00 <- mu$mu00
  # mu10 <- mu$mu10

  ## compute mu11, sd11
  mu$mu11 <- mu$mu10 + (mu$mu01 - mu$mu00) / (sd$sd00 / sd$sd10)
  sd$sd11 <- sd$sd10 * sd$sd01 / sd$sd00

  ## "true effect"
  prY1 <- prY0 <- rep(NA, n_choice)
  for (j in 1:n_choice) {
    prY0[j] <- pnorm(kappaJ[j+1], mean = mu$mu11, sd = sd$sd11) -
                pnorm(kappaJ[j], mean = mu$mu11, sd = sd$sd11)
    prY1[j] <- pnorm(kappaJ[j+1], mean = mu_obs, sd = sd_obs) -
                pnorm(kappaJ[j], mean = mu_obs, sd = sd_obs)
  }

  ## simulate data
  lu <- list()
  lu$Y00 <- latent_sim(link, mu = mu$mu00, sd = sd$sd00, n0)
  lu$Y01 <- latent_sim(link, mu = mu$mu01, sd = sd$sd01, n0)
  lu$Y10 <- latent_sim(link, mu = mu$mu10, sd = sd$sd10, n1)
  lu$Y11 <- latent_sim(link, mu = mu$mu11, sd = sd$sd11, n1)
  lu$Yobs  <- latent_sim(link, mu = mu_obs, sd = sd_obs, n1)


  ## generate categorical data
  ord <- list()

  ord$Y00 <- as.numeric(cut(lu$Y00, breaks = kappaJ))
  ord$Y01 <- as.numeric(cut(lu$Y01, breaks = kappaJ))
  ord$Y10 <- as.numeric(cut(lu$Y10, breaks = kappaJ))
  ord$Y11 <- as.numeric(cut(lu$Y11, breaks = kappaJ))
  ord$Yobs <- as.numeric(cut(lu$Yobs, breaks = kappaJ))

  dat <- data.frame(Ynew = c(ord$Yobs, ord$Y01), Yold = c(ord$Y10, ord$Y00),
  treat = treat)

  return(list(ord=ord, lu=lu, kappa=kappa, kappaJ=kappaJ, dat = dat,
    mu=mu, sd=sd, n1=n1, n0=n0, prY0 = prY0, prY1 = prY1))
}
