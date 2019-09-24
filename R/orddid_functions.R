# log-likelihood
log_like_probit <- function(par, Y, cut) {
  mu <- par[1]; ss <- exp(par[2])
  j_min <- min(Y)
  j_max <- max(Y)
  item <- 1
  ll <- 0


  ## evalute the likelihood
  nj <- sum(Y == j_min)
  ll <- nj * log(pnorm((cut[item] - mu)/ss))
  item <- item + 1

  for (j in (j_min+1):(j_max-1)) {
    nj <- sum(Y == j)
    ll <- ll + nj * log(
      pnorm((cut[item] - mu) / ss) -
      pnorm((cut[item-1] - mu) / ss)
    )
    item <- item + 1
  }

  nj <- sum(Y == j_max)
  ll <- ll + nj * log(1 - pnorm((max(cut) - mu) / ss))
  return(-ll)
}

# ordered probit
fit_ord_probit <- function(Y, init = NULL, cut) {
  if (is.null(init)) par_init <- c(0.1, 1)
  fit <- optim(par = par_init, fn = log_like_probit,
  Y = Y, cut = cut, method = 'BFGS')
  return(list(mu = fit$par[1], sd = exp(fit$par[2]), ll = fit$value))
}



## function
ord_did_run <- function(Ynew, Yold, treat, cut) {

  ## fit four probit
  ## Y[d = 0, t = 0]
  fit00 <- fit_ord_probit( Y = Yold[treat == 0], cut = cut )
  ## Y[d = 1, t = 0]
  fit01 <- fit_ord_probit( Y = Yold[treat == 1], cut = cut )
  ## Y[d = 0, t = 1]
  fit10 <- fit_ord_probit( Y = Ynew[treat == 0], cut = cut )

  # observed Y[t=1,d=1]
  fit_tr <- fit_ord_probit( Y = Ynew[treat == 1], cut = cut )

  ## estimate counter factual distribution
  mu11 <- fit01$mu + (fit10$mu - fit00$mu) *
          (fit01$sd / fit00$sd)
  ss11 <- fit10$sd * fit01$sd / fit00$sd

  ## compute the probability
  cut_new <- c(-Inf, cut, Inf)
  Ypred <- rep(NA, (length(cut)+1))
  for (j in 1:(length(cut)+1)) {
    Ypred[j] <- pnorm(cut_new[j+1], mean = mu11, sd = ss11)  -
                pnorm(cut_new[j], mean = mu11, sd = ss11)
  }

  ## compute observed Pr(Y(1) = j | D)
  Yobs <- as.vector(prop.table(table(Ynew[treat==1])))


  return(list("Y1" = Yobs, "Y0" = Ypred, 'mu11' = mu11, 'ss11' = ss11,
    fit00 = fit00, fit01 = fit01, fit10 = fit10, fit_tr = fit_tr
  ))
}




#
#
# variance_delta <- function(kappa, mu_point, sd_point, mu_vec, sd_vec) {
#   VCOV <- diag(c(var(mu_vec), var(sd_vec)))
#   VCOV[VCOV==0] <- cov(mu_vec, sd_vec)
#   zj <- sqrt(2) * (kappa - mu_point) / sd_point
#   grad <- c(
#     sqrt(2) / sqrt(pi) * exp(-zj^2) / sd_point,
#     exp(-zj^2) * zj / (sqrt(pi) * sd_point)
#   )
#
#   varj <- as.vector(grad %*% VCOV %*% grad)
#   return(varj)
#
# }
