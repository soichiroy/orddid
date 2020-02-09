##
## Probit functions
##


## --------------------------------------------------------------- ##
##                       Probit Regression                         ##
## --------------------------------------------------------------- ##

#' Log-likelihood Function
#' @keywords internal
log_like_probit <- function(par, Y, cut) {
  ## prep parameters
  mu    <- par[1]
  ss    <- exp(par[2])
  j_min <- min(Y)
  j_max <- max(Y)
  item  <- 1
  ll    <- 0

  ## if the number of categories is more than 4,
  ## we estimate cutoff
  if (length(par[-c(1,2)]) > 0) {
    cut <- c(cut, cut[length(cut)] + cumsum(exp(par[-c(1, 2)])))
  }

  ## evalute the likelihood
  nj   <- sum(Y == j_min)
  ll   <- nj * log(pnorm((cut[item] - mu)/ss))
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


#' Probit Regression
#'
#' Fitting ordered probit regression
#' @param Y A vector of outcome
#' @keywords internal
fit_ord_probit <- function(Y, init = NULL, cut) {
  ## input check
  ord_probit_check_cat(Y)

  ## initial value setup
  if (is.null(init)) {
    par_init <- c(0.1, 1)
    ## initialize cutoff
    n_cat <- (max(Y) - min(Y)) + 1
    if (n_cat > 3) {
      n_cutoffs <- n_cat - 1
      ## first two cutoffs are fixed: default is c(0, 1)
      par_init <- c(par_init, rep(-0.5, n_cutoffs-2))
    }
  }

  ## fit ordered probit
  fit <- optim(
    par = par_init, fn = log_like_probit,
    Y = Y, cut = cut, method = 'BFGS'
  )

  ## return cutoff
  cutoff <- c(cut, cut[length(cut)] + cumsum(exp(fit$par[-c(1,2)])))
  return(list(mu = fit$par[1], sd = exp(fit$par[2]), ll = fit$value, cutoff = cutoff))
}



## --------------------------------------------------------------- ##
##               Probit Regression with Group Dummies              ##
## --------------------------------------------------------------- ##


#' Log-likelihood Function
#'
#' @param par A vector of parameters.
#'  means: \code{par[1:3]}; scale: \code{par[4:6]}; cutoff: \code{par[7:]} (only when J > 3)
#' @param Y A vector of categorical outcome.
#' @param cut A set of two cutoffs fixed for identification. Default is \code{cut = c(0, 1)}.
#' @param id_group A vector of group indicator.
#' @keywords internal
log_like_probit_group <- function(par, Y, cut, id_group) {
  ## prep parameters
  n_group <- length(unique(id_group))
  mu    <- par[1:n_group]
  ss    <- exp(par[(n_group+1):(n_group*2)])
  j_min <- min(Y)
  j_max <- max(Y)

  ## if the number of categories is more than 4,
  ## we estimate cutoff
  if (length(par[-c(1:(n_group*2))]) > 0) {
    cut <- c(cut, cut[length(cut)] + cumsum(exp(par[-c(1:(n_group*2))])))
  }

  ## evalute the likelihood
  ll    <- 0
  for (g in 1:n_group) {
    ## get group params
    mu_g <- mu[g];
    ss_g <- ss[g]
    item  <- 1

    ## evaluate ll function
    nj   <- sum(Y[id_group == g] == j_min)
    ll   <- ll + nj * log(pnorm((cut[item] - mu_g) / ss_g))
    item <- item + 1

    for (j in (j_min+1):(j_max-1)) {
      nj <- sum(Y[id_group == g] == j)
      ll <- ll + nj * log(
        pnorm((cut[item] - mu_g) / ss_g) -
        pnorm((cut[item-1] - mu_g) / ss_g)
      )
      item <- item + 1
    }

    nj <- sum(Y[id_group == g] == j_max)
    ll <- ll + nj * log(1 - pnorm((max(cut) - mu_g) / ss_g))
  }

  return(-ll)
}

#' Probit regression with group dummies
#'
#' @param Y A vector of the ordinal outcome.
#' @param id_group A vector of group indicator.
#' @keywords internal
fit_ord_probit_gr <- function(Y, id_group, init = NULL, cut) {
  ## input check
  ord_probit_check_cat(Y)
  n_group <- length(unique(id_group))

  ## initial value setup
  if (is.null(init)) {
    par_init <- c(rep(0.1, n_group), rep(0.5, n_group))
    ## initialize cutoff
    n_cat <- (max(Y) - min(Y)) + 1
    if (n_cat > 3) {
      n_cutoffs <- n_cat - 1
      ## first two cutoffs are fixed: default is c(0, 1)
      par_init <- c(par_init, rep(-0.5, n_cutoffs-2))
    }
  }

  ## fit ordered probit
  fit <- optim(
    par = par_init, fn = log_like_probit_group,
    Y = Y, cut = cut, id_group = id_group,
    method = 'BFGS'
  )

  ## return object
  cutoff <- c(cut, cut[length(cut)] + cumsum(exp(fit$par[-c(1:(n_group*2))])))
  mu_vec <- fit$par[1:n_group]
  sd_vec <- fit$par[(n_group+1):(2*n_group)]

  return(list(mu = mu_vec, sd = exp(sd_vec), cutoff = cutoff, ll = fit$value))
}
