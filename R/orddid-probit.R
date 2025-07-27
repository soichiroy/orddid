#' Log-likelihood Function for ordered Probit
#' @keywords internal
#' @importFrom purrr map_dbl
.LogLikelihoodProbit <- function(par, Y, cutoff_val = NULL) {
  # Extract parameters
  # * The first element of par is mean: mu_dt
  # * The second element is standard deviation (log-scale): sd_dt
  mu_y <- par[1]

  # Fix the variance to 1 if the cutoffs need to be estimated.
  sd_y <- ifelse(is.null(cutoff_val), 1, exp(par[2]))

  if (is.null(cutoff_val) && length(par) > 2) {
    # Estimate cutoff if not provided
    # Fix the first element as 0
    cut <- 0
    cut <- c(cut, cumsum(exp(par[-c(1, 2)])))
  } else {
    # Use the provided cutoff values
    cut <- cutoff_val
  }
  # Add end points
  cut <- c(-Inf, cut, Inf)

  # Data summary
  j_min <- min(Y)
  j_max <- max(Y)
  item <- 1
  ll <- 0

  # Evaluate the likelihood
  loglikelihood <- purrr::map_dbl(1:j_max, function(j) {
    nj <- sum(Y == j)
    ll <- nj *
      log(pnorm((cut[j + 1] - mu_y) / sd_y) - pnorm((cut[j] - mu_y) / sd_y))
    return(ll)
  })

  # Return the negative log-likelihood for minimization
  return(-sum(loglikelihood))
}


#' Fit Probit Regression
#'
#' Fitting ordered probit regression
#' @param Y A vector of outcome
#' @keywords internal
.FitOrderedProbit <- function(Y, init = NULL, cutoff_val = NULL) {
  ## input check
  ord_probit_check_cat(Y)

  # Set initial values
  if (is.null(init)) {
    par_init <- c(0, 1)  # mu and sd
  }

  if (is.null(cutoff_val)) {
    par_init <- c(par_init, runif(length(unique(Y)) - 2))
  }

  ## fit ordered probit
  fit <- optim(
    par = par_init,
    fn = .LogLikelihoodProbit,
    Y = Y,
    cutoff_val = cutoff_val,
    method = 'BFGS'
  )

  # calculate the cutoff values
  if (is.null(cutoff_val)) {
    cutoff <- c(0, cumsum(exp(fit$par[-c(1, 2)])))
  } else {
    cutoff <- cutoff_val
  }

  return(list(
    mu = fit$par[1],
    sd = ifelse(is.null(cutoff_val), 1, exp(fit$par[2])),
    cutoff = cutoff,
    ll = fit$value
  ))
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
log_like_probit_group <- function(
  par,
  Y,
  cut,
  id_group,
  n_group,
  n_gj,
  j_min,
  j_max
) {
  ## prep parameters
  # n_group <- length(unique(id_group))
  mu <- par[1:n_group]
  ss <- exp(par[(n_group + 1):(n_group * 2)])

  ## if the number of categories is more than 4,
  ## we estimate cutoff
  if (length(par[-c(1:(n_group * 2))]) > 0) {
    cut <- c(cut, cut[length(cut)] + cumsum(exp(par[-c(1:(n_group * 2))])))
  }

  ## evalute the likelihood
  ll <- 0
  for (g in 1:n_group) {
    ## get group params
    mu_g <- mu[g]
    ss_g <- ss[g]
    item <- 1

    ## evaluate ll function
    nj <- n_gj[g, item]
    ll <- ll + nj * log(pnorm(cut[item], mean = mu_g, sd = ss_g))
    item <- item + 1

    for (j in (j_min + 1):(j_max - 1)) {
      nj <- n_gj[g, item]
      ll <- ll +
        nj *
          log(
            pnorm(cut[item], mean = mu_g, sd = ss_g) -
              pnorm(cut[item - 1], mean = mu_g, sd = ss_g)
          )
      item <- item + 1
    }

    nj <- n_gj[g, item]
    ll <- ll + nj * log(1 - pnorm(max(cut), mean = mu_g, sd = ss_g))
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
  j_min <- min(Y)
  j_max <- max(Y)
  n_cat <- (j_max - j_min) + 1

  ## initial value setup ----------------------------------
  if (is.null(init)) {
    par_init <- c(rep(0.1, n_group), rep(0.5, n_group))
    ## initialize cutoff
    if (n_cat > 3) {
      n_cutoffs <- n_cat - 1
      ## first two cutoffs are fixed: default is c(0, 1)
      par_init <- c(par_init, rep(-0.5, n_cutoffs - 2))
    }
  }

  ## intput prep; count statistics ------------------------
  n_gj <- count_nj(Y, id_group, n_group, n_cat, j_max, j_min)

  ## fit ordered probit -----------------------------------
  fit <- optim(
    par = par_init,
    fn = log_like_probit_group,
    Y = Y,
    cut = cut,
    id_group = id_group,
    n_group = n_group,
    n_gj = n_gj,
    j_min = j_min,
    j_max = j_max,
    method = 'BFGS'
  )

  ## return object ----------------------------------------
  cutoff <- c(cut, cut[length(cut)] + cumsum(exp(fit$par[-c(1:(n_group * 2))])))
  mu_vec <- fit$par[1:n_group]
  sd_vec <- fit$par[(n_group + 1):(2 * n_group)]

  return(list(mu = mu_vec, sd = exp(sd_vec), cutoff = cutoff, ll = fit$value))
}


#' Counting category selections by groups
#' @keywords internal
count_nj <- function(Y, id_group, n_group, n_cat, j_max, j_min) {
  n_gj <- matrix(0, nrow = n_group, ncol = n_cat)
  for (g in 1:n_group) {
    item <- 1
    Yg <- Y[id_group == g]
    n_gj[g, item] <- sum(Yg == j_min)
    item <- item + 1
    for (j in (j_min + 1):(j_max - 1)) {
      n_gj[g, item] <- sum(Yg == j)
      item <- item + 1
    }
    n_gj[g, item] <- sum(Yg == j_max)
  }

  return(n_gj)
}
