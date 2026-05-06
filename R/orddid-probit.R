#' Base CDF dispatcher for the ordinal latent variable
#'
#' Returns the CDF F(.) of the standard base distribution implied by `link`.
#' `"probit"` uses the standard normal; `"logit"` uses the standard logistic.
#' Used by both the MLE objective and the predicted-probability calculation.
#' @keywords internal
#' @noRd
.BaseCDF <- function(x, link = c("probit", "logit")) {
  link <- match.arg(link)
  switch(link,
    probit = stats::pnorm(x),
    logit  = stats::plogis(x)
  )
}

#' Fit Ordinal Probit/Logit Regression via MLE
#'
#' @param Y A vector of outcome
#' @param init A vector of initial values. The first element is for the mean,
#'   and the second element is for the standard deviation. If NULL, the default
#'   values are used: 0 for mean, 1 for standard deviation.
#' @param cutoff_val A vector of cutoff values. If NULL, the cutoffs are
#'   estimated fix the first cutoff is fixed at 1.
#' @param link Character; `"probit"` (default, standard normal F) or `"logit"`
#'   (standard logistic F).
#' @importFrom stats runif optim
#' @keywords internal
.FitOrderedProbit <- function(Y, init = NULL, cutoff_val = NULL,
                              link = c("probit", "logit")) {
  link <- match.arg(link)

  # Set initial values
  if (is.null(init)) {
    par_init <- c(0, 1) # mu and sd
  }

  if (is.null(cutoff_val)) {
    par_init <- c(par_init, runif(length(unique(Y)) - 2))
  }

  ## fit ordered probit / logit
  fit <- optim(
    par = par_init,
    fn = .LogLikelihoodProbit,
    Y = Y,
    cutoff_val = cutoff_val,
    link = link,
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

#' Log-likelihood Function for ordered Probit/Logit
#' @keywords internal
#' @importFrom purrr map_dbl
#' @noRd
.LogLikelihoodProbit <- function(par, Y, cutoff_val = NULL,
                                 link = c("probit", "logit")) {
  link <- match.arg(link)

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
  cut <- .PadInf(cut)

  # Evaluate the likelihood over sorted categories
  cats <- sort(unique(Y))
  loglikelihood <- purrr::map_dbl(seq_along(cats), function(k) {
    nj <- sum(Y == cats[k])
    nj * log(.BaseCDF((cut[k + 1] - mu_y) / sd_y, link = link) -
             .BaseCDF((cut[k]     - mu_y) / sd_y, link = link))
  })

  # Return the negative log-likelihood for minimization
  return(-sum(loglikelihood))
}
