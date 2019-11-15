#
# title: functions for diagnostic test
#
# author: Soichiro Yamauchi
#
# note:
#   steps for the test
#    1. Estimate all parameters using ordered probit
#    2. Obtain variance of parameters using bootstrap
#    3. Compute t(v) using qd() function
#    4. Compute Var(t(v))
#


## --------------------------------------------------------------- ##
##                      Main Functions                             ##
## --------------------------------------------------------------- ##


#' Conduct an equivalence test of q1(v) = q0(v)
#'
#' @param object a fitted object from \code{orddid}.
#' @param alpha level of a test. Should be between 0 and 1. Default is 0.05.
#' @param threshold an equivalance threshold.
#' @return a list of outputs.
#' @importFrom Matrix bdiag
#' @export
equivalence_test <- function(object, alpha = 0.05, threshold = NULL) {
  ## check input
  if (!("orddid.fit" %in% class(object))) {
    stop("object should be the output of orddid function!")
  }

  ## compute delta if not provided
  if( is.null(threshold) ) {
    threshold <- calc_threshold(object)
  }

  ## extract information
  theta <- object$fit$theta
  nobs  <- attr(object, 'n') * 2
  vcov  <- cov(object$boot_params)

  ## zero out under independent unit assumption
  v0 <- vcov[1:4, 1:4]; v1 <- vcov[5:8, 5:8]
  vcov <- as.matrix(bdiag(list(v0, v1)))


  ## compute t(v)
  v_range      <- seq(0.001, 0.999, by = 0.01)
  tv <- tv_var <- rep(NA, length(v_range))
  Uv <- Lv     <- rep(NA, length(v_range))
  zscore       <- rep(NA, length(v_range))
  Upvalue      <- Lpvalue <- rep(NA, length(v_range))

  for (v in 1:length(v_range)) {
    ## compute q1(v) - q0(v)
    tv[v] <- qd(v_range[v], mu1 = theta$mu11, mu0 = theta$mu10,
                s1 = theta$sd11, s0 = theta$sd10) -
             qd(v_range[v], mu1 = theta$mu01, mu0 = theta$mu00,
                s1 = theta$sd01, s0 = theta$sd00)

    ## compute Var(t(v))
    grad_v    <- tv_gradient(v_range[v], theta)
    tv_var[v] <- as.vector(t(grad_v) %*% vcov %*% grad_v)

    ## compute zscore
    zscore[v] <- tv[v] / sqrt(tv_var[v])

    ## compute point-wise confidence interval with alpha level
    Uv[v] <- tv[v] + qnorm(1 - alpha) * sqrt(tv_var[v])
    Lv[v] <- tv[v] - qnorm(1 - alpha) * sqrt(tv_var[v])

    ## compute point-wise pvalues
    Upvalue[v] <- 1 - pnorm((threshold - tv[v]) / sqrt(tv_var[v]))
    Lpvalue[v] <- 1 - pnorm((tv[v] + threshold) / sqrt(tv_var[v]))
  }

  ## some summary quantities
  Umax <- max(Uv); Lmin <- min(Lv)
  tmax <- max(abs(tv))

  ## decide if we want to reject the null or not
  pvalue <- max(Upvalue, Lpvalue)
  reject <- (Umax < threshold) & (Lmin > -threshold)


  ## return
  return_list <- list(
    tv = tv, tv_var = tv_var, v_range = v_range,
    Uv = Uv, Lv = Lv,
    Umax = Umax, Lmin = Lmin, tmax = tmax,
    Upvalue = Upvalue, Lpvalue = Lpvalue, pvalue = pvalue,
    zscore = zscore, reject = reject
  )

  attr(return_list, "threshold") <- threshold
  class(return_list)             <- c("orddid", "orddid.test")
  return(return_list)
}



#' Selecting Equivalence Threshold
#'
#' A function to compute the data-dependent threshold for the equivalence test.
#' @param object a fitted object from \code{orddid}.
#' @return a scalar value of threshold.
#' @export
calc_threshold <- function(object, omega = 0.05) {
  if (!("orddid.fit" %in% class(object))) {
    stop("object should be the output of orddid function!")
  }

  ## check the sample size
  n  <- attr(object, "n")
  n1 <- attr(object, "n1")
  n0 <- n - n1

  const   <- sqrt(-log(omega) / 2)
  n_scale <- sqrt((n1 + n0) / (n1 * n0))
  delta   <- const * n_scale
  if (delta > 1.0) {
    warning("threshold is too large due to smalle sample size.")
    delta <- 1
  }

  return(delta)
}

## --------------------------------------------------------------- ##
##                     Auxiliary Functions                         ##
## --------------------------------------------------------------- ##

#' Compute qd(v)
#'
#' @return a scalar of qd(v) evaluated at u.
#' @importFrom VGAM erf
#' @keywords internal
qd <- function(u, mu1, mu0, s1, s0) {
  B <- erf(2 * u - 1, inverse = TRUE) / (s0 / s1)
  A <- erf((mu1 - mu0) / (s0 * sqrt(2)) + B)
  return((1 + A) / 2)
}


#' Compute tvmax
#'
#' Compute max_v | q1(v) - q0(v) | over 0 < v < 1.
#'
#' @param theta a list of parameters.
#' @return a scalar value of max_v |q1(v) - q0(v)|.
#' @keywords internal
tv_max <- function(theta, grid = 0.01) {
  v_range <- seq(0.001, 0.999, by = grid)
  tv <- rep(NA, length(v_range))
  for (v in 1:length(v_range)) {
    tv[v] <- qd(v_range[v], mu1 = theta$mu11, mu0 = theta$mu10,
               s1 = theta$sd11, s0 = theta$sd10) -
            qd(v_range[v], mu1 = theta$mu01, mu0 = theta$mu00,
                       s1 = theta$sd01, s0 = theta$sd00)
  }

  tv_max <- max(tv)
  attr(tv_max, "tv") <- tv
  return(tv_max)
}


#' Compute zd(v)
#'
#' @importFrom VGAM erf
#' @keywords internal
zd <- function(v, mu1, mu0, s1, s0) {
  A <- (mu1 - mu0) / (s0 * sqrt(2))
  B <- erf(2 * v - 1, inverse = TRUE) / (s0 / s1)
  return(A + B)
}



#' Compute gradient of t(v; theta)
#' @param theta a list of parameters.
#' @return a gradient vector.
#' @importFrom VGAM erf
#' @keywords internal
tv_gradient <- function(v, theta) {
  ## compute intermediate quantities
  z0 <- zd(v, mu1 = theta$mu01, mu0 = theta$mu00, s1 = theta$sd01, s0 = theta$sd00)
  z1 <- zd(v, mu1 = theta$mu11, mu0 = theta$mu10, s1 = theta$sd11, s0 = theta$sd10)
  exp_z0   <- exp(-z0^2)
  exp_z1   <- exp(-z1^2)
  erf_1    <- erf(2 * v - 1, inverse = TRUE)
  sqrt_pi  <- sqrt(pi)
  sqrt_2pi <- sqrt(2 * pi)

  ## compute the gradient
  grad <- c(
    ## derivative wrt theta 0
    exp_z0 / (sqrt_2pi * theta$sd00),
    exp_z0 * z0  / (sqrt_pi * theta$sd00),
    -exp_z0 / (sqrt_2pi * theta$sd00),
    -exp_z0 * erf_1 / (sqrt_pi * theta$sd00),
    ## wrt theta 1
    -exp_z1 / (sqrt_2pi * theta$sd10),
    -exp_z1 * z1  / (sqrt_pi * theta$sd10),
    exp_z1 / (sqrt_2pi * theta$sd10),
    exp_z1 * erf_1 / (sqrt_pi * theta$sd10)
  )

  return(grad)
}
