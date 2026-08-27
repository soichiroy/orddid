## ============================================================
## orddid-semiparametric-npmle.R
## Interval-censoring NPMLE for the three-category model
## ============================================================

#' Compute the joint NPMLE using all three ordinal categories
#'
#' Under the internal convention `b = -beta`, the error intervals are
#' `(-Inf, Xb]`, `(Xb, Xb + alpha]`, and `(Xb + alpha, Inf)`.
#'
#' @param Y00 integer outcomes coded 0, 1, 2
#' @param X00 design matrix
#' @param b_00 normalized index coefficients
#' @param alpha positive second cutoff
#' @return list with sorted `knots`, monotone `values`, and the fitted object
#' @keywords internal
.compute_joint_npmle <- function(Y00, X00, b_00, alpha) {
  Y00 <- as.integer(Y00)
  X00 <- as.matrix(X00)
  if (anyNA(Y00) || !identical(sort(unique(Y00)), 0:2)) {
    stop("Y00 must contain exactly the three categories 0, 1, and 2.")
  }
  if (nrow(X00) != length(Y00) || ncol(X00) != length(b_00) ||
      any(!is.finite(X00)) || any(!is.finite(b_00))) {
    stop("Y00, X00, and b_00 have incompatible dimensions.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0) {
    stop("alpha must be one positive finite cutoff.")
  }

  u <- as.numeric(X00 %*% b_00)
  lower <- ifelse(Y00 == 0L, -Inf, ifelse(Y00 == 1L, u, u + alpha))
  upper <- ifelse(Y00 == 0L, u, ifelse(Y00 == 1L, u + alpha, Inf))
  interval_data <- data.frame(lower = lower, upper = upper)

  fit <- tryCatch(
    icenReg::ic_np(cbind(lower, upper) ~ 0, data = interval_data),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  curves <- icenReg::getSCurves(fit)
  raw <- list(
    knots = curves$Tbull_ints[, 1L],
    values = 1 - curves$S_curves$baseline
  )
  step <- tryCatch(.prepare_step_cdf(raw), error = function(e) NULL)
  if (is.null(step)) return(NULL)

  Fhat_inside <- stats::approxfun(
    step$knots,
    step$values,
    method = "constant",
    rule = 2,
    f = 0,
    ties = max
  )
  Fhat <- function(u) {
    value <- Fhat_inside(pmin(pmax(u, min(step$knots)), max(step$knots)))
    value[u < min(step$knots)] <- 0
    value[u > max(step$knots)] <- 1
    value
  }
  list(
    Ftilde = Fhat,
    knots = step$knots,
    values = step$values,
    fit = fit
  )
}
