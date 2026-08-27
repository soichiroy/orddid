## ============================================================
## orddid-semiparametric-smoothing.R
## Regularization of the interval-censoring NPMLE
## ============================================================

#' Prepare a step-function CDF for regularization
#'
#' @param Fhat list with numeric components `knots` and `values`
#' @return list with sorted, unique knots and monotone values
#' @keywords internal
.prepare_step_cdf <- function(Fhat) {
  if (is.null(Fhat$knots) || is.null(Fhat$values)) {
    stop("Fhat must contain 'knots' and 'values'.")
  }

  knots <- as.numeric(Fhat$knots)
  values <- as.numeric(Fhat$values)
  keep <- is.finite(knots) & is.finite(values)
  knots <- knots[keep]
  values <- values[keep]
  if (length(knots) < 2L) {
    stop("The NPMLE must contain at least two finite support points.")
  }

  ord <- order(knots, values)
  knots <- knots[ord]
  values <- values[ord]
  groups <- cumsum(c(TRUE, diff(knots) != 0))
  split_knots <- split(knots, groups)
  split_values <- split(values, groups)
  knots <- vapply(split_knots, function(x) x[1L], numeric(1))
  values <- vapply(split_values, max, numeric(1))
  values <- cummax(pmin(pmax(values, 0), 1))

  if (length(knots) < 2L || diff(range(knots)) <= sqrt(.Machine$double.eps)) {
    stop("The NPMLE support is degenerate and cannot be regularized.")
  }
  list(knots = knots, values = values)
}

#' Smooth a step-function NPMLE using an I-spline
#'
#' The returned functions evaluate the spline basis directly. This avoids the
#' optimizer-dependent approximation introduced by interpolating a fixed grid.
#'
#' @param Fhat list with numeric components `knots` and `values`
#' @param n0 sample size of the control-pre cell
#' @param K number of requested internal knots; defaults to floor(n0^(1/3))
#' @param degree spline degree
#' @return list containing CDF and density functions and fit details
#' @keywords internal
.smooth_Fhat <- function(Fhat, n0, K = NULL, degree = 3) {
  step <- .prepare_step_cdf(Fhat)
  grid_u <- step$knots
  grid_F <- step$values
  u_lo <- min(grid_u)
  u_hi <- max(grid_u)
  F_lo <- grid_F[1L]

  if (is.null(K)) K <- max(3L, floor(n0^(1 / 3)))
  probs <- seq(0, 1, length.out = K + 2L)[-c(1L, K + 2L)]
  int_knots <- unique(as.numeric(stats::quantile(
    grid_u, probs = probs, names = FALSE, type = 7
  )))
  tolerance <- sqrt(.Machine$double.eps) * max(1, abs(c(u_lo, u_hi)))
  int_knots <- int_knots[
    int_knots > u_lo + tolerance & int_knots < u_hi - tolerance
  ]
  spline_knots <- if (length(int_knots)) int_knots else NULL

  basis <- splines2::iSpline(
    grid_u,
    knots = spline_knots,
    degree = degree,
    Boundary.knots = c(u_lo, u_hi),
    intercept = TRUE
  )
  n_basis <- ncol(basis)
  target <- grid_F - F_lo
  Dmat <- crossprod(basis)
  diag(Dmat) <- diag(Dmat) + 1e-8
  upper_basis <- as.numeric(splines2::iSpline(
    u_hi,
    knots = spline_knots,
    degree = degree,
    Boundary.knots = c(u_lo, u_hi),
    intercept = TRUE
  ))
  solution <- quadprog::solve.QP(
    Dmat = Dmat,
    dvec = as.numeric(crossprod(basis, target)),
    Amat = cbind(diag(n_basis), -upper_basis),
    bvec = c(rep(0, n_basis), -(1 - F_lo)),
    meq = 0
  )
  coefficients <- pmax(solution$solution, 0)
  fitted_increment <- sum(upper_basis * coefficients)
  allowed_increment <- max(0, 1 - F_lo)
  if (fitted_increment > allowed_increment && fitted_increment > 0) {
    coefficients <- coefficients * allowed_increment / fitted_increment
  }

  Ftilde <- function(u) {
    if (!length(u)) return(numeric(0))
    below <- u < u_lo
    above <- u > u_hi
    u_eval <- pmin(pmax(u, u_lo), u_hi)
    new_basis <- splines2::iSpline(
      u_eval,
      knots = spline_knots,
      degree = degree,
      Boundary.knots = c(u_lo, u_hi),
      intercept = TRUE
    )
    value <- F_lo + as.numeric(new_basis %*% coefficients)
    value[below] <- 0
    value[above] <- 1
    pmin(pmax(value, 0), 1)
  }

  ftilde <- function(u) {
    if (!length(u)) return(numeric(0))
    inside <- is.finite(u) & u >= u_lo & u <= u_hi
    value <- numeric(length(u))
    if (any(inside)) {
      new_basis <- splines2::mSpline(
        u[inside],
        knots = spline_knots,
        degree = degree,
        Boundary.knots = c(u_lo, u_hi),
        intercept = TRUE
      )
      value[inside] <- pmax(as.numeric(new_basis %*% coefficients), 0)
    }
    value
  }

  list(
    Ftilde = Ftilde,
    ftilde = ftilde,
    coef = coefficients,
    knot_seq = int_knots,
    u_range = c(u_lo, u_hi),
    F_range = range(Ftilde(c(-Inf, Inf))),
    method = "ispline"
  )
}

#' Apply the paper's vanishing-ridge regularization to an NPMLE
#'
#' @param Fhat list with numeric components `knots` and `values`
#' @param regularization_n number of independent units defining epsilon = 1/n
#' @return list containing the regularized CDF and fit details
#' @keywords internal
.smooth_Fhat_ridge <- function(Fhat, regularization_n) {
  step <- .prepare_step_cdf(Fhat)
  if (length(regularization_n) != 1L || !is.finite(regularization_n) ||
      regularization_n <= 2) {
    stop("regularization_n must be a finite number greater than two.")
  }

  u_lo <- min(step$knots)
  u_hi <- max(step$knots)
  u_span <- u_hi - u_lo
  epsilon <- 1 / regularization_n
  F_step_inside <- stats::approxfun(
    step$knots,
    step$values,
    method = "constant",
    rule = 2,
    f = 0,
    ties = max
  )
  F_step <- function(u) {
    value <- F_step_inside(pmin(pmax(u, u_lo), u_hi))
    value[u < u_lo] <- 0
    value[u > u_hi] <- 1
    value
  }

  Ftilde <- function(u) {
    ramp <- pmin(pmax((u - u_lo) / u_span, 0), 1)
    value <- (1 - 2 * epsilon) * F_step(u) + epsilon + epsilon * ramp
    pmin(pmax(value, epsilon), 1)
  }

  list(
    Ftilde = Ftilde,
    ftilde = NULL,
    eps = epsilon,
    u_range = c(u_lo, u_hi),
    F_range = range(Ftilde(c(-Inf, Inf))),
    method = "ridge"
  )
}

#' Select the second-stage NPMLE regularization
#' @keywords internal
.regularize_Fhat <- function(Fhat, n0, method = c("ispline", "ridge"),
                             regularization_n = n0) {
  method <- match.arg(method)
  if (method == "ispline") {
    .smooth_Fhat(Fhat, n0 = n0)
  } else {
    .smooth_Fhat_ridge(Fhat, regularization_n = regularization_n)
  }
}
