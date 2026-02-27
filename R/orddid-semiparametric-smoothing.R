## ============================================================
## orddid-semiparametric-smoothing.R
## I-spline smoothing of step-function CDF from PAVA or NPMLE
## ============================================================

#' Smooth an isotonic step-function CDF using I-splines
#'
#' @param Fhat  FhatIsotonic object (or list with $knots, $values)
#' @param n0    sample size of the (D=0,t=0) cell (for choosing K)
#' @param K     number of internal knots; if NULL, uses floor(n0^(1/3))
#' @param degree spline degree (default 3, giving C^2 smoothness)
#' @return list with Ftilde (smoothed CDF), ftilde (density), coef, knot_seq
#' @keywords internal
.smooth_Fhat <- function(Fhat, n0, K = NULL, degree = 3) {

  if (!requireNamespace("splines2", quietly = TRUE)) {
    stop("Package 'splines2' is required for semiparametric estimation. ",
         "Install it with install.packages('splines2').")
  }
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("Package 'quadprog' is required for semiparametric estimation. ",
         "Install it with install.packages('quadprog').")
  }

  grid_u <- Fhat$knots
  grid_F <- Fhat$values
  G      <- length(grid_u)

  if (is.null(K)) K <- max(3L, floor(n0^(1/3)))

  probs     <- seq(0, 1, length.out = K + 2)[-c(1, K + 2)]
  int_knots <- stats::quantile(grid_u, probs = probs)

  u_lo <- min(grid_u)
  u_hi <- max(grid_u)
  F_lo <- grid_F[1]
  F_hi <- grid_F[G]

  ## Degenerate case: single knot or no range
  if (G < 2 || (u_hi - u_lo) < 1e-12) {
    Ftilde <- function(u) rep(F_lo, length(u))
    ftilde <- function(u) rep(0, length(u))
    return(list(
      Ftilde   = Ftilde,
      ftilde   = ftilde,
      coef     = numeric(0),
      knot_seq = numeric(0),
      u_range  = c(u_lo, u_hi),
      F_range  = c(F_lo, F_hi)
    ))
  }

  B <- splines2::iSpline(grid_u, knots = int_knots,
                          degree = degree,
                          Boundary.knots = c(u_lo, u_hi),
                          intercept = TRUE)
  n_basis <- ncol(B)

  target <- grid_F - F_lo

  Dmat <- crossprod(B)
  dvec <- crossprod(B, target)
  diag(Dmat) <- diag(Dmat) + 1e-8

  Amat <- diag(n_basis)
  bvec <- rep(0, n_basis)

  sol   <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 0)
  c_hat <- pmax(sol$solution, 0)

  ## Precompute CDF and density on a fine grid for fast O(log n) lookup
  ## via approxfun, avoiding repeated spline basis construction
  n_grid <- 10000L
  grid_u <- seq(u_lo, u_hi, length.out = n_grid)

  B_grid <- splines2::iSpline(grid_u, knots = int_knots,
                               degree = degree,
                               Boundary.knots = c(u_lo, u_hi),
                               intercept = TRUE)
  F_grid <- F_lo + as.numeric(B_grid %*% c_hat)
  F_grid <- pmin(pmax(F_grid, 0), 1)
  ## Enforce monotonicity on grid
  for (k in 2:n_grid) F_grid[k] <- max(F_grid[k], F_grid[k - 1])

  M_grid <- splines2::mSpline(grid_u, knots = int_knots,
                               degree = degree,
                               Boundary.knots = c(u_lo, u_hi),
                               intercept = TRUE)
  f_grid <- pmax(as.numeric(M_grid %*% c_hat), 0)

  Ftilde <- stats::approxfun(grid_u, F_grid, rule = 2, ties = mean)
  ftilde <- stats::approxfun(grid_u, f_grid, rule = 2, ties = mean)

  list(
    Ftilde   = Ftilde,
    ftilde   = ftilde,
    coef     = c_hat,
    knot_seq = int_knots,
    u_range  = c(u_lo, u_hi),
    F_range  = c(F_lo, F_hi)
  )
}
