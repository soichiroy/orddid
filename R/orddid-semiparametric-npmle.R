## ============================================================
## orddid-semiparametric-npmle.R
## Interval-censoring NPMLE for error distribution + smooth alpha
## ============================================================

#' Estimate cutoff alpha using the smoothed Ftilde
#'
#' Solves: mean(1 - 1\{Y = max(Y)\} - Ftilde(X'b + alpha)) = 0
#'
#' @param Y00    outcomes for (D=0, t=0) cell
#' @param X00    design matrix
#' @param b_00   Step 1 index coefficients
#' @param Ftilde smoothed CDF function
#' @return estimated alpha, or NA if root-finding fails
#' @keywords internal
.estimate_alpha_smooth <- function(Y00, X00, b_00, Ftilde) {
  u     <- as.numeric(X00 %*% b_00)
  y_max <- max(Y00)
  y3    <- as.integer(Y00 == y_max)

  psi <- function(a) mean(1 - y3 - Ftilde(u + a))

  lo <- 0
  hi <- 2 * diff(range(u))

  for (i in 1:20) {
    if (psi(hi) <= 0) break
    hi <- hi * 2
  }
  if (psi(hi) > 0) return(NA)

  tryCatch(
    stats::uniroot(psi, lower = lo, upper = hi, tol = 1e-8)$root,
    error = function(e) NA
  )
}

#' Compute joint NPMLE for error distribution using all ordinal categories
#'
#' Treats ordered response as interval-censored data:
#'   Y = min: epsilon in (-Inf, u_i]
#'   Y = middle: epsilon in (u_i, u_i+alpha]
#'   Y = max: epsilon in (u_i+alpha, Inf)
#'
#' @param Y00   outcomes for (D=0, t=0) cell
#' @param X00   design matrix
#' @param b_00  Step 1 index coefficients
#' @param alpha cutoff parameter
#' @return list with Ftilde, knots, values; or NULL on failure
#' @keywords internal
.compute_joint_npmle <- function(Y00, X00, b_00, alpha) {

  if (!requireNamespace("icenReg", quietly = TRUE)) {
    stop("Package 'icenReg' is required for NPMLE estimation. ",
         "Install it with install.packages('icenReg').")
  }

  u     <- as.numeric(X00 %*% b_00)
  n     <- length(Y00)
  y_min <- min(Y00)
  y_max <- max(Y00)

  L <- R <- numeric(n)
  for (i in seq_len(n)) {
    if (Y00[i] == y_min) {
      L[i] <- -Inf; R[i] <- u[i]
    } else if (Y00[i] == y_max) {
      L[i] <- u[i] + alpha; R[i] <- Inf
    } else {
      L[i] <- u[i]; R[i] <- u[i] + alpha
    }
  }

  df  <- data.frame(L = L, R = R)
  fit <- tryCatch(icenReg::ic_np(cbind(L, R) ~ 0, data = df),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  sc        <- icenReg::getSCurves(fit)
  t_vals    <- sc$Tbull_ints[, 1]
  S_vals    <- sc$S_curves$baseline
  finite_idx <- is.finite(t_vals)
  t_fin     <- t_vals[finite_idx]
  F_fin     <- pmin(pmax((1 - S_vals)[finite_idx], 0), 1)
  if (length(t_fin) < 2) return(NULL)

  Ftilde <- stats::approxfun(t_fin, F_fin, method = "constant",
                              rule = 2, f = 0, ties = mean)
  list(Ftilde = Ftilde, knots = t_fin, values = F_fin)
}
