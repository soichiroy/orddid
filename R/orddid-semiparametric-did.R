## ============================================================
## orddid-semiparametric-did.R
## Canonical three-category semiparametric DID estimator
## ============================================================

#' Full semiparametric DID estimation for one sample
#'
#' The first stage uses an all-category ordered-probit initializer followed by
#' two interval-censoring NPMLE, coefficient-moment, and cutoff-moment updates,
#' as recommended for the implementation in the paper. A direct I-spline is
#' used as a common numerical approximation when
#' solving the first-stage moments, so the first-stage parameters are shared by
#' both reported regularizations. A final raw NPMLE is recomputed at the updated
#' parameters and then regularized by either an I-spline or the paper's
#' vanishing ridge before the two profile-likelihood fits.
#'
#' @param Y00,X00 outcomes and covariates for control pre
#' @param Y01,X01 outcomes and covariates for control post
#' @param Y10,X10 outcomes and covariates for treated pre
#' @param Y11,X11 outcomes and covariates for treated post
#' @param fixed_index coefficient normalized to one in the internal `b=-beta`
#'   convention
#' @param verbose print first-stage progress
#' @param tau_assumption relative-effect bound assumptions
#' @param semiparametric_regularization final NPMLE regularization
#' @param regularization_n independent-unit count used for ridge epsilon=1/n
#' @return estimates, effects, and first/second-stage diagnostics
#' @keywords internal
.semiparametric_did_once <- function(
    Y00, X00, Y01, X01, Y10, X10, Y11, X11,
    fixed_index = 2L, verbose = FALSE, tau_assumption = "none",
    semiparametric_regularization = c("ispline", "ridge"),
    regularization_n = NULL) {
  semiparametric_regularization <- match.arg(
    semiparametric_regularization
  )
  matrices <- lapply(list(X00, X01, X10, X11), as.matrix)
  X00 <- matrices[[1L]]
  X01 <- matrices[[2L]]
  X10 <- matrices[[3L]]
  X11 <- matrices[[4L]]
  outcomes <- lapply(list(Y00, Y01, Y10, Y11), as.integer)
  Y00 <- outcomes[[1L]]
  Y01 <- outcomes[[2L]]
  Y10 <- outcomes[[3L]]
  Y11 <- outcomes[[4L]]

  if (any(vapply(outcomes, anyNA, logical(1))) ||
      !identical(sort(unique(c(Y00, Y01, Y10, Y11))), 0:2) ||
      !identical(sort(unique(Y00)), 0:2)) {
    stop("Semiparametric DID requires categories 0:2, including all three in Y00.")
  }
  p <- ncol(X00)
  if (p < 1L || any(vapply(matrices, ncol, integer(1)) != p) ||
      any(vapply(seq_along(matrices), function(i) {
        nrow(matrices[[i]]) != length(outcomes[[i]]) ||
          any(!is.finite(matrices[[i]]))
      }, logical(1)))) {
    stop("All design matrices must be finite, have common columns, and match Y.")
  }
  if (is.null(regularization_n)) {
    regularization_n <- length(Y00) + length(Y10)
  }

  first_stage <- .run_semiparametric_first_stage(
    Y00 = Y00,
    X00 = X00,
    fixed_index = fixed_index,
    max_iter = 2L,
    verbose = verbose
  )
  if (first_stage$stop_reason == "max_iter") {
    first_stage$stop_reason <- "two_updates"
  }
  beta_00 <- first_stage$beta
  alpha <- first_stage$alpha
  kappa <- c(-Inf, 0, alpha, Inf)
  regularized <- .regularize_Fhat(
    first_stage$npmle,
    n0 = length(Y00),
    method = semiparametric_regularization,
    regularization_n = regularization_n
  )

  ## Both regularizers can be nonsmooth: ridge retains the NPMLE jumps and the
  ## direct I-spline uses finite 0/1 endpoint tails. Numeric BFGS matches the
  ## paper simulation implementation and avoids an invalid global derivative.
  ## Since BFGS and Nelder-Mead reach different basins for some NPMLE samples,
  ## evaluate both from the same deterministic starts and retain the best
  ## converged likelihood.
  optimizer <- "both"
  profile_density <- NULL
  fit01 <- .step2_profile_likelihood(
    Y01, X01, regularized$Ftilde, kappa,
    beta_init = beta_00,
    ftilde = profile_density,
    optimizer = optimizer
  )
  fit10 <- .step2_profile_likelihood(
    Y10, X10, regularized$Ftilde, kappa,
    beta_init = beta_00,
    ftilde = profile_density,
    optimizer = optimizer
  )

  effects <- .compute_semiparametric_effects(
    beta_00 = beta_00,
    beta_01 = fit01$beta,
    xi_01 = fit01$xi,
    beta_10 = fit10$beta,
    xi_10 = fit10$xi,
    kappa = kappa,
    Ftilde = regularized$Ftilde,
    X11 = X11,
    Y11 = Y11,
    tau_assumption = tau_assumption
  )

  list(
    beta_00 = beta_00,
    alpha = alpha,
    beta_01 = fit01$beta,
    xi_01 = fit01$xi,
    beta_10 = fit10$beta,
    xi_10 = fit10$xi,
    zeta = effects$zeta,
    tau = effects$tau,
    obs_probs = effects$obs_probs,
    cf_probs_avg = effects$cf_probs_avg,
    convergence = c(
      first_stage = if (first_stage$completed) 0L else 1L,
      fit01 = fit01$convergence,
      fit10 = fit10$convergence
    ),
    semiparametric_regularization = semiparametric_regularization,
    regularization = regularized,
    first_stage = first_stage,
    profile_fits = list(control_post = fit01, treated_pre = fit10)
  )
}
