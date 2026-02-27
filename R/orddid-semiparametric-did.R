## ============================================================
## orddid-semiparametric-did.R
## Main semiparametric DID workhorse (NPMLE method)
## ============================================================

#' Full semiparametric DID estimation (single run, NPMLE method)
#'
#' Pipeline: probit init -> PAVA -> interval-censoring NPMLE ->
#'           I-spline smooth -> profile likelihood -> counterfactual
#'
#' @param Y00,X00 outcome and covariates for control pre
#' @param Y01,X01 outcome and covariates for control post
#' @param Y10,X10 outcome and covariates for treated pre
#' @param Y11,X11 outcome and covariates for treated post
#' @param fixed_index which coefficient to normalize in Step 1
#' @param verbose print progress
#' @return list with all estimates and treatment effects
#' @keywords internal
.semiparametric_did_once <- function(Y00, X00, Y01, X01,
                                     Y10, X10, Y11, X11,
                                     fixed_index = 2,
                                     verbose = FALSE) {

  n0 <- length(Y00)

  ## ---- Step 1: probit init + PAVA + NPMLE ----
  b_00 <- .init_beta_probit(Y00, X00, fixed_index)
  if (verbose) cat("Step 1: b_00 =", round(b_00, 3), "\n")

  Fhat    <- stage1i_isotonic(Y00, X00, b_00)
  beta_00 <- -b_00

  ## Initial smooth + alpha from PAVA
  sm    <- .smooth_Fhat(Fhat, n0 = n0)
  alpha <- .estimate_alpha_smooth(Y00, X00, b_00, sm$Ftilde)
  if (is.na(alpha)) {
    alpha_fit <- EstimateCutoff(Y00, X00, b_00, Fhat = Fhat, verbose = FALSE)
    alpha     <- alpha_fit$alpha
  }

  ## Iterate: NPMLE using all categories, then re-estimate alpha
  Ftilde  <- sm$Ftilde
  sm_last <- sm
  for (iter in 1:2) {
    npmle <- .compute_joint_npmle(Y00, X00, b_00, alpha)
    if (is.null(npmle)) break

    npmle_fhat <- structure(
      list(knots = npmle$knots, values = npmle$values),
      class = "FhatIsotonic"
    )
    sm_last  <- .smooth_Fhat(npmle_fhat, n0 = n0)
    Ftilde   <- sm_last$Ftilde

    alpha_new <- .estimate_alpha_smooth(Y00, X00, b_00, Ftilde)
    if (!is.na(alpha_new)) alpha <- alpha_new
  }

  if (verbose) {
    cat("Step 1: beta_00 =", round(beta_00, 3),
        ", alpha =", round(alpha, 3), "\n")
  }

  ## Build kappa vector
  kappa <- c(-Inf, 0, alpha, Inf)

  ## ---- Step 2: profile likelihood for (D=0,t=1) and (D=1,t=0) ----
  ## NOTE: analytical gradient disabled — the piecewise-linear approxfun
  ## density causes BFGS to converge to worse local optima than

  ## finite-difference BFGS.  Use ftilde = NULL for now.
  fit01 <- .step2_profile_likelihood(Y01, X01, Ftilde, kappa,
                                     beta_init = beta_00,
                                     ftilde = NULL)
  fit10 <- .step2_profile_likelihood(Y10, X10, Ftilde, kappa,
                                     beta_init = beta_00,
                                     ftilde = NULL)

  ## ---- Step 3: counterfactual and treatment effects ----
  te <- .compute_semiparametric_effects(
    beta_00  = beta_00,
    beta_01  = fit01$beta,
    xi_01    = fit01$xi,
    beta_10  = fit10$beta,
    xi_10    = fit10$xi,
    kappa    = kappa,
    Ftilde   = Ftilde,
    X11      = X11,
    Y11      = Y11
  )

  list(
    beta_00 = beta_00,
    alpha   = alpha,
    beta_01 = fit01$beta,
    xi_01   = fit01$xi,
    beta_10 = fit10$beta,
    xi_10   = fit10$xi,
    zeta    = te$zeta,
    tau     = te$tau,
    obs_probs    = te$obs_probs,
    cf_probs_avg = te$cf_probs_avg,
    convergence  = c(fit01 = fit01$convergence,
                     fit10 = fit10$convergence)
  )
}
