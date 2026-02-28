## ============================================================
## orddid-semiparametric-counterfactual.R
## Counterfactual distribution and treatment effects
## ============================================================

#' Compute counterfactual probabilities and treatment effects
#'
#' Uses identification result:
#'   mu_11(x)    = x'beta_10 + exp(x'xi_10) * x'(beta_01 - beta_00)
#'   sigma_11(x) = exp(x'(xi_10 + xi_01))
#'
#' @param beta_00 coefficient vector from Step 1 (control pre)
#' @param beta_01 coefficient vector from Step 2 (control post)
#' @param xi_01   scale coefficient vector (control post)
#' @param beta_10 coefficient vector from Step 2 (treated pre)
#' @param xi_10   scale coefficient vector (treated pre)
#' @param kappa   cutoff vector c(-Inf, 0, alpha, Inf)
#' @param Ftilde  smoothed CDF function
#' @param X11     design matrix for treated post-treatment units
#' @param Y11     observed outcomes for treated post-treatment units
#' @return list with zeta, tau, obs_probs, cf_probs_avg
#' @keywords internal
.compute_semiparametric_effects <- function(beta_00, beta_01, xi_01,
                                            beta_10, xi_10,
                                            kappa, Ftilde,
                                            X11, Y11) {

  n1 <- nrow(X11)
  J  <- length(kappa) - 1

  mu_11    <- as.numeric(X11 %*% beta_10) +
              exp(as.numeric(X11 %*% xi_10)) *
              as.numeric(X11 %*% (beta_01 - beta_00))

  sigma_11 <- exp(as.numeric(X11 %*% (xi_10 + xi_01)))

  cf_probs <- matrix(0, n1, J)
  for (j in seq_len(J)) {
    u_hi <- (kappa[j + 1] - mu_11) / sigma_11
    u_lo <- (kappa[j]     - mu_11) / sigma_11

    F_hi <- ifelse(is.finite(u_hi), Ftilde(u_hi), ifelse(u_hi > 0, 1, 0))
    F_lo <- ifelse(is.finite(u_lo), Ftilde(u_lo), ifelse(u_lo > 0, 1, 0))

    cf_probs[, j] <- pmax(F_hi - F_lo, 0)
  }
  row_sums <- rowSums(cf_probs)
  cf_probs <- cf_probs / row_sums

  obs_probs <- numeric(J)
  for (j in seq_len(J)) {
    obs_probs[j] <- mean(Y11 == (j - 1))
  }

  cf_avg <- colMeans(cf_probs)
  zeta   <- obs_probs - cf_avg

  tau <- .ComputeRelativeEffect(obs_probs, cf_avg)

  list(
    zeta          = zeta,
    tau           = tau,
    obs_probs     = obs_probs,
    cf_probs_avg  = cf_avg,
    cf_probs_unit = cf_probs,
    mu_11         = mu_11,
    sigma_11      = sigma_11
  )
}
