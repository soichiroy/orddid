## ============================================================
## orddid-parametric-covariates.R
## Parametric (ordered probit) DID with covariates
##
## Uses F = Phi (known normal CDF) throughout.
## Step 1: ordered probit MLE on (D=0,t=0) -> beta_00, kappa.
## Step 2: profile likelihood with F=Phi for other cells.
## Step 3: counterfactual via identification formula with F=Phi.
## ============================================================

#' Ordered probit MLE with covariates
#'
#' Fits P(Y <= j | X) = Phi(kappa_j - X'beta) with sigma = 1.
#' Normalization: kappa_1 = 0, remaining cutoff gaps on log scale.
#'
#' @param Y integer outcome vector (0, 1, ..., J-1)
#' @param X n x p design matrix (with intercept column)
#' @return list with beta, kappa, alpha, loglik, convergence
#' @keywords internal
.ordered_probit_mle <- function(Y, X) {

  n <- length(Y)
  p <- ncol(X)
  cats <- sort(unique(Y))
  J <- length(cats)

  ## Category indicator matrix
  ind <- matrix(0L, n, J)
  for (j in seq_len(J)) ind[, j] <- as.integer(Y == cats[j])

  ## Free cutoff parameters: kappa = (-Inf, 0, cumsum(exp(delta)), Inf)
  n_delta <- J - 2L

  ## Initialize beta from binary probit (Y=0 vs Y>0)
  y_bin <- as.integer(Y == cats[1])
  beta_init <- tryCatch({
    fit <- stats::glm(y_bin ~ X - 1, family = stats::binomial(link = "probit"))
    -as.numeric(stats::coef(fit))
  }, error = function(e) rep(0, p))

  delta_init <- if (n_delta > 0) rep(0, n_delta) else numeric(0)
  theta_init <- c(beta_init, delta_init)

  ## Cache for combined fn+gradient evaluation
  .cache <- new.env(parent = emptyenv())
  .cache$theta <- NULL

  .eval_at <- function(theta) {
    if (identical(theta, .cache$theta)) return(invisible(NULL))

    beta <- theta[1:p]
    if (n_delta > 0) {
      delta <- theta[(p + 1):(p + n_delta)]
      kappa_inner <- cumsum(exp(delta))
    } else {
      kappa_inner <- numeric(0)
      delta <- numeric(0)
    }
    kappa <- c(-Inf, 0, kappa_inner, Inf)

    mu <- as.numeric(X %*% beta)
    ll <- 0
    dll_dbeta <- rep(0, p)
    dll_dkappa <- if (n_delta > 0) rep(0, n_delta) else numeric(0)

    for (j in seq_len(J)) {
      if (is.finite(kappa[j + 1])) {
        F_hi <- stats::pnorm(kappa[j + 1] - mu)
        f_hi <- stats::dnorm(kappa[j + 1] - mu)
      } else {
        F_hi <- if (kappa[j + 1] > 0) rep(1, n) else rep(0, n)
        f_hi <- rep(0, n)
      }
      if (is.finite(kappa[j])) {
        F_lo <- stats::pnorm(kappa[j] - mu)
        f_lo <- stats::dnorm(kappa[j] - mu)
      } else {
        F_lo <- if (kappa[j] < 0) rep(0, n) else rep(1, n)
        f_lo <- rep(0, n)
      }

      prob_j <- pmax(F_hi - F_lo, 1e-12)
      ll <- ll + sum(ind[, j] * log(prob_j))

      w_j <- ind[, j] / prob_j

      ## Gradient w.r.t. beta
      dll_dbeta <- dll_dbeta +
        as.numeric(crossprod(X, w_j * (f_lo - f_hi)))

      ## Gradient w.r.t. interior cutoffs (kappa_inner[m] = kappa[m+2])
      if (n_delta > 0) {
        ## kappa_inner[m] is upper bound for category j = m+1
        m_hi <- j - 1L
        if (m_hi >= 1 && m_hi <= n_delta) {
          dll_dkappa[m_hi] <- dll_dkappa[m_hi] + sum(w_j * f_hi)
        }
        ## kappa_inner[m] is lower bound for category j = m+2
        m_lo <- j - 2L
        if (m_lo >= 1 && m_lo <= n_delta) {
          dll_dkappa[m_lo] <- dll_dkappa[m_lo] - sum(w_j * f_lo)
        }
      }
    }

    ## Chain rule: kappa_inner[m] = sum(exp(delta[1:m]))
    ## d(ll)/d(delta[k]) = exp(delta[k]) * sum_{m >= k} d(ll)/d(kappa_inner[m])
    if (n_delta > 0) {
      cum_dk <- rev(cumsum(rev(dll_dkappa)))
      dll_ddelta <- exp(delta) * cum_dk
    } else {
      dll_ddelta <- numeric(0)
    }

    .cache$theta <- theta
    .cache$val   <- -ll
    .cache$grad  <- -c(dll_dbeta, dll_ddelta)
  }

  negloglik <- function(theta) { .eval_at(theta); .cache$val }
  grad_fn   <- function(theta) { .eval_at(theta); .cache$grad }

  fit <- stats::optim(theta_init, negloglik, gr = grad_fn,
                      method = "BFGS",
                      control = list(maxit = 500, reltol = 1e-8))

  beta <- fit$par[1:p]
  if (n_delta > 0) {
    kappa_inner <- cumsum(exp(fit$par[(p + 1):(p + n_delta)]))
  } else {
    kappa_inner <- numeric(0)
  }
  kappa <- c(-Inf, 0, kappa_inner, Inf)

  list(
    beta        = beta,
    kappa       = kappa,
    alpha       = if (n_delta > 0) kappa_inner[1] else NA_real_,
    loglik      = -fit$value,
    convergence = fit$convergence
  )
}


#' Full parametric DID with covariates (single run)
#'
#' Pipeline: ordered probit MLE -> profile likelihood (F=Phi) -> counterfactual.
#' Reuses .step2_profile_likelihood() and .compute_semiparametric_effects()
#' with Ftilde = pnorm and ftilde = dnorm.
#'
#' @param Y00,X00 outcome and covariates for control pre
#' @param Y01,X01 outcome and covariates for control post
#' @param Y10,X10 outcome and covariates for treated pre
#' @param Y11,X11 outcome and covariates for treated post
#' @param fixed_index ignored (included for interface compatibility)
#' @param verbose print progress
#' @return list with all estimates and treatment effects
#' @keywords internal
.parametric_covariates_did_once <- function(Y00, X00, Y01, X01,
                                            Y10, X10, Y11, X11,
                                            fixed_index = 2,
                                            verbose = FALSE) {

  ## Step 1: ordered probit MLE on (D=0,t=0) cell
  fit00 <- .ordered_probit_mle(Y00, X00)
  beta_00 <- fit00$beta
  kappa   <- fit00$kappa

  if (verbose) {
    cat("Step 1: beta_00 =", round(beta_00, 3),
        ", alpha =", round(fit00$alpha, 3), "\n")
  }

  ## Step 2: profile likelihood with F = Phi for other cells
  fit01 <- .step2_profile_likelihood(Y01, X01, stats::pnorm, kappa,
                                     beta_init = beta_00,
                                     ftilde = stats::dnorm)
  fit10 <- .step2_profile_likelihood(Y10, X10, stats::pnorm, kappa,
                                     beta_init = beta_00,
                                     ftilde = stats::dnorm)

  if (verbose) {
    cat("Step 2 (01): beta =", round(fit01$beta, 3),
        ", xi =", round(fit01$xi, 3), "\n")
    cat("Step 2 (10): beta =", round(fit10$beta, 3),
        ", xi =", round(fit10$xi, 3), "\n")
  }

  ## Step 3: counterfactual and treatment effects with F = Phi
  te <- .compute_semiparametric_effects(
    beta_00  = beta_00,
    beta_01  = fit01$beta,
    xi_01    = fit01$xi,
    beta_10  = fit10$beta,
    xi_10    = fit10$xi,
    kappa    = kappa,
    Ftilde   = stats::pnorm,
    X11      = X11,
    Y11      = Y11
  )

  list(
    beta_00 = beta_00,
    alpha   = fit00$alpha,
    beta_01 = fit01$beta,
    xi_01   = fit01$xi,
    beta_10 = fit10$beta,
    xi_10   = fit10$xi,
    zeta    = te$zeta,
    tau     = te$tau,
    obs_probs    = te$obs_probs,
    cf_probs_avg = te$cf_probs_avg,
    convergence  = c(fit00 = fit00$convergence,
                     fit01 = fit01$convergence,
                     fit10 = fit10$convergence)
  )
}
